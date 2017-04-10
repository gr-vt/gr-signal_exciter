

#include "signal_ofdm.hpp"
#include <stdio.h>
#include <stdexcept>

Signal_OFDM::Signal_OFDM(size_t fftsize, size_t cp_len, size_t active_carriers, size_t syms_per_frame, 
                         bool pilot_per_frame, size_t pilot_count, size_t* pilot_locations, float backoff,
                         int mod_type, int mod_order, float mod_offset, int seed, bool add_sync,
                         float* symbol_taper, size_t sample_overlap, float* interp_taps, size_t tap_len, int interp,
                         bool enable, size_t buff_size, size_t min_notify)
  : d_fftsize(fftsize),
    d_cp_len(cp_len),
    d_active(active_carriers),
    d_spf(syms_per_frame),
    d_ppf(pilot_per_frame),
    d_npilots(pilot_count),
    d_backoff(backoff),//this is in dB
    d_sync(add_sync),
    d_protected(0),
    d_symbol_idx(0),
    d_first_pass(true),
    d_enable(enable),
    d_buffer_size(buff_size),
    d_interp(interp),
    d_samp_offset(0),
    d_branch_offset(0),
    d_frame_offset(0),
    d_samp_overlap(sample_overlap)
{
  boost::mutex::scoped_lock scoped_lock(fftw_lock());

  d_active_list = std::vector<size_t>(active_carriers,0);
  d_pilot_list = std::vector<size_t>(pilot_count,0);
  
  for(size_t idx = 0; idx < pilot_count; idx++){
    d_pilot_list[idx] = pilot_locations[idx];
  }

  sort_subcarriers();

  float pulseshape[1] = {1.};
  if(mod_type == gr::signal_exciter::QAM){
    d_mod = new Signal_QAM(mod_order,mod_offset,1,pulseshape,1,seed,enable,buff_size,min_notify);
  }
  else if(mod_type == gr::signal_exciter::PSK){
    d_mod = new Signal_PSK(mod_order,mod_offset,1,pulseshape,1,seed,enable,buff_size,min_notify);
  }
  else if(mod_type == gr::signal_exciter::PAM){
    d_mod = new Signal_PAM(mod_order,mod_offset,1,pulseshape,1,seed,enable,buff_size,min_notify);
  }
  else if(mod_type == gr::signal_exciter::ASK){
    d_mod = new Signal_PAM(mod_order,mod_offset,1,pulseshape,1,seed,enable,buff_size,min_notify);
  }
  else{
    throw_runtime("Unknown Base Modulation Type, choose ('PSK','QAM','PAM','ASK').\n");
  }

  d_fft_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*d_fftsize);
  d_fft_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*d_fftsize);
  memset( d_fft_in, 0, sizeof(fftwf_complex)*d_fftsize );
  memset( d_fft_out, 0, sizeof(fftwf_complex)*d_fftsize );
  d_fft = fftwf_plan_dft_1d(d_fftsize, d_fft_in, d_fft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

  generate_pilots();

  d_backoff = pow(10.,d_backoff/10.);
  d_sym_counter = 0;
  items_written = 0;

  if(tap_len==0){
    d_interp_taps = std::vector<float>(1,1.);
  }
  else{
    double power_check = 0.;
    d_interp_taps = std::vector<float>(tap_len);
    for(size_t idx = 0; idx < tap_len; idx++){
      d_interp_taps[idx] = interp_taps[idx];
      power_check += interp_taps[idx]*interp_taps[idx];
    }
    double normalizer = sqrt(double(interp)/power_check);
    for(size_t idx = 0; idx < tap_len; idx++){
      d_interp_taps[idx] *= normalizer;
    }
  }
  if(d_samp_overlap){
    d_taper = std::vector<float>(d_samp_overlap);
    for(size_t idx = 0; idx < d_samp_overlap; idx++){
      d_taper[idx] = symbol_taper[idx];
    }
  }
  else{
    d_taper = std::vector<float>(0);
  }

  load_firs();
  d_symbol_length = d_fftsize+d_cp_len;

  d_align = volk_get_alignment();

  if(d_enable){
    d_running = true;
  }
}

Signal_OFDM::~Signal_OFDM()
{
  boost::mutex::scoped_lock scoped_lock(s_mutex_fftw);
  d_running = false;
  fftwf_free(d_fft_in);
  fftwf_free(d_fft_out);
  fftwf_destroy_plan(d_fft);

  for(size_t idx = 0; idx < d_interp; idx++){
    delete d_firs[idx];
  }

  delete d_mod;
}


void
Signal_OFDM::generate_symbols(complexf* output, size_t symbol_count)
{
  std::vector<complexf> raw_sym(d_active, complexf(0.,0.));
  d_mod->generate_symbols( &raw_sym[0], d_active );
  for(size_t idx = 0; idx < d_active; idx++){
    output[d_active_list[idx]] = raw_sym[idx];
  }
  d_symb_gen_count++;
}

void
Signal_OFDM::generate_signal(complexf* output, size_t sample_count)
{
  size_t oo(0);
  d_samp_gen_count++;

  if(d_first_pass){
    if(d_sync){
      //if using a 'sync' then increment the symbols per frame to accomodate
      d_spf += 2;
    }
    //how long is a frame
    d_frame_length = d_spf*(d_symbol_length);
    float frame_counter = float(d_hist)/float(d_frame_length);
    d_frames_needed = ceil(frame_counter);

    //make a frame at a time
    d_frame = std::vector<complexf>(d_frame_length+d_samp_overlap, complexf(0.,0.));
    d_past = std::vector<complexf>(d_hist);
    
    //space for the individual ofdm symbols
    d_prefft_cache = std::vector<complexf>(d_fftsize, complexf(0.,0.));
    d_postfft_cache = std::vector<complexf>(d_fftsize, complexf(0.,0.));

    //pseudo sequences for the 'sync' symbols
    d_pn1 = std::vector<size_t>(d_fftsize/2, 0);
    d_pn2 = std::vector<size_t>(d_fftsize/2, 0);
    d_pn2[0] = 3;
    for(size_t idx = 1; idx < d_fftsize/2; idx++){
      d_pn1[idx] = (d_pn1[idx-1]+1)%4;
      d_pn2[idx] = (d_pn2[idx-1]+3)%4;
    }

    size_t left_to_fill = d_hist;
    size_t mincpy;

    //fill the frame_cache
    for(size_t fidx = 0; fidx < d_frames_needed; fidx++){
      //generate frame
      generate_frame();

      //determine how much to copy of the frame
      mincpy = std::min(d_frame_length+d_samp_overlap, left_to_fill);
      if(mincpy > 0){
        if(mincpy > d_frame_length){//need to blend the overlap of the start of the frame
          for(size_t idx = (d_frame_length+d_samp_overlap)-mincpy; idx < d_samp_overlap; idx++){
            d_past[left_to_fill-mincpy] = d_past[left_to_fill-mincpy]*d_taper[d_samp_overlap-1-idx]
                                        + d_frame[d_frame_length-mincpy]*d_taper[idx];
            mincpy--;//remove the overlap from the copy count
          }
        }
        memcpy( &d_past[left_to_fill-mincpy], &d_frame[d_frame_length+d_samp_overlap-mincpy], mincpy*sizeof(complexf) );
        left_to_fill -= mincpy;
      }
    }

    //print_buffer_complexf("Prehistory", &d_past[0], d_past.size());

    d_first_pass = false;
  }

  filter(sample_count, output);

}



void
Signal_OFDM::sort_subcarriers()
{
  size_t half = d_active/2;
  size_t offset = d_fftsize/2+(d_fftsize/2-half);
  offset = (d_active%2) ? offset-1 : offset;
  for(size_t idx = 0; idx < half+(d_active%2); idx++){
    d_active_list[idx] = offset++;
  }
  offset = (d_fftsize!=d_active);
  for(size_t idx = d_active-half; idx < d_active; idx++){
    d_active_list[idx] = offset++;
  }
}


void
Signal_OFDM::generate_sync()
{
  std::vector<complexf> first_sym(d_fftsize/2, complexf(0.,0.));
  d_mod->generate_symbols(&first_sym[0], d_fftsize/2);
  memset( &d_prefft_cache[0], 0, sizeof(complexf)*d_fftsize );
  for(size_t idx = 0; idx < d_fftsize/2; idx++){
    d_prefft_cache[idx*2] = first_sym[idx]*complexf(sqrt(2.),0.);
  }
  do_fft();
  memcpy( &d_frame[0], &d_postfft_cache[d_fftsize-(d_cp_len+d_samp_overlap)], sizeof(complexf)*(d_cp_len+d_samp_overlap) );
  memcpy( &d_frame[d_cp_len+d_samp_overlap], &d_postfft_cache[0], sizeof(complexf)*d_fftsize );
  for(size_t idx = 0; idx < d_fftsize/2; idx++){
    d_prefft_cache[idx*2] = first_sym[idx] * exp(complexf(0., 2*M_PI*d_pn1[idx]/4.));
    d_prefft_cache[idx*2+1] = d_prefft_cache[idx*2] * exp(complexf(0., 2*M_PI*d_pn2[idx]/4.));
  }
  do_fft();
  for(size_t idx = 0; idx < d_samp_overlap; idx++){
    d_frame[d_symbol_length+idx] = d_frame[d_symbol_length+idx]*d_taper[d_samp_overlap-1-idx] + d_postfft_cache[d_fftsize-(d_cp_len+d_samp_overlap)+idx]*d_taper[idx];
  }
  memcpy( &d_frame[d_symbol_length+d_samp_overlap], &d_postfft_cache[d_fftsize-d_cp_len], sizeof(complexf)*d_cp_len );
  memcpy( &d_frame[d_symbol_length+d_cp_len+d_samp_overlap], &d_postfft_cache[0], sizeof(complexf)*d_fftsize );
}


void
Signal_OFDM::generate_pilots()
{
  d_pilots = std::vector<complexf>(d_npilots,complexf(0.,0.));
  d_mod->generate_symbols(&d_pilots[0], d_npilots);
}


void
Signal_OFDM::do_fft()
{
  memset( &d_fft_in[0], 0, sizeof(fftwf_complex)*d_fftsize );
  memcpy( &d_fft_in[0], &d_prefft_cache[0], sizeof(complexf)*d_fftsize );
  //print_buffer_complex( "Input", &d_fft_in[0], d_fftsize );
  memset( &d_fft_out[0], 0, sizeof(fftwf_complex)*d_fftsize );
  fftwf_execute( d_fft );
  //print_buffer_complex( "Output", &d_fft_out[0], d_fftsize );
  float scale = 1/sqrt(float(d_fftsize));
  for(size_t idx = 0; idx < d_fftsize; idx++){
    d_postfft_cache[idx] = complexf(d_fft_out[idx][0]*scale,d_fft_out[idx][1]*scale);
  }
}

void
Signal_OFDM::print_buffer_complexf(std::string str, complexf* buff, size_t length)
{
  printf("%s Buffer (%lu):\n[(%1.3f,%1.3f)",str.c_str(),length,buff[0].real(),buff[0].imag());
  for(size_t idx=1; idx<length; idx++){
    printf(", (%1.3f,%1.3f)",buff[idx].real(),buff[idx].imag());
  }
  printf("]\n");
}
void
Signal_OFDM::print_buffer_complex(std::string str, fftwf_complex* buff, size_t length)
{
  printf("%s Buffer (%lu):\n[(%1.3f,%1.3f)",str.c_str(),length,buff[0][0],buff[0][1]);
  for(size_t idx=1; idx<length; idx++){
    printf(", (%1.3f,%1.3f)",buff[idx][0],buff[idx][1]);
  }
  printf("]\n");
}





void 
Signal_OFDM::auto_fill_symbols()
{}

void 
Signal_OFDM::auto_fill_signal()
{}



void
Signal_OFDM::filter( size_t nout, complexf* out )
{
  size_t total_samps = d_past.size()*d_interp;
  size_t used_samps = d_branch_offset;//how many samples are used from total_samps
  size_t part_samps = (d_interp-(d_branch_offset%d_interp))%d_interp;//frac samps left of first symbol (part_samps/d_interp) < d_interp;
  
  //subtract out the samples that went last time
  total_samps = total_samps - used_samps;

  //subtract out the symbols needed for the filter
  total_samps = total_samps - (d_hist*d_interp);

  if(total_samps > d_past.size()*d_interp){
    printf("This is a logic error.\n");
  }

  //printf("There are %lu samples supposedly left in the past.\n",total_samps);

  size_t samps_needed;
  if(nout < total_samps){
    //printf("This is highly unlikely.. hopefully.\n");
    samps_needed = 0;
  }
  else{
    samps_needed = nout - total_samps;
  }
  float fractional_N = float(samps_needed);
  float fractional_D = float(d_interp);
  float fractional = fractional_N/fractional_D;
  size_t symbs_needed = ceil(fractional);
  fractional_N = float(symbs_needed);
  fractional_D = float(d_frame_length);
  fractional = fractional_N/fractional_D;
  size_t frames_needed = ceil(fractional);
  //printf("Need %lu symbs, %lu samps, from %lu frames.\n", symbs_needed, samps_needed, frames_needed);

  symbs_needed = frames_needed * d_frame_length;

  size_t total_input_len = symbs_needed + d_past.size() + d_samp_overlap;

  size_t min_backtrack;
  d_filt_in = (gr_complex*) volk_malloc( total_input_len*sizeof(complexf), d_align );
  memcpy( &d_filt_in[0], &d_past[0], d_past.size()*sizeof(complexf) );
  for(size_t ff = 0; ff < frames_needed; ff++){
    generate_frame();

    min_backtrack = std::min(ff*d_frame_length + d_past.size(), d_samp_overlap);
    if(min_backtrack){
      for(size_t idx = (d_samp_overlap)-min_backtrack; idx < d_samp_overlap; idx++){
        d_filt_in[ff*d_frame_length + d_past.size() - d_samp_overlap + idx] = d_filt_in[ff*d_frame_length + d_past.size() - d_samp_overlap + idx] * d_taper[d_samp_overlap-1-idx]
                                                                            + d_frame[idx]*d_taper[idx];
      }
      memcpy( &d_filt_in[ff*d_frame_length + d_past.size()], &d_frame[d_samp_overlap], d_frame_length*sizeof(complexf) );
    }
    else{//there might still be overlap needed
      if(d_samp_overlap){//still need overlap, whelp doing it...
        for(size_t idx = 0; idx < d_samp_overlap; idx++){
          d_filt_in[idx] = d_filt_in[idx]*d_taper[d_samp_overlap-1-idx] + d_frame[idx]*d_taper[idx];
        }
        memcpy( &d_filt_in[ff*d_frame_length + d_past.size() + d_samp_overlap], &d_frame[d_samp_overlap], d_frame_length*sizeof(complexf) );
      }
      else{//no need for overlap
        memcpy( &d_filt_in[ff*d_frame_length + d_past.size()], &d_frame[0], d_frame_length*sizeof(complexf) );
      }
    }
  }

  size_t oo = 0, ii = 0;

  while(oo<nout){
    out[oo] = d_firs[d_branch_offset]->filter( &d_filt_in[ii] );
    d_branch_offset = (d_branch_offset+1)%d_interp;
    if(d_branch_offset == 0){
      ii++;
    }
    oo++;
  }
  size_t remaining = total_input_len - ii;
  if(remaining < d_hist){
    printf("OFDM - There isn't enough left in the buffer!!!\n");
  }
  d_past = std::vector<complexf>( &d_filt_in[ii], &d_filt_in[total_input_len] );
  
  volk_free(d_filt_in);

}

void
Signal_OFDM::load_firs()
{
  //printf("OFDM::Attempting to load filters.\n");
  size_t intp = d_interp;
  size_t ts = d_interp_taps.size() / intp;
  if((d_interp_taps.size() % intp)){
    ts++;
  }
  //printf("OFDM:: Intp(%lu), ts(%lu), sps(%d), d_overlap(%lu).\n",intp,ts,d_sps,d_overlap);

  //std::vector< std::vector<float> > xtaps(intp);
  d_taps = std::vector< std::vector<float> >(intp);

  //printf("OFDM:: taps_size(%lu).\n",d_taps.size());

  for(size_t idx = 0; idx < intp; idx++){
    d_taps[idx].resize(ts);
    memset( &d_taps[idx][0], 0., ts*sizeof(float) );
  }
  //printf("OFDM:: taps init 0.\n");

  for(size_t idx = 0; idx < d_interp_taps.size(); idx++){
    d_taps[idx % intp][idx / intp] = d_interp_taps[idx];
  }
  //printf("OFDM:: taps filled.\n");

  d_firs = std::vector< gr::filter::kernel::fir_filter_ccf *>(intp);
  std::vector<float> dummy_taps(1,0.);
  for(size_t idx = 0; idx < intp; idx++){
    d_firs[idx] = new gr::filter::kernel::fir_filter_ccf(1,dummy_taps);
  }
  //printf("OFDM:: filters made.\n");
  for(size_t idx = 0; idx < intp; idx++){
    d_firs[idx]->set_taps(d_taps[idx]);
  }
  //printf("OFDM:: taps loaded.\n");
  d_hist = ts-1;
}


void
Signal_OFDM::generate_frame()
{
  d_symbol_idx = 0;
  //printf("Generating Frame\n");
  if(d_sync){
    //printf("Inserting sync symbols (0,1)\n");
    // if sync is to be added, is it time to add?
    generate_sync();
    d_symbol_idx += 2;
  }

  size_t pidx = 0;
  while(d_symbol_idx < d_spf){
    memset( &d_prefft_cache[0], 0, sizeof(complexf)*d_fftsize );
    memset( &d_postfft_cache[0], 0, sizeof(complexf)*d_fftsize );
    generate_symbols( &d_prefft_cache[0], 1 );

    //add pilots here
    if(d_ppf){
      //pilots are given per frame
      //if 'sync' is used, start on second symbol
      //printf("Adding frame based pilots on symbol (%lu)\n",d_symbol_idx);
      size_t poffset = d_sync ? 2 : 0;
      while((d_pilot_list[pidx]-d_active*(d_symbol_idx-poffset))<d_active){
        //while within a single frame
        //load the framed pilots
        d_prefft_cache[d_active_list[d_pilot_list[pidx]-d_active*(d_symbol_idx-poffset)]] = d_pilots[pidx];
        pidx++;
        if(pidx >= d_npilots){
          //all pilots have been loaded
          break;
        }
      }
    }
    else{
      //pilots are given per symbol
      //printf("Adding symbol based pilots on symbol (%lu)\n",d_symbol_idx);
      while(pidx < d_npilots){
        //load pilots into symbol
        d_prefft_cache[d_active_list[d_pilot_list[pidx]]] = d_pilots[pidx];
        pidx++;
      }
      //pilots are reset every symbol
      pidx = 0;
    }
    do_fft();

    for(size_t idx = 0; idx < d_samp_overlap; idx++){
      d_frame[(d_symbol_length)*(d_symbol_idx)+idx] = d_frame[(d_symbol_length)*(d_symbol_idx)+idx]*d_taper[d_samp_overlap-1-idx]
                                                    + d_postfft_cache[d_fftsize-(d_cp_len+d_samp_overlap)+idx]*d_taper[idx];
    }
    memcpy( &d_frame[(d_symbol_length)*(d_symbol_idx)+d_samp_overlap],
            &d_postfft_cache[d_fftsize-d_cp_len], sizeof(complexf)*d_cp_len );
    memcpy( &d_frame[(d_symbol_length)*(d_symbol_idx)+d_cp_len+d_samp_overlap],
            &d_postfft_cache[0], sizeof(complexf)*d_fftsize );
    d_symbol_idx++;
    d_sym_counter++;
  }
  d_symbol_idx = 0;
  complexf backoff(1/d_backoff,0.);
  for(size_t idx = 0; idx < d_frame.size(); idx++){
    d_frame[idx] = d_frame[idx]*backoff;
  }
}





void
Signal_OFDM::throw_runtime(std::string err)
{
  d_running = false;
  fftwf_free(d_fft_in);
  fftwf_free(d_fft_out);
  fftwf_destroy_plan(d_fft);

  for(size_t idx = 0; idx < d_interp; idx++){
    delete d_firs[idx];
  }

  delete d_mod;

  std::string err_msg = "Signal_OFDM::";
  err_msg += err;
  throw std::runtime_error(err_msg.c_str());
}







