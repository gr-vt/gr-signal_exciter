

#include "signal_usb.hpp"
#include <stdio.h>//////////////////////////////////

Signal_USB::Signal_USB(float mod_idx, size_t components, float* mu,
                      float* sigma, float* weight, float samp_rate,
                      size_t tap_count, int seed, bool norm, float fso,
                      bool enable, size_t buff_size, size_t min_notify)
  : d_mod_idx(mod_idx),
    d_tap_count(tap_count),
    d_enable(enable),
    d_buffer_size(buff_size),
    d_notify_size(min_notify),
    d_agc(),
    d_norm(norm)
{
  set_seed(seed);
  boost::mutex::scoped_lock scoped_lock(fftw_lock());
  d_gmm_tap_gen.set_params(components, mu, sigma, weight, samp_rate, tap_count);
  generate_taps();

  d_rng = new gr::random(d_seed, 0, 1);

  d_burn = 20;
  d_agc.set_rate(5e-4);
  if(d_enable){
    d_running = true;
    auto_fill_symbols();
    auto_fill_signal();
  }
  if(d_norm){
    std::vector<complexf> burn_buff(50);
    for(size_t count = 0; count < d_burn; count++){
      generate_signal( &burn_buff[0], d_burn );
    }
  }
  d_first_pass = true;


  printf("usb: fso: %0.3e\n",fso);
  d_fso = fso;

  d_align = volk_get_alignment();
  // Generate and load the GNURadio FIR Filters with the pulse shape.
  load_fir();
  //printf("DSB::Loaded FIR filters.\n");
}

Signal_USB::~Signal_USB()
{
  delete d_fir;
  delete d_frac_filt;
  if(d_enable){
    d_running = false;
    d_TGroup.join_all();
    delete d_Sy;
  }
}

void
Signal_USB::generate_symbols(complexf* output, size_t symbol_count)
{
  if(d_enable){
    size_t filled(0);
    while((filled < symbol_count) && d_running){
      filled += d_Sy->bmemcpy( &output[filled], symbol_count-filled, false );
    }
  }
  else{
    //printf("Generating message.\n");
    std::vector<float> message(symbol_count,0.);
    float scale = 1./3.;
    for(size_t idx = 0; idx < symbol_count; idx++){
      message[idx] = scale * d_rng->gasdev();
      output[idx] = complexf(d_mod_idx*message[idx],0.);
    }
    //printf("Generated message.\n");
    //printf("[%1.3e",message[0]);
    //for(size_t idx = 1; idx<symbol_count; idx++){
    //  printf(", %1.3e",message[idx]);
    //}
    //printf("]\n");
    //printf("[%1.3e",output[0].real());
    //for(size_t idx = 1; idx<symbol_count; idx++){
    //  printf(", %1.3e",output[idx].real());
    //}
    //printf("]\n");
  }
}

void
Signal_USB::generate_signal(complexf* output, size_t sample_count)
{
  if(d_first_pass){
    d_past = std::vector<float>(d_hist);
    d_symbol_cache = std::vector<complexf>(d_hist);
    generate_symbols( &d_symbol_cache[0], d_hist );
    for(size_t idx = 0; idx < d_hist; idx++){
      d_past[idx] = d_symbol_cache[idx].real();
    }
    //fso needs
    filter( d_frac_cache.size(), &d_frac_cache[0] );

    d_first_pass = false;
  }

  d_time_shift_in = (complexf *) volk_malloc(
                      (sample_count+d_frac_cache.size())*sizeof(complexf),
                      d_align);
  memcpy( &d_time_shift_in[0], &d_frac_cache[0],
          d_frac_cache.size()*sizeof(complexf) );
  filter( sample_count, &d_time_shift_in[d_frac_cache.size()] );
  d_frac_filt->filterN( &output[0], &d_time_shift_in[0], sample_count );
  memcpy( &d_frac_cache[0], &d_time_shift_in[sample_count],
          d_frac_cache.size()*sizeof(complexf) );
  volk_free(d_time_shift_in);

  if(d_norm){
    d_agc.scaleN( output, output, sample_count );
  }
}

void
Signal_USB::filter( size_t nout, complexf* out )
{
  size_t total_samps = d_past.size();

  size_t symbs_needed = nout;

  size_t total_input_len = symbs_needed + d_past.size();
  d_filt_in = (float*) volk_malloc( total_input_len*sizeof(float), d_align );
  memcpy( &d_filt_in[0], &d_past[0], d_past.size()*sizeof(float) );
  d_symbol_cache = std::vector<complexf>(symbs_needed);
  generate_symbols( &d_symbol_cache[0], symbs_needed );
  for(size_t idx = 0; idx < symbs_needed; idx++){
    d_filt_in[total_samps+idx] = d_symbol_cache[idx].real();
  }

  size_t oo = 0, ii = 0;
  d_fm = std::vector<float>(nout,0.);
  d_output_cache = std::vector<complexf>(nout,complexf(0.,0.));

  while( oo < nout ){
    d_fm[oo++] = d_fir->filter( &d_filt_in[ii++] );
  }

  size_t remaining = total_input_len - ii;
  if(remaining < d_hist){
    fprintf(stderr,"USB: nout(%lu), til(%lu), ii(%lu), past(%lu)\n",nout,total_input_len,ii,d_past.size());
    fprintf(stderr,"USB - There isn't enough left in the buffer!!! (%lu,%lu)\n",remaining,d_hist);
  }
  d_past = std::vector<float>( &d_filt_in[ii], &d_filt_in[total_input_len] );

  boost::mutex::scoped_lock scoped_lock(fftw_lock());
  d_fft_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nout);
  d_fft_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nout);
  d_fft = fftwf_plan_dft_r2c_1d( nout, &d_fm[0], d_fft_in, FFTW_ESTIMATE );
  d_ifft = fftwf_plan_dft_1d( nout, d_fft_in, d_fft_out, FFTW_BACKWARD, FFTW_ESTIMATE );

  do_fft(nout);
  for(size_t idx = 1; idx < nout/2; idx++){
    d_fft_in[idx][0] *= 2.;
    d_fft_in[idx][1] *= 2.;
  }
  for(size_t idx = nout/2+1; idx <= nout; idx++){
    d_fft_in[idx%nout][0] = 0.;
    d_fft_in[idx%nout][1] = 0.;
  }







  do_ifft(nout);//output stored in d_output_cache

  fftwf_free( d_fft_in );
  fftwf_free( d_fft_out );
  fftwf_destroy_plan( d_fft );
  fftwf_destroy_plan( d_ifft );

  memcpy( &out[0], &d_output_cache[0], sizeof(complexf)*nout);

  volk_free(d_filt_in);

}


void
Signal_USB::generate_taps()
{
  d_gmm_tap_gen.get_taps(d_taps);
  d_window = std::vector<float>(d_tap_count,0.);
  for(size_t idx = 0; idx < d_tap_count; idx++)
  {
    if(!((idx==0)||(idx==d_tap_count-1)))
    {
      d_window[idx] = 0.42 - 0.5*cos((2*M_PI*float(idx))/float(d_tap_count-1)) + 0.08*cos((4*M_PI*float(idx))/float(d_tap_count-1));
    }
    d_taps[idx] = d_taps[idx]*d_window[idx];
  }
  d_hist = d_tap_count-1;

  d_proto_taps = gr::filter::firdes::low_pass_2(1,1,0.5,0.05,61,
                        gr::filter::firdes::WIN_BLACKMAN_hARRIS);
}

void
Signal_USB::load_fir()
{
  std::vector<float> dummy_taps;

  d_fir = new gr::filter::kernel::fir_filter_fff(1, dummy_taps);

  d_fir->set_taps(d_taps);

  d_frac_filt = new gr::filter::kernel::fir_filter_ccf(1, dummy_taps);
  std::vector<float> shifted_taps;
  time_offset(shifted_taps, d_proto_taps, d_fso);
  d_frac_filt->set_taps(shifted_taps);
  d_frac_cache = std::vector<complexf>(shifted_taps.size()-1);
}


void
Signal_USB::do_fft(size_t samp_count)
{
  memset( &d_fft_in[0], 0, sizeof(fftwf_complex)*samp_count );
  fftwf_execute( d_fft );
}

void
Signal_USB::do_ifft(size_t samp_count)
{
  memset( &d_fft_out[0], 0, sizeof(fftwf_complex)*samp_count );
  fftwf_execute( d_ifft );
  complexf scale = complexf(1/float(samp_count),0.);
  for(size_t idx = 0; idx < samp_count; idx++){
    d_output_cache[idx] = complexf(d_fft_out[idx][0],d_fft_out[idx][1])*scale;
  }
}

void
Signal_USB::auto_fill_symbols()
{
  d_Sy = new signal_threaded_buffer<complexf>(d_buffer_size,d_notify_size);

  d_TGroup.create_thread( boost::bind(&Signal_USB::auto_gen_GM, this) );

}

void
Signal_USB::auto_fill_signal()
{}



void
Signal_USB::auto_gen_GM()
{
  size_t buff_size(0), buff_pnt(0);
  float scale = 1./3.;
  std::vector<float> buffer(d_buffer_size,0.);
  std::vector<complexf> message(d_buffer_size,complexf(0.,0.));
  while(d_running){
    for(size_t idx = 0; idx < buffer.size(); idx++){
      buffer[idx] = scale * d_rng->gasdev();
      message[idx] = complexf(d_mod_idx*buffer[idx],0.);
    }
    while((buff_pnt < buffer.size()) && d_running){
      buff_pnt += d_Sy->bmemcpy( &message[buff_pnt], message.size()-buff_pnt, true );
    }
    buff_pnt = 0;
  }
}
