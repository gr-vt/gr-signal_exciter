

#include "signal_qam.hpp"
#include <stdio.h>//////////////////////////////////
#include <stdexcept>
#include <algorithm>

Signal_QAM::Signal_QAM(int order, float offset, int sps, float* pulse_shape, size_t length, int seed,
                        bool enable_fso, float fso, bool enable, size_t buff_size, size_t min_notify)
  : d_order(order),
    d_offset(offset),
    d_sps(sps),
    d_samp_offset(0),
    d_protected(0),
    d_enable(enable),
    d_buffer_size(buff_size),
    d_notify_size(min_notify),
    d_sample_count(0),
    d_symbol_count(0)
{
  //printf("Init.\n");
  get_indicator();
  set_seed(seed);
  //printf("Seeded.\n");

  double power_check = 0.;
  d_pulse_shape = std::vector<float>(length, 0.);
  for(size_t idx = 0; idx < length; idx++){
    d_pulse_shape[idx] = pulse_shape[idx];
    power_check += pulse_shape[idx]*pulse_shape[idx];
  }

  double normalizer = sqrt(sps/power_check);
  for(size_t idx = 0; idx<length; idx++){
    d_pulse_shape[idx] *= normalizer;
  }

  //printf("Shaped.\n");
  float check = float(length)/float(sps);
  d_overlap = size_t(int(floor(check)));
  if(floor(check) < check)
  {
    d_overlap++;
  }
  //printf("Overlaped.\n");

  d_first_pass = true;

  int new_order=0;
  int k=0;
  while(order>0){
    order = order>>1;
    if(!new_order) new_order = 1;
    else new_order = new_order*2;
    k++;
  }
  //printf("New Order.\n");

  if(new_order != d_order){
    if(new_order > 0) d_order = new_order;
    else d_order = 2;
  }
  //printf("Ok.\n");

  create_symbol_list(k-1);

  d_rng = new gr::random(d_seed, 0, d_order);


  if(d_sps > d_pulse_shape.size()){
    printf("The pulse_shape is shorter than sps, this will crash.\n");
  }


  enable_fractional_offsets(enable_fso, fso);

  d_align = volk_get_alignment();

  load_firs();
  //printf("QAM::Loaded FIR filters.\n");


  if(d_enable){
    d_running = true;
    auto_fill_symbols();
    auto_fill_signal();
  }
  //printf("QAM::Fully made\n");
}

Signal_QAM::~Signal_QAM()
{
  if(d_enable){
    d_running = false;/*
    std::vector<complexf> null_buff(d_buffer_size,complexf(0.,0.));
    d_Sy->bmemcpy( &null_buff[0], d_buffer_size, true );
    d_Sy->bmemcpy( &null_buff[0], d_buffer_size, false );
    d_Sy->cleanup();*/
    d_TGroup.join_all();
    delete d_Sy;
  }
  for(size_t idx = 0; idx < d_sps; idx++){
    delete d_firs[idx];
  }
  delete d_rng;

}

void
Signal_QAM::generate_symbols(complexf* output, size_t symbol_count)
{
  if(d_enable){
    size_t filled(0);
    while((filled < symbol_count) && d_running){
      filled += d_Sy->bmemcpy( &output[filled], symbol_count-filled, false );
    }
    // All requested symbols have been made

    size_t zeros(0),infs(0),nans(0);
    for(size_t idx = 0; idx<symbol_count; idx++){
      if(std::isnan(output[idx].real())||std::isnan(output[idx].imag()))                        nans++;
      if(std::isinf(output[idx].real())||std::isinf(output[idx].imag()))                        infs++;
      if((output[idx].real()*output[idx].real() + output[idx].imag()*output[idx].imag()) == 0.) zeros++;
    }
    if(zeros | infs | nans){
      printf("Fault found in generate_symbols call #%lu.\n",d_symb_gen_count);
      printf("Consumer requested %lu symbols, generated %lu samples.\n",symbol_count,filled);
      printf("Within generated samples found: %lu zeros, %lu NANs, %lu Infs.\n",zeros,infs,nans);
    }
  }
  else{
    for(size_t idx = 0; idx<symbol_count; idx++){
      //int data = rand()%d_order;
      int data = d_rng->ran_int();
      output[idx] = d_symbol_list[data];
    }
  }
}

void
Signal_QAM::generate_signal(complexf* output, size_t sample_count)
{
  size_t oo = 0;
  size_t interp = size_t(d_sps);
  d_samp_gen_count++;

  // If using threading
  if(d_first_pass){
    //printf("First pass: doc(%lu)\n",d_hist);
    d_past = std::vector<complexf>(d_hist);
    generate_symbols( &d_past[0], d_hist );
    d_first_pass = false;
  }

  filter(sample_count, output);

/*/////////////////DEBUGGING
  size_t zeros(0),infs(0),nans(0);
  std::vector<size_t> index_locations;
  for(size_t idx = 0; idx<sample_count; idx++){
    if(std::isnan(output[idx].real())||std::isnan(output[idx].imag()))                        nans++;
    if(std::isinf(output[idx].real())||std::isinf(output[idx].imag()))                        infs++;
    if((output[idx].real()*output[idx].real() + output[idx].imag()*output[idx].imag()) == 0.){
      zeros++;
      index_locations.push_back(idx);
    }
  }
  if(zeros | infs | nans){
    printf("Fault found in generate_signal call #%lu.\n",d_samp_gen_count);
    printf("Consumer requested %lu samples, generated %lu samples.\n",sample_count,oo);
    printf("Within generated samples found: %lu zeros, %lu NANs, %lu Infs.\n",zeros,infs,nans);
    for(size_t locs = 0; locs < index_locations.size(); locs++){
      printf("%lu, ",index_locations[locs]);
    }
    printf("\n");
  }
////////////////////////////*/
  d_sample_count += sample_count;
}



void
Signal_QAM::create_symbol_list(int k)
{
  //printf("Creating the symbol list\n");
  float spillage = float(k)/2.0;
  int r(floor(spillage)),i(floor(spillage));
  if(floor(spillage) < spillage){
    r++;
  }
  int R(1),I(1);
  for(int idx = 0; idx < i; idx++){
    R = R*2;
    I = I*2;
  }
  if(r>i) R = R*2;

  std::vector<float> reals(R,0.);
  std::vector<float> imags(I,0.);

  for(int idx = 0; idx < R; idx++){
    reals[idx] = idx*2-R+1;
  }
  if(I > 1){
    for(int idx = 0; idx < I; idx++){
      imags[idx] = idx*2-I+1;
    }
  }

  d_symbol_list = std::vector<complexf>(d_order,complexf(0.,0.));
  for(int data = 0; data < d_order; data++){
    r = data / I;
    i = data % I;
    d_symbol_list[data] = complexf(reals[r],imags[i]);
  }

  float power_check = 0.;
  float fr(0.),fi(0.);
  for(int idx = 0; idx < d_order; idx++){
    fr = d_symbol_list[idx].real();
    fi = d_symbol_list[idx].imag();
    power_check += (fr*fr + fi*fi)/float(d_order);
  }
//  printf("The average power seen is %1.3f\n",power_check);///////////////////

  complexf scale_shift = complexf(sqrt(1/power_check),0.)*exp(complexf(0.,d_offset));

//  printf("The scale_shift is (%1.3f,%1.3f)\n",scale_shift.real(),scale_shift.imag());//////////////////

  for(int idx = 0; idx < d_order; idx++){
    d_symbol_list[idx] = d_symbol_list[idx]*scale_shift;
  }
/*/////////////////////////////////////////////////////////////////////
  power_check = 0.;
  for(int idx = 0; idx < d_order; idx++){
    fr = d_symbol_list[idx].real();
    fi = d_symbol_list[idx].imag();
    power_check += (fr*fr + fi*fi)/float(d_order);
  }
  printf("The average power seen is %1.3f\n",power_check);
/////////////////////////////////////////////////////////////////////*/

/*/////////////////////////////////////////////////////////////////////
  printf("\033[38;5;21mThe symbols are:[(%1.3f,%1.3f)",d_symbol_list[0].real(),d_symbol_list[0].imag());
  for(int idx = 1; idx < d_order; idx++){
    printf(", (%1.3f,%1.3f)",d_symbol_list[idx].real(),d_symbol_list[idx].imag());
  }
  printf("]\n\033[0m");
/////////////////////////////////////////////////////////////////////*/
  //printf("The bit count is %d\n",k);
  //printf("Created the symbol list\n");

}



void
Signal_QAM::filter( size_t nout, complexf* out )
{
/*  size_t oo(0),ii(0),cached(0),interp(d_sps),sidx(0);
  float fractional(0);
  size_t frac_fl(0),frac_cl(0);
  if(d_samp_offset) cached = interp - d_samp_offset;

  fractional = (float(nout)-float(cached))/float(interp);
  if(nout > cached){
    frac_fl = floor(fractional);
    frac_cl = ceil(fractional);
    d_output_cache = std::vector<complexf>(frac_cl);
    generate_symbols( &d_output_cache[0], frac_cl );
    ii += frac_cl;
  }
  else{
    d_output_cache = std::vector<complexf>(0);
  }

  d_filt_in = (gr_complex*) volk_malloc((d_symbol_cache.size()+frac_cl)*sizeof(complexf),d_align);
  memcpy( &d_filt_in[0], &d_symbol_cache[0], d_symbol_cache.size()*sizeof(complexf) );
  memcpy( &d_filt_in[d_symbol_cache.size()], &d_output_cache[0], d_output_cache.size()*sizeof(complexf) );

  if(cached > nout){
    for(size_t branch = d_samp_offset; (branch<interp)&&(oo<nout); branch++){
      out[oo] = d_firs[branch]->filter( &d_filt_in[0] );
      oo++;
      d_samp_offset = (branch+1)%interp;
    }

    if(!d_samp_offset){
      std::vector<complexf> temp(d_symbol_cache.begin()+1,d_symbol_cache.end());
      d_symbol_cache = std::vector<complexf>(temp.begin(),temp.end());
      sidx++;
    }

  }
  else{
    if(cached > 0){
      for(size_t branch = d_samp_offset; (branch<interp)&&(oo<nout); branch++){
        out[oo] = d_firs[branch]->filter( &d_filt_in[0] );
        oo++;
        d_samp_offset = (branch+1)%interp;
      }

      if(!d_samp_offset){
        std::vector<complexf> temp(d_symbol_cache.begin()+1,d_symbol_cache.end());
        d_symbol_cache = std::vector<complexf>(temp.begin(),temp.end());
        sidx++;
      }

    }

    size_t remaining(0);
    remaining = (nout>oo) ? nout - oo : 0;
    if(remaining){
      cached = oo;
      while(oo-cached < frac_fl*interp){
        for(size_t branch = 0; branch < interp; branch++){
          out[oo] = d_firs[branch]->filter( &d_filt_in[sidx] );
          oo++;
        }

        sidx++;
      }

      remaining = nout-oo;
    }
    //Fractional space remaining
    if(remaining){
      for(size_t branch = 0; branch < remaining; branch++){
        out[oo] = d_firs[branch]->filter( &d_filt_in[sidx] );
        oo++;
      }

      d_samp_offset = remaining;
    }

    if(d_samp_offset){
      d_symbol_cache.resize(d_hist+1);
      if(frac_cl > d_hist+1){
        memcpy( &d_symbol_cache[0], &d_output_cache[frac_cl-(d_hist+1)], (d_hist+1)*sizeof(complexf) );
      }
      else{
        std::vector<complexf> temp(d_symbol_cache.begin()+frac_cl-1,d_symbol_cache.end());
        memcpy( &d_symbol_cache[0], &temp[0], temp.size()*sizeof(complexf) );
        memcpy( &d_symbol_cache[(d_hist+1)-frac_cl], &d_output_cache[0], frac_cl*sizeof(complexf) );
      }

    }
    else{
      if(frac_cl >= d_hist){
        memcpy( &d_symbol_cache[0], &d_output_cache[frac_cl-d_hist], d_hist*sizeof(complexf) );
      }
      else{
        std::vector<complexf> temp(d_symbol_cache.begin()+frac_cl,d_symbol_cache.end());
        memcpy( &d_symbol_cache[0], &temp[0], temp.size()*sizeof(complexf) );
        memcpy( &d_symbol_cache[d_hist-frac_cl], &d_output_cache[0], frac_cl*sizeof(complexf) );
      }

    }

  }

  volk_free(d_filt_in);


/////////////////////DEBUGGING
  size_t zeros(0),infs(0),nans(0);
  std::vector<size_t> index_locations;
  for(size_t idx = 0; idx<nout; idx++){
    if(std::isnan(out[idx].real())||std::isnan(out[idx].imag()))                        nans++;
    if(std::isinf(out[idx].real())||std::isinf(out[idx].imag()))                        infs++;
    if((out[idx].real()*out[idx].real() + out[idx].imag()*out[idx].imag()) == 0.){
      zeros++;
      index_locations.push_back(idx);
    }
  }
  if(zeros | infs | nans){
    printf("Fault found in generate_signal call #%lu.\n",d_samp_gen_count);
    printf("Consumer requested %lu samples, generated %lu samples.\n",nout,oo);
    printf("Within generated samples found: %lu zeros, %lu NANs, %lu Infs.\n",zeros,nans,infs);
    //for(size_t locs = 0; locs < index_locations.size(); locs++){
    //  printf("%lu, ",index_locations[locs]);
    //}
    //printf("\n");
    //if(debug_need)
    //  printf("symbols_in = [%1.1f, %1.1f",debug_val,d_symbol_cache[0].real());
    //else
      printf("symbols_in = [%1.1f",d_symbol_cache[0].real());
    for(size_t idx = 1; idx < d_symbol_cache.size(); idx++){
      printf(", %1.1f",d_symbol_cache[idx].real());
    }
    for(size_t idx = 0; idx < d_output_cache.size(); idx++){
      printf(", %1.1f",d_output_cache[idx].real());
    }
    printf("];\n");
  }
/////////////////////////////*/
  size_t total_samps = d_past.size()*d_sps;
  size_t used_samps = d_samp_offset;
  size_t part_samps = (d_sps-(d_samp_offset%d_sps))%d_sps;

  total_samps = total_samps - used_samps;

  total_samps = total_samps - (d_hist*d_sps);

  if(total_samps > d_past.size()*d_sps){
    printf("This is a logic error.\n");
  }

  size_t samps_needed;
  if(nout < total_samps){
    samps_needed = 0;
  }
  else{
    samps_needed = nout - total_samps;
  }
  float fractional_N = float(samps_needed);
  float fractional_D = float(d_sps);
  float fractional = fractional_N/fractional_D;
  size_t symbs_needed = ceil(fractional);

  size_t total_input_len = symbs_needed + d_past.size();

  d_filt_in = (gr_complex*) volk_malloc( total_input_len*sizeof(complexf), d_align );
  memcpy( &d_filt_in[0], &d_past[0], d_past.size()*sizeof(complexf) );
  generate_symbols( &d_filt_in[d_past.size()], symbs_needed );

  size_t oo = 0, ii = 0;

  while( oo < nout ){
    out[oo] = d_firs[d_samp_offset]->filter( &d_filt_in[ii] );
    d_samp_offset = (d_samp_offset+1)%d_sps;
    if(d_samp_offset == 0){
      ii++;
    }
    oo++;
  }
  size_t remaining = total_input_len - ii;
  if(remaining < d_hist){
    fprintf(stderr,"PSK - There isn't enough left in the buffer!!!\n");
  }
  d_past = std::vector<complexf>( &d_filt_in[ii], &d_filt_in[total_input_len] );

  volk_free(d_filt_in);

}




void
Signal_QAM::auto_fill_symbols()
{
  d_Sy = new signal_threaded_buffer<complexf>(d_buffer_size,d_notify_size);
  d_TGroup.create_thread( boost::bind(&Signal_QAM::auto_gen_SYMS, this) );
}

void
Signal_QAM::auto_fill_signal()
{}


void
Signal_QAM::load_firs()
{
  //printf("QAM::Attempting to load filters.\n");
  size_t intp = d_sps;
  //printf("QAM:: Intp(%lu), ts(%lu), sps(%d), d_overlap(%lu).\n",intp,ts,d_sps,d_overlap);

  d_firs = std::vector< gr::filter::kernel::fir_filter_ccf *>(intp);
  std::vector<float> dummy_taps;
  for(size_t idx = 0; idx < intp; idx++){
    d_firs[idx] = new gr::filter::kernel::fir_filter_ccf(1,dummy_taps);
  }

  //std::cout << d_indicator << ": FSO: " << d_fso << "\n";

  std::vector<float> shifted_taps;
  prototype_augemnt_fractional_delay(d_sps, 1., d_pulse_shape, d_fso, shifted_taps);

  //size_t leftover = (intp - (d_pulse_shape.size() % intp))%intp;
  //d_proto_taps = std::vector<float>(d_pulse_shape.size() + leftover, 0.);

  size_t leftover = (intp - (shifted_taps.size() % intp))%intp;
  d_proto_taps = std::vector<float>(shifted_taps.size() + leftover, 0.);

  memcpy( &d_proto_taps[0], &shifted_taps[0],
          shifted_taps.size()*sizeof(float) );

  if( d_proto_taps.size() % intp ){
    throw_runtime("signal_qam: error setting pulse shaping taps.\n");
  }

  //std::vector<float> shifted_taps = d_proto_taps;

  d_taps = std::vector< std::vector<float> >(intp);

  size_t ts = d_proto_taps.size() / intp;
  for(size_t idx = 0; idx < intp; idx++){
    d_taps[idx].resize(ts);
  }
  //printf("QAM:: taps init 0.\n");

  for(size_t idx = 0; idx < d_proto_taps.size(); idx++){
    d_taps[idx % intp][idx / intp] = d_proto_taps[idx];
  }
  //printf("QAM:: taps filled.\n");

  //printf("QAM:: filters made.\n");
  for(size_t idx = 0; idx < intp; idx++){
    d_firs[idx]->set_taps(d_taps[idx]);
  }
  //printf("QAM:: taps loaded.\n");
  d_hist = ts-1;
}


void
Signal_QAM::auto_gen_SYMS()
{
  size_t buff_size(0), buff_pnt(0);
  std::vector<complexf> buffer(d_buffer_size,complexf(0.,0.));
  while(d_running){
    for(size_t idx = 0; idx < d_buffer_size; idx++){
      //int data = rand()%d_order;
      int data = d_rng->ran_int();
      buffer[idx] = d_symbol_list[data];
      /*d_symbol_count++;
      if(buffer[idx].real()*buffer[idx].real()+buffer[idx].imag()*buffer[idx].imag() > 20.){
        printf("Auto_gen: Problem @ symbol %lu\n",d_symbol_count);
      }*/
    }

    while((buff_pnt < buffer.size()) && d_running){
      buff_pnt += d_Sy->bmemcpy( &buffer[buff_pnt], buffer.size()-buff_pnt, true );
    }
    buff_pnt = 0;
  }
}

void
Signal_QAM::throw_runtime(std::string err, size_t sc, size_t os, size_t si)
{
  //printf("Error occured while trying to make %lu samples, on sample %lu, symbol_idx %lu\n",sc,os, si);
  //printf("Error occured (QAM): A= %lu, B= %lu, C= %lu\n",sc,os, si);
  d_running = false;
  d_TGroup.join_all();
  delete d_Sy;
  for(size_t idx = 0; idx < d_sps; idx++){
    delete d_firs[idx];
  }
  printf("Threads ended.\n");
  throw std::runtime_error(err.c_str());
}
