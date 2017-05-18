

#include "signal_cpm.hpp"
#include <stdio.h>//////////////////////////
#include <stdexcept>
#include <algorithm>

Signal_CPM::Signal_CPM(int order, gr::analog::cpm::cpm_type phase_type,
                      int sps, int overlap, float mod_idx, int seed,
                      double beta, float* phase_shape,
                      size_t phase_shape_length, bool enable_fso, float fso, bool enable,
                      size_t buff_size, size_t min_notify)
  : d_order(order),
    d_sps(sps),
    d_L(overlap),
    d_h(mod_idx),
    d_beta(beta),
    d_samp_offset(0),
    d_phase_acm(0.),
    d_enable(enable),
    d_buffer_size(buff_size),
    d_notify_size(min_notify)
{
  //printf("CPM INIT\n");
  get_indicator();
  set_seed(seed);

  d_first_pass = true;

  int new_order=0;
  while(order>0){
    order = order>>1;
    if(!new_order) new_order = 1;
    else new_order = new_order*2;
  }

  if(new_order != d_order){
    if(new_order > 0) d_order = new_order;
    else d_order = 2;
  }

  //printf("CPM ORDER\n");

  enable_fractional_offsets(enable_fso, fso);

  create_symbol_list();

  d_rng = new gr::random(d_seed, 0, d_order);

  //printf("CPM SYMBOLS\n");

  if(phase_type != gr::analog::cpm::GENERIC){
    d_phase_shape = gr::analog::cpm::phase_response(phase_type, d_sps, d_L, d_beta);
  }
  else{
    d_phase_shape = std::vector<float>(phase_shape_length,0.);
    for(size_t idx = 0; idx < phase_shape_length; idx++){
      d_phase_shape[idx] = phase_shape[idx];
    }
    //d_phase_shape = std::vector<float>(phase_shape, phase_shape+phase_shape_length);
  }
  d_cpm_type = phase_type;

  float weighting = 0.;
  for(size_t idx = 0; idx < d_phase_shape.size(); idx++){
    weighting += d_phase_shape[idx];
  }
  if(weighting){
    weighting = 1./weighting;
  }
  else{
    weighting = 1.;
  }
  for(size_t idx = 0; idx < d_phase_shape.size(); idx++){
    d_phase_shape[idx] *= weighting;
  }

  //printf("CPM PHASE SHAPE\n");

  d_align = volk_get_alignment();

  load_firs();

  //printf("CPM FIRS\n");

  if(d_enable){
    d_running = true;
    auto_fill_symbols();
    auto_fill_signal();
  }

  //printf("CPM END\n");
}

Signal_CPM::~Signal_CPM()
{
  if(d_enable){
    d_running = false;
    d_TGroup.join_all();
    delete d_Sy;
    //for(size_t idx = 0; idx < d_sps; idx++){
    //  delete d_firs[idx];
    //}
  }
  delete d_rng;
  delete d_frac_filt;
}

void
Signal_CPM::generate_symbols(complexf* output, size_t symbol_count)
{
  //printf("CPM:: SYM, COUNT %lu\n", symbol_count);
  if(d_enable){
    size_t filled(0);
    while((filled < symbol_count) && d_running){
      filled += d_Sy->bmemcpy( &output[filled], symbol_count-filled, false );
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
Signal_CPM::generate_signal(complexf* output, size_t sample_count)
{
  if(d_first_pass){
    d_past = std::vector<float>(d_hist);
    std::vector<complexf> temp(d_hist);
    generate_symbols( &temp[0], d_hist );
    for(size_t idx = 0; idx < d_hist; idx++){
      d_past[idx] = temp[idx].real();
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

}

void
Signal_CPM::create_symbol_list()
{
  //printf("Creating the symbol list\n");
  int R(d_order);


  std::vector<float> reals(R,0.);

  for(int idx = 0; idx < R; idx++){
    reals[idx] = idx*2-R+1;
  }

  d_symbol_amps = std::vector<float>(reals.begin(),reals.end());
  d_symbol_list = std::vector<complexf>(R);
  for(size_t idx = 0; idx < d_symbol_amps.size(); idx++){
    d_symbol_list[idx] = complexf(d_symbol_amps[idx],0.);
  }

  //if(d_enable_fractional){
    size_t tap_len = (d_sps%2) ? 11*d_sps : 11*d_sps+1;
    d_proto_taps = std::vector<float>(tap_len,0.);
    d_proto_taps[(tap_len-1)/2] = 1.;
  //}
  //else{
  //  d_proto_taps = std::vector<float>(1,1.);
  //}
}

void
Signal_CPM::filter( size_t nout, complexf* out )
{
  //printf("CPM:: FILTER, %lu\n", nout);

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
  //printf("CPM:: INPUT LEN, %lu\n", total_input_len);

  d_filt_in = (float*) volk_malloc( total_input_len*sizeof(float), d_align );
  d_filt_out = (float*) volk_malloc( d_sps*total_input_len*sizeof(float), d_align );
  memcpy( &d_filt_in[0], &d_past[0], d_past.size()*sizeof(float) );

  complexf* symb_buff = (complexf*) volk_malloc( symbs_needed*sizeof(complexf), d_align );
  //generate_symbols( &d_filt_in[d_past.size()], symbs_needed );
  generate_symbols( &symb_buff[0], symbs_needed );
  //printf("CPM:: SYMS_GEN'd\n");
  volk_32fc_deinterleave_real_32f( &d_filt_in[d_past.size()], &symb_buff[0], symbs_needed );
  //printf("CPM:: DEINTERLEAVED\n");
  volk_free(symb_buff);
  //printf("CPM:: VOLK FREED 2\n");

  size_t oo = 0, ii = 0;

  double oi,oq;

  //printf("[");
  //fflush(stdout);

  while( oo < nout ){
    d_filt_out[oo] = d_firs[d_samp_offset]->filter( &d_filt_in[ii] );
    //printf(" %1.3e", d_filt_out[oo]);
    //fflush(stdout);
    //pulled from gr::analog::fm
    d_phase_acm += d_filt_out[oo] * M_PI * d_h;
    //d_phase_acm = std::fmod( d_phase_acm + M_PI, M_2_PIl ) - M_PI;
    gr::sincos(d_phase_acm, &oi, &oq);
    out[oo] = complexf(oi,oq);


    d_samp_offset = (d_samp_offset+1)%d_sps;
    if(d_samp_offset == 0){
      ii++;
    }
    oo++;
  }
  //printf("]\n");
  size_t remaining = total_input_len - ii;
  if((remaining < d_hist) || (total_input_len < ii)){
    fprintf(stderr,"CPM - There isn't enough left in the buffer!!!\n");
  }
  //printf("CPM:: FILTERED, REMAINS %lu\n", remaining);
  d_past = std::vector<float>( &d_filt_in[ii], &d_filt_in[total_input_len] );
  //printf("CPM:: PAST %lu\n",d_past.size());

  volk_free(d_filt_in);
  //printf("CPM:: VOLK FREED 0\n");
  volk_free(d_filt_out);
  //printf("CPM:: VOLK FREED 1\n");
}

void
Signal_CPM::auto_fill_symbols()
{
  //
  d_Sy = new signal_threaded_buffer<complexf>(d_buffer_size, d_notify_size);
  // Start the symbol generation thread
  d_TGroup.create_thread( boost::bind(&Signal_CPM::auto_gen_SYMS, this) );
}

void
Signal_CPM::auto_fill_signal()
{}

void
Signal_CPM::load_firs()
{
  size_t intp = d_sps;
  size_t ts = d_phase_shape.size() / intp;
  if((d_phase_shape.size() % intp)){
    ts++;
  }

  d_taps = std::vector< std::vector<float> >(intp);


  for(size_t idx = 0; idx < intp; idx++){
    d_taps[idx].resize(ts);
    memset( &d_taps[idx][0], 0., ts*sizeof(float) );
  }

  for(size_t idx = 0; idx < d_phase_shape.size(); idx++){
    d_taps[idx % intp][idx / intp] = d_phase_shape[idx];
  }

  d_firs = std::vector< gr::filter::kernel::fir_filter_fff *>(intp);
  std::vector<float> dummy_taps(1,0.);
  for(size_t idx = 0; idx < intp; idx++){
    d_firs[idx] = new gr::filter::kernel::fir_filter_fff(1,dummy_taps);
  }
  for(size_t idx = 0; idx < intp; idx++){
    d_firs[idx]->set_taps(d_taps[idx]);
  }
  d_hist = ts-1;

  //std::cout << d_indicator << ": FSO: " << d_fso << "\n";

  d_frac_filt = new gr::filter::kernel::fir_filter_ccf(1, dummy_taps);
  std::vector<float> shifted_taps;
  prototype_augemnt_fractional_delay(1., d_sps, d_proto_taps, d_fso, shifted_taps);
  //std::vector<float> shifted_taps = d_proto_taps;
  d_frac_filt->set_taps(shifted_taps);
  d_frac_cache = std::vector<complexf>(shifted_taps.size()-1);
}

void
Signal_CPM::auto_gen_SYMS()
{
  size_t buff_size(0), buff_pnt(0);
  std::vector<complexf> buffer(d_buffer_size,complexf(0.,0.));
  while(d_running){
    for(size_t idx = 0; idx < d_buffer_size; idx++){
      //int data = rand()%d_order;
      int data = d_rng->ran_int();
      buffer[idx] = d_symbol_list[data];
    }

    while((buff_pnt < buffer.size()) && d_running){
      buff_pnt += d_Sy->bmemcpy( &buffer[buff_pnt], buffer.size()-buff_pnt, true );
    }
    buff_pnt = 0;
  }
}

void
Signal_CPM::throw_runtime(std::string err)
{
  //printf("Error occured while trying to make %lu samples, on sample %lu, symbol_idx %lu\n",sc,os, si);
  //printf("Error occured (PSK): A= %lu, B= %lu, C= %lu\n",sc,os, si);
  d_running = false;
  d_TGroup.join_all();
  delete d_Sy;
  //for(size_t idx = 0; idx < d_sps; idx++){
  //  delete d_firs[idx];
  //}
  //printf("Threads ended.\n");
  throw std::runtime_error(err.c_str());
}
