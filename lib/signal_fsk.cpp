

#include "signal_fsk.hpp"
#include <stdio.h>//////////////////////////
#include <stdexcept>
#include <algorithm>

Signal_FSK::Signal_FSK(int order, int sps, float mod_idx, int seed,
                      bool enable_fso, float fso, bool enable,
                      size_t buff_size, size_t min_notify)
  : d_order(order),
    d_sps(sps),
    d_h(mod_idx),
    d_n(0),
    d_samp_offset(0),
    d_enable(enable),
    d_buffer_size(buff_size),
    d_notify_size(min_notify)
{
  //printf("FSK INIT\n");
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

  //printf("FSK ORDER\n");

  enable_fractional_offsets(enable_fso, fso);

  d_f = (d_h/float(d_sps))/2.;

  create_symbol_list();

  d_rng = new gr::random(d_seed, 0, d_order);

  //printf("FSK SYMBOLS\n");

  d_align = volk_get_alignment();

  load_firs();

  //printf("FSK FIRS\n");


  //printf("Order(%d)\nSPS(%d)\nh(%1.3e)\nd_f(%1.3e)\n",d_order,d_sps,d_h,d_f);

  if(d_enable){
    d_running = true;
    auto_fill_symbols();
    auto_fill_signal();
  }

  //printf("FSK END\n");
}

Signal_FSK::~Signal_FSK()
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
Signal_FSK::generate_symbols(complexf* output, size_t symbol_count)
{
  //printf("FSK:: SYM, COUNT %lu\n", symbol_count);
  if(d_enable){
    size_t filled(0);
    while((filled < symbol_count) && d_running){
      filled += d_Sy->bmemcpy( &output[filled], symbol_count-filled, false );
    }
  }
  else{
    for(size_t idx = 0; idx < symbol_count; idx++){
      int data = d_rng->ran_int();
      output[idx] = complexf(d_symbol_amps[data],0.);
    }
  }
}

void
Signal_FSK::generate_signal(complexf* output, size_t sample_count)
{
  if(d_first_pass){
    d_symbol_cache = std::vector<complexf>(0);
    int needed = int(std::ceil(float(d_frac_cache.size())/float(d_sps)));
    std::vector<complexf> temp(needed);
    d_n = -needed*d_sps;
    filter( needed, &temp[0] );
    memcpy( &d_frac_cache[0], &temp[needed*d_sps-d_frac_cache.size()],
            d_frac_cache.size()*sizeof(complexf) );
    d_n = 0;
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
Signal_FSK::create_symbol_list()
{
  //printf("Creating the symbol list\n");
  int R(d_order);


  std::vector<float> reals(R,0.);

  for(int idx = 0; idx < R; idx++){
    reals[idx] = d_f*(idx*2-R+1);
  }
  /*printf("Symbs = [");
  for(size_t idx = 0; idx < R; idx++){
    printf(" %1.3e",reals[idx]);
  }
  printf("];\n");*/

  d_symbol_amps = std::vector<float>(reals.begin(),reals.end());

  size_t tap_len = (d_sps%2) ? 11*d_sps : 11*d_sps+1;
  d_proto_taps = std::vector<float>(tap_len,0.);
  d_proto_taps[(tap_len-1)/2] = 1.;
}

void
Signal_FSK::filter( size_t nout, complexf* out )
{
  double angle = 0., oi, oq;
  int oo(0);
  int remain = (d_sps - d_samp_offset)%d_sps;
  while((d_samp_offset != 0)&&(oo<nout)){
    angle = (2.*M_PI*d_n)*(d_symbol_cache[0].real());
    gr::sincos(angle, &oq, &oi);
    out[oo] = complexf(oi,oq);
    oo++;
    d_samp_offset = (d_samp_offset+1)%d_sps;
  }

  int needed = int(std::ceil(float(nout-remain)/float(d_sps)));
  std::vector<complexf> temp(needed);
  generate_symbols( &temp[0], needed );

  size_t idx(0);
  while(oo < nout){
    angle = (2.*M_PI*d_n)*(temp[idx].real());
    gr::sincos(angle, &oq, &oi);
    out[oo] = complexf(oi,oq);
    d_n++;
    oo++;
    d_samp_offset = (d_samp_offset+1)%d_sps;
    if(d_samp_offset == 0){
      idx++;
    }
  }
  if(d_samp_offset){
    d_symbol_cache = std::vector<complexf>( &temp[idx], &temp[idx+1] );
  }
  else{
    d_symbol_cache = std::vector<complexf>(0);
  }
}

void
Signal_FSK::auto_fill_symbols()
{
  //
  d_Sy = new signal_threaded_buffer<complexf>(d_buffer_size, d_notify_size);
  // Start the symbol generation thread
  d_TGroup.create_thread( boost::bind(&Signal_FSK::auto_gen_SYMS, this) );
}

void
Signal_FSK::auto_fill_signal()
{}

void
Signal_FSK::load_firs()
{
  std::vector<float> dummy_taps(0);
  //std::cout << d_indicator << ": FSO: " << d_fso << "\n";

  d_frac_filt = new gr::filter::kernel::fir_filter_ccf(1, dummy_taps);
  std::vector<float> shifted_taps;
  prototype_augemnt_fractional_delay(1., d_sps, d_proto_taps, d_fso, shifted_taps);
  //std::vector<float> shifted_taps = d_proto_taps;
  d_frac_filt->set_taps(shifted_taps);
  d_frac_cache = std::vector<complexf>(shifted_taps.size()-1);
}

void
Signal_FSK::auto_gen_SYMS()
{
  size_t buff_size(0), buff_pnt(0);
  std::vector<complexf> buffer(d_buffer_size,complexf(0.,0.));
  while(d_running){
    for(size_t idx = 0; idx < d_buffer_size; idx++){
      int data = d_rng->ran_int();
      buffer[idx] = complexf(d_symbol_amps[data],0.);
    }

    while((buff_pnt < buffer.size()) && d_running){
      buff_pnt += d_Sy->bmemcpy( &buffer[buff_pnt], buffer.size()-buff_pnt, true );
    }
    buff_pnt = 0;
  }
}

void
Signal_FSK::throw_runtime(std::string err)
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
