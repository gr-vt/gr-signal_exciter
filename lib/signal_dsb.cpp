

#include "signal_dsb.hpp"
#include <stdio.h>//////////////////////////////////

Signal_DSB::Signal_DSB(float mod_idx, size_t components, float* mu,
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
  boost::mutex::scoped_lock scoped_lock(s_mutex_fftw);
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


  d_fso = fso;

  d_align = volk_get_alignment();
  // Generate and load the GNURadio FIR Filters with the pulse shape.
  load_fir();
  //printf("DSB::Loaded FIR filters.\n");

}

Signal_DSB::~Signal_DSB()
{
  delete d_fir;
  if(d_enable){
    d_running = false;
    d_TGroup.join_all();
    delete d_Sy;
  }
}

void
Signal_DSB::generate_symbols(complexf* output, size_t symbol_count)
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
Signal_DSB::generate_signal(complexf* output, size_t sample_count)
{
/*  complexf carrier(1.,0.);
  d_symbol_cache = std::vector<complexf>(sample_count, complexf(0.,0.));
  generate_symbols( &d_symbol_cache[0], sample_count );
  for(size_t idx = 0; idx < sample_count; idx++){
    output[idx] = d_symbol_cache[idx] + carrier;
  }
  if(d_norm){
    d_agc.scaleN( output, output, sample_count );
  }*/

  if(d_first_pass){
    d_past = std::vector<float>(d_hist);
    d_symbol_cache = std::vector<complexf>(d_hist);
    generate_symbols( &d_symbol_cache[0], d_hist );
    for(size_t idx = 0; idx < d_hist; idx++){
      d_past[idx] = d_symbol_cache[idx].real();
    }
    d_first_pass = false;
  }

  filter( sample_count, output );

  if(d_norm){
    d_agc.scaleN( output, output, sample_count );
  }
  //printf("Output[0] = (%1.3e, %1.3e)\n",output[0].real(), output[0].imag());
}

void
Signal_DSB::filter( size_t nout, complexf* out )
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
  //complexf carrier(1.,0.);
  float outsamp;

  while( oo < nout ){
    outsamp = d_fir->filter( &d_filt_in[ii++] );
    //out[oo++] = complexf( outsamp, 0. ) + carrier;
    out[oo++] = complexf( 1. + outsamp, 0. );
  }

  size_t remaining = total_input_len - ii;
  if(remaining < d_hist){
    fprintf(stderr,"DSB: nout(%lu), til(%lu), ii(%lu), past(%lu)\n",nout,total_input_len,ii,d_past.size());
    fprintf(stderr,"DSB - There isn't enough left in the buffer!!! (%lu,%lu)\n",remaining,d_hist);
  }
  d_past = std::vector<float>( &d_filt_in[ii], &d_filt_in[total_input_len] );

  volk_free(d_filt_in);

}


void
Signal_DSB::generate_taps()
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
}

void
Signal_DSB::load_fir()
{
  std::vector<float> dummy_taps;

  d_fir = new gr::filter::kernel::fir_filter_fff(1, dummy_taps);

  d_proto_taps = std::vector<float>( d_taps.begin(), d_taps.end() );
  time_offset(d_taps, d_proto_taps, d_fso);

  d_fir->set_taps(d_taps);
}

void
Signal_DSB::auto_fill_symbols()
{
  d_Sy = new signal_threaded_buffer<complexf>(d_buffer_size,d_notify_size);

  d_TGroup.create_thread( boost::bind(&Signal_DSB::auto_gen_GM, this) );

}

void
Signal_DSB::auto_fill_signal()
{}


void
Signal_DSB::auto_gen_GM()
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
