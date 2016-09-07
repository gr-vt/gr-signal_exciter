

#include "signal_fm.hpp"
#include <stdio.h>//////////////////////////////////

Signal_FM::Signal_FM(float mod_idx, float f_max, float var1, float var2, float thresh, int seed,
                      bool enable, size_t buff_size, size_t min_notify)
  : d_mod_idx(mod_idx),
    d_fmax(f_max),
    d_var1(var1),
    d_var2(var2),
    d_thresh(thresh),
    d_cum(0.),
    d_enable(enable),
    d_buffer_size(buff_size),
    d_notify_size(min_notify)
{
  set_seed(seed);
  generate_taps();
  d_gm = Gaussian_Mixture(sqrt(d_var1),sqrt(d_var2),d_fmax,d_thresh,d_seed);
  if(d_enable){
    d_running = true;
    auto_fill_symbols();
    auto_fill_signal();
  }
}

Signal_FM::~Signal_FM()
{
  if(d_enable){
    d_running = false;
    d_TGroup.join_all();
    delete d_Sy;
  }
}

void
Signal_FM::generate_symbols(complexf* output, size_t symbol_count)
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
    d_gm.get_message( &message[0], symbol_count );

    for(size_t idx = 0; idx < symbol_count; idx++){
      output[idx] = complexf(message[idx],0.);
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
Signal_FM::generate_signal(complexf* output, size_t sample_count)
{
  d_symbol_cache = std::vector<complexf>(sample_count, complexf(0.,0.));
  generate_symbols( &d_symbol_cache[0], sample_count );
  double mod = M_2_PIl*(d_mod_idx*d_fmax);
  for(size_t idx = 0; idx < sample_count; idx++){
    d_cum += d_symbol_cache[idx].real();
    output[idx] = exp(complexf(0.,mod*d_cum));
  }
}


void
Signal_FM::generate_taps()
{
  d_taps = std::vector<float>(51,0.);
  d_window = std::vector<float>(51,0.);
  for(size_t idx = 0; idx < 51; idx++)
  {
    if((idx-25) == 0)
    {
      d_taps[idx] = d_fmax;
    }
    else
    {
      d_taps[idx] = d_fmax*sin(M_PI*d_fmax*(float(idx)-25))/(M_PI*d_fmax*(float(idx)-25));
    }
    if(!((idx==0)||(idx==50)))
    {
      d_window[idx] = 0.42 - 0.5*cos((2*M_PI*float(idx))/50) + 0.08*cos((4*M_PI*float(idx))/50);
    }
    d_taps[idx] = d_taps[idx]*d_window[idx];
  }
}


void
Signal_FM::filter( float* n, float* m, size_t length )
{
  /*printf("d_m = [ %f",m[0]);
  for(size_t idx = 1; idx < length; idx++){
    printf(", %f",m[idx]);
  }
  printf("];\n");*/
  for(size_t idx = 0; idx < length; idx++){
    n[idx] = 0.;
    for(int ind_h(0), ind_m(idx); (ind_m >= 0) && (ind_h < (int)d_taps.size()); ind_m--, ind_h++){
      n[idx] += d_taps[ind_h]*m[ind_m];
      //printf("d_n(%lu,%lu) = %f;\n",idx+1,ind_h+1,d_taps[ind_h]*m[ind_m]);
      //printf("n(%lu), h(%d), m(%d)\n",idx,ind_h,ind_m);
    }
  }
  /*printf("d_n = [ %f",n[0]);
  for(size_t idx = 1; idx < length; idx++){
    printf(", %f",n[idx]);
  }
  printf("];\n");*/
}


void
Signal_FM::do_fft(size_t samp_count)
{
  memset( &d_fft_in[0], 0, sizeof(fftwf_complex)*samp_count );
  fftwf_execute( d_fft );
}

void
Signal_FM::do_ifft(size_t samp_count)
{
  memset( &d_fft_out[0], 0, sizeof(fftwf_complex)*samp_count );
  fftwf_execute( d_ifft );
  complexf scale = complexf(1/float(samp_count),0.);
  for(size_t idx = 0; idx < samp_count; idx++){
    d_output_cache[idx] = complexf(d_fft_out[idx][0],d_fft_out[idx][1])*scale;
  }
}

void 
Signal_FM::auto_fill_symbols()
{
  d_Sy = new signal_threaded_buffer<complexf>(d_buffer_size,d_notify_size);
  
  d_TGroup.create_thread( boost::bind(&Signal_FM::auto_gen_GM, this) );
  
}

void 
Signal_FM::auto_fill_signal()
{}


void
Signal_FM::auto_gen_GM()
{
  size_t buff_size(0), buff_pnt(0);
  std::vector<float> buffer(d_buffer_size,0.);
  std::vector<complexf> message(d_buffer_size,complexf(0.,0.));
  while(d_running){
    d_gm.get_message( &buffer[0], buffer.size() );
    for(size_t idx = 0.; idx < buffer.size(); idx++){
      message[idx] = complexf(buffer[idx],0.);
    }
    while((buff_pnt < buffer.size()) && d_running){
      buff_pnt += d_Sy->bmemcpy( &message[buff_pnt], message.size()-buff_pnt, true );
    }
    buff_pnt = 0;
  }
}





 
