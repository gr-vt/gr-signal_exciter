

#include "signal_lsb.hpp"
#include <stdio.h>//////////////////////////////////

Signal_LSB::Signal_LSB(float mod_idx, float f_max, float var1, float var2, float thresh, int seed,
                        bool norm, bool enable, size_t buff_size, size_t min_notify)
  : d_mod_idx(mod_idx),
    d_fmax(f_max),
    d_var1(var1),
    d_var2(var2),
    d_thresh(thresh),
    d_enable(enable),
    d_buffer_size(buff_size),
    d_notify_size(min_notify),
    d_agc(),
    d_norm(norm)
{
  set_seed(seed);
  generate_taps();
  d_gm = Gaussian_Mixture(sqrt(d_var1),sqrt(d_var2),d_fmax,d_thresh,d_seed);
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
}

Signal_LSB::~Signal_LSB()
{
  if(d_enable){
    d_running = false;
    d_TGroup.join_all();
    delete d_Sy;
  }
}

void
Signal_LSB::generate_symbols(complexf* output, size_t symbol_count)
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
Signal_LSB::generate_signal(complexf* output, size_t sample_count)
{
  boost::mutex::scoped_lock scoped_lock(s_mutex_fftw);
  d_symbol_cache = std::vector<complexf>(sample_count, complexf(0.,0.));
  d_output_cache = std::vector<complexf>(sample_count, complexf(0.,0.));
  d_fm = std::vector<float>(sample_count,0.);
  generate_symbols( &d_symbol_cache[0], sample_count );
  for(size_t idx = 0; idx < sample_count; idx++){
    d_fm[idx] = d_symbol_cache[idx].real();
  }
  /*printf("Signal LSB = [ %f",d_fm[0]);
  for(size_t idx = 1; idx < sample_count; idx++){
    printf(", %f",d_fm[idx]);
  }
  printf("];\n");*/

  d_fft_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*sample_count);
  d_fft_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*sample_count);
  d_fft = fftwf_plan_dft_r2c_1d( sample_count, &d_fm[0], d_fft_in, FFTW_ESTIMATE );
  d_ifft = fftwf_plan_dft_1d( sample_count, d_fft_in, d_fft_out, FFTW_BACKWARD, FFTW_ESTIMATE );

  do_fft(sample_count);
  for(size_t idx = 1; idx < sample_count/2; idx++){
    d_fft_in[idx][0] *= 2.;
    d_fft_in[idx][1] *= 2.;
  }
  for(size_t idx = sample_count/2+1; idx < sample_count; idx++){
    d_fft_in[idx][0] = 0.;
    d_fft_in[idx][1] = 0.;
  }
  for(size_t idx = 1, lidx=sample_count-1; idx < sample_count/2; idx++, lidx--){
    d_fft_in[lidx][0] = d_fft_in[idx][0];
    d_fft_in[lidx][1] = d_fft_in[idx][1];
    d_fft_in[idx][0] = 0.;
    d_fft_in[idx][1] = 0.;
  }
  d_fft_in[0][0] = 0.; d_fft_in[0][1] = 0.;
  do_ifft(sample_count);

  fftwf_free( d_fft_in );
  fftwf_free( d_fft_out );
  fftwf_destroy_plan( d_fft );
  fftwf_destroy_plan( d_ifft );

  memcpy( &output[0], &d_output_cache[0], sizeof(complexf)*sample_count);
  if(d_norm){
    d_agc.scaleN( output, output, sample_count );
  }
}


void
Signal_LSB::generate_taps()
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
      float x = (float(idx)-25)*d_fmax;
      d_taps[idx] = d_fmax*sin(M_PI*x)/(M_PI*x);
    }
    if(!((idx==0)||(idx==50)))
    {
      d_window[idx] = 0.42 - 0.5*cos((2*M_PI*float(idx))/50) + 0.08*cos((4*M_PI*float(idx))/50);
    }
    d_taps[idx] = d_taps[idx]*d_window[idx];
  }
}


void
Signal_LSB::filter( float* n, float* m, size_t length )
{
  for(size_t idx = 0; idx < length; idx++){
    n[idx] = 0.;
    for(int ind_h(0), ind_m(idx); (ind_m >= 0) && (ind_h < (int)d_taps.size()); ind_m--, ind_h++){
      n[idx] += d_taps[ind_h]*m[ind_m];
    }
  }
}


void
Signal_LSB::do_fft(size_t samp_count)
{
  memset( &d_fft_in[0], 0, sizeof(fftwf_complex)*samp_count );
  fftwf_execute( d_fft );
}

void
Signal_LSB::do_ifft(size_t samp_count)
{
  memset( &d_fft_out[0], 0, sizeof(fftwf_complex)*samp_count );
  fftwf_execute( d_ifft );
  complexf scale = complexf(1/float(samp_count),0.);
  for(size_t idx = 0; idx < samp_count; idx++){
    d_output_cache[idx] = complexf(d_fft_out[idx][0],d_fft_out[idx][1])*scale;
  }
}

void 
Signal_LSB::auto_fill_symbols()
{
  d_Sy = new signal_threaded_buffer<complexf>(d_buffer_size,d_notify_size);
  
  d_TGroup.create_thread( boost::bind(&Signal_LSB::auto_gen_GM, this) );
  
}

void 
Signal_LSB::auto_fill_signal()
{}



void
Signal_LSB::auto_gen_GM()
{
  size_t buff_size(0), buff_pnt(0);
  std::vector<float> buffer(d_buffer_size,0.);
  std::vector<complexf> message(d_buffer_size,complexf(0.,0.));
  while(d_running){
    d_gm.get_message( &buffer[0], buffer.size() );
    for(size_t idx = 0.; idx < buffer.size(); idx++){
      message[idx] = complexf(d_mod_idx*buffer[idx],0.);
    }
    while((buff_pnt < buffer.size()) && d_running){
      buff_pnt += d_Sy->bmemcpy( &message[buff_pnt], message.size()-buff_pnt, true );
    }
    buff_pnt = 0;
  }
}



  
