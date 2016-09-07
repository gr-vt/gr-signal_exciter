

#include "signal_psk.hpp"
#include <stdio.h>//////////////////////////////////
#include <stdexcept>
#include <algorithm>

Signal_PSK::Signal_PSK(int order, float offset, int sps, float* pulse_shape, size_t length, int seed,
                        bool enable, size_t buff_size, size_t min_notify)
  : d_order(order),
    d_offset(offset),
    d_sps(sps),
    d_samp_offset(0),
    d_protected(0),
    d_enable(enable),
    d_buffer_size(buff_size),
    d_notify_size(min_notify),
    d_sample_count(0),
    d_symbol_count(0),
    d_samp_gen_count(0)
{
  //printf("PSK::Init.\n");
  set_seed(seed);
  //printf("PSK::Seeded.\n");

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
  //printf("PSK::Filled Pulse Shape.\n");

  // Checking for how many symbols will overlap
  float check = float(length)/float(sps);
  d_overlap = size_t(int(floor(check)));
  // If there is a factional overlap, the maximum overlap is +1
  if(floor(check) < check)
  {
    d_overlap++;
  }
  //printf("Overlaped.\n");

  d_first_pass = true;

  // Force a proper modulation order
  int new_order=0;

  while(order>0){
    order = order>>1;
    // Init new order
    if(!new_order) new_order = 1;
    // Grow new order
    else new_order = new_order*2;
  }
  //printf("New Order.\n");
  
  // If not equal, impose the new_order
  if(new_order != d_order){
    if(new_order > 0) d_order = new_order;
    else d_order = 2;
  }
  //printf("Ok.\n");

  // Generate the symbol look up table
  create_symbol_list();
  //printf("PSK::Symbols Created.\n");

  d_rng = new gr::random(d_seed, 0, d_order);

  // Check that there is a tap for the number of samples per symbol.
  if(d_sps > d_pulse_shape.size()){
    printf("The pulse_shape is shorter than sps, this will crash.\n");
  }

  // Enable background threads for signal generation parallel processing.
  if(d_enable){
    d_align = volk_get_alignment();
    // Generate and load the GNURadio FIR Filters with the pulse shape.
    load_firs();
    //printf("PSK::Loaded FIR filters.\n");

    // Start the running flag for the threads to exit on shutdown
    d_running = true;
    auto_fill_symbols();
    auto_fill_signal();
  }
  //printf("PSK::Fully made\n");

  // Everything has started, now just need a consumer downstream
}

Signal_PSK::~Signal_PSK()
{
  // Threads are alive, kill them and clean up
  if(d_enable){
    d_running = false;/*
    std::vector<complexf> null_buff(d_buffer_size,complexf(0.,0.));
    d_Sy->bmemcpy( &null_buff[0], d_buffer_size, true );
    d_Sy->bmemcpy( &null_buff[0], d_buffer_size, false );
    d_Sy->cleanup();*/
    d_TGroup.join_all();
    delete d_Sy;
    for(size_t idx = 0; idx < d_sps; idx++){
      delete d_firs[idx];
    }
  }
  delete d_rng;
}

void
Signal_PSK::generate_symbols(complexf* output, size_t symbol_count)
{
  d_symb_gen_count++;
  // If using threading
  if(d_enable){
    // Init how many symbols have been consumed
    size_t filled(0);
    // Until all requested symbols are made pull more symbols, conditioned on the system is running
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
    // Fall back method, generate all symbols on demand (slower execution, less load)
    for(size_t idx = 0; idx<symbol_count; idx++){
      //int data = rand()%d_order;
      int data = d_rng->ran_int();
      output[idx] = d_symbol_list[data];
    }
  }
}

void
Signal_PSK::generate_signal(complexf* output, size_t sample_count)
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
Signal_PSK::create_symbol_list()
{
  //printf("Creating the symbol list\n");
  // The symbol look up table is length of the modulation order
  d_symbol_list = std::vector<complexf>(d_order);
  for(int data = 0; data < d_order; data++){
    d_symbol_list[data] = exp( complexf( 0., 2*M_PI*float(data)/d_order + d_offset ) );
  }
  //printf("Created the symbol list\n");
}



void
Signal_PSK::filter( size_t nout, complexf* out )
{
/*////////////////////DEBUGGING
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
Signal_PSK::auto_fill_symbols()
{
  // Create a buffer for the generation thread to fill, and consumer to consumer from
  d_Sy = new signal_threaded_buffer<complexf>(d_buffer_size,d_notify_size);
  // Start the symbol generation thread
  d_TGroup.create_thread( boost::bind(&Signal_PSK::auto_gen_SYMS, this) );
}

void 
Signal_PSK::auto_fill_signal()
{//Currently not doing anything
}


void
Signal_PSK::load_firs()
{
  size_t interp = size_t(d_sps);
  
  std::vector<float> dummy_taps;

  d_firs = std::vector< gr::filter::kernel::fir_filter_ccf* >(interp);

  for(size_t samp_offset = 0; samp_offset < interp; samp_offset++){
    d_firs[samp_offset] = new gr::filter::kernel::fir_filter_ccf(1, dummy_taps);
  }

  std::vector<float> stretch = d_pulse_shape;
  size_t increase = stretch.size() % interp;
  if(increase > 0){
    increase = interp-increase;
    while(increase-- > 0){
      stretch.insert( stretch.end(), 0. );
    }
  }

  if( stretch.size() % interp != 0 ){
    throw_runtime("signal_psk: error setting pulse shaping taps.\n");
  }

  d_taps = std::vector< std::vector<float> >(interp);

  size_t taps_per_branch = stretch.size()/interp;
  for(size_t branch = 0; branch < interp; branch++){
    d_taps[branch].resize(taps_per_branch);
  }

  for(size_t tap_idx = 0; tap_idx < stretch.size(); tap_idx++){
    d_taps[tap_idx % interp][tap_idx/interp] = stretch[tap_idx];
  }

  for(size_t branch = 0; branch < interp; branch++){
    d_firs[branch]->set_taps(d_taps[branch]);
  }

  d_hist = taps_per_branch-1;

}


void
Signal_PSK::auto_gen_SYMS()
{
  // Thread for generating symbols and pushing them to d_Sy buffer
  size_t buff_size(0), buff_pnt(0);
  // Using local buffer of max size d_buffer_size
  std::vector<complexf> buffer(d_buffer_size,complexf(0.,0.));
  while(d_running){
    // While the system is running
    for(size_t idx = 0; idx < d_buffer_size; idx++){
      // Generate a integer for the symbol buffer to generate
      //int data = rand()%d_order;
      int data = d_rng->ran_int();
      // Load the local buffer
      buffer[idx] = d_symbol_list[data];

      /*d_symbol_count++;
      if(buffer[idx].real()*buffer[idx].real()+buffer[idx].imag()*buffer[idx].imag() > 20.){
        printf("Auto_gen: Problem @ symbol %lu\n",d_symbol_count);
      }*/
    }
    // Local buffer is full, fill the shared buffer

    // While local buffer still has data, fill the shared buffer, conditioned on system still running
    while((buff_pnt < buffer.size()) && d_running){
      // Copy over the local buffer, increment by how much was actually copied
      buff_pnt += d_Sy->bmemcpy( &buffer[buff_pnt], buffer.size()-buff_pnt, true );
    }
    // Reset local pointer
    buff_pnt = 0;
  }
}

void
Signal_PSK::throw_runtime(std::string err, size_t sc, size_t os, size_t si)
{
  //printf("Error occured while trying to make %lu samples, on sample %lu, symbol_idx %lu\n",sc,os, si);
  //printf("Error occured (PSK): A= %lu, B= %lu, C= %lu\n",sc,os, si);
  d_running = false;
  d_TGroup.join_all();
  delete d_Sy;
  for(size_t idx = 0; idx < d_sps; idx++){
    delete d_firs[idx];
  }
  printf("Threads ended.\n");
  throw std::runtime_error(err.c_str());
}






  







