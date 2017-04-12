//#################################################################
//# GNU Radio C++ Module
//# Title: Signal Exciter File
//#################################################################
/* 
 * Copyright 2016 Bill Clark.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


// General Includes
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <exception>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <libconfig.h>

// GNU Radio Includes
#include <gnuradio/top_block.h>
#include <gnuradio/filter/firdes.h>
#include <gnuradio/filter/pm_remez.h>
#include <gnuradio/filter/interp_fir_filter_ccf.h>
#include <gnuradio/filter/pfb_arb_resampler_ccf.h>
#include <gnuradio/filter/freq_xlating_fir_filter_ccf.h>
#include <gnuradio/filter/pfb_channelizer_ccf.h>
#include <gnuradio/filter/pfb_synthesizer_ccf.h>
#include <gnuradio/filter/fir_filter_ccf.h>
#include <gnuradio/blocks/multiply_const_cc.h>
#include <gnuradio/blocks/add_cc.h>
#include <gnuradio/blocks/stream_to_streams.h>
#include <gnuradio/blocks/null_source.h>
#include <gnuradio/blocks/file_sink.h>
#include <gnuradio/blocks/head.h>
#include <gnuradio/blocks/skiphead.h>
#include <gnuradio/blocks/null_sink.h>
#include <gnuradio/blocks/throttle.h>
#include <gnuradio/blocks/multiply_cc.h>
#include <gnuradio/channels/channel_model.h>
#include <gnuradio/analog/cpm.h>
#include <gnuradio/analog/sig_source_c.h>

// OOTM Includes
#include <signal_exciter/random_signal_config.h>
#include <signal_exciter/random_signal.h>
#include <signal_exciter/one_pass_gate.h>
#include <signal_exciter/periodic_gate.h>
#include <signal_exciter/random_gate.h>

// Boost Includes
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/sinc.hpp>

std::vector<float> gen_arb_interp_taps(double rate);
std::vector<float> gen_int_interp_taps(double rate);
std::vector<float> gen_synthz_taps(int channels, float samp_rate = 1.0);
std::vector<float> gen_channl_taps(int channels, float samp_rate = 1.0);
std::vector<double> kaiser_window(size_t ntaps, double beta);

#define SYNTH_SIZE 64
#define TAP_SWITCH 1
#define SIGNAL_MAX 20

//Used if TAP_SWITCH 0
//Do not use just yet
#define TAPS_PER_CHAN 10
#define SYNTH_BETA 9.5
#define CHANN_BETA 9.5

#define OUTPUT_ENABLE_PARSER    1
#define OUTPUT_ENABLE_GENERATE  1
#define OUTPUT_ENABLE_RUN       1

#if OUTPUT_ENABLE_PARSER
#define Pprintf(s, c, ...) printf("\033[38;5;%lum" s "\033[0m\n", c, __VA_ARGS__)
#else
#define Pprintf(s, c, ...)
#endif
#if OUTPUT_ENABLE_GENERATE
#define Cprintf(s, c, ...) printf("\033[38;5;%lum" s "\033[0m\n", c, __VA_ARGS__)
#else
#define Cprintf(s, c, ...)
#endif
#if OUTPUT_ENABLE_RUN
#define Rprintf(s, c, ...) printf("\033[38;5;%lum" s "\033[0m\n", c, __VA_ARGS__)
#else
#define Rprintf(s, c, ...)
#endif


//////////////////////////// DEBUGGING ///////////////////////////
#define DEBUG_OUTPUT 0
//////////////////////////////////////////////////////////////////

struct system_var
{
  bool EnableNoise        ;
  int NumSignals          ;
  int Seed                ;
  double SysBW            ;
  double CenterFreq       ;
  double RunTime          ;
  double NoiseFloor       ;
  std::string outputFile  ;
};


int parse_config(std::string config_file, system_var* system_container, gr::signal_exciter::sig_params* signal_container, size_t max_signals=SIGNAL_MAX)
{
  config_t cfg, *cfp;

  cfp = &cfg;
  config_init(cfp);

  if(!config_read_file(cfp, config_file.c_str())){
    fprintf(stderr, "Failed to read the config file: %s\n",config_file.c_str());
    config_destroy(cfp);
    return 1;
  }

  // Check for the seed value
  if(!config_lookup_int(cfp, "Seed", &(*system_container).Seed)){
    //Failed to find keyword
    system_container->Seed = -1;
  }
  Pprintf("Seed: %d", size_t(75),system_container->Seed);

  // Check for the number of signals to expect
  if(!config_lookup_int(cfp, "Nsignals", &(*system_container).NumSignals)){
    //Failed to find keyword
    system_container->NumSignals = 0;
  }
  Pprintf("Nsignals: %d", size_t(75),system_container->NumSignals);
  if(system_container->NumSignals > max_signals){
    config_destroy(cfp);
    fprintf(stderr, "The config file requests more signals than allowed.\n\tRemake with higher allowance, or cut back on the number of signals.\n");
    return 1;
  }

  // Check for the center frequency to be emulated at
  if(!config_lookup_float(cfp, "CenterFreq", &(*system_container).CenterFreq)){
    //Failed to find keyword
    system_container->CenterFreq = 0;
  }
  Pprintf("CenterFreq: %1.3e", size_t(75),system_container->CenterFreq);
  // Check for the system bandwidth to operate at
  if(!config_lookup_float(cfp, "SysBW", &(*system_container).SysBW)){
    //Failed to find keyword
    system_container->SysBW = 0;
  }
  Pprintf("SysBW: %1.3e", size_t(75),system_container->SysBW);
  // Check for how long the file should be in seconds
  if(!config_lookup_float(cfp, "RunTime", &(*system_container).RunTime)){
    //Failed to find keyword
    system_container->RunTime = 0;
  }
  Pprintf("RunTime: %1.3e", size_t(75),system_container->RunTime);
  // Check for the noise floor level in dBm/Hz
  if(!config_lookup_float(cfp, "NoiseFloor", &(*system_container).NoiseFloor)){
    //Failed to find keyword
    system_container->NoiseFloor = -100.;
  }
  Pprintf("NoiseFloor: %1.3e", size_t(75),system_container->NoiseFloor);

  // Check if noise should be added to the file
  int EnableNoise;
  if(!config_lookup_bool(cfp, "NoiseEnable", &EnableNoise)){
    //Failed to find keyword
    system_container->EnableNoise = false;
  }
  else{
    //Found Keyword
    system_container->EnableNoise = EnableNoise;
  }
  Pprintf("NoiseEnable: %d", size_t(75),system_container->EnableNoise);

  // Check for the output file name
  const char* outputFile;
  if(!config_lookup_string(cfp, "outfile", &outputFile)){
    //Failed to find keyword
    system_container->outputFile = "SignalExciterGenerated.fc32";
  }
  else{
    //Found Keyword
    system_container->outputFile = outputFile;
  }
  Pprintf("outfile: %s", size_t(75),system_container->outputFile.c_str());


  for(int sig = 1; sig <= system_container->NumSignals; sig++){
    char signum[15];
    sprintf(signum, "sig%d", sig);
    const char* signumber = signum;

    char str[50];

    // Checking the modulation type
    strcpy(str,signumber);
    strcat(str,".type");
    const char* type;
    if(!config_lookup_string(cfp, str, &type)){
      //Failed to find keyword
      signal_container[sig-1].type = gr::signal_exciter::PSK;
    }
    else{
      //Found Keyword
      std::string type_lower = type;
      std::transform( type_lower.begin(), type_lower.end(), type_lower.begin(), ::tolower );
      if(strcmp(type_lower.c_str(),"psk")==0){
        signal_container[sig-1].type = gr::signal_exciter::PSK;
      }
      else if(strcmp(type_lower.c_str(),"qam")==0){
        signal_container[sig-1].type = gr::signal_exciter::QAM;
      }
      else if(strcmp(type_lower.c_str(),"pam")==0){
        signal_container[sig-1].type = gr::signal_exciter::PAM;
      }
      else if(strcmp(type_lower.c_str(),"ask")==0){
        signal_container[sig-1].type = gr::signal_exciter::ASK;
      }
      else if(strcmp(type_lower.c_str(),"fsk")==0){
        signal_container[sig-1].type = gr::signal_exciter::FSK;
      }
      else if(strcmp(type_lower.c_str(),"cwmorse")==0){
        signal_container[sig-1].type = gr::signal_exciter::CWMORSE;
      }
      else if(strcmp(type_lower.c_str(),"dsb")==0){
        signal_container[sig-1].type = gr::signal_exciter::DSB;
        signal_container[sig-1].am_norm = true;
      }
      else if(strcmp(type_lower.c_str(),"dsbsc")==0){
        signal_container[sig-1].type = gr::signal_exciter::DSBSC;
        signal_container[sig-1].am_norm = true;
      }
      else if(strcmp(type_lower.c_str(),"usb")==0){
        signal_container[sig-1].type = gr::signal_exciter::USB;
        signal_container[sig-1].am_norm = true;
      }
      else if(strcmp(type_lower.c_str(),"lsb")==0){
        signal_container[sig-1].type = gr::signal_exciter::LSB;
        signal_container[sig-1].am_norm = true;
      }
      else if(strcmp(type_lower.c_str(),"fm")==0){
        signal_container[sig-1].type = gr::signal_exciter::FM;
      }
      else if(strcmp(type_lower.c_str(),"cpm")==0){
        signal_container[sig-1].type = gr::signal_exciter::CPM;
      }
      else if(strcmp(type_lower.c_str(),"gmsk")==0){
        signal_container[sig-1].type = gr::signal_exciter::GMSK;
      }
      else if(strcmp(type_lower.c_str(),"msk")==0){
        signal_container[sig-1].type = gr::signal_exciter::MSK;
      }
      else if(strcmp(type_lower.c_str(),"gfsk")==0){
        signal_container[sig-1].type = gr::signal_exciter::GFSK;
      }
      else if(strcmp(type_lower.c_str(),"ofdm")==0){
        signal_container[sig-1].type = gr::signal_exciter::OFDM;
      }
      else{
        config_destroy(cfp);
        fprintf(stderr, "Unknown 'type' : signal %d : %s\n", sig, type);
        return 1;
      }
    }
    Pprintf("Signal %d, type: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].type);

    // Checking the secondary modulation type (OFDM primary)
    strcpy(str,signumber);
    strcat(str,".mod");
    const char* mod;
    if(!config_lookup_string(cfp, str, &mod)){
      //Failed to find keyword
      signal_container[sig-1].mod = gr::signal_exciter::PSK;
    }
    else{
      //Found Keyword
      std::string mod_lower = mod;
      std::transform( mod_lower.begin(), mod_lower.end(), mod_lower.begin(), ::tolower );
      if(strcmp(mod_lower.c_str(),"psk")==0){
        signal_container[sig-1].mod = gr::signal_exciter::PSK;
      }
      else if(strcmp(mod_lower.c_str(),"qam")==0){
        signal_container[sig-1].mod = gr::signal_exciter::QAM;
      }
      else if(strcmp(mod_lower.c_str(),"pam")==0){
        signal_container[sig-1].mod = gr::signal_exciter::PAM;
      }
      else if(strcmp(mod_lower.c_str(),"ask")==0){
        signal_container[sig-1].mod = gr::signal_exciter::ASK;
      }
      else{
        config_destroy(cfp);
        fprintf(stderr, "Unknown 'mod' : signal %d : %s\n", sig, mod);
        return 1;
      }
    }
    Pprintf("Signal %d, mod: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].mod);

    // Checking the phase shaping type (CPM...)
    strcpy(str,signumber);
    strcat(str,".ptype");
    const char* ptype;
    if(!config_lookup_string(cfp, str, &ptype)){
      //Failed to find keyword
      signal_container[sig-1].phase_type = gr::analog::cpm::LREC;
    }
    else{
      //Found Keyword
      std::string ptype_lower = ptype;
      std::transform( ptype_lower.begin(), ptype_lower.end(), ptype_lower.begin(), ::tolower );
      if(strcmp(ptype_lower.c_str(),"lrec")==0){
        signal_container[sig-1].phase_type = gr::analog::cpm::LREC;
      }
      else if(strcmp(ptype_lower.c_str(),"lrc")==0){
        signal_container[sig-1].phase_type = gr::analog::cpm::LRC;
      }
      else if(strcmp(ptype_lower.c_str(),"lsrc")==0){
        signal_container[sig-1].phase_type = gr::analog::cpm::LSRC;
      }
      else if(strcmp(ptype_lower.c_str(),"tfm")==0){
        signal_container[sig-1].phase_type = gr::analog::cpm::TFM;
      }
      else if(strcmp(ptype_lower.c_str(),"gaussian")==0){
        signal_container[sig-1].phase_type = gr::analog::cpm::GAUSSIAN;
      }
      else{
        config_destroy(cfp);
        fprintf(stderr, "Unknown 'ptype' : signal %d : %s\n", sig, ptype);
        return 1;
      }
    }
    Pprintf("Signal %d, ptype: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].phase_type);

    // Checking for the number of pilot symbols (OFDM)
    strcpy(str,signumber);
    strcat(str,".Npilot");
    int pilot_count;
    if(!config_lookup_int(cfp, str, &pilot_count)){
      //Failed to find keyword
      signal_container[sig-1].pilot_count = 0;
    }
    else{
      //Found keyword
      signal_container[sig-1].pilot_count = pilot_count;
    }
    Pprintf("Signal %d, Npilot: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].pilot_count);

    // Checking for the number of taps in the pulse shape (PSK/PAM/QAM/OFDM...)
    strcpy(str,signumber);
    strcat(str,".pulse_length");
    int pulse_len;
    if(!config_lookup_int(cfp, str, &pulse_len)){
      //Failed to find keyword
      signal_container[sig-1].pulse_len = 0;
    }
    else{
      //Found keyword
      signal_container[sig-1].pulse_len = pulse_len;
    }
    Pprintf("Signal %d, pulse_length: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].pulse_len);

    // Checking pilot location vector (OFDM)
    const config_setting_t *ref;
    strcpy(str,signumber);
    strcat(str,".Ploc");
    ref = config_lookup(cfp, str);
    if(ref){
      // Found Keyword
      std::vector<int> Ploc(signal_container[sig-1].pilot_count, 0);
      std::vector<size_t> ploc(signal_container[sig-1].pilot_count, 0);
      for(int idx = 0; idx < signal_container[sig-1].pilot_count; idx++){
        Ploc[idx] = config_setting_get_int_elem(ref,idx);
        ploc[idx] = size_t(Ploc[idx]);
      }
      signal_container[sig-1].pilot_locations = std::vector<size_t>(ploc.begin(),ploc.end());
    }
    Pprintf("Signal %d, Ploc size: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].pilot_locations.size());

    // Checking pulse shaping vector (PSK/PAM/QAM...) / interpolation taps (OFDM)
    strcpy(str,signumber);
    strcat(str,".pulse_shape");
    ref = config_lookup(cfp, str);
    if(ref){
      // Found Keyword
      std::vector<double> pulse_shape(signal_container[sig-1].pulse_len, 1.);
      std::vector<float> ps(signal_container[sig-1].pulse_len, 1.);
      for(int idx = 0; idx < signal_container[sig-1].pulse_len; idx++){
        pulse_shape[idx] = config_setting_get_float_elem(ref,idx);
        ps[idx] = float(pulse_shape[idx]);
      }
      signal_container[sig-1].pulse_shape = std::vector<float>(ps.begin(),ps.end());
    }
    else{
      signal_container[sig-1].pulse_shape = std::vector<float>(1,1.);
      signal_container[sig-1].pulse_len = 1;
    }
    Pprintf("Signal %d, pulse_shape size: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].pulse_shape.size());

    // Checking for the modulation order
    strcpy(str,signumber);
    strcat(str,".order");
    if(!config_lookup_int(cfp, str, &signal_container[sig-1].order)){
      //Failed to find keyword
      signal_container[sig-1].order = 2;
    }
    Pprintf("Signal %d, order: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].order);

    // Checking for the number of subcarriers (OFDM)
    strcpy(str,signumber);
    strcat(str,".Nsc");
    int active_carriers;
    if(!config_lookup_int(cfp, str, &active_carriers)){
      //Failed to find keyword
      signal_container[sig-1].active_carriers = 0;
    }
    else{
      //Found keyword
      signal_container[sig-1].active_carriers = active_carriers;
    }
    Pprintf("Signal %d, Nsc: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].active_carriers);

    // Checking for the fft size (OFDM)
    strcpy(str,signumber);
    strcat(str,".NFFT");
    int fftsize;
    if(!config_lookup_int(cfp, str, &fftsize)){
      //Failed to find keyword
      signal_container[sig-1].fftsize = 256;
    }
    else{
      //Found keyword
      signal_container[sig-1].fftsize = fftsize;
    }
    Pprintf("Signal %d, NFFT: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].fftsize);

    // Checking for the cyclic prefix length (OFDM)
    strcpy(str,signumber);
    strcat(str,".CPL");
    int cp_len;
    if(!config_lookup_int(cfp, str, &cp_len)){
      //Failed to find keyword
      signal_container[sig-1].cp_len = 0;
    }
    else{
      //Found keyword
      signal_container[sig-1].cp_len = cp_len;
    }
    Pprintf("Signal %d, CPL: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].cp_len);

    // Checking for the number of OFDM symbols in an OFDM frame (OFDM)
    strcpy(str,signumber);
    strcat(str,".SPF");
    int syms_per_frame;
    if(!config_lookup_int(cfp, str, &syms_per_frame)){
      //Failed to find keyword
      signal_container[sig-1].syms_per_frame = 1;
    }
    else{
      //Found keyword
      signal_container[sig-1].syms_per_frame = syms_per_frame;
    }
    Pprintf("Signal %d, SPF: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].syms_per_frame);

    // Checking for the number of OFDM sample overlapping in an OFDM symbols (OFDM)
    strcpy(str,signumber);
    strcat(str,".samp_overlap");
    int samp_overlap;
    if(!config_lookup_int(cfp, str, &samp_overlap)){
      //Failed to find keyword
      signal_container[sig-1].samp_overlap = 0;
    }
    else{
      //Found keyword
      signal_container[sig-1].samp_overlap = samp_overlap;
    }
    Pprintf("Signal %d, samp_overlap: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].samp_overlap);

    // Checking taper vector (OFDM...)
    strcpy(str,signumber);
    strcat(str,".taper");
    ref = config_lookup(cfp, str);
    if(ref){
      // Found Keyword
      std::vector<double> tapering(signal_container[sig-1].samp_overlap, 1.);
      std::vector<float> taper(signal_container[sig-1].samp_overlap, 1.);
      for(int idx = 0; idx < signal_container[sig-1].samp_overlap; idx++){
        tapering[idx] = config_setting_get_float_elem(ref,idx);
        taper[idx] = float(tapering[idx]);
      }
      signal_container[sig-1].taper = std::vector<float>(taper.begin(),taper.end());
    }
    else{
      signal_container[sig-1].taper = std::vector<float>(0);
    }
    Pprintf("Signal %d, taper size: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].taper.size());

    // Checking for the number of samples per symbol (PSK/PAM/QAM...) / integer interpolation (OFDM)
    strcpy(str,signumber);
    strcat(str,".SPS");
    if(!config_lookup_int(cfp, str, &signal_container[sig-1].sps)){
      //Failed to find keyword
      strcpy(str,signumber);
      strcat(str,".sps");//I keep forgetting this one
      if(!config_lookup_int(cfp, str, &signal_container[sig-1].sps)){
        //Failed to find both keywords
        signal_container[sig-1].sps = 1;
      }
    }
    Pprintf("Signal %d, SPS: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].sps);

    // Checking for the number characters per word (CWMORSE)
    strcpy(str,signumber);
    strcat(str,".char_per_word");
    if(!config_lookup_int(cfp, str, &signal_container[sig-1].char_per_word)){
      //Failed to find keyword
      signal_container[sig-1].char_per_word = 4;
    }
    Pprintf("Signal %d, char_per_word: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].char_per_word);

    // Checking for the number overlapping symbols with phase shaping (CPM...)
    strcpy(str,signumber);
    strcat(str,".L");
    int L;
    if(!config_lookup_int(cfp, str, &L)){
      //Failed to find keyword
      signal_container[sig-1].L = 4;
    }
    else{
      //Found keyword
      signal_container[sig-1].L = L;
    }
    Pprintf("Signal %d, L: %u", size_t((sig-1)%8+1), sig, signal_container[sig-1].L);

    // Checking for the sample rate for generation of the signal
    strcpy(str,signumber);
    strcat(str,".Fs");
    double fs;
    if(!config_lookup_float(cfp, str, &fs)){
      //Failed to find keyword
      signal_container[sig-1].fs = 1.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].fs = fs;
    }
    Pprintf("Signal %d, Fs: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].fs);

    // Checking for the center frequency of the signal
    strcpy(str,signumber);
    strcat(str,".CenterFreq");
    double fc;
    if(!config_lookup_float(cfp, str, &fc)){
      //Failed to find keyword
      signal_container[sig-1].fc = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].fc = fc;
    }
    Pprintf("Signal %d, CenterFreq: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].fc);

    // Checking for the gain (dB) of the signal
    strcpy(str,signumber);
    strcat(str,".Gain");
    double gain;
    if(!config_lookup_float(cfp, str, &gain)){
      //Failed to find keyword
      signal_container[sig-1].gain = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].gain = gain;
    }
    Pprintf("Signal %d, Gain: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].gain);

    // Checking for the offset of the signal (phase ..)
    strcpy(str,signumber);
    strcat(str,".offset");
    double offset;
    if(!config_lookup_float(cfp, str, &offset)){
      //Failed to find keyword
      signal_container[sig-1].offset = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].offset = offset;
    }
    Pprintf("Signal %d, offset: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].offset);

    // Checking for the digital backoff (dB) (OFDM)
    strcpy(str,signumber);
    strcat(str,".backoff");
    double backoff;
    if(!config_lookup_float(cfp, str, &backoff)){
      //Failed to find keyword
      signal_container[sig-1].backoff = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].backoff = backoff;
    }
    Pprintf("Signal %d, backoff: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].backoff);

    // Checking for the analog component count (Analog...)
    strcpy(str,signumber);
    strcat(str,".components");
    int components;
    if(!config_lookup_int(cfp, str, &components)){
      //Failed to find keyword
      signal_container[sig-1].components = 0;
    }
    else{
      //Found keyword
      signal_container[sig-1].components = components;
    }
    Pprintf("Signal %d, components: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].components);

    // Checking mu vector (Analog...)
    strcpy(str,signumber);
    strcat(str,".mu");
    ref = config_lookup(cfp, str);
    if(ref){
      // Found Keyword
      std::vector<double> mus(signal_container[sig-1].components, 0.);
      std::vector<float> mu(signal_container[sig-1].components, 0.);
      for(int idx = 0; idx < signal_container[sig-1].components; idx++){
        mus[idx] = config_setting_get_float_elem(ref,idx);
        mu[idx] = float(mus[idx]);
      }
      signal_container[sig-1].mu = std::vector<float>(mu.begin(),mu.end());
    }
    else{
      signal_container[sig-1].mu = std::vector<float>(0);
    }
    Pprintf("Signal %d, mu size: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].mu.size());

    // Checking sigma vector (Analog...)
    strcpy(str,signumber);
    strcat(str,".sigma");
    ref = config_lookup(cfp, str);
    if(ref){
      // Found Keyword
      std::vector<double> sigs(signal_container[sig-1].components, 0.);
      std::vector<float> sigm(signal_container[sig-1].components, 0.);
      for(int idx = 0; idx < signal_container[sig-1].components; idx++){
        sigs[idx] = config_setting_get_float_elem(ref,idx);
        sigm[idx] = float(sigs[idx]);
      }
      signal_container[sig-1].sigma = std::vector<float>(sigm.begin(),sigm.end());
    }
    else{
      signal_container[sig-1].sigma = std::vector<float>(0);
    }
    Pprintf("Signal %d, sigma size: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].sigma.size());

    // Checking weight vector (Analog...)
    strcpy(str,signumber);
    strcat(str,".weight");
    ref = config_lookup(cfp, str);
    if(ref){
      // Found Keyword
      std::vector<double> weights(signal_container[sig-1].components, 0.);
      std::vector<float> weight(signal_container[sig-1].components, 0.);
      for(int idx = 0; idx < signal_container[sig-1].components; idx++){
        weights[idx] = config_setting_get_float_elem(ref,idx);
        weight[idx] = float(weights[idx]);
      }
      signal_container[sig-1].weight = std::vector<float>(weight.begin(),weight.end());
    }
    else{
      signal_container[sig-1].weight = std::vector<float>(0);
    }
    Pprintf("Signal %d, weight size: %lu", size_t((sig-1)%8+1), sig, signal_container[sig-1].weight.size());

    // Checking for the modulation index
    strcpy(str,signumber);
    strcat(str,".mod_idx");
    double mod_idx;
    if(!config_lookup_float(cfp, str, &mod_idx)){
      //Failed to find keyword
      signal_container[sig-1].mod_idx = 0.5;
    }
    else{
      //Found keyword
      signal_container[sig-1].mod_idx = mod_idx;
    }
    Pprintf("Signal %d, mod_idx: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].mod_idx);

    // Checking for the number of words per minute (CWMORSE)
    strcpy(str,signumber);
    strcat(str,".words_per_min");
    double words_per_minute;
    if(!config_lookup_float(cfp, str, &words_per_minute)){
      //Failed to find keyword
      signal_container[sig-1].words_per_minute = 120;
    }
    else{
      //Found keyword
      signal_container[sig-1].words_per_minute = words_per_minute;
    }
    Pprintf("Signal %d, words_per_minute: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].words_per_minute);

    // Checking for the pulse shaping parameter
    strcpy(str,signumber);
    strcat(str,".beta");
    double beta;
    if(!config_lookup_float(cfp, str, &beta)){
      //Failed to find keyword
      signal_container[sig-1].beta = 0.35;
    }
    else{
      //Found keyword
      signal_container[sig-1].beta = beta;
    }
    Pprintf("Signal %d, beta: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].beta);

    // Checking for the base word (CWMORSE)
    strcpy(str,signumber);
    strcat(str,".base_word");
    int base_word;
    if(!config_lookup_bool(cfp, str, &base_word)){
      //Failed to find keyword
      signal_container[sig-1].base_word = false;
    }
    else{
      //Found keyword
      signal_container[sig-1].base_word = base_word;
    }
    Pprintf("Signal %d, base_word: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].base_word);

    // Checking for the Pilots Per Frame flag (OFDM)
    strcpy(str,signumber);
    strcat(str,".ppf");
    int ppf;
    if(!config_lookup_bool(cfp, str, &ppf)){
      //Failed to find keyword
      signal_container[sig-1].pilot_per_frame = false;
    }
    else{
      //Found keyword
      signal_container[sig-1].pilot_per_frame = ppf;
    }
    Pprintf("Signal %d, ppf: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].pilot_per_frame);

    // Checking for the add sync symbols flag (OFDM)
    strcpy(str,signumber);
    strcat(str,".addSync");
    int addSync;
    if(!config_lookup_bool(cfp, str, &addSync)){
      //Failed to find keyword
      signal_container[sig-1].add_sync = false;
    }
    else{
      //Found keyword
      signal_container[sig-1].add_sync = addSync;
    }
    Pprintf("Signal %d, addSync: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].add_sync);

    // Checking for the opg flag
    strcpy(str,signumber);
    strcat(str,".ops_gate");
    int opsg;
    if(!config_lookup_bool(cfp, str, &opsg)){
      //Failed to find keyword
      signal_container[sig-1].ops_gate = false;
    }
    else{
      //Found keyword
      signal_container[sig-1].ops_gate = opsg;
    }
    Pprintf("Signal %d, ops_gate: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].ops_gate);

    // Checking for the per_gate flag
    strcpy(str,signumber);
    strcat(str,".per_gate");
    int perg;
    if(!config_lookup_bool(cfp, str, &perg)){
      //Failed to find keyword
      signal_container[sig-1].per_gate = false;
    }
    else{
      //Found keyword
      signal_container[sig-1].per_gate = perg;
    }
    Pprintf("Signal %d, per_gate: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].per_gate);

    // Checking for the rnd_gate flag
    strcpy(str,signumber);
    strcat(str,".rnd_gate");
    int rndg;
    if(!config_lookup_bool(cfp, str, &rndg)){
      //Failed to find keyword
      signal_container[sig-1].rnd_gate = false;
    }
    else{
      //Found keyword
      signal_container[sig-1].rnd_gate = rndg;
    }
    Pprintf("Signal %d, rnd_gate: %d", size_t((sig-1)%8+1), sig, signal_container[sig-1].rnd_gate);

    // Checking for the per gate on duration
    strcpy(str,signumber);
    strcat(str,".per_gate_on");
    double pgon;
    if(!config_lookup_float(cfp, str, &pgon)){
      //Failed to find keyword
      signal_container[sig-1].per_gate_on = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].per_gate_on = pgon;
    }
    Pprintf("Signal %d, per_gate_on: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].per_gate_on);

    // Checking for the per gate off duration
    strcpy(str,signumber);
    strcat(str,".per_gate_off");
    double pgoff;
    if(!config_lookup_float(cfp, str, &pgoff)){
      //Failed to find keyword
      signal_container[sig-1].per_gate_off = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].per_gate_off = pgoff;
    }
    Pprintf("Signal %d, per_gate_off: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].per_gate_off);

    // Checking for the per gate periodic offset
    strcpy(str,signumber);
    strcat(str,".per_gate_offset");
    double pgoffset;
    if(!config_lookup_float(cfp, str, &pgoffset)){
      //Failed to find keyword
      signal_container[sig-1].per_gate_offset = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].per_gate_offset = pgoffset;
    }
    Pprintf("Signal %d, per_gate_offset: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].per_gate_offset);

    // Checking for the one pass gate on duration
    strcpy(str,signumber);
    strcat(str,".ops_gate_on");
    double ogon;
    if(!config_lookup_float(cfp, str, &ogon)){
      //Failed to find keyword
      signal_container[sig-1].ops_gate_on = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].ops_gate_on = ogon;
    }
    Pprintf("Signal %d, ops_gate_on: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].ops_gate_on);

    // Checking for the one pas gate off duration
    strcpy(str,signumber);
    strcat(str,".ops_gate_off");
    double ogoff;
    if(!config_lookup_float(cfp, str, &ogoff)){
      //Failed to find keyword
      signal_container[sig-1].ops_gate_off = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].ops_gate_off = ogoff;
    }
    Pprintf("Signal %d, ops_gate_off: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].ops_gate_off);

    // Checking for the random gate on duration min
    strcpy(str,signumber);
    strcat(str,".rnd_gate_on_min");
    double rgonmin;
    if(!config_lookup_float(cfp, str, &rgonmin)){
      //Failed to find keyword
      signal_container[sig-1].rnd_gate_on_min = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].rnd_gate_on_min = rgonmin;
    }
    Pprintf("Signal %d, rnd_gate_on_min: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].rnd_gate_on_min);

    // Checking for the random gate on duration max
    strcpy(str,signumber);
    strcat(str,".rnd_gate_on_max");
    double rgonmax;
    if(!config_lookup_float(cfp, str, &rgonmax)){
      //Failed to find keyword
      signal_container[sig-1].rnd_gate_on_max = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].rnd_gate_on_max = rgonmax;
    }
    Pprintf("Signal %d, rnd_gate_on_max: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].rnd_gate_on_max);

    // Checking for the random gate off duration min
    strcpy(str,signumber);
    strcat(str,".rnd_gate_off_min");
    double rgoffmin;
    if(!config_lookup_float(cfp, str, &rgoffmin)){
      //Failed to find keyword
      signal_container[sig-1].rnd_gate_off_min = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].rnd_gate_off_min = rgoffmin;
    }
    Pprintf("Signal %d, rnd_gate_off_min: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].rnd_gate_off_min);

    // Checking for the random gate off duration max
    strcpy(str,signumber);
    strcat(str,".rnd_gate_off_max");
    double rgoffmax;
    if(!config_lookup_float(cfp, str, &rgoffmax)){
      //Failed to find keyword
      signal_container[sig-1].rnd_gate_off_max = 0.0;
    }
    else{
      //Found keyword
      signal_container[sig-1].rnd_gate_off_max = rgoffmax;
    }
    Pprintf("Signal %d, rnd_gate_off_max: %1.3e", size_t((sig-1)%8+1), sig, signal_container[sig-1].rnd_gate_off_max);


    //Value tweaks
    if(signal_container[sig-1].type == gr::signal_exciter::CWMORSE){
      if(signal_container[sig-1].base_word)
        signal_container[sig-1].fs = signal_container[sig-1].words_per_minute;
      else
        signal_container[sig-1].fs = signal_container[sig-1].words_per_minute/1.2;
    }

  }

  config_destroy(cfp);
  return 0;
}

int GenerateFlowGraph(system_var system_config, gr::signal_exciter::sig_params* signal_list, gr::top_block_sptr& tb, size_t& byteCount)
{
  size_t SynthSize = SYNTH_SIZE;
  size_t SignalCount = system_config.NumSignals;
  size_t SampleCount = size_t(ceil(system_config.SysBW * system_config.RunTime));
  byteCount = SampleCount*8;
  double NoiseAmplitude = sqrt( 50. * pow(10., ( system_config.NoiseFloor + 10.*log10( system_config.SysBW ) - 30 ) / 10. ) );

  Cprintf("This config file should generate %lu signals, with %lu samples.\nA total of %lu bytes will be written.", size_t(190), SignalCount, SampleCount, byteCount);

  Cprintf("Init gnuradio blocks%s", size_t(190), ".");

  // The signal source blocks
  gr::signal_exciter::random_signal::sptr           sig_sources[SignalCount];

  // The gain multipliers
  gr::blocks::multiply_const_cc::sptr               multipliers[SignalCount];

  // Interpolation taps
  std::vector< std::vector<float> >                 arb_resamp_taps;
  std::vector< std::vector< std::vector<float> > >  int_resamp_taps;

  // Channelizer Taps
  std::vector< std::vector<float> >                 channelizer_taps;

  // Synthesizer Taps
  std::vector<float>                                synthesizer_taps;
  std::vector<float>                                synth_dec_taps;

  // Arbitrary Resamplers
  gr::filter::pfb_arb_resampler_ccf::sptr           arb_resamplers[SignalCount];
  gr::blocks::skiphead::sptr                        arb_skip[SignalCount];
  // Integer Resamplers - > Allows for max rate increase of 100 per resamp
  std::vector< gr::filter::interp_fir_filter_ccf::sptr >  int_resamplers[SignalCount];
  std::vector< gr::blocks::skiphead::sptr >               int_skip[SignalCount];

  // Frequency Translation
  gr::analog::sig_source_c::sptr                    xlators[SignalCount];
  gr::blocks::multiply_cc::sptr                     xlate_m[SignalCount];

  // Signal Gates
  gr::signal_exciter::one_pass_gate::sptr           one_pass_gates[SignalCount];
  gr::signal_exciter::periodic_gate::sptr           periodic_gates[SignalCount];
  gr::signal_exciter::random_gate::sptr             random_gates[SignalCount];

  // Channelizer
  gr::blocks::stream_to_streams::sptr               chan_streamer[SignalCount];
  gr::filter::pfb_channelizer_ccf::sptr             channelizers[SignalCount];
  std::vector< gr::blocks::skiphead::sptr >         chann_skip[SignalCount];

  // Add up the channel connections
  gr::blocks::add_cc::sptr                          adders[SynthSize-1];
  std::vector<size_t>                               connections(SynthSize-1,0);

  // Synthesizer
  gr::filter::pfb_synthesizer_ccf::sptr             synthesizer;
  gr::filter::fir_filter_ccf::sptr                  synth_decim;
  gr::blocks::skiphead::sptr                        synth_skip;
  gr::blocks::skiphead::sptr                        syndec_skip;

  // Noise Source
  gr::channels::channel_model::sptr                 noise_source;

  // File saving
  gr::blocks::head::sptr                            head_clip;
  gr::blocks::file_sink::sptr                       file_saver;

#if DEBUG_OUTPUT
  gr::blocks::file_sink::sptr                       debug_sources[SignalCount];
  gr::blocks::file_sink::sptr                       debug_gains[SignalCount];
  gr::blocks::file_sink::sptr                       debug_arb[SignalCount];
  gr::blocks::file_sink::sptr                       debug_int[SignalCount];
  gr::blocks::file_sink::sptr                       debug_xlate[SignalCount];
  gr::blocks::file_sink::sptr                       debug_opsg[SignalCount];
  gr::blocks::file_sink::sptr                       debug_perg[SignalCount];
  gr::blocks::file_sink::sptr                       debug_rndg[SignalCount];
  gr::blocks::file_sink::sptr                       debug_synth;
  char debug_head[50];
  char str[50];
#endif

  // Optimize Connection Routine
  std::vector<gr::basic_block_sptr>                 block_list[SignalCount];

  Cprintf("Building up signal chains%s", size_t(190), ".");

  // Bandwidth of a single synth channel?
  double synth_chan = double(system_config.SysBW)/double(SynthSize);

  // What are the limitations of the synthesizer? (Ignoring Wrap around channel)
  std::vector<double> synth_chan_boundaries(SynthSize,0.);
  std::vector<double> synth_chan_centers(SynthSize-1,0.);

  synth_chan_boundaries[0] = -double(system_config.SysBW)/2. + synth_chan/2.;
  for(size_t idx = 0; idx < SynthSize-1; idx++){
    //Define boundaries and centers
    synth_chan_boundaries[idx+1] = synth_chan_boundaries[idx] + synth_chan;
    synth_chan_centers[idx] = (synth_chan_boundaries[idx] + synth_chan_boundaries[idx+1])/2.;

    // Initialize the adder blocks
    adders[idx] = gr::blocks::add_cc::make(1);
  }

  Cprintf("The synthesizer is operating over a %1.3e bandwidth, with %lu channels, and each channel having a bandwidth of %1.3e.", size_t(190), system_config.SysBW, SynthSize, synth_chan);

  for(size_t sig_idx = 0; sig_idx < SignalCount; sig_idx++){
    //Initialize block list for this chain
    block_list[sig_idx] = std::vector<gr::basic_block_sptr>(0);

    // Track the signal type for special signal type tweaks
    gr::signal_exciter::sig_type_t sig_type = signal_list[sig_idx].type;

    // Space needed for this signal
    double rrate = double(signal_list[sig_idx].fs)/synth_chan;
    size_t nbins = int(ceil(rrate));
    size_t crit_bins = nbins + (nbins==SynthSize ? 0 : (nbins%2 ? 0 : 1));//This is the number of bins needed to do as requested
    nbins = crit_bins;

    Cprintf("Signal %lu fits into %1.3e bins, and has a critical bin size of %lu.", sig_idx%8+1, sig_idx+1, rrate, crit_bins);

    bool good_config = false; //Run until a good configuration is found, or unobtainable
    double xlate_freq = 0.;   //How much translation is needed?
    size_t synth_chan_center; //Which synth channel is it centered in?
    size_t chann_chan_center; //Which chann channel is it centered in? (due to channel map)
    size_t synth_start_connection; //Which synth channel should it start connecting at?
    double target_fs;
    while(!good_config){
      Cprintf("Signal %lu aiming for %lu channels.", sig_idx%8+1, sig_idx+1, nbins);
      if(nbins>SynthSize){
        fprintf(stderr, "Signal %lu requires a channelizer larger than synthesizer. This is not currently allowed.\n", sig_idx+1);
        return 1;
      }

      xlate_freq = signal_list[sig_idx].fc - system_config.CenterFreq;
      if((xlate_freq <= synth_chan_boundaries[0]) || (xlate_freq >= synth_chan_boundaries[SynthSize-1])){
        fprintf(stderr, "Signal %lu is centered outside the system's bandwidth. This is not curently allowed.\n", sig_idx+1);
        return 1;
      }

      Cprintf("Signal %lu needs to translate by %1.3e.", sig_idx%8+1, sig_idx+1, xlate_freq);

      synth_chan_center = 0;
      while((xlate_freq > synth_chan_boundaries[synth_chan_center]) && (synth_chan_center < SynthSize-1)){
        if(xlate_freq > synth_chan_boundaries[synth_chan_center+1]){
          // Centerd in at least the next channel
          synth_chan_center++;
        }
        else{
          // Definately in this channel
          break;
        }
      }

      Cprintf("Signal %lu is centered in channel %lu with a channel center %1.3e.", sig_idx%8+1, sig_idx+1, synth_chan_center, synth_chan_centers[synth_chan_center]);

      xlate_freq = xlate_freq - synth_chan_centers[synth_chan_center];

      Cprintf("Signal %lu is requires a relative offset of %1.3e for that synth channel.", sig_idx%8+1, sig_idx+1, xlate_freq);

      double upper_limit = xlate_freq + double(signal_list[sig_idx].fs)/2. + synth_chan_centers[synth_chan_center];
      double lower_limit = xlate_freq - double(signal_list[sig_idx].fs)/2. + synth_chan_centers[synth_chan_center];
      if(nbins == 1){
        //Check if the xlate causes the bin size to increase
        if((lower_limit < synth_chan_boundaries[synth_chan_center]) || (upper_limit > synth_chan_boundaries[synth_chan_center+1])){
          // The translation incurred channel requried count
          nbins = 2;
          Cprintf("\tNeed to increase the channel count.\n\t\t[Lower limit(%1.3e), Lower bound(%1.3e)]\n\t\t[Upper limit(%1.3e), Upper bound(%1.3e)]", sig_idx%8+1, lower_limit, synth_chan_boundaries[synth_chan_center], upper_limit, synth_chan_boundaries[synth_chan_center+1]);
        }
        else{
          // Signal is fully within a single synth channel
          target_fs = 2*synth_chan;//2x oversampled synthesizer
          chann_chan_center = 0;
          synth_start_connection = synth_chan_center;
          good_config = true;
        }
      }
      else if(nbins==SynthSize){
        //Signal fits exactly into desired spectrum... special case
        if((lower_limit < -system_config.SysBW/2.) || (upper_limit > system_config.SysBW/2.)){
          nbins = nbins*2;
          Cprintf("\tNeed to increase the channel count.\n\t\t[Lower limit(%1.3e), Lower bound(%1.3e)]\n\t\t[Upper limit(%1.3e), Upper bound(%1.3e)]", sig_idx%8+1, lower_limit, synth_chan_boundaries[synth_chan_center], upper_limit, synth_chan_boundaries[synth_chan_center+1]);
        }
        else{
          target_fs = system_config.SysBW;
          chann_chan_center = nbins/2;
          synth_start_connection = 0;
          crit_bins = SynthSize-1;
          good_config = true;
        }
      }
      else{
        //Needs at least 2 bins, check and make sure all bins fit.
        //Hard coding to power of 2 channelizer
        //////////////////////////////////////////////////////////
        int chan_counter(1);
        float chan_temp(nbins);
        while(chan_temp > 1.){
          chan_temp = chan_temp/2.;
          chan_counter *= 2;
        }
        nbins = chan_counter;
        //////////////////////////////////////////////////////////
        //There are an even number of bins, but the wrap around bin will be tossed
        int edge_room = (nbins-1)/2;//room upper and lower
        if((int(synth_chan_center) - edge_room < 0) || (int(synth_chan_center)+edge_room > SynthSize-1)){
          fprintf(stderr, "Signal %lu has critical bandwidth outside acceptable parameters.\n",sig_idx+1);
          return 1;
        }
        if((lower_limit < synth_chan_boundaries[synth_chan_center-edge_room]) || (upper_limit > synth_chan_boundaries[synth_chan_center+1+edge_room])){
          // True bandwidth escaped the critical boundary
          nbins *= 2;
          Cprintf("\tNeed to increase the channel count.\n\t\t[Lower limit(%1.3e), Lower bound(%1.3e)]\n\t\t[Upper limit(%1.3e), Upper bound(%1.3e)]", sig_idx%8+1, lower_limit, synth_chan_boundaries[synth_chan_center-edge_room], upper_limit, synth_chan_boundaries[synth_chan_center+1+edge_room]);
        }
        else{
          // The critical bandwidth is acceptable
          target_fs = synth_chan * nbins; // channelizer/synth both at 2x
          chann_chan_center = nbins/2;// [-nbins/2, nbins/2-1]
          synth_start_connection = synth_chan_center - (crit_bins-1)/2;
          good_config = true;
        }
      }
    }

    Cprintf("Signal %lu will use a %lu channel channelizer and only connect the center %lu channels.\n\tThe starting connection will be on synthesizer channel %lu.", sig_idx%8+1, sig_idx+1, nbins, crit_bins, synth_start_connection);

    //Initialize the signal source
    sig_sources[sig_idx] = gr::signal_exciter::random_signal::make(signal_list[sig_idx], system_config.Seed++);
    block_list[sig_idx].push_back(sig_sources[sig_idx]);
#if DEBUG_OUTPUT
    sprintf(debug_head, "debug_sig%02lu_source.fc32",sig_idx+1);
    strcpy(str,debug_head);
    debug_sources[sig_idx] = gr::blocks::file_sink::make(sizeof(gr_complex), str, false);
    tb->connect(sig_sources[sig_idx], 0, debug_sources[sig_idx], 0);
#endif

    double upsamp = target_fs/signal_list[sig_idx].fs;

    Cprintf("Signal %lu will require %1.3e upsampling.", sig_idx%8+1, sig_idx+1, upsamp);

    int int_filts = 0;
    float arb_rem = upsamp;
    std::vector<int> upsamp_factors(0);
    while(arb_rem >= 2.0){
      if(arb_rem > 100.0){
        upsamp_factors.push_back(100);
        arb_rem = arb_rem/100.0;
      }
      else{
        int rem = int(arb_rem);
        arb_rem = arb_rem/double(rem);
        upsamp_factors.push_back(rem);
      }
      int_filts++;
    }

    /*Cprintf("Signal %lu breaks up the resampling as: [", sig_idx%8+1, sig_idx+1);
    for( size_t idx=0; idx < int_filts; idx++ ){
      Cprintf("\t%1.3e", sig_idx%8+1, upsamp_factors[idx]);
    }
    Cprintf("\t%1.3e]", sig_idx%8+1, arb_rem);*/

    Cprintf("Signal %lu will use %d integer upsamplers and an arbitrary filter with %1.3e.", sig_idx%8+1, sig_idx+1, int_filts, arb_rem);

    // Check for necessary blocks
    bool use_gan = (signal_list[sig_idx].gain != 0.);
    bool use_arb = (arb_rem != 1.0);
    bool use_int = (int_filts > 0);
    bool use_xlt = (xlate_freq != 0.);
    bool use_per = (signal_list[sig_idx].per_gate);
    bool use_rnd = (signal_list[sig_idx].rnd_gate);
    bool use_opg = (signal_list[sig_idx].ops_gate);
    bool use_chn = (nbins > 1);

    Cprintf("Signal %lu, use gain? %d, use arb_resamp? %d, use int_resamp? %d, use xlator? %d, use one pass gate? %d, use periodic gate? %d, use random gate? %d, use channelizer? %d.", sig_idx%8+1, sig_idx+1, use_gan, use_arb, use_int, use_xlt, use_opg, use_per, use_rnd, use_chn);

    if(use_gan){
      Cprintf("Signal %lu, making gain.",sig_idx%8+1,sig_idx+1);
      double gain = pow(10., signal_list[sig_idx].gain/10.);
      multipliers[sig_idx] = gr::blocks::multiply_const_cc::make(gr_complex(gain,0.),1);
      block_list[sig_idx].push_back(multipliers[sig_idx]);
#if DEBUG_OUTPUT
      sprintf(debug_head, "debug_sig%02lu_multi.fc32",sig_idx+1);
      strcpy(str,debug_head);
      debug_gains[sig_idx] = gr::blocks::file_sink::make(sizeof(gr_complex), str, false);
      tb->connect(multipliers[sig_idx], 0, debug_gains[sig_idx], 0);
#endif
    }
    if(use_arb){
      Cprintf("Signal %lu, making arb.",sig_idx%8+1,sig_idx+1);
      arb_resamp_taps.push_back(gen_arb_interp_taps(arb_rem));
      arb_resamplers[sig_idx] = gr::filter::pfb_arb_resampler_ccf::make(arb_rem, arb_resamp_taps[sig_idx], 32);
      block_list[sig_idx].push_back(arb_resamplers[sig_idx]);
      size_t artl = (arb_resamp_taps[sig_idx].size()-1)/2;
      arb_skip[sig_idx] = gr::blocks::skiphead::make(sizeof(gr_complex), artl);
      block_list[sig_idx].push_back(arb_skip[sig_idx]);
#if DEBUG_OUTPUT
      sprintf(debug_head, "debug_sig%02lu_arb.fc32",sig_idx+1);
      strcpy(str,debug_head);
      debug_arb[sig_idx] = gr::blocks::file_sink::make(sizeof(gr_complex), str, false);
      tb->connect(arb_resamplers[sig_idx], 0, debug_arb[sig_idx], 0);
#endif
    }
    else{
      arb_resamp_taps.push_back(std::vector<float>(0));
    }
    if(use_int){
      Cprintf("Signal %lu, making int.",sig_idx%8+1,sig_idx+1);
      int_resamp_taps.push_back(std::vector< std::vector<float> >(0));
      for(size_t filtidx = 0; filtidx < int_filts; filtidx++){
        int_resamp_taps[sig_idx].push_back(gen_int_interp_taps(upsamp_factors[filtidx]));
        int_resamplers[sig_idx].push_back(gr::filter::interp_fir_filter_ccf::make(upsamp_factors[filtidx], int_resamp_taps[sig_idx][filtidx]));
        block_list[sig_idx].push_back(int_resamplers[sig_idx][filtidx]);
        size_t irtl = (int_resamp_taps[sig_idx][filtidx].size()-1)/2;
        int_skip[sig_idx].push_back(gr::blocks::skiphead::make(sizeof(gr_complex), irtl));
        block_list[sig_idx].push_back(int_skip[sig_idx][filtidx]);
      }
#if DEBUG_OUTPUT
      sprintf(debug_head, "debug_sig%02lu_int.fc32",sig_idx+1);
      strcpy(str,debug_head);
      debug_int[sig_idx] = gr::blocks::file_sink::make(sizeof(gr_complex), str, false);
      tb->connect(int_resamplers[sig_idx][int_filts-1], 0, debug_int[sig_idx], 0);
#endif
    }
    else{
      int_resamp_taps.push_back(std::vector< std::vector<float> >(0));
    }
    if(use_xlt){
      Cprintf("Signal %lu, making xlate.",sig_idx%8+1,sig_idx+1);
      xlators[sig_idx] = gr::analog::sig_source_c::make(target_fs, gr::analog::GR_COS_WAVE, xlate_freq, 1., 0.);
      xlate_m[sig_idx] = gr::blocks::multiply_cc::make(1);
      tb->connect(xlators[sig_idx], 0, xlate_m[sig_idx], 1);
      block_list[sig_idx].push_back(xlate_m[sig_idx]);
#if DEBUG_OUTPUT
      sprintf(debug_head, "debug_sig%02lu_xlate.fc32",sig_idx+1);
      strcpy(str,debug_head);
      debug_xlate[sig_idx] = gr::blocks::file_sink::make(sizeof(gr_complex), str, false);
      tb->connect(xlate_m[sig_idx], 0, debug_xlate[sig_idx], 0);
#endif
    }
    if(use_per){
      Cprintf("Signal %lu, making per.",sig_idx%8+1,sig_idx+1);
      periodic_gates[sig_idx] = gr::signal_exciter::periodic_gate::make(target_fs, signal_list[sig_idx].per_gate_off, signal_list[sig_idx].per_gate_on, signal_list[sig_idx].per_gate_offset);
      block_list[sig_idx].push_back(periodic_gates[sig_idx]);
#if DEBUG_OUTPUT
      sprintf(debug_head, "debug_sig%02lu_perg.fc32",sig_idx+1);
      strcpy(str,debug_head);
      debug_perg[sig_idx] = gr::blocks::file_sink::make(sizeof(gr_complex), str, false);
      tb->connect(periodic_gates[sig_idx], 0, debug_perg[sig_idx], 0);
#endif
    }
    if(use_rnd){
      Cprintf("Signal %lu, making rnd.",sig_idx%8+1,sig_idx+1);
      random_gates[sig_idx] = gr::signal_exciter::random_gate::make(target_fs, signal_list[sig_idx].rnd_gate_off_min, signal_list[sig_idx].rnd_gate_off_max, signal_list[sig_idx].rnd_gate_on_min, signal_list[sig_idx].rnd_gate_on_max,system_config.Seed++);
      block_list[sig_idx].push_back(random_gates[sig_idx]);
#if DEBUG_OUTPUT
      sprintf(debug_head, "debug_sig%02lu_rndg.fc32",sig_idx+1);
      strcpy(str,debug_head);
      debug_rndg[sig_idx] = gr::blocks::file_sink::make(sizeof(gr_complex), str, false);
      tb->connect(random_gates[sig_idx], 0, debug_rndg[sig_idx], 0);
#endif
    }
    if(use_opg){
      Cprintf("Signal %lu, making opg.",sig_idx%8+1,sig_idx+1);
      one_pass_gates[sig_idx] = gr::signal_exciter::one_pass_gate::make(target_fs, signal_list[sig_idx].ops_gate_off, signal_list[sig_idx].ops_gate_on, false);
      block_list[sig_idx].push_back(one_pass_gates[sig_idx]);
#if DEBUG_OUTPUT
      sprintf(debug_head, "debug_sig%02lu_opsg.fc32",sig_idx+1);
      strcpy(str,debug_head);
      debug_opsg[sig_idx] = gr::blocks::file_sink::make(sizeof(gr_complex), str, false);
      tb->connect(one_pass_gates[sig_idx], 0, debug_opsg[sig_idx], 0);
#endif
    }
    if(use_chn){
      Cprintf("Signal %lu, making chn.",sig_idx%8+1,sig_idx+1);
      channelizer_taps.push_back(gen_channl_taps(nbins));
      chan_streamer[sig_idx] = gr::blocks::stream_to_streams::make(sizeof(gr_complex), nbins);
      channelizers[sig_idx] = gr::filter::pfb_channelizer_ccf::make(nbins, channelizer_taps[sig_idx], 2.0);
      channelizers[sig_idx]->set_tag_propagation_policy(gr::block::TPP_DONT);
      std::vector< std::vector<float> > taps;
      taps = channelizers[sig_idx]->taps();
      size_t ctl = (taps[0].size()-1)/2;
      std::vector<int> channl_map(crit_bins);
      for(int mapper = 0; mapper < crit_bins; mapper++){
        channl_map[mapper] = (mapper - (crit_bins-1)/2)%nbins;
        chann_skip[sig_idx].push_back(gr::blocks::skiphead::make(sizeof(gr_complex), ctl));
      }
      channelizers[sig_idx]->set_channel_map(channl_map);
      for(size_t cidx = 0; cidx < nbins; cidx++){
        tb->connect(chan_streamer[sig_idx], cidx, channelizers[sig_idx], cidx);
      }
      for(size_t cidx = 0; cidx < crit_bins; cidx++){
        tb->connect(channelizers[sig_idx], cidx, chann_skip[sig_idx][cidx], 0);
        tb->connect(chann_skip[sig_idx][cidx], 0, adders[synth_start_connection+cidx], connections[synth_start_connection+cidx]++);
      }
      block_list[sig_idx].push_back(chan_streamer[sig_idx]);
    }
    else{
      channelizer_taps.push_back(std::vector<float>(0));
      tb->connect( block_list[sig_idx][block_list[sig_idx].size()-1], 0, adders[synth_start_connection], connections[synth_start_connection]++ );
    }
    Cprintf("Signal %lu, all made.",sig_idx%8+1,sig_idx+1);

    for(size_t bidx = 0; bidx < block_list[sig_idx].size()-1; bidx++){
      tb->connect( block_list[sig_idx][bidx], 0, block_list[sig_idx][bidx+1], 0 );
    }
    Cprintf("Signal %lu is being connected from source through adders.", sig_idx%8+1, sig_idx+1);
  }

  Cprintf("All %lu signals have been connected from source through adders. Now connecting to synthesier and remainders.", size_t(190), SignalCount);

  size_t connections_needed = 0;
  std::vector<size_t> connection_map;
  std::vector<int> synth_mapping;
  for(size_t sidx = 0; sidx < SynthSize-1; sidx++){
    if(connections[sidx]>0){
      connections_needed++;
      connection_map.push_back(sidx);
      synth_mapping.push_back( (sidx+1+(SynthSize/2))%SynthSize );
      if(synth_mapping[synth_mapping.size()-1] >= SynthSize/2){
        synth_mapping[synth_mapping.size()-1] += SynthSize;
      }
    }
  }

  Cprintf("The synthesizer needs %lu connections in total.", size_t(190), connections_needed);

  synthesizer_taps = gen_synthz_taps(SynthSize);
  synthesizer = gr::filter::pfb_synthesizer_ccf::make(SynthSize, synthesizer_taps, true);
  synthesizer->set_tag_propagation_policy(gr::block::TPP_DONT);

    

  synthesizer->set_channel_map(synth_mapping);

  Cprintf("The synthesizer has been created (%lu, %lu, %lu).", size_t(190), SynthSize, connections_needed, synthesizer_taps.size());

  for(size_t sidx = 0; sidx < connections_needed; sidx++){
    tb->connect( adders[connection_map[sidx]], 0, synthesizer, sidx );
  }
#if DEBUG_OUTPUT
  sprintf(debug_head, "debug_synth.fc32");
  strcpy(str,debug_head);
  debug_synth = gr::blocks::file_sink::make(sizeof(gr_complex), str, false);
  tb->connect(synthesizer, 0, debug_synth, 0);
#endif

  Cprintf("The synthesizer has been connected to %lu inputs.", size_t(190), connections_needed);

  size_t stl = (synthesizer_taps.size()-1)/2;
  synth_skip = gr::blocks::skiphead::make(sizeof(gr_complex), stl);
  tb->connect( synthesizer, 0, synth_skip, 0);

  synth_dec_taps = gr::filter::firdes::low_pass_2(2,2,(double(SynthSize)-1.)/double(SynthSize), 1./double(SynthSize), 100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
  synth_decim = gr::filter::fir_filter_ccf::make(2, synth_dec_taps);
  size_t sdtl = (synth_dec_taps.size()-1)/2;
  syndec_skip = gr::blocks::skiphead::make(sizeof(gr_complex), sdtl);

  tb->connect( synth_skip, 0, synth_decim, 0 );
  tb->connect( synth_decim, 0, syndec_skip, 0 );

  head_clip = gr::blocks::head::make(sizeof(gr_complex), SampleCount);
  file_saver = gr::blocks::file_sink::make(sizeof(gr_complex), system_config.outputFile.c_str());

  if(system_config.EnableNoise){
    noise_source = gr::channels::channel_model::make(NoiseAmplitude,0.,1.,std::vector<gr_complex>(1,gr_complex(1.,0.)));
    tb->connect(syndec_skip, 0, noise_source, 0);
    tb->connect(noise_source, 0, head_clip, 0);
  }
  else{
    tb->connect(syndec_skip, 0, head_clip, 0);
  }

  tb->connect(head_clip, 0, file_saver, 0);

  Cprintf("The TB is fully connected%s", size_t(190), ".");

  return 0;
}

int RunFlowGraph(gr::top_block_sptr& tb, size_t byteCount)
{

  Rprintf("TB starting, expecting generation of %lu bytes.", size_t(206), byteCount);
  tb->start();

  printf("Press Enter to end:");
  int c = getchar();

  tb->stop();
  tb->wait();
  return 0;
}

std::vector<float> gen_arb_interp_taps(double rate)
{
  std::vector<float> taps;
  if(rate != 1.){
    if(rate < 5./6.){
      double cutoff = 0.4*rate;
      double transb = 0.2*rate;
      taps = gr::filter::firdes::low_pass_2(32.0,32.0,cutoff,transb,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
    }
    else if(rate > 6./5.){
      double cutoff = 0.5/rate;
      double transb = 0.1/rate;
      taps = gr::filter::firdes::low_pass_2(32.0,32.0,cutoff,transb,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
    }
    else if(rate < 1.){
      double cutoff = 0.4*rate;
      double transb = 0.5 - cutoff;
      taps = gr::filter::firdes::low_pass_2(32.0,32.0,cutoff,transb,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
    }
    else{
      double cutoff = 0.4/rate;
      double transb = 0.5 - cutoff;
      taps = gr::filter::firdes::low_pass_2(32.0,32.0,cutoff,transb,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
    }
  }
  else{
    taps = std::vector<float>(1,1.);
  }
  return taps;
}

std::vector<float> gen_int_interp_taps(double rate)
{
  std::vector<float> taps;
  if(rate != 1.){
    taps = gr::filter::firdes::low_pass_2(1,1,0.5/rate,0.5/rate*0.05,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
  }
  else{
    taps = std::vector<float>(1,1.);
  }
  return taps;
}

std::vector<float> gen_synthz_taps(int channels, float samp_rate)
{
  std::vector<float> taps;
  if(TAP_SWITCH==0){
    size_t ntaps = channels*TAPS_PER_CHAN;
    taps = std::vector<float>(ntaps,0.);
    ntaps = (ntaps%2) ? ntaps : ntaps-1;

    std::vector<double> window = kaiser_window(ntaps,SYNTH_BETA);
    std::vector<double> sinc_w(ntaps,0.);
    double s(0.),y(0.);
    for(size_t idx = 0; idx<ntaps; idx++){
      y = (double(idx)*2.-(double(ntaps)-1))/double(channels);
      s = boost::math::sinc_pi(M_PI * y);
      sinc_w[idx] = s;
    }

    for(size_t idx = 0; idx<ntaps; idx++){
      taps[idx] = window[idx]*sinc_w[idx];
    }
  }
  else if(TAP_SWITCH==1){
    float num_channels(channels);
    taps = gr::filter::firdes::low_pass_2(num_channels/2.,num_channels,0.55,0.005,120, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
  }

  return taps;
}

std::vector<float> gen_channl_taps(int channels, float samp_rate)
{
  std::vector<float> taps;
  if(TAP_SWITCH==0){
    size_t ntaps = channels*TAPS_PER_CHAN;
    taps = std::vector<float>(ntaps,0.);
    ntaps = ntaps-1;

    std::vector<double> window = kaiser_window(ntaps,CHANN_BETA);
    std::vector<double> sinc_w(ntaps,0.);
    double s(0.),y(0.);
    for(size_t idx = 0; idx<ntaps; idx++){
      y = (double(idx)-(double(ntaps)-1)/2.)/double(channels);
      s = boost::math::sinc_pi(M_PI * y);
      sinc_w[idx] = s;
    }

    for(size_t idx = 0; idx<ntaps; idx++){
      taps[idx] = window[idx]*sinc_w[idx];
    }
  }
  else if(TAP_SWITCH==1){
    float num_channels(channels);
    //taps = gr::filter::firdes::low_pass_2(num_channels/2.,2.*samp_rate,samp_rate/num_channels,2.*samp_rate/(num_channels*5.),80, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
    switch(channels){
      case 2:
      {
        taps = gr::filter::firdes::low_pass_2(num_channels/2.,num_channels,0.5,0.232,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
        break;
      }
      case 4:
      {
        taps = gr::filter::firdes::low_pass_2(num_channels/2.,num_channels,0.5,0.232,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
        break;
      }
      case 8:
      {
        taps = gr::filter::firdes::low_pass_2(num_channels/2.,num_channels,0.5,0.232,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
        break;
      }
      case 16:
      {
        taps = gr::filter::firdes::low_pass_2(num_channels/2.,num_channels,0.5,0.232,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
        break;
      }
      case 32:
      {
        taps = gr::filter::firdes::low_pass_2(num_channels/2.,num_channels,0.5,0.232,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
        break;
      }
      case 64:
      {
        taps = gr::filter::firdes::low_pass_2(num_channels/2.,num_channels,0.5,0.232,100, gr::filter::firdes::WIN_BLACKMAN_HARRIS);
        break;
      }
    }
  }

  return taps;
}

std::vector<double> kaiser_window(size_t ntaps, double beta)
{
  std::vector<double> window(ntaps,0.);


  double r(0.),x(0.),y(0.),d(0.);

  d = boost::math::cyl_bessel_i(0.,beta);

  for(size_t idx = 0; idx < ntaps; idx++){
    y = (double(idx)-(double(ntaps)-1)/2.)/((double(ntaps)-1)/2.);
    x = beta*sqrt(1-y*y);
    r = boost::math::cyl_bessel_i(0.,x);
    window[idx] = r/d;
  }

  return window;
}




int main(int argc, char** argv)
{

  if(argc < 2){
    std::cerr << "Usage: " << argv[0] << " <config_file>.cfg" << std::endl;
    return 0;
  }

  std::string config_file(argv[1]);


  size_t max_signals = SIGNAL_MAX;
  system_var system_config;
  gr::signal_exciter::sig_params signal_list[max_signals];

  std::string run_name = "Signal Exciter " + config_file;

  gr::top_block_sptr tb = gr::make_top_block(run_name.c_str());

  if(parse_config(config_file, &system_config, signal_list, max_signals)){
    fprintf(stderr, "Parser failure.\n");
    return 1;
  }

  printf("The parser returned, building flowgraph.\n");

  size_t byteCount;

  if(GenerateFlowGraph(system_config, signal_list, tb, byteCount)){
    fprintf(stderr, "Flowgraph generation failure.\n");
    return 2;
  }

  printf("The flowgraph is built, executing flowgraph.\n");

  if(RunFlowGraph(tb, byteCount)){
    fprintf(stderr, "Flowgraph execution failure.\n");
    return 3;
  }

  return 0;
}
















