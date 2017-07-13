/* -*- c++ -*- */
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


#ifndef INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_CONFIG_H
#define INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_CONFIG_H

#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <signal_exciter/api.h>
#include <vector>
#include <gnuradio/analog/cpm.h>

namespace gr {
  namespace signal_exciter {

    enum sig_type_t {
      NONE    =-1,
      CWMORSE = 0,
      DSB     = 1,
      DSBSC   = 2,
      USB     = 3,
      LSB     = 4,
      FM      = 5,
      PSK     = 6,
      QAM     = 7,
      PAM     = 8,
      OFDM    = 9,
      FSK     = 10,
      CPM     = 11,
      GMSK    = 12,
      MSK     = 13,
      GFSK    = 14,
      ASK     = 15,
      NCFSK   = 16
    };

    struct SIGNAL_EXCITER_API sig_params
    {
     public:
      // Common to all signals
      sig_type_t type;                      // Signal Type
      float fc;                             // Signal center frequency
      float gain;                           // Signal Gain
      float fs;                             // Signal Sample Rate (note: not always used directly)
      bool ops_gate;                        // Use one pass gate (OPG)?
      bool per_gate;                        // Use periodic gate (PG)?
      bool rnd_gate;                        // Use random gate   (RG)?
      float ops_gate_off;                   // OPG inital off duration
      float ops_gate_on;                    // OPG on duration
      float per_gate_off;                   // PG off duration
      float per_gate_on;                    // PG on duration
      float per_gate_offset;                // PG periodic offset (off/on -> period)
      float rnd_gate_off_min;               // RG minimum off duration
      float rnd_gate_off_max;               // RG maximum off duration
      float rnd_gate_on_min;                // RG minimum on duration
      float rnd_gate_on_max;                // RG maximum on duartion
      // Unique to signal sets
      sig_type_t mod;                       // Base signal type (OFDM)
      int order;                            // Digital Modulation Order
      float offset;                         // Rotate the symbols?
      int sps;                              // Samples per symbol (matched with Pulse Shape/Len)
      std::vector<float> pulse_shape;       // The pulse shaping
      size_t pulse_len;                     // The length of the pulse shaping
      float mod_idx;                        // Modulation index (analog mods)
      size_t fftsize;                       // OFDM fft size
      size_t cp_len;                        // Cyclic Prefix length of OFDM
      size_t active_carriers;               // Number of carriers to use from the fftsize
      size_t syms_per_frame;                // How many OFDM symbols are in a frame
      bool pilot_per_frame;                 // True: Pilots are defined per frame/ False: per symbol
      size_t pilot_count;                   // Total number of pilots
      std::vector<size_t> pilot_locations;  // Index of the pilot locations
      float backoff;                        // OFDM digital backoff to avoid clipping (dB)
      bool add_sync;                        // Add full fftsize Schmidl & Cox sync preamble per frame
      int char_per_word;                    // CWMORSE - how many characters in a 'codeword'
      float words_per_minute;               // How many words (based on base_word) to send per minute
      bool base_word;                       // base_word ? PARIS : CODEX -> set's the Fs of signal
      unsigned L;                           // number of overlapping symbols in CPM mods
      float beta;                           // beta for gaussian & spectral raised cosine
      int phase_type;                       // gr::analog::cpm::cpm_type

      bool am_norm;                         // Enable the agc on the AM signals

      std::vector<float> taper;             // OFDM Specific
      size_t samp_overlap;                  // OFDM Specific

      //reworking analog right now
      //old
      //float f_max;                          // I don't remember (analog mods)
      //float var1;                           // Variance 1 of Gaussian mixture (analog mods)
      //float var2;                           // Variance 2 of Gaussian mixture (analog mods)
      //float thresh;                         // randU > thresh ? dist1 : dist2 (analog mods)
      //new
      std::vector<float> mu;                // mean frequency of each spectral guassian component (analog mods)
      std::vector<float> sigma;             // std of each spectral guassian component (analog mods)
      std::vector<float> weight;            // weights of each spectral guassian component (analog mods)
      size_t components;                    // number of spectral guassian components (analog mods)
      size_t spectral_len;                  // number of points/taps in the spectral tap generation
      float f_max;                          // Maximum frequency at generation (pre-interp)


    };

    class SIGNAL_EXCITER_API signal_parameters
    {
     protected:
      sig_type_t mod;                       // Base signal type (OFDM)
      sig_type_t type;                      // Signal Type
      bool ops_gate;                        // Use one pass gate (OPG)?
      bool per_gate;                        // Use periodic gate (PG)?
      bool rnd_gate;                        // Use random gate   (RG)?
      bool pilot_per_frame;                 // True: Pilots are defined per frame/ False: per symbol
      bool add_sync;                        // Add full fftsize Schmidl & Cox sync preamble per frame
      bool base_word;                       // base_word ? PARIS : CODEX -> set's the Fs of signal
      bool am_norm;                         // Enable the agc on the AM signals
      float fc;                             // Signal center frequency
      float gain;                           // Signal Gain
      float fs;                             // Signal Sample Rate (note: not always used directly)
      float ops_gate_off;                   // OPG inital off duration
      float ops_gate_on;                    // OPG on duration
      float per_gate_off;                   // PG off duration
      float per_gate_on;                    // PG on duration
      float per_gate_offset;                // PG periodic offset (off/on -> period)
      float rnd_gate_off_min;               // RG minimum off duration
      float rnd_gate_off_max;               // RG maximum off duration
      float rnd_gate_on_min;                // RG minimum on duration
      float rnd_gate_on_max;                // RG maximum on duartion
      float offset;                         // Rotate the symbols?
      float mod_idx;                        // Modulation index (analog mods)
      float backoff;                        // OFDM digital backoff to avoid clipping (dB)
      float words_per_minute;               // How many words (based on base_word) to send per
      float beta;                           // beta for gaussian & spectral raised cosine
      float f_max;                          // Maximum frequency at generation (pre-interp)
      int order;                            // Digital Modulation Order
      int sps;                              // Samples per symbol (matched with Pulse Shape/Len)
      int char_per_word;                    // CWMORSE - how many characters in a 'codeword'
      int phase_type;                       // gr::analog::cpm::cpm_type
      unsigned L;                           // number of overlapping symbols in CPM mods
      size_t pulse_len;                     // The length of the pulse shaping
      size_t fftsize;                       // OFDM fft size
      size_t cp_len;                        // Cyclic Prefix length of OFDM
      size_t active_carriers;               // Number of carriers to use from the fftsize
      size_t syms_per_frame;                // How many OFDM symbols are in a frame
      size_t samp_overlap;                  // OFDM Specific
      size_t pilot_count;                   // Total number of pilots
      size_t components;                    // number of spectral guassian components (analog mods)
      size_t spectral_len;                  // number of points/taps in the spectral tap generation
      std::vector<size_t> pilot_locations;  // Index of the pilot locations minute
      std::vector<float> pulse_shape;       // The pulse shaping
      std::vector<float> taper;             // OFDM Specific
      std::vector<float> mu;                // mean frequency of each spectral guassian component (analog mods)
      std::vector<float> sigma;             // std of each spectral guassian component (analog mods)
      std::vector<float> weight;            // weights of each spectral guassian component (analog mods)

     public:
      signal_parameters()
      {
        mod = NONE;
        type = NONE;
        ops_gate = false;
        per_gate = false;
        rnd_gate = false;
        pilot_per_frame = false;
        add_sync = false;
        base_word = false;
        am_norm = false;
        fc = 0.;
        gain = 0.;
        fs = 0.;
        ops_gate_off = 0.;
        ops_gate_on = 0.;
        per_gate_off = 0.;
        per_gate_on = 0.;
        per_gate_offset = 0.;
        rnd_gate_off_min = 0.;
        rnd_gate_off_max = 0.;
        rnd_gate_on_min = 0.;
        rnd_gate_on_max = 0.;
        offset = 0.;
        mod_idx = 0.;
        backoff = 0.;
        words_per_minute = 0.;
        beta = 0.;
        f_max = 0.;
        order = 0;
        sps = 0;
        char_per_word = 0;
        phase_type = 0;
        L = 0;
        pulse_len = 0;
        fftsize = 0;
        cp_len = 0;
        active_carriers = 0;
        syms_per_frame = 0;
        samp_overlap = 0;
        pilot_count = 0;
        components = 0;
        spectral_len = 0;
        pilot_locations = std::vector<size_t>(0);
        pulse_shape = std::vector<float>(0);
        taper = std::vector<float>(0);
        mu = std::vector<float>(0);
        sigma = std::vector<float>(0);
        weight = std::vector<float>(0);
      }

      signal_parameters(const signal_parameters &rhs)
      {
        mod = rhs.mod;
        type = rhs.type;
        ops_gate = rhs.ops_gate;
        per_gate = rhs.per_gate;
        rnd_gate = rhs.rnd_gate;
        pilot_per_frame = rhs.pilot_per_frame;
        add_sync = rhs.add_sync;
        base_word = rhs.base_word;
        am_norm = rhs.am_norm;
        fc = rhs.fc;
        gain = rhs.gain;
        fs = rhs.fs;
        ops_gate_off = rhs.ops_gate_off;
        ops_gate_on = rhs.ops_gate_on;
        per_gate_off = rhs.per_gate_off;
        per_gate_on = rhs.per_gate_on;
        per_gate_offset = rhs.per_gate_offset;
        rnd_gate_off_min = rhs.rnd_gate_off_min;
        rnd_gate_off_max = rhs.rnd_gate_off_max;
        rnd_gate_on_min = rhs.rnd_gate_on_min;
        rnd_gate_on_max = rhs.rnd_gate_on_max;
        offset = rhs.offset;
        mod_idx = rhs.mod_idx;
        backoff = rhs.backoff;
        words_per_minute = rhs.words_per_minute;
        beta = rhs.beta;
        f_max = rhs.f_max;
        order = rhs.order;
        sps = rhs.sps;
        char_per_word = rhs.char_per_word;
        phase_type = rhs.phase_type;
        L = rhs.L;
        pulse_len = rhs.pulse_len;
        fftsize = rhs.fftsize;
        cp_len = rhs.cp_len;
        active_carriers = rhs.active_carriers;
        syms_per_frame = rhs.syms_per_frame;
        samp_overlap = rhs.samp_overlap;
        pilot_count = rhs.pilot_count;
        components = rhs.components;
        spectral_len = rhs.spectral_len;
        pilot_locations = std::vector<size_t>(rhs.pilot_locations.begin(), rhs.pilot_locations.end());
        pulse_shape = std::vector<float>(rhs.pulse_shape.begin(), rhs.pulse_shape.end());
        taper = std::vector<float>(rhs.taper.begin(), rhs.taper.end());
        mu = std::vector<float>(rhs.mu.begin(), rhs.mu.end());
        sigma = std::vector<float>(rhs.sigma.begin(), rhs.sigma.end());
        weight = std::vector<float>(rhs.weight.begin(), rhs.weight.end());
      }

      signal_parameters& operator = (const signal_parameters &rhs)
      {
        mod = rhs.mod;
        type = rhs.type;
        ops_gate = rhs.ops_gate;
        per_gate = rhs.per_gate;
        rnd_gate = rhs.rnd_gate;
        pilot_per_frame = rhs.pilot_per_frame;
        add_sync = rhs.add_sync;
        base_word = rhs.base_word;
        am_norm = rhs.am_norm;
        fc = rhs.fc;
        gain = rhs.gain;
        fs = rhs.fs;
        ops_gate_off = rhs.ops_gate_off;
        ops_gate_on = rhs.ops_gate_on;
        per_gate_off = rhs.per_gate_off;
        per_gate_on = rhs.per_gate_on;
        per_gate_offset = rhs.per_gate_offset;
        rnd_gate_off_min = rhs.rnd_gate_off_min;
        rnd_gate_off_max = rhs.rnd_gate_off_max;
        rnd_gate_on_min = rhs.rnd_gate_on_min;
        rnd_gate_on_max = rhs.rnd_gate_on_max;
        offset = rhs.offset;
        mod_idx = rhs.mod_idx;
        backoff = rhs.backoff;
        words_per_minute = rhs.words_per_minute;
        beta = rhs.beta;
        f_max = rhs.f_max;
        order = rhs.order;
        sps = rhs.sps;
        char_per_word = rhs.char_per_word;
        phase_type = rhs.phase_type;
        L = rhs.L;
        pulse_len = rhs.pulse_len;
        fftsize = rhs.fftsize;
        cp_len = rhs.cp_len;
        active_carriers = rhs.active_carriers;
        syms_per_frame = rhs.syms_per_frame;
        samp_overlap = rhs.samp_overlap;
        pilot_count = rhs.pilot_count;
        components = rhs.components;
        spectral_len = rhs.spectral_len;
        pilot_locations = std::vector<size_t>(rhs.pilot_locations.begin(), rhs.pilot_locations.end());
        pulse_shape = std::vector<float>(rhs.pulse_shape.begin(), rhs.pulse_shape.end());
        taper = std::vector<float>(rhs.taper.begin(), rhs.taper.end());
        mu = std::vector<float>(rhs.mu.begin(), rhs.mu.end());
        sigma = std::vector<float>(rhs.sigma.begin(), rhs.sigma.end());
        weight = std::vector<float>(rhs.weight.begin(), rhs.weight.end());
        return *this;
      }

      ~signal_parameters()
      {
      }

      void set_mod(sig_type_t new_mod){ mod = new_mod;}
      void set_type(sig_type_t new_type){ type = new_type;}
      void set_ops_gate(bool new_ops_gate){ ops_gate = new_ops_gate;}
      void set_per_gate(bool new_per_gate){ per_gate = new_per_gate;}
      void set_rnd_gate(bool new_rnd_gate){ rnd_gate = new_rnd_gate;}
      void set_pilot_per_frame(bool new_pilot_per_frame){ pilot_per_frame = new_pilot_per_frame;}
      void set_add_sync(bool new_add_sync){ add_sync = new_add_sync;}
      void set_base_word(bool new_base_word){ base_word = new_base_word;}
      void set_am_norm(bool new_am_norm){ am_norm = new_am_norm;}
      void set_fc(float new_fc){ fc = new_fc;}
      void set_gain(float new_gain){ gain = new_gain;}
      void set_fs(float new_fs){ fs = new_fs;}
      void set_ops_gate_off(float new_ops_gate_off){ ops_gate_off = new_ops_gate_off;}
      void set_ops_gate_on(float new_ops_gate_on){ ops_gate_on = new_ops_gate_on;}
      void set_per_gate_off(float new_per_gate_off){ per_gate_off = new_per_gate_off;}
      void set_per_gate_on(float new_per_gate_on){ per_gate_on = new_per_gate_on;}
      void set_per_gate_offset(float new_per_gate_offset){ per_gate_offset = new_per_gate_offset;}
      void set_rnd_gate_off_min(float new_rnd_gate_off_min){ rnd_gate_off_min = new_rnd_gate_off_min;}
      void set_rnd_gate_off_max(float new_rnd_gate_off_max){ rnd_gate_off_max = new_rnd_gate_off_max;}
      void set_rnd_gate_on_min(float new_rnd_gate_on_min){ rnd_gate_on_min = new_rnd_gate_on_min;}
      void set_rnd_gate_on_max(float new_rnd_gate_on_max){ rnd_gate_on_max = new_rnd_gate_on_max;}
      void set_offset(float new_offset){ offset = new_offset;}
      void set_mod_idx(float new_mod_idx){ mod_idx = new_mod_idx;}
      void set_backoff(float new_backoff){ backoff = new_backoff;}
      void set_words_per_minute(float new_words_per_minute){ words_per_minute = new_words_per_minute;}
      void set_beta(float new_beta){ beta = new_beta;}
      void set_f_max(float new_f_max){ f_max = new_f_max;}
      void set_order(int new_order){ order = new_order;}
      void set_sps(int new_sps){ sps = new_sps;}
      void set_char_per_word(int new_char_per_word){ char_per_word = new_char_per_word;}
      void set_phase_type(int new_phase_type){ phase_type = new_phase_type;}
      void set_L(unsigned new_L){ L = new_L;}
      void set_pulse_len(size_t new_pulse_len){ pulse_len = new_pulse_len;}
      void set_fftsize(size_t new_fftsize){ fftsize = new_fftsize;}
      void set_cp_len(size_t new_cp_len){ cp_len = new_cp_len;}
      void set_active_carriers(size_t new_active_carriers){ active_carriers = new_active_carriers;}
      void set_syms_per_frame(size_t new_syms_per_frame){ syms_per_frame = new_syms_per_frame;}
      void set_samp_overlap(size_t new_samp_overlap){ samp_overlap = new_samp_overlap;}
      void set_pilot_count(size_t new_pilot_count){ pilot_count = new_pilot_count;}
      void set_components(size_t new_components){ components = new_components;}
      void set_spectral_len(size_t new_spectral_len){ spectral_len = new_spectral_len;}
      void set_pilot_locations(std::vector<size_t> &new_pilot_locations){ pilot_locations = std::vector<size_t>(new_pilot_locations.begin(), new_pilot_locations.end());}
      void set_pulse_shape(std::vector<float> &new_pulse_shape){ pulse_shape = std::vector<float>(new_pulse_shape.begin(),new_pulse_shape.end());}
      void set_taper(std::vector<float> &new_taper){ taper = std::vector<float>(new_taper.begin(),new_taper.end());}
      void set_mu(std::vector<float> &new_mu){ mu = std::vector<float>(new_mu.begin(),new_mu.end());}
      void set_sigma(std::vector<float> &new_sigma){ sigma = std::vector<float>(new_sigma.begin(),new_sigma.end());}
      void set_weight(std::vector<float> &new_weight){ weight = std::vector<float>(new_weight.begin(),new_weight.end());}
      void set_pilot_locations(const std::vector<size_t> &new_pilot_locations){ pilot_locations = std::vector<size_t>(new_pilot_locations.begin(), new_pilot_locations.end());}
      void set_pulse_shape(const std::vector<float> &new_pulse_shape){ pulse_shape = std::vector<float>(new_pulse_shape.begin(),new_pulse_shape.end());}
      void set_taper(const std::vector<float> &new_taper){ taper = std::vector<float>(new_taper.begin(),new_taper.end());}
      void set_mu(const std::vector<float> &new_mu){ mu = std::vector<float>(new_mu.begin(),new_mu.end());}
      void set_sigma(const std::vector<float> &new_sigma){ sigma = std::vector<float>(new_sigma.begin(),new_sigma.end());}
      void set_weight(const std::vector<float> &new_weight){ weight = std::vector<float>(new_weight.begin(),new_weight.end());}

      sig_type_t get_mod(){ return mod; }
      sig_type_t get_type(){ return type; }
      bool get_ops_gate(){ return ops_gate; }
      bool get_per_gate(){ return per_gate; }
      bool get_rnd_gate(){ return rnd_gate; }
      bool get_pilot_per_frame(){ return pilot_per_frame; }
      bool get_add_sync(){ return add_sync; }
      bool get_base_word(){ return base_word; }
      bool get_am_norm(){ return am_norm; }
      float get_fc(){ return fc; }
      float get_gain(){ return gain; }
      float get_fs(){ return fs; }
      float get_ops_gate_off(){ return ops_gate_off; }
      float get_ops_gate_on(){ return ops_gate_on; }
      float get_per_gate_off(){ return per_gate_off; }
      float get_per_gate_on(){ return per_gate_on; }
      float get_per_gate_offset(){ return per_gate_offset; }
      float get_rnd_gate_off_min(){ return rnd_gate_off_min; }
      float get_rnd_gate_off_max(){ return rnd_gate_off_max; }
      float get_rnd_gate_on_min(){ return rnd_gate_on_min; }
      float get_rnd_gate_on_max(){ return rnd_gate_on_max; }
      float get_offset(){ return offset; }
      float get_mod_idx(){ return mod_idx; }
      float get_backoff(){ return backoff; }
      float get_words_per_minute(){ return words_per_minute; }
      float get_beta(){ return beta; }
      float get_f_max(){ return f_max; }
      int get_order(){ return order; }
      int get_sps(){ return sps; }
      int get_char_per_word(){ return char_per_word; }
      int get_phase_type(){ return phase_type; }
      unsigned get_L(){ return L; }
      size_t get_pulse_len(){ return pulse_len; }
      size_t get_fftsize(){ return fftsize; }
      size_t get_cp_len(){ return cp_len; }
      size_t get_active_carriers(){ return active_carriers; }
      size_t get_syms_per_frame(){ return syms_per_frame; }
      size_t get_samp_overlap(){ return samp_overlap; }
      size_t get_pilot_count(){ return pilot_count; }
      size_t get_components(){ return components; }
      size_t get_spectral_len(){ return spectral_len; }
      std::vector<size_t> get_pilot_locations(){ return pilot_locations; }
      std::vector<float> get_pulse_shape(){ return pulse_shape; }
      std::vector<float> get_taper(){ return taper; }
      std::vector<float> get_mu(){ return mu; }
      std::vector<float> get_sigma(){ return sigma; }
      std::vector<float> get_weight(){ return weight; }
      size_t* get_pilot_locations_ptr(){ return &pilot_locations[0]; }
      float* get_pulse_shape_ptr(){ return &pulse_shape[0]; }
      float* get_taper_ptr(){ return &taper[0]; }
      float* get_mu_ptr(){ return &mu[0]; }
      float* get_sigma_ptr(){ return &sigma[0]; }
      float* get_weight_ptr(){ return &weight[0]; }
      size_t get_pilot_locations_size(){ return pilot_locations.size(); }
      size_t get_pulse_shape_size(){ return pulse_shape.size(); }
      size_t get_taper_size(){ return taper.size(); }
      size_t get_mu_size(){ return mu.size(); }
      size_t get_sigma_size(){ return sigma.size(); }
      size_t get_weight_size(){ return weight.size(); }

      void print_params() const
      {
        printf("mod = %d\n",mod);
        printf("type = %d\n",type);
        printf("ops_gate = %d\n",ops_gate);
        printf("per_gate = %d\n",per_gate);
        printf("rnd_gate = %d\n",rnd_gate);
        printf("pilot_per_frame = %d\n",pilot_per_frame);
        printf("add_sync = %d\n",add_sync);
        printf("base_word = %d\n",base_word);
        printf("am_norm = %d\n",am_norm);
        printf("fc = %1.8e\n",fc);
        printf("gain = %1.8e\n",gain);
        printf("fs = %1.8e\n",fs);
        printf("ops_gate_off = %1.8e\n",ops_gate_off);
        printf("ops_gate_on = %1.8e\n",ops_gate_on);
        printf("per_gate_off = %1.8e\n",per_gate_off);
        printf("per_gate_on = %1.8e\n",per_gate_on);
        printf("per_gate_offset = %1.8e\n",per_gate_offset);
        printf("rnd_gate_off_min = %1.8e\n",rnd_gate_off_min);
        printf("rnd_gate_off_max = %1.8e\n",rnd_gate_off_max);
        printf("rnd_gate_on_min = %1.8e\n",rnd_gate_on_min);
        printf("rnd_gate_on_max = %1.8e\n",rnd_gate_on_max);
        printf("offset = %1.8e\n",offset);
        printf("mod_idx = %1.8e\n",mod_idx);
        printf("backoff = %1.8e\n",backoff);
        printf("words_per_minute = %1.8e\n",words_per_minute);
        printf("beta = %1.8e\n",beta);
        printf("f_max = %1.8e\n",f_max);
        printf("order = %d\n",order);
        printf("sps = %d\n",sps);
        printf("char_per_word = %d\n",char_per_word);
        printf("phase_type = %d\n",phase_type);
        printf("L = %u\n",L);
        printf("pulse_len = %lu\n",pulse_len);
        printf("fftsize = %lu\n",fftsize);
        printf("cp_len = %lu\n",cp_len);
        printf("active_carriers = %lu\n",active_carriers);
        printf("syms_per_frame = %lu\n",syms_per_frame);
        printf("samp_overlap = %lu\n",samp_overlap);
        printf("pilot_count = %lu\n",pilot_count);
        printf("components = %lu\n",components);
        printf("spectral_len = %lu\n",spectral_len);
        printf("pilot_locations = [");
        for(size_t idx = 0; idx < pilot_locations.size(); idx++){
          printf(" %lu",pilot_locations[idx]);
        }
        printf("]\n");
        printf("pulse_shape = [");
        for(size_t idx = 0; idx < pulse_shape.size(); idx++){
          printf(" %1.8e",pulse_shape[idx]);
        }
        printf("]\n");
        printf("taper = [");
        for(size_t idx = 0; idx < taper.size(); idx++){
          printf(" %1.8e",taper[idx]);
        }
        printf("]\n");
        printf("mu = [");
        for(size_t idx = 0; idx < mu.size(); idx++){
          printf(" %1.8e",mu[idx]);
        }
        printf("]\n");
        printf("sigma = [");
        for(size_t idx = 0; idx < sigma.size(); idx++){
          printf(" %1.8e",sigma[idx]);
        }
        printf("]\n");
        printf("weight = [");
        for(size_t idx = 0; idx < weight.size(); idx++){
          printf(" %1.8e",weight[idx]);
        }
        printf("]\n");
      }

    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_CONFIG_H */
