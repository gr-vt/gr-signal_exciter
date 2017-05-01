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

#include <signal_exciter/api.h>
#include <vector>
#include <gnuradio/analog/cpm.h>

namespace gr {
  namespace signal_exciter {

    enum sig_type_t {
      CWMORSE = 0,
      DSB     = 1,
      DSBSC   = 2,
      USB     = 3,
      LSB     = 4,
      FM      = 5,
      PSK     = 6,
      QAM     = 7,
      PAM     = 8,
      ASK     = 8,
      OFDM    = 9,
      FSK     = 10,
      CPM     = 11,
      GMSK    = 12,
      MSK     = 13,
      GFSK    = 14
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
      float frac_symb_offset;               // Fractional Symbol Offset (pairs with interp)
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

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_CONFIG_H */
