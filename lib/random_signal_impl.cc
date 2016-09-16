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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "random_signal_impl.h"
#include <gnuradio/filter/firdes.h>
#include <stdio.h>

namespace gr {
  namespace signal_exciter {

    random_signal::sptr
    //random_signal::make(int seed)
    random_signal::make(sig_params sig, int seed)
    {
      return gnuradio::get_initial_sptr
        //(new random_signal_impl(seed));
        (new random_signal_impl(sig, seed));
    }

    /*
     * The private constructor
     */
    //random_signal_impl::random_signal_impl(int seed)
    random_signal_impl::random_signal_impl(sig_params sig, int seed)
      : gr::sync_block("random_signal",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(1, 1, sizeof(gr_complex)))
    {
      /*printf("Got to the creation\n");
      roundone = true;*/

      d_params = sig;
      sig_type_t mod_type = d_params.type;
      if(mod_type == FM){
        d_mod = new Signal_FM(d_params.mod_idx,d_params.f_max,d_params.var1,d_params.var2,d_params.thresh,seed);
      }
      else if(mod_type == DSB){
        d_mod = new Signal_DSB(d_params.mod_idx,d_params.f_max,d_params.var1,d_params.var2,d_params.thresh,seed,d_params.am_norm);
      }
      else if(mod_type == DSBSC){
        d_mod = new Signal_DSBSC(d_params.mod_idx,d_params.f_max,d_params.var1,d_params.var2,d_params.thresh,seed,d_params.am_norm);
      }
      else if(mod_type == USB){
        d_mod = new Signal_USB(d_params.mod_idx,d_params.f_max,d_params.var1,d_params.var2,d_params.thresh,seed,d_params.am_norm);
      }
      else if(mod_type == LSB){
        d_mod = new Signal_LSB(d_params.mod_idx,d_params.f_max,d_params.var1,d_params.var2,d_params.thresh,seed,d_params.am_norm);
      }
      else if(mod_type == PSK){
        d_mod = new Signal_PSK(d_params.order,d_params.offset,d_params.sps,&d_params.pulse_shape[0],d_params.pulse_len,seed);
      }
      else if(mod_type == QAM){
        d_mod = new Signal_QAM(d_params.order,d_params.offset,d_params.sps,&d_params.pulse_shape[0],d_params.pulse_len,seed);
      }
      else if(mod_type == PAM){
        d_mod = new Signal_PAM(d_params.order,d_params.offset,d_params.sps,&d_params.pulse_shape[0],d_params.pulse_len,seed);
      }
      else if(mod_type == OFDM){
        d_mod = new Signal_OFDM(d_params.fftsize,d_params.cp_len,d_params.active_carriers,d_params.syms_per_frame,
          d_params.pilot_per_frame, d_params.pilot_count, &d_params.pilot_locations[0], d_params.backoff,
          d_params.mod, d_params.order,d_params.offset, seed, d_params.add_sync,&d_params.pulse_shape[0],d_params.pulse_len,int(d_params.sps));
      }
      else if(mod_type == CWMORSE){
        //printf("cpw = %d\nwpm = %0.0f\nbw = %u\nsr = %lf",sig.char_per_word,sig.words_per_minute,sig.base_word,(double(sig.words_per_minute*(sig.base_word ? 60 : 50))/60.));
        d_mod = new Signal_CWMORSE(d_params.char_per_word,d_params.words_per_minute,d_params.base_word,seed);
      }
      else if(mod_type == MSK){
        d_mod = new Signal_CPM(2,gr::analog::cpm::LREC,d_params.sps,1,0.5,seed);
      }
      else if(mod_type == GMSK){
        d_mod = new Signal_CPM(2,gr::analog::cpm::GAUSSIAN,d_params.sps,d_params.L,0.5,seed,d_params.beta);
      }
      else if(mod_type == FSK){
        d_mod = new Signal_CPM(d_params.order,gr::analog::cpm::LREC,d_params.sps,1,d_params.mod_idx,seed);
      }
      else if(mod_type == GFSK){
        std::vector<float> gt = gr::filter::firdes::gaussian(1,d_params.sps,d_params.beta,d_params.L*d_params.sps);
        std::vector<float> rt(d_params.sps,1.);
        std::vector<float> taps(gt.size()+rt.size()-1);
        for(size_t idx = 0; idx < taps.size(); idx++){
          size_t count = std::min(std::min(idx,rt.size()),gt.size());

          float summer = 0.;
          for(size_t dp = 0; dp < count; dp++){
            if(idx-dp < gt.size())
              summer += gt[idx-dp]*rt[dp];
          }
          taps[idx] = summer;
        }
        d_mod = new Signal_CPM(d_params.order,gr::analog::cpm::GENERIC,d_params.sps,d_params.L,d_params.mod_idx,seed,d_params.beta,&taps[0],taps.size());
      }
      else if(mod_type == CPM){
        gr::analog::cpm::cpm_type ptype = gr::analog::cpm::cpm_type(d_params.phase_type);
        d_mod = new Signal_CPM(d_params.order,ptype,d_params.sps,d_params.L,d_params.mod_idx,seed,d_params.beta);
      }
      else{
        printf("UNKNOWN.\n");
      }
    }

    /*
     * Our virtual destructor.
     */
    random_signal_impl::~random_signal_impl()
    {
      delete d_mod;
    }

    int
    random_signal_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      gr_complex *out = (gr_complex *) output_items[0];

      // Do <+signal processing+>

      d_mod->generate_signal(&out[0], noutput_items);
      //d_mod->generate_symbols(&out[0], noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace signal_exciter */
} /* namespace gr */

