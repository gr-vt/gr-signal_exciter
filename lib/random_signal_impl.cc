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
#include "pulseshapedes.h"
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

    random_signal::sptr
    //random_signal::make(int seed)
    random_signal::make(const signal_parameters &sig, int seed)
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
#ifndef SIGNAL_EXCITER_USING_GNURADIO_FFT_MUTEX
#define SIGNAL_EXCITER_USING_GNURADIO_FFT_MUTEX
      Signal_Base::set_fftw_mutex(&(gr::fft::planner::mutex()));
#endif //SIGNAL_EXCITER_USING_GNURADIO_FFT_MUTEX
      /*printf("Got to the creation\n");
      roundone = true;*/
      //set_output_multiple(4);

      d_params = sig;
      sig_type_t mod_type = d_params.type;
      if(mod_type == FM){
        d_mod = new Signal_FM(d_params.mod_idx,d_params.components,
                              &d_params.mu[0],&d_params.sigma[0],
                              &d_params.weight[0],d_params.f_max,
                              d_params.spectral_len,seed,
                              &d_params.pulse_shape[0],d_params.pulse_len,
                              d_params.sps);
      }
      else if(mod_type == DSB){
        d_mod = new Signal_DSB(d_params.mod_idx,d_params.components,
                              &d_params.mu[0],&d_params.sigma[0],
                              &d_params.weight[0],d_params.f_max,
                              d_params.spectral_len,seed,d_params.am_norm,
                              &d_params.pulse_shape[0],d_params.pulse_len,
                              d_params.sps);
      }
      else if(mod_type == DSBSC){
        d_mod = new Signal_DSBSC(d_params.mod_idx,d_params.components,
                                &d_params.mu[0],&d_params.sigma[0],
                                &d_params.weight[0],d_params.f_max,
                                d_params.spectral_len,seed,d_params.am_norm,
                                &d_params.pulse_shape[0],d_params.pulse_len,
                                d_params.sps);
      }
      else if(mod_type == USB){
        d_mod = new Signal_USB(d_params.mod_idx,d_params.components,
                              &d_params.mu[0],&d_params.sigma[0],
                              &d_params.weight[0],d_params.f_max,
                              d_params.spectral_len,seed,d_params.am_norm,
                              &d_params.pulse_shape[0],d_params.pulse_len,
                              d_params.sps);
      }
      else if(mod_type == LSB){
        d_mod = new Signal_LSB(d_params.mod_idx,d_params.components,
                              &d_params.mu[0],&d_params.sigma[0],
                              &d_params.weight[0],d_params.f_max,
                              d_params.spectral_len,seed,d_params.am_norm,
                              &d_params.pulse_shape[0],d_params.pulse_len,
                              d_params.sps);
      }
      else if(mod_type == PSK){
        d_mod = new Signal_PSK(d_params.order,d_params.offset,d_params.sps,
                              &d_params.pulse_shape[0],d_params.pulse_len,
                              seed);
      }
      else if(mod_type == QAM){
        d_mod = new Signal_QAM(d_params.order,d_params.offset,d_params.sps,
                              &d_params.pulse_shape[0],d_params.pulse_len,
                              seed);
      }
      else if(mod_type == PAM){
        d_mod = new Signal_PAM(d_params.order,d_params.offset,d_params.sps,
                              &d_params.pulse_shape[0],d_params.pulse_len,
                              seed);
      }
      else if(mod_type == ASK){
        d_mod = new Signal_ASK(d_params.order,d_params.offset,d_params.sps,
                              &d_params.pulse_shape[0],d_params.pulse_len,
                              seed);
      }
      else if(mod_type == OFDM){
        d_mod = new Signal_OFDM(d_params.fftsize,d_params.cp_len,
                                d_params.active_carriers,
                                d_params.syms_per_frame,
                                d_params.pilot_per_frame,
                                d_params.pilot_count,
                                &d_params.pilot_locations[0],
                                d_params.backoff,
                                d_params.mod, d_params.order,
                                d_params.offset, seed, d_params.add_sync,
                                &d_params.taper[0], d_params.samp_overlap,
                                &d_params.pulse_shape[0], d_params.pulse_len,
                                int(d_params.sps));
      }
      else if(mod_type == CWMORSE){
        //printf("cpw = %d\nwpm = %0.0f\nbw = %u\nsr = %lf",sig.char_per_word,sig.words_per_minute,sig.base_word,(double(sig.words_per_minute*(sig.base_word ? 60 : 50))/60.));
        d_mod = new Signal_CWMORSE(d_params.char_per_word,
                                  d_params.words_per_minute,
                                  d_params.base_word,seed,
                                  &d_params.pulse_shape[0],
                                  d_params.pulse_len,d_params.sps);
      }
      else if(mod_type == MSK){
        if(d_params.pulse_len==0){
          d_mod = new Signal_CPM(2,gr::analog::cpm::LREC,
                                  d_params.sps,1,0.5,seed,0.,NULL,0,NULL,0);
        }
        else{
          d_mod = new Signal_CPM(2,gr::analog::cpm::LREC,
                                  d_params.sps,1,0.5,seed,0.,NULL,0,
                                  &d_params.pulse_shape[0],d_params.pulse_len);
        }
      }
      else if(mod_type == GMSK){
        if(d_params.pulse_len==0){
          d_mod = new Signal_CPM(2,gr::analog::cpm::GAUSSIAN,
                                  d_params.sps,d_params.L,0.5,
                                  seed,d_params.beta,NULL,0,NULL,0);
        }
        else{
          d_mod = new Signal_CPM(2,gr::analog::cpm::GAUSSIAN,
                                  d_params.sps,d_params.L,0.5,
                                  seed,d_params.beta,NULL,0,
                                  &d_params.pulse_shape[0],d_params.pulse_len);
        }
      }
      else if(mod_type == FSK){
        if(d_params.pulse_len==0){
          d_mod = new Signal_CPM(d_params.order,gr::analog::cpm::LREC,
                                  d_params.sps,1,d_params.mod_idx,seed,
                                  0.,NULL,0,NULL,0);
        }
        else{
          d_mod = new Signal_CPM(d_params.order,gr::analog::cpm::LREC,
                                  d_params.sps,1,d_params.mod_idx,seed,
                                  0.,NULL,0,&d_params.pulse_shape[0],d_params.pulse_len);
        }
      }
      else if(mod_type == GFSK){
        std::vector<float> gt = gr::filter::firdes::gaussian(1,d_params.sps,d_params.beta,d_params.L*d_params.sps);
        std::vector<float> rt(d_params.sps,1.);
        std::vector<float> taps;
        conv(taps,gt,rt);
        if(d_params.pulse_len==0){
          d_mod = new Signal_CPM(d_params.order,gr::analog::cpm::GENERIC,
                                d_params.sps,d_params.L,d_params.mod_idx,
                                seed,d_params.beta,&taps[0],taps.size(),NULL,0);
        }
        else{
          d_mod = new Signal_CPM(d_params.order,gr::analog::cpm::GENERIC,
                                d_params.sps,d_params.L,d_params.mod_idx,
                                seed,d_params.beta,&taps[0],taps.size(),
                                &d_params.pulse_shape[0],d_params.pulse_len);
        }
      }
      else if(mod_type == C4FM){
        std::vector<float> rct;
        pulseshapedes::rcos_pulse( rct, float(d_params.L), d_params.sps, .2f );
        std::vector<float> rt(d_params.sps,1.);
        std::vector<float> taps;
        conv(taps,rct,rt);
        if(d_params.pulse_len==0){
          d_mod = new Signal_CPM(4,gr::analog::cpm::GENERIC,
                                d_params.sps,d_params.L,d_params.mod_idx,
                                seed,0.2,&taps[0],taps.size(),NULL,0);
        }
        else{
          d_mod = new Signal_CPM(4,gr::analog::cpm::GENERIC,
                                d_params.sps,d_params.L,d_params.mod_idx,
                                seed,0.2,&taps[0],taps.size(),
                                &d_params.pulse_shape[0],d_params.pulse_len);
        }
      }
      else if(mod_type == CPM){
        gr::analog::cpm::cpm_type ptype =
                        gr::analog::cpm::cpm_type(d_params.phase_type);
        if(d_params.pulse_len==0){
          d_mod = new Signal_CPM(d_params.order,ptype,d_params.sps,
                                d_params.L,d_params.mod_idx,seed,
                                d_params.beta,NULL,0,NULL,0);
        }
        else{
          d_mod = new Signal_CPM(d_params.order,ptype,d_params.sps,
                                d_params.L,d_params.mod_idx,seed,
                                d_params.beta,NULL,0,
                                &d_params.pulse_shape[0],d_params.pulse_len);
        }
      }
      else if(mod_type == NCFSK){
        if(d_params.pulse_len==0){
          d_mod = new Signal_FSK(d_params.order,d_params.sps,d_params.mod_idx,
                                  seed,NULL,0);
        }
        else{
          d_mod = new Signal_FSK(d_params.order,d_params.sps,d_params.mod_idx,
                                  seed,&d_params.pulse_shape[0],d_params.pulse_len);
        }
      }
      else{
        printf("UNKNOWN.\n");
      }
    }

    /*
     * The private constructor
     */
    //random_signal_impl::random_signal_impl(int seed)
    random_signal_impl::random_signal_impl(const signal_parameters &sig, int seed)
      : gr::sync_block("random_signal",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(1, 1, sizeof(gr_complex)))
    {
#ifndef SIGNAL_EXCITER_USING_GNURADIO_FFT_MUTEX
#define SIGNAL_EXCITER_USING_GNURADIO_FFT_MUTEX
      Signal_Base::set_fftw_mutex(&(gr::fft::planner::mutex()));
#endif //SIGNAL_EXCITER_USING_GNURADIO_FFT_MUTEX
      /*printf("Got to the creation\n");
      roundone = true;*/
      //set_output_multiple(4);

      d_parameters = sig;

      sig_type_t mod_type = d_parameters.get_type();
      if(mod_type == FM){
        d_mod = new Signal_FM(d_parameters.get_mod_idx(),d_parameters.get_components(),
                              d_parameters.get_mu_ptr(),d_parameters.get_sigma_ptr(),
                              d_parameters.get_weight_ptr(),d_parameters.get_f_max(),
                              d_parameters.get_spectral_len(),seed,
                              d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len(),
                              d_parameters.get_sps());
      }
      else if(mod_type == DSB){
        d_mod = new Signal_DSB(d_parameters.get_mod_idx(),d_parameters.get_components(),
                              d_parameters.get_mu_ptr(),d_parameters.get_sigma_ptr(),
                              d_parameters.get_weight_ptr(),d_parameters.get_f_max(),
                              d_parameters.get_spectral_len(),seed,d_parameters.get_am_norm(),
                              d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len(),
                              d_parameters.get_sps());
      }
      else if(mod_type == DSBSC){
        d_mod = new Signal_DSBSC(d_parameters.get_mod_idx(),d_parameters.get_components(),
                                d_parameters.get_mu_ptr(),d_parameters.get_sigma_ptr(),
                                d_parameters.get_weight_ptr(),d_parameters.get_f_max(),
                                d_parameters.get_spectral_len(),seed,d_parameters.get_am_norm(),
                                d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len(),
                                d_parameters.get_sps());
      }
      else if(mod_type == USB){
        d_mod = new Signal_USB(d_parameters.get_mod_idx(),d_parameters.get_components(),
                              d_parameters.get_mu_ptr(),d_parameters.get_sigma_ptr(),
                              d_parameters.get_weight_ptr(),d_parameters.get_f_max(),
                              d_parameters.get_spectral_len(),seed,d_parameters.get_am_norm(),
                              d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len(),
                              d_parameters.get_sps());
      }
      else if(mod_type == LSB){
        d_mod = new Signal_LSB(d_parameters.get_mod_idx(),d_parameters.get_components(),
                              d_parameters.get_mu_ptr(),d_parameters.get_sigma_ptr(),
                              d_parameters.get_weight_ptr(),d_parameters.get_f_max(),
                              d_parameters.get_spectral_len(),seed,d_parameters.get_am_norm(),
                              d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len(),
                              d_parameters.get_sps());
      }
      else if(mod_type == PSK){
        //printf("rs: psk: fso: %1.3e\n",d_parameters.get_frac_offset(),d_parameters.get_frac_symb_offset());
        d_mod = new Signal_PSK(d_parameters.get_order(),d_parameters.get_offset(),d_parameters.get_sps(),
                              d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len(),
                              seed);
      }
      else if(mod_type == QAM){
        d_mod = new Signal_QAM(d_parameters.get_order(),d_parameters.get_offset(),d_parameters.get_sps(),
                              d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len(),
                              seed);
      }
      else if(mod_type == PAM){
        d_mod = new Signal_PAM(d_parameters.get_order(),d_parameters.get_offset(),d_parameters.get_sps(),
                              d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len(),
                              seed);
      }
      else if(mod_type == ASK){
        d_mod = new Signal_ASK(d_parameters.get_order(),d_parameters.get_offset(),d_parameters.get_sps(),
                              d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len(),
                              seed);
      }
      else if(mod_type == OFDM){
        d_mod = new Signal_OFDM(d_parameters.get_fftsize(),d_parameters.get_cp_len(),
                                d_parameters.get_active_carriers(),
                                d_parameters.get_syms_per_frame(),
                                d_parameters.get_pilot_per_frame(),
                                d_parameters.get_pilot_count(),
                                d_parameters.get_pilot_locations_ptr(),
                                d_parameters.get_backoff(),
                                d_parameters.get_mod(), d_parameters.get_order(),
                                d_parameters.get_offset(), seed, d_parameters.get_add_sync(),
                                d_parameters.get_taper_ptr(), d_parameters.get_samp_overlap(),
                                d_parameters.get_pulse_shape_ptr(), d_parameters.get_pulse_len(),
                                int(d_parameters.get_sps()));
      }
      else if(mod_type == CWMORSE){
        //printf("cpw = %d\nwpm = %0.0f\nbw = %u\nsr = %lf",sig.char_per_word,sig.words_per_minute,sig.base_word,(double(sig.words_per_minute*(sig.base_word ? 60 : 50))/60.));
        d_mod = new Signal_CWMORSE(d_parameters.get_char_per_word(),
                                  d_parameters.get_words_per_minute(),
                                  d_parameters.get_base_word(),seed,
                                  d_parameters.get_pulse_shape_ptr(),
                                  d_parameters.get_pulse_len(),d_parameters.get_sps());
      }
      else if(mod_type == MSK){
        if(d_parameters.get_pulse_len()==0){
          d_mod = new Signal_CPM(2,gr::analog::cpm::LREC,
                                  d_parameters.get_sps(),1,0.5,seed,0.,NULL,0,NULL,0);
        }
        else{
          d_mod = new Signal_CPM(2,gr::analog::cpm::LREC,
                                  d_parameters.get_sps(),1,0.5,seed,0.,NULL,0,
                                  d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len());
        }
      }
      else if(mod_type == GMSK){
        if(d_parameters.get_pulse_len()==0){
          d_mod = new Signal_CPM(2,gr::analog::cpm::GAUSSIAN,
                                  d_parameters.get_sps(),d_parameters.get_L(),0.5,
                                  seed,d_parameters.get_beta(),NULL,0,NULL,0);
        }
        else{
          d_mod = new Signal_CPM(2,gr::analog::cpm::GAUSSIAN,
                                  d_parameters.get_sps(),d_parameters.get_L(),0.5,
                                  seed,d_parameters.get_beta(),NULL,0,
                                  d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len());
        }
      }
      else if(mod_type == FSK){
        if(d_parameters.get_pulse_len()==0){
          d_mod = new Signal_CPM(d_parameters.get_order(),gr::analog::cpm::LREC,
                                  d_parameters.get_sps(),1,d_parameters.get_mod_idx(),seed,
                                  0.,NULL,0,NULL,0);
        }
        else{
          d_mod = new Signal_CPM(d_parameters.get_order(),gr::analog::cpm::LREC,
                                  d_parameters.get_sps(),1,d_parameters.get_mod_idx(),seed,
                                  0.,NULL,0,d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len());
        }
      }
      else if(mod_type == GFSK){
        std::vector<float> gt = gr::filter::firdes::gaussian(1,d_parameters.get_sps(),d_parameters.get_beta(),d_parameters.get_L()*d_parameters.get_sps());
        std::vector<float> rt(d_parameters.get_sps(),1.);
        std::vector<float> taps;
        conv(taps,gt,rt);
        if(d_parameters.get_pulse_len()==0){
          d_mod = new Signal_CPM(d_parameters.get_order(),gr::analog::cpm::GENERIC,
                                d_parameters.get_sps(),d_parameters.get_L(),d_parameters.get_mod_idx(),
                                seed,d_parameters.get_beta(),&taps[0],taps.size(),NULL,0);
        }
        else{
          d_mod = new Signal_CPM(d_parameters.get_order(),gr::analog::cpm::GENERIC,
                                d_parameters.get_sps(),d_parameters.get_L(),d_parameters.get_mod_idx(),
                                seed,d_parameters.get_beta(),&taps[0],taps.size(),
                                d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len());
        }
      }
      else if(mod_type == C4FM){
        std::vector<float> rct;
        pulseshapedes::rcos_pulse( rct, float(d_parameters.get_L()), d_parameters.get_sps(), 0.2f );
        std::vector<float> rt(d_parameters.get_sps(),1.);
        std::vector<float> taps;
        conv(taps,rct,rt);
        if(d_parameters.get_pulse_len()==0){
          d_mod = new Signal_CPM(4,gr::analog::cpm::GENERIC,
                                d_parameters.get_sps(),d_parameters.get_L(),d_parameters.get_mod_idx(),
                                seed,0.2,&taps[0],taps.size(),NULL,0);
        }
        else{
          d_mod = new Signal_CPM(4,gr::analog::cpm::GENERIC,
                                d_parameters.get_sps(),d_parameters.get_L(),d_parameters.get_mod_idx(),
                                seed,0.2,&taps[0],taps.size(),
                                d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len());
        }
      }
      else if(mod_type == CPM){
        gr::analog::cpm::cpm_type ptype =
                        gr::analog::cpm::cpm_type(d_parameters.get_phase_type());
        if(d_parameters.get_pulse_len()==0){
          d_mod = new Signal_CPM(d_parameters.get_order(),ptype,d_parameters.get_sps(),
                                d_parameters.get_L(),d_parameters.get_mod_idx(),seed,
                                d_parameters.get_beta(),NULL,0,NULL,0);
        }
        else{
          d_mod = new Signal_CPM(d_parameters.get_order(),ptype,d_parameters.get_sps(),
                                d_parameters.get_L(),d_parameters.get_mod_idx(),seed,
                                d_parameters.get_beta(),NULL,0,
                                d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len());
        }
      }
      else if(mod_type == NCFSK){
        if(d_parameters.get_pulse_len()==0){
          d_mod = new Signal_FSK(d_parameters.get_order(),d_parameters.get_sps(),d_parameters.get_mod_idx(),
                                  seed,NULL,0);
        }
        else{
          d_mod = new Signal_FSK(d_parameters.get_order(),d_parameters.get_sps(),d_parameters.get_mod_idx(),
                                  seed,d_parameters.get_pulse_shape_ptr(),d_parameters.get_pulse_len());
        }
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
