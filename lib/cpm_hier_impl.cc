/* -*- c++ -*- */
/* 
 * Copyright 2016 <+YOU OR YOUR COMPANY+>.
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
#include "cpm_hier_impl.h"
#include <gnuradio/analog/cpm.h>
#include <stdexcept>
#include <stdio.h>
#include <math.h>
#include <algorithm>

namespace gr {
  namespace signal_exciter {

    cpm_hier::sptr
    cpm_hier::make(sig_params sig, int seed)
    {
      return gnuradio::get_initial_sptr
        (new cpm_hier_impl(sig, seed));
    }

    /*
     * The private constructor
     */
    cpm_hier_impl::cpm_hier_impl(sig_params sig, int seed)
      : gr::hier_block2("cpm_hier",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(1, 1, sizeof(gr_complex)))
    {
      //printf("CPM HIER STARTING.\n");
      float syms[2] = {-1.,1.};
      d_symbols = std::vector<float>(syms,syms+sizeof(syms)/sizeof(float));

      if(seed<0) d_seed = (time(NULL));
      else d_seed = seed;

      d_src = gr::analog::random_uniform_source_b::make(0,256,d_seed);
      d_upk = gr::blocks::packed_to_unpacked_bb::make(1,gr::GR_MSB_FIRST);
      d_sym = gr::digital::chunks_to_symbols_bf::make(d_symbols,1);

      switch(sig.type){
        case GFSK :{
          d_phase_type = gr::analog::cpm::GAUSSIAN;
          d_sym_overlap = sig.cpm_params.L;
          d_sensitivity = sig.cpm_params.fm_sense;
          d_beta = sig.cpm_params.beta;

          std::vector<float> gt = gr::filter::firdes::gaussian(1,sig.sps,d_beta,d_sym_overlap*sig.sps);
          std::vector<float> rt(sig.sps,1.);
          d_taps = std::vector<float>(gt.size()+rt.size()-1);
          for(size_t idx = 0; idx < d_taps.size(); idx++){
            size_t count = std::min(std::min(idx,rt.size()),gt.size());

            float summer = 0.;
            for(size_t dp = 0; dp <= count; dp++){
              if(idx-dp < gt.size())
                summer += gt[idx-dp]*rt[dp];
            }
            d_taps[idx] = summer;
          }
          //d_taps = gr::filter::firdes::gaussian(1,sig.sps,d_beta,d_sym_overlap*sig.sps);
          //printf("(%1.14e",d_taps[0]);
          //for(size_t idx = 1; idx < d_taps.size(); idx++){
          //  printf(", %1.14e",d_taps[idx]);
          //}
          //printf(")\n");
          break;
        }
        case CPM :{
          d_phase_type = gr::analog::cpm::cpm_type(sig.cpm_params.phase_type);
          d_sym_overlap = sig.cpm_params.L;
          d_sensitivity = sig.cpm_params.fm_sense;
          d_beta = sig.cpm_params.beta;

          d_taps = gr::analog::cpm::phase_response(d_phase_type, sig.sps, d_sym_overlap, d_beta);
          break;
        }
        case MSK :{
          d_phase_type = gr::analog::cpm::LREC;
          d_sym_overlap = 1;
          d_sensitivity = 0.5;
          d_beta = 0;

          d_taps = gr::analog::cpm::phase_response(d_phase_type, sig.sps, d_sym_overlap, d_beta);
          break;
        }
        case GMSK :{
          d_phase_type = gr::analog::cpm::GAUSSIAN;
          d_sym_overlap = sig.cpm_params.L;
          d_sensitivity = M_PI/2.;
          d_beta = sig.cpm_params.beta;

          d_taps = gr::analog::cpm::phase_response(d_phase_type, sig.sps, d_sym_overlap, d_beta);
          break;
        }
        case FSK :{
          d_phase_type = gr::analog::cpm::LREC;
          d_sym_overlap = 1;
          d_sensitivity = sig.cpm_params.fm_sense;
          d_beta = 0;

          //d_taps = std::vector<float>(sig.sps,1.);
          d_taps = gr::analog::cpm::phase_response(d_phase_type, sig.sps, d_sym_overlap, d_beta);
          //printf("(%1.14e",d_taps[0]);
          //for(size_t idx = 1; idx < d_taps.size(); idx++){
          //  printf(", %1.14e",d_taps[idx]);
          //}
          //printf(")\n");
          break;
        }
        default :{
          throw std::runtime_error("Unknown type with CPM_HIER\n");
          break;
        }
      }
      d_fir = gr::filter::interp_fir_filter_fff::make(sig.sps, d_taps);
      d_fm = gr::analog::frequency_modulator_fc::make(d_sensitivity);


      connect(d_src, 0, d_upk, 0);
      connect(d_upk, 0, d_sym, 0);
      //connect(d_sym, 0,self(), 0);
      connect(d_sym, 0, d_fir, 0);
      connect(d_fir, 0, d_fm , 0);
      connect(d_fm , 0, self(), 0);
      
    }

    /*
     * Our virtual destructor.
     */
    cpm_hier_impl::~cpm_hier_impl()
    {
    }


  } /* namespace signal_exciter */
} /* namespace gr */

