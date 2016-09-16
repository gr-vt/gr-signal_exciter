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
    cpm_hier::make(sig_params sig, bool gray)
    {
      return gnuradio::get_initial_sptr
        (new cpm_hier_impl(sig, gray));
    }

    /*
     * The private constructor
     */
    cpm_hier_impl::cpm_hier_impl(sig_params sig, bool gray)
      : gr::hier_block2("cpm_hier",
              gr::io_signature::make(1, 1, sizeof(unsigned char)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_gray(gray)
    {
      //printf("CPM HIER STARTING.\n");
      d_order = sig.order;
      if(d_gray){
        int order = sig.order;
        int new_order=0;
        while(order>0){
          order = order>>1;
          if(!new_order) new_order = 1;
          else new_order = new_order*2;
        }
        if((new_order != d_order) || (d_order == 0)){
          if(new_order > 0) d_order = new_order;
          else d_order = 2;
        }
      }
      if(d_order <= 1) d_order = 2;

      create_symbol_list();

      //d_src = gr::analog::random_uniform_source_b::make(0,256,d_seed);
      //d_upk = gr::blocks::packed_to_unpacked_bb::make(1,gr::GR_MSB_FIRST);
      d_sym = gr::digital::chunks_to_symbols_bf::make(d_symbols,1);

      switch(sig.type){
        case GFSK :{
          d_phase_type = gr::analog::cpm::GAUSSIAN;
          d_sym_overlap = sig.L;
          d_mod_idx = sig.mod_idx;
          d_beta = sig.beta;

          std::vector<float> gt = gr::filter::firdes::gaussian(1,sig.sps,d_beta,d_sym_overlap*sig.sps);

          std::vector<float> rt(sig.sps,1.);
          d_taps = std::vector<float>(gt.size()+rt.size()-1);
          for(size_t idx = 0; idx < d_taps.size(); idx++){
            size_t count = std::min(std::min(idx,rt.size()),gt.size());

            float summer = 0.;
            for(size_t dp = 0; dp < count; dp++){
              if(idx-dp < gt.size())
                summer += gt[idx-dp]*rt[dp];
            }
            d_taps[idx] = summer;
          }
          break;
        }
        case CPM :{
          d_phase_type = gr::analog::cpm::cpm_type(sig.phase_type);
          d_sym_overlap = sig.L;
          d_mod_idx = sig.mod_idx;
          d_beta = sig.beta;

          d_taps = gr::analog::cpm::phase_response(d_phase_type, sig.sps, d_sym_overlap, d_beta);
          break;
        }
        case MSK :{
          d_phase_type = gr::analog::cpm::LREC;
          d_sym_overlap = 1;
          d_mod_idx = 0.5;
          d_beta = 0;

          d_taps = gr::analog::cpm::phase_response(d_phase_type, sig.sps, d_sym_overlap, d_beta);
          break;
        }
        case GMSK :{
          d_phase_type = gr::analog::cpm::GAUSSIAN;
          d_sym_overlap = sig.L;
          d_mod_idx = 0.5;
          d_beta = sig.beta;

          d_taps = gr::analog::cpm::phase_response(d_phase_type, sig.sps, d_sym_overlap, d_beta);
          break;
        }
        case FSK :{
          d_phase_type = gr::analog::cpm::LREC;
          d_sym_overlap = 1;
          d_mod_idx = sig.mod_idx;
          d_beta = 0;

          d_taps = gr::analog::cpm::phase_response(d_phase_type, sig.sps, d_sym_overlap, d_beta);
          break;
        }
        default :{
          throw std::runtime_error("Unknown type with CPM_HIER\n");
          break;
        }
      }

      //printf("sps = %d;taps = (%1.14e",sig.sps,d_taps[0]);
      //for(size_t idx = 1; idx < d_taps.size(); idx++){
      //  printf(", %1.14e",d_taps[idx]);
      //}
      //printf(")\n");

      d_fir = gr::filter::interp_fir_filter_fff::make(sig.sps, d_taps);
      d_fm = gr::analog::frequency_modulator_fc::make(M_PI * d_mod_idx);


      //connect(d_src, 0, d_upk, 0);
      //connect(d_upk, 0, d_sym, 0);
      //connect(d_sym, 0,self(), 0);
      connect(self(), 0, d_sym, 0);
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

    void
    cpm_hier_impl::create_symbol_list()
    {
      //gray coding pulled from gr-digital/python/digital/utils/gray_code.py

      //printf("Making symbols.\n");
      d_symbols = std::vector<float>(d_order,0.);
      if(d_gray){
        std::vector<int> sym_base(2,0);
        sym_base[1] = 1;
        //printf("sym_base = [%d,%d];\n",sym_base[0],sym_base[1]);
        int lp2 = 2;
        int np2 = 4;
        int ind = 2;
        int result;
        while(sym_base.size() < d_order){
          if(ind == lp2)
            result = ind + ind/2;
          else
            result = sym_base[2*lp2-1-ind] + lp2;
          sym_base.push_back(result);
          ind++;
          if(ind == np2){
            lp2 = ind;
            np2 = ind*2;
          }
        }
        for(int idx = 0; idx < d_order; idx++){
        //  printf("%d ", sym_base[idx]*2-int(d_order)+1);
          d_symbols[idx] = float(sym_base[idx]*2-int(d_order)+1);
        }
        //printf("\n");
      }
      else{
        for(int idx = 0; idx < d_order; idx++){
          d_symbols[idx] = float(idx*2-int(d_order)+1);
        }
      }

      //printf("Symbols: [%1.0f",d_symbols[0]);
      //for(int idx = 1; idx < d_order; idx++){
      //  printf(", %1.0f", d_symbols[idx]);
      //}
      //printf("]\n");
    }


  } /* namespace signal_exciter */
} /* namespace gr */

