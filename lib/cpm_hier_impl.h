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

#ifndef INCLUDED_SIGNAL_EXCITER_CPM_HIER_IMPL_H
#define INCLUDED_SIGNAL_EXCITER_CPM_HIER_IMPL_H

#include <signal_exciter/cpm_hier.h>
#include <gnuradio/analog/cpm.h>
#include <gnuradio/analog/random_uniform_source_b.h>
#include <gnuradio/filter/interp_fir_filter_fff.h>
#include <gnuradio/filter/firdes.h>
#include <gnuradio/analog/frequency_modulator_fc.h>
#include <gnuradio/blocks/packed_to_unpacked_bb.h>
#include <gnuradio/digital/chunks_to_symbols_bf.h>

#include <time.h>
#include <signal_exciter/signal_base.hpp>

namespace gr {
  namespace signal_exciter {

    class cpm_hier_impl : public cpm_hier
    {
     private:
      gr::analog::cpm::cpm_type d_phase_type;
      size_t d_order;
      size_t d_sym_overlap;
      double d_mod_idx;
      double d_beta;
      //size_t d_seed;
      bool d_gray;
      bool d_conv;

      std::vector<float> d_symbols;
      std::vector<float> d_taps;
      
      gr::analog::random_uniform_source_b::sptr d_src;
      gr::analog::frequency_modulator_fc::sptr d_fm;
      gr::filter::interp_fir_filter_fff::sptr d_fir;
      gr::blocks::packed_to_unpacked_bb::sptr d_upk;
      gr::digital::chunks_to_symbols_bf::sptr d_sym;

      void create_symbol_list();

     public:
      cpm_hier_impl(sig_params sig, bool gray=false);
      ~cpm_hier_impl();

    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_CPM_HIER_IMPL_H */

