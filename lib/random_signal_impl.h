/* -*- c++ -*- */
/*
 * Copyright 2015 <+YOU OR YOUR COMPANY+>.
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

#ifndef INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_IMPL_H
#define INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_IMPL_H

#include <signal_exciter/random_signal.h>
#include <gnuradio/fft/fft.h>//gain access to a shared fftw mutex
#include <signal_exciter/signal_base.hpp>
#include "signal_fm.hpp"
#include "signal_dsb.hpp"
#include "signal_dsbsc.hpp"
#include "signal_usb.hpp"
#include "signal_lsb.hpp"
#include "signal_psk.hpp"
#include "signal_qam.hpp"
#include "signal_pam.hpp"
#include "signal_ask.hpp"
#include "signal_ofdm.hpp"
#include "signal_cwmorse.hpp"
#include "signal_cpm.hpp"
#include "signal_fsk.hpp"

namespace gr {
  namespace signal_exciter {

    class random_signal_impl : public random_signal
    {
     private:
      sig_params d_params;
      Signal_Base* d_mod;
      bool roundone;

     public:
      //random_signal_impl(int seed=-1);
      random_signal_impl(sig_params sig, int seed=-1);
      ~random_signal_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_IMPL_H */
