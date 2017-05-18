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

#ifndef INCLUDED_SIGNAL_EXCITER_RANDOM_GATE_IMPL_H
#define INCLUDED_SIGNAL_EXCITER_RANDOM_GATE_IMPL_H

#include <signal_exciter/random_gate.h>
#include <gnuradio/random.h>
#include <boost/random/random_device.hpp>
#include <stdio.h>

namespace gr {
  namespace signal_exciter {

    class random_gate_impl : public random_gate
    {
     private:
      float d_samp_rate;
      float d_min_on_dur;
      float d_max_on_dur;
      float d_min_off_dur;
      float d_max_off_dur;

      int d_seed;
      size_t d_on_count;
      size_t d_off_count;
      size_t d_on_counter;
      size_t d_off_counter;

      size_t d_cycle_count;
      size_t d_cycle_counter;

      boost::random_device d_rd;
      gr::random* d_rng;

      size_t rand_on();
      size_t rand_off();

     public:
      random_gate_impl(float samp_rate, float min_off_dur, float max_off_dur, float min_on_dur, float max_on_dur, int seed);
      ~random_gate_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);

      void set_seed(int seed)
      {
        d_seed = seed;
        if(d_seed < 0) d_seed = d_rd();
      }
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_RANDOM_GATE_IMPL_H */

