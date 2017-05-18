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

#ifndef INCLUDED_SIGNAL_EXCITER_ONE_PASS_GATE_IMPL_H
#define INCLUDED_SIGNAL_EXCITER_ONE_PASS_GATE_IMPL_H

#include <signal_exciter/one_pass_gate.h>
#include <stdio.h>

namespace gr {
  namespace signal_exciter {

    class one_pass_gate_impl : public one_pass_gate
    {
     private:
      float d_samp_rate;
      float d_on_dur;
      float d_off_dur;
      bool d_consume;

      size_t d_on_count;
      size_t d_off_count;
      size_t d_on_counter;
      size_t d_off_counter;

     public:
      one_pass_gate_impl(float samp_rate, float off_duration, float on_duration, bool consume);
      ~one_pass_gate_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_ONE_PASS_GATE_IMPL_H */

