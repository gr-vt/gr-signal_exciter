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

#ifndef INCLUDED_SIGNAL_EXCITER_ZERO_COUNTER_IMPL_H
#define INCLUDED_SIGNAL_EXCITER_ZERO_COUNTER_IMPL_H

#include <signal_exciter/zero_counter.h>

namespace gr {
  namespace signal_exciter {

    class zero_counter_impl : public zero_counter
    {
     private:
      int d_itemsize;
      size_t d_faults;
      size_t d_nan;
      size_t d_zero;
      size_t d_inf;

     public:
      zero_counter_impl(int item_size);
      ~zero_counter_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);

      int get_faults();
      int get_nans();
      int get_infs();
      int get_zeros();
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_ZERO_COUNTER_IMPL_H */

