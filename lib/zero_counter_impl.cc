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
#include "zero_counter_impl.h"
#include <stdexcept>
#include <cmath>
#include <stdio.h>

namespace gr {
  namespace signal_exciter {

    zero_counter::sptr
    zero_counter::make(int item_size)
    {
      return gnuradio::get_initial_sptr
        (new zero_counter_impl(item_size));
    }

    /*
     * The private constructor
     */
    zero_counter_impl::zero_counter_impl(int item_size)
      : gr::sync_block("zero_counter",
              gr::io_signature::make(1, 1, item_size),
              gr::io_signature::make(0, 0, 0)),
        d_faults(0),
        d_itemsize(item_size),
        d_zero(0),
        d_nan(0),
        d_inf(0)
    {
      if(!((d_itemsize==4)||(d_itemsize==8))){
        throw std::runtime_error("zero_counter: Invalid item size\n");
      }
    }

    /*
     * Our virtual destructor.
     */
    zero_counter_impl::~zero_counter_impl()
    {
      printf("Faults detected: %lu\n",d_faults);
      printf("Zeros  detected: %lu\n",d_zero);
      printf("NANs   detected: %lu\n",d_nan);
      printf("Infs   detected: %lu\n",d_inf);
    }

    int
    zero_counter_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      if(d_itemsize == 4){
        const float *in = (const float *) input_items[0];
        for(size_t ii = 0; ii < noutput_items; ii++){
          if(std::isnan(in[ii])) d_nan++;
          if(std::isinf(in[ii])) d_inf++;
          if(in[ii] == 0.)  d_zero++;
        }
      }
      else if(d_itemsize == 8){
        const gr_complex *in = (const gr_complex *) input_items[0];
        for(size_t ii = 0; ii < noutput_items; ii++){
          if(std::isnan(in[ii].real())||std::isnan(in[ii].imag()))  d_nan++;
          if(std::isinf(in[ii].real())||std::isinf(in[ii].imag()))  d_inf++;
          if((in[ii].real() == 0.)&&(in[ii].imag() == 0.))          d_zero++;
        }
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

    int
    zero_counter_impl::get_faults()
    {
      return d_zero+d_nan+d_inf;
    }

    int
    zero_counter_impl::get_nans()
    {
      return d_nan;
    }

    int
    zero_counter_impl::get_infs()
    {
      return d_inf;
    }

    int
    zero_counter_impl::get_zeros()
    {
      return d_zero;
    }

  } /* namespace signal_exciter */
} /* namespace gr */

