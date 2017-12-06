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
#include "one_pass_gate_impl.h"
#include <stdexcept>
#include <string.h>
#include <stdio.h>

namespace gr {
  namespace signal_exciter {

    one_pass_gate::sptr
    one_pass_gate::make(float samp_rate, float off_duration, float on_duration, bool consume)
    {
      return gnuradio::get_initial_sptr
        (new one_pass_gate_impl(samp_rate, off_duration, on_duration, consume));
    }

    /*
     * The private constructor
     */
    one_pass_gate_impl::one_pass_gate_impl(float samp_rate, float off_duration, float on_duration, bool consume)
      : gr::block("one_pass_gate",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_samp_rate(samp_rate),
        d_on_dur(on_duration),
        d_off_dur(off_duration),
        d_consume(consume)
    {
      if((d_samp_rate<0)||(d_on_dur<0)||(d_off_dur<0)){
        throw std::runtime_error("One Pass Gate requires all inputs to be positive.\n");
      }
      if(d_samp_rate==0.){
        throw std::runtime_error("One Pass Gate requires sample rate > 0.\n");
      }
      
      d_on_count = int(d_samp_rate * d_on_dur);
      d_off_count = int(d_samp_rate * d_off_dur);

      d_on_counter = 0;
      d_off_counter = 0;

    }

    /*
     * Our virtual destructor.
     */
    one_pass_gate_impl::~one_pass_gate_impl()
    {
    }

    int
    one_pass_gate_impl::general_work(int noutput_items,
        gr_vector_int& ninput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      memset( &out[0], 0, noutput_items*sizeof(gr_complex) );

      if(d_on_counter < d_on_count){
        size_t oo(0),min_cpy(0);
        while((oo < noutput_items)&&(d_on_counter < d_on_count)){
          if(d_off_counter < d_off_count){
            min_cpy = (noutput_items < (d_off_count - d_off_counter)) ? noutput_items : (d_off_count - d_off_counter);
            oo += min_cpy;
            d_off_counter += min_cpy;
          }
          else if(d_on_counter < d_on_count){
            min_cpy = ((noutput_items-oo) < (d_on_count-d_on_counter)) ? (noutput_items-oo) : (d_on_count-d_on_counter);
            memcpy( &out[oo], &in[oo], min_cpy*sizeof(gr_complex) );
            oo += min_cpy;
            d_on_counter += min_cpy;
          }
        }
        if(!d_consume) consume_each (noutput_items);
      }

      // Tell runtime system how many input items we consumed on
      // each input stream.
      if(d_consume) consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace signal_exciter */
} /* namespace gr */

