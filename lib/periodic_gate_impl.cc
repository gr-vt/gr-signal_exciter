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
#include "periodic_gate_impl.h"

namespace gr {
  namespace signal_exciter {

    periodic_gate::sptr
    periodic_gate::make(float samp_rate, float off_duration, float on_duration, float period_offset)
    {
      return gnuradio::get_initial_sptr
        (new periodic_gate_impl(samp_rate, off_duration, on_duration, period_offset));
    }

    /*
     * The private constructor
     */
    periodic_gate_impl::periodic_gate_impl(float samp_rate, float off_duration, float on_duration, float period_offset)
      : gr::sync_block("periodic_gate",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_samp_rate(samp_rate),
        d_on_dur(on_duration),
        d_off_dur(off_duration),
        d_period_offset(period_offset)
    {
      if((d_samp_rate<0)||(d_on_dur<0)||(d_off_dur<0)||(d_period_offset<0)){
        throw std::runtime_error("Periodic Gate requires all inputs to be positive.\n");
      }
      if(d_samp_rate==0.){
        throw std::runtime_error("Periodic Gate requires sample rate > 0.\n");
      }

      d_on_counter = 0;
      d_off_counter = 0;

      d_period = d_on_dur+d_off_dur;
      
      d_on_count = int(d_samp_rate * d_on_dur);
      d_off_count = int(d_samp_rate * d_off_dur);
      d_cycle_count = int(d_samp_rate * d_period);

      if(d_cycle_count != d_on_count+d_off_count){
        //throw std::runtime_error("Periodic Gate: This isn't a true error, but I haven't decided how to actually handle this, or if it should be handled at all.\n");
        printf("Difference found in the periodic gate: Found as %lu, but going to be %lu\n",d_cycle_count,(d_on_count+d_off_count));
        d_cycle_count = d_on_count+d_off_count;
      }

      d_cycle_counter = int(d_samp_rate * fmod(d_period_offset,d_period));
      if(d_cycle_counter > d_off_count){
        d_off_counter = d_off_count;
        d_on_counter = d_cycle_counter - d_off_counter;
      }
      else{
        d_off_counter = d_cycle_counter;
      }
    }

    /*
     * Our virtual destructor.
     */
    periodic_gate_impl::~periodic_gate_impl()
    {
    }

    int
    periodic_gate_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      memset( &out[0], 0, noutput_items*sizeof(gr_complex) );

      size_t oo(0),min_cpy(0);

      while(oo < noutput_items){
        //printf("noutput_items(%d), oo(%lu), d_off_counter(%lu), d_on_counter(%lu), d_cycle_counter(%lu), d_off_count(%lu), d_on_count(%lu), d_cycle_count(%lu)\n",noutput_items,oo, d_off_counter, d_on_counter, d_cycle_counter, d_off_count, d_on_count, d_cycle_count);
        if(d_cycle_counter < d_off_count){
          min_cpy = ((noutput_items-oo) < (d_off_count-d_off_counter)) ? (noutput_items-oo) : (d_off_count-d_off_counter);
          oo += min_cpy;
          d_off_counter += min_cpy;
          d_cycle_counter += min_cpy;
        }
        else if(d_cycle_counter < d_cycle_count){
          min_cpy = ((noutput_items-oo) < (d_on_count-d_on_counter)) ? (noutput_items-oo) : (d_on_count-d_on_counter);
          memcpy( &out[oo], &in[oo], min_cpy*sizeof(gr_complex) );
          oo += min_cpy;
          d_on_counter += min_cpy;
          d_cycle_counter += min_cpy;
        }
        if(d_cycle_counter == d_cycle_count){
          d_on_counter = 0;
          d_off_counter = 0;
          d_cycle_counter = 0;
        }
      }


      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace signal_exciter */
} /* namespace gr */

