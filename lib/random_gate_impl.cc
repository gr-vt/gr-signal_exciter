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
#include "random_gate_impl.h"
#include <stdio.h>

namespace gr {
  namespace signal_exciter {

    random_gate::sptr
    random_gate::make(float samp_rate, float min_off_dur, float max_off_dur, float min_on_dur, float max_on_dur, int seed)
    {
      return gnuradio::get_initial_sptr
        (new random_gate_impl(samp_rate, min_off_dur, max_off_dur, min_on_dur, max_on_dur, seed));
    }

    /*
     * The private constructor
     */
    random_gate_impl::random_gate_impl(float samp_rate, float min_off_dur, float max_off_dur, float min_on_dur, float max_on_dur, int seed)
      : gr::sync_block("random_gate",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
        d_samp_rate(samp_rate),
        d_min_off_dur(min_off_dur),
        d_max_off_dur(max_off_dur),
        d_min_on_dur(min_on_dur),
        d_max_on_dur(max_on_dur)
    {
      if((d_samp_rate<0)||(d_min_on_dur<0)||(d_max_on_dur<0)||(d_min_off_dur<0)||(d_max_off_dur<0)){
        throw std::runtime_error("Random Gate requires all inputs to be positive, sans seed.\n");
      }
      if(d_samp_rate==0.){
        throw std::runtime_error("Random Gate requires sample rate > 0.\n");
      }

      set_seed(seed);

      d_rng = new gr::random(d_seed,0,1337);

      d_on_counter = 0;
      d_off_counter = 0;
      d_cycle_counter = 0;

      d_on_count = rand_on();
      d_off_count = rand_off();
      d_cycle_count = d_on_count + d_off_count;


      //printf("Inital settings: OFF(%lu) ON(%lu)\n",d_off_count, d_on_count);
    }

    /*
     * Our virtual destructor.
     */
    random_gate_impl::~random_gate_impl()
    {
      delete d_rng;
    }

    int
    random_gate_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];

      memset( &out[0], 0, noutput_items*sizeof(gr_complex) );

      size_t oo(0),min_cpy(0);

      while(oo < noutput_items){
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

          d_on_count = rand_on();
          d_off_count = rand_off();
          d_cycle_count = d_on_count + d_off_count;

          //printf("Next settings: OFF(%lu) ON(%lu)\n",d_off_count, d_on_count);
        }
      }


      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

    size_t
    random_gate_impl::rand_on()
    {
      //float u = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
      float u = d_rng->ran1();
      return int((u*(d_max_on_dur-d_min_on_dur) + d_min_on_dur)*d_samp_rate);
    }

    size_t
    random_gate_impl::rand_off()
    {
      //float u = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
      float u = d_rng->ran1();
      return int((u*(d_max_off_dur-d_min_off_dur) + d_min_off_dur)*d_samp_rate);
    }

  } /* namespace signal_exciter */
} /* namespace gr */

