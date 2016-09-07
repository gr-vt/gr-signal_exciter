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
#include "gmm_impl.h"

namespace gr {
  namespace signal_exciter {

    gmm::sptr
    gmm::make(float sigma1, float sigma2, float fmax, float thresh, int seed)
    {
      return gnuradio::get_initial_sptr
        (new gmm_impl(sigma1, sigma2, fmax, thresh, seed));
    }

    /*
     * The private constructor
     */
    gmm_impl::gmm_impl(float sigma1, float sigma2, float fmax, float thresh, int seed)
      : gr::sync_block("gmm",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(1, 1, sizeof(float))),
        d_sig1(sigma1),
        d_sig2(sigma2),
        d_trsh(thresh),
        d_fmax(fmax),
        d_seed(seed),
        d_gm(d_sig1,d_sig2,d_fmax,d_trsh,d_seed)
    {}

    /*
     * Our virtual destructor.
     */
    gmm_impl::~gmm_impl()
    {
    }

    int
    gmm_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      float *out = (float *) output_items[0];

      d_gm.get_message( &out[0], noutput_items );

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace signal_exciter */
} /* namespace gr */

