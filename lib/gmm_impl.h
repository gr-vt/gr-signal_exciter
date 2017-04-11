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

#ifndef INCLUDED_SIGNAL_EXCITER_GMM_IMPL_H
#define INCLUDED_SIGNAL_EXCITER_GMM_IMPL_H

#include <signal_exciter/gmm.h>
#include <signal_exciter/gaussian_mixture.h>

namespace gr {
  namespace signal_exciter {

    class gmm_impl : public gmm
    {
     private:
      float d_sig1;
      float d_sig2;
      float d_fmax;
      float d_trsh;
      int   d_seed;
      Gaussian_Mixture d_gm;

     public:
      gmm_impl(float sigma1, float sigma2, float fmax, float thresh, int seed);
      ~gmm_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);

    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_GMM_IMPL_H */

