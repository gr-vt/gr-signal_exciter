/* -*- c++ -*- */
/*
 *     gr-signal_exciter  Copyright (C) 2015-2017  Bill Clark
 *     This program comes with ABSOLUTELY NO WARRANTY;
 *     This is free software, and you are welcome to redistribute it
 *     under certain conditions;
 */

#ifndef INCLUDED_SIGNAL_EXCITER_WHITEN_AND_COMPRESS_BLOCK_IMPL_H
#define INCLUDED_SIGNAL_EXCITER_WHITEN_AND_COMPRESS_BLOCK_IMPL_H

#include <signal_exciter/whiten_and_compress_block.h>
#include <math.h>
#include <volk/volk.h>
#include <gnuradio/random.h>
#include <boost/random/random_device.hpp>

namespace gr {
  namespace signal_exciter {

    typedef std::complex<float> complexf;

    class whiten_and_compress_block_impl : public whiten_and_compress_block
    {
     private:
      int d_seed;
      size_t d_block_size;
      double d_snr_db;
      double d_noise_var;
      double d_scale_base;

      gr::random *d_rng;
      boost::random_device *d_rd;

      complexf* d_sig;
      complexf* d_nse;

     public:
      whiten_and_compress_block_impl(size_t blocksize, double snr_db, int seed);
      ~whiten_and_compress_block_impl();

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);

      void set_seed(int seed)
      {
        d_seed = seed;
        if(d_seed < 0){
          d_rd = new boost::random_device();
          d_seed = (*d_rd)();
          delete d_rd;
        }
      }
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_WHITEN_AND_COMPRESS_BLOCK_IMPL_H */

