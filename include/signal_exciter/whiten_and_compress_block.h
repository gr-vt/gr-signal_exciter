/* -*- c++ -*- */
/*
 *     gr-signal_exciter  Copyright (C) 2015-2017  Bill Clark
 *     This program comes with ABSOLUTELY NO WARRANTY;
 *     This is free software, and you are welcome to redistribute it
 *     under certain conditions;
 */


#ifndef INCLUDED_SIGNAL_EXCITER_WHITEN_AND_COMPRESS_BLOCK_H
#define INCLUDED_SIGNAL_EXCITER_WHITEN_AND_COMPRESS_BLOCK_H

#include <signal_exciter/api.h>
#include <gnuradio/sync_interpolator.h>
#include <complex>

namespace gr {
  namespace signal_exciter {

    /*!
     * \brief Whiten and Compress: Take a clean signal, and combine with AWGN
     * then compress from complex float to complex int16
     * \ingroup signal_exciter
     *
     */
    class SIGNAL_EXCITER_API whiten_and_compress_block : virtual public gr::sync_interpolator
    {
     public:
      typedef boost::shared_ptr<whiten_and_compress_block> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of signal_exciter::whiten_and_compress_block.
       *
       * To avoid accidental use of raw pointers, signal_exciter::whiten_and_compress_block's
       * constructor is in a private implementation
       * class. signal_exciter::whiten_and_compress_block::make is the public interface for
       * creating new instances.
       */
      static sptr make(size_t blocksize, double snr_db, int seed=-1);
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_WHITEN_AND_COMPRESS_BLOCK_H */

