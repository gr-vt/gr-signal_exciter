/* -*- c++ -*- */
/* 
 * Copyright 2015 <+YOU OR YOUR COMPANY+>.
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


#ifndef INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_H
#define INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_H

#include <signal_exciter/api.h>
#include <gnuradio/sync_block.h>
#include <signal_exciter/random_signal_config.h>



namespace gr {
  namespace signal_exciter {

    /*!
     * \brief <+description of block+>
     * \ingroup signal_exciter
     *
     */
    class SIGNAL_EXCITER_API random_signal : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<random_signal> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of signal_exciter::random_signal.
       *
       * To avoid accidental use of raw pointers, signal_exciter::random_signal's
       * constructor is in a private implementation
       * class. signal_exciter::random_signal::make is the public interface for
       * creating new instances.
       */
      //static sptr make(int seed=-1);
      static sptr make(sig_params sig, int seed=-1);
      static sptr make(const signal_parameters &sig, int seed=-1);
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_RANDOM_SIGNAL_H */

