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


#ifndef INCLUDED_SIGNAL_EXCITER_ZERO_COUNTER_H
#define INCLUDED_SIGNAL_EXCITER_ZERO_COUNTER_H

#include <signal_exciter/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace signal_exciter {

    /*!
     * \brief <+description of block+>
     * \ingroup signal_exciter
     *
     */
    class SIGNAL_EXCITER_API zero_counter : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<zero_counter> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of signal_exciter::zero_counter.
       *
       * To avoid accidental use of raw pointers, signal_exciter::zero_counter's
       * constructor is in a private implementation
       * class. signal_exciter::zero_counter::make is the public interface for
       * creating new instances.
       */
      static sptr make(int item_size);

      virtual int get_faults() = 0;
      virtual int get_nans() = 0;
      virtual int get_infs() = 0;
      virtual int get_zeros() = 0;
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_ZERO_COUNTER_H */

