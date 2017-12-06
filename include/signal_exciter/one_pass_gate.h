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


#ifndef INCLUDED_SIGNAL_EXCITER_ONE_PASS_GATE_H
#define INCLUDED_SIGNAL_EXCITER_ONE_PASS_GATE_H

#include <signal_exciter/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace signal_exciter {

    /*!
     * \brief <+description of block+>
     * \ingroup signal_exciter
     *
     */
    class SIGNAL_EXCITER_API one_pass_gate : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<one_pass_gate> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of signal_exciter::one_pass_gate.
       *
       * To avoid accidental use of raw pointers, signal_exciter::one_pass_gate's
       * constructor is in a private implementation
       * class. signal_exciter::one_pass_gate::make is the public interface for
       * creating new instances.
       */
      static sptr make(float samp_rate, float off_duration, float on_duration, bool consume);
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_ONE_PASS_GATE_H */

