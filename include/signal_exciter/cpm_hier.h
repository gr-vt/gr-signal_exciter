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


#ifndef INCLUDED_SIGNAL_EXCITER_CPM_HIER_H
#define INCLUDED_SIGNAL_EXCITER_CPM_HIER_H

#include <signal_exciter/api.h>
#include <gnuradio/hier_block2.h>
#include <signal_exciter/random_signal_config.h>

namespace gr {
  namespace signal_exciter {

    /*!
     * \brief <+description of block+>
     * \ingroup signal_exciter
     *
     */
    class SIGNAL_EXCITER_API cpm_hier : virtual public gr::hier_block2
    {
     public:
      typedef boost::shared_ptr<cpm_hier> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of signal_exciter::cpm_hier.
       *
       * To avoid accidental use of raw pointers, signal_exciter::cpm_hier's
       * constructor is in a private implementation
       * class. signal_exciter::cpm_hier::make is the public interface for
       * creating new instances.
       */
      static sptr make(sig_params sig, bool gray=false);
    };

  } // namespace signal_exciter
} // namespace gr

#endif /* INCLUDED_SIGNAL_EXCITER_CPM_HIER_H */

