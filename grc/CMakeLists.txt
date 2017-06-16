# Copyright 2011 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
install(FILES
    signal_exciter_random_signal.xml
    signal_exciter_zero_counter.xml
    signal_exciter_cpm_hier.xml
    signal_exciter_periodic_gate.xml
    signal_exciter_random_gate.xml
    signal_exciter_one_pass_gate.xml
    signal_exciter_gmm.xml
    signal_exciter_whiten_and_compress_block.xml DESTINATION share/gnuradio/grc/blocks
)
