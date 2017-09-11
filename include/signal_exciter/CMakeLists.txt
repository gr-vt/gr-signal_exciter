# Copyright 2011,2012 Free Software Foundation, Inc.
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

########################################################################
# Install public header files
########################################################################
install(FILES
    api.h
    signal_base.hpp
    random_signal_config.h
    random_signal.h
    signal_threaded_buffer.h
    gaussian_mixture.h
    zero_counter.h
    cpm_hier.h
    periodic_gate.h
    random_gate.h
    one_pass_gate.h
    gmm.h
    gmm_spectral_taps.h
    whiten_and_compress_block.h DESTINATION include/signal_exciter
)
