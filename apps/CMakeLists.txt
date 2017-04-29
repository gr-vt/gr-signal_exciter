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

include(GrPlatform) #define LIB_SUFFIX

include_directories(${CMAKE_INSTALL_PREFIX}/include/)
link_directories(${CMAKE_INSTALL_PREFIX}/lib/)

add_executable(Signal_Exciter_File Signal_Exciter_File.cpp)
target_link_libraries(Signal_Exciter_File
	${GNURADIO_RUNTIME_LIBRARIES} ${GNURADIO_BLOCKS_LIBRARIES}
	${GNURADIO_FILTER_LIBRARIES} ${GNURADIO_DIGITAL_LIBRARIES}
	${GNURADIO_CHANNELS_LIBRARIES}
	gnuradio-signal_exciter ${Boost_LIBRARIES} fftw3 fftw3f config)

include(GrPython)

install(
    TARGETS
    Signal_Exciter_File
    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)
