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
# Setup library
########################################################################
include(GrPlatform) #define LIB_SUFFIX

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
list(APPEND signal_exciter_sources
    gmm_impl.cc
    signal_base.cpp
    signal_threaded_buffer.cpp
    signal_cwmorse.cpp
    signal_fm.cpp
    signal_dsb.cpp
    signal_dsbsc.cpp
    signal_usb.cpp
    signal_lsb.cpp
    signal_psk.cpp
    signal_qam.cpp
    signal_pam.cpp
    signal_cpm.cpp
    signal_ofdm.cpp
    random_signal_impl.cc
    zero_counter_impl.cc
    cpm_hier_impl.cc
    periodic_gate_impl.cc
    random_gate_impl.cc
    one_pass_gate_impl.cc )

set(signal_exciter_sources "${signal_exciter_sources}" PARENT_SCOPE)
if(NOT signal_exciter_sources)
	MESSAGE(STATUS "No C++ sources... skipping lib/")
	return()
endif(NOT signal_exciter_sources)

add_library(gnuradio-signal_exciter SHARED ${signal_exciter_sources})
target_link_libraries(gnuradio-signal_exciter ${Boost_LIBRARIES} ${GNURADIO_ALL_LIBRARIES}
	  fftw3f)
set_target_properties(gnuradio-signal_exciter PROPERTIES DEFINE_SYMBOL "gnuradio_signal_exciter_EXPORTS")
add_definitions(-DEIGEN_DONT_PARALLELIZE)

if(APPLE)
    set_target_properties(gnuradio-signal_exciter PROPERTIES
        INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    )
endif(APPLE)

########################################################################
# Install built library files
########################################################################
install(TARGETS gnuradio-signal_exciter
    LIBRARY DESTINATION lib${LIB_SUFFIX} # .so/.dylib file
    ARCHIVE DESTINATION lib${LIB_SUFFIX} # .lib file
    RUNTIME DESTINATION bin              # .dll file
)

########################################################################
# Build and register unit test
########################################################################
include(GrTest)

include_directories(${CPPUNIT_INCLUDE_DIRS})

list(APPEND test_signal_exciter_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/test_signal_exciter.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_signal_exciter.cc
)

add_executable(test-signal_exciter ${test_signal_exciter_sources})

target_link_libraries(
  test-signal_exciter
  ${GNURADIO_RUNTIME_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CPPUNIT_LIBRARIES}
  gnuradio-signal_exciter
)

GR_ADD_TEST(test_signal_exciter test-signal_exciter)
