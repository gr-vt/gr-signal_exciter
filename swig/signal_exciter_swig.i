/* -*- c++ -*- */

#define SIGNAL_EXCITER_API

%include "gnuradio.i"			// the common stuff
%include "std_vector.i"
//%include "carrays.i"
//
//%array_class(size_t, size_tArray);
//%array_class(float, floatArray);

//# define SWIGPY_SLICE_ARG(obj) ((PySliceObject*) (obj))

//load generated python docstrings
%include "signal_exciter_swig_doc.i"

%{
#define SWIGPY_SLICE_ARG(obj) ((PySliceObject*) (obj))
#include "signal_exciter/random_signal_config.h"
#include "signal_exciter/random_signal.h"
#include "signal_exciter/zero_counter.h"
#include "signal_exciter/cpm_hier.h"
#include "signal_exciter/periodic_gate.h"
#include "signal_exciter/random_gate.h"
#include "signal_exciter/one_pass_gate.h"
#include "signal_exciter/gmm.h"
#include "signal_exciter/whiten_and_compress_block.h"
%}

%rename(__assign__) *::operator=;

namespace std{
  %template(FloatVector) vector<float>;
  %template(SizeVector) vector<size_t>;
}


%include "signal_exciter/random_signal_config.h"
%include "signal_exciter/random_signal.h"
GR_SWIG_BLOCK_MAGIC2(signal_exciter, random_signal);
%include "signal_exciter/zero_counter.h"
GR_SWIG_BLOCK_MAGIC2(signal_exciter, zero_counter);
%include "signal_exciter/cpm_hier.h"
GR_SWIG_BLOCK_MAGIC2(signal_exciter, cpm_hier);
%include "signal_exciter/periodic_gate.h"
GR_SWIG_BLOCK_MAGIC2(signal_exciter, periodic_gate);
%include "signal_exciter/random_gate.h"
GR_SWIG_BLOCK_MAGIC2(signal_exciter, random_gate);
%include "signal_exciter/one_pass_gate.h"
GR_SWIG_BLOCK_MAGIC2(signal_exciter, one_pass_gate);
%include "signal_exciter/gmm.h"
GR_SWIG_BLOCK_MAGIC2(signal_exciter, gmm);
%include "signal_exciter/whiten_and_compress_block.h"
GR_SWIG_BLOCK_MAGIC2(signal_exciter, whiten_and_compress_block);

