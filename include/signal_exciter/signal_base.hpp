#ifndef INCLUDED_SIGNAL_BASE_HPP
#define INCLUDED_SIGNAL_BASE_HPP


#include <vector>
#include <complex>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <signal_exciter/signal_threaded_buffer.h>
#include <boost/thread/mutex.hpp>
#include <gnuradio/random.h>


typedef std::complex<float> complexf;

#ifndef M_2_PIl
#define M_2_PIl (2*M_PIl)
#endif


class Signal_Base
{
  private:
//    static bool d_seeded;
//    static int d_seed;

    virtual void auto_fill_symbols() = 0;
    virtual void auto_fill_signal() = 0;

  protected:
    gr::random *d_rng;
    //For thread safe create of fftw objects.
    static boost::mutex s_mutex_fftw;

    virtual boost::mutex& fftw_lock() const {return s_mutex_fftw;}

//    void set_seed(int seed=-1);

  public:
    virtual ~Signal_Base() = 0;

    virtual void generate_signal(complexf* output, size_t sample_count) = 0;
    virtual void generate_symbols(complexf* output, size_t symbol_count) = 0;

//    int get_seed();
};

inline Signal_Base::~Signal_Base()
{}

#endif /* INCLUDED_SIGNAL_BASE_HPP */
