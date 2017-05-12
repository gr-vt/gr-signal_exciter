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

//boost stuff
#include <boost/thread/mutex.hpp>
#include <boost/random/random_device.hpp>
#include <boost/math/special_functions/sinc.hpp>

//gr stuff
#include <gnuradio/random.h>
#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include <gnuradio/filter/firdes.h>

//eigen stuff
#include <eigen3/Eigen/Dense>


typedef std::complex<float> complexf;
typedef std::complex<double> complexd;

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

    //For thread safe create of fftw objects.
    static boost::mutex s_mutex_fftw;
    static boost::mutex* d_mutex_ptr;

  protected:
    boost::random_device d_rd;
    gr::random *d_rng;

    virtual boost::mutex& fftw_lock() const {return *d_mutex_ptr;}

    void prototype_augemnt_fractional_delay(
          double interp, std::vector<float> &proto,
          double frac_delay, std::vector<float> &taps,
          std::vector<float> &extended_proto);

  public:
    virtual ~Signal_Base() = 0;

    virtual void generate_signal(complexf* output, size_t sample_count) = 0;
    virtual void generate_symbols(complexf* output, size_t symbol_count) = 0;

    static void set_mutex_pointer(boost::mutex* ext_mutex){ d_mutex_ptr = ext_mutex; }
//    int get_seed();
};

inline Signal_Base::~Signal_Base()
{}

#endif /* INCLUDED_SIGNAL_BASE_HPP */
