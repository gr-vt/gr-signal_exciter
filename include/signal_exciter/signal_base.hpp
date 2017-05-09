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
#include <boost/random/random_device.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <gnuradio/fft/window.h>


typedef std::complex<float> complexf;

#ifndef M_2_PIl
#define M_2_PIl (2*M_PIl)
#endif

#define ENABLE_FRAC_DELAY_IN_BLOCK false

class Signal_Base
{
  private:
//    static bool d_seeded;
//    static int d_seed;

    virtual void auto_fill_symbols() = 0;
    virtual void auto_fill_signal() = 0;

  protected:
    boost::random_device d_rd;
    gr::random *d_rng;
    //For thread safe create of fftw objects.
    static boost::mutex s_mutex_fftw;

    virtual boost::mutex& fftw_lock() const {return s_mutex_fftw;}

//    void set_seed(int seed=-1);
    virtual void time_offset(std::vector<float> &taps,
                             std::vector<float> &proto,
                             float offset);

  public:
    virtual ~Signal_Base() = 0;

    virtual void generate_signal(complexf* output, size_t sample_count) = 0;
    virtual void generate_symbols(complexf* output, size_t symbol_count) = 0;

//    int get_seed();
};

inline Signal_Base::~Signal_Base()
{}

inline void Signal_Base::time_offset(std::vector<float> &taps,
                                     std::vector<float> &proto,
                                     float offset)
{
  if(ENABLE_FRAC_DELAY_IN_BLOCK){
    std::vector<double> proto2(proto.begin(),proto.end());
    std::vector<double> taps2(proto.size(),0.);
    std::vector<float> window = gr::fft::window::blackman_harris(proto.size());
    std::vector<double> window2(window.begin(), window.end());

    for(int idx = 0; idx < proto.size(); idx++){
      for(int ind = 0; ind < proto.size(); ind++){
        taps2[idx] += proto2[ind]*window2[ind]*
            boost::math::sinc_pi(M_PI*(double(idx)-double(ind) - double(offset)));
      }
      taps2[idx] *= window2[idx];
    }
    taps = std::vector<float>(taps2.begin(),taps2.end());
  }
  else{
    taps = proto;
  }
}

#endif /* INCLUDED_SIGNAL_BASE_HPP */
