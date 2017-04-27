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
  if((offset > - 1.19e-07)&&(offset<1.19e-07)){
    taps = std::vector<float>(proto.begin(), proto.end());
  }
  else{
    taps = std::vector<float>(proto.size(),0.);

    for(int idx = 0; idx < proto.size(); idx++){
      for(int ind = 0; ind < proto.size(); ind++){
        taps[idx] += proto[ind]*
            boost::math::sinc_pi(M_PI*(float(idx)-float(ind) - offset));
      }
    }
  }

}

#endif /* INCLUDED_SIGNAL_BASE_HPP */
