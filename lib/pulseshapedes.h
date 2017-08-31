#ifndef INCLUDED_SIGNAL_EXCITER_PULSESHAPEDES_H
#define INCLUDED_SIGNAL_EXCITER_PULSESHAPEDES_H

#include <signal_exciter/api.h>
#include <cmath>
#include <boost/math/special_functions/sinc.hpp>

#ifndef M_PIl
#define M_PIl          3.141592653589793238462643383279502884L
#endif

namespace gr {
namespace signal_exciter{
 /* Taken from https://stackoverflow.com/questions/24518989/how-to-perform-1-dimensional-valid-convolution
  * User: 101010
  * Modified
  */
  template<typename T>
  inline std::vector<T>
  conv(std::vector<T> &out, std::vector<T> const &f, std::vector<T> const &g) {
    int const nf = f.size();
    int const ng = g.size();
    int const n  = nf + ng - 1;
    out = std::vector<T>(n, T());
    for(size_t i(0); i < n; ++i) {
      int const jmn = (i >= ng - 1)? i - (ng - 1) : 0;
      int const jmx = (i <  nf - 1)? i            : nf - 1;
      for(size_t j(jmn); j <= jmx; ++j) {
        out[i] += (f[j] * g[i - j]);
      }
    }
    return out;
  }

  class pulseshapedes {
   public:
    enum pulse_type {
      PS_NONE = -1,
      PS_RECT = 0,
      PS_RCOS = 1,
    };

    template <class FPprec>
    static void
    rect_pulse( std::vector<FPprec> &pulse_shape, FPprec span, int sps )
    {
      if(span < 1.) span = 1.;
      size_t psl = size_t(std::ceil(span*FPprec(sps)));
      psl = psl + (1-psl%2);//want odd
      pulse_shape = std::vector<FPprec>(psl,FPprec(0.));
      size_t remainder = psl-sps;
      size_t head = size_t(std::ceil(FPprec(remainder)/FPprec(2.)));
      FPprec val = FPprec(1.)/FPprec(sps);
      for(size_t idx = head; idx < head+sps; idx++){
        pulse_shape[idx] = val;
      }
    }

    template <class FPprec>
    static void
    rcos_pulse( std::vector<FPprec> &pulse_shape, FPprec span, int sps, FPprec beta )
    {
      if(span < 1.) span = 1.;

      FPprec lpi = FPprec(M_PIl);
      FPprec lsps = FPprec(sps);

      size_t psl = size_t(std::ceil(span*lsps));
      psl = psl + (1-psl%2);//want odd
      pulse_shape = std::vector<FPprec>(psl,FPprec(0.));

      FPprec t = -(FPprec(psl-1)/FPprec(2.));
      for(size_t idx = 0; idx < pulse_shape.size(); idx++, t+=FPprec(1.)){
        if(std::abs(t) == (lsps/(FPprec(2.)*beta))){
          pulse_shape[idx] = (FPprec(2.)*beta)*std::sin(lpi/(FPprec(2.)*beta))/(FPprec(4.)*sps);
        }
        else{
          pulse_shape[idx] = ( boost::math::sinc_pi(lpi*t/lsps) )
                * std::cos(lpi * beta * t / lsps)
                / (1 - ((FPprec(2.) * beta * t) / lsps) 
                  * ((FPprec(2.) * beta * t) / lsps)) / lsps;
        }
      }
    }


  };

} /* signal_exciter */
} /* gr */

#endif /* INCLUDED_SIGNAL_EXCITER_PULSESHAPEDES_H */
