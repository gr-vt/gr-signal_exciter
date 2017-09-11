
#include <signal_exciter/signal_base.hpp>
#include <iostream>

boost::mutex Signal_Base::s_mutex_fftw;
//boost::mutex Signal_Base::s_mutex_fso;
boost::mutex* Signal_Base::d_mutex_ptr = &s_mutex_fftw;
size_t Signal_Base::s_indicator = 0;

void
Signal_Base::norm_f(std::vector<float> &to_norm, float scale)
{
  double wt = 0.;
  for(size_t idx = 0; idx < to_norm.size(); idx++){
    wt += to_norm[idx]*to_norm[idx];
  }
  if(wt) wt = std::sqrt(scale/wt);
  else wt = 1.;
  for(size_t idx = 0; idx < to_norm.size(); idx++){
    to_norm[idx] *= wt;
  }
}

void
Signal_Base::norm_d(std::vector<double> &to_norm, double scale)
{
  double wt = 0.;
  for(size_t idx = 0; idx < to_norm.size(); idx++){
    wt += to_norm[idx]*to_norm[idx];
  }
  if(wt) wt = std::sqrt(scale/wt);
  else wt = 1.;
  for(size_t idx = 0; idx < to_norm.size(); idx++){
    to_norm[idx] *= wt;
  }
}

void
Signal_Base::print_f(std::string Name, float* vec, size_t len)
{
  std::cout << d_indicator << ": ";
  if(len){
    std::cout << Name << " = [" << vec[0];
    for(size_t idx = 1; idx < len; idx++){
      std::cout << ", " << vec[idx];
    }
    std::cout << "];\n";
  }
  else{
    std::cout << Name << " = [];\n";
  }
}

void
Signal_Base::print_fcr(std::string Name, complexf* vec, size_t len)
{
  std::cout << d_indicator << ": ";
  if(len){
    std::cout << Name << " = [" << vec[0].real();
    for(size_t idx = 1; idx < len; idx++){
      std::cout << ", " << vec[idx].real();
    }
    std::cout << "];\n";
  }
  else{
    std::cout << Name << " = [];\n";
  }
}

void
Signal_Base::print_d(std::string Name, double* vec, size_t len)
{
  std::cout << d_indicator << ": ";
  if(len){
    std::cout << Name << " = [" << vec[0];
    for(size_t idx = 1; idx < len; idx++){
      std::cout << ", " << vec[idx];
    }
    std::cout << "];\n";
  }
  else{
    std::cout << Name << " = [];\n";
  }
}
