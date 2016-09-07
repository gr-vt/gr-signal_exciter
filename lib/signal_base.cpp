
#include <signal_exciter/signal_base.hpp>


boost::mutex Signal_Base::s_mutex_fftw;
/*bool Signal_Base::d_seeded = false;
int Signal_Base::d_seed = -1;


void Signal_Base::set_seed(int seed){
  if(!d_seeded){
    d_seed = seed;
    if(d_seed<0){
      d_seed = time(NULL);
    }
    srand(d_seed);
    d_seeded = true;
  }
}

int Signal_Base::get_seed(){
  return d_seed;
}*/


