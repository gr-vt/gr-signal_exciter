#ifndef INCLUDED_SIGNAL_DSB_HPP
#define INCLUDED_SIGNAL_DSB_HPP


#include <signal_exciter/signal_base.hpp>
#include <signal_exciter/gmm_spectral_taps.h>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/analog/agc.h>
#include <volk/volk.h>


class Signal_DSB : public Signal_Base
{
  private:
    int d_seed;
    float d_mod_idx;
    bool d_first_pass;

    size_t d_tap_count;
    std::vector<float> d_taps;
    std::vector<float> d_window;
    GMM_Spectral_Taps d_gmm_tap_gen;
    void generate_taps();

    gr::filter::kernel::fir_filter_fff* d_fir;
    void load_fir();
    size_t d_hist;
    std::vector<float> d_past;
    void filter( size_t nout, complexf* output );

    std::vector<complexf> d_symbol_cache;
    std::vector<complexf> d_output_cache;

    bool d_norm;
    gr::analog::kernel::agc_cc d_agc;
    size_t d_burn;

    //void filter(float* n, float* m, size_t length);

    size_t d_enable;
    size_t d_buffer_size;
    size_t d_notify_size;
    bool d_running;
    void auto_fill_symbols();
    void auto_fill_signal();
/*    signal_threaded_buffer<float>*      d_G1;
    signal_threaded_buffer<float>*      d_G2;
    signal_threaded_buffer<float>*      d_Mx;
    signal_threaded_buffer<float>*      d_Fn;*/
    signal_threaded_buffer<complexf>*   d_Sy;

    boost::thread_group d_TGroup;

/*    void gen_gaussian(signal_threaded_buffer<float>* buffer, float var);
    void gaussian_mix(signal_threaded_buffer<float>* obuff, signal_threaded_buffer<float>* ibuff1, signal_threaded_buffer<float>* ibuff2, float thresh);
    void auto_filter(signal_threaded_buffer<float>* obuff, signal_threaded_buffer<float>* ibuff);
    void auto_load_symbols(signal_threaded_buffer<complexf>* obuff, signal_threaded_buffer<float>* ibuff);*/

    void auto_gen_GM();

    //volk things
    int d_align;
    float* d_filt_in;

  public:
//    Signal_DSB(float mod_idx, float f_max, float var1, float var2, float thresh, int seed,
//                bool norm=false, bool enable=true, size_t buff_size=8192, size_t min_notify=512);
    Signal_DSB(float mod_idx, size_t components, float* mu, float* sigma, float* weight, float samp_rate, size_t tap_count, int seed,
                bool norm = false, bool enable=true, size_t buff_size=8192, size_t min_notify=512);
    ~Signal_DSB();

    void generate_signal(complexf* output, size_t sample_count);
    void generate_symbols(complexf* output, size_t symbol_count);

    void reset(){ d_first_pass = true; }

    void set_seed(int seed)
    {
      d_seed = seed;
      if(d_seed < 0) d_seed = d_rd();
    }

};

#endif /* INCLUDED_SIGNAL_DSB_HPP */
