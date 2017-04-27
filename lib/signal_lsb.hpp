#ifndef INCLUDED_SIGNAL_LSB_HPP
#define INCLUDED_SIGNAL_LSB_HPP


#include <signal_exciter/signal_base.hpp>
#include <signal_exciter/gmm_spectral_taps.h>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/filter/firdes.h>
#include <gnuradio/analog/agc.h>
#include <volk/volk.h>
#include <fftw3.h>


class Signal_LSB : public Signal_Base
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

    std::vector<float> d_fm;
    fftwf_plan d_fft;
    fftwf_plan d_ifft;
    fftwf_complex* d_fft_in;
    fftwf_complex* d_fft_out;
    void do_fft(size_t samp_count);
    void do_ifft(size_t samp_count);

    size_t d_enable;
    size_t d_buffer_size;
    size_t d_notify_size;
    bool d_running;
    void auto_fill_symbols();
    void auto_fill_signal();
    signal_threaded_buffer<complexf>*   d_Sy;

    boost::thread_group d_TGroup;

    float d_fso;
    std::vector<float> d_proto_taps;

    void auto_gen_GM();

    //volk things
    int d_align;
    float* d_filt_in;
    complexf* d_time_shift_in;
    gr::filter::kernel::fir_filter_ccf* d_frac_filt;
    std::vector<complexf> d_frac_cache;

  public:
    Signal_LSB(float mod_idx, size_t components, float* mu, float* sigma,
                float* weight, float samp_rate, size_t tap_count, int seed,
                bool norm=false, float fso=0., bool enable=true,
                size_t buff_size=8192, size_t min_notify=512);
    ~Signal_LSB();

    void generate_signal(complexf* output, size_t sample_count);
    void generate_symbols(complexf* output, size_t symbol_count);

    void reset(){ d_first_pass = true; }

    void set_seed(int seed)
    {
      d_seed = seed;
      if(d_seed < 0) d_seed = d_rd();
    }

};

#endif /* INCLUDED_SIGNAL_LSB_HPP */
