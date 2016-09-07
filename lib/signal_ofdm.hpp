#ifndef INCLUDED_SIGNAL_OFDM_HPP
#define INCLUDED_SIGNAL_OFDM_HPP

#include <signal_exciter/signal_base.hpp>
#include <signal_exciter/random_signal_config.h>
#include "signal_psk.hpp"
#include "signal_qam.hpp"
#include "signal_pam.hpp"
#include <fftw3.h>
#include <gnuradio/filter/fir_filter.h>
#include <volk/volk.h>


class Signal_OFDM : public Signal_Base
{
  private:
    std::vector<complexf> d_output_cache;
    std::vector<complexf> d_prefft_cache;
    std::vector<complexf> d_postfft_cache;

    size_t d_fftsize;
    size_t d_cp_len;
    size_t d_active;
    size_t d_spf;
    bool d_ppf;
    size_t d_npilots;
    std::vector<size_t> d_active_list;
    std::vector<size_t> d_pilot_list;
    std::vector<complexf> d_pilots;
    float d_backoff;//this is in dB
    bool d_sync;

    size_t d_protected;
    size_t d_symbol_idx;
    bool d_first_pass;

    std::vector<size_t> d_pn1;
    std::vector<size_t> d_pn2;

    Signal_Base* d_mod;

    fftwf_plan d_fft;
    fftwf_complex* d_fft_in;
    fftwf_complex* d_fft_out;
    void do_fft();

    void sort_subcarriers();
    void generate_sync();
    void generate_pilots();

    void print_buffer_complex(std::string str, fftwf_complex* buff, size_t length);
    void print_buffer_complexf(std::string str, complexf* buff, size_t length);

    size_t d_sym_counter;
    size_t items_written;

    size_t d_enable;
    size_t d_buffer_size;
    size_t d_notify_size;
    bool d_running;
    void auto_fill_symbols();
    void auto_fill_signal();
    signal_threaded_buffer<float>*      d_G1;
    signal_threaded_buffer<float>*      d_G2;
    signal_threaded_buffer<float>*      d_Mx;
    signal_threaded_buffer<float>*      d_Fn;
    signal_threaded_buffer<complexf>*   d_Sy;

    boost::thread_group d_TGroup;


    size_t d_samp_gen_count;
    size_t d_symb_gen_count;

    int d_interp;
    size_t d_branch_offset;//starting branch in interp filter
    size_t d_samp_offset;//starting sample in frame
    size_t d_frame_offset;//starting frame in cache
    std::vector<float> d_window;
    std::vector< gr::filter::kernel::fir_filter_ccf* > d_firs;
    std::vector< std::vector<float> > d_taps;
    size_t d_hist;
    std::vector< complexf > d_past;

    //volk things
    int d_align;
    complexf* d_filt_in;

    void load_firs();
    void filter( size_t nout, complexf* out);

    void generate_frame();
    std::vector<complexf> d_frame;
    //std::vector< std::vector<complexf> > d_frame_cache;
    std::vector<complexf> d_frame_cache;
    bool d_frame_cache_enlarged;
    size_t d_symbol_length;
    size_t d_frame_length;
    size_t d_frames_needed;
    size_t d_cached;

    void throw_runtime(std::string err);

  public:
    Signal_OFDM(size_t fftsize, size_t cp_len, size_t active_carriers, size_t syms_per_frame, 
                bool pilot_per_frame, size_t pilot_count, size_t* pilot_locations, float backoff,
                int mod_type, int mod_order, float mod_offset, int seed, bool add_sync=false,
                float* window=NULL, size_t length=0, int interp=1,
                bool enable=true, size_t buff_size=8192, size_t min_notify=512);
    ~Signal_OFDM();

    void generate_signal(complexf* output, size_t sample_count);
    void generate_symbols(complexf* output, size_t symbol_count);

};

#endif /* INCLUDED_SIGNAL_OFDM_HPP */
