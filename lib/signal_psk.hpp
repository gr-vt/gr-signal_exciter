#ifndef INCLUDED_SIGNAL_PSK_HPP
#define INCLUDED_SIGNAL_PSK_HPP


#include <signal_exciter/signal_base.hpp>
#include <gnuradio/filter/fir_filter.h>
#include <volk/volk.h>


class Signal_PSK : public Signal_Base
{
  private:
    int d_seed;
    int d_order;
    float d_offset;
    size_t d_overlap;
    std::vector<float> d_pulse_shape;
    int d_sps;
    bool d_first_pass;

    std::vector<float> d_proto_taps;

    std::vector<complexf> d_symbol_list;

    std::vector<complexf> d_symbol_cache;

    size_t d_samp_offset;
    //std::vector< int > d_data_cache;
    std::vector<complexf> d_output_cache;
    size_t d_protected;

    void create_symbol_list();

    size_t d_enable;
    size_t d_buffer_size;
    size_t d_notify_size;
    bool d_running;
    void auto_fill_symbols();
    void auto_fill_signal();
    signal_threaded_buffer<complexf>*   d_Sy;

    boost::thread_group d_TGroup;

    void auto_gen_SYMS();

    std::vector<gr::filter::kernel::fir_filter_ccf *> d_firs;
    void load_firs();
    std::vector< std::vector<float> > d_taps;
    size_t d_hist;

    //volk things
    int d_align;
    complexf* d_filt_in;

    void throw_runtime(std::string err, size_t sc=0, size_t os=0, size_t si=0);

    size_t d_sample_count;
    size_t d_symbol_count;

    size_t d_samp_gen_count;
    size_t d_symb_gen_count;

    std::vector<complexf> d_past;
    void filter( size_t nout, complexf* out);

  public:
    Signal_PSK(int order, float offset, int sps, float* pusle_shape,
                size_t length, int seed, bool enable=true,
                size_t buff_size=8192, size_t min_notify=512);
    ~Signal_PSK();

    void generate_signal(complexf* output, size_t sample_count);
    void generate_symbols(complexf* output, size_t symbol_count);

    void reset(){ d_first_pass = true; }

    void set_seed(int seed)
    {
      d_seed = seed;
      if(d_seed < 0) d_seed = (*d_rd)();
    }

};

#endif /* INCLUDED_SIGNAL_PSK_HPP */
