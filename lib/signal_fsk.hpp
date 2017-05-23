#ifndef INCLUDED_Signal_FSK_HPP
#define INCLUDED_Signal_FSK_HPP

#include <signal_exciter/signal_base.hpp>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/filter/firdes.h>
#include <gnuradio/sincos.h>
#include <volk/volk.h>

class Signal_FSK : public Signal_Base
{
  private:
    int                         d_seed;
    int                         d_order;
    int                         d_sps;
    float                       d_h;
    float                       d_f;
    bool                        d_first_pass;

    std::vector<float>          d_symbol_amps;
    std::vector<complexf>       d_symbol_cache;
    long                        d_n;
    size_t                      d_samp_offset;

    void throw_runtime(std::string err);

    void create_symbol_list();

    // Auto gen?
    size_t d_enable;
    size_t d_buffer_size;
    size_t d_notify_size;
    bool d_running;
    void auto_fill_symbols();
    void auto_fill_signal();
    signal_threaded_buffer<complexf>*   d_Sy;

    boost::thread_group d_TGroup;

    void auto_gen_SYMS();

    // Phase Filtering
    //volk things
    int                               d_align;

    std::vector<float> d_proto_taps;

    //volk things
    complexf* d_time_shift_in;
    gr::filter::kernel::fir_filter_ccf* d_frac_filt;
    std::vector<complexf> d_frac_cache;

    void load_firs();
    void filter( size_t nout, complexf* out );

  public:
    Signal_FSK(int order, int sps, float mod_idx, int seed,
              bool enable_fso = false, float fso=0., bool enable=true, size_t buff_size=8192,
              size_t min_notify=512);
    ~Signal_FSK();

    void generate_signal(complexf* output, size_t sample_count);
    void generate_symbols(complexf* output, size_t symbol_count);

    void reset(){ d_first_pass = true; }

    void set_seed(int seed)
    {
      d_seed = seed;
      if(d_seed < 0) d_seed = d_rd();
    }
};
#endif /* INCLUDED_Signal_FSK_HPP */
