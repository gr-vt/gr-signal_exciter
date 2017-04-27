#ifndef INCLUDED_SIGNAL_CPM_HPP
#define INCLUDED_SIGNAL_CPM_HPP

#include <signal_exciter/signal_base.hpp>
#include <gnuradio/analog/cpm.h>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/sincos.h>
#include <volk/volk.h>

class Signal_CPM : public Signal_Base
{
  private:
    int                         d_seed;
    int                         d_order;
    int                         d_sps;
    int                         d_L;
    float                       d_h;
    double                      d_beta;
    bool                        d_first_pass;

    double                      d_phase_acm;

    gr::analog::cpm::cpm_type   d_cpm_type;

    std::vector<float>          d_phase_shape;

    std::vector<float>          d_symbol_amps;
    std::vector<complexf>       d_symbol_list;

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
    size_t                            d_hist;
    size_t                            d_samp_offset;
    float*                            d_filt_in;
    float*                            d_filt_out;
    std::vector<float>                d_past;
    std::vector< std::vector<float> > d_taps;
    std::vector<gr::filter::kernel::fir_filter_fff *> d_firs;

    float d_fso;
    std::vector<float> d_proto_taps;

    void load_firs();
    void filter( size_t nout, complexf* out );

  public:
    Signal_CPM(int order, gr::analog::cpm::cpm_type phase_type, int sps,
              int overlap, float mod_idx, int seed, double beta=0.3,
              float* phase_shape=NULL, size_t phase_shape_length=0,
              float fso=0., bool enable=true, size_t buff_size=8192,
              size_t min_notify=512);
    ~Signal_CPM();

    void generate_signal(complexf* output, size_t sample_count);
    void generate_symbols(complexf* output, size_t symbol_count);

    void reset(){ d_first_pass = true; }

    void set_seed(int seed)
    {
      d_seed = seed;
      if(d_seed < 0) d_seed = d_rd();
    }
};
#endif /* INCLUDED_SIGNAL_CPM_HPP */
