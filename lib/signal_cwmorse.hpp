#ifndef INCLUDED_SIGNAL_CWMORSE_HPP
#define INCLUDED_SIGNAL_CWMORSE_HPP


#include <signal_exciter/signal_base.hpp>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/filter/firdes.h>
#include <volk/volk.h>
#include <algorithm>


class Signal_CWMORSE : public Signal_Base
{
  private:
    int d_seed;
    int d_cpw;
    float d_wpm;
    bool d_word;
    bool d_first_pass;
    size_t d_letter_count;
    int d_state;
    std::vector< std::vector<size_t> > d_alphabet;
    std::vector< char > d_chars;
    std::vector< size_t > d_syllable;

    std::vector< std::vector<complexf> > d_symbol_list;
    std::vector<complexf> d_char_end;
    std::vector<complexf> d_word_end;

    std::vector<complexf> d_symbol_cache;

    int d_samp_offset;
    //std::vector< int > d_data_cache;
    std::vector<complexf> d_output_cache;
    size_t d_protected;

    size_t d_char_in_word;
    std::vector<complexf> d_leftover_symbol;
    size_t d_leftover_count;

    size_t d_burn;


    void create_symbol_list();
    void print_symbol_list();

    size_t d_enable;
    size_t d_buffer_size;
    size_t d_notify_size;
    bool d_running;
    void auto_fill_symbols();
    void auto_fill_signal();
    signal_threaded_buffer<complexf>*   d_Sy;

    boost::thread_group d_TGroup;

    void auto_gen_SYMS();

    float d_fso;
    std::vector<float> d_proto_taps;

    //volk things
    int d_align;

    // adding in interp option
    size_t d_interp;
    size_t d_branch_offset;
    std::vector<float> d_interp_taps;
    std::vector< gr::filter::kernel::fir_filter_ccf* > d_firs;
    std::vector< std::vector<float> > d_taps;
    std::vector<complexf> d_past;
    size_t d_hist;
    complexf* d_filt_in;

    void load_firs();
    void filter( size_t nout, complexf* out );

  public:
    Signal_CWMORSE(int d_char_per_word, float words_per_minute, bool base_word,
                    int seed, float* interp_taps=NULL, size_t tap_len=0,
                    int interp=1, float fso=0., bool enable=true,
                    size_t buff_size=8192, size_t min_notify=512);
    ~Signal_CWMORSE();

    void generate_signal(complexf* output, size_t sample_count);
    void generate_symbols(complexf* output, size_t symbol_count);

    void reset(){ d_first_pass = true; }

    void set_seed(int seed)
    {
      d_seed = seed;
      if(d_seed < 0) d_seed = d_rd();
    }

};

#endif /* INCLUDED_SIGNAL_CWMORSE_HPP */
