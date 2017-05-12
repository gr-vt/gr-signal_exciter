

#include "signal_cwmorse.hpp"
#include <stdio.h>//////////////////////////////////

Signal_CWMORSE::Signal_CWMORSE(int char_per_word, float words_per_minute,
                            bool base_word, int seed, float* interp_taps,
                            size_t tap_len, int interp, float fso, bool enable,
                            size_t buff_size, size_t min_notify)
  : d_cpw(char_per_word),
    d_wpm(words_per_minute),
    d_word(base_word),
    d_leftover_count(0),
    d_char_in_word(0),
    d_samp_offset(0),
    d_state(0),
    d_protected(0),
    d_interp(interp),
    d_branch_offset(0),
    d_enable(enable),
    d_buffer_size(buff_size)
{
  //printf("Init.\n");
  set_seed(seed);
  //printf("Seeded.\n");
  d_burn = buff_size;

  d_first_pass = true;

  if(tap_len){
    double power_check = 0.;
    d_interp_taps = std::vector<float>(tap_len);
    for(size_t idx = 0; idx < tap_len; idx++){
      d_interp_taps[idx] = interp_taps[idx];
      power_check += interp_taps[idx]*interp_taps[idx];
    }
    double normalizer = sqrt(double(interp)/power_check);
    for(size_t idx = 0; idx < tap_len; idx++){
      d_interp_taps[idx] *= normalizer;
    }
  }
  else{
    d_interp = 1;
    d_interp_taps = gr::filter::firdes::low_pass_2(1,1,0.5,0.05,61,
                          gr::filter::firdes::WIN_BLACKMAN_hARRIS);
  }

  d_fso = fso;

  d_align = volk_get_alignment();

  load_firs();
  create_symbol_list();

  d_rng = new gr::random(d_seed, 0, d_letter_count);

  if(d_enable){
    d_running = true;
    auto_fill_symbols();
    auto_fill_signal();
  }
}

Signal_CWMORSE::~Signal_CWMORSE()
{
  if(d_enable){
    d_running = false;
    d_TGroup.join_all();
    delete d_Sy;
  }
  delete d_rng;
  for(size_t idx = 0; idx < d_interp; idx++){
    delete d_firs[idx];
  }
}

void
Signal_CWMORSE::generate_symbols(complexf* output, size_t symbol_count)
{
  if(d_enable){
    size_t filled(0);
    while((filled < symbol_count) && d_running){
      filled += d_Sy->bmemcpy( &output[filled], symbol_count-filled, false );
    }
  }
  else{
    size_t oo(0), to_fill(0);
    if(d_leftover_symbol.size()){
      to_fill = std::min(symbol_count-oo,d_leftover_symbol.size());
      memcpy( &output[oo], &d_leftover_symbol[0], sizeof(complexf)*to_fill );
      oo += to_fill;
      if(to_fill < d_leftover_symbol.size()){
        d_leftover_symbol = std::vector<complexf>(d_leftover_symbol.begin()+to_fill, d_leftover_symbol.end());
      }
    }
    while(oo < symbol_count){
      //int data = rand()%d_letter_count;
      if(d_state == 0){
        int data = d_rng->ran_int();
        to_fill = std::min(symbol_count-oo, d_symbol_list[data].size());
        memcpy( &output[oo], &d_symbol_list[data][0], sizeof(complexf)*to_fill );
        oo += to_fill;
        if(to_fill < d_symbol_list[data].size()){
          d_leftover_symbol = std::vector<complexf>(d_symbol_list[data].begin()+to_fill, d_symbol_list[data].end());
        }
        d_char_in_word++;
        if(d_char_in_word < d_cpw){
          d_state = 1;
        }
        else{
          d_state = 2;
        }
      }
      else if(d_state = 1){
        to_fill = std::min(symbol_count-oo, d_char_end.size());
        memcpy( &output[oo], &d_char_end[0], sizeof(complexf)*to_fill );
        oo += to_fill;
        if(to_fill < d_char_end.size()){
          d_leftover_symbol = std::vector<complexf>(d_char_end.begin()+to_fill, d_char_end.end());
        }
        d_state = 0;
      }
      else{
        to_fill = std::min(symbol_count-oo, d_word_end.size());
        memcpy( &output[oo], &d_word_end[0], sizeof(complexf)*to_fill );
        oo += to_fill;
        if(to_fill < d_word_end.size()){
          d_leftover_symbol = std::vector<complexf>(d_word_end.begin()+to_fill, d_word_end.end());
        }
        d_state = 0;
        d_char_in_word = 0;
      }
    }
  }
}

void
Signal_CWMORSE::generate_signal(complexf* output, size_t sample_count)
{
  if(d_first_pass){
    d_past = std::vector<complexf>(d_hist);
    generate_symbols( &d_past[0], d_past.size() );
    d_first_pass = false;
  }
  filter( sample_count, output );

}

void
Signal_CWMORSE::filter( size_t nout, complexf* out )
{
  size_t total_samps = d_past.size()*d_interp;
  size_t used_samps = d_branch_offset;
  size_t part_samps = (d_interp-(d_branch_offset%d_interp))%d_interp;

  total_samps -= used_samps;
  total_samps -= d_hist*d_interp;

  size_t samps_needed;
  if(nout < total_samps){
    samps_needed = 0;
  }
  else{
    samps_needed = nout - total_samps;
  }
  float fractional_N = float(samps_needed);
  float fractional_D = float(d_interp);
  float fractional = fractional_N/fractional_D;
  size_t symbs_needed = ceil(fractional);

  size_t total_input_len = symbs_needed + d_past.size();

  d_filt_in = (complexf*)volk_malloc( total_input_len*sizeof(complexf), d_align);
  memcpy( &d_filt_in[0], &d_past[0], d_past.size()*sizeof(complexf) );
  generate_symbols( &d_filt_in[d_past.size()], symbs_needed );

  size_t oo(0),ii(0);
  while( oo < nout ){
    out[oo] = d_firs[d_branch_offset]->filter( &d_filt_in[ii] );
    d_branch_offset = (d_branch_offset+1)%d_interp;
    if(!d_branch_offset){
      ii++;
    }
    oo++;
  }
  size_t remaining = total_input_len - ii;
  d_past = std::vector<complexf>( &d_filt_in[ii], &d_filt_in[total_input_len] );
  volk_free(d_filt_in);
}



void
Signal_CWMORSE::create_symbol_list()
{
/*
A: 10111
B: 111010101
C: 11101011101
D: 1110101
E: 1
F: 101011101
G: 111011101
H: 1010101
I: 101
J: 1011101110111
K: 111010111
L: 101110101
M: 1110111
N: 11101
O: 11101110111
P: 10111011101
Q: 1110111010111
R: 1011101
S: 10101
T: 111
U: 1010111
V: 101010111
W: 101110111
X: 11101010111
Y: 1110101110111
Z: 11101110101
0: 11101110111011101111
1: 10111011101110111
2: 101011101110111
3: 1010101110111
4: 10101010111
5: 101010101
6: 11101010101
7: 1110111010101
8: 111011101110101
9: 11101110111011101
-: 111010101010111
(: 1110101110111010111
): 111010111011101
.: 10111010111010111
/: 1110101011101
': 1011101110111011101
:: 11101110111010101
,: 1110111010101110111
?: 101011101110101
!: 1110101110101110111
;: 11101011101011101
=: 1110101010111
+: 1011101011101
": 101110101011101
*/
  size_t syllable_count[] = {2,4,4,3,1,4,3,4,2,4,3,4,2,2,3,4,4,3,3,1,3,4,3,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,5,6,5,6,6,6,6,6,6,5,5,6};
  d_syllable = std::vector<size_t>(syllable_count, syllable_count+sizeof(syllable_count)/sizeof(size_t));
  size_t syllab[] = {1,3,3,1,1,1,3,1,3,1,3,1,1,1,1,1,3,1,3,3,1,1,1,1,1,1,1,1,3,3,3,3,1,3,1,3,1,1,3,3,3,1,3,3,3,1,3,3,1,3,3,1,3,
                     1,3,1,1,1,1,3,1,1,3,1,1,1,3,1,3,3,3,1,1,3,3,1,3,3,3,3,1,1,3,3,3,3,3,1,3,3,3,3,1,1,3,3,3,1,1,1,3,3,1,1,1,1,
                     3,1,1,1,1,1,3,1,1,1,1,3,3,1,1,1,3,3,3,1,1,3,3,3,3,1,3,1,1,1,1,3,3,1,3,3,1,3,3,1,3,3,1,1,3,1,3,1,3,3,1,1,3,
                     1,1,3,3,3,3,1,3,3,3,1,1,1,3,3,1,1,3,3,1,1,3,3,1,1,3,1,3,1,3,3,3,1,3,1,3,1,3,1,1,1,3,1,3,1,3,1,1,3,1,1,3,1};
  std::vector<size_t> syllables(syllab, syllab + sizeof(syllab)/sizeof(size_t));
  char chars[] = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','0',
'1','2','3','4','5','6','7','8','9','-','(',')','.','/','\'',':',',','?','!',';','=','+','"'};
  d_chars = std::vector<char>(chars, chars+sizeof(chars)/sizeof(char));
  d_letter_count = d_chars.size();
  size_t syb_pointer = 0;
  size_t one_count = 0;
  size_t zero_count = 0;
  for(size_t chr = 0; chr < d_syllable.size(); chr++){
    //size_t syb_start = syb_pointer;
    size_t chr_space = 0;
    for(size_t syb = 0; syb < d_syllable[chr]; syb++){
      chr_space += syllables[syb_pointer+syb];
      one_count += syllables[syb_pointer+syb];
    }
    chr_space += d_syllable[chr]-1;
    zero_count += chr_space;
    d_alphabet.push_back(std::vector<size_t>(chr_space,0));
    d_symbol_list.push_back(std::vector<complexf>(chr_space,complexf(0.,0.)));
    size_t chr_pointer = 0;
    for(size_t syb = 0; syb < d_syllable[chr]; syb++){
      for(size_t cnt = 0; cnt < syllables[syb_pointer+syb]; cnt++){
        d_symbol_list[chr][chr_pointer] = complexf(1.,0.);
        d_alphabet[chr][chr_pointer++] = 1;
      }
      if(syb < d_syllable[chr]-1){
        chr_pointer++;
      }
    }
    syb_pointer += d_syllable[chr];
  }
  zero_count = zero_count - one_count + 3+7;

  //printf("1: %lu, 0: %lu\n",one_count,zero_count);

  d_char_end = std::vector<complexf>(3, complexf(0.,0.));
  d_word_end = std::vector<complexf>(7, complexf(0.,0.));


  //print_symbol_list();
}

void
Signal_CWMORSE::print_symbol_list(void)
{
  for(size_t idx1 = 0; idx1 < d_alphabet.size(); idx1++){
    printf("%c:\t%lu",d_chars[idx1],d_alphabet[idx1][0]);
    for(size_t idx2 = 1; idx2 < d_alphabet[idx1].size(); idx2++){
      printf(", %lu",d_alphabet[idx1][idx2]);
    }
    printf("\n");
  }
  for(size_t idx1 = 0; idx1 < d_symbol_list.size(); idx1++){
    printf("%c:\t%0.0f",d_chars[idx1],d_symbol_list[idx1][0].real());
    for(size_t idx2 = 1; idx2 < d_symbol_list[idx1].size(); idx2++){
      printf(", %0.0f",d_symbol_list[idx1][idx2].real());
    }
    printf("\n");
  }
  printf("char end:\t0, 0, 0\n");
  printf("word end:\t0, 0, 0, 0, 0, 0, 0\n");
}

void
Signal_CWMORSE::load_firs()
{
  size_t intp = d_interp;

  d_firs = std::vector< gr::filter::kernel::fir_filter_ccf *>(intp);
  std::vector<float> dummy_taps;
  for(size_t idx = 0; idx < intp; idx++){
    d_firs[idx] = new gr::filter::kernel::fir_filter_ccf(1,dummy_taps);
  }

  size_t leftover = (intp - (d_interp_taps.size() % intp))%intp;
  d_proto_taps = std::vector<float>(d_interp_taps.size() + leftover, 0.);

  memcpy( &d_proto_taps[0], &d_interp_taps[0],
          d_interp_taps.size()*sizeof(float) );


  //std::vector<float> shifted_taps;
  //time_offset(shifted_taps, d_proto_taps, d_interp*d_fso);
  std::vector<float> shifted_taps = d_proto_taps;

  //std::vector< std::vector<float> > xtaps(intp);
  d_taps = std::vector< std::vector<float> >(intp);

  size_t ts = shifted_taps.size() / intp;
  for(size_t idx = 0; idx < intp; idx++){
    d_taps[idx].resize(ts);
  }
  //printf("OFDM:: taps init 0.\n");

  for(size_t idx = 0; idx < d_interp_taps.size(); idx++){
    d_taps[idx % intp][idx / intp] = shifted_taps[idx];
  }
  //printf("OFDM:: taps filled.\n");

  //printf("OFDM:: filters made.\n");
  for(size_t idx = 0; idx < intp; idx++){
    d_firs[idx]->set_taps(d_taps[idx]);
  }
  //printf("OFDM:: taps loaded.\n");
  d_hist = ts-1;
}


void
Signal_CWMORSE::auto_fill_symbols()
{
  d_Sy = new signal_threaded_buffer<complexf>(d_buffer_size, d_notify_size);
  // Start the symbol generation thread
  d_TGroup.create_thread( boost::bind(&Signal_CWMORSE::auto_gen_SYMS, this) );
}

void
Signal_CWMORSE::auto_fill_signal()
{}


void
Signal_CWMORSE::auto_gen_SYMS()
{
  size_t buff_size(d_buffer_size), buff_pnt(0);
  std::vector<complexf> buffer(d_buffer_size,complexf(0.,0.));
  std::vector<complexf> leftover(0);
  size_t to_fill(0);
  while(d_running){
    if(leftover.size()){
      to_fill = leftover.size();
      memcpy( &buffer[buff_pnt], &leftover[0], sizeof(complexf)*to_fill );
      buff_pnt += to_fill;
      leftover = std::vector<complexf>(0);
    }
    while(buff_pnt < buff_size){
      //int data = rand()%d_letter_count;
      if(d_state == 0){
        int data = d_rng->ran_int();
        to_fill = std::min(buff_size-buff_pnt, d_symbol_list[data].size());
        memcpy( &buffer[buff_pnt], &d_symbol_list[data][0], sizeof(complexf)*to_fill );
        buff_pnt += to_fill;
        if(to_fill < d_symbol_list[data].size()){
          leftover = std::vector<complexf>(d_symbol_list[data].begin()+to_fill, d_symbol_list[data].end());
        }
        d_char_in_word++;
        if(d_char_in_word < d_cpw){
          d_state = 1;
        }
        else{
          d_state = 2;
        }
      }
      else if(d_state == 1){
        to_fill = std::min(buff_size-buff_pnt, d_char_end.size());
        memcpy( &buffer[buff_pnt], &d_char_end[0], sizeof(complexf)*to_fill );
        buff_pnt += to_fill;
        if(to_fill < d_char_end.size()){
          leftover = std::vector<complexf>(d_char_end.begin()+to_fill, d_char_end.end());
        }
        d_state = 0;
      }
      else{
        to_fill = std::min(buff_size-buff_pnt, d_word_end.size());
        memcpy( &buffer[buff_pnt], &d_word_end[0], sizeof(complexf)*to_fill );
        buff_pnt += to_fill;
        if(to_fill < d_word_end.size()){
          leftover = std::vector<complexf>(d_word_end.begin()+to_fill, d_word_end.end());
        }
        d_state = 0;
      }
    }
    buff_pnt = 0;
    while((buff_pnt < buff_size) && d_running){
      buff_pnt += d_Sy->bmemcpy( &buffer[buff_pnt], buff_size-buff_pnt, true );
    }
    buff_pnt = 0;
  }
}
