/* -*- c++ -*- */
/*
 *     gr-signal_exciter  Copyright (C) 2015-2017  Bill Clark
 *     This program comes with ABSOLUTELY NO WARRANTY;
 *     This is free software, and you are welcome to redistribute it
 *     under certain conditions;
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "whiten_and_compress_block_impl.h"
#include <stdio.h>

namespace gr {
  namespace signal_exciter {

    whiten_and_compress_block::sptr
    whiten_and_compress_block::make(size_t blocksize, double snr_db, int seed)
    {
      return gnuradio::get_initial_sptr
        (new whiten_and_compress_block_impl(blocksize, snr_db, seed));
    }

    /*
     * The private constructor
     */
    whiten_and_compress_block_impl::whiten_and_compress_block_impl(size_t blocksize, double snr_db, int seed)
      : gr::sync_interpolator("whiten_and_compress_block",
              gr::io_signature::make(1, 1, sizeof(complexf)),
              gr::io_signature::make(1, 1, sizeof(int16_t)),
              2),
        d_block_size(blocksize),
        d_snr_db(snr_db)
    {
      set_seed(seed);
      d_rng = new gr::random(d_seed, 0, 1);

      d_noise_var = std::pow(10.,-(d_snr_db/10.));

      set_output_multiple(2*d_block_size);
      const int align = volk_get_alignment();
      d_sig = (complexf *) volk_malloc(d_block_size*sizeof(complexf), align);
      d_nse = (complexf *) volk_malloc(d_block_size*sizeof(complexf), align);

      d_scale_base = 32767.; //max val of signed 16bit (2^15-1)

    }

    /*
     * Our virtual destructor.
     */
    whiten_and_compress_block_impl::~whiten_and_compress_block_impl()
    {
      delete d_rng;
      volk_free(d_sig);
      volk_free(d_nse);
    }

    int
    whiten_and_compress_block_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const complexf *in = (const complexf *) input_items[0];
      int16_t *out = (int16_t *) output_items[0];
      float *sig = (float *) d_sig;
      float *nse = (float *) d_nse;

      size_t ninput_items = noutput_items/interpolation();
      size_t iters = ninput_items/(d_block_size);
      size_t oo = 0;

      float span = float(d_block_size);
      float scale;

      //printf("Generating %d items, from %lu items, over %lu iterations\n",noutput_items,ninput_items,iters);
      //printf("The noise variance is given as: %1.3e\n",d_noise_var);

      for(size_t cnt = 0; cnt < iters; cnt++){
        memcpy( &d_sig[0], &in[d_block_size*cnt], d_block_size*sizeof(complexf) );
        for(size_t idx = 0; idx < d_block_size; idx++){
          d_nse[idx] = d_rng->rayleigh_complex();
        }
        float sig_w(0.),nse_w(0.);
        for(size_t idx = 0; idx < d_block_size; idx++){
          sig_w += (d_sig[idx] * std::conj(d_sig[idx])).real();
          nse_w += (d_nse[idx] * std::conj(d_nse[idx])).real();
        }

        scale = sig_w/span;
        if(scale) scale = 1./std::sqrt(scale);//might be all zeros, don't div by 0
        volk_32f_s32f_multiply_32f( sig, sig, scale, 2*d_block_size );


        scale = nse_w/span;
        scale = std::sqrt(d_noise_var/scale);
        volk_32f_s32f_multiply_32f( nse, nse, scale, 2*d_block_size );

/*        sig_w = 0;
        nse_w = 0;
        for(size_t idx = 0; idx < d_block_size; idx++){
          sig_w += (d_sig[idx] * std::conj(d_sig[idx])).real();
          nse_w += (d_nse[idx] * std::conj(d_nse[idx])).real();
        }
        sig_w = sig_w/span;
        nse_w = nse_w/span;
        printf("SIG(%1.5e), NOISE(%1.5e)\t",sig_w,nse_w);
        printf("SNR(%1.5e), SNRdB(%02.1f)\n",sig_w/nse_w,10.*log(sig_w/nse_w)/log(10.));*/


        volk_32f_x2_add_32f( sig, sig, nse, 2*d_block_size );

        float max_val(0.);
        for(size_t idx = 0; idx < 2*d_block_size; idx++){
          if(std::abs(sig[idx]) > max_val) max_val = std::abs(sig[idx]);
        }

        scale = float(d_scale_base)/max_val;
        //printf("The scale is: %1.16f\n",scale);
        volk_32f_s32f_multiply_32f( sig, sig, scale, 2*d_block_size );
        
        for(size_t idx = 0; idx < 2*d_block_size; idx++){
          out[oo+idx] = (int16_t) roundf(sig[idx]);
        }

        /*printf("The in/out relationship\ntemp = [ (%1.3e+%1.3ej), (%1.3e+%1.3ej), %05d, %05d",
          in[d_block_size*cnt].real(),
          in[d_block_size*cnt].imag(),
          sig[0],sig[1],
          out[oo],out[oo+1]);
        for(size_t idx = 1; idx < d_block_size; idx++){
          printf("\n  (%1.3e+%1.3ej), (%1.3e+%1.3ej), %05d, %05d",
            in[d_block_size*cnt+idx].real(),
            in[d_block_size*cnt+idx].imag(),
            sig[2*idx],sig[2*idx+1],
            out[oo+2*idx],out[oo+2*idx+1]);
        }
        printf("];\n");*/

        oo += 2*d_block_size;
      }

      //printf("oo(%lu)\n",oo);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace signal_exciter */
} /* namespace gr */

