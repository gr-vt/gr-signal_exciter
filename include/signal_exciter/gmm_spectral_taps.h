/* 
 * Copyright 2017 Bill Clark.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */
#ifndef INCLUDED_GMM_SPECTRAL_TAPS_H
#define INCLUDED_GMM_SPECTRAL_TAPS_H

#include <stdlib.h>
#include <vector>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <fftw3.h>

/*
Working prototype model
fs = 8000;
mu = [200,450,1500,2500,3500]/fs;
sig = [50,75,750,250,250]/fs;
rho = [.29,.60,.08,.01,.02];
*/

class GMM_Spectral_Taps
{
 private:
  size_t              d_tap_count;
  float               d_samp_rate;
  bool                d_made;
  std::vector<float>  d_mu;
  std::vector<float>  d_sigma;
  std::vector<float>  d_weight;

  std::vector<float>  d_freq;

  fftwf_plan d_ifft;
  fftwf_complex* d_fft_in;
  fftwf_complex* d_fft_out;

  void destroy_fft()
  {
    if(d_made){
      fftwf_free( d_fft_in );
      fftwf_free( d_fft_out );
      fftwf_destroy_plan( d_ifft );
      d_made = false;
    }
  }

  void build_fft()
  {
    d_fft_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*d_tap_count);
    d_fft_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*d_tap_count);
    d_ifft = fftwf_plan_dft_1d( d_tap_count, d_fft_in, d_fft_out, FFTW_BACKWARD, FFTW_ESTIMATE );
    d_made = true;
  }

 public:
  GMM_Spectral_Taps()
  {
    d_tap_count = 0;
    d_mu = std::vector<float>(0);
    d_sigma = std::vector<float>(0);
    d_weight = std::vector<float>(0);
    d_samp_rate = 1.;
    d_freq = std::vector<float>(0);
    d_made = false;
  }

  GMM_Spectral_Taps(size_t components, float* mus, float* sigmas, float* weights, float samp_rate, size_t tap_count)
  {
    d_mu      = std::vector<float>(mus, mus+components);
    d_sigma   = std::vector<float>(sigmas, sigmas+components);
    d_weight  = std::vector<float>(weights, weights+components);
    d_tap_count = tap_count;

    build_fft();

    float o = (1./float(d_tap_count))*float(d_tap_count%2);
    float p = std::pow(2.,float(d_tap_count%2));
    d_freq = std::vector<float>(d_tap_count);

    float a = -.5+o/p;
    float b = .5-1./float(d_tap_count)/p;
    float c = (b-a)/(float(d_tap_count-1));
    for(size_t idx = 0; idx < d_tap_count; idx++){
      d_freq[idx] = a + idx*c;
    }
  }

  ~GMM_Spectral_Taps()
  {
    if(d_made){
      destroy_fft();
    }
  }

  void get_taps(std::vector<float>& taps)
  {
    if(d_made){
      taps = std::vector<float>(d_tap_count,0.);
      memset( d_fft_in, 0, sizeof(fftwf_complex)*d_tap_count );
      memset( d_fft_out, 0, sizeof(fftwf_complex)*d_tap_count );
      size_t offset = (size_t)std::ceil(float(d_tap_count)/2.);//ifftshift
      
      for(size_t kk = 0; kk < d_mu.size(); kk++){
        for(size_t idx = 0; idx < d_tap_count; idx++){
          d_fft_in[(idx-offset)%d_tap_count][0] += d_weight[kk]/2. * (1./(std::sqrt(2*M_PI)*(d_sigma[kk]/d_samp_rate))) * std::exp(-std::pow((d_freq[idx]-(d_mu[idx]/d_samp_rate)),2.)/(2.*d_sigma[kk]*d_sigma[kk]/d_samp_rate/d_samp_rate));
          d_fft_in[(idx-offset)%d_tap_count][0] += d_weight[kk]/2. * (1./(std::sqrt(2*M_PI)*(d_sigma[kk]/d_samp_rate))) * std::exp(-std::pow((d_freq[idx]+(d_mu[idx]/d_samp_rate)),2.)/(2.*d_sigma[kk]*d_sigma[kk]/d_samp_rate/d_samp_rate));
        }
      }
      fftwf_execute( d_ifft );

      offset = (size_t)std::floor(float(d_tap_count)/2.);//fftshift
      for(size_t idx = 0; idx < d_tap_count; idx++){
        taps[idx] = d_fft_out[(idx-offset)%d_tap_count][0];
      }

      float pwr_chk = 0;
      for(size_t idx = 0; idx < d_tap_count; idx++){
        pwr_chk += taps[idx]*taps[idx];
      }
      float scale = 1/std::sqrt(pwr_chk);
      for(size_t idx = 0; idx < d_tap_count; idx++){
        taps[idx] *= scale;
      }
    }
    else{
      taps = std::vector<float>(0);
    }
  }

  void set_params(size_t components, float* mus, float* sigmas, float* weights, float samp_rate, size_t tap_count)
  {
    d_mu      = std::vector<float>(mus, mus+components);
    d_sigma   = std::vector<float>(sigmas, sigmas+components);
    d_weight  = std::vector<float>(weights, weights+components);
    d_tap_count = tap_count;

    build_fft();

    float o = (1./float(d_tap_count))*float(d_tap_count%2);
    float p = std::pow(2.,float(d_tap_count%2));
    d_freq = std::vector<float>(d_tap_count);

    float a = -.5+o/p;
    float b = .5-1./float(d_tap_count)/p;
    float c = (b-a)/(float(d_tap_count-1));
    for(size_t idx = 0; idx < d_tap_count; idx++){
      d_freq[idx] = a + idx*c;
    }
  }
};


#endif //INCLUDED_GMM_SPECTRAL_TAPS_H
