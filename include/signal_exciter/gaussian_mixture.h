/* 
 * Copyright 2016 Bill Clark.
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
#ifndef INCLUDED_GAUSSIAN_MIXTURE_H
#define INCLUDED_GAUSSIAN_MIXTURE_H

#include <stdlib.h>
#include <vector>
#include <string.h>
#include <math.h>
#include <time.h>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/random_device.hpp>
#include <stdio.h>
#include <gnuradio/filter/fir_filter.h>

class Gaussian_Mixture
{
  private:
    int d_seed;

    boost::random_device d_rd;
    boost::mt19937 d_rng;

  public:

    void set_seed(int seed)
    {
      //printf("Setting the seed\n");
      d_seed = seed;
      if(d_seed < 0) d_seed = d_rd();
      d_rng.seed(d_seed);
    }

  private:
    float d_sigma1;
    float d_sigma2;
    float d_fmax;
    float d_thresh;
    
    std::vector<float> d_1;
    std::vector<float> d_2;
    std::vector<float> d_3;
    std::vector<float> d_4;

    std::vector<float> d_hist;
    std::vector<float> d_taps;
    float d_max_found;

    size_t d_thr;

    gr::filter::kernel::fir_filter_fff* d_fir;


    //boost::mt19937 *d_rng;
    //boost::normal_distribution<> d_dist1;
    //boost::normal_distribution<> d_dist2;
    //boost::variate_generator< boost::mt19937, boost::normal_distribution<> >* d_gen1;
    //boost::variate_generator< boost::mt19937, boost::normal_distribution<> >* d_gen2;
    boost::random::variate_generator< boost::mt19937&, boost::normal_distribution<float> > d_gen1;
    boost::random::variate_generator< boost::mt19937&, boost::normal_distribution<float> > d_gen2;

    void gen_gauss1(float* output, size_t count)
    {
      for(size_t idx = 0; idx < count; idx++){
        //output[idx] = (*d_gen1)();
        output[idx] = d_sigma1 * d_gen1();
      }
    }

    void gen_gauss2(float* output, size_t count)
    {
      for(size_t idx = 0; idx < count; idx++){
        //output[idx] = (*d_gen2)();
        output[idx] = d_sigma2 * d_gen2();
      }
    }

    void sel_gauss(float* output, size_t count)
    {
      size_t count1(0),count2(0),c1(0),c2(0);
      for(size_t idx = 0; idx < count; idx++){
        if(rand() > d_thr){
          count1++;
          output[idx] = 1.;
        }
      }
      count2 = count-count1;
      d_1 = std::vector<float>(count1,0.);
      d_2 = std::vector<float>(count2,0.);
      gen_gauss1( &d_1[0], count1 );
      gen_gauss2( &d_2[0], count2 );
      for(size_t idx = 0; idx < count; idx++){
        if(output[idx]!=0.){
          output[idx] = d_1[c1++];
        }
        else{
          output[idx] = d_2[c2++];
        }
      }
    }

    void filt_gauss(float* output, size_t count)
    {
      d_3 = std::vector<float>(count+d_hist.size(),0.);

      sel_gauss( &d_3[d_hist.size()], count );
      memcpy( &d_3[0], &d_hist[0], d_hist.size()*sizeof(float) );
      memcpy( &d_hist[0], &d_3[d_3.size()-d_hist.size()], d_hist.size()*sizeof(float) );

      memset( &output[0], 0., count*sizeof(float) );


      /*for(size_t idx = 0; idx < count; idx++){
        output[idx] = 0.;
        for(int ind_h(0), ind_m(d_hist.size()+idx); (ind_m >= 0) && (ind_h < (int)d_taps.size()); ind_m--, ind_h++){
          output[idx] += d_taps[ind_h]*d_3[ind_m];
        }
        if(output[idx] > d_max_found) d_max_found = output[idx];
        if(-output[idx] > d_max_found) d_max_found = -output[idx];
      }*/

      d_fir->filterN( &output[0], &d_3[0], count );
      for(size_t idx = 0; idx < count; idx++){
      //  if(output[idx] > d_max_found) d_max_found = output[idx];
      //  if(-output[idx] > d_max_found) d_max_found = output[idx];
        if(output[idx] > 1.) output[idx] = 1;
        if(-output[idx] > 1.) output[idx] = -1;
      }


      //for(size_t idx = 0; idx < count; idx++){
      //  output[idx] = output[idx]/d_max_found;
      //  if(isnan(output[idx])) output[idx] = 0.;
      //}
    }

    void generate_taps()
    {
      d_taps = std::vector<float>(51,0.);
      d_hist = std::vector<float>(50,0.);

      d_max_found = 0.;
      std::vector<float> d_window(51,0.);
      for(size_t idx = 0; idx < 51; idx++)
      {
        if((idx-25) == 0)
        {
          d_taps[idx] = d_fmax;
        }
        else
        {
          float x = (float(idx)-25)*d_fmax;
          d_taps[idx] = d_fmax*sin(M_PI*x)/(M_PI*x);
        }
        if(!((idx==0)||(idx==50)))
        {
          d_window[idx] = 0.42 - 0.5*cos((2*M_PI*float(idx))/50) + 0.08*cos((4*M_PI*float(idx))/50);
        }
        d_taps[idx] = d_taps[idx]*d_window[idx];
      }

      d_fir = new gr::filter::kernel::fir_filter_fff(1, d_taps);
    }

  public:
    Gaussian_Mixture()
      : d_rng(),
        d_gen1(d_rng, boost::normal_distribution<float>(0.,1.0)),
        d_gen2(d_rng, boost::normal_distribution<float>(0.,1.0))
    {
    }

    Gaussian_Mixture(float sigma_1, float sigma_2, float fmax, float thresh, int seed=-1)
      : d_sigma1(sigma_1),
        d_sigma2(sigma_2),
        d_fmax(fmax),
        d_thresh(thresh),
        d_rng(),
        d_gen1(d_rng, boost::normal_distribution<float>(0.,1.0)),
        d_gen2(d_rng, boost::normal_distribution<float>(0.,1.0))
    {
      //printf("GM: Making object.\n");
      set_seed(seed);
      //printf("GM: Set seed.\n");
      if(d_sigma1 < 0) d_sigma1 = -d_sigma1;
      if(d_sigma2 < 0) d_sigma2 = -d_sigma2;
      if(d_thresh < 0) d_thresh = 0.;
      if(d_thresh > 1) d_thresh = 1.;
      //printf("GM: Set params.\n");

      d_thr = size_t(float(RAND_MAX)*d_thresh);

      generate_taps();
      //printf("GM: Made taps.\n");
    }

    ~Gaussian_Mixture()
    {
      delete d_fir;
    }

    Gaussian_Mixture & operator=(const Gaussian_Mixture& rhs){
      if(this != &rhs){
        set_seed(rhs.d_seed);
        d_sigma1 = rhs.d_sigma1;
        d_sigma2 = rhs.d_sigma2;
        d_thresh = rhs.d_thresh;
        d_fmax = rhs.d_fmax;

        d_thr = size_t(float(RAND_MAX)*d_thresh);
        generate_taps();
        d_max_found = rhs.d_max_found;
      }
      return *this;
    }

    void get_message(float* output, size_t message_count)
    {
      //printf("GM: Making Message.\n");
      filt_gauss(output,message_count);
      /*printf("GM: Made Message.\n[%1.3f",output[0]);
      for(size_t idx=1; idx < message_count; idx++){
        printf(", %1.3e",output[idx]);
      }
      printf("]\n%1.3e, %1.3e\n",d_sigma1,d_sigma2);*/
    }

};
















#endif //INCLUDED_GAUSSIAN_MIXTURE_H
