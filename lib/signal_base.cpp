
#include <signal_exciter/signal_base.hpp>
#include <iostream>

boost::mutex Signal_Base::s_mutex_fftw;
//boost::mutex Signal_Base::s_mutex_fso;
boost::mutex* Signal_Base::d_mutex_ptr = &s_mutex_fftw;
size_t Signal_Base::s_indicator = 0;


void
Signal_Base::prototype_augemnt_fractional_delay(
      double interp, double interpd, std::vector<float> &proto,
      double frac_delay, std::vector<float> &taps)
{
  if(d_enable_fractional){
    size_t ntaps = proto.size();
    if(ntaps){
      double fd = std::fmod(frac_delay,1.0)*(2.*interp*interpd);//negative is maintained
      std::vector<float> base_taps = gr::filter::firdes::low_pass_2(1,1,.25,.05,
            120,static_cast<gr::filter::firdes::win_type>(5));
      double alpha = 0.3;
      double center;

      std::vector<float> extended_proto;
      norm_f(proto);
      norm_f(base_taps);
      convolve_2x_1(extended_proto,proto,base_taps);
      norm_f(extended_proto);

      //print_f("Proto", &proto[0], proto.size());
      //print_f("Base", &base_taps[0], base_taps.size());
      //print_f("Extend", &extended_proto[0], extended_proto.size());

      //now have the true prototype
      int L = extended_proto.size();
      int N = L-1;
      if(N%2){//odd
        center = double(N-1)/2;
      }
      else{//even
        center = double(N)/2.;
      }
      double D = center + fd;
      //std::cout << "L(" << L <<"), N("<<N<<"), center("<<center<<"), D("<<D<<")\n";

      int nbins = 2*L;
      gr::fft::fft_complex* to_spec = new gr::fft::fft_complex(nbins,true);
      complexf *ib1 = to_spec->get_inbuf();
      complexf *ob1 = to_spec->get_outbuf();
      memset( ib1, 0, nbins*sizeof(complexf) );
      memset( ob1, 0, nbins*sizeof(complexf) );
      std::vector<double> mag_spec0(nbins);
      std::vector<double> mag_specA(nbins);
      for(size_t idx = 0; idx < L; idx++){
        ib1[idx] = extended_proto[idx];
      }
      //print_fcr("Input_Buff", ib1, nbins);
      to_spec->execute();
      for(size_t idx = 0; idx < nbins; idx++){
        mag_spec0[idx] = double((ob1[idx]*std::conj(ob1[idx])).real());
      }
      //print_d("Mag0", &mag_spec0[0], nbins);
      delete to_spec;
      size_t pointA = int(std::ceil(float(nbins)/2.));
      memcpy( &mag_specA[0], &mag_spec0[pointA], (nbins-pointA)*sizeof(double) );
      memcpy( &mag_specA[pointA], &mag_spec0[0], (pointA)*sizeof(double) );
      //print_d("MagA", &mag_specA[0], nbins);

      double incrmt = (2*M_PI - 1./double(nbins))/(double(nbins-1));
      std::vector<double> omegaA(nbins);
      for(size_t idx = 0; idx < nbins; idx++){
        omegaA[idx] = -M_PI + double(idx)*incrmt;
      }
      //print_d("OmegaA", &omegaA[0], nbins);

      //std::cout << "nbins("<<nbins<<"), alpha("<<alpha<<"), lower("
      //          << (1-2*alpha)/2.<<"),upper("<<(1-(1-2*alpha)/2.)<<")\n";

      int start_at = int( std::max(0., round(double(nbins)*(1-2*alpha)/2.)) );
      int end_at = int( std::min(double(nbins), round(double(nbins)*(1-(1-2*alpha)/2.))) );

      std::vector<double> omegaB( omegaA.begin()+start_at, omegaA.begin()+end_at );
      std::vector<double> mag_specB( mag_specA.begin()+start_at, mag_specA.begin()+end_at );

      int chk = omegaB.size();
      //print_d("MagB", &mag_specB[0], chk);
      //print_d("OmegaB", &omegaB[0], chk);
      std::vector<complexd> Hid(chk);
      std::vector<complexd> ek(chk);
      std::vector<complexd> el(chk);
      std::vector<double> c(chk);
      std::vector<double> s(chk);
      for(size_t idx = 0; idx < chk; idx++){
        Hid[idx] = std::exp(complexd(0.,-omegaB[idx]*D));
      }

      Eigen::MatrixXd P(L,L);
      Eigen::VectorXd p(L);

      double Ckl(0.),pk(0.),scale_w(1./double(chk));
      for(int k = 0; k < L; k++){
        pk = 0.;
        for(size_t w = 0; w < chk; w++){
          ek[w] = std::exp(complexd(0.,-k*omegaB[w]));
          c[w] = std::cos(k*omegaB[w]);
          s[w] = std::sin(k*omegaB[w]);
          pk += mag_specB[w]*(Hid[w].real()*c[w] - Hid[w].imag()*s[w]);
        }
        p(k) = pk*scale_w;
        for(int l = 0; l < L; l++){
          Ckl = 0.;
          for(size_t w = 0; w < chk; w++){
            el[w] = std::exp(complexd(0.,l*omegaB[w]));//conjugate compared to ek
            Ckl += mag_specB[w]*(ek[w]*el[w]).real();
          }
          P(k,l) = Ckl*scale_w;
        }
      }

      char conditions = Eigen::ComputeThinU | Eigen::ComputeThinV;
      Eigen::VectorXd aug = P.jacobiSvd(conditions).solve(p);

      std::vector<double> aug_raw( aug.data(), aug.data() + aug.size() );
      norm_d(aug_raw);
      std::vector<float> aug_raw2( aug_raw.begin(), aug_raw.end() );
      //print_f("Aug", &aug_raw2[0], aug_raw2.size());

      std::vector<float> interp_taps;
      convolve(interp_taps, extended_proto, aug_raw2);
      //print_f("Interp", &interp_taps[0], interp_taps.size());
      //interp_taps should now be a new filter with desired fractional delay
      // @ 2x oversampled

      // going to simplify the taps so that the interp->decim is self-contained.

      std::vector<float> full_taps;
      convolve(full_taps, interp_taps, base_taps);
      //print_f("Full", &full_taps[0], full_taps.size());

      size_t count = (full_taps.size()%2) ? (full_taps.size()+1)/2 : full_taps.size()/2;
      std::vector<float> decim_taps(count,0.);
      for(size_t idx = 0; idx < count; idx++){
        decim_taps[idx] = full_taps[2*idx];
      }
      //print_f("Decim", &decim_taps[0], decim_taps.size());

      size_t half = (count%2) ? (count-1)/2 : count/2;
      size_t htap = (ntaps%2) ? (ntaps-1)/2 : ntaps/2;
      //taps_fd_2[(len(taps_fd_2)-1)/2-(clip_len-1)/2:(len(taps_fd_2)-1)/2+(clip_len-1)/2+1]
      taps = std::vector<float>(ntaps,0.);
      for(size_t idx = 0; idx < ntaps; idx++){
        taps[idx] = decim_taps[half-htap+idx];
      }
      norm_f(taps, interp);
      //print_f("Taps", &taps[0], taps.size());

    }
    else{
      taps = std::vector<float>(0);
    }
  }
  else{
    taps = std::vector<float>( proto.begin(), proto.end() );
    norm_f(taps,interp);
  }
}


void
Signal_Base::convolve_2x_1( std::vector<float> &out,
                            std::vector<float> &in2x,
                            std::vector<float> &keep )
{
  int nt(in2x.size()),nt2(2*nt-1),bt(keep.size()),ct(nt2+bt-1);
  out = std::vector<float>(ct,0.);
  int llidx,lhidx,bplidx,bbhidx,ipidx,uidx,pidx,bidx;
  for(size_t idx = 0; idx < ct; idx++){
    //here we're convolving an interpolated proto with base
    // in2x => [0:nt-1]
    // iin2x[0:2:nt2] = in2x; 0 o.w.
    // keep => [0:bt-1]

    llidx = (idx < std::min(nt2,bt)) ? 0 : idx - (std::min(nt2,bt)-1);
    lhidx = (idx < std::max(nt2,bt)) ? idx+1 : std::max(nt2,bt);

    bplidx = std::max(0,int(idx)-int(bt-1));
    //bphidx = std::min(nt2-1,idx);
    //bblidx = std::min(std::max(0,int(idx)-int(nt2-1)),std::max(bt-1,nt2-1));
    bbhidx = std::min(int(idx),bt-1);

    for(size_t ind = 0; ind < (lhidx-llidx); ind++){
      ipidx = bplidx + ind;//interpolated prototype index
      uidx = 1-ipidx%2;//used index?
      if(uidx){
        pidx = ipidx/2;//prototype index
        bidx = bbhidx - ind;//base index
        out[idx] += in2x[pidx] * keep[bidx];
      }
    }
  }
}

void
Signal_Base::convolve(std::vector<float> &out,
                      std::vector<float> &a,
                      std::vector<float> &b)
{
  int at(a.size()),bt(b.size()),ct(at+bt-1);
  out = std::vector<float>(ct,0.);
  int llidx,lhidx,bplidx,bbhidx,pidx,bidx;
  for(size_t idx = 0; idx < ct; idx++){

    llidx = (idx < std::min(at,bt)) ? 0 : idx - (std::min(at,bt)-1);
    lhidx = (idx < std::max(at,bt)) ? idx+1 : std::max(at,bt);

    bplidx = std::max(0,int(idx)-int(bt-1));
    //bphidx = std::min(nt2-1,idx);
    //bblidx = std::min(std::max(0,int(idx)-int(at-1)),std::max(bt-1,at-1));
    bbhidx = std::min(int(idx),bt-1);

    for(size_t ind = 0; ind < (lhidx-llidx); ind++){
      pidx = bplidx + ind;//prototype index
      bidx = bbhidx - ind;//base index
      out[idx] += a[pidx] * b[bidx];
    }
  }
}

void
Signal_Base::norm_f(std::vector<float> &to_norm, float scale)
{
  double wt = 0.;
  for(size_t idx = 0; idx < to_norm.size(); idx++){
    wt += to_norm[idx]*to_norm[idx];
  }
  if(wt) wt = std::sqrt(scale/wt);
  else wt = 1.;
  for(size_t idx = 0; idx < to_norm.size(); idx++){
    to_norm[idx] *= wt;
  }
}

void
Signal_Base::norm_d(std::vector<double> &to_norm, double scale)
{
  double wt = 0.;
  for(size_t idx = 0; idx < to_norm.size(); idx++){
    wt += to_norm[idx]*to_norm[idx];
  }
  if(wt) wt = std::sqrt(scale/wt);
  else wt = 1.;
  for(size_t idx = 0; idx < to_norm.size(); idx++){
    to_norm[idx] *= wt;
  }
}

void
Signal_Base::print_f(std::string Name, float* vec, size_t len)
{
  std::cout << d_indicator << ": ";
  if(len){
    std::cout << Name << " = [" << vec[0];
    for(size_t idx = 1; idx < len; idx++){
      std::cout << ", " << vec[idx];
    }
    std::cout << "];\n";
  }
  else{
    std::cout << Name << " = [];\n";
  }
}

void
Signal_Base::print_fcr(std::string Name, complexf* vec, size_t len)
{
  std::cout << d_indicator << ": ";
  if(len){
    std::cout << Name << " = [" << vec[0].real();
    for(size_t idx = 1; idx < len; idx++){
      std::cout << ", " << vec[idx].real();
    }
    std::cout << "];\n";
  }
  else{
    std::cout << Name << " = [];\n";
  }
}

void
Signal_Base::print_d(std::string Name, double* vec, size_t len)
{
  std::cout << d_indicator << ": ";
  if(len){
    std::cout << Name << " = [" << vec[0];
    for(size_t idx = 1; idx < len; idx++){
      std::cout << ", " << vec[idx];
    }
    std::cout << "];\n";
  }
  else{
    std::cout << Name << " = [];\n";
  }
}
