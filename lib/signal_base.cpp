
#include <signal_exciter/signal_base.hpp>


boost::mutex Signal_Base::s_mutex_fftw;
boost::mutex* Signal_Base::d_mutex_ptr = &s_mutex_fftw;


void
Signal_Base::prototype_augemnt_fractional_delay(
      double interp, std::vector<float> &proto,
      double frac_delay, std::vector<float> &taps,
      std::vector<float> &extended_proto)
{
  size_t ntaps = proto.size();
  if(ntaps){
    double alpha = 0.3;
    double center;
    int offset;
    double fd = std::fmod(frac_delay,1.0)*(2.*interp);//negative is maintained
    std::vector<float> base_taps = gr::filter::firdes::low_pass_2(1,1,.25,.05,
          120,static_cast<gr::filter::firdes::win_type>(5));
    double wpt(0.),wbt(0.);
    for(size_t idx = 0; idx < ntaps; idx++){
      wpt += proto[idx]*proto[idx];
    }
    for(size_t idx = 0; idx < base_taps.size(); idx++){
      wbt += base_taps[idx]*base_taps[idx];
    }

    if(wpt) wpt = 1./std::sqrt(wpt);
    else wpt = 1;
    if(wbt) wbt = 1./std::sqrt(wbt);
    else wbt = 1;

    for(size_t idx = 0; idx < ntaps; idx++){
      proto[idx] *= wpt;
    }
    for(size_t idx = 0; idx < base_taps.size(); idx++){
      base_taps[idx] *= wbt;
    }

    size_t bplidx,bphidx,bblidx,bbhidx,llidx,lhidx;
    size_t ipidx,pidx,bidx,uidx;
    size_t nt2(2*ntaps-1),bt(base_taps.size());
    extended_proto = std::vector<float>(nt2+bt-1,0.);
    for(size_t idx = 0; idx < extended_proto.size(); idx++){
      //here we're convolving an interpolated proto with base
      // proto => [0:ntaps-1]
      // iproto[0:2:2*ntaps-1] = proto; 0 o.w.
      // base => [0:base_taps.size()]

      llidx = (idx < std::min(nt2,bt)) ? 0 : idx - (std::min(nt2,bt)-1);
      lhidx = (idx < std::max(nt2,bt)) ? idx+1 : std::max(nt2,bt);

      bplidx = std::max(0,int(idx)-int(bt-1));
      //bphidx = std::min(nt2-1,idx);
      //bblidx = std::min(std::max(0,int(idx)-int(nt2-1)),std::max(bt-1,nt2-1));
      bbhidx = std::min(idx,bt-1);

      for(size_t ind = 0; ind < (lhidx-llidx); ind++){
        ipidx = bplidx + ind;//interpolated prototype index
        uidx = 1-ipidx%2;//used index?
        if(uidx){
          pidx = ipidx/2;//prototype index
          bidx = bbhidx - ind;//base index
          extended_proto[idx] += proto[pidx] * base_taps[bidx];
        }
      }
    }

    double wet(0.);
    for(size_t idx = 0; idx < extended_proto.size(); idx++){
      wet += extended_proto[idx]*extended_proto[idx];
    }
    if(wet) wet = 1./std::sqrt(wet);
    else wet = 1.;
    for(size_t idx = 0; idx < extended_proto.size(); idx++){
      extended_proto[idx] *= wet;
    }
    //now have the true prototype
    int L = extended_proto.size();
    int N = L-1;
    if(N%2){//odd
      center = double(N-1)/2;
      offset = (N-1)/2;
    }
    else{//even
      center = double(N)/2.;
      offset = N/2;
    }
    double D = center + fd;

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
    to_spec->execute();
    for(size_t idx = 0; idx < nbins; idx++){
      mag_spec0[idx] = double((ob1[idx]*std::conj(ob1[idx])).real());
    }
    size_t pointA = int(std::ceil(float(nbins)/2.));
    memcpy( &mag_specA[0], &mag_spec0[pointA], (nbins-pointA)*sizeof(double) );
    memcpy( &mag_specA[pointA], &mag_spec0[0], (pointA)*sizeof(double) );

    double incrmt = (2*M_PI - 1./double(nbins))/(double(nbins-1));
    std::vector<double> omegaA(nbins);
    for(size_t idx = 0; idx < nbins; idx++){
      omegaA[idx] = -M_PI + double(idx)*incrmt;
    }

    int start_at = int( std::max(0., round(double(nbins)*(1-2*alpha)/2.)) );
    int end_at = int( std::min(double(nbins), round(double(nbins)*(1-2*alpha)/2.)) );

    std::vector<double> omegaB( omegaA.begin()+start_at, omegaA.begin()+end_at );
    std::vector<double> mag_specB( mag_specA.begin()+start_at, mag_specA.begin()+end_at );

    int chk = omegaB.size();
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
    double wat(0.);
    for(size_t idx = 0; idx < aug_raw.size(); idx++){
      wat += aug_raw[idx]*aug_raw[idx];
    }
    if(wat) wat = 1./std::sqrt(wat);
    else wat = 1.;
    for(size_t idx = 0; idx < aug_raw.size(); idx++){
      aug_raw[idx] *= wat;
    }

    taps = std::vector<float>(L,0.);
    for(size_t cidx = 0; cidx < L; cidx++){
      //here we're convolving an extended_proto with LMS augmentor
      size_t idx = cidx + offset;//only compute center set of taps

      llidx = (idx < L) ? 0 : idx - N;
      lhidx = (idx < L) ? idx+1 : L;

      bplidx = std::max(0,int(idx)-N);
      //bphidx = std::min(N,idx);
      //bblidx = std::min(std::max(0,int(idx)-N),N);
      bbhidx = std::min(int(idx),N);

      for(size_t ind = 0; ind < (lhidx-llidx); ind++){
        pidx = bplidx + ind;//extended prototype index
        bidx = bbhidx - ind;//augmentor index
        taps[idx] += extended_proto[pidx] * aug_raw[bidx];
      }
    }
    //taps should now be a new filter with desired fractional delay
  }
  else{
    taps = std::vector<float>(0);
  }
}
