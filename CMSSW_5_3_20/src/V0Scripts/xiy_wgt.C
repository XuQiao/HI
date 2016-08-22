// weight_choice < 10 are 7 TeV.  weight > choice > 10 are 900 GeV
// weight_choice 1,11 are D6T, 2,12 are pythia8, 3,13 are P0
double xiy_wgt( double y, int weight_choice ) {

  if( !weight_choice ) return 1.;

  else if( weight_choice == 1 ) {//D6T
    const int nbins = 11;
    const Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.3};
    const Double_t y_weight[nbins] = {1.146013, 0.9239685, 1.111753, 0.9692743, 1.088724, 1.04348, 0.9885751, 0.884877, 1.123645, 0.9240211, 0.8222993};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( y < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return y_weight[ndx];
  }

  else if( weight_choice == 2 ) {//pythia8
    return 1;
  }

  else if( weight_choice == 3 ) {//P0
    return 1;
  }

  else if( weight_choice == 11 ) {//900 GeV D6T
    const int nbins = 11;
    const Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.3};
    const Double_t y_weight[nbins] = {1.261148, 1.346277, 2.059355, 1.133999, 0.779956, 1.967053, 1.033872, 0.7615587, 0.7697292, 0.7707394, 0.8239518};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( y < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return y_weight[ndx];
  }


}
