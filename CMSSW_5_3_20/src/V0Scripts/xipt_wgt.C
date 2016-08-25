// weight_choice < 10 are 7 TeV.  weight > choice > 10 are 900 GeV
// weight_choice 1,11 are D6T, 2,12 are pythia8, 3,13 are P0
double xipt_wgt( double pt, int weight_choice ) {

  if( !weight_choice ) return 1.;

  else if( weight_choice == 1 ) {//D6T
    const int nbins = 10;
    const Double_t xAxis1[nbins+1] = {1, 1.4, 1.6, 1.8, 2, 2.2, 2.6, 3, 3.4, 4, 6};
    const Double_t pt_weight[nbins] = {1.660826, 1.464409, 1.329119, 1.406882, 1.122656, 0.9803696, 0.7671991, 0.7260265, 0.561658, 0.5607015};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( pt < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return pt_weight[ndx];
  }

  else if( weight_choice == 2 ) {//pythia8
    return 1;
  }

  else if( weight_choice == 3 ) {//P0
    return 1;
  }

  else if( weight_choice == 11 ) {//900 GeV D6T
    const int nbins = 10;
    const Double_t xAxis1[nbins+1] = {1, 1.4, 1.6, 1.8, 2, 2.2, 2.6, 3, 3.4, 4, 6};
    const Double_t pt_weight[nbins] = {1.246124, 0.9699563, 0.9253977, 1.141979, 0.734954, 1.162812, 0.9290018, 2.069333, 0.8191425, 0.2775709};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( pt < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return pt_weight[ndx];
  }


}
