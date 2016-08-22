double lamy_wgt( double y, int weight_choice ) {

  if( weight_choice == 1 ) {//D6T
    const int nbins = 23;
    const Double_t y_weight[nbins] = {1.051177, 1.049749, 1.06767, 1.022956, 1.054922, 1.049431, 1.027359, 1.021834, 1.02841, 1.037216, 1.058857, 1.026179, 1.057839, 1.013954, 0.9920087, 0.9752364, 0.9693014, 0.9585418, 0.9669884, 0.8867812, 0.8287872, 0.7438514, 0.710441};
    int ndx = (int) 10.0*y;
    if (ndx >= nbins) ndx = nbins-1;
    return y_weight[ndx];
  }

  else if( weight_choice == 2 ) {//pythia8
    const int nbins = 23;
    const Double_t y_weight[nbins] = {1.028947, 1.031015, 1.032598, 0.9730167, 1.018098, 1.014717, 1.002446, 0.9856325, 1.035051, 1.009572, 1.020354, 1.001895, 1.024646, 0.98586, 0.9812227, 0.9871389, 1.004868, 0.9943691, 1.010525, 0.9549695, 0.9412455, 0.9044774, 0.8537855};
    int ndx = (int) 10.0*y;
    if (ndx >= nbins) ndx = nbins-1;
    return y_weight[ndx];
  }

  else if( weight_choice == 3 ) {//P0
    const int nbins = 23;
    const Double_t y_weight[nbins] = {1.008809, 0.9989802, 0.9748529, 0.976692, 1.004819, 0.9957665, 1.000636, 0.9928182, 1.016604, 0.9792475, 1.037794, 0.9859791, 1.040842, 0.9764901, 0.9958794, 1.003388, 1.006224, 1.012714, 1.036967, 1.011207, 1.00944, 0.9400557, 0.9380598};
    int ndx = (int) 10.0*y;
    if (ndx >= nbins) ndx = nbins-1;
    return y_weight[ndx];
  }

  else if( weight_choice == 11 ) {//900 GeV D6T
    const int nbins = 23;
    const Double_t y_weight[nbins] = {1.037445, 1.117442, 1.035941, 1.010867, 1.050708, 1.019198, 1.013151, 1.028558, 0.9920126, 1.045949, 0.9945654, 0.953168, 0.9643439, 0.9838724, 0.9611377, 0.977955, 0.97627, 0.9681172, 0.923717, 0.9609616, 0.9235286, 0.8530482, 0.9754884};
    int ndx = (int) 10.0*y;
    if (ndx >= nbins) ndx = nbins-1;
    return y_weight[ndx];
  }

  else if( weight_choice == 1002 ) {//pythia8 Lambda (not Lambdabar)
    const int nbins = 23;
    const Double_t y_weight[nbins] = {1.22193, 1.191453, 1.210649, 1.133128, 1.230681, 1.19997, 1.174719, 1.206048, 1.241441, 1.272074, 1.242612, 1.217095, 1.245536, 1.087139, 1.124356, 1.097888, 1.100562, 1.084934, 1.089828, 1.057444, 0.9827095, 0.9870367, 0.9032874};
    int ndx = (int) 10.0*y;
    if (ndx >= nbins) ndx = nbins-1;
    return y_weight[ndx];
  }

  else if( weight_choice == 2002 ) {//pythia8 Lambdabar (not Lambda)
    const int nbins = 23;
    const Double_t y_weight[nbins] = {1.036643, 1.034983, 1.045942, 0.9411083, 1.018528, 1.04661, 0.9916014, 0.9688672, 1.039244, 1.029095, 1.041136, 0.9956324, 1.024812, 0.9827136, 0.9698727, 0.9982474, 1.015968, 1.008592, 1.016995, 0.9294032, 0.8798677, 0.9357154, 0.8509101};
    int ndx = (int) 10.0*y;
    if (ndx >= nbins) ndx = nbins-1;
    return y_weight[ndx];
  }

  return 1.0;

}
