const int ncent = 1;
const int nhar = 3;
const int nbbcz = 10;
const int nsub = 3;
const int npt = 20;
const int ncorr = 3;
const float ptbin[npt+1] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0};
const int nsample = 10;

//for QCumulant
const double ptmin = 0.3;
const double ptmax = 3.0; 
const double ptbinv[]={0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
const double etabinv[]={-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4};
const int nptv = 9;
const int netav = 12;
const int nsamples = 50;
const TString method = "recurrence";
const int maxH = 8;
const bool doLoops = false;
