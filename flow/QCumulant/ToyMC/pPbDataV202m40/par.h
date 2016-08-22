#include "TF1.h"
TF1 *V2vsPt = new TF1("V2vsPt","((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);
const int nFileAll = 300;
const TString dir = "/lio/lfs/cms/store/user/qixu/flow/NewSTEG/pPbDataV202m40";
const TString outdir = "/scratch/xuq7/flow/QCumulant/pPbDataV202m40";
const double ptmin = 0.3;
const double ptmax = 6.0; 
const double ptbinv[]={0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
const double etabinv[]={-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4};
const int nptv= 9;
const int netav= 12;
const int trkbin[]={41,39};
const int nbin = 1;
const int nsamples=50;
const double inV2=0.0502178*0.4;
const TString method = "recurrence";
const int maxN = 8;
const bool useWeight = false;
