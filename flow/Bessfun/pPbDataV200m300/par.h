#include "TF1.h"

//TF1 *PtDistr = new TF1("PtDistr","[0]*x*TMath::Power(1+(sqrt(x*x+[1]*[1])-[1]/[2]),-[3])",0.2,10);
//	PtDistr->SetParameters(118.836,-0.335972,0.759243,118.836);     //Real data fit with Tsallis
//TF1 *V2vsPt = new TF1("V2vsPt","((x/5)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);    //Toy V~0.06
  TF1 *V2vsPt = new TF1("V2vsPt","((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);
const int nFileAll = 300;
const TString dir = "/lio/lfs/cms/store/user/qixu/flow/STEG/pPbDataV200m300";
const double ptmin = 0.3;
const double ptmax = 6.0; 
//const double ptbinV[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
//const double ptbinV[]={0.3,6.0};
//const int nptV= 1;
const double ptbinv[]={0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
const double etabinv[]={-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4};
const int nptv= 9;
const int netav= 12;
const int trkbin[]={301,299};
const int nbin = 1;
//const int trkbin[]={185,180,175,170,165,160,155,150};
const int nn=2;
const int mm=1;
const int isSum=1;
const double j01=2.404826;
const int ntheta = 1;
const int nstepr=200;
const bool isSimple=1;
const int nsamples=50;
const double inV2=0;
