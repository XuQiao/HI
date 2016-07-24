#include "TF1.h"

//TF1 *PtDistr = new TF1("PtDistr","[0]*x*TMath::Power(1+(sqrt(x*x+[1]*[1])-[1]/[2]),-[3])",0.2,10);
//	PtDistr->SetParameters(118.836,-0.335972,0.759243,118.836);     //Real data fit with Tsallis
TF1 *V2vsPt = new TF1("V2vsPt","((x/5)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);    //Toy V~0.06
const int nFileAll = 300;
const TString dir = "/lio/lfs/cms/store/user/qixu/flow/STEG/ToyV206m150/";
const double ptmin = 0.2;
const double ptmax = 10.0; 
//const double ptbinV[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
const double ptbinV[]={0.2,10.0};
const int nptV= 1;
const double ptbinv[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
const int nptv= 13;
const int trkbin[]={160,140};
const int nbin = 1;
//const int trkbin[]={185,180,175,170,165,160,155,150};
const int nn=2;
const int mm=1;
const int isSum=0;
const double j01=2.404826;
const int ntheta = 5;
const int nstepr=200;
const bool isSimple=0;
