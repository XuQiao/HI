#include "TF1.h"

//TF1 *PtDistr  = new TF1("PtDistr","exp (-(x/0.90))+0.15*exp (-(x/15))", 0.2,10);	//V~=0.06
//TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.2,10);
//V2vsPt->SetParameters(4.7,1.8,3.0,1.8,0.8);
TF1 *V2vsPt = new TF1("V2vsPt","((x/4.7)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);    //Toy V~0.06

const int nFileAll = 51;
const TString dir = "/lio/lfs/cms/store/user/qixu/flow/STEG/ToyV206m300/";
const double ptmin = 0.2;
const double ptmax = 10.0; 
const double ptbin[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
//const double ptbin[]={0.2,5.0,10.0};
const int npt= 13;
const int trkbin[]={310,290};
const int nbin = 1;
//const int trkbin[]={185,180,175,170,165,160,155,150};
const int nn=2;
const int mm=1;
const double j01=2.404826;
const int ntheta = 5;
const int nstepr=200;
const bool isSimple = 1;
