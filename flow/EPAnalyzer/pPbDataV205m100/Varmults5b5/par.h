#include "TF1.h"

//TF1 *PtDistr  = new TF1("PtDistr","exp (-(x/0.90))+0.15*exp (-(x/15))", 0.2,10);	//V~=0.06
//TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.2,10);
//V2vsPt->SetParameters(4.7,1.8,3.0,1.8,0.8);
//TF1 *V2vsPt = new TF1("V2vsPt","((x/4.7)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);    //Toy V~0.06
TF1 *V2vsPt = new TF1("V2vsPt","((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);

const int nFileAll = 20;
const TString dir = "/cms/store/user/qixu/flow/STEGwithnf/pPbDataV205m100/Varmults5b5";
const double ptmin = 0.3;
const double ptmax = 6.0; 
const double etamax = 4.0; 
//const double ptbin[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
const double ptbin[]={0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
const int npt= 9;
const int trkbin[]={150,50};
const int nbin = 1;
const int neta=5;
const double etap[neta]={0,1.0,1.5,2,3};
const double etan[neta]={0,-1.0,-1.5,-2,-3};
const double etac = 0.35;
//const int trkbin[]={185,180,175,170,165,160,155,150};
const int nn=2;
const int mm=1;
