#ifndef constants__h__
#define constants__h__

#include <TROOT.h>
#include <TMath.h>
#include <cmath>
#include <vector>
using namespace std;

const double Pi_ = 3.14159265358979312;
const double denom = 2.50662827463100024;//sqrt(2.0*Pi_);

double ksMassXmin = 0.340;
double ksMassXmax = 0.740;
int ksMassNbins = 400;
double ksMassBinWidth = (ksMassXmax - ksMassXmin) / ksMassNbins;
//double ksMassBinWidth_DG = ksMassBinWidth / sqrt(2.0*Pi_);//TMath::Pi());
double ksMassBinWidth_DG = ksMassBinWidth / denom;
//extern double ksMassBinWidth_DG;

double ksMassXmin_narrow = 0.440;
double ksMassXmax_narrow = 0.640;
int ksMassNbins_narrow = 200;

double lamMassXmin = 1.088;
double lamMassXmax = 1.188;
int lamMassNbins = 200;
double lamMassBinWidth = (lamMassXmax - lamMassXmin) / lamMassNbins;
//double lamMassBinWidth_DG = lamMassBinWidth / sqrt(2.0*Pi_);//TMath::Pi());
double lamMassBinWidth_DG = lamMassBinWidth / denom;
//extern lamMassBinWidth_DG;

double lamMassXmin_narrow = 1.09;
double lamMassXmax_narrow = 1.15;
int lamMassNbins_narrow = 120;

double ksyXmin = 0.;//-2.5;
double ksyXmax = 2.3;
int ksyNbins = 23;
double ksyBinWidth = (ksyXmax - ksyXmin) / ksyNbins;

double ksyXmin_signed = -2.3;//-2.5;
double ksyXmax_signed = 2.3;
int ksyNbins_signed = 46;
double ksyBinWidth_signed = (ksyXmax_signed - ksyXmin_signed) / ksyNbins_signed;

double lamyXmin = 0.;
double lamyXmax = 2.3;
int lamyNbins = 23;
double lamyBinWidth = (lamyXmax - lamyXmin) / lamyNbins;

double lamyXmin_signed = -2.3;
double lamyXmax_signed = 2.3;
int lamyNbins_signed = 46;
double lamyBinWidth_signed = (lamyXmax_signed - lamyXmin_signed) / lamyNbins;
//44
double ksptarray[44] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0,7.0,8.0,10.0};
int ksptNbins = sizeof(ksptarray)/sizeof(double) -1;
//cout << "pT " << ksptNbins << endl;
extern vector<double> ksptBinMinima;
//std::cout << ksptBinMinima.size() << endl;
/*for ( int ndx = 0; ndx <= ksptNbins; ndx++ ) {
  ksptBinMinima.push_back( ksptarray[ndx] );
  }*/
double kspt_default_binWidth = 0.1;
//41
double lamptarray[41] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0,7.0,8.0,10.0};
int lamptNbins = sizeof(lamptarray)/sizeof(double) -1;
//cout << "lpT " << lamptNbins << endl;
extern vector<double> lamptBinMinima;
/*for ( int ndx = 0; ndx <= lamptNbins; ndx++ ) {
  lamptBinMinima.push_back( lamptarray[ndx] );
  }*/
double lampt_default_binWidth = 0.1;
//44
double ksparray[44] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.4, 5.8, 6.2, 6.6, 7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 16.0, 18.0, 20.0};
int kspNbins = sizeof(ksparray)/sizeof(double) -1;
//cout << "p " << kspNbins << endl;
extern vector<double> kspBinMinima;
/*for ( int ndx = 0; ndx <= kspNbins; ndx++ ) {
  kspBinMinima.push_back( ksparray[ndx] );
  }*/
double ksp_default_binWidth = 0.2;
//41
double lamparray[42] = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.4,5.8,6.2,6.6,7.0,7.4,7.8,8.2,8.6,9.0,10.0,11.0,12.0,13.0,14.0,16.0,18.0,20.0};
int lampNbins = sizeof(lamparray)/sizeof(double) -1;
//cout << "lp " << lampNbins << endl;
extern vector<double> lampBinMinima;
/*for ( int ndx = 0; ndx <= lampNbins; ndx++ ) {
  lampBinMinima.push_back( lamparray[ndx] );
  }*/
double lamp_default_binWidth = 0.2;
//61
double kslifearray[61] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.4,6.8,7.2,7.6,8.0,8.8,9.6,10.4,11.2,12.0};
int kslifeNbins = sizeof(kslifearray)/sizeof(double) -1;
//cout << "life " << kslifeNbins << endl;
extern vector<double> kslifeBinMinima;
/*for ( int ndx = 0; ndx <= kslifeNbins; ndx++ ) {
  kslifeBinMinima.push_back( kslifearray[ndx] );
  }*/
double kslife_default_binWidth = 0.1;
//41
double lamlifearray[42] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0, 8.8, 9.6, 10.4, 11.2, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};
int lamlifeNbins = sizeof(lamlifearray)/sizeof(double) -1;
//cout << "llife " << lamlifeNbins << endl;
extern vector<double> lamlifeBinMinima;
/*for ( int ndx = 0; ndx <= lamlifeNbins; ndx++ ) {
  lamlifeBinMinima.push_back( lamlifearray[ndx] );
  }*/
double lamlife_default_binWidth = 0.1;

double ntrkarray[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,84,88,92,96,100,110,120,140,160,200};
int ntrkNbins = sizeof(ntrkarray)/sizeof(double) -1;
//cout << "ntrk " << ntrkNbins << endl;
extern vector<double> ntrkBinMinima;
/*for ( int ndx = 0; ndx <= ntrkNbins; ndx++ ) {
  ntrkBinMinima.push_back( ntrkarray[ndx] );
  }*/
double ntrk_default_binWidth = 0.1;
#endif

