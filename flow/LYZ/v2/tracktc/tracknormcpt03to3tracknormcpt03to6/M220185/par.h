#include "../par.h"
const double ptmin = 0.3;
const double ptmax = 3.0; 
const int nFileAll = 544;
//const double ptbinV[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
//const double ptbinV[]={0.3,6.0};
//const double etabinV[]={-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4};
//const int nptV= 1;
//const int netaV= 24;
//const double ptbinv[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
//const double ptbinv[]={0.3,0.4,0.7,1.2,1.7,2.2,2.7,3.4,4.4,5.4,6.0};
const double ptbinv[]={0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
const double etabinv[]={-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4};
//const double etabinv[]={-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4};
const int nptv= 9;
const int netav= 24;//12
const int trkbin[]={220,185};
const int nbin = 1;
//const int trkbin[]={185,180,175,170,165,160,155,150};
const int nn=2;
const int mm=1;
const double j01=2.404826;
const int ntheta = 5;
const int nstepr=200;
const bool isSimple=0;
