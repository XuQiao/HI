#include "fstream.h"
#include "SimplifyLife.C"

void makesysratio(string t_="ptccentc"){
TString type = t_.c_str();
if(type.Contains("pt") && type != "ptccentc"){
const int ncent = 1;
double centbin[ncent+1] = {0,5}; 
double centmidp[ncent] = {50.3333}; //mb
if(type=="ptIn"){//1.0-3.0
const int npt = 1;
double ptbin[npt+1] = {1.0,3.0};
double ptmean[npt] = {1.36878};
}
else if(type=="ptIn25_4"){
const int npt = 1;
double ptbin[npt+1] = {2.5,4.0};
double ptmean[npt] = {2.92489};
}
else if(type=="ptcoarser"){
const int npt = 4;
double ptbin[npt+1] = {0.2,1.0,2.0,3.0,5.0};
const double ptmean[npt] = {0.519639, 1.29345, 2.32523, 3.51803};
}
else if(type=="ptfiner"){
const int npt = 10;
const double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
const double ptmean[npt] = {0.360943, 0.691833, 1.1911, 1.69654, 2.20117, 2.70571, 3.2097, 3.71372, 4.21814, 4.72014};
}
}
else if(type=="centIn"){
const int ncent = 1;
const double centbin[ncent+1] = {0,100};
const double centmidp[ncent] = {50};
const int npt = 10;
const double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
const double ptmean[npt] = {0.360943, 0.691833, 1.1911, 1.69654, 2.20117, 2.70571, 3.2097, 3.71372, 4.21814, 4.72014};
}
else if(type=="ptccentc"){
const int ncent = 2;
const double centbin[ncent+1] = {0,5,20};
const double centmidp[ncent] = {2.5,12.5};
const int npt = 5;
const int npt = 10;
const double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
const double ptmean[npt] = {0.360943, 0.691833, 1.1911, 1.69654, 2.20117, 2.70571, 3.2097, 3.71372, 4.21814, 4.72014};
}
double c1north[ncent][npt];
double c1south[ncent][npt];
double c1northerr[ncent][npt];
double c1southerr[ncent][npt];
double c2north[ncent][npt];
double c2south[ncent][npt];
double c2northerr[ncent][npt];
double c2southerr[ncent][npt];
double c3north[ncent][npt];
double c3south[ncent][npt];
double c3northerr[ncent][npt];
double c3southerr[ncent][npt];
double c1north0[ncent][npt];
double c1south0[ncent][npt];
double c1northerr0[ncent][npt];
double c1southerr0[ncent][npt];
double c2north0[ncent][npt];
double c2south0[ncent][npt];
double c2northerr0[ncent][npt];
double c2southerr0[ncent][npt];
double c3north0[ncent][npt];
double c3south0[ncent][npt];
double c3northerr0[ncent][npt];
double c3southerr0[ncent][npt];
double c1north1[ncent][npt];
double c1south1[ncent][npt];
double c1northerr1[ncent][npt];
double c1southerr1[ncent][npt];
double c2north1[ncent][npt];
double c2south1[ncent][npt];
double c2northerr1[ncent][npt];
double c2southerr1[ncent][npt];
double c3north1[ncent][npt];
double c3south1[ncent][npt];
double c3northerr1[ncent][npt];
double c3southerr1[ncent][npt];
double c1northratio0[ncent][npt];
double c1southratio0[ncent][npt];
double c1northerrratio0[ncent][npt];
double c1southerrratio0[ncent][npt];
double c2northratio0[ncent][npt];
double c2southratio0[ncent][npt];
double c2northerrratio0[ncent][npt];
double c2southerrratio0[ncent][npt];
double c3northratio0[ncent][npt];
double c3southratio0[ncent][npt];
double c3northerrratio0[ncent][npt];
double c3southerrratio0[ncent][npt];
double c1northratio1[ncent][npt];
double c1southratio1[ncent][npt];
double c1northerrratio1[ncent][npt];
double c1southerrratio1[ncent][npt];
double c2northratio1[ncent][npt];
double c2southratio1[ncent][npt];
double c2northerrratio1[ncent][npt];
double c2southerrratio1[ncent][npt];
double c3northratio1[ncent][npt];
double c3southratio1[ncent][npt];
double c3northerrratio1[ncent][npt];
double c3southerrratio1[ncent][npt];
double c1north2[ncent][npt];
double c1south2[ncent][npt];
double c1northerr2[ncent][npt];
double c1southerr2[ncent][npt];
double c2north2[ncent][npt];
double c2south2[ncent][npt];
double c2northerr2[ncent][npt];
double c2southerr2[ncent][npt];
double c3north2[ncent][npt];
double c3south2[ncent][npt];
double c3northerr2[ncent][npt];
double c3southerr2[ncent][npt];
double c1northratio2[ncent][npt];
double c1southratio2[ncent][npt];
double c1northerrratio2[ncent][npt];
double c1southerrratio2[ncent][npt];
double c2northratio2[ncent][npt];
double c2southratio2[ncent][npt];
double c2northerrratio2[ncent][npt];
double c2southerrratio2[ncent][npt];
double c3northratio2[ncent][npt];
double c3southratio2[ncent][npt];
double c3northerrratio2[ncent][npt];
double c3southerrratio2[ncent][npt];
ifstream fnorth;
ifstream fsouth;
ifstream fnorth0;
ifstream fsouth0;
ifstream fnorth1;
ifstream fsouth1;
ifstream fnorth2;
ifstream fsouth2;
fnorth.open(Form("c1_c2_central_%s_north.dat",type.Data()));
fsouth.open(Form("c1_c2_West_PC3Sigma2_central_%s_south.dat",type.Data()));
fnorth0.open(Form("c1_c2_central_%s_north.dat",type.Data()));
fsouth0.open(Form("c1_c2_West_PC3Sigma2_central_%s_south.dat",type.Data()));
//fnorth1.open(Form("c1_c2_PUlow_central_%s_north.dat",type.Data()));
//fsouth1.open(Form("c1_c2_PUlow_PC3Sigma2_central_%s_south.dat",type.Data()));
//fnorth2.open(Form("c1_c2_PUhigh_central_%s_north.dat",type.Data()));
//fsouth2.open(Form("c1_c2_PUhigh_PC3Sigma2_central_%s_south.dat",type.Data()));
//fnorth1.open(Form("c1_c2_PUlow_central_%s_north.dat",type.Data()));
//fsouth1.open(Form("c1_c2_west_PC3Sigma2_central_%s_south.dat",type.Data()));
//fnorth2.open(Form("c1_c2_PUhigh_central_%s_north.dat",type.Data()));
//fsouth2.open(Form("c1_c2_West_PC3Sigma2_central_%s_south.dat",type.Data()));
fnorth1.open(Form("c1_c2_bbcr2_central_%s_north.dat","ptfiner"));
fsouth1.open(Form("c1_c2_west_bbcr2_central_%s_south.dat","ptfiner"));
ofstream fnorthout0("/dev/null");
ofstream fsouthout0("/dev/null");
ofstream fnorthout1("/dev/null");
ofstream fsouthout1;
ofstream fnorthout2("/dev/null");
ofstream fsouthout2;
//fnorthout0.open(Form("c1_c2_PC3Sigma2ratio_central_%s_north.dat",type.Data()));
//fsouthout0.open(Form("c1_c2_PC3Sigma2ratio_central_%s_south.dat",type.Data()));
//fnorthout1.open(Form("c1_c2_PUlowratio_central_%s_north.dat",type.Data()));
//fsouthout1.open(Form("c1_c2_PUlowratio_PC3Sigma2_central_%s_south.dat",type.Data()));
//fsouthout1.open(Form("c1_c2_Westratio_PC3Sigma2_central_%s_south.dat",type.Data()));
fsouthout1.open(Form("c1_c2_west_bbcr2ratio_PC3Sigma2_central_%s_south.dat",type.Data()));
//fnorthout2.open(Form("c1_c2_PUhighratio_central_%s_north.dat",type.Data()));
//fsouthout2.open(Form("c1_c2_PUhighratio_PC3Sigma2_central_%s_south.dat",type.Data()));
//fsouthout2.open(Form("c1_c2_westratio_PC3Sigma2_central_%s_south.dat",type.Data()));
//fsouthout2.open(Form("c1_c2_bbcr2ratio_PC3Sigma2_central_%s_south.dat",type.Data()));
for(int icent=0;icent<ncent;icent++){
for(int ipt=0;ipt<npt;ipt++){
fnorth>>c1north[icent][ipt];
fnorth>>c1northerr[icent][ipt];
fnorth>>c2north[icent][ipt];
fnorth>>c2northerr[icent][ipt];
fnorth>>c3north[icent][ipt];
fnorth>>c3northerr[icent][ipt];
fsouth>>c1south[icent][ipt];
fsouth>>c1southerr[icent][ipt];
fsouth>>c2south[icent][ipt];
fsouth>>c2southerr[icent][ipt];
fsouth>>c3south[icent][ipt];
fsouth>>c3southerr[icent][ipt];
fnorth0>>c1north0[icent][ipt];
fnorth0>>c1northerr0[icent][ipt];
fnorth0>>c2north0[icent][ipt];
fnorth0>>c2northerr0[icent][ipt];
fnorth0>>c3north0[icent][ipt];
fnorth0>>c3northerr0[icent][ipt];
fsouth0>>c1south0[icent][ipt];
fsouth0>>c1southerr0[icent][ipt];
fsouth0>>c2south0[icent][ipt];
fsouth0>>c2southerr0[icent][ipt];
fsouth0>>c3south0[icent][ipt];
fsouth0>>c3southerr0[icent][ipt];
fnorth1>>c1north1[icent][ipt];
fnorth1>>c1northerr1[icent][ipt];
fnorth1>>c2north1[icent][ipt];
fnorth1>>c2northerr1[icent][ipt];
fnorth1>>c3north1[icent][ipt];
fnorth1>>c3northerr1[icent][ipt];
fsouth1>>c1south1[icent][ipt];
fsouth1>>c1southerr1[icent][ipt];
fsouth1>>c2south1[icent][ipt];
fsouth1>>c2southerr1[icent][ipt];
fsouth1>>c3south1[icent][ipt];
fsouth1>>c3southerr1[icent][ipt];
fnorth2>>c1north2[icent][ipt];
fnorth2>>c1northerr2[icent][ipt];
fnorth2>>c2north2[icent][ipt];
fnorth2>>c2northerr2[icent][ipt];
fnorth2>>c3north2[icent][ipt];
fnorth2>>c3northerr2[icent][ipt];
fsouth2>>c1south2[icent][ipt];
fsouth2>>c1southerr2[icent][ipt];
fsouth2>>c2south2[icent][ipt];
fsouth2>>c2southerr2[icent][ipt];
fsouth2>>c3south2[icent][ipt];
fsouth2>>c3southerr2[icent][ipt];

c1northratio0[icent][ipt] = c1north0[icent][ipt]/c1north[icent][ipt];
c1northerrratio0[icent][ipt]=fsys(c1north0[icent][ipt],c1northerr0[icent][ipt],c1north[icent][ipt],c1northerr[icent][ipt]);
c2northratio0[icent][ipt] = c2north0[icent][ipt]/c2north[icent][ipt];
c2northerrratio0[icent][ipt]=fsys(c2north0[icent][ipt],c2northerr0[icent][ipt],c2north[icent][ipt],c2northerr[icent][ipt]);
c3northratio0[icent][ipt] = c3north0[icent][ipt]/c3north[icent][ipt];
c3northerrratio0[icent][ipt]=fsys(c3north0[icent][ipt],c3northerr0[icent][ipt],c3north[icent][ipt],c3northerr[icent][ipt]);
c1southratio0[icent][ipt] = c1south0[icent][ipt]/c1south[icent][ipt];
c1southerrratio0[icent][ipt]=fsys(c1south0[icent][ipt],c1southerr0[icent][ipt],c1south[icent][ipt],c1southerr[icent][ipt]);
c2southratio0[icent][ipt] = c2south0[icent][ipt]/c2south[icent][ipt];
c2southerrratio0[icent][ipt]=fsys(c2south0[icent][ipt],c2southerr0[icent][ipt],c2south[icent][ipt],c2southerr[icent][ipt]);
c3southratio0[icent][ipt] = c3south0[icent][ipt]/c3south[icent][ipt];
c3southerrratio0[icent][ipt]=fsys(c3south0[icent][ipt],c3southerr0[icent][ipt],c3south[icent][ipt],c3southerr[icent][ipt]);

c1northratio1[icent][ipt] = c1north1[icent][ipt]/c1north[icent][ipt];
c1northerrratio1[icent][ipt]=fsys(c1north1[icent][ipt],c1northerr1[icent][ipt],c1north[icent][ipt],c1northerr[icent][ipt]);
c2northratio1[icent][ipt] = c2north1[icent][ipt]/c2north[icent][ipt];
c2northerrratio1[icent][ipt]=fsys(c2north1[icent][ipt],c2northerr1[icent][ipt],c2north[icent][ipt],c2northerr[icent][ipt]);
c3northratio1[icent][ipt] = c3north1[icent][ipt]/c3north[icent][ipt];
c3northerrratio1[icent][ipt]=fsys(c3north1[icent][ipt],c3northerr1[icent][ipt],c3north[icent][ipt],c3northerr[icent][ipt]);
c1southratio1[icent][ipt] = c1south1[icent][ipt]/c1south[icent][ipt];
c1southerrratio1[icent][ipt]=fsys(c1south1[icent][ipt],c1southerr1[icent][ipt],c1south[icent][ipt],c1southerr[icent][ipt]);
c2southratio1[icent][ipt] = c2south1[icent][ipt]/c2south[icent][ipt];
c2southerrratio1[icent][ipt]=fsys(c2south1[icent][ipt],c2southerr1[icent][ipt],c2south[icent][ipt],c2southerr[icent][ipt]);
c3southratio1[icent][ipt] = c3south1[icent][ipt]/c3south[icent][ipt];
c3southerrratio1[icent][ipt]=fsys(c3south1[icent][ipt],c3southerr1[icent][ipt],c3south[icent][ipt],c3southerr[icent][ipt]);

c1northratio2[icent][ipt] = c1north2[icent][ipt]/c1north[icent][ipt];
c1northerrratio2[icent][ipt]=fsys(c1north2[icent][ipt],c1northerr1[icent][ipt],c1north[icent][ipt],c1northerr[icent][ipt]);
c2northratio2[icent][ipt] = c2north2[icent][ipt]/c2north[icent][ipt];
c2northerrratio2[icent][ipt]=fsys(c2north2[icent][ipt],c2northerr2[icent][ipt],c2north[icent][ipt],c2northerr[icent][ipt]);
c3northratio2[icent][ipt] = c3north2[icent][ipt]/c3north[icent][ipt];
c3northerrratio2[icent][ipt]=fsys(c3north2[icent][ipt],c3northerr2[icent][ipt],c3north[icent][ipt],c3northerr[icent][ipt]);
c1southratio2[icent][ipt] = c1south2[icent][ipt]/c1south[icent][ipt];
c1southerrratio2[icent][ipt]=fsys(c1south2[icent][ipt],c1southerr2[icent][ipt],c1south[icent][ipt],c1southerr[icent][ipt]);
c2southratio2[icent][ipt] = c2south2[icent][ipt]/c2south[icent][ipt];
c2southerrratio2[icent][ipt]=fsys(c2south2[icent][ipt],c2southerr2[icent][ipt],c2south[icent][ipt],c2southerr[icent][ipt]);
c3southratio2[icent][ipt] = c3south2[icent][ipt]/c3south[icent][ipt];
c3southerrratio2[icent][ipt]=fsys(c3south2[icent][ipt],c3southerr2[icent][ipt],c3south[icent][ipt],c3southerr[icent][ipt]);

fnorthout0<<c1northratio0[icent][ipt];
fnorthout0<<"\t";
fnorthout0<<c1northerrratio0[icent][ipt];
fnorthout0<<"\t";
fnorthout0<<c2northratio0[icent][ipt];
fnorthout0<<"\t";
fnorthout0<<c2northerrratio0[icent][ipt];
fnorthout0<<"\t";
fnorthout0<<c3northratio0[icent][ipt];
fnorthout0<<"\t";
fnorthout0<<c3northerrratio0[icent][ipt];
fnorthout0<<"\t"<<endl;
fsouthout0<<c1southratio0[icent][ipt];
fsouthout0<<"\t";
fsouthout0<<c1southerrratio0[icent][ipt];
fsouthout0<<"\t";
fsouthout0<<c2southratio0[icent][ipt];
fsouthout0<<"\t";
fsouthout0<<c2southerrratio0[icent][ipt];
fsouthout0<<"\t";
fsouthout0<<c3southratio0[icent][ipt];
fsouthout0<<"\t";
fsouthout0<<c3southerrratio0[icent][ipt];
fsouthout0<<"\t"<<endl;

fnorthout1<<c1northratio1[icent][ipt];
fnorthout1<<"\t";
fnorthout1<<c1northerrratio1[icent][ipt];
fnorthout1<<"\t";
fnorthout1<<c2northratio1[icent][ipt];
fnorthout1<<"\t";
fnorthout1<<c2northerrratio1[icent][ipt];
fnorthout1<<"\t";
fnorthout1<<c3northratio1[icent][ipt];
fnorthout1<<"\t";
fnorthout1<<c3northerrratio1[icent][ipt];
fnorthout1<<"\t"<<endl;
fsouthout1<<c1southratio1[icent][ipt];
fsouthout1<<"\t";
fsouthout1<<c1southerrratio1[icent][ipt];
fsouthout1<<"\t";
fsouthout1<<c2southratio1[icent][ipt];
fsouthout1<<"\t";
fsouthout1<<c2southerrratio1[icent][ipt];
fsouthout1<<"\t";
fsouthout1<<c3southratio1[icent][ipt];
fsouthout1<<"\t";
fsouthout1<<c3southerrratio1[icent][ipt];
fsouthout1<<"\t"<<endl;

fnorthout2<<c1northratio2[icent][ipt];
fnorthout2<<"\t";
fnorthout2<<c1northerrratio2[icent][ipt];
fnorthout2<<"\t";
fnorthout2<<c2northratio2[icent][ipt];
fnorthout2<<"\t";
fnorthout2<<c2northerrratio2[icent][ipt];
fnorthout2<<"\t";
fnorthout2<<c3northratio2[icent][ipt];
fnorthout2<<"\t";
fnorthout2<<c3northerrratio2[icent][ipt];
fnorthout2<<"\t"<<endl;
fsouthout2<<c1southratio2[icent][ipt];
fsouthout2<<"\t";
fsouthout2<<c1southerrratio2[icent][ipt];
fsouthout2<<"\t";
fsouthout2<<c2southratio2[icent][ipt];
fsouthout2<<"\t";
fsouthout2<<c2southerrratio2[icent][ipt];
fsouthout2<<"\t";
fsouthout2<<c3southratio2[icent][ipt];
fsouthout2<<"\t";
fsouthout2<<c3southerrratio2[icent][ipt];
fsouthout2<<"\t"<<endl;
}
}

}

double fsys(double a, double ea, double b, double eb){
    double c = a/b;
    return c*sqrt((ea/a)**2+(eb/b)**2);
}
