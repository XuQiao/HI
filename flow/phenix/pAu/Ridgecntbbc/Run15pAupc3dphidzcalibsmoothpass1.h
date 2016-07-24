#include "TF1.h"
double calcsdphi(double dphi, int arm, int ch, double mom);
double calcsdphi(double dphi, int arm, int ch, double mom);


double calcsdphi(double dphi, int arm, int ch, double mom){
    TF1 *fm = new TF1("fm","[0]+[1]*x+[2]/x+[3]/sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0,10);
    TF1 *fs = new TF1("fs","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]/TMath::Sqrt(x)+[7]/x/x",0,10);
    int iarm, ich;
        if(arm==0) iarm=0;
        else if(arm==1) iarm=1;
        else return -9999;

        if(ch>0) ich = 0;
        else if(ch<0) ich = 1;
        else return -9999;
if(iarm==0 && ich==0)
fm->SetParameters(0.000390595,0.000162176,0.000251133,0.000396356,-2.60345e-05,-7.99175e-05,-6.25229e-05);
else if(iarm==0 && ich ==1)
fm->SetParameters(0.000223158,7.94167e-05,0.000227081,0.000265483,8.77645e-05,2.71799e-05,8.84937e-06);
else if(iarm==1 && ich ==0)
fm->SetParameters(-0.00011239,2.14409e-05,-0.000291798,-0.000241901,-0.00014483,-2.42003e-05,1.42447e-05);
else if(iarm==1 && ich ==1)
fm->SetParameters(-0.00042893,-0.000197551,-0.000157779,-0.000371043,9.79158e-05,7.84645e-05,3.25075e-05);

double dphimean = fm->Eval(mom);

if(iarm==0 && ich==0)
fs->SetParameters(0.000510613,0.000100915,1.44829e-05,1.74355e-06,1.33778e-07,-1.45307e-08,0.000831428,0.000731629);
else if(iarm==0 && ich ==1)
fs->SetParameters(0.000521444,0.000106384,1.45324e-05,1.44706e-06,1.39689e-08,-4.98669e-08,0.000813302,0.000595076);
else if(iarm==1 && ich ==0)
fs->SetParameters(0.000625007,0.000126489,1.72462e-05,1.74888e-06,3.53801e-08,-5.19499e-08,0.000979725,0.000727521);
else if(iarm==1 && ich ==1)
fs->SetParameters(0.000642572,0.000132943,1.88101e-05,2.10006e-06,1.13522e-07,-3.4543e-08,0.000997731,0.00072021);
double dphisigma = fs->Eval(mom);
delete fm;
delete fs;
return (dphi-dphimean)/dphisigma;
 }

double calcsdz(double dz, int arm, int ch, double mom){
    TF1 *fm = new TF1("fm","[0]+[1]*x+[2]/x+[3]/sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0,10);
    TF1 *fs = new TF1("fs","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]/TMath::Sqrt(x)+[7]/x/x",0,10);
    int iarm, ich;
        if(arm==0) iarm=0;
        else if(arm==1) iarm=1;
        else return -9999;

        if(ch>0) ich = 0;
        else if(ch<0) ich = 1;
        else return -9999;
if(iarm==0 && ich==0)
fm->SetParameters(-1.36626,0.123143,0.436233,1.88276,-2.09071,1.32394,-0.283277);
else if(iarm==0 && ich ==1)
fm->SetParameters(0.140301,-0.0127923,-0.107547,0.0139815,-0.000816549,0.0373167,-0.0129189);
else if(iarm==1 && ich ==0)
fm->SetParameters(0.550266,0.0857204,0.0695367,1.9055,-1.60001,1.00307,-0.211375);
else if(iarm==1 && ich ==1)
fm->SetParameters(0.441096,0.113657,0.24521,1.92236,-1.55237,0.940731,-0.200312);
double dzmean = fm->Eval(mom);
if(iarm==0 && ich==0)
fs->SetParameters(0.933913,0.128747,-0.030867,0.00308107,0.0035487,-0.000753025,0.579085,0.103329);
else if(iarm==0 && ich ==1)
fs->SetParameters(0.316736,0.192084,0.0250001,-0.00634049,-0.0012941,0.000249124,1.13892,0.0390641);
else if(iarm==1 && ich ==0)
fs->SetParameters(0.327402,0.141654,0.0203898,-0.00486675,-0.000895856,0.000203072,1.36016,0.0153886);
else if(iarm==1 && ich ==1)
fs->SetParameters(0.337825,0.161366,0.0214142,-0.00533378,-0.00100516,0.000196641,1.27023,0.0226192);
double dzsigma = fs->Eval(mom);
delete fm;
delete fs;
return (dz-dzmean)/dzsigma;
}
