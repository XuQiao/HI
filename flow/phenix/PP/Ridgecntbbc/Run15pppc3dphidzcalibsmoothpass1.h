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
fm->SetParameters(0.000430828,0.000184655,0.000255305,0.000424993,-4.54449e-05,-9.50015e-05,-7.1242e-05);
else if(iarm==0 && ich ==1)
fm->SetParameters(0.000225081,8.27579e-05,0.00022704,0.000264965,9.43863e-05,3.53653e-05,1.50746e-05);
else if(iarm==1 && ich ==0)
fm->SetParameters(-0.000141353,3.34851e-05,-0.000349821,-0.000300917,-0.000153129,-1.51126e-05,2.18853e-05);
else if(iarm==1 && ich ==1)
fm->SetParameters(-0.00045831,-0.000164498,-0.000225325,-0.000443789,9.92433e-05,9.66958e-05,4.55827e-05);

double dphimean = fm->Eval(mom);

if(iarm==0 && ich==0)
fs->SetParameters(0.000521002,0.000109566,1.78554e-05,2.77799e-06,4.18492e-07,5.93107e-08,0.000835199,0.000723165);
else if(iarm==0 && ich ==1)
fs->SetParameters(0.000505937,0.000103122,1.57483e-05,2.18814e-06,2.66865e-07,2.21593e-08,0.000815376,0.000701127);
else if(iarm==1 && ich ==0)
fs->SetParameters(0.000606298,0.00012682,1.9435e-05,2.67133e-06,3.1842e-07,2.48182e-08,0.000955767,0.000757969);
else if(iarm==1 && ich ==1)
fs->SetParameters(0.000661511,0.000146384,2.15976e-05,2.55679e-06,1.87326e-07,-2.12939e-08,0.000975135,0.00055282);
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
fm->SetParameters(0.0235701,0.00629276,-0.0161878,0.0121738,-0.0310031,-0.00672129,0.00732975);
else if(iarm==0 && ich ==1)
fm->SetParameters(0.00506442,0.0231584,-0.00957774,-0.0100336,0.00666027,0.00305449,-0.00359048);
else if(iarm==1 && ich ==0)
fm->SetParameters(1.09544,0.0473387,-0.119676,0.750498,-0.263252,0.114951,-0.0154389);
else if(iarm==1 && ich ==1)
fm->SetParameters(1.11282,0.0376028,-0.101171,0.779789,-0.272143,0.111515,-0.0110738);

double dzmean = fm->Eval(mom);
if(iarm==0 && ich==0)
fs->SetParameters(0.801104,0.200279,-0.00278408,-0.00722532,-0.000753957,0.000513378,0.690238,0.10851);
else if(iarm==0 && ich ==1)
fs->SetParameters(0.627932,0.117278,0.0169538,0.00198363,-0.000214271,-0.000294563,0.96669,0.0450957);
else if(iarm==1 && ich ==0)
fs->SetParameters(1.60725,0.0941545,-0.0942928,0.00776694,0.00816979,-0.00147169,0.118279,0.173903);
else if(iarm==1 && ich ==1)
fs->SetParameters(1.12485,0.192289,-0.0683635,0.000692004,0.00674185,-0.00112982,0.467116,0.134625);
double dzsigma = fs->Eval(mom);
delete fm;
delete fs;
return (dz-dzmean)/dzsigma;
}
