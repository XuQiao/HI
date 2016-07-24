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
fm->SetParameters(0.00041587,0.000175554,0.00025773,0.000416249,-3.35605e-05,-8.62726e-05,-6.62246e-05);
else if(iarm==0 && ich ==1)
fm->SetParameters(0.00022411,8.16704e-05,0.000221541,0.000262646,8.49391e-05,2.7523e-05,9.96506e-06);
else if(iarm==1 && ich ==0)
fm->SetParameters(-0.000109552,4.0495e-05,-0.000330074,-0.000265666,-0.00016253,-2.38534e-05,1.89108e-05);
else if(iarm==1 && ich ==1)
fm->SetParameters(-0.000495161,-0.000211503,-0.000184715,-0.000436749,0.000121629,9.10395e-05,3.25172e-05);

double dphimean = fm->Eval(mom);

if(iarm==0 && ich==0)
fs->SetParameters(0.000513172,0.000104684,1.62529e-05,2.33205e-06,3.01278e-07,2.93183e-08,0.000831069,0.000730834);
else if(iarm==0 && ich ==1)
fs->SetParameters(0.000494769,9.73738e-05,1.39591e-05,1.70024e-06,1.39795e-07,-1.01996e-08,0.000806399,0.000705748);
else if(iarm==1 && ich ==0)
fs->SetParameters(0.000619904,0.000127482,1.93093e-05,2.61998e-06,3.02616e-07,1.99548e-08,0.00098767,0.00080852);
else if(iarm==1 && ich ==1)
fs->SetParameters(0.000745113,0.000157073,2.01111e-05,1.41393e-06,-2.27421e-07,-1.4613e-07,0.00110136,0.000571697);
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
fm->SetParameters(-0.221457,0.0123548,0.116753,-0.0682329,-0.0205868,-0.0789914,0.0306131);
else if(iarm==0 && ich ==1)
fm->SetParameters(0.0511468,-0.0359837,-0.200837,-0.106342,0.0427844,0.0745268,-0.0346661);
else if(iarm==1 && ich ==0)
fm->SetParameters(1.28166,0.00436434,-0.336366,0.70848,-0.245502,0.17595,-0.0394824);
else if(iarm==1 && ich ==1)
fm->SetParameters(1.1991,0.0135342,-0.157248,0.766315,-0.237425,0.113147,-0.0219663);

double dzmean = fm->Eval(mom);
if(iarm==0 && ich==0)
fs->SetParameters(-0.402972,0.218426,0.102248,-0.0158504,-0.00844338,0.00162026,1.76882,-0.0340538);
else if(iarm==0 && ich ==1)
fs->SetParameters(1.15546,0.0699654,-0.0392842,0.00658073,0.00333088,-0.00073815,0.294522,0.131076);
else if(iarm==1 && ich ==0)
fs->SetParameters(1.15707,0.0154272,-0.046935,0.00730642,0.00377216,-0.000682076,0.499964,0.109393);
else if(iarm==1 && ich ==1)
fs->SetParameters(1.10562,0.0126302,-0.0384047,0.00974606,0.00383479,-0.000904214,0.504238,0.0994903);
double dzsigma = fs->Eval(mom);
delete fm;
delete fs;
return (dz-dzmean)/dzsigma;
}
