  double calcsdphi(double dphi, int arm, int ch, double mom);
  double calcsdz(double dphi, int arm, int ch, double mom);
  double calcsdphi(double dphi, int arm, int ch, double mom){
     int iarm, ich;
     if(arm==0) iarm=0;
     else if(arm==1) iarm=1;
     else return -9999;
     if(ch>0) ich = 0;
         else if(ch<0) ich = 1;
         else return -9999;
     double dphimean=-9999;
     double dphisigma=-9999;
if(iarm==0 && ich==0){
dphimean = 0.000386696+0.000205701*mom+0.000233234/mom+0.000365682/sqrt(mom)+-1.64513e-05/mom/mom+-7.44536e-05/mom/mom/mom+-6.09302e-05/mom/mom/mom/mom;
}
if(iarm==0 && ich==1){
dphimean = 0.000254517+0.000114193*mom+0.000242441/mom+0.00028339/sqrt(mom)+0.000102609/mom/mom+3.42829e-05/mom/mom/mom+1.17273e-05/mom/mom/mom/mom;
}
if(iarm==1 && ich==0){
dphimean = -0.000182701+-1.97934e-05*mom+-0.000286152/mom+-0.000275695/sqrt(mom)+-0.000131197/mom/mom+-2.11361e-05/mom/mom/mom+1.28311e-05/mom/mom/mom/mom;
}
if(iarm==1 && ich==1){
dphimean = -0.000412149+-0.000207506*mom+-0.000215095/mom+-0.000378615/sqrt(mom)+5.74654e-05/mom/mom+8.79539e-05/mom/mom/mom+5.35159e-05/mom/mom/mom/mom;
}
if(iarm==0 && ich==0){
dphisigma = 0.000489527+0.000116515*mom+1.85468e-05*mom*mom+1.81262e-06*mom*mom*mom+-2.39962e-07*mom*mom*mom*mom+-2.36303e-07*mom*mom*mom*mom*mom+0.00074777/sqrt(mom)+0.00069108/mom/mom;
}
if(iarm==0 && ich==1){
dphisigma = 0.000479732+0.000113641*mom+1.75746e-05*mom*mom+1.48511e-06*mom*mom*mom+-3.42453e-07*mom*mom*mom*mom+-2.66497e-07*mom*mom*mom*mom*mom+0.000731169/sqrt(mom)+0.000666555/mom/mom;
}
if(iarm==1 && ich==0){
dphisigma = 0.000603644+0.000149171*mom+2.42196e-05*mom*mom+2.32824e-06*mom*mom*mom+-3.61916e-07*mom*mom*mom*mom+-3.34432e-07*mom*mom*mom*mom*mom+0.000900971/sqrt(mom)+0.000776876/mom/mom;
}
if(iarm==1 && ich==1){
dphisigma = 0.000601728+0.000155385*mom+2.59688e-05*mom*mom+2.62575e-06*mom*mom*mom+-3.3234e-07*mom*mom*mom*mom+-3.3801e-07*mom*mom*mom*mom*mom+0.000871266/sqrt(mom)+0.000671047/mom/mom;
}
return (dphi-dphimean)/dphisigma;
}
  double calcsdz(double dz, int arm, int ch, double mom){
     int iarm, ich;
     if(arm==0) iarm=0;
     else if(arm==1) iarm=1;
     else return -9999;
     if(ch>0) ich = 0;
         else if(ch<0) ich = 1;
         else return -9999;
     double dzmean=-9999;
     double dzsigma=-9999;
if(iarm==0 && ich==0){
dzmean = 1.72817+0.0163229*mom+-0.917444/mom+0.44419/sqrt(mom)+-0.0645308/mom/mom+0.503528/mom/mom/mom+-0.184579/mom/mom/mom/mom;
}
if(iarm==0 && ich==1){
dzmean = 1.34462+0.0818882*mom+-0.484148/mom+0.550664/sqrt(mom)+-0.166213/mom/mom+0.2978/mom/mom/mom+-0.0958029/mom/mom/mom/mom;
}
if(iarm==1 && ich==0){
dzmean = 2.06963+-0.0214897*mom+-0.439009/mom+1.09432/sqrt(mom)+-0.391103/mom/mom+0.296423/mom/mom/mom+-0.0642472/mom/mom/mom/mom;
}
if(iarm==1 && ich==1){
dzmean = 1.80124+0.0584395*mom+-0.188296/mom+1.07296/sqrt(mom)+-0.354665/mom/mom+0.181106/mom/mom/mom+-0.0319963/mom/mom/mom/mom;
}
if(iarm==0 && ich==0){
dzsigma = -1.10337+0.234011*mom+0.301749*mom*mom+-0.0438853*mom*mom*mom+-0.0396016*mom*mom*mom*mom+0.00916634*mom*mom*mom*mom*mom+2.42644/sqrt(mom)+-0.0813971/mom/mom;
}
if(iarm==0 && ich==1){
dzsigma = 1.4073+0.121023*mom+-0.110989*mom*mom+0.0141683*mom*mom*mom+0.015246*mom*mom*mom*mom+-0.0036355*mom*mom*mom*mom*mom+0.0867441/sqrt(mom)+0.155359/mom/mom;
}
if(iarm==1 && ich==0){
dzsigma = 1.40226+0.0422227*mom+-0.116503*mom*mom+0.0200999*mom*mom*mom+0.0166596*mom*mom*mom*mom+-0.00435766*mom*mom*mom*mom*mom+0.303847/sqrt(mom)+0.131408/mom/mom;
}
if(iarm==1 && ich==1){
dzsigma = 1.4749+0.127885*mom+-0.126459*mom*mom+0.00892987*mom*mom*mom+0.0151614*mom*mom*mom*mom+-0.00300786*mom*mom*mom*mom*mom+0.0944681/sqrt(mom)+0.163548/mom/mom;
}
return (dz-dzmean)/dzsigma;
}
