  double calcsdphi(double dphi, int arm, int ch, double mom, int RunNumber);
  double calcsdz(double dphi, int arm, int ch, double mom, int RunNumber);
  double calcsdphi(double dphi, int arm, int ch, double mom, int RunNumber){
     int iarm, ich;
     if(arm==0) iarm=0;
     else if(arm==1) iarm=1;
     else return -9999;
     if(ch>0) ich = 0;
         else if(ch<0) ich = 1;
         else return -9999;
     double dphimean=-9999;
     double dphisigma=-9999;
 if(RunNumber>=454774 && RunNumber<=455639){ //200 dAu
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
}
 if(RunNumber>=455792 && RunNumber<=456283){ //62 dAu
if(iarm==0 && ich==0){
dphimean = 0.000387409+0.000206806*mom+0.000233481/mom+0.000365928/sqrt(mom)+-1.52045e-05/mom/mom+-7.2903e-05/mom/mom/mom+-5.97065e-05/mom/mom/mom/mom;
}
if(iarm==0 && ich==1){
dphimean = 0.000249061+0.000112395*mom+0.000235651/mom+0.000276396/sqrt(mom)+9.91323e-05/mom/mom+3.30341e-05/mom/mom/mom+1.13401e-05/mom/mom/mom/mom;
}
if(iarm==1 && ich==0){
dphimean = -0.000185517+1.88245e-05*mom+-0.000356329/mom+-0.000323626/sqrt(mom)+-0.000161804/mom/mom+-1.40352e-05/mom/mom/mom+2.87395e-05/mom/mom/mom/mom;
}
if(iarm==1 && ich==1){
dphimean = -0.000429468+-0.000224325*mom+-0.000209137/mom+-0.000384583/sqrt(mom)+6.36625e-05/mom/mom+8.39601e-05/mom/mom/mom+4.55467e-05/mom/mom/mom/mom;
}
if(iarm==0 && ich==0){
dphisigma = 0.000486074+0.00011793*mom+1.94442e-05*mom*mom+2.12472e-06*mom*mom*mom+-1.50274e-07*mom*mom*mom*mom+-2.13014e-07*mom*mom*mom*mom*mom+0.000737748/sqrt(mom)+0.000674239/mom/mom;
}
if(iarm==0 && ich==1){
dphisigma = 0.000473043+0.000117384*mom+1.98721e-05*mom*mom+2.35315e-06*mom*mom*mom+-5.49284e-08*mom*mom*mom*mom+-1.76102e-07*mom*mom*mom*mom*mom+0.000709081/sqrt(mom)+0.000624218/mom/mom;
}
if(iarm==1 && ich==0){
dphisigma = 0.000607082+0.000155927*mom+2.73802e-05*mom*mom+3.49415e-06*mom*mom*mom+3.68985e-08*mom*mom*mom*mom+-2.02983e-07*mom*mom*mom*mom*mom+0.00089355/sqrt(mom)+0.000746462/mom/mom;
}
if(iarm==1 && ich==1){
dphisigma = 0.000615041+0.000163793*mom+2.8824e-05*mom*mom+3.55073e-06*mom*mom*mom+-2.01926e-08*mom*mom*mom*mom+-2.31095e-07*mom*mom*mom*mom*mom+0.000876461/sqrt(mom)+0.000636966/mom/mom;
}
}
 if(RunNumber>=457634 && RunNumber<=458167){ //39 dAu
if(iarm==0 && ich==0){
dphimean = 0.000406997+0.000221303*mom+0.000240543/mom+0.000380802/sqrt(mom)+-1.63473e-05/mom/mom+-7.41776e-05/mom/mom/mom+-6.02209e-05/mom/mom/mom/mom;
}
if(iarm==0 && ich==1){
dphimean = 0.000210701+8.96879e-05*mom+0.00021994/mom+0.00024439/sqrt(mom)+0.000106031/mom/mom+4.17751e-05/mom/mom/mom+1.72042e-05/mom/mom/mom/mom;
}
if(iarm==1 && ich==0){
dphimean = -0.00019166+7.1178e-06*mom+-0.00032683/mom+-0.000311094/sqrt(mom)+-0.000137977/mom/mom+-9.97279e-06/mom/mom/mom+2.422e-05/mom/mom/mom/mom;
}
if(iarm==1 && ich==1){
dphimean = -0.000446483+-0.000229112*mom+-0.000235813/mom+-0.000408847/sqrt(mom)+4.78238e-05/mom/mom+8.07889e-05/mom/mom/mom+4.84334e-05/mom/mom/mom/mom;
}
if(iarm==0 && ich==0){
dphisigma = 0.000495064+0.000123498*mom+2.18714e-05*mom*mom+3.08965e-06*mom*mom*mom+2.10078e-07*mom*mom*mom*mom+-8.53014e-08*mom*mom*mom*mom*mom+0.000745083/sqrt(mom)+0.000667545/mom/mom;
}
if(iarm==0 && ich==1){
dphisigma = 0.000486628+0.000127741*mom+2.41616e-05*mom*mom+3.96244e-06*mom*mom*mom+5.29461e-07*mom*mom*mom*mom+2.96802e-08*mom*mom*mom*mom*mom+0.000713247/sqrt(mom)+0.000589567/mom/mom;
}
if(iarm==1 && ich==0){
dphisigma = 0.000601272+0.000154738*mom+2.77539e-05*mom*mom+4.03304e-06*mom*mom*mom+3.62383e-07*mom*mom*mom*mom+-5.31903e-08*mom*mom*mom*mom*mom+0.000882439/sqrt(mom)+0.000717924/mom/mom;
}
if(iarm==1 && ich==1){
dphisigma = 0.000586283+0.000148054*mom+2.52194e-05*mom*mom+3.16782e-06*mom*mom*mom+8.58217e-08*mom*mom*mom*mom+-1.38318e-07*mom*mom*mom*mom*mom+0.000862509/sqrt(mom)+0.000695773/mom/mom;
}
}
 if(RunNumber>=456652 && RunNumber<=457298){ //20 dAu
if(iarm==0 && ich==0){
dphimean = 0.000494108+0.000282654*mom+0.000257806/mom+0.000442862/sqrt(mom)+-4.98834e-05/mom/mom+-0.000105136/mom/mom/mom+-7.9712e-05/mom/mom/mom/mom;
}
if(iarm==0 && ich==1){
dphimean = 0.000174598+7.36485e-05*mom+0.000204678/mom+0.000210871/sqrt(mom)+0.000126382/mom/mom+6.68864e-05/mom/mom/mom+3.53613e-05/mom/mom/mom/mom;
}
if(iarm==1 && ich==0){
dphimean = -8.78753e-05+0.000264206*mom+-0.000542035/mom+-0.000409182/sqrt(mom)+-0.000235003/mom/mom+2.27997e-05/mom/mom/mom+8.39267e-05/mom/mom/mom/mom;
}
if(iarm==1 && ich==1){
dphimean = -0.000622742+-0.000120075*mom+-0.000444449/mom+-0.000684558/sqrt(mom)+8.28602e-05/mom/mom+0.000146324/mom/mom/mom+7.54122e-05/mom/mom/mom/mom;
}
if(iarm==0 && ich==0){
dphisigma = 0.000564799+0.000147437*mom+2.84773e-05*mom*mom+4.68664e-06*mom*mom*mom+5.29861e-07*mom*mom*mom*mom+-4.51868e-08*mom*mom*mom*mom*mom+0.000842234/sqrt(mom)+0.000759544/mom/mom;
}
if(iarm==0 && ich==1){
dphisigma = 0.000533504+0.000148673*mom+3.1678e-05*mom*mom+6.46836e-06*mom*mom*mom+1.32531e-06*mom*mom*mom*mom+2.72348e-07*mom*mom*mom*mom*mom+0.000770083/sqrt(mom)+0.000627898/mom/mom;
}
if(iarm==1 && ich==0){
dphisigma = 0.000710837+0.000178106*mom+2.97078e-05*mom*mom+3.25434e-06*mom*mom*mom+-2.18964e-07*mom*mom*mom*mom+-3.16682e-07*mom*mom*mom*mom*mom+0.00105287/sqrt(mom)+0.000881485/mom/mom;
}
if(iarm==1 && ich==1){
dphisigma = 0.000936985+0.00027252*mom+3.74523e-05*mom*mom+-3.38057e-06*mom*mom*mom+-4.93978e-06*mom*mom*mom*mom+-2.49566e-06*mom*mom*mom*mom*mom+0.00118588/sqrt(mom)+0.00039032/mom/mom;
}
}
return (dphi-dphimean)/dphisigma;
}
  double calcsdz(double dz, int arm, int ch, double mom, int RunNumber){
     int iarm, ich;
     if(arm==0) iarm=0;
     else if(arm==1) iarm=1;
     else return -9999;
     if(ch>0) ich = 0;
         else if(ch<0) ich = 1;
         else return -9999;
     double dzmean=-9999;
     double dzsigma=-9999;
 if(RunNumber>=454774 && RunNumber<=455639){ //200 dAu
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
}
 if(RunNumber>=455792 && RunNumber<=456283){ //62 dAu
if(iarm==0 && ich==0){
dzmean = 0.976533+0.0657801*mom+-0.411341/mom+0.336518/sqrt(mom)+-0.0926394/mom/mom+0.225118/mom/mom/mom+-0.079196/mom/mom/mom/mom;
}
if(iarm==0 && ich==1){
dzmean = 1.27852+-0.00405691*mom+-0.727664/mom+0.319641/sqrt(mom)+-0.0973805/mom/mom+0.375486/mom/mom/mom+-0.130673/mom/mom/mom/mom;
}
if(iarm==1 && ich==0){
dzmean = 1.3102+0.0466208*mom+0.0793568/mom+0.97011/sqrt(mom)+-0.373307/mom/mom+0.0220716/mom/mom/mom+0.0311447/mom/mom/mom/mom;
}
if(iarm==1 && ich==1){
dzmean = 1.42052+0.0950962*mom+0.292751/mom+0.246439/sqrt(mom)+0.723935/mom/mom+-1.04107/mom/mom/mom+0.322048/mom/mom/mom/mom;
}
if(iarm==0 && ich==0){
dzsigma = 0.510622+0.0688104*mom+0.0243863*mom*mom+0.00834748*mom*mom*mom+0.000573259*mom*mom*mom*mom+-0.00103288*mom*mom*mom*mom*mom+0.965039/sqrt(mom)+0.0528075/mom/mom;
}
if(iarm==0 && ich==1){
dzsigma = 0.623527+0.201441*mom+0.00908589*mom*mom+-0.0105551*mom*mom*mom+-0.00250359*mom*mom*mom*mom+0.00117454*mom*mom*mom*mom*mom+0.665861/sqrt(mom)+0.106809/mom/mom;
}
if(iarm==1 && ich==0){
dzsigma = 0.523201+0.0155795*mom+0.0140036*mom*mom+0.00923089*mom*mom*mom+0.00132981*mom*mom*mom*mom+-0.000941257*mom*mom*mom*mom*mom+1.12587/sqrt(mom)+0.0333357/mom/mom;
}
if(iarm==1 && ich==1){
dzsigma = 0.590526+0.100946*mom+0.00642141*mom*mom+-0.000863441*mom*mom*mom+-0.000316806*mom*mom*mom*mom+-5.09038e-06*mom*mom*mom*mom*mom+0.903611/sqrt(mom)+0.0647382/mom/mom;
}
}
 if(RunNumber>=457634 && RunNumber<=458167){ //39 dAu
if(iarm==0 && ich==0){
dzmean = 0.561138+0.105224*mom+5.94803/mom+-2.46042/sqrt(mom)+-6.60173/mom/mom+3.97303/mom/mom/mom+-0.882612/mom/mom/mom/mom;
}
if(iarm==0 && ich==1){
dzmean = 1.74212+-0.145844*mom+-1.49443/mom+0.0444759/sqrt(mom)+0.0506072/mom/mom+0.734216/mom/mom/mom+-0.279868/mom/mom/mom/mom;
}
if(iarm==1 && ich==0){
dzmean = 1.91251+-0.118987*mom+-0.707666/mom+0.686354/sqrt(mom)+-0.190515/mom/mom+0.30761/mom/mom/mom+-0.0875273/mom/mom/mom/mom;
}
if(iarm==1 && ich==1){
dzmean = 2.22972+-0.100269*mom+-0.808292/mom+-0.0759343/sqrt(mom)+0.936153/mom/mom+-0.524453/mom/mom/mom+0.10447/mom/mom/mom/mom;
}
if(iarm==0 && ich==0){
dzsigma = 1.33641+-0.107382*mom+-0.184095*mom*mom+0.0506931*mom*mom*mom+0.036273*mom*mom*mom*mom+-0.010279*mom*mom*mom*mom*mom+0.0920524/sqrt(mom)+0.128164/mom/mom;
}
if(iarm==0 && ich==1){
dzsigma = 1.55255+0.108122*mom+-0.220042*mom*mom+0.0177923*mom*mom*mom+0.0318554*mom*mom*mom*mom+-0.00641271*mom*mom*mom*mom*mom+-0.416948/sqrt(mom)+0.221174/mom/mom;
}
if(iarm==1 && ich==0){
dzsigma = -0.978535+0.0959833*mom+0.225944*mom*mom+-0.0359598*mom*mom*mom+-0.0327582*mom*mom*mom*mom+0.00814467*mom*mom*mom*mom*mom+2.24418/sqrt(mom)+-0.0608286/mom/mom;
}
if(iarm==1 && ich==1){
dzsigma = 1.7225+-0.0187753*mom+-0.203394*mom*mom+0.0249642*mom*mom*mom+0.0227213*mom*mom*mom*mom+-0.00490385*mom*mom*mom*mom*mom+-0.326446/sqrt(mom)+0.191603/mom/mom;
}
}
 if(RunNumber>=456652 && RunNumber<=457298){ //20 dAu
if(iarm==0 && ich==0){
dzmean = -0.977605+0.501572*mom+7.57686/mom+-2.27109/sqrt(mom)+-6.62342/mom/mom+3.25547/mom/mom/mom+-0.63435/mom/mom/mom/mom;
}
if(iarm==0 && ich==1){
dzmean = -0.899973+0.489294*mom+1.27438/mom+0.337394/sqrt(mom)+0.0786158/mom/mom+-0.644518/mom/mom/mom+0.207994/mom/mom/mom/mom;
}
if(iarm==1 && ich==0){
dzmean = 3.5645+-0.613853*mom+-1.78594/mom+0.93803/sqrt(mom)+-0.636048/mom/mom+1.05689/mom/mom/mom+-0.291635/mom/mom/mom/mom;
}
if(iarm==1 && ich==1){
dzmean = -0.221801+0.460199*mom+1.71792/mom+0.609929/sqrt(mom)+0.934703/mom/mom+-1.90396/mom/mom/mom+0.598957/mom/mom/mom/mom;
}
if(iarm==0 && ich==0){
dzsigma = 0.343255+-0.307807*mom+-0.0138482*mom*mom+0.0762099*mom*mom*mom+0.0274379*mom*mom*mom*mom+-0.0141321*mom*mom*mom*mom*mom+1.24671/sqrt(mom)+-0.0386834/mom/mom;
}
if(iarm==0 && ich==1){
dzsigma = 0.299536+-0.293147*mom+0.00209439*mom*mom+0.0812177*mom*mom*mom+0.0268495*mom*mom*mom*mom+-0.0160664*mom*mom*mom*mom*mom+1.17065/sqrt(mom)+-0.0151954/mom/mom;
}
if(iarm==1 && ich==0){
dzsigma = 0.538422+0.139061*mom+-0.0156019*mom*mom+-0.0195822*mom*mom*mom+-0.00381578*mom*mom*mom*mom+0.00263508*mom*mom*mom*mom*mom+0.657148/sqrt(mom)+0.117537/mom/mom;
}
if(iarm==1 && ich==1){
dzsigma = 0.279068+-0.451931*mom+0.00970137*mom*mom+0.131759*mom*mom*mom+0.0423344*mom*mom*mom*mom+-0.0293924*mom*mom*mom*mom*mom+1.39832/sqrt(mom)+-0.036058/mom/mom;
}
}
return (dz-dzmean)/dzsigma;
}
