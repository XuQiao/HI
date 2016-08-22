/*
const double *ptbin = rTotalEff3D->GetYaxis()->GetXbins()->GetArray();
const int npt = gFakPt->GetN();
const double *etabin = rTotalEff3D->GetXaxis()->GetXbins()->GetArray();
const int neta = rTotalEff3D->GetNbinsX();
*/

const int npt = 41;
const double ptbin[npt+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.4,1.6,1.8,2,2.2,2.4,3.2,4,4.8,5.6,6.4,7.2,9.6,12,14.4,19.2,24,28.8,35.2,41.6,48,60.8,73.6,86.4,103.6,120.8,138,155.2,172.4};
const int neta = 30;
const double etabin[neta+1] = {-3,-2.8,-2.6,-2.4,-2.2,-2,-1.8,-1.6,-1.4,-1.2,-1,-0.8,-0.6,-0.4,-0.2,4.47035e-08,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3};
const int nMult=14;
const double Multbin[nMult+1] = {0,1,2,3,4,5,6,7,8,9,10,25,50,200,1000};
const double ptmin = 0.1;
//const int nMulttot=200;
//int Multbintot[nMulttot];
//for(int iMult=0;iMult<nMulttot;iMult++)
//      Multbintot[iMult]=iMult;
