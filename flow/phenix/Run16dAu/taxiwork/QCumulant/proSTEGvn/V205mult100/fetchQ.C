#include "SimplifyLife.C"
#include "RpPar.h"
#include <vector>

double getsigma(vector<double> s){
    double sum = 0;
    double sum2 = 0;
    int n = s.size();
    for(int i=0;i<n;i++){
        sum += s[i];
        sum2 += s[i]*s[i];
    }
    return sqrt(sum2/n-sum/n*sum/n);
}

double fetchsub(TFile *fin, int icent, int isub, int ihar, int icorr, int ipt, int choosedcv, int iscorr){
    TH1D *hwd0 = (TH1D*)fin->Get(Form("hwd_%d_%d_%d_%d",icent,isub,ihar,0));
    TH1D *hwd1 = (TH1D*)fin->Get(Form("hwd_%d_%d_%d_%d",icent,isub,ihar,1));
    TH1D *hwd2 = (TH1D*)fin->Get(Form("hwd_%d_%d_%d_%d",icent,isub,ihar,2));
    TH1D *hwcos1 = (TH1D*)fin->Get(Form("hwcos1_%d_%d_%d",icent,isub,ihar));
    TH1D *hwsin1 = (TH1D*)fin->Get(Form("hwsin1_%d_%d_%d",icent,isub,ihar));
    TH1D *hwcos1p2 = (TH1D*)fin->Get(Form("hwcos1p2_%d_%d_%d",icent,isub,ihar));
    TH1D *hwsin1p2 = (TH1D*)fin->Get(Form("hwsin1p2_%d_%d_%d",icent,isub,ihar));
    TH1D *hwcos1m2m3 = (TH1D*)fin->Get(Form("hwcos1m2m3_%d_%d_%d",icent,isub,ihar));
    TH1D *hwsin1m2m3 = (TH1D*)fin->Get(Form("hwsin1m2m3_%d_%d_%d",icent,isub,ihar));
    
    if(ipt>=0){
    TH1D *hwdpr0 = (TH1D*)fin->Get(Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,0));
    TH1D *hwdpr1 = (TH1D*)fin->Get(Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,1));
    TH1D *hwdpr2 = (TH1D*)fin->Get(Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,2));
    TH1D *hwcos1pr = (TH1D*)fin->Get(Form("hwcos1pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
    TH1D *hwsin1pr = (TH1D*)fin->Get(Form("hwsin1pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
    TH1D *hwcos1p2pr = (TH1D*)fin->Get(Form("hwcos1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
    TH1D *hwsin1p2pr = (TH1D*)fin->Get(Form("hwsin1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
    TH1D *hwcos1p2m3pr = (TH1D*)fin->Get(Form("hwcos1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
    TH1D *hwsin1p2m3pr = (TH1D*)fin->Get(Form("hwsin1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
    TH1D *hwcos1m2m3pr = (TH1D*)fin->Get(Form("hwcos1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
    TH1D *hwsin1m2m3pr = (TH1D*)fin->Get(Form("hwsin1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
    }
    
    double wd0 = hwd0->GetMean();
    double wd1 = hwd1->GetMean();
    double wd2 = hwd2->GetMean();
    
    double wcos1 = hwcos1->GetMean();
    double wsin1 = hwsin1->GetMean();
    double wcos1p2 = hwcos1p2->GetMean();
    double wsin1p2 = hwsin1p2->GetMean();
    double wcos1m2m3 = hwcos1m2m3->GetMean();
    double wsin1m2m3 = hwsin1m2m3->GetMean();
    
    double wc0 = wd0;
    double wc1 = wd1-2*wd0*wd0;
    double wc2 = wd2-9*wd0*wd1+12*wd0*wd0*wd0;
    double wv0 = wc0>=0?sqrt(wc0):-9999;
    double wv1 = wc1<=0?sqrt(sqrt(-wc1)):-9999;
    double wv2 = wc2>=0?sqrt(sqrt(sqrt(1./4*wc2))):-9999;

    double wdcorr0 = wd0;
    double wdcorr1 = wd1;
    double wdcorr2 = wd2;
    double wccorr0 = wdcorr0-wcos1*wcos1-wsin1*wsin1;
    double wccorr1 = wdcorr1-2*wdcorr0*wdcorr0-4*wcos1*wcos1m2m3+4*wsin1*wsin1m2m3-wcos1p2*wcos1p2-wsin1p2*wsin1p2+4*wcos1p2*(wcos1*wcos1-wsin1*wsin1)+8*wsin1p2*wsin1*wcos1+8*wd0*(wcos1*wcos1+wsin1*wsin1)-6*(wcos1*wcos1+wsin1*wsin1)*(wcos1*wcos1+wsin1*wsin1);
    double wccorr2 = wdcorr2-9*wdcorr0*wdcorr1+12*wdcorr0*wdcorr0*wdcorr0;
    double wvcorr0 = wccorr0>=0?sqrt(wccorr0):-9999;
    double wvcorr1 = wccorr1<=0?sqrt(sqrt(-wccorr1)):-9999;
    double wvcorr2 = wccorr2>=0?sqrt(sqrt(sqrt(1./4*wccorr2))):-9999;

    if(ipt>=0){
    double wdpr0 = hwdpr0->GetMean();
    double wdpr1 = hwdpr1->GetMean();
    double wdpr2 = hwdpr2->GetMean();
    
    double wcos1pr = hwcos1pr->GetMean();
    double wsin1pr = hwsin1pr->GetMean();
    double wcos1p2pr = hwcos1p2pr->GetMean();
    double wsin1p2pr = hwsin1p2pr->GetMean();
    double wcos1p2m3pr = hwcos1p2m3pr->GetMean();
    double wsin1p2m3pr = hwsin1p2m3pr->GetMean();
    double wcos1m2m3pr = hwcos1m2m3pr->GetMean();
    double wsin1m2m3pr = hwsin1m2m3pr->GetMean();

    double wcpr0 = wdpr0;
    double wcpr1 = wdpr1-2*wdpr0*wdpr0;
    double wcpr2 = wdpr2-9*wdpr0*wdpr1+12*wdpr0*wdpr0*wdpr0;
    double wvpr0 = wcpr0/wv0;
    double wvpr1 = -wcpr1/wv1/wv1/wv1;
    double wvpr2 = -9999;
    
    double wdcorrpr0 = wdpr0;
    double wdcorrpr1 = wdpr1;
    double wdcorrpr2 = wdpr2;
    double wccorrpr0 = wdcorrpr0-wcos1pr*wcos1-wsin1pr*wsin1;
    double wccorrpr1 = wdcorrpr1-2*wdcorrpr0*wdcorr0 
-wcos1pr*wcos1m2m3
+wsin1pr*wsin1m2m3
-wcos1*wcos1m2m3pr
+wsin1*wsin1m2m3pr
-2*wcos1*wcos1p2m3pr
-2*wsin1*wsin1p2m3pr
-wcos1p2pr*wcos1p2
-wsin1p2pr*wsin1p2
+2*wcos1p2*(wcos1pr*wcos1-wsin1pr*wsin1)+2*wsin1p2*(wcos1pr*wsin1+wsin1pr*wcos1)
+4*wdcorr0*(wcos1pr*wcos1+wsin1pr*wsin1)
+2*wcos1p2pr*(wcos1*wcos1-wsin1*wsin1)
+4*wsin1p2pr*wcos1*wsin1
+4*wdcorrpr0*(wcos1*wcos1+wsin1*wsin1)
-6*(wcos1*wcos1-wsin1*wsin1)*
 (wcos1pr*wcos1-wsin1pr*wsin1)
-12*wcos1*wsin1*
 (wsin1pr*wcos1-wcos1pr*wsin1);
    double wccorrpr2 = wdcorrpr2-9*wdcorrpr0*wdcorrpr1+12*wdcorrpr0*wdcorrpr0*wdcorrpr0;
    double wvcorrpr0 = wccorrpr0/wvcorr0;
    double wvcorrpr1 = -wccorrpr1/wvcorr1/wvcorr1/wvcorr1;
    double wvcorrpr2 = -9999;
    }
 //   fin->Close();
 if(!iscorr){
    if(ipt<0){
        if(icorr==0 && choosedcv == 0)
            return wd0;
        else if(icorr==1 && choosedcv == 0)
            return wd1;
        else if(icorr==2 && choosedcv == 0)
            return wd2;
        else if(icorr==0  && choosedcv == 1)
            return wc0;
        else if(icorr==1  && choosedcv == 1)
            return wc1;
        else if(icorr==2  && choosedcv == 1)
            return wc2;
        else if(icorr==0  && choosedcv == 2)
            return wv0;
        else if(icorr==1  && choosedcv == 2)
            return wv1;
        else if(icorr==2  && choosedcv == 2)
            return wv2;
    }
    else{
        if(icorr==0 && choosedcv == 0)
            return wdpr0;
        else if(icorr==1 && choosedcv == 0)
            return wdpr1;
        else if(icorr==2 && choosedcv == 0)
            return wdpr2;
        else if(icorr==0  && choosedcv == 1)
            return wcpr0;
        else if(icorr==1  && choosedcv == 1)
            return wcpr1;
        else if(icorr==2  && choosedcv == 1)
            return wcpr2;
        else if(icorr==0  && choosedcv == 2)
            return wvpr0;
        else if(icorr==1  && choosedcv == 2)
            return wvpr1;
        else if(icorr==2  && choosedcv == 2)
            return wvpr2;
    }
 }
 else{
    if(ipt<0){
        if(icorr==0 && choosedcv == 0)
            return wdcorr0;
        else if(icorr==1 && choosedcv == 0)
            return wdcorr1;
        else if(icorr==2 && choosedcv == 0)
            return wdcorr2;
        else if(icorr==0  && choosedcv == 1)
            return wccorr0;
        else if(icorr==1  && choosedcv == 1)
            return wccorr1;
        else if(icorr==2  && choosedcv == 1)
            return wccorr2;
        else if(icorr==0  && choosedcv == 2)
            return wvcorr0;
        else if(icorr==1  && choosedcv == 2)
            return wvcorr1;
        else if(icorr==2  && choosedcv == 2)
            return wvcorr2;
    }
    else{
        if(icorr==0 && choosedcv == 0)
            return wdcorrpr0;
        else if(icorr==1 && choosedcv == 0)
            return wdcorrpr1;
        else if(icorr==2 && choosedcv == 0)
            return wdcorrpr2;
        else if(icorr==0  && choosedcv == 1)
            return wccorrpr0;
        else if(icorr==1  && choosedcv == 1)
            return wccorrpr1;
        else if(icorr==2  && choosedcv == 1)
            return wccorrpr2;
        else if(icorr==0  && choosedcv == 2)
            return wvcorrpr0;
        else if(icorr==1  && choosedcv == 2)
            return wvcorrpr1;
        else if(icorr==2  && choosedcv == 2)
            return wvcorrpr2;
    }
 }
    return -9999;
}

TString choosesub(int isub){
    TString str;
        if(isub==0)
         str = "FVTX1LS";
        else if(isub==1)
         str = "FVTX2LS";
        else if(isub==2)
         str = "FVTX3LS";
        else if(isub==3)
         str = "FVTX4LS";
        else if(isub==5)
         str = "FVTX1S";
        else if(isub==4)
         str = "BBCS";
        else if(isub==6)
         str = "CNT";
        else
         str = "ABORT";
    return str;
}

void fetchQ(){
  gStyle->SetOptStat(kFALSE);
  TFile *fin = TFile::Open("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/work/200GeV/MCAna_mult100_fullacc.root","ReadOnly");
  TFile *fins[nsample];
  for(int isample=0;isample<nsample;isample++){
    fins[isample] = TFile::Open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/work/200GeV/treeout/MCAna_mult100_fullacc_%d.root",isample),"ReadOnly");
  }
  TFile *fout = new TFile("outP.root","recreate");

  TH1D* hd[ncent][nsub][nhar][ncorr]; 
  TH1D* hc[ncent][nsub][nhar][ncorr];
  TH1D* hv[ncent][nsub][nhar][ncorr];
  
  TH1D* hwd[ncent][nsub][nhar][ncorr];
  TH1D* hwc[ncent][nsub][nhar][ncorr];
  TH1D* hwv[ncent][nsub][nhar][ncorr];

  double wd[ncent][nsub][nhar][ncorr];
  double wc[ncent][nsub][nhar][ncorr];
  double wv[ncent][nsub][nhar][ncorr];
  double wderr[ncent][nsub][nhar][ncorr];
  double wcerr[ncent][nsub][nhar][ncorr];
  double wverr[ncent][nsub][nhar][ncorr];
  
  TH1D* hwcos1[ncent][nsub][nhar]; 
  TH1D* hwsin1[ncent][nsub][nhar];
  TH1D* hwcos1p2[ncent][nsub][nhar]; 
  TH1D* hwsin1p2[ncent][nsub][nhar];
  TH1D* hwcos1m2m3[ncent][nsub][nhar]; 
  TH1D* hwsin1m2m3[ncent][nsub][nhar];

  double wcos1[ncent][nsub][nhar]; 
  double wsin1[ncent][nsub][nhar];
  double wcos1p2[ncent][nsub][nhar]; 
  double wsin1p2[ncent][nsub][nhar];
  double wcos1m2m3[ncent][nsub][nhar]; 
  double wsin1m2m3[ncent][nsub][nhar];

  double wdcorr[ncent][nsub][nhar][ncorr];
  double wccorr[ncent][nsub][nhar][ncorr];
  double wvcorr[ncent][nsub][nhar][ncorr];
  double wdcorrerr[ncent][nsub][nhar][ncorr];
  double wccorrerr[ncent][nsub][nhar][ncorr];
  double wvcorrerr[ncent][nsub][nhar][ncorr];

  TH1D* hdpr[ncent][nsub][nhar][npt][ncorr];
  TH1D* hcpr[ncent][nsub][nhar][npt][ncorr];
  TH1D* hvpr[ncent][nsub][nhar][npt][ncorr];
  
  TH1D* hwdpr[ncent][nsub][nhar][npt][ncorr];
  TH1D* hwcpr[ncent][nsub][nhar][npt][ncorr];
  TH1D* hwvpr[ncent][nsub][nhar][npt][ncorr];
  
  double wdpr[ncent][nsub][nhar][ncorr][npt];
  double wcpr[ncent][nsub][nhar][ncorr][npt];
  double wvpr[ncent][nsub][nhar][ncorr][npt];
  double wdprerr[ncent][nsub][nhar][ncorr][npt];
  double wcprerr[ncent][nsub][nhar][ncorr][npt];
  double wvprerr[ncent][nsub][nhar][ncorr][npt];

  TH1D* hwcos1pr[ncent][nsub][nhar][npt];
  TH1D* hwsin1pr[ncent][nsub][nhar][npt];
  TH1D* hwcos1p2pr[ncent][nsub][nhar][npt]; 
  TH1D* hwsin1p2pr[ncent][nsub][nhar][npt];
  TH1D* hwcos1p2m3pr[ncent][nsub][nhar][npt]; 
  TH1D* hwsin1p2m3pr[ncent][nsub][nhar][npt];
  TH1D* hwcos1m2m3pr[ncent][nsub][nhar][npt]; 
  TH1D* hwsin1m2m3pr[ncent][nsub][nhar][npt]; 
  
  double wcos1pr[ncent][nsub][nhar][npt]; 
  double wsin1pr[ncent][nsub][nhar][npt];
  double wcos1p2pr[ncent][nsub][nhar][npt]; 
  double wsin1p2pr[ncent][nsub][nhar][npt];
  double wcos1p2m3pr[ncent][nsub][nhar][npt]; 
  double wsin1p2m3pr[ncent][nsub][nhar][npt];
  double wcos1m2m3pr[ncent][nsub][nhar][npt]; 
  double wsin1m2m3pr[ncent][nsub][nhar][npt];
  
  double wdcorrpr[ncent][nsub][nhar][ncorr][npt];
  double wccorrpr[ncent][nsub][nhar][ncorr][npt];
  double wvcorrpr[ncent][nsub][nhar][ncorr][npt];
  double wdcorrprerr[ncent][nsub][nhar][ncorr][npt];
  double wccorrprerr[ncent][nsub][nhar][ncorr][npt];
  double wvcorrprerr[ncent][nsub][nhar][ncorr][npt];

  TGraphErrors* grwd[ncent][nsub][nhar];
  TGraphErrors* grwc[ncent][nsub][nhar];
  TGraphErrors* grwv[ncent][nsub][nhar];
  
  TGraphErrors* grwdcorr[ncent][nsub][nhar];
  TGraphErrors* grwccorr[ncent][nsub][nhar];
  TGraphErrors* grwvcorr[ncent][nsub][nhar];
  
  TGraphErrors* grwdpr[ncent][nsub][nhar][ncorr];
  TGraphErrors* grwcpr[ncent][nsub][nhar][ncorr];
  TGraphErrors* grwvpr[ncent][nsub][nhar][ncorr];
  
  TGraphErrors* grwdcorrpr[ncent][nsub][nhar][ncorr];
  TGraphErrors* grwccorrpr[ncent][nsub][nhar][ncorr];
  TGraphErrors* grwvcorrpr[ncent][nsub][nhar][ncorr];
  

  int iharE=0;
  if(nhar==1 || nhar==2) iharE = 1.0;

  for(int icent=0;icent<ncent;icent++){
    for(int ihar=0;ihar<nhar;ihar++){
        int n = ihar + 1 + iharE;
        ofstream f(Form("Result/vraw%d_cent%d.dat",n,icent));
        ofstream fdet(Form("Result/vdet%d_cent%d.dat",n,icent));
      for(int isub=0;isub<nsub;isub++){
          cout<<"processing icent = "<<icent<<" ihar = "<<ihar<<" isub = "<<isub<<endl;
        fdet<<"isub = "<<isub<<endl;
        hwcos1[icent][isub][ihar]=(TH1D*)fin->Get(Form("hwcos1_%d_%d_%d",icent,isub,ihar));
        hwsin1[icent][isub][ihar]=(TH1D*)fin->Get(Form("hwsin1_%d_%d_%d",icent,isub,ihar));
        hwcos1p2[icent][isub][ihar]=(TH1D*)fin->Get(Form("hwcos1p2_%d_%d_%d",icent,isub,ihar));
        hwsin1p2[icent][isub][ihar]=(TH1D*)fin->Get(Form("hwsin1p2_%d_%d_%d",icent,isub,ihar));
        hwcos1m2m3[icent][isub][ihar]=(TH1D*)fin->Get(Form("hwcos1m2m3_%d_%d_%d",icent,isub,ihar));
        hwsin1m2m3[icent][isub][ihar]=(TH1D*)fin->Get(Form("hwsin1m2m3_%d_%d_%d",icent,isub,ihar));
        fdet<<" <cos1> = " << hwcos1[icent][isub][ihar]->GetMean();
        fdet<<" <sin1> = " << hwcos1[icent][isub][ihar]->GetMean();
        fdet<<" <cos1p2> = " << hwcos1p2[icent][isub][ihar]->GetMean();
        fdet<<" <sin1p2> = " << hwsin1p2[icent][isub][ihar]->GetMean();
        fdet<<" <cos1m2m3> = " << hwcos1m2m3[icent][isub][ihar]->GetMean();
        fdet<<" <sin1m2m3> = " << hwsin1m2m3[icent][isub][ihar]->GetMean()<<endl;
        for(int ipt=0;ipt<npt;ipt++){
        fdet<<"isub = "<<isub<<endl;
        hwcos1pr[icent][isub][ihar][ipt]=(TH1D*)fin->Get(Form("hwcos1pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
        hwsin1pr[icent][isub][ihar][ipt]=(TH1D*)fin->Get(Form("hwsin1pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
        hwcos1p2pr[icent][isub][ihar][ipt]=(TH1D*)fin->Get(Form("hwcos1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
        hwsin1p2pr[icent][isub][ihar][ipt]=(TH1D*)fin->Get(Form("hwsin1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
        hwcos1p2m3pr[icent][isub][ihar][ipt]=(TH1D*)fin->Get(Form("hwcos1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
        hwsin1p2m3pr[icent][isub][ihar][ipt]=(TH1D*)fin->Get(Form("hwsin1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
        hwcos1m2m3pr[icent][isub][ihar][ipt]=(TH1D*)fin->Get(Form("hwcos1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
        hwsin1m2m3pr[icent][isub][ihar][ipt]=(TH1D*)fin->Get(Form("hwsin1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt));
        fdet<<" <cos1pr> = " << hwcos1pr[icent][isub][ihar][ipt]->GetMean();
        fdet<<" <sin1pr> = " << hwsin1pr[icent][isub][ihar][ipt]->GetMean();
        fdet<<" <cos1p2pr> = " << hwcos1p2pr[icent][isub][ihar][ipt]->GetMean();
        fdet<<" <sin1p2pr> = " << hwsin1p2pr[icent][isub][ihar][ipt]->GetMean();
        fdet<<" <cos1p2m3> = " << hwcos1p2m3pr[icent][isub][ihar][ipt]->GetMean();
        fdet<<" <sin1p2m3> = " << hwsin1p2m3pr[icent][isub][ihar][ipt]->GetMean();
        fdet<<" <cos1m2m3> = " << hwcos1m2m3pr[icent][isub][ihar][ipt]->GetMean();
        fdet<<" <sin1m2m3> = " << hwsin1m2m3pr[icent][isub][ihar][ipt]->GetMean();
        }
        for(int icorr=0;icorr<ncorr;icorr++){
        f<<"isub = "<<isub<< " " << 2*(icorr+1) << " - particle correlation:"<<endl;
        hd[icent][isub][ihar][icorr] = (TH1D*)fin->Get(Form("hd_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hc[icent][isub][ihar][icorr] = (TH1D*)fin->Get(Form("hc_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hv[icent][isub][ihar][icorr] = (TH1D*)fin->Get(Form("hv_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hwd[icent][isub][ihar][icorr] = (TH1D*)fin->Get(Form("hwd_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hwc[icent][isub][ihar][icorr] = (TH1D*)fin->Get(Form("hwc_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hwv[icent][isub][ihar][icorr] = (TH1D*)fin->Get(Form("hwv_%d_%d_%d_%d",icent,isub,ihar,icorr));
        f<<" <d> = " << hd[icent][isub][ihar][icorr]->GetMean() <<" <c> = " << hc[icent][isub][ihar][icorr]->GetMean() <<" <v> = " << hv[icent][isub][ihar][icorr]->GetMean() << endl;
        f<<" <wd> = " << hwd[icent][isub][ihar][icorr]->GetMean() <<" <wc> = " << hwc[icent][isub][ihar][icorr]->GetMean() <<" <wv> = " << hwv[icent][isub][ihar][icorr]->GetMean() << endl;
        for(int ipt=0;ipt<npt;ipt++){
            f<<"ipt = "<<ipt<< ", " << ptbin[ipt] << " <= pt < "<< ptbin[ipt+1] << endl;
            hdpr[icent][isub][ihar][ipt][icorr] = (TH1D*)fin->Get(Form("hdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hcpr[icent][isub][ihar][ipt][icorr] = (TH1D*)fin->Get(Form("hcpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hvpr[icent][isub][ihar][ipt][icorr] = (TH1D*)fin->Get(Form("hvpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hwdpr[icent][isub][ihar][ipt][icorr] = (TH1D*)fin->Get(Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hwcpr[icent][isub][ihar][ipt][icorr] = (TH1D*)fin->Get(Form("hwcpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hwvpr[icent][isub][ihar][ipt][icorr] = (TH1D*)fin->Get(Form("hwvpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
        f<<" <d> = " << hdpr[icent][isub][ihar][ipt][icorr]->GetMean() <<" <c> = " << hcpr[icent][isub][ihar][ipt][icorr]->GetMean() <<" <v> = " << hvpr[icent][isub][ihar][ipt][icorr]->GetMean() << endl;
        f<<" <wd> = " << hwdpr[icent][isub][ihar][ipt][icorr]->GetMean() <<" <wc> = " << hwcpr[icent][isub][ihar][ipt][icorr]->GetMean() <<" <wv> = " << hwvpr[icent][isub][ihar][ipt][icorr]->GetMean() << endl;
        }
    }
    }
    }
  }
  double corrbin[ncorr] = {2,4,6};
  double ptmean[npt] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9};

  for(int icent=0;icent<ncent;icent++){
    for(int ihar=0;ihar<nhar;ihar++){
        int n = ihar + 1 + iharE;
        ofstream fV(Form("Result/V%d_cent%d.dat",n,icent));
        ofstream fv(Form("Result/v%d_cent%d.dat",n,icent));
        ofstream fVcorr(Form("Result/Vcorr%d_cent%d.dat",n,icent));
        ofstream fvcorr(Form("Result/vcorr%d_cent%d.dat",n,icent));
      for(int isub=0;isub<nsub;isub++){
        wcos1[icent][isub][ihar] = hwcos1[icent][isub][ihar]->GetMean();
        wsin1[icent][isub][ihar] = hwsin1[icent][isub][ihar]->GetMean();
        wcos1p2[icent][isub][ihar] = hwcos1p2[icent][isub][ihar]->GetMean();
        wsin1p2[icent][isub][ihar] = hwsin1p2[icent][isub][ihar]->GetMean();
        wcos1m2m3[icent][isub][ihar] = hwcos1m2m3[icent][isub][ihar]->GetMean();
        wsin1m2m3[icent][isub][ihar] = hwsin1m2m3[icent][isub][ihar]->GetMean();
        for(int ipt=0;ipt<npt;ipt++){
            wcos1pr[icent][isub][ihar][ipt] = hwcos1pr[icent][isub][ihar][ipt]->GetMean();
            wsin1pr[icent][isub][ihar][ipt] = hwsin1pr[icent][isub][ihar][ipt]->GetMean();
            wcos1p2pr[icent][isub][ihar][ipt] = hwcos1p2pr[icent][isub][ihar][ipt]->GetMean();
            wsin1p2pr[icent][isub][ihar][ipt] = hwsin1p2pr[icent][isub][ihar][ipt]->GetMean();
            wcos1p2m3pr[icent][isub][ihar][ipt] = hwcos1p2m3pr[icent][isub][ihar][ipt]->GetMean();
            wsin1p2m3pr[icent][isub][ihar][ipt] = hwsin1p2m3pr[icent][isub][ihar][ipt]->GetMean();
            wcos1m2m3pr[icent][isub][ihar][ipt] = hwcos1m2m3pr[icent][isub][ihar][ipt]->GetMean();
            wsin1m2m3pr[icent][isub][ihar][ipt] = hwsin1m2m3pr[icent][isub][ihar][ipt]->GetMean();
        }
        for(int icorr=0;icorr<ncorr;icorr++){
        fV<<"isub = "<<isub<< " " << 2*(icorr+1) << " - particle correlation"<<endl;
        fVcorr<<"isub = "<<isub<< " " << 2*(icorr+1) << " - particle correlation corrected:"<<endl;
        wd[icent][isub][ihar][icorr] = hwd[icent][isub][ihar][icorr]->GetMean();

        wc[icent][isub][ihar][0] = wd[icent][isub][ihar][0];
        wc[icent][isub][ihar][1] = wd[icent][isub][ihar][1]-2*wd[icent][isub][ihar][0]*wd[icent][isub][ihar][0];
        wc[icent][isub][ihar][2] = wd[icent][isub][ihar][2]-9*wd[icent][isub][ihar][0]*wd[icent][isub][ihar][1]+12*wd[icent][isub][ihar][0]*wd[icent][isub][ihar][0]*wd[icent][isub][ihar][0];
        wv[icent][isub][ihar][0] = wc[icent][isub][ihar][0]>=0?sqrt(wc[icent][isub][ihar][0]):-9999;
        wv[icent][isub][ihar][1] = wc[icent][isub][ihar][1]<=0?sqrt(sqrt(-wc[icent][isub][ihar][1])):-9999;
        wv[icent][isub][ihar][2] = wc[icent][isub][ihar][2]>=0?sqrt(sqrt(sqrt(1./4*wc[icent][isub][ihar][2]))):-9999;
        
        wdcorr[icent][isub][ihar][icorr] = wd[icent][isub][ihar][icorr];
        
        wccorr[icent][isub][ihar][0] = wdcorr[icent][isub][ihar][0]-wcos1[icent][isub][ihar]*wcos1[icent][isub][ihar]-wsin1[icent][isub][ihar]*wcos1[icent][isub][ihar];
        wccorr[icent][isub][ihar][1] = wdcorr[icent][isub][ihar][1]-2*wdcorr[icent][isub][ihar][0]*wdcorr[icent][isub][ihar][0]-4*wcos1[icent][isub][ihar]*wcos1m2m3[icent][isub][ihar]+4*wsin1[icent][isub][ihar]*wsin1m2m3[icent][isub][ihar]-wcos1p2[icent][isub][ihar]*wcos1p2[icent][isub][ihar]-wsin1p2[icent][isub][ihar]*wsin1p2[icent][isub][ihar]+4*wcos1p2[icent][isub][ihar]*(wcos1[icent][isub][ihar]*wcos1[icent][isub][ihar]-wsin1[icent][isub][ihar]*wsin1[icent][isub][ihar])+8*wsin1p2[icent][isub][ihar]*wsin1[icent][isub][ihar]*wcos1[icent][isub][ihar]+8*wd[icent][isub][ihar][0]*(wcos1[icent][isub][ihar]*wcos1[icent][isub][ihar]+wsin1[icent][isub][ihar]*wsin1[icent][isub][ihar])-6*(wcos1[icent][isub][ihar]*wcos1[icent][isub][ihar]+wsin1[icent][isub][ihar]*wsin1[icent][isub][ihar])*(wcos1[icent][isub][ihar]*wcos1[icent][isub][ihar]+wsin1[icent][isub][ihar]*wsin1[icent][isub][ihar]);
        wccorr[icent][isub][ihar][2] = wdcorr[icent][isub][ihar][2]-9*wdcorr[icent][isub][ihar][0]*wdcorr[icent][isub][ihar][1]+12*wdcorr[icent][isub][ihar][0]*wdcorr[icent][isub][ihar][0]*wdcorr[icent][isub][ihar][0];
        wvcorr[icent][isub][ihar][0] = wccorr[icent][isub][ihar][0]>=0?sqrt(wccorr[icent][isub][ihar][0]):-9999;
        wvcorr[icent][isub][ihar][1] = wccorr[icent][isub][ihar][1]<=0?sqrt(sqrt(-wccorr[icent][isub][ihar][1])):-9999;
        wvcorr[icent][isub][ihar][2] = wccorr[icent][isub][ihar][2]>=0?sqrt(sqrt(sqrt(1./4*wccorr[icent][isub][ihar][2]))):-9999;

        vector<double> vecd;
        vector<double> vecc;
        vector<double> vecv;
        vector<double> vecdcorr;
        vector<double> vecccorr;
        vector<double> vecvcorr;
        for(int isample=0;isample<nsample;isample++){
            vecd.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, -1, 0, 0));
            vecc.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, -1, 1, 0));
            vecv.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, -1, 2, 0));
            
            vecdcorr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, -1, 0, 1));
            vecccorr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, -1, 1, 1));
            vecvcorr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, -1, 2, 1));
        }
        wderr[icent][isub][ihar][icorr] = getsigma(vecd);
        wcerr[icent][isub][ihar][icorr] = getsigma(vecc);
        wverr[icent][isub][ihar][icorr] = getsigma(vecv);
        
        wdcorrerr[icent][isub][ihar][icorr] = getsigma(vecdcorr);
        wccorrerr[icent][isub][ihar][icorr] = getsigma(vecccorr);
        wvcorrerr[icent][isub][ihar][icorr] = getsigma(vecvcorr);

        fV<<" <d> = " << wd[icent][isub][ihar][icorr] <<" <c> = " << wc[icent][isub][ihar][icorr] <<" <v> = " << wv[icent][isub][ihar][icorr] << endl;
        fVcorr<<" <dcorr> = " << wdcorr[icent][isub][ihar][icorr] <<" <ccorr> = " << wccorr[icent][isub][ihar][icorr] <<" <vcorr> = " << wvcorr[icent][isub][ihar][icorr] << endl;

        for(int ipt=0;ipt<npt;ipt++){
        fv<<"isub = "<<isub<< " " << 2*(icorr+1) << " - particle correlation:"<<endl;
        fvcorr<<"isub = "<<isub<< " " << 2*(icorr+1) << " - particle correlation corrected:"<<endl;
        wdpr[icent][isub][ihar][icorr][ipt] = hwdpr[icent][isub][ihar][ipt][icorr]->GetMean();

        wcpr[icent][isub][ihar][0][ipt] = wdpr[icent][isub][ihar][0][ipt];
        wcpr[icent][isub][ihar][1][ipt] = wdpr[icent][isub][ihar][1][ipt]-2*wdpr[icent][isub][ihar][0][ipt]*wd[icent][isub][ihar][0];
        wvpr[icent][isub][ihar][0][ipt] = wcpr[icent][isub][ihar][0][ipt]/wv[icent][isub][ihar][0];
        wvpr[icent][isub][ihar][1][ipt] = -wcpr[icent][isub][ihar][1][ipt]/wv[icent][isub][ihar][1]/wv[icent][isub][ihar][1]/wv[icent][isub][ihar][1];
        
        wdcorrpr[icent][isub][ihar][icorr][ipt] = wdpr[icent][isub][ihar][icorr][ipt];

        wccorrpr[icent][isub][ihar][0][ipt] = wdcorrpr[icent][isub][ihar][0][ipt]-wcos1pr[icent][isub][ihar][ipt]*wcos1[icent][isub][ihar]-wsin1pr[icent][isub][ihar][ipt]*wsin1[icent][isub][ihar];
        wccorrpr[icent][isub][ihar][1][ipt] = wdcorrpr[icent][isub][ihar][1][ipt]-2*wdcorrpr[icent][isub][ihar][0][ipt]*wdcorr[icent][isub][ihar][0] 
-wcos1pr[icent][isub][ihar][ipt]*wcos1m2m3[icent][isub][ihar]
+wsin1pr[icent][isub][ihar][ipt]*wsin1m2m3[icent][isub][ihar]
-wcos1[icent][isub][ihar]*wcos1m2m3pr[icent][isub][ihar][ipt]
+wsin1[icent][isub][ihar]*wsin1m2m3pr[icent][isub][ihar][ipt]
-2*wcos1[icent][isub][ihar]*wcos1p2m3pr[icent][isub][ihar][ipt]
-2*wsin1[icent][isub][ihar]*wsin1p2m3pr[icent][isub][ihar][ipt]
-wcos1p2pr[icent][isub][ihar][ipt]*wcos1p2[icent][isub][ihar]
-wsin1p2pr[icent][isub][ihar][ipt]*wsin1p2[icent][isub][ihar]
+2*wcos1p2[icent][isub][ihar]*(wcos1pr[icent][isub][ihar][ipt]*wcos1[icent][isub][ihar]-wsin1pr[icent][isub][ihar][ipt]*wsin1[icent][isub][ihar])+2*wsin1p2[icent][isub][ihar]*(wcos1pr[icent][isub][ihar][ipt]*wsin1[icent][isub][ihar]+wsin1pr[icent][isub][ihar][ipt]*wcos1[icent][isub][ihar])
+4*wdcorr[icent][isub][ihar][0]*(wcos1pr[icent][isub][ihar][ipt]*wcos1[icent][isub][ihar]+wsin1pr[icent][isub][ihar][ipt]*wsin1[icent][isub][ihar])
+2*wcos1p2pr[icent][isub][ihar][ipt]*(wcos1[icent][isub][ihar]*wcos1[icent][isub][ihar]-wsin1[icent][isub][ihar]*wsin1[icent][isub][ihar])
+4*wsin1p2pr[icent][isub][ihar][ipt]*wcos1[icent][isub][ihar]*wsin1[icent][isub][ihar]
+4*wdpr[icent][isub][ihar][0][ipt]*(wcos1[icent][isub][ihar]*wcos1[icent][isub][ihar]+wsin1[icent][isub][ihar]*wsin1[icent][isub][ihar])
-6*(wcos1[icent][isub][ihar]*wcos1[icent][isub][ihar]-wsin1[icent][isub][ihar]*wsin1[icent][isub][ihar])*
 (wcos1pr[icent][isub][ihar][ipt]*wcos1[icent][isub][ihar]-wsin1pr[icent][isub][ihar][ipt]*wsin1[icent][isub][ihar])
-12*wcos1[icent][isub][ihar]*wsin1[icent][isub][ihar]*
 (wsin1pr[icent][isub][ihar][ipt]*wcos1[icent][isub][ihar]-wcos1pr[icent][isub][ihar][ipt]*wsin1[icent][isub][ihar]);
        wvcorrpr[icent][isub][ihar][0][ipt] = wccorrpr[icent][isub][ihar][0][ipt]/wvcorr[icent][isub][ihar][0];
        wvcorrpr[icent][isub][ihar][1][ipt] = -wccorrpr[icent][isub][ihar][1][ipt]/wvcorr[icent][isub][ihar][1]/wvcorr[icent][isub][ihar][1]/wvcorr[icent][isub][ihar][1];

        vector<double> vecdpr;
        vector<double> veccpr;
        vector<double> vecvpr;
        vector<double> vecdcorrpr;
        vector<double> vecccorrpr;
        vector<double> vecvcorrpr;
        for(int isample=0;isample<nsample;isample++){
            vecdpr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, ipt, 0, 0));
            veccpr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, ipt, 1, 0));
            vecvpr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, ipt, 2, 0));
            
            vecdcorrpr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, ipt, 0, 1));
            vecccorrpr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, ipt, 1, 1));
            vecvcorrpr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, ipt, 2, 1));
        }
        wdprerr[icent][isub][ihar][icorr][ipt] = getsigma(vecdpr);
        wcprerr[icent][isub][ihar][icorr][ipt] = getsigma(veccpr);
        wvprerr[icent][isub][ihar][icorr][ipt] = getsigma(vecvpr);

        wdcorrprerr[icent][isub][ihar][icorr][ipt] = getsigma(vecdcorrpr);
        wccorrprerr[icent][isub][ihar][icorr][ipt] = getsigma(vecccorrpr);
        wvcorrprerr[icent][isub][ihar][icorr][ipt] = getsigma(vecvcorrpr);

        fv<<" <d> = " << wdpr[icent][isub][ihar][icorr][ipt] <<" <c> = " << wcpr[icent][isub][ihar][icorr][ipt] <<" <v> = " << wvpr[icent][isub][ihar][icorr][ipt] << endl;
        fvcorr<<" <dcorr> = " << wdcorrpr[icent][isub][ihar][icorr][ipt] <<" <ccorr> = " << wccorrpr[icent][isub][ihar][icorr][ipt] <<" <vcorr> = " << wvcorrpr[icent][isub][ihar][icorr][ipt] << endl;
        }
        fout->cd();
        grwdpr[icent][isub][ihar][icorr] = new TGraphErrors(npt,ptmean,wdpr[icent][isub][ihar][icorr],0,wdprerr[icent][isub][ihar][icorr]);
        grwdpr[icent][isub][ihar][icorr]->SetName(Form("grwdpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        grwdpr[icent][isub][ihar][icorr]->Write();
        grwcpr[icent][isub][ihar][icorr] = new TGraphErrors(npt,ptmean,wcpr[icent][isub][ihar][icorr],0,wcprerr[icent][isub][ihar][icorr]);
        grwcpr[icent][isub][ihar][icorr]->SetName(Form("grwcpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        grwcpr[icent][isub][ihar][icorr]->Write();
        grwvpr[icent][isub][ihar][icorr] = new TGraphErrors(npt,ptmean,wvpr[icent][isub][ihar][icorr],0,wvprerr[icent][isub][ihar][icorr]);
        grwvpr[icent][isub][ihar][icorr]->SetName(Form("grwvpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        grwvpr[icent][isub][ihar][icorr]->Write();
        
        grwdcorrpr[icent][isub][ihar][icorr] = new TGraphErrors(npt,ptmean,wdcorrpr[icent][isub][ihar][icorr],0,wdcorrprerr[icent][isub][ihar][icorr]);
        grwdcorrpr[icent][isub][ihar][icorr]->SetName(Form("grwdcorrpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        grwdcorrpr[icent][isub][ihar][icorr]->Write();
        grwccorrpr[icent][isub][ihar][icorr] = new TGraphErrors(npt,ptmean,wccorrpr[icent][isub][ihar][icorr],0,wccorrprerr[icent][isub][ihar][icorr]);
        grwccorrpr[icent][isub][ihar][icorr]->SetName(Form("grwccorrpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        grwccorrpr[icent][isub][ihar][icorr]->Write();
        grwvcorrpr[icent][isub][ihar][icorr] = new TGraphErrors(npt,ptmean,wvcorrpr[icent][isub][ihar][icorr],0,wvcorrprerr[icent][isub][ihar][icorr]);
        grwvcorrpr[icent][isub][ihar][icorr]->SetName(Form("grwvcorrpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        grwvcorrpr[icent][isub][ihar][icorr]->Write();
        }
        fout->cd();
        grwd[icent][isub][ihar] = new TGraphErrors(ncorr,corrbin,wd[icent][isub][ihar],0,wderr[icent][isub][ihar]);
        grwd[icent][isub][ihar]->SetName(Form("grwd_%d_%d_%d",icent,isub,ihar));
        grwd[icent][isub][ihar]->Write();
        grwc[icent][isub][ihar] = new TGraphErrors(ncorr,corrbin,wc[icent][isub][ihar],0,wcerr[icent][isub][ihar]);
        grwc[icent][isub][ihar]->SetName(Form("grwc_%d_%d_%d",icent,isub,ihar));
        grwc[icent][isub][ihar]->Write();
        grwv[icent][isub][ihar] = new TGraphErrors(ncorr,corrbin,wv[icent][isub][ihar],0,wverr[icent][isub][ihar]);
        grwv[icent][isub][ihar]->SetName(Form("grwv_%d_%d_%d",icent,isub,ihar));
        grwv[icent][isub][ihar]->Write();
        
        grwdcorr[icent][isub][ihar] = new TGraphErrors(ncorr,corrbin,wdcorr[icent][isub][ihar],0,wdcorrerr[icent][isub][ihar]);
        grwdcorr[icent][isub][ihar]->SetName(Form("grwdcorr_%d_%d_%d",icent,isub,ihar));
        grwdcorr[icent][isub][ihar]->Write();
        grwccorr[icent][isub][ihar] = new TGraphErrors(ncorr,corrbin,wccorr[icent][isub][ihar],0,wccorrerr[icent][isub][ihar]);
        grwccorr[icent][isub][ihar]->SetName(Form("grwccorr_%d_%d_%d",icent,isub,ihar));
        grwccorr[icent][isub][ihar]->Write();
        grwvcorr[icent][isub][ihar] = new TGraphErrors(ncorr,corrbin,wvcorr[icent][isub][ihar],0,wvcorrerr[icent][isub][ihar]);
        grwvcorr[icent][isub][ihar]->SetName(Form("grwvcorr_%d_%d_%d",icent,isub,ihar));
        grwvcorr[icent][isub][ihar]->Write();
      }
    }
  }

}
