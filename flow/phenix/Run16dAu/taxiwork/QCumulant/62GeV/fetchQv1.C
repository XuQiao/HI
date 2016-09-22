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

double fetchsub(TFile *fin, int icent, int isub, int ihar, int icorr, int ipt, int choosedcv){
    TH1F *hwd0 = (TH1F*)fin->Get(Form("hwd_%d_%d_%d_%d",icent,isub,ihar,0));
    TH1F *hwd1 = (TH1F*)fin->Get(Form("hwd_%d_%d_%d_%d",icent,isub,ihar,1));
    TH1F *hwd2 = (TH1F*)fin->Get(Form("hwd_%d_%d_%d_%d",icent,isub,ihar,2));
    if(ipt>=0){
    TH1F *hwdpr0 = (TH1F*)fin->Get(Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,0));
    TH1F *hwdpr1 = (TH1F*)fin->Get(Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,1));
    TH1F *hwdpr2 = (TH1F*)fin->Get(Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,2));
    }
    double wd0 = hwd0->GetMean();
    double wd1 = hwd1->GetMean();
    double wd2 = hwd2->GetMean();
    double wc0 = wd0;
    double wc1 = wd1-2*wd0*wd0;
    double wc2 = wd2-9*wd0*wd1+12*wd0*wd0*wd0;
    double wv0 = wc0>=0?sqrt(wc0):-9999;
    double wv1 = wc1<=0?sqrt(sqrt(-wc1)):-9999;
    double wv2 = wc2>=0?sqrt(sqrt(sqrt(1./4*wc2))):-9999;
    if(ipt>=0){
    double wdpr0 = hwdpr0->GetMean();
    double wdpr1 = hwdpr1->GetMean();
    double wdpr2 = hwdpr2->GetMean();
    double wcpr0 = wdpr0;
    double wcpr1 = wdpr1-2*wdpr0*wdpr0;
    double wcpr2 = wdpr2-9*wdpr0*wdpr1+12*wdpr0*wdpr0*wdpr0;
    double wvpr0 = wcpr0/wv0;
    double wvpr1 = -wcpr1/wv1/wv1/wv1;
    double wvpr2 = -9999;
    }
 //   fin->Close();
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
  TFile *fin = TFile::Open("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/work/62GeV/output_Qcumu.root");
  TFile *fins[nsample];
  for(int isample=0;isample<nsample;isample++){
    fins[isample] = TFile::Open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/work/62GeV/output_Qcumu_%d.root",isample));
  }

  TH1F* hd[ncent][nsub][nhar][ncorr]; 
  TH1F* hc[ncent][nsub][nhar][ncorr];
  TH1F* hv[ncent][nsub][nhar][ncorr];
  
  TH1F* hwd[ncent][nsub][nhar][ncorr];
  TH1F* hwc[ncent][nsub][nhar][ncorr];
  TH1F* hwv[ncent][nsub][nhar][ncorr];

  double wd[ncent][nsub][nhar][ncorr];
  double wc[ncent][nsub][nhar][ncorr];
  double wv[ncent][nsub][nhar][ncorr];
  double wderr[ncent][nsub][nhar][ncorr];
  double wcerr[ncent][nsub][nhar][ncorr];
  double wverr[ncent][nsub][nhar][ncorr];

  TH1F* hdpr[ncent][nsub][nhar][npt][ncorr];
  TH1F* hcpr[ncent][nsub][nhar][npt][ncorr];
  TH1F* hvpr[ncent][nsub][nhar][npt][ncorr];
  
  TH1F* hwdpr[ncent][nsub][nhar][npt][ncorr];
  TH1F* hwcpr[ncent][nsub][nhar][npt][ncorr];
  TH1F* hwvpr[ncent][nsub][nhar][npt][ncorr];
  
  double wdpr[ncent][nsub][nhar][ncorr][npt];
  double wcpr[ncent][nsub][nhar][ncorr][npt];
  double wvpr[ncent][nsub][nhar][ncorr][npt];
  double wdprerr[ncent][nsub][nhar][ncorr][npt];
  double wcprerr[ncent][nsub][nhar][ncorr][npt];
  double wvprerr[ncent][nsub][nhar][ncorr][npt];
  
  TGraphErrors* grwd[ncent][nsub][nhar];
  TGraphErrors* grwc[ncent][nsub][nhar];
  TGraphErrors* grwv[ncent][nsub][nhar];
  
  TGraphErrors* grwdpr[ncent][nsub][nhar][ncorr];
  TGraphErrors* grwcpr[ncent][nsub][nhar][ncorr];
  TGraphErrors* grwvpr[ncent][nsub][nhar][ncorr];

  int iharE=0;
  if(nhar==1 || nhar==2) iharE = 1.0;

  for(int icent=0;icent<ncent;icent++){
    for(int ihar=0;ihar<nhar;ihar++){
        int n = ihar + 1 + iharE;
        ofstream f(Form("Result/vraw%d_cent%d.dat",n,icent));
      for(int isub=0;isub<nsub;isub++){
        for(int icorr=0;icorr<ncorr;icorr++){
        f<<"isub = "<<isub<< " " << 2*(icorr+1) << " - particle correlation:"<<endl;
        hd[icent][isub][ihar][icorr] = (TH1F*)fin->Get(Form("hd_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hc[icent][isub][ihar][icorr] = (TH1F*)fin->Get(Form("hc_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hv[icent][isub][ihar][icorr] = (TH1F*)fin->Get(Form("hv_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hwd[icent][isub][ihar][icorr] = (TH1F*)fin->Get(Form("hwd_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hwc[icent][isub][ihar][icorr] = (TH1F*)fin->Get(Form("hwc_%d_%d_%d_%d",icent,isub,ihar,icorr));
        hwv[icent][isub][ihar][icorr] = (TH1F*)fin->Get(Form("hwv_%d_%d_%d_%d",icent,isub,ihar,icorr));
        f<<" <d> = " << hd[icent][isub][ihar][icorr]->GetMean() <<" <c> = " << hc[icent][isub][ihar][icorr]->GetMean() <<" <v> = " << hv[icent][isub][ihar][icorr]->GetMean() << endl;
        f<<" <wd> = " << hwd[icent][isub][ihar][icorr]->GetMean() <<" <wc> = " << hwc[icent][isub][ihar][icorr]->GetMean() <<" <wv> = " << hwv[icent][isub][ihar][icorr]->GetMean() << endl;
        for(int ipt=0;ipt<npt;ipt++){
            f<<"ipt = "<<ipt<< ", " << ptbin[ipt] << " <= pt < "<< ptbin[ipt+1] << endl;
            hdpr[icent][isub][ihar][ipt][icorr] = (TH1F*)fin->Get(Form("hdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hcpr[icent][isub][ihar][ipt][icorr] = (TH1F*)fin->Get(Form("hcpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hvpr[icent][isub][ihar][ipt][icorr] = (TH1F*)fin->Get(Form("hvpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hwdpr[icent][isub][ihar][ipt][icorr] = (TH1F*)fin->Get(Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hwcpr[icent][isub][ihar][ipt][icorr] = (TH1F*)fin->Get(Form("hwcpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
            hwvpr[icent][isub][ihar][ipt][icorr] = (TH1F*)fin->Get(Form("hwvpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr));
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
        ofstream f(Form("Result/V%d_cent%d.dat",n,icent));
        ofstream fv(Form("Result/v%d_cent%d.dat",n,icent));
      for(int isub=0;isub<nsub;isub++){
        for(int icorr=0;icorr<ncorr;icorr++){
        f<<"isub = "<<isub<< " " << 2*(icorr+1) << " - particle correlation:"<<endl;
        wd[icent][isub][ihar][icorr] = hwd[icent][isub][ihar][icorr]->GetMean();

        wc[icent][isub][ihar][0] = wd[icent][isub][ihar][0];
        wc[icent][isub][ihar][1] = wd[icent][isub][ihar][1]-2*wd[icent][isub][ihar][0]*wd[icent][isub][ihar][0];
        wc[icent][isub][ihar][2] = wd[icent][isub][ihar][2]-9*wd[icent][isub][ihar][0]*wd[icent][isub][ihar][1]+12*wd[icent][isub][ihar][0]*wd[icent][isub][ihar][0]*wd[icent][isub][ihar][0];
        wv[icent][isub][ihar][0] = wc[icent][isub][ihar][0]>=0?sqrt(wc[icent][isub][ihar][0]):-9999;
        wv[icent][isub][ihar][1] = wc[icent][isub][ihar][1]<=0?sqrt(sqrt(-wc[icent][isub][ihar][1])):-9999;
        wv[icent][isub][ihar][2] = wc[icent][isub][ihar][2]>=0?sqrt(sqrt(sqrt(1./4*wc[icent][isub][ihar][2]))):-9999;


        vector<double> vecd;
        vector<double> vecc;
        vector<double> vecv;
        for(int isample=0;isample<nsample;isample++){
            vecd.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, -1, 0));
            vecc.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, -1, 1));
            vecv.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, -1, 2));
        }
        wderr[icent][isub][ihar][icorr] = getsigma(vecd);
        wcerr[icent][isub][ihar][icorr] = getsigma(vecc);
        wverr[icent][isub][ihar][icorr] = getsigma(vecv);

        f<<" <d> = " << wd[icent][isub][ihar][icorr] <<" <c> = " << wc[icent][isub][ihar][icorr] <<" <v> = " << wv[icent][isub][ihar][icorr] << endl;
        for(int ipt=0;ipt<npt;ipt++){
        fv<<"isub = "<<isub<< " " << 2*(icorr+1) << " - particle correlation:"<<endl;
        wdpr[icent][isub][ihar][icorr][ipt] = hwdpr[icent][isub][ihar][ipt][icorr]->GetMean();

        wcpr[icent][isub][ihar][0][ipt] = wdpr[icent][isub][ihar][0][ipt];
        wcpr[icent][isub][ihar][1][ipt] = wdpr[icent][isub][ihar][1][ipt]-2*wdpr[icent][isub][ihar][0][ipt]*wd[icent][isub][ihar][0];
        wvpr[icent][isub][ihar][0][ipt] = wcpr[icent][isub][ihar][0][ipt]/wv[icent][isub][ihar][0];
        wvpr[icent][isub][ihar][1][ipt] = -wcpr[icent][isub][ihar][1][ipt]/wv[icent][isub][ihar][1]/wv[icent][isub][ihar][1]/wv[icent][isub][ihar][1];
        
        vector<double> vecdpr;
        vector<double> veccpr;
        vector<double> vecvpr;
        for(int isample=0;isample<nsample;isample++){
            vecdpr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, ipt, 0));
            veccpr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, ipt, 1));
            vecvpr.push_back(fetchsub(fins[isample], icent, isub, ihar, icorr, ipt, 2));
        }
        wdprerr[icent][isub][ihar][icorr][ipt] = getsigma(vecdpr);
        wcprerr[icent][isub][ihar][icorr][ipt] = getsigma(veccpr);
        wvprerr[icent][isub][ihar][icorr][ipt] = getsigma(vecvpr);


        fv<<" <d> = " << wdpr[icent][isub][ihar][icorr][ipt] <<" <c> = " << wcpr[icent][isub][ihar][icorr][ipt] <<" <v> = " << wvpr[icent][isub][ihar][icorr][ipt] << endl;
        }
        grwdpr[icent][isub][ihar][icorr] = new TGraphErrors(npt,ptmean,wdpr[icent][isub][ihar][icorr],0,wdprerr[icent][isub][ihar][icorr]);
        grwdpr[icent][isub][ihar][icorr]->SetName(Form("grwdpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        grwcpr[icent][isub][ihar][icorr] = new TGraphErrors(npt,ptmean,wcpr[icent][isub][ihar][icorr],0,wcprerr[icent][isub][ihar][icorr]);
        grwcpr[icent][isub][ihar][icorr]->SetName(Form("grwcpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        grwvpr[icent][isub][ihar][icorr] = new TGraphErrors(npt,ptmean,wvpr[icent][isub][ihar][icorr],0,wvprerr[icent][isub][ihar][icorr]);
        grwvpr[icent][isub][ihar][icorr]->SetName(Form("grwvpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        }
        grwd[icent][isub][ihar] = new TGraphErrors(ncorr,corrbin,wd[icent][isub][ihar],0,wderr[icent][isub][ihar]);
        grwd[icent][isub][ihar]->SetName(Form("grwd_%d_%d_%d",icent,isub,ihar));
        grwc[icent][isub][ihar] = new TGraphErrors(ncorr,corrbin,wc[icent][isub][ihar],0,wcerr[icent][isub][ihar]);
        grwc[icent][isub][ihar]->SetName(Form("grwc_%d_%d_%d",icent,isub,ihar));
        grwv[icent][isub][ihar] = new TGraphErrors(ncorr,corrbin,wv[icent][isub][ihar],0,wverr[icent][isub][ihar]);
        grwv[icent][isub][ihar]->SetName(Form("grwv_%d_%d_%d",icent,isub,ihar));
      }
    }
  }

  fin->Close();
  for(int isample=0;isample<nsample;isample++){
      fins[isample]->Close();
  }
  TH1F* h = new TH1F("","",100,0,10);
  int color[12] = {1,2,4,7,8,5,1,2,4,7,8,5};
  int style[12] = {20,21,24,25,26,27,29,30,31,32,33,34};
  TCanvas *c1 = new TCanvas();
  c1->cd();
  int xcent = 0;
  int xhar = 1;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"n-particle correlation", Form("<d%d>",n),"");
  SetXRange(*h,0,8);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwd[xcent][isub][xhar],1.2,color[isub],style[isub]);
    grwd[xcent][isub][xhar]->Draw("Psame");
    leg->AddEntry(grwd[xcent][isub][xhar],Form("<d%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c1->Print(Form("d%d_cent%d.png",n,xcent));

  TCanvas *c2 = new TCanvas();
  c2->cd();
  int xcent = 0;
  int xhar = 1;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"n-particle correlation", Form("<c%d>",n),"");
  SetXRange(*h,0,8);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwc[xcent][isub][xhar],1.2,color[isub],style[isub]);
    grwc[xcent][isub][xhar]->Draw("Psame");
    leg->AddEntry(grwd[xcent][isub][xhar],Form("<c%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c2->Print(Form("c%d_cent%d.png",n,xcent));

  TCanvas *c3 = new TCanvas();
  c3->cd();
  int xcent = 0;
  int xhar = 1;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"n-particle correlation", Form("<v%d>",n),"");
  SetXRange(*h,0,8);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwv[xcent][isub][xhar],1.2,color[isub],style[isub]);
    grwv[xcent][isub][xhar]->Draw("Psame");
    leg->AddEntry(grwv[xcent][isub][xhar],Form("<v%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c3->Print(Form("v%d_cent%d.png",n,xcent));

  TCanvas *c4 = new TCanvas();
  c4->cd();
  int xcent = 0;
  int xhar = 1;
  int xcorr = 1;
  int n = xhar + 1 + iharE;
  TH1F* h = new TH1F("","",50,0,5);
  SetTitle(*h,"p_{T} (GeV/c)", Form("<d%d>",n),"");
  SetXRange(*h,0,5);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  TLatex *tl = new TLatex();
  tl->SetNDC();
  tl->DrawLatex(0.4,0.7,Form("%d-particle correlation differential",2*(xcorr+1)));
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwdpr[xcent][isub][xhar][xcorr],1.2,color[isub],style[isub]);
    grwdpr[xcent][isub][xhar][xcorr]->Draw("Psame");
    leg->AddEntry(grwdpr[xcent][isub][xhar][xcorr],Form("<d%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c4->Print(Form("dpr%d_cent%d.png",n,xcent));
  
  TCanvas *c5 = new TCanvas();
  c5->cd();
  int xcent = 0;
  int xhar = 1;
  int xcorr = 1;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"p_{T} (GeV/c)", Form("<c%d>",n),"");
  SetXRange(*h,0,5);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  TLatex *tl = new TLatex();
  tl->SetNDC();
  tl->DrawLatex(0.4,0.7,Form("%d-particle correlation differential",2*(xcorr+1)));
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwcpr[xcent][isub][xhar][xcorr],1.2,color[isub],style[isub]);
    grwcpr[xcent][isub][xhar][xcorr]->Draw("Psame");
    leg->AddEntry(grwcpr[xcent][isub][xhar][xcorr],Form("<c%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c5->Print(Form("cpr%d_cent%d.png",n,xcent));
  
  TCanvas *c6 = new TCanvas();
  c6->cd();
  int xcent = 0;
  int xhar = 1;
  int xcorr = 1;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"p_{T} (GeV/c)", Form("<v%d>",n),"");
  SetXRange(*h,0,5);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.1,0.6,0.45,0.9);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  TLatex *tl = new TLatex();
  tl->SetTextSize(0.04);
  tl->SetNDC();
  tl->DrawLatex(0.5,0.75,Form("cent: %d-%d",0,5));
  tl->DrawLatex(0.5,0.8,Form("%d-particle correlation differential",2*(xcorr+1)));
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwvpr[xcent][isub][xhar][xcorr],1.2,color[isub],style[isub]);
    grwvpr[xcent][isub][xhar][xcorr]->Draw("Psame");
    leg->AddEntry(grwvpr[xcent][isub][xhar][xcorr],Form("<v%d>, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c6->Print(Form("vpr%d_cent%d.png",n,xcent));
}
