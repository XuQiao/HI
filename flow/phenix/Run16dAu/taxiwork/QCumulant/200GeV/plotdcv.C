#include "SimplifyLife.C"
#include "RpPar.h"
#include <vector>

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

void plotdcv(){
  TGraphErrors* grwd[ncent][nsub][nhar];
  TGraphErrors* grwc[ncent][nsub][nhar];
  TGraphErrors* grwv[ncent][nsub][nhar];
  
  TGraphErrors* grwdpr[ncent][nsub][nhar][ncorr];
  TGraphErrors* grwcpr[ncent][nsub][nhar][ncorr];
  TGraphErrors* grwvpr[ncent][nsub][nhar][ncorr];

  TFile *fin = TFile::Open("outwP.root");
  for(int icent=0;icent<ncent;icent++){
    for(int ihar=0;ihar<nhar;ihar++){
      for(int isub=0;isub<nsub;isub++){
  grwd[icent][isub][ihar] = (TGraphErrors*)fin->Get(Form("grwd_%d_%d_%d",icent,isub,ihar));
  grwc[icent][isub][ihar] = (TGraphErrors*)fin->Get(Form("grwc_%d_%d_%d",icent,isub,ihar));
  grwv[icent][isub][ihar] = (TGraphErrors*)fin->Get(Form("grwv_%d_%d_%d",icent,isub,ihar));
        for(int icorr=0;icorr<ncorr;icorr++){
  grwdpr[icent][isub][ihar][icorr] = (TGraphErrors*)fin->Get(Form("grwdpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
  grwcpr[icent][isub][ihar][icorr] = (TGraphErrors*)fin->Get(Form("grwcpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
  grwvpr[icent][isub][ihar][icorr] = (TGraphErrors*)fin->Get(Form("grwvpr_%d_%d_%d_%d",icent,isub,ihar,icorr));
        }
      }
    }
  }
  
  int iharE=0;
  if(nhar==1 || nhar==2) iharE = 1.0;

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
    leg->AddEntry(grwd[xcent][isub][xhar],Form("<c%d> cent:%d-%d%%, %s",n,0,5,choosesub(isub).Data()),"P");
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
    leg->AddEntry(grwv[xcent][isub][xhar],Form("<v%d> cent:%d-%d%%, %s",n,0,5,choosesub(isub).Data()),"P");
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
  tl->SetTextSize(0.035);
  tl->SetNDC();
  tl->DrawLatex(0.4,0.7,Form("%d-particle correlation differential",2*(xcorr+1)));
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwdpr[xcent][isub][xhar][xcorr],1.2,color[isub],style[isub]);
    grwdpr[xcent][isub][xhar][xcorr]->Draw("Psame");
    leg->AddEntry(grwdpr[xcent][isub][xhar][xcorr],Form("<d%d> cent:%d-%d%%, %s",n,0,5,choosesub(isub).Data()),"P");
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
  tl->SetTextSize(0.035);
  tl->SetNDC();
  tl->DrawLatex(0.4,0.7,Form("%d-particle correlation differential",2*(xcorr+1)));
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwcpr[xcent][isub][xhar][xcorr],1.2,color[isub],style[isub]);
    grwcpr[xcent][isub][xhar][xcorr]->Draw("Psame");
    leg->AddEntry(grwcpr[xcent][isub][xhar][xcorr],Form("<c%d> cent:%d-%d%%, %s",n,0,5,choosesub(isub).Data()),"P");
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
  tl->SetTextSize(0.035);
  tl->SetNDC();
  tl->DrawLatex(0.5,0.8,Form("%d-particle correlation differential",2*(xcorr+1)));
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwvpr[xcent][isub][xhar][xcorr],1.2,color[isub],style[isub]);
    grwvpr[xcent][isub][xhar][xcorr]->Draw("Psame");
    leg->AddEntry(grwvpr[xcent][isub][xhar][xcorr],Form("<v%d>, cent: %d-%d%%, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c6->Print(Form("vpr%d_cent%d.png",n,xcent));
}
