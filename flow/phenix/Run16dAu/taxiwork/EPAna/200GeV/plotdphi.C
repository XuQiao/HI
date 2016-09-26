#include "RpPar.h"
void plotdphi(){
    gStyle->SetOptStat(kFALSE);
    int ihar = 2;
    int icent = 0;
    int isub = 5;
    int iharE = 0;
    float ptmin  = 1.0;
    float ptmax  = 3.0;
     if(nhar==1 || nhar ==2) iharE=1;
     int n = ihar+1.0+iharE;
    TH1::SetDefaultSumw2();
    //const int np = 13;
    if(n==2){
    const int np = 26;
    float start = 1.56;
    }
    if(n==3){
    const int np = 13;
    float start = 1.04;
    }
    double xbin[np+1+2];
    //float start = 1.04;
    // {-1.56,-1.-1.2,-1.05,-0.9,-0.75,-0.6,-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5};
    for(int i=1;i<np+2;i++){
        xbin[i]=-start + (i-1) * (start * 2) / np;
    }
    xbin[0] = -4./n;
    xbin[1] = -start - 1e-5;
    xbin[np+1] = start + 1e-5;
    xbin[np+2] = 4./n;
    TFile *fin = TFile::Open(Form("dphiv%d.root",n));
  //  for(int isub=0;isub<nsub;isub++){
        TH2F* hpt = (TH2F*)fin->Get(Form("hdphinall_00_%d_%d_%d",icent,ihar,isub));
        int xbinmin = hpt->GetXaxis()->FindBin(ptmin);
        int xbinmax = hpt->GetXaxis()->FindBin(ptmax);
        TH1F* h1 = (TH1F*)hpt->ProjectionY("h1",xbinmin,xbinmax);
        //h = remake(h1);
        //h->Rebin(5);
        h1->Rebin(np+2,"h",xbin);
        int ybinmin = h->GetXaxis()->FindBin(-TMath::Pi()/n)+1;
        int ybinmax = h->GetXaxis()->FindBin(TMath::Pi()/n)-1;
        h->GetXaxis()->SetRange(ybinmin,ybinmax);
        h->SetBinContent(ybinmin-1,0);
        h->SetBinContent(ybinmax+1,0);
        cout<<ybinmin<<" "<<ybinmax<<endl;
        float Norm = h->Integral(ybinmin,ybinmax)/(ybinmax-ybinmin+1);
        cout<<Norm<<endl;
        h->Scale(1./Norm);
   // }
    TF1 *fun = new TF1("fun",Form("[0]+[1]*cos(%d*(x-[3]))+[2]*cos(2*%d*(x-[3]))",n,n),-TMath::Pi()/n,TMath::Pi()/n);
    fun->FixParameter(3,0);
    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    h->Fit("fun","R");
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.04);
    h->SetMarkerSize(1.2);
    h->SetMarkerStyle(20);
    h->SetTitle("");
    float min = h->GetMinimum();
    float max = h->GetMaximum();
    if(n==3)
    h->GetYaxis()->SetRangeUser(0.99, 1.01);
    if(n==2)
    h->GetYaxis()->SetRangeUser(0.8, 1.2);
    //h->GetYaxis()->SetRangeUser(0., 2);
    h->GetXaxis()->SetRangeUser(-4./n,4./n);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetTitle(Form("#phi-#Psi_{%d}",n));
    h->GetYaxis()->SetTitle("Arbitrary scale");
    h->Draw("PE");
    t.DrawLatex(0.18,0.85,"Run 16 d+Au 200 GeV, 0-5\% Central");
    t.DrawLatex(0.18,0.75,"|#eta| < 0.35");
    t.DrawLatex(0.18,0.70,"1.0 < p_{T} < 3.0 GeV/c");
    t.SetTextColor(2);
    if(n==2)
    t.DrawLatex(0.55,0.72,Form("Secord order plane"));
    if(n==3)
    t.DrawLatex(0.55,0.72,Form("third order plane"));
    TLegend *leg = new TLegend(0.12,0.15,0.8,0.3);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->SetFillColor(0);
    leg->AddEntry(fun,Form("fit function: [0]+[1]*cos(%dx)+[2]*cos(%dx)",n,2*n),"l");
    leg->Draw("same");
    c1->Print(Form("dphin_%d_%d_%d.png",n,icent,isub));
}

TH1F* remake(TH1F *h){
    int n = 0;
    int n1 = 0;
    for(int i=1;i<h->GetNbinsX();i++){
        if(h->GetBinContent(i)!=0)  n++;
        if(n==1) n1 = i;
    }
    cout<<"n1 = " << n1 <<endl;
     n = n - 2;
    TH1F *h1 = new TH1F(Form("%s",h->GetName()),h->GetTitle(),n,h->GetBinLowEdge(n1+1),h->GetBinLowEdge(n1+n+1));
    for(int i=1;i<h1->GetNbinsX();i++){
        h1->SetBinContent(i,h->GetBinContent(n1+i));
    }
    return h1;
}

