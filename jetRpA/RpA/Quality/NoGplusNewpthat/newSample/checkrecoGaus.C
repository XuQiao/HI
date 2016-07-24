#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"

static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};

const double pPbLumi = 15.78 ;

void checkrecoGaus(){

    TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.12);
//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

const int N=8;
const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12",""};
const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};

TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,1.2,1.6);
hFrame->SetTitle("");
//hFrame->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
hFrame->GetXaxis()->SetTitle("p_{T}^{reco} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetXaxis()->SetLimits(10,100);
hFrame->GetYaxis()->SetRangeUser(0,0.1);
hFrame->DrawCopy();
TLegend* leg1=new TLegend(0.7,0.65,0.88,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.03);
TLatex *T1=new TLatex();
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->DrawLatex(0.2,0.88,"PYTHIA+HIJING");
T1->DrawLatex(0.2,0.84,"CMS simulation");
T1->DrawLatex(0.6,0.45,"47<=p_{T}^{ref}<55 (GeV/c)");

TH2F* histo2F[N]; TH1F* histo[N];	
for(int i=N-1;i<N;i++){
if(i==N-1)
histo2F[i]=(TH2F*)fPPb->Get(Form("refptpt_0-100%%"));
else
histo2F[i]=(TH2F*)fPPb->Get(Form("refptptEtaBin%s_Cen0-100%%",etabinname[i].Data()));

histo[i]=(TH1F*)histo2F[i]->ProjectionY(Form("histo_%d",i),histo2F[i]->GetXaxis()->FindBin(47),histo2F[i]->GetXaxis()->FindBin(55));
//histo[i]->Rebin(2);
histo[i]->Scale(1./histo[i]->Integral());
histo[i]->SetMarkerColor(colorCode[i]);
histo[i]->SetLineColor(colorCode[i]);
histo[i]->SetMarkerSize(1.2);
histo[i]->SetMarkerStyle(markerCode[i]);
histo[i]->Draw("same");
double mean=histo[i]->GetMean();
TH1F* histocp=(TH1F*)histo[i]->Clone("histocp");
for(int ibin=1;ibin<histo[i]->GetNbinsX();ibin++){
histocp->SetBinContent(104-ibin,histo[i]->GetBinContent(ibin));}

leg1->AddEntry(histo[i],etastring[i],"lp");
cout<<mean<<endl;

cout<<100*(histo[i]->Integral(0,30)/histo[i]->Integral(0,100))<<"%"<<endl;
cout<<100*(histocp->Integral(0,30)/histocp->Integral(0,100))<<"%"<<endl;
histo[i]->Fit("gaus","","same",mean,100);
histocp->SetMarkerSize(1.2);
histocp->SetMarkerColor(2);
histocp->SetMarkerStyle(24);
histocp->Draw("same");
//gaus->Draw("same");
}
leg1->Draw("same");
T1->Draw("same");
gPad->RedrawAxis();
TString outname="checkrecodis";
c1->Print(Form("%sPPb.png",outname.Data()));
c1->Print(Form("%sPPb.pdf",outname.Data()));

}
