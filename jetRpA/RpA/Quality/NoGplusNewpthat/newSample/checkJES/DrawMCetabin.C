#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"
void DrawMCetabin(){
gStyle->SetOptStat(kFALSE);
gStyle->SetOptFit(kFALSE);
gStyle->SetErrorX(0);
static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};

const int N=8;
const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12",""};
//const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
//ofstream fstr("JESsysetabin.txt");
TH2F* hrawptJES[N];	 TProfile* hProrawptJES[N]; 
TH2F* hrefptJES[N];	 TProfile* hProrefptJES[N]; 
c1 = new TCanvas("c1","c1",600,600);
c2 = new TCanvas("c2","c2",600,600);
c2->cd(1)->SetLogx();
makeMultiPanelCanvas(c1,1,1,-0.12,0,0.14,0.12,0.03);
makeMultiPanelCanvas(c2,1,1,-0.12,0,0.14,0.12,0.03);
TLegend *leg1=new TLegend(0.72,0.55,0.88,0.9);
TLegend *leg2=new TLegend(0.72,0.55,0.88,0.9);
leg1->SetTextSize(0.035);
leg2->SetTextSize(0.035);
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(22.1,599.9);
hFrame->SetTitle("");
hFrame->GetYaxis()->SetRangeUser(0.95,1.3);
fixedFontHist(hFrame,1.2,1.6);
TLatex *T1=new TLatex;
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);
TLine *l=new TLine(31,1,599.9,1);
l->SetLineStyle(2);
l->SetLineColor(1);
for(int i=0;i<N;i++){
if(i==N-1){
hrawptJES[i]=(TH2F*)fPPb->Get("rawptJES_0-100%");
hrefptJES[i]=(TH2F*)fPPb->Get("refptJES_0-100%");
}
else{
hrawptJES[i]=(TH2F*)fPPb->Get(Form("rawptJESEtaBin%s_Cen0-100%%",etabinname[i].Data()));
hrefptJES[i]=(TH2F*)fPPb->Get(Form("refptJESEtaBin%s_Cen0-100%%",etabinname[i].Data()));
}
hrawptJES[i]->SetName(Form("hrawptJES_%d",i));
hrefptJES[i]->SetName(Form("hrefptJES_%d",i));
hProrawptJES[i]=(TProfile*)hrawptJES[i]->ProfileX();
hProrefptJES[i]=(TProfile*)hrefptJES[i]->ProfileX();
hProrawptJES[i]=(TProfile*)hProrawptJES[i]->Rebin(Nbin_pt,hProrawptJES[i]->GetName(),binbound_pt);
hProrefptJES[i]=(TProfile*)hProrefptJES[i]->Rebin(Nbin_pt,hProrefptJES[i]->GetName(),binbound_pt);
}

c1->cd(1);
hFrame->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
hFrame->GetYaxis()->SetTitle("p_{T}^{jet}/p_{T}^{raw}");
hFrame->DrawCopy();
for(int i=0;i<N;i++){
hProrawptJES[i]->SetMarkerStyle(markerCode[i]);
hProrawptJES[i]->SetMarkerColor(colorCode[i]);
hProrawptJES[i]->SetLineColor(colorCode[i]);
hProrawptJES[i]->SetMarkerSize(1.4);
leg1->AddEntry(hProrawptJES[i],etastring[i],"lp");
hProrawptJES[i]->DrawCopy("same");
}
for(int i=0;i<N-1;i++)	if(i!=3) hProrawptJES[i]->DrawCopy("same");

T1->DrawLatex(0.18,0.90,"PYTHIA+HIJING");
T1->DrawLatex(0.18,0.85,"CMS simulation");
T1->DrawLatex(0.18,0.80,"Anti-k_{T} Particle Flow Jets R=0.3");
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->Draw("same");
//drawCMS(0.4,0.9);
c1->Print(Form("OldrawptJES_etabin_%s.pdf",algo.Data()));
c2->cd(1)->SetLogx();
hFrame->GetXaxis()->SetNdivisions(1020);
TF1 *fitfun[N];
hFrame->GetYaxis()->SetRangeUser(0.98,1.1);
hFrame->GetXaxis()->SetRangeUser(29,600);
//TF1* fitfun=new TF1("fit","sqrt(pow([0],2)+pow([1]/sqrt(x),2)+pow([2]/x,2))",0,1000);
//fitfun->SetParameters(1,1,100,1000);
hFrame->GetXaxis()->SetTitle("p_{T}^{ref} (GeV/c)");
hFrame->GetYaxis()->SetTitle("p_{T}^{jet}/p_{T}^{ref}");
hFrame->DrawCopy();
//fstr<<"fitfun[i]=new TF1(Form(\"fit_\%d\",i),\"[0] + [1]/x + [2]/x/x + [3]/x/x/x\",1,1000);"<<endl;
for(int i=0;i<N;i++){
fitfun[i]=new TF1(Form("fit_%d",i),"[0] + [1]/x + [2]/x/x + [3]/x/x/x",1,1000);
fitfun[i]->SetParameters(1,-0.893019,-23.1071,4217.43);
fitfun[i]->FixParameter(0,1);
fitfun[i]->SetLineStyle(2);
fitfun[i]->SetLineColor(colorCode[i]);
hProrefptJES[i]->SetMarkerStyle(markerCode[i]);
hProrefptJES[i]->SetMarkerColor(colorCode[i]);
hProrefptJES[i]->SetLineColor(colorCode[i]);
hProrefptJES[i]->SetMarkerSize(1.4);
leg2->AddEntry(hProrefptJES[i],etastring[i],"lp");
/*hProrefptJES[i]->Fit(fitfun[i],"Q","same",20,600);
fstr<<etastring[i]<<endl;
fstr<<fitfun[i]->GetParameter(0)<<"\t";
fstr<<fitfun[i]->GetParameter(1)<<'\t';
fstr<<fitfun[i]->GetParameter(2)<<'\t';
fstr<<fitfun[i]->GetParameter(3)<<endl<<endl;
cout<<"chis="<<fitfun[i]->GetChisquare()<<endl;
cout<<"ndf="<<fitfun[i]->GetNDF()<<endl;*/
l->Draw("same");
hProrefptJES[i]->DrawCopy("same");
}
for(int i=0;i<N-1;i++) hProrefptJES[i]->DrawCopy("same");
//fitfun->Draw("same");
//leg2->AddEntry(fitfun,"fit function","l");
//for(int i=0;i<N-1;i++)  hProrefptJES[i]->DrawCopy("same");
T1->DrawLatex(0.18,0.90,"PYTHIA+HIJING");
T1->DrawLatex(0.18,0.85,"CMS simulation");
T1->DrawLatex(0.18,0.80,"Anti-k_{T} Particle Flow Jets R=0.3");
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->Draw("same");
c2->Print(Form("OldrefptJES_etabin_%s.pdf",algo.Data()));
//drawCMS(0.4,0.9);

}
