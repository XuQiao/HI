#include "../../Quality/root_setting.h"
#include "../produceandcheck/file.h"
#include <iomanip>
#include <iostream>
#include "extract.C"

ofstream fstr[Neta];
int ilist0=21, ilist1=6, ilist2=14;

TH1D* makehisto(int ilist, int i,double cut, TString coll){
TString JetIDName = JetIDNameList[ilist];
if(ilist== 6)  {    JetID = "Max p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0]=0.05; JetIDcut[1]=1; }
else if(ilist== 7)  {    JetID = "#Sum p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "Max p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "#Sum p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10)  {   JetID = "Max p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "#Sum p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{jet}}";JetIDcut[0]=0; JetIDcut[1]=cut;}
else if(ilist== 13) {    JetID = "#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}/p_{T}^{raw}"; JetIDcut[0]=0;JetIDcut[1]=1.13;}
else if(ilist== 14)   {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})";JetIDcut[0]=0;JetIDcut[1]=0.975;}
else if(ilist== 21)   {  JetID = "PP cut Tighter";JetIDcut[0]=0;JetIDcut[1]=0;}

if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
TH2F* h2F;
if(coll.Contains("Data")){
	TFile *fPPbJetID = fDataPPbJetID;
	TFile *fPbPJetID = fDataPbPJetID;
	TString jetq = "";
}
else{
	TFile *fPPbJetID = fMCPPb;
	TFile *fPbPJetID = fMCPbP;
	if(coll.Contains("fake"))
		TString jetq = "_fake1";
	else if(coll.Contains("real"))
		TString jetq = "_real1";
	else
		TString jetq = "";
}

if(i==7){
    TString histonameID=Form("jetpt%s%s",JetIDName.Data(),jetq.Data());
}
else{
    TString histonameID=Form("jetpt%sEtaBin%s%s",JetIDName.Data(),etabinnameswap[i].Data(),jetq.Data());
}
    if(coll.Contains("PPb")){
    h2F=(TH2F*)fPPbJetID->Get(histonameID);
    }
    else{
    h2F=(TH2F*)fPbPJetID->Get(histonameID);
    }
xrange_JetIDcut[0]=JetIDcut[0]+1e-4;
xrange_JetIDcut[1]=JetIDcut[1]-1e-4;
TH1D* histo=new TH1D("","",1000,0.,1000.);
histo->Sumw2();
if(ilist==21) {
	for(int ibin=1;ibin<h2F->GetYaxis()->GetNbins();ibin++){
		double binxcenter = h2F->GetYaxis()->GetBinCenter(ibin);
		TVectorD* vecx2 = (TVectorD*)extract(binxcenter);
		TH1D* histo_=(TH1D*)h2F->ProjectionX(Form("histo_%d",ibin),ibin,ibin);
		if((*vecx2)[(int)cut]!=0)
		histo->Add(histo_);
	}
}
else histo=(TH1D*)h2F->ProjectionX(Form("histo"),h2F->GetYaxis()->FindBin(xrange_JetIDcut[0]),h2F->GetYaxis()->FindBin(xrange_JetIDcut[1]));
rehisto=(TH1D*)histo->Rebin(Nbin_pt,Form("rehisto"),binbound_pt);
normalizeByBinWidth(rehisto);
    for(int ibin = 1; ibin <= rehisto->GetNbinsX(); ibin++) {
        rehisto->SetBinError(ibin,1e-18);
    }
rehisto->SetMarkerSize(1.4);
h2F=NULL;
return rehisto;
}

void DrawsysJetIDcut_Etabin(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
//gStyle->SetOptFit(1);
//gStyle->SetOptTitle(0);
//c1 = new TCanvas("c1","",600,1000);
//makeMultiPanelCanvas(c1,1,2,0.03,0.03,0.1,0.12,0.03);
c1 = new TCanvas("c1"," ",1200,600);
makeMultiPanelCanvas(c1,4,2,0,0,0.25,0.2,0.03);
TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,2.0,3.0);
hFrame->SetTitle("");
//hFrame->GetYaxis()->SetTitle("Yield Ratio");
hFrame->GetXaxis()->SetLimits(27.5,600);
hFrame->GetYaxis()->SetRangeUser(0.98,1.01);
TLatex T;

for(int i=0;i<Neta;i++){
    c1->cd(canvas[i]+1)->SetGridx();
    c1->cd(canvas[i]+1)->SetGridy();
    if(canvas[i]==0  || canvas[i]==4)
        hFrame->GetYaxis()->SetTitle("Yield Ratio");
    else
        hFrame->GetYaxis()->SetTitle("");
    if(canvas[i]>=4)
        hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
    else
        hFrame->GetXaxis()->SetTitle("");

    hFrame->DrawCopy();
//fstr[i].open(Form("jetIDsys%s.txt",etabinname[i].Data()));
//fstr[i]<<setprecision(4)<<fixed;
//TH1D* histo0 = makehisto(ilist0,i,1.01);
TString coll="DataPPb";
TH1D* histo1PPb = makehisto(ilist0,i,0,coll);
TH1D* histo2PPb = makehisto(ilist0,i,1,coll);
TH1D* histo3PPb = makehisto(ilist0,i,2,coll);
TH1D* histo4PPb = makehisto(ilist0,i,3,coll);
TH1D* hratioPPb2 = (TH1D*)histo2PPb->Clone();
TH1D* hratioPPb3 = (TH1D*)histo3PPb->Clone();
TH1D* hratioPPb4 = (TH1D*)histo4PPb->Clone();
hratioPPb2->Divide(histo1PPb);
hratioPPb3->Divide(histo1PPb);
hratioPPb4->Divide(histo1PPb);
hratioPPb2->SetMarkerSize(1.2);
hratioPPb2->SetMarkerStyle(20);
hratioPPb2->SetMarkerColor(1);
hratioPPb2->SetLineColor(1);
hratioPPb3->SetMarkerSize(1.2);
hratioPPb3->SetMarkerStyle(24);
hratioPPb3->SetMarkerColor(2);
hratioPPb3->SetLineColor(2);
hratioPPb4->SetMarkerSize(1.2);
hratioPPb4->SetMarkerStyle(32);
hratioPPb4->SetMarkerColor(4);
hratioPPb4->SetLineColor(4);
/*
histo1PbP = makehisto(ilist1,i,0.99,"PbP");
histo2PbP = makehisto(ilist2,i,1.03,"PbP");
hratioPbP = (TH1D*)histo2PbP->Clone();
hratioPbP->Divide(histo1PbP);
hratioPbP->SetMarkerSize(1.2);
hratioPbP->SetMarkerStyle(27);
hratioPbP->SetMarkerColor(6);
hratioPbP->SetLineColor(6);
*/
if(canvas[i]!=4){   
hratioPPb2->Draw("same");
hratioPPb3->Draw("same");
hratioPPb4->Draw("same");
//hratioPbP->Draw("same");
}
/*for(int ibin=1;ibin<histo1PPb->GetNbinsX();ibin++){
if(histo1PPb->GetBinContent(ibin)!=0 && histo1PPb->GetBinCenter(ibin)>25 && histo1PPb->GetBinCenter(ibin)<=600){
fstr[i]<<histo1PPb->GetBinCenter(ibin)<<'\t';
fstr[i]<<100*(TMath::Abs(histo2PPb->GetBinContent(ibin)/histo1PPb->GetBinContent(ibin)-1)+TMath::Abs(histo2PbP->GetBinContent(ibin)/histo1PbP->GetBinContent(ibin)-1))/2<<endl;
}
}
*/
/*TString JetID0 = histo0->GetTitle();
JetID0=Form("JetID systematics");
JetID1.Prepend("Cut1: ");
JetID2.Prepend("Cut2: ");
histo0->SetMarkerStyle(20);
histo0->SetMarkerColor(1);
histo0->SetLineColor(1);
histo1->SetMarkerStyle(24);
histo1->SetMarkerColor(2);
histo1->SetLineColor(2);
histo2->SetMarkerStyle(29);
histo2->SetMarkerColor(4);
histo2->SetLineColor(4);
histo0->Draw("E1same");
histo1->Draw("E1same");
histo2->Draw("E1same");
TLegend *leg1=new TLegend(0.20,0.68,0.90,0.90);
TLegend *leg2=new TLegend(0.70,0.80,0.85,0.90);
leg1->SetHeader(JetID0);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetTextSize(0.04);
leg1->AddEntry(hFrame,"","");
leg1->AddEntry(histo0,"Original: Var < 1.01","p");
leg1->AddEntry(histo1,"Cut1: Var < 0.99","p");
leg1->AddEntry(histo2,"Cut2: Var < 1.03","p");
leg1->Draw("same");
TLatex *T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.04);
T->SetTextColor(1);
T->SetTextFont(42);
c1->cd(2);
TH1D* ratio1=(TH1D*)histo1->Clone("ratio1");
TH1D* ratio2=(TH1D*)histo2->Clone("ratio2");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle(Form("Ratio"));
hFrame->GetYaxis()->SetRangeUser(0.8,1.19);
hFrame->DrawCopy();
ratio1->Divide(histo0);
ratio2->Divide(histo0);
ratio1->Draw("E1same");
ratio2->Draw("E1same");
leg2->AddEntry(ratio1,"Cut1/Original","lp");
leg2->AddEntry(ratio2,"Cut2/Original","lp");
leg2->Draw("same");
*/
TLine *l =new TLine(30,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
if(canvas[i]!=4)
l->Draw("same");
T.SetNDC();
T.SetTextAlign(12);
T.SetTextSize(0.04);
T.SetTextColor(1);
T.SetTextFont(42);
T.SetTextSize(0.065);
if(canvas[i]<4)
T.DrawLatex(0.4,0.18,etastring[i]);
else if(canvas[i]>4)
T.DrawLatex(0.4,0.30,etastring[i]);
else{
    TLegend *leg1=new TLegend(0.28,0.25,0.85,0.45);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetTextSize(0.050);
    leg1->AddEntry(hratioPPb2,"nF<0.7,eF<0.8/nF<0.8,eF<0.8","p");
    leg1->AddEntry(hratioPPb3,"nF<0.6,eF<0.8/nF<0.8,eF<0.8","p");
    leg1->AddEntry(hratioPPb4,"nF<0.7,eF<0.7/nF<0.8,eF<0.8","p");
    leg1->Draw("same");
    T.SetTextSize(0.065);
    if(coll.Contains("PPb"))
    T.DrawLatex(0.35,0.72,"PbP(Pb going postive)");
    else
    T.DrawLatex(0.35,0.72,"PPb(Proton going postive)");
    if(coll.Contains("Data"))
    T.DrawLatex(0.35,0.80,"CMS Preliminary");
    else{
    T.DrawLatex(0.35,0.80,"CMS Simulation");
	if(coll.Contains("fake"))
    T.DrawLatex(0.35,0.6,"Fake jets");
	else if(coll.Contains("real"))
    T.DrawLatex(0.35,0.66,"Real jets");
	else
    T.DrawLatex(0.35,0.66,"Inclusive jets");
    }
}
}
c1->Print(Form("pic/JetIDcutsys_%s_Etabin.gif",coll.Data()));
c1->Print(Form("pic/JetIDcutsys_%s_Etabin.pdf",coll.Data()));

}
