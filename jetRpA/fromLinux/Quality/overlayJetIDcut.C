#include "root_setting.h"
#include "file.h"

void overlayJetIDcut(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
double cut;
if(ilist== 6){      JetID = "chargedMax/p_{T}";	cut=0.3;}
else if(ilist== 7){      JetID = "chargedSum/p_{T}";	cut=0.6;}
else if(ilist== 8){     JetID = "neutralMax/p_{T}";	cut=0.08;}
else if(ilist== 9){     JetID = "neutralSum/p_{T}";	cut=0.15;}
//else if(ilist== 10){     JetID = "photonMax/p_{T}";	}
else if(ilist== 11){     JetID = "photonSum/p_{T}";	cut=0.3;}
else if(ilist== 12){     JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{jet}";cut =1.01;}
else if(ilist== 13){     JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{raw}";cut =1.13;}
else{   exit();}

if(ilist==12 || ilist==13 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}


double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};


int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;


TString filename1="NoResidual/newsubmitted/CombinedJetTrigakPu3PFJetAllHist.root";

TString filename2=Form("NoResidual/newsubmitted/%sCombinedJetTrigakPu3PFJetAllHist.root",JetIDName.Data());

TString dir1="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID";
TString dir2="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID";
TString histoname="jetpt_0-100%";
TString histoname1="jetptEta_0-100%";
//TString histoname1=Form("jetpt%s_0-100%%",JetIDName.Data());

TFile *file1=TFile::Open(Form("%s/%s",dir1.Data(),filename1.Data()));
TFile *file2=TFile::Open(Form("%s/%s",dir2.Data(),filename2.Data()));

TH2F* histo12F=(TH2F*)file1->Get(histoname1);
TH2F* histo22F=(TH2F*)file2->Get(histoname1);
//TH1D* histo1=(TH1D*)histo12F->ProjectionY("histo1",histo12F->GetXaxis()->FindBin(binbound_pt[5]),histo12F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
//TH1D* histo2=(TH1D*)histo22F->ProjectionY("histo2",histo22F->GetXaxis()->FindBin(binbound_pt[5]),histo22F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
TH1D* histo1=(TH1D*)histo12F->ProjectionX("histo1");
TH1D* histo2=(TH1D*)histo22F->ProjectionX("histo2");

//TH1D* histo1=(TH1D*)file1->Get(histoname);
//TH1D* histo2=(TH1D*)file2->Get(histoname);

rehisto1=(TH1D*)histo1->Clone("rehisto1");
rehisto2=(TH1D*)histo2->Clone("rehisto2");
rehisto1=(TH1D*)rehisto1->Rebin(Nbin_pt,"rehisto1",binbound_pt);
rehisto2=(TH1D*)rehisto2->Rebin(Nbin_pt,"rehisto2",binbound_pt);
rehisto1->GetXaxis()->SetRangeUser(xrange_pt[0],xrange_pt[1]);
rehisto2->GetXaxis()->SetRangeUser(xrange_pt[0],xrange_pt[1]);
normalizeByBinWidth(rehisto1);
normalizeByBinWidth(rehisto2);
//rehisto1->Scale(1/rehisto1->Integral());
//rehisto2->Scale(1/rehisto2->Integral());
rehisto1->SetMarkerStyle(24);
rehisto1->SetMarkerSize(1.4);
rehisto1->SetMarkerColor(1);
rehisto1->GetYaxis()->SetTitle("Number of Events");
rehisto2->SetMarkerStyle(20);
rehisto2->SetMarkerSize(1.2);
rehisto2->SetMarkerColor(2);
fixedFontHist(rehisto1,1.8,1.5);
fixedFontHist(rehisto2,1.8,1.5);

c1 = new TCanvas("c1","",600,600);
makeMultiPanelCanvas(c1,1,2,0.03,0.03,0.1,0.12,0.03);

c1->cd(1)->SetLogy();
rehisto1->SetTitle("");
rehisto1->Draw("E1");
rehisto2->Draw("E1same");
TLegend *leg1=new TLegend(0.60,0.75,0.85,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.06);
leg1->AddEntry(rehisto1,"Before cut","lp");
leg1->AddEntry(rehisto2,"After cut","lp");
leg1->Draw("same");
TLatex *T=new TLatex(0.35,0.30,Form("Cut: %s<%.2f",JetID.Data(),cut));
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.06);
T->SetTextColor(1);
T->SetTextFont(42);
c1->cd(2);
TH1D* ratio=(TH1D*)rehisto2->Clone("ratio");
ratio->Divide(rehisto1);
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.4);
ratio->SetMarkerColor(1);
ratio->SetMinimum(ratio->GetMinimum()/1.2);
ratio->SetMaximum(1.09);
ratio->SetTitle("");
ratio->GetYaxis()->SetTitle("Ratio: After cut/Before cut");
ratio->Draw("E1");
TLine *l =new TLine(binbound_pt[5],1,binbound_pt[Nbin_pt-1],1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
T->Draw("same");

c1->Print(Form("pic/JetIDcut/%s/overlay_JetIDcut.png",JetIDName.Data()));
}
