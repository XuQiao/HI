#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void CompareHong(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;


TString filename1="Datacombined.root";
TString filename2="AllEtaLumiCombinedSpectraAllEtaBin.root";

TString dir1="/scratch/xuq7/RpA/TreeAna";
TString dir2="/scratch/xuq7/RpA/TreeAna";

TString histoname1="jetptCombinedSpectraInEtaBin-3_3";
TString histoname2="CombinedSpectraInEtaBin-3_3";
//TString histoname1=Form("jetpt%s_0-100%%",JetIDName.Data());

TFile *file1=TFile::Open(Form("%s/%s",dir1.Data(),filename1.Data()));
TFile *file2=TFile::Open(Form("%s/%s",dir2.Data(),filename2.Data()));

TH1F* histo1=(TH1F*)file1->Get(histoname1);
TH1F* histo2=(TH1F*)file2->Get(histoname2);

histo1=(TH1F*)histo1->Rebin(Nbin_pt,"histo1",binbound_pt);
histo2=(TH1F*)histo2->Rebin(Nbin_pt,"histo2",binbound_pt);
normalizeByBinWidth(histo1);
normalizeByBinWidth(histo2);
//histo1->Scale(1/histo1->Integral());
//histo2->Scale(1/histo2->Integral());
histo1->SetMarkerStyle(24);
histo1->SetMarkerSize(1.2);
histo1->SetMarkerColor(1);
histo2->SetMarkerStyle(20);
histo2->SetMarkerSize(1.2);
histo2->SetMarkerColor(2);

c1 = new TCanvas("c1"," ",400,600);
makeMultiPanelCanvas(c1,1,2,-0.08,0,0.14,0.16,0.03);

TH1F* hFrame=new TH1F("","",20000,-1000,1000);
fixedFontHist(hFrame,2,2.2);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetTitle("Event Probability");
hFrame->GetXaxis()->SetLimits(0,100);
//hFrame->GetXaxis()->SetLimits(-3,3);
//hFrame->GetYaxis()->SetRangeUser(1.01e-10,1e-1);
hFrame->GetYaxis()->SetRangeUser(1,1e8);
c1->cd(1)->SetLogy();
hFrame->DrawCopy();
histo1->Draw("E1same");
histo2->Draw("E1same");
TLegend *leg1=new TLegend(0.60,0.75,0.85,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.06);
leg1->AddEntry(histo1,"Qiao","lp");
leg1->AddEntry(histo2,"Hong","lp");
leg1->Draw("same");
c1->cd(2);
TLatex *T1=new TLatex(0.5,0.92,"|#eta|<1  Data");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
TH1D* ratio=(TH1D*)histo2->Clone("ratio");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Ratio: Hong/Qiao");
hFrame->GetYaxis()->SetRangeUser(0.6,1.39);
hFrame->DrawCopy();
ratio->Divide(histo1);
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.2);
ratio->SetMarkerColor(1);
ratio->SetTitle("");
ratio->Draw("E1same");
TLine *l =new TLine(0,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
T1->Draw("same");
}
