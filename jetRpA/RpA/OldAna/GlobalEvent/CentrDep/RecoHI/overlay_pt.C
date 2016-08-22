#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void overlay_pt(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

TString filename1="NoGplusCombinedJetTrigak3PFJetAllHist.root";
TString filename2="HFCorrNoGplusCombinedJetTrigak3PFJetAllHist.root";
TString filename3="NtrkCorrNoGplusCombinedJetTrigak3PFJetAllHist.root";
TString dir1="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/NoGplusNewpthat";
TString dir2="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/GlobalEvent/CentrDep";
TString dir3="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/GlobalEvent/CentrDep";

TString histoname="jetpt_0-100%";
TString histoname1="jetptEta_0-100%";
//TString histoname1=Form("jetpt%s_0-100%%",JetIDName.Data());

TFile *file1=TFile::Open(Form("%s/%s",dir1.Data(),filename1.Data()));
TFile *file2=TFile::Open(Form("%s/%s",dir2.Data(),filename2.Data()));
TFile *file3=TFile::Open(Form("%s/%s",dir3.Data(),filename3.Data()));

TH2F* histo12F=(TH2F*)file1->Get(histoname1);
TH2F* histo22F=(TH2F*)file2->Get(histoname1);
TH2F* histo32F=(TH2F*)file3->Get(histoname1);
//TH1D* histo1=(TH1D*)histo12F->ProjectionY("histo1",histo12F->GetXaxis()->FindBin(binbound_pt[5]),histo12F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
//TH1D* histo2=(TH1D*)histo22F->ProjectionY("histo2",histo22F->GetXaxis()->FindBin(binbound_pt[5]),histo22F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
TH1D* histo1=(TH1D*)histo12F->ProjectionX("histo1");
TH1D* histo2=(TH1D*)histo22F->ProjectionX("histo2");
TH1D* histo3=(TH1D*)histo32F->ProjectionX("histo3");

TH1D* histo1=(TH1D*)file1->Get(histoname);
TH1D* histo2=(TH1D*)file2->Get(histoname);
TH1D* histo3=(TH1D*)file3->Get(histoname);

rehisto1=(TH1D*)histo1->Clone("rehisto1");
rehisto2=(TH1D*)histo2->Clone("rehisto2");
rehisto3=(TH1D*)histo3->Clone("rehisto3");
rehisto1=(TH1D*)rehisto1->Rebin(Nbin_pt,"rehisto1",binbound_pt);
rehisto2=(TH1D*)rehisto2->Rebin(Nbin_pt,"rehisto2",binbound_pt);
rehisto3=(TH1D*)rehisto3->Rebin(Nbin_pt,"rehisto3",binbound_pt);
normalizeByBinWidth(rehisto1);
normalizeByBinWidth(rehisto2);
normalizeByBinWidth(rehisto3);
//rehisto1->Scale(1/rehisto1->Integral());
//rehisto2->Scale(1/rehisto2->Integral());
rehisto1->SetMarkerStyle(24);
rehisto1->SetMarkerSize(1.2);
rehisto1->SetMarkerColor(1);
rehisto1->SetLineColor(1);
rehisto2->SetMarkerStyle(20);
rehisto2->SetMarkerSize(1.2);
rehisto2->SetMarkerColor(2);
rehisto2->SetLineColor(2);
rehisto3->SetMarkerStyle(28);
rehisto3->SetMarkerSize(1.2);
rehisto3->SetMarkerColor(4);
rehisto3->SetLineColor(4);

c1 = new TCanvas("c1"," ",600,800);
makeMultiPanelCanvas(c1,1,2,-0.05,0,0.12,0.14,0.03);
TH1F* hFrame=new TH1F("hFrame","",10000,0,1000);
hFrame->GetYaxis()->SetTitle("1/N_{evt}dN/dp_{T}");
fixedFontHist(hFrame,2.2,2);

c1->cd(1)->SetLogy();
hFrame->GetXaxis()->SetRangeUser(30,600);
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetRangeUser(1.01e-11,2e-3);
hFrame->DrawCopy();
rehisto1->Draw("E1same");
rehisto2->Draw("E1same");
rehisto3->Draw("E1same");
TLegend *leg1=new TLegend(0.60,0.70,0.85,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.06);
leg1->AddEntry(rehisto1,"AN-12-377","lp");
leg1->AddEntry(rehisto2,"Use trackMult","lp");
leg1->AddEntry(rehisto3,"Use HF #eta>4","lp");
leg1->Draw("same");
TLatex *T1=new TLatex(0.5,0.90,"|#eta| < 1");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.06);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
c1->cd(2);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Ratio");
hFrame->GetYaxis()->SetRangeUser(0.8,1.19);
hFrame->DrawCopy();
TH1D* ratio1=(TH1D*)rehisto2->Clone("ratio1");
TH1D* ratio2=(TH1D*)rehisto3->Clone("ratio2");
ratio1->Divide(rehisto1);
ratio2->Divide(rehisto1);
TLegend *leg2=new TLegend(0.55,0.75,0.70,0.92);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetTextSize(0.05);
leg2->AddEntry(ratio1,"TrackMult/AN-12-377","lp");
leg2->AddEntry(ratio2,"HF #eta>4/AN-12-377","lp");
ratio1->Draw("E1same");
ratio2->Draw("E1same");
leg2->Draw("same");
TLine *l =new TLine(30,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
c1->Print("Comparept.png");
}
