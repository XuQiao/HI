#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"
void CompareOldNew(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

TString histoname="jetptEtaBin-10_10";
TString histoname1="jetpt";
TString histoname2F="jetptEta";
TFile *fMCPPbYX = new TFile("/home/xuq/Documents/HI/RpA/output/TreeAna/MCPPbakPu3PFJetSpectraCombineTriggerJetPtNoJetIDCutNoResidualTrkCorrFile0_8.root");
TFile *fDataPPbYX = new TFile("/home/xuq/Documents/HI/RpA/output/TreeAna/DATAPPbakPu3PFJetSpectraCombineTriggerJetPtNoJetIDCutNoResidualTrkCorrFile0_0.root");

TH2F* histo12F=(TH2F*)fMCPPb->Get(histoname2F);
TH2F* histo22F=(TH2F*)fMCPPbYX->Get(histoname2F);
//TH1D* histo1=(TH1D*)histo12F->ProjectionY("histo1",histo12F->GetXaxis()->FindBin(binbound_pt[5]),histo12F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
//TH1D* histo2=(TH1D*)histo22F->ProjectionY("histo2",histo22F->GetXaxis()->FindBin(binbound_pt[5]),histo22F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
//TH1D* histo1=(TH1D*)histo12F->ProjectionX("histo1");
//TH1D* histo2=(TH1D*)histo22F->ProjectionX("histo2");

TH1D* histo1=(TH1D*)fDataPPb->Get(histoname1);
TH1D* histo2=(TH1D*)fDataPPbYX->Get(histoname);

histo1=(TH1D*)histo1->Rebin(Nbin_pt,"histo1",binbound_pt);
histo2=(TH1D*)histo2->Rebin(Nbin_pt,"histo2",binbound_pt);
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

c1 = new TCanvas("c1"," ",500,800);
makeMultiPanelCanvas(c1,1,2,-0.14,0,0.2,0.14,0.03);
TLatex *T1 = new TLatex();
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
TH1F* hFrame=new TH1F("","",20000,-1000,1000);
fixedFontHist(hFrame,2,3.1);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{T}d#eta}");
hFrame->GetXaxis()->SetLimits(0,1000);
//hFrame->GetYaxis()->SetRangeUser(1.01e-11,1e-1);
hFrame->GetYaxis()->SetRangeUser(1.01e-1,1e6);
c1->cd(1)->SetLogy();
hFrame->DrawCopy();
histo1->Draw("E1same");
histo2->Draw("E1same");
TLegend *leg1=new TLegend(0.70,0.75,0.9,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.06);
leg1->AddEntry(histo1,"New","lp");
leg1->AddEntry(histo2,"Old","lp");
leg1->Draw("same");
T1->Draw("same");
for(int ibin=0;ibin<histo1->GetNbinsX();ibin++)
if(histo1->GetBinContent(ibin)!=0)
cout<<(histo1->GetBinContent(ibin)-histo2->GetBinContent(ibin))/histo1->GetBinContent(ibin)<<endl;
c1->cd(2);
TLatex *T1=new TLatex(0.5,0.92,"|#eta|<1");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
TH1D* ratio=(TH1D*)histo2->Clone("ratio");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Ratio: New/Old");
hFrame->GetYaxis()->SetRangeUser(0.6,1.39);
hFrame->DrawCopy();
ratio->Divide(histo1);
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.2);
ratio->SetMarkerColor(1);
ratio->SetTitle("");
ratio->Draw("E1same");
TLine *l =new TLine(30,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
T1->Draw("same");
}

