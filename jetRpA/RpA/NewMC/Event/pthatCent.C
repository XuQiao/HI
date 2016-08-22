#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
#include "/home/xuq7/HI/jetRpA/RpA/NewMC/produceandcheck/file.h"
#include <iomanip>
void pthatCent(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TString coll="PPb";
TString MCOPPbfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPbakPu3PFskimUniv.root";
TFile *fOMCPPb=TFile::Open(MCOPPbfile);
TString MCPPbfilepthat="/scratch/xuq7/RpA/NewMC/MCPPbakPu3PF_useskimMLast2pthat.root";
TFile *fMCPPbpthat = TFile::Open(MCPPbfilepthat);
bool Save=kTRUE;
TF1 * fCen = new TF1("fCen","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x+[5]*x*x*x*x+[6]*x*x*x*x*x)", 0., 100.);
if(coll=="PPb"){
//fCen->SetParameters(1.20916e-02, 5.02157e+00, -3.38300e-02, 1.87647e-03, -6.76442e-05, 9.08602e-07, -4.01536e-09);//parameterize on 05.03 after approval
fCen->SetParameters(7.92204e-03, 4.52005e+00, 9.77340e-02, -5.00362e-03, 9.74735e-05, -8.93897e-07, 3.39375e-09);//parameterize on new official MC after QM on 12/03/14
TH1F* histo0=(TH1F*)fMCPPbpthat->Get("Cent");
TH1F* histo1=(TH1F*)fOMCPPb->Get("Cent");
//TH1F* histo2=(TH1F*)fNMCPPb->Get("Cent");
TH1F* histo2=(TH1F*)fMCPPb->Get("Cent");
TString data="Pb going positive side";
}
//fCen->SetParameters(8.68073e-03, 5.09356e+00, -1.33053e-02, 1.46904e-03, -6.99681e-05, 1.06721e-06, -5.21398e-09); //original
else if(coll=="PbP"){
//fCen->SetParameters(1.05408e-02, 5.27477e+00, -8.03382e-02, 3.51669e-03, -8.85332e-05, 1.08917e-06, -4.90091e-09);
//fCen->SetParameters(1.14851e-02, 5.31172e+00, -8.52366e-02, 3.00268e-03, -6.04667e-05, 6.24105e-07, -2.43580e-09);
fCen->SetParameters(2.89263e-02, 3.43643e+00, 5.62562e-02, -1.86555e-03, 1.97924e-06, 3.15416e-07, -1.97946e-09);//parameterize on new official MC after QM on 12/03/14
TH1F* histo0=(TH1F*)fMCPbPpthat->Get("Cent");
TH1F* histo1=(TH1F*)fOMCPb->Get("Cent");
//TH1F* histo2=(TH1F*)fNMCPbP->Get("Cent");
TH1F* histo2=(TH1F*)fMCPbP->Get("Cent");
TString data="Proton going positive side";
}
cout<<"pthat: "<<histo0->GetEntries()<<endl;
cout<<"Old: "<<histo1->GetEntries()<<endl;
cout<<"All: "<<histo2->GetEntries()<<endl;
histo1->SetName(Form("%sOldMC_unweighted",coll.Data()));
histo2->SetName(Form("%sNewMC_unweighted",coll.Data()));
histo0->SetName(Form("%spthat",coll.Data()));
//histo1=(TH1F*)histo1->Rebin(Nbin_vz,"histo1",binbound_vz);
//histo2=(TH1F*)histo2->Rebin(Nbin_vz,"histo2",binbound_vz);
histo1->Scale(1/histo1->Integral());
histo2->Scale(1/histo2->Integral());
histo0->Scale(1/histo0->Integral());
histo1->SetMarkerStyle(24);
histo1->SetMarkerSize(1.2);
histo1->SetMarkerColor(2);
histo1->SetLineColor(2);
histo2->SetMarkerStyle(29);
histo2->SetMarkerSize(1.2);
//histo2->SetFillStyle(3004);
//histo2->SetFillColor(2);
histo2->SetLineColor(4);
histo2->SetMarkerColor(4);
histo2->SetMarkerColor(4);
histo0->SetMarkerStyle(20);
histo0->SetMarkerSize(1.2);
histo0->SetMarkerColor(1);
histo0->SetLineColor(1);

TCanvas* c1 = new TCanvas("c1"," ",500,500);
TCanvas* c2 = new TCanvas("c2"," ",500,500);
makeMultiPanelCanvas(c1,1,1,-0.16,0,0.16,0.14,0.03);
makeMultiPanelCanvas(c2,1,1,-0.16,0,0.16,0.14,0.03);
TH1F* hFrame=new TH1F("","",20000,-1000,1000);
fixedFontHist(hFrame,1.2,1.7);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitleSize(0.05);
hFrame->GetYaxis()->SetTitleSize(0.05);
hFrame->GetXaxis()->SetTitle("Centrality");
hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetXaxis()->SetLimits(0,100);
//hFrame->GetXaxis()->SetLimits(-3,3);
hFrame->GetYaxis()->SetRangeUser(0,0.07);
c1->cd(1);
hFrame->DrawCopy();
histo0->Draw("same");
histo1->Draw("same");
histo2->Draw("same");
TLegend *leg1=new TLegend(0.20,0.80,0.85,0.92);
TLegend *leg2=new TLegend(0.20,0.80,0.85,0.92);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetTextSize(0.04);
leg1->AddEntry(histo1,"Old MC Before CW","lp");
leg1->AddEntry(histo0,"New MC Before CW Last 2 pthat","lp");
leg1->AddEntry(histo2,"New MC Before CW all","lp");
leg1->Draw("same");
TLatex *T1=new TLatex(0.25,0.92,"");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
c2->cd(1);
hFrame->GetYaxis()->SetTitle("Ratio");
hFrame->GetXaxis()->SetLimits(0,100);
hFrame->GetYaxis()->SetRangeUser(0,2);
hFrame->DrawCopy();
TH1F* ratio = (TH1F*)histo0->Clone(Form("%sratio",coll.Data()));
ratio->SetTitle("");
ratio->Divide(histo2);

ratio->Draw();
leg2->AddEntry(ratio,"Last 2 pthat/all","lp");
leg2->Draw("same");

c1->Print(Form("Centpthat%s.png",coll.Data()));

}
