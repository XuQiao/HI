#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
void Compare(){
    const Double_t jetPtBin[]={3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};

    const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;
    bool isMC=1;
    if(isMC){
    TFile *f1 = new TFile("/cms/store/user/qixu/jetRpA/RpA/NewMC/MCPPbakPu3PF_useskim.root");
    TFile *f2 = new TFile("/home/maoy/Public/PPb_UnfoPriorGen_akPu3PFOfficialMCNoIDCut_MCJECv8_jtpt20_EtaBin-10_10_Inc_v6.root");
    TFile *f3 = new TFile("/cms/store/user/qixu/jetRpA/RpA/fromLinux/TreeAna/MCPPbakPu3PFskimUniv.root");
    TH1D* h1 = (TH1D*)f1->Get("jetpt");
    TH1D* h2 = (TH1D*)f2->Get("hMeas0");
    TH1D* h3 = (TH1D*)f3->Get("jetpt");
    }
    else{
    TFile *f1 = new TFile("/cms/store/user/qixu/jetRpA/RpA/NewMC/DATAPPbakPu3PFskimUniv.root");
    TFile *f2 = new TFile("/home/maoy/Public/DATAPPbakPu3PFJetSpectraKurtCombineJetTriggerEtaWeight7EtabinJetIDCutRecoPt.root");
    TFile *f3 = new TFile("/cms/store/user/qixu/jetRpA/RpA/fromLinux/TreeAna/DATAPPbakPu3PFskimUniv.root");
    TH1D* h1 = (TH1D*)f1->Get("jetpt");
    TH1D* h2 = (TH1D*)f2->Get("jetptEtaBin-10_10");
    TH1D* h3 = (TH1D*)f3->Get("jetpt");
   }
    TCanvas *c1 = new TCanvas("c1","c1",400,780);
    c1->Divide(1,2);
    c1->cd(1)->SetLogy();
    TH1D* h1_ = (TH1D*)h1->Rebin(nJetPtBin,"h1_",jetPtBin);
    TH1D* h2_ = (TH1D*)h2->Rebin(nJetPtBin,"h2_",jetPtBin);
    TH1D* h3_ = (TH1D*)h3->Rebin(nJetPtBin,"h3_",jetPtBin);
    normalizeByBinWidth(h1_);
    normalizeByBinWidth(h2_);
    normalizeByBinWidth(h3_);
    h1_->SetMarkerColor(2);
    h1_->SetLineColor(2);
    h1_->SetMarkerSize(1.2);
    h1_->SetMarkerStyle(20);
    h2_->SetMarkerColor(1);
    h2_->SetMarkerSize(1.2);
    h2_->SetLineColor(1);
    h2_->SetMarkerStyle(24);
    h3_->SetMarkerColor(1);
    h3_->SetMarkerSize(1.2);
    h3_->SetLineColor(1);
    h3_->SetMarkerStyle(24);

    h1_->Draw("P");
   // h2_->Draw("Psame");
    h3_->Draw("Psame");
    c1->cd(2);
   // TH1D* hratio = (TH1D*)h2_->Clone();
    TH1D* hratio = (TH1D*)h3_->Clone();
    hratio->Divide(h1_);
    hratio->GetYaxis()->SetRangeUser(0.8,1.2);
    hratio->Draw();
    TLine *l = new TLine(h2_->GetXaxis()->GetXmin(),1,h2_->GetXaxis()->GetXmax(),1);
    l->SetLineStyle(2);
    l->Draw("same");
}
    
