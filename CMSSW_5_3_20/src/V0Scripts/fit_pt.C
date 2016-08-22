{
#include <tdrstyle.C>
#include "tsallis.C"


bool isMC = false;
bool kshort = true;


  // Style and other ROOT commands
setTDRStyle();
gStyle->SetOptStat(1);
TGaxis::SetMaxDigits(3);
//gROOT->SetBatch(kTRUE);

TF1 *tsfunc = new TF1("tsfunc",tsallis,0.0,10.0,3);
tsfunc->SetParNames("Amplitude","T","n");
tsfunc->SetParameters(1000000.0.,0.2,7.2); 

TString inFileName;
TString outFileName;

    outFileName = "ks_pythia8/plots/kspt_cor_data_pythia8_final.root";
    inFileName = "ks_pythia8/plots/kspt_cor_data_pythia8.root";
    TFile* fin = TFile::Open(inFileName, "READ");
    TCanvas *c1 = (TCanvas*)fin.Get("c1");
    TH1F* v0_dNdpT_cor_data = (TH1F*)c1->GetPrimitive("ks_dNdpT_scale");
    fin->Close();
    delete c1;
    delete fin;
//if( !isMC ) {
//  if ( kshort ) {
//  }
//  else {
//    inFileName = "lam_pythia8/plots/lampt_cor_data_pythia8.root";
//    outFileName = "lam_pythia8/plots/lampt_cor_data_pythia8_final.root";
//  }
//}
//else {
//  if ( kshort ) {
//    inFileName = "ks_pythia8/plots/kspt_cor_mc_pythia8.root";
//    outFileName = "ks_pythia8/plots/kspt_cor_mc_pythia8_final.root";
//  }
//  else {
//    inFileName = "lam_pythia8/plots/lampt_cor_mc_pythia8.root";
//    outFileName = "lam_pythia8/plots/lampt_cor_mc_pythia8_final.root";
//  }
//}
//double dataentries = v0_dNdpT_cor_data->Integral();

TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

v0_dNdpT_cor_data->Fit("tsfunc","IEM0");
tsfunc->SetNpx(50000);
v0_dNdpT_cor_data->Draw();
tsfunc->Draw("same"); 

c1->SaveAs(outFileName);

}
