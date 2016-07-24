{
#include <TH1F.h>
#include <TFile.h>
#include "tsallis.C"
  gROOT->SetBatch(0);

  TFile* fin = TFile::Open("plots/lampt_cor_data_pythia8.root", "READ");

  TCanvas* c1 = (TCanvas*) fin->Get("c1");
  TH1F* lam_dNdpT_in = (TH1F*) c1->GetPrimitive("lam_dNdpT_scale");

  double lamptarray[] = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0,7.0,8.0,10.0};
  int lamptNbins = sizeof(lamptarray)/sizeof(double) -1;
  vector<Double_t> lamptBinMinima;
  for ( int ndx = 0; ndx <= lamptNbins; ndx++ ) {
    lamptBinMinima.push_back( lamptarray[ndx] );
  }
  Double_t lampt_default_binWidth = 0.1;

  Double_t lamptOrigNbins = lamptNbins - 3;
  
  //TH1F* lam_ptscale = new TH1F("ks_ptscale", "pT scaling factor", lamptNbins, *lamptBinMinima[0]);
  TH1F* lam_dNdpT = new TH1F("ks_dNdpT_scale", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);

  for( int ndx = 1; ndx <= lamptNbins; ndx++ ) {
    if( ndx == 1 || ndx == 2 || ndx == 3) continue;
    int origNdx = ndx - 3;
    lam_dNdpT->SetBinContent(ndx, lam_dNdpT_in->GetBinContent(origNdx));
  }

  TF1* tsfunc = new TF1("tsfunc", tsallis, 0.3, 10., 3);
  TF1* tsfunc2 = new TF1("tsfunc2", tsallis, 0., 10., 3);

  tsfunc->SetParNames("Amplitude", "T", "n");
  tsfunc->SetParameters(1000000., 0.2, 7.2);
  tsfunc2->SetParNames("Amplitude","T", "n");
  tsfunc2->SetParameters(1000000., 0.2, 7.2);

  lam_dNdpT_in->SetMarkerColor(kRed);

  TCanvas* can = new TCanvas("c1","V0 pT distribution");

  lam_dNdpT_in->Fit("tsfunc", "IEMO");
  tsfunc->SetNpx(50000);
  tsfunc2->SetParameters( tsfunc->GetParameters() );
  tsfunc2->SetNpx(50000);

  lam_dNdpT->Draw();
  lam_dNdpT->SetMarkerStyle(1);
  tsfunc2->Draw("same");
  can->SaveAs("test.png");
}
