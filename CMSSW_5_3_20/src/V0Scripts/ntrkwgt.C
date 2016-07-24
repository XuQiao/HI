{
#include <tdrstyle.C>
#include <sstream>
  //#include <iostream>
#include <string>

  gROOT->Reset();
  setTDRStyle();
  //gStyle->SetLabelOffset(0.017, "X");
  gStyle->SetStatFormat("5.5g");
  gStyle->SetStatY(0.94);
  gStyle->SetOptFit(111);
  gStyle->SetFitFormat("5.5g");
  gStyle->SetStatFontSize(0.03);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.3);
  gStyle->SetTextSize(0.038);
  gStyle->SetTextFont(42);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gStyle->SetStatFormat("NELU",AutoPrecision(3));

  double ntrkarray[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,84,88,92,96,100,110,120,140,160,200};
  int ntrkNbins = sizeof(ntrkarray)/sizeof(double) -1;
  vector<double> ntrkBinMinima;
  for ( int ndx = 0; ndx <= ntrkNbins; ndx++ ) {
    ntrkBinMinima.push_back( ntrkarray[ndx] );
  }
  Double_t ntrk_default_binWidth = 0.1;

  TString ntrk_data = "ntrk.root";
  //TString ntrk_mc = "ntrk_MC_pythia8.root";
  TString ntrk_mc = "ntrk_MC_tuneD6T.root";

  TFile* infile = TFile::Open(ntrk_data, "READ");
  TCanvas *c1 = (TCanvas*)infile.Get("c1");
  TH1F* ntrk_rec_data = (TH1F*)c1->GetPrimitive("ntrk_scale");
  double datarecentries = ntrk_rec_data->GetSumOfWeights();
  infile->Close();
  delete c1;
  delete infile;

  TFile* infile = TFile::Open(ntrk_mc, "READ");
  TCanvas *c1 = (TCanvas*)infile.Get("c1");
  TH1F* ntrk_rec_mc = (TH1F*)c1->GetPrimitive("ntrk_scale");
  double mcrecentries = ntrk_rec_mc->GetSumOfWeights();
  infile->Close();
  delete c1;
  delete infile;

  double scaleFactor = datarecentries/mcrecentries;

  TH1F *ntrk_weight = (TH1F*)ntrk_rec_data->Clone("ntrk_weight");
  TH1F *scale_ntrk = (TH1F*)ntrk_rec_data->Clone("scale_ntrk");
  for( int ndx = 0; ndx <= ntrkNbins+1; ndx++ ) {
    scale_ntrk->SetBinContent( ndx, scaleFactor );
    scale_ntrk->SetBinError( ndx, 0.0 );
  }
  ntrk_rec_mc->Multiply(scale_ntrk);
  ntrk_weight->Divide(ntrk_rec_data,ntrk_rec_mc);

  cout << "data reconstructed entries = " << datarecentries << endl;
  cout << "mc reconstructed entries = " << mcrecentries << endl;


  TCanvas* c1 = new TCanvas("c1", "Track multiplicity");

  ntrk_weight->SetAxisRange(0.4, 6.0, "Y");
  ntrk_weight->GetXaxis()->SetNdivisions(508);
  ntrk_weight->SetLineColor(1);
  ntrk_weight->SetLineWidth(1);
  ntrk_weight->GetXaxis()->SetTitle("Track multiplicity");
  ntrk_weight->GetYaxis()->SetTitle("Weight (data/MC)");
  ntrk_weight->SetMarkerSize(0.8);
  ntrk_weight->Draw();

  //TLatex *cmsprelim = new TLatex();
  //cmsprelim->SetNDC();
  //cmsprelim->SetTextSize(txtsize1);
  //cmsprelim->SetTextAlign(22);
  //cmsprelim->DrawLatex(0.40,0.91,"CMS Preliminary");
  //TLatex *cmsdata = new TLatex();
  //cmsdata->SetNDC();
  //cmsdata->SetTextSize(txtsize1);
  //cmsdata->SetTextAlign(22);
  //cmsdata->DrawLatex(0.40,0.87,"#sqrt{s} = 7 TeV");
  //
  //TLegend *leg = new TLegend(0.48,0.41,0.73,0.50);
  //leg->SetTextSize(0.040);
  //leg->SetFillColor(10);
  //leg->SetBorderSize(1); // no frame, no shadow
  //leg->AddEntry(ntrk_gen_mc, "Generated" , "p" );
  //leg->AddEntry(ntrk_cor_mc, "Corrected" , "p" );
  //leg->Draw("same");

  c1->SaveAs("ntrk_weight.png"); 
  c1->SaveAs("ntrk_weight.eps"); 
  c1->SaveAs("ntrk_weight.C");

}
