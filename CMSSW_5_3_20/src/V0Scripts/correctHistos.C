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

  // algo determines which efficiency to use to correct the histogram.
  //  1: D6T
  //  2: pythia8
  //  3: P0
  int algo = 11;
  int isWeighted = 3;
  // be CAREFUL about this.  Are you weighting tracks and checking MC?  Don't do it.
  bool isMC = false;
  bool isKshort = true;

  int weight = 0;
  ostringstream weightoss;
  if( isWeighted ) {
    weight = algo + 100;
  }
  weightoss << weight;
  //if(!isKshort) weightoss << "_subgeant2";

  TString energy = "7TeV";
  TString energy2 = "7000";
if ( algo >= 10 ) {
  energy = "900GeV";
  energy2 = "900";
}
  TString algostr;
  if( algo == 1 || algo == 11 ) algostr = "tuneD6T";
  if( algo == 2 ) algostr = "pythia8";
  if( algo == 3 ) algostr = "tuneP0";

  //TString effSuffixAlgo1 = "tuneD6T_" + energy + "_weight" + weightoss.str().c_str() + ".root";
  //TString effSuffixAlgo2 = "pythia8_" + energy + "_weight" + weightoss.str().c_str() + ".root";
  //TString effSuffixAlgo3 = "tuneP0_" + energy + "_weight" + weightoss.str().c_str() + ".root";

  // Fit functions, constants
  double piMass = 0.13957018;
  double protonMass = 0.93827203;
  double ksMassConst = 0.497614;
  double lamMassConst = 1.115683;
  double cConst = 2.99792458e-2; // in cm/ps
  double ksPDGlife = 89.53;//ps
  double lamPDGlife = 263.1;//ps

  ostringstream ksctoss;
  ksctoss << "exp([0])*exp(-x/(" << cConst << "*[1]))";

  ostringstream lamctoss;
  lamctoss << "exp([0])*exp(-x/(" << cConst << "*[1]))";

  ksctFit = new TF1("ksctFit", ksctoss.str().c_str(), 0., 12.);
  ksctFit->SetParName(1, "t [ps]");

  lamctFit = new TF1("lamctFit", lamctoss.str().c_str(), 0., 20.);
  lamctFit->SetParName(1, "t [ps]");

  if( isKshort ) {

    if( isMC ) {
      TString kslife_rec = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdlife_MC_" + algostr + ".root";
      TString ksy_rec = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdy_MC_" + algostr + ".root";
      TString ksy_pos_rec = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdy_pos_MC_" + algostr + ".root";
      TString ksy_neg_rec = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdy_neg_MC_" + algostr + ".root";
      TString kspt_rec = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdpT_MC_" + algostr + ".root";
    }
    else {
      kslife_rec = "ks_data_" + energy2 + "/plots/ks_dNdlife.root";
      ksy_rec = "ks_data_" + energy2 + "/plots/ks_dNdy.root";
      ksy_pos_rec = "ks_data_" + energy2 + "/plots/ks_dNdy_pos.root";
      ksy_neg_rec = "ks_data_" + energy2 + "/plots/ks_dNdy_neg.root";
      kspt_rec = "ks_data_" + energy2 + "/plots/ks_dNdpT.root";
    }

    kslife_eff = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdlife_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + ".root";
    ksy_eff = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdy_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + ".root";
    ksy_pos_eff = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdy_pos_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + ".root";
    ksy_neg_eff = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdy_neg_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + ".root";
    kspt_eff = "ks_" + algostr + "_" + energy2 + "/plots/ks_dNdpT_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + ".root";

    TString kscorLtimeFname;
    TString kscoryFname;
    TString kscorpTFname;
    TString kscorLtimeFname_data;
    TString kscoryFname_data;
    TString kscorpTFname_data;
    kscorLtimeFname = "plots/kslife_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";
    kscoryFname = "plots/ksy_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";
    kscory_negFname = "plots/ksy_neg_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";
    kscory_posFname = "plots/ksy_pos_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";
    kscorpTFname = "plots/kspt_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";
    kscorLtimeFname_data = "plots/kslife_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";
    kscoryFname_data = "plots/ksy_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";
    kscory_negFname_data = "plots/ksy_neg_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";
    kscory_posFname_data = "plots/ksy_pos_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";
    kscorpTFname_data = "plots/kspt_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "";

    /*
      cout << kslife_rec << endl;
      cout << ksy_rec << endl;
      cout << kspt_rec << endl;
      cout << kslife_eff << endl;
      cout << ksy_eff << endl;
      cout << kspt_eff << endl;
    */

    // Get histos to correct

    TFile* fin = TFile::Open(kslife_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdlife_rec_data = (TH1F*)c1->GetPrimitive("ks_dNdlife_scale");
    TH1F* ks_dNdlife_cor_data = (TH1F*)c1->GetPrimitive("ks_dNdlife_scale");
    double reclifeentries = ks_dNdlife_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(ksy_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdy_rec_data = (TH1F*)c1->GetPrimitive("ks_dNdy");
    TH1F* ks_dNdy_cor_data = (TH1F*)c1->GetPrimitive("ks_dNdy");
    double recyentries = ks_dNdy_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(ksy_neg_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdy_neg_rec_data = (TH1F*)c1->GetPrimitive("ks_dNdy_neg");
    TH1F* ks_dNdy_neg_cor_data = (TH1F*)c1->GetPrimitive("ks_dNdy_neg");
    double recynegentries = ks_dNdy_neg_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(ksy_pos_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdy_pos_rec_data = (TH1F*)c1->GetPrimitive("ks_dNdy_pos");
    TH1F* ks_dNdy_pos_cor_data = (TH1F*)c1->GetPrimitive("ks_dNdy_pos");
    double recyposentries = ks_dNdy_pos_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(kspt_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdpT_rec_data = (TH1F*)c1->GetPrimitive("ks_dNdpT_scale");
    TH1F* ks_dNdpT_cor_data = (TH1F*)c1->GetPrimitive("ks_dNdpT_scale");
    double recptentries = ks_dNdpT_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    // Get efficiency histos

    TFile* fin = TFile::Open(kslife_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdlife_eff_mc = (TH1F*)c1->GetPrimitive("ks_dNdlife_eff");
    double mcefflifeentries = ks_dNdlife_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(ksy_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdy_eff_mc = (TH1F*)c1->GetPrimitive("ks_dNdy_eff");
    double mceffyentries = ks_dNdy_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(ksy_pos_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdy_pos_eff_mc = (TH1F*)c1->GetPrimitive("ks_dNdy_pos_eff");
    double mceffy_posentries = ks_dNdy_pos_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(ksy_neg_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdy_neg_eff_mc = (TH1F*)c1->GetPrimitive("ks_dNdy_neg_eff");
    double mceffy_negentries = ks_dNdy_neg_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(kspt_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* ks_dNdpT_eff_mc = (TH1F*)c1->GetPrimitive("ks_dNdpT_eff");
    double mceffptentries = ks_dNdpT_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    cout << "reconstructed lifetime entries: " << reclifeentries << endl;
    cout << "reconstructed y entries: " << recyentries << endl;
    cout << "reconstructed pT entries: " << recptentries << endl;

    ks_dNdlife_cor_data->Divide(ks_dNdlife_rec_data, ks_dNdlife_eff_mc);
    ks_dNdy_cor_data->Divide(ks_dNdy_rec_data, ks_dNdy_eff_mc);
    ks_dNdy_pos_cor_data->Divide(ks_dNdy_pos_rec_data, ks_dNdy_pos_eff_mc);
    ks_dNdy_neg_cor_data->Divide(ks_dNdy_neg_rec_data, ks_dNdy_neg_eff_mc);
    ks_dNdpT_cor_data->Divide(ks_dNdpT_rec_data, ks_dNdpT_eff_mc);

    TCanvas* c1 = new TCanvas("c1", "K0S corrected yields");

    Double_t ksp0 = log( ks_dNdlife_cor_data->GetMaximum() );

    ksctFit->SetParameter(0, ksp0);
    ksctFit->SetParameter(1, ksPDGlife);

    ks_dNdlife_cor_data->Fit("ksctFit", "EIM");

    //ks_dNdlife_cor_data->SetAxisRange(0., 500000., "Y");
    //ks_dNdlife_cor_data->SetMinimum(0.);
    ks_dNdlife_cor_data->GetXaxis()->SetNdivisions(508);
    ks_dNdlife_cor_data->SetLineColor(1);
    ks_dNdlife_cor_data->SetLineWidth(1);
    ks_dNdlife_cor_data->GetXaxis()->SetTitle("K_{S}^{0}  ct  [cm]");
    ks_dNdlife_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1 cm");
    ks_dNdlife_cor_data->SetMarkerSize(0.8);
    ks_dNdlife_cor_data->Draw();
    //gPad->SetLogy(1);

    if( !isMC ) {
      c1->SaveAs(kscorLtimeFname_data + ".png");
      c1->SaveAs(kscorLtimeFname_data + ".eps");
      c1->SaveAs(kscorLtimeFname_data + ".root");
    }
    else {
      c1->SaveAs(kscorLtimeFname + ".png");
      c1->SaveAs(kscorLtimeFname + ".eps");
      c1->SaveAs(kscorLtimeFname + ".root");
    }
    //gPad->SetLogy(0);

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "K0S corrected yields");

    ks_dNdy_cor_data->SetAxisRange(0., ks_dNdy_cor_data->GetMaximum()*1.4, "Y");
    ks_dNdy_cor_data->GetXaxis()->SetNdivisions(508);
    ks_dNdy_cor_data->SetLineColor(1);
    ks_dNdy_cor_data->SetLineWidth(1);
    ks_dNdy_cor_data->GetXaxis()->SetTitle("K_{S}^{0}  rapidity");
    ks_dNdy_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1");
    ks_dNdy_cor_data->SetMarkerSize(0.8);
    ks_dNdy_cor_data->Draw();

    if( !isMC ) {
      c1->SaveAs(kscoryFname_data + ".png");
      c1->SaveAs(kscoryFname_data + ".eps");
      c1->SaveAs(kscoryFname_data + ".root");
    }
    else {
      c1->SaveAs(kscoryFname + ".png");
      c1->SaveAs(kscoryFname + ".eps");
      c1->SaveAs(kscoryFname + ".root");
    }

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "K0S corrected yields");

    ks_dNdy_neg_cor_data->SetAxisRange(0., ks_dNdy_neg_cor_data->GetMaximum()*1.4, "Y");
    ks_dNdy_neg_cor_data->GetXaxis()->SetNdivisions(508);
    ks_dNdy_neg_cor_data->SetLineColor(1);
    ks_dNdy_neg_cor_data->SetLineWidth(1);
    ks_dNdy_neg_cor_data->GetXaxis()->SetTitle("K_{S}^{0}  rapidity");
    ks_dNdy_neg_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1");
    ks_dNdy_neg_cor_data->SetMarkerSize(0.8);
    ks_dNdy_neg_cor_data->Draw();

    if( !isMC ) {
      c1->SaveAs(kscory_negFname_data + ".png");
      c1->SaveAs(kscory_negFname_data + ".eps");
      c1->SaveAs(kscory_negFname_data + ".root");
    }
    else {
      c1->SaveAs(kscory_negFname + ".png");
      c1->SaveAs(kscory_negFname + ".eps");
      c1->SaveAs(kscory_negFname + ".root");
    }

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "K0S corrected yields");

    ks_dNdy_pos_cor_data->SetAxisRange(0., ks_dNdy_pos_cor_data->GetMaximum()*1.4, "Y");
    ks_dNdy_pos_cor_data->GetXaxis()->SetNdivisions(508);
    ks_dNdy_pos_cor_data->SetLineColor(1);
    ks_dNdy_pos_cor_data->SetLineWidth(1);
    ks_dNdy_pos_cor_data->GetXaxis()->SetTitle("K_{S}^{0}  rapidity");
    ks_dNdy_pos_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1");
    ks_dNdy_pos_cor_data->SetMarkerSize(0.8);
    ks_dNdy_pos_cor_data->Draw();

    if( !isMC ) {
      c1->SaveAs(kscory_posFname_data + ".png");
      c1->SaveAs(kscory_posFname_data + ".eps");
      c1->SaveAs(kscory_posFname_data + ".root");
    }
    else {
      c1->SaveAs(kscory_posFname + ".png");
      c1->SaveAs(kscory_posFname + ".eps");
      c1->SaveAs(kscory_posFname + ".root");
    }

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "K0S corrected yields");

    //ks_dNdpT_cor_data->SetAxisRange(0., ks_dNdpT_cor_data->GetMaximum()*1.05, "Y");
    ks_dNdpT_cor_data->SetMaximum(ks_dNdpT_cor_data->GetMaximum()*1.05);
    ks_dNdpT_cor_data->GetXaxis()->SetNdivisions(508);
    ks_dNdpT_cor_data->SetLineColor(1);
    ks_dNdpT_cor_data->SetLineWidth(1);
    ks_dNdpT_cor_data->GetXaxis()->SetTitle("K_{S}^{0}  p_{T}  [GeV/c]");
    ks_dNdpT_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1 GeV/c");
    ks_dNdpT_cor_data->SetMarkerSize(0.8);
    ks_dNdpT_cor_data->Draw();

    if( !isMC ) {
      c1->SaveAs(kscorpTFname_data + ".png");
      c1->SaveAs(kscorpTFname_data + ".eps");
      c1->SaveAs(kscorpTFname_data + ".root");
    }
    else {
      c1->SaveAs(kscorpTFname + ".png");
      c1->SaveAs(kscorpTFname + ".eps");
      c1->SaveAs(kscorpTFname + ".root");
    }

    delete c1;

  }

  else {
//Lambda

    if( isMC ) {
      TString lamlife_rec = "lam_" + algostr + "_" + energy2 + "/plots/lam_dNdlife_MC_" + algostr + ".root";
      TString lamy_rec = "lam_" + algostr + "_" + energy2 + "/plots/lam_dNdy_MC_" + algostr + ".root";
      TString lamRegy_rec = "lam_" + algostr + "_" + energy2 + "/plots/lamReg_dNdy_MC_" + algostr + ".root";
      TString lamBary_rec = "lam_" + algostr + "_" + energy2 + "/plots/lamBar_dNdy_MC_" + algostr + ".root";
      TString lampt_rec = "lam_" + algostr + "_" + energy2 + "/plots/lam_dNdpT_MC_" + algostr + ".root";
      TString lamRegpt_rec = "lam_" + algostr + "_" + energy2 + "/plots/lamReg_dNdpT_MC_" + algostr + ".root";
      TString lamBarpt_rec = "lam_" + algostr + "_" + energy2 + "/plots/lamBar_dNdpT_MC_" + algostr + ".root";
    }
    else {
      lamlife_rec = "lam_data_" + energy2 + "/plots/lam_dNdlife.root";
      lamy_rec = "lam_data_" + energy2 + "/plots/lam_dNdy.root";
      lamRegy_rec = "lam_data_" + energy2 + "/plots/lamReg_dNdy.root";
      lamBary_rec = "lam_data_" + energy2 + "/plots/lamBar_dNdy.root";
      lampt_rec = "lam_data_" + energy2 + "/plots/lam_dNdpT.root";
      lamRegpt_rec = "lam_data_" + energy2 + "/plots/lamReg_dNdpT.root";
      lamBarpt_rec = "lam_data_" + energy2 + "/plots/lamBar_dNdpT.root";
    }

    lamlife_eff = "lam_" + algostr + "_" + energy2 + "/plots/lam_dNdlife_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2.root";
    lamy_eff = "lam_" + algostr + "_" + energy2 + "/plots/lam_dNdy_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2.root";
    lamRegy_eff = "lam_" + algostr + "_" + energy2 + "/plots/lamReg_dNdy_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2.root";
    lamBary_eff = "lam_" + algostr + "_" + energy2 + "/plots/lamBar_dNdy_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2.root";
    lampt_eff = "lam_" + algostr + "_" + energy2 + "/plots/lam_dNdpT_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2.root";
    lamRegpt_eff = "lam_" + algostr + "_" + energy2 + "/plots/lamReg_dNdpT_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2.root";
    lamBarpt_eff = "lam_" + algostr + "_" + energy2 + "/plots/lamBar_dNdpT_eff_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2.root";

    TString lamcorLtimeFname;
    TString lamcoryFname;
    TString lamRegcoryFname;
    TString lamBarcoryFname;
    TString lamcorpTFname;
    TString lamRegcorpTFname;
    TString lamBarcorpTFname;
    TString lamcorLtimeFname_data;
    TString lamcoryFname_data;
    TString lamRegcoryFname_data;
    TString lamBarcoryFname_data;
    TString lamcorpTFname_data;
    TString lamRegcorpTFname_data;
    TString lamBarcorpTFname_data;
    lamcorLtimeFname = "plots/lamlife_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamcoryFname = "plots/lamy_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamRegcoryFname = "plots/lamRegy_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamBarcoryFname = "plots/lamBary_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";

    lamcorpTFname = "plots/lampt_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamRegcorpTFname = "plots/lamRegpt_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamBarcorpTFname = "plots/lamBarpt_cor_mc_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamcorLtimeFname_data = "plots/lamlife_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamcoryFname_data = "plots/lamy_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamcorpTFname_data = "plots/lampt_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamRegcoryFname_data = "plots/lamRegy_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamRegcorpTFname_data = "plots/lamRegpt_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamBarcoryFname_data = "plots/lamBary_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";
    lamBarcorpTFname_data = "plots/lamBarpt_cor_data_" + algostr + "_" + energy + "_weight" + weightoss.str().c_str() + "_subgeant2";

    /*
      cout << lamlife_rec << endl;
      cout << lamy_rec << endl;
      cout << lampt_rec << endl;
      cout << lamlife_eff << endl;
      cout << lamy_eff << endl;
      cout << lampt_eff << endl;
    */

    // Get histos to correct

    TFile* fin = TFile::Open(lamlife_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lam_dNdlife_rec_data = (TH1F*)c1->GetPrimitive("lam_dNdlife_scale");
    TH1F* lam_dNdlife_cor_data = (TH1F*)c1->GetPrimitive("lam_dNdlife_scale");
    double reclifeentries = lam_dNdlife_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamy_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lam_dNdy_rec_data = (TH1F*)c1->GetPrimitive("lam_dNdy");
    TH1F* lam_dNdy_cor_data = (TH1F*)c1->GetPrimitive("lam_dNdy");
    double recyentries = lam_dNdy_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamRegy_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lamReg_dNdy_rec_data = (TH1F*)c1->GetPrimitive("lamReg_dNdy");
    TH1F* lamReg_dNdy_cor_data = (TH1F*)c1->GetPrimitive("lamReg_dNdy");
    double recRyentries = lamReg_dNdy_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamBary_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lamBar_dNdy_rec_data = (TH1F*)c1->GetPrimitive("lamBar_dNdy");
    TH1F* lamBar_dNdy_cor_data = (TH1F*)c1->GetPrimitive("lamBar_dNdy");
    double recByentries = lamBar_dNdy_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lampt_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lam_dNdpT_rec_data = (TH1F*)c1->GetPrimitive("lam_dNdpT_scale");
    TH1F* lam_dNdpT_cor_data = (TH1F*)c1->GetPrimitive("lam_dNdpT_scale");
    double recptentries = lam_dNdpT_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamRegpt_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lamReg_dNdpT_rec_data = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_scale");
    TH1F* lamReg_dNdpT_cor_data = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_scale");
    double recptentries = lamReg_dNdpT_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamBarpt_rec, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lamBar_dNdpT_rec_data = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_scale");
    TH1F* lamBar_dNdpT_cor_data = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_scale");
    double recptentries = lamBar_dNdpT_rec_data->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    // Get efficiency histos

    TFile* fin = TFile::Open(lamlife_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lam_dNdlife_eff_mc = (TH1F*)c1->GetPrimitive("lam_dNdlife_eff");
    double mcefflifeentries = lam_dNdlife_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamy_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lam_dNdy_eff_mc = (TH1F*)c1->GetPrimitive("lam_dNdy_eff");
    double mceffyentries = lam_dNdy_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamRegy_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lamReg_dNdy_eff_mc = (TH1F*)c1->GetPrimitive("lamReg_dNdy_eff");
    double mceffRyentries = lamReg_dNdy_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamBary_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lamBar_dNdy_eff_mc = (TH1F*)c1->GetPrimitive("lamBar_dNdy_eff");
    double mceffByentries = lamBar_dNdy_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lampt_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lam_dNdpT_eff_mc = (TH1F*)c1->GetPrimitive("lam_dNdpT_eff");
    double mceffptentries = lam_dNdpT_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamRegpt_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lamReg_dNdpT_eff_mc = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_eff");
    double mceffRptentries = lamReg_dNdpT_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    TFile* fin = TFile::Open(lamBarpt_eff, "READ");
    TCanvas* c1 = (TCanvas*)fin.Get("c1");
    TH1F* lamBar_dNdpT_eff_mc = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_eff");
    double mceffBptentries = lamBar_dNdpT_eff_mc->GetSumOfWeights();
    fin->Close();
    delete c1;
    delete fin;

    cout << "reconstructed lifetime entries: " << reclifeentries << endl;
    cout << "reconstructed y entries: " << recyentries << endl;
    cout << "reconstructed pT entries: " << recptentries << endl;

    lam_dNdlife_cor_data->Divide(lam_dNdlife_rec_data, lam_dNdlife_eff_mc);
    lam_dNdy_cor_data->Divide(lam_dNdy_rec_data, lam_dNdy_eff_mc);
    lamReg_dNdy_cor_data->Divide(lamReg_dNdy_rec_data, lamReg_dNdy_eff_mc);
    lamBar_dNdy_cor_data->Divide(lamBar_dNdy_rec_data, lamBar_dNdy_eff_mc);
    lam_dNdpT_cor_data->Divide(lam_dNdpT_rec_data, lam_dNdpT_eff_mc);
    lamReg_dNdpT_cor_data->Divide(lamReg_dNdpT_rec_data, lamReg_dNdpT_eff_mc);
    lamBar_dNdpT_cor_data->Divide(lamBar_dNdpT_rec_data, lamBar_dNdpT_eff_mc);

    TCanvas* c1 = new TCanvas("c1", "Lambda corrected yields");

    Double_t lamp0 = log( lam_dNdlife_cor_data->GetMaximum() );

    lamctFit->SetParameter(0, lamp0);
    lamctFit->SetParameter(1, lamPDGlife);

    lam_dNdlife_cor_data->Fit("lamctFit", "EIM");

    //lam_dNdlife_cor_data->SetAxisRange(0., 500000., "Y");
    //lam_dNdlife_cor_data->SetMinimum(0.);
    lam_dNdlife_cor_data->GetXaxis()->SetNdivisions(508);
    lam_dNdlife_cor_data->SetLineColor(1);
    lam_dNdlife_cor_data->SetLineWidth(1);
    lam_dNdlife_cor_data->GetXaxis()->SetTitle("#Lambda^{0}  ct  [cm]");
    lam_dNdlife_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1 cm");
    lam_dNdlife_cor_data->SetMarkerSize(0.8);
    lam_dNdlife_cor_data->Draw();
    //gPad->SetLogy(1);

    if( !isMC ) {
      c1->SaveAs(lamcorLtimeFname_data + ".png");
      c1->SaveAs(lamcorLtimeFname_data + ".eps");
      c1->SaveAs(lamcorLtimeFname_data + ".root");
    }
    else {
      c1->SaveAs(lamcorLtimeFname + ".png");
      c1->SaveAs(lamcorLtimeFname + ".eps");
      c1->SaveAs(lamcorLtimeFname + ".root");
    }
    //gPad->SetLogy(0);

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "Lambda corrected yields");

    lam_dNdy_cor_data->SetAxisRange(0., lam_dNdy_cor_data->GetMaximum()*1.4, "Y");
    lam_dNdy_cor_data->GetXaxis()->SetNdivisions(508);
    lam_dNdy_cor_data->SetLineColor(1);
    lam_dNdy_cor_data->SetLineWidth(1);
    lam_dNdy_cor_data->GetXaxis()->SetTitle("#Lambda^{0}  rapidity");
    lam_dNdy_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1");
    lam_dNdy_cor_data->SetMarkerSize(0.8);
    lam_dNdy_cor_data->Draw();

    if( !isMC ) {
      c1->SaveAs(lamcoryFname_data + ".png");
      c1->SaveAs(lamcoryFname_data + ".eps");
      c1->SaveAs(lamcoryFname_data + ".root");
    }
    else {
      c1->SaveAs(lamcoryFname + ".png");
      c1->SaveAs(lamcoryFname + ".eps");
      c1->SaveAs(lamcoryFname + ".root");
    }

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "Lambda corrected yields");

    lamReg_dNdy_cor_data->SetAxisRange(0., lamReg_dNdy_cor_data->GetMaximum()*1.4, "Y");
    lamReg_dNdy_cor_data->GetXaxis()->SetNdivisions(508);
    lamReg_dNdy_cor_data->SetLineColor(1);
    lamReg_dNdy_cor_data->SetLineWidth(1);
    lamReg_dNdy_cor_data->GetXaxis()->SetTitle("#Lambda^{0}  rapidity");
    lamReg_dNdy_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1");
    lamReg_dNdy_cor_data->SetMarkerSize(0.8);
    lamReg_dNdy_cor_data->Draw();

    if( !isMC ) {
      c1->SaveAs(lamRegcoryFname_data + ".png");
      c1->SaveAs(lamRegcoryFname_data + ".eps");
      c1->SaveAs(lamRegcoryFname_data + ".root");
    }
    else {
      c1->SaveAs(lamRegcoryFname + ".png");
      c1->SaveAs(lamRegcoryFname + ".eps");
      c1->SaveAs(lamRegcoryFname + ".root");
    }

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "LambdaBar corrected yields");

    lamBar_dNdy_cor_data->SetAxisRange(0., lamBar_dNdy_cor_data->GetMaximum()*1.4, "Y");
    lamBar_dNdy_cor_data->GetXaxis()->SetNdivisions(508);
    lamBar_dNdy_cor_data->SetLineColor(1);
    lamBar_dNdy_cor_data->SetLineWidth(1);
    lamBar_dNdy_cor_data->GetXaxis()->SetTitle("#Lambda^{0}  rapidity");
    lamBar_dNdy_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1");
    lamBar_dNdy_cor_data->SetMarkerSize(0.8);
    lamBar_dNdy_cor_data->Draw();

    if( !isMC ) {
      c1->SaveAs(lamBarcoryFname_data + ".png");
      c1->SaveAs(lamBarcoryFname_data + ".eps");
      c1->SaveAs(lamBarcoryFname_data + ".root");
    }
    else {
      c1->SaveAs(lamBarcoryFname + ".png");
      c1->SaveAs(lamBarcoryFname + ".eps");
      c1->SaveAs(lamBarcoryFname + ".root");
    }

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "Lambda corrected yields");

    //lam_dNdpT_cor_data->SetAxisRange(0., lam_dNdpT_cor_data->GetMaximum()*1.05, "Y");
    lam_dNdpT_cor_data->SetMaximum(lam_dNdpT_cor_data->GetMaximum()*1.05);
    lam_dNdpT_cor_data->GetXaxis()->SetNdivisions(508);
    lam_dNdpT_cor_data->SetLineColor(1);
    lam_dNdpT_cor_data->SetLineWidth(1);
    lam_dNdpT_cor_data->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
    lam_dNdpT_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1 GeV/c");
    lam_dNdpT_cor_data->SetMarkerSize(0.8);
    lam_dNdpT_cor_data->Draw();

    if( !isMC ) {
      c1->SaveAs(lamcorpTFname_data + ".png");
      c1->SaveAs(lamcorpTFname_data + ".eps");
      c1->SaveAs(lamcorpTFname_data + ".root");
    }
    else {
      c1->SaveAs(lamcorpTFname + ".png");
      c1->SaveAs(lamcorpTFname + ".eps");
      c1->SaveAs(lamcorpTFname + ".root");
    }

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "Lambda corrected yields");

    //lam_dNdpT_cor_data->SetAxisRange(0., lam_dNdpT_cor_data->GetMaximum()*1.05, "Y");
    lamReg_dNdpT_cor_data->SetMaximum(lamReg_dNdpT_cor_data->GetMaximum()*1.05);
    lamReg_dNdpT_cor_data->GetXaxis()->SetNdivisions(508);
    lamReg_dNdpT_cor_data->SetLineColor(1);
    lamReg_dNdpT_cor_data->SetLineWidth(1);
    lamReg_dNdpT_cor_data->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
    lamReg_dNdpT_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1 GeV/c");
    lamReg_dNdpT_cor_data->SetMarkerSize(0.8);
    lamReg_dNdpT_cor_data->Draw();

    if( !isMC ) {
      c1->SaveAs(lamRegcorpTFname_data + ".png");
      c1->SaveAs(lamRegcorpTFname_data + ".eps");
      c1->SaveAs(lamRegcorpTFname_data + ".root");
    }
    else {
      c1->SaveAs(lamRegcorpTFname + ".png");
      c1->SaveAs(lamRegcorpTFname + ".eps");
      c1->SaveAs(lamRegcorpTFname + ".root");
    }

    delete c1;

    TCanvas* c1 = new TCanvas("c1", "LambdaBar corrected yields");

    //lam_dNdpT_cor_data->SetAxisRange(0., lam_dNdpT_cor_data->GetMaximum()*1.05, "Y");
    lamBar_dNdpT_cor_data->SetMaximum(lamBar_dNdpT_cor_data->GetMaximum()*1.05);
    lamBar_dNdpT_cor_data->GetXaxis()->SetNdivisions(508);
    lamBar_dNdpT_cor_data->SetLineColor(1);
    lamBar_dNdpT_cor_data->SetLineWidth(1);
    lamBar_dNdpT_cor_data->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
    lamBar_dNdpT_cor_data->GetYaxis()->SetTitle("Corrected candidates / 0.1 GeV/c");
    lamBar_dNdpT_cor_data->SetMarkerSize(0.8);
    lamBar_dNdpT_cor_data->Draw();

    if( !isMC ) {
      cout << lamBarcorpTFname_data << endl;
      c1->SaveAs(lamBarcorpTFname_data + ".png");
      c1->SaveAs(lamBarcorpTFname_data + ".eps");
      c1->SaveAs(lamBarcorpTFname_data + ".root");
    }
    else {
      c1->SaveAs(lamBarcorpTFname + ".png");
      c1->SaveAs(lamBarcorpTFname + ".eps");
      c1->SaveAs(lamBarcorpTFname + ".root");
    }

    delete c1;

  }


}
