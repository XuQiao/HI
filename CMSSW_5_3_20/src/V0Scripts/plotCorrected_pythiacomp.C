{
#include <tdrstyle.C>
#include "tsallis_keith.C"
//#include "tsallis.C"
//#include "tsallismean.C"
//#include "tsallisint.C"
#include "calc_integral.C"
// basic file operations
using namespace std;

  // Style and other ROOT commands
setTDRStyle();
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
TGaxis::SetMaxDigits(3);
//gROOT->SetBatch(kTRUE);
gStyle->SetErrorX(0.5);
gStyle->SetTitleSize(0.05, "XYZ");
gStyle->SetLabelSize(0.045, "XYZ");
gStyle->SetTitleXOffset(1.1);
gStyle->SetTitleYOffset(1.5);
std::cout.precision(3);

Double_t mks = 0.497614;
Double_t mlam = 1.115683;

// /nfs/data37/cms/v0ntuple/mc_7TeV_tuneD6T_50M_START3X_V26B_june3fixes_PARTIAL.root
//double mcevt_d6t_7000 = 0.8081*37425190;
// /nfs/data37/cms/v0ntuple/mc_7TeV_tuneD6T_50M_START3X_V26B_june4fixesGood.root
double mcevt_d6t_7000 = 30610982;
// mc_7TeV_pythia8_START3X_V26B_june3fixes_PARTIAL.root
//double mcevt_py8_7000 = 0.8083*9814844;
// /nfs/data37/cms/v0ntuple/mc_7TeV_pythia8_START3X_V26B_june4fixesGood.root
double mcevt_py8_7000 = 8701075;
// /nfs/data37/cms/v0ntuple/mc_7TeV_tuneP0_START3X_V26B_june3fixes.root
//double mcevt_p0_7000 = 0.8083*10932196;
// /nfs/data37/cms/v0ntuple/mc_7TeV_tuneP0_START3X_V26B_june4fixesGood.root
double mcevt_p0_7000 = 8837092;
// /nfs/data37/cms/v0ntuple/mc_900GeV_START3X_V22A_june4fixes.root
//double mcevt_d6t_900 = 0.7762*10101170;
// /nfs/data37/cms/v0ntuple/mc_900GeV_START3X_V22A_june4fixesGood_PARTIAL.root
double mcevt_d6t_900 = 7840577;


double trig7000 = 15540204;
double trig900 = 12325647;
double lum7000 = 291850.0; // inverse millibarns
double lum900 = 321970.0; // inverse millibarns
double ptbinwid = 0.1;
double ybinwid = 0.2; // 0.2 since it is a 0.1 bin on each side of 0.
double lifebinwid = 0.1;

TF1 *tsfunc900 = new TF1("tsfunc900",tsallis,0.0,10.0,3);
tsfunc900->SetParNames("Amplitude","T","n");
tsfunc900->SetParameters(1.0.,0.2,7.2); 
TF1 *tsfunc900d6t = new TF1("tsfunc900d6t",tsallis,0.0,10.0,3);
tsfunc900d6t->SetParNames("Amplitude","T","n");
tsfunc900d6t->SetParameters(1.0.,0.2,7.2); 
TF1 *tsfunc7000 = new TF1("tsfunc7000",tsallis,0.0,10.0,3);
tsfunc7000->SetParNames("Amplitude","T","n");
tsfunc7000->SetParameters(1.0.,0.2,7.2); 
TF1 *tsfunc7000d6t = new TF1("tsfunc7000d6t",tsallis,0.0,10.0,3);
tsfunc7000d6t->SetParNames("Amplitude","T","n");
tsfunc7000d6t->SetParameters(1.0.,0.2,7.2); 
TF1 *tsfunc7000py8 = new TF1("tsfunc7000py8",tsallis,0.0,10.0,3);
tsfunc7000py8->SetParNames("Amplitude","T","n");
tsfunc7000py8->SetParameters(1.0.,0.2,7.2); 
TF1 *tsfunc7000p0 = new TF1("tsfunc7000p0",tsallis,0.0,10.0,3);
tsfunc7000p0->SetParNames("Amplitude","T","n");
tsfunc7000p0->SetParameters(1.0.,0.2,7.2); 

TF1 *tsmean900 = new TF1("tsmean900",tsallismean,0.0,10.0,3);
tsmean900->SetParNames("Amplitude","T","n");
tsmean900->SetParameters(1.0.,0.9,7.2); 
TF1 *tsmean900d6t = new TF1("tsmean900d6t",tsallismean,0.0,10.0,3);
tsmean900d6t->SetParNames("Amplitude","T","n");
tsmean900d6t->SetParameters(1.0.,0.9,7.2); 
TF1 *tsmean7000 = new TF1("tsmean7000",tsallismean,0.0,10.0,3);
tsmean7000->SetParNames("Amplitude","T","n");
tsmean7000->SetParameters(1.0.,0.9,7.2); 
TF1 *tsmean7000d6t = new TF1("tsmean7000d6t",tsallismean,0.0,10.0,3);
tsmean7000d6t->SetParNames("Amplitude","T","n");
tsmean7000d6t->SetParameters(1.0.,0.9,7.2); 
TF1 *tsmean7000py8 = new TF1("tsmean7000py8",tsallismean,0.0,10.0,3);
tsmean7000py8->SetParNames("Amplitude","T","n");
tsmean7000py8->SetParameters(1.0.,0.9,7.2); 
TF1 *tsmean7000p0 = new TF1("tsmean7000p0",tsallismean,0.0,10.0,3);
tsmean7000p0->SetParNames("Amplitude","T","n");
tsmean7000p0->SetParameters(1.0.,0.9,7.2); 

TLatex *cmsprelim = new TLatex();
cmsprelim->SetNDC();
cmsprelim->SetTextSize(0.045);
cmsprelim->SetTextAlign(22);

TString inFile;
TString outFile_pt_evt = "final_plots/kspt_cor_datamc_perevt";

Double_t ptent[43];

TString outFile_y_evt = "final_plots/ksy_cor_datamc_perevt";

Double_t yent[23];


// Pt for Kshort per event

inFile = "plots/kspt_cor_data_pythia8_7TeV_weight102_nobinomial2.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdpT_cor_data_7000 = (TH1F*)c1->GetPrimitive("ks_dNdpT_scale");
fin->Close();
delete c1;
delete fin;
ks_dNdpT_cor_data_7000->Scale(1.0/(trig7000*ptbinwid));
double ks_yld_7000_pt = calc_integral(ks_dNdpT_cor_data_7000, 1.0);
ks_dNdpT_cor_data_7000->GetXaxis()->SetTitle("K_{S}^{0}  p_{T}  [GeV/c]");
ks_dNdpT_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dp_{T}  (GeV/c)^{-1}");
ks_dNdpT_cor_data_7000->SetLineColor(4);
ks_dNdpT_cor_data_7000->SetMarkerColor(4);
ks_dNdpT_cor_data_7000->SetMarkerStyle(20);
ks_dNdpT_cor_data_7000->SetMarkerSize(0.8);
TH1F *ks_dNdpT_cor_data_7000_sys = (TH1F*)ks_dNdpT_cor_data_7000->Clone("ks_dNdpT_cor_data_7000_sys");
cout << "ks pt 7000 integral / event: " << ks_yld_7000_pt << endl;
for (int i = 1; i <= ks_dNdpT_cor_data_7000_sys->GetXaxis()->GetNbins(); ++i) {
  double value = ks_dNdpT_cor_data_7000_sys->GetBinContent(i);
  double error = ks_dNdpT_cor_data_7000_sys->GetBinError(i);
  double syserr = 0.02*value;
  ptent[i-1] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdpT_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i-1] << endl;
}
ks_dNdpT_cor_data_7000_sys->SetError(ptent);
ks_dNdpT_cor_data_7000->Fit("tsmean7000","EM0");
ks_dNdpT_cor_data_7000->Fit("tsfunc7000","EM0");
tsfunc7000->SetNpx(2000);
tsfunc7000->SetLineColor(4);


inFile = "ks_tuneD6T_7000_new2/plots/ks_dNdpT_gen_tuneD6T_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdpT_cor_d6t_7000 = (TH1F*)c1->GetPrimitive("ks_dNdpT_gen_scale");
fin->Close();
delete c1;
delete fin;
ks_dNdpT_cor_d6t_7000->Scale(1.0/(mcevt_d6t_7000*ptbinwid));
double ks_yld_d6t_7000_pt = calc_integral(ks_dNdpT_cor_d6t_7000, 1.0);
ks_dNdpT_cor_d6t_7000->GetXaxis()->SetTitle("K_{S}^{0}  p_{T}  [GeV/c]");
ks_dNdpT_cor_d6t_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dp_{T}  (GeV/c)^{-1}");
ks_dNdpT_cor_d6t_7000->SetLineColor(4);
ks_dNdpT_cor_d6t_7000->SetMarkerStyle(1);
ks_dNdpT_cor_d6t_7000->SetLineStyle(1);
cout << "ks D6T pt 7000 integral / event: " << ks_yld_d6t_7000_pt << endl;
ks_dNdpT_cor_d6t_7000->Fit("tsmean7000d6t","EM0");
ks_dNdpT_cor_d6t_7000->Fit("tsfunc7000d6t","EM0");
tsfunc7000d6t->SetNpx(2000);
tsfunc7000d6t->SetLineColor(4);
tsfunc7000d6t->SetLineStyle(1);
TH1F *ks_dNdpT_data_d6t_7000 = (TH1F*)ks_dNdpT_cor_d6t_7000->Clone("ks_dNdpT_data_d6t_7000");
ks_dNdpT_data_d6t_7000->Divide(ks_dNdpT_cor_data_7000_sys);

inFile = "ks_pythia8_7000_new2/plots/ks_dNdpT_gen_pythia8_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdpT_cor_py8_7000 = (TH1F*)c1->GetPrimitive("ks_dNdpT_gen_scale");
fin->Close();
delete c1;
delete fin;
ks_dNdpT_cor_py8_7000->Scale(1.0/(mcevt_py8_7000*ptbinwid));
double ks_yld_py8_7000_pt = calc_integral(ks_dNdpT_cor_py8_7000, 1.0);
ks_dNdpT_cor_py8_7000->GetXaxis()->SetTitle("K_{S}^{0}  p_{T}  [GeV/c]");
ks_dNdpT_cor_py8_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dp_{T}  (GeV/c)^{-1}");
ks_dNdpT_cor_py8_7000->SetLineColor(4);
ks_dNdpT_cor_py8_7000->SetMarkerStyle(1);
ks_dNdpT_cor_py8_7000->SetLineStyle(2);
cout << "ks PY8 pt 7000 integral / event: " << ks_yld_py8_7000_pt << endl;
ks_dNdpT_cor_py8_7000->Fit("tsmean7000py8","EM0");
ks_dNdpT_cor_py8_7000->Fit("tsfunc7000py8","EM0");
tsfunc7000py8->SetNpx(2000);
tsfunc7000py8->SetLineColor(4);
tsfunc7000py8->SetLineStyle(2);
TH1F *ks_dNdpT_data_py8_7000 = (TH1F*)ks_dNdpT_cor_py8_7000->Clone("ks_dNdpT_data_py8_7000");
ks_dNdpT_data_py8_7000->Divide(ks_dNdpT_cor_data_7000_sys);

inFile = "ks_tuneP0_7000_new2/plots/ks_dNdpT_gen_tuneP0_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdpT_cor_p0_7000 = (TH1F*)c1->GetPrimitive("ks_dNdpT_gen_scale");
fin->Close();
delete c1;
delete fin;
ks_dNdpT_cor_p0_7000->Scale(1.0/(mcevt_p0_7000*ptbinwid));
double ks_yld_p0_7000_pt = calc_integral(ks_dNdpT_cor_p0_7000, 1.0);
ks_dNdpT_cor_p0_7000->GetXaxis()->SetTitle("K_{S}^{0}  p_{T}  [GeV/c]");
ks_dNdpT_cor_p0_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dp_{T}  (GeV/c)^{-1}");
ks_dNdpT_cor_p0_7000->SetLineColor(4);
ks_dNdpT_cor_p0_7000->SetMarkerStyle(1);
ks_dNdpT_cor_p0_7000->SetLineStyle(3);
cout << "ks P0 pt 7000 integral / event: " << ks_yld_p0_7000_pt << endl;
ks_dNdpT_cor_p0_7000->Fit("tsmean7000p0","EM0");
ks_dNdpT_cor_p0_7000->Fit("tsfunc7000p0","EM0");
tsfunc7000p0->SetNpx(2000);
tsfunc7000p0->SetLineColor(4);
tsfunc7000p0->SetLineStyle(3);
TH1F *ks_dNdpT_data_p0_7000 = (TH1F*)ks_dNdpT_cor_p0_7000->Clone("ks_dNdpT_data_p0_7000");
ks_dNdpT_data_p0_7000->Divide(ks_dNdpT_cor_data_7000_sys);



inFile = "plots/kspt_cor_data_tuneD6T_900GeV_weight111_nobinomial2.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdpT_cor_data_900 = (TH1F*)c1->GetPrimitive("ks_dNdpT_scale");
fin->Close();
delete c1;
delete fin;
double yld_900_pt = ks_dNdpT_cor_data_900->Integral();
ks_dNdpT_cor_data_900->Scale(1.0/(trig900*ptbinwid));
double ks_yld_900_pt = calc_integral(ks_dNdpT_cor_data_900, 1.0);
ks_dNdpT_cor_data_900->GetXaxis()->SetTitle("K_{S}^{0}  p_{T}  [GeV/c]");
ks_dNdpT_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dp_{T}  (GeV/c)^{-1}");
ks_dNdpT_cor_data_900->SetLineColor(2);
ks_dNdpT_cor_data_900->SetMarkerColor(2);
ks_dNdpT_cor_data_900->SetMarkerStyle(21);
ks_dNdpT_cor_data_900->SetMarkerSize(0.8);
TH1F *ks_dNdpT_cor_data_900_sys = (TH1F*)ks_dNdpT_cor_data_900->Clone("ks_dNdpT_cor_data_7000_sys");
//cout << "ks pt 900 integral / event: " << ks_yld_900_pt << endl;
for (int i = 1; i <= ks_dNdpT_cor_data_900_sys->GetXaxis()->GetNbins(); ++i) {
  double value = ks_dNdpT_cor_data_900_sys->GetBinContent(i);
  double error = ks_dNdpT_cor_data_900_sys->GetBinError(i);
  double syserr = 0.02*value;
  ptent[i-1] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdpT_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i-1] << endl;
}
ks_dNdpT_cor_data_900_sys->SetError(ptent);
ks_dNdpT_cor_data_900->Fit("tsmean900","EM0");
ks_dNdpT_cor_data_900->Fit("tsfunc900","EM0");
tsfunc900->SetNpx(2000);
tsfunc900->SetLineColor(2);

inFile = "ks_tuneD6T_900_new2/plots/ks_dNdpT_gen_tuneD6T_900GeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdpT_cor_d6t_900 = (TH1F*)c1->GetPrimitive("ks_dNdpT_gen_scale");
fin->Close();
delete c1;
delete fin;
ks_dNdpT_cor_d6t_900->Scale(1.0/(mcevt_d6t_900*ptbinwid));
double ks_yld_d6t_900_pt = calc_integral(ks_dNdpT_cor_d6t_900, 1.0);
ks_dNdpT_cor_d6t_900->GetXaxis()->SetTitle("K_{S}^{0}  p_{T}  [GeV/c]");
ks_dNdpT_cor_d6t_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dp_{T}  (GeV/c)^{-1}");
ks_dNdpT_cor_d6t_900->SetLineColor(2);
ks_dNdpT_cor_d6t_900->SetMarkerStyle(1);
ks_dNdpT_cor_d6t_900->SetLineStyle(1);
cout << "ks D6T pt 900 integral / event: " << ks_yld_d6t_900_pt << endl;
ks_dNdpT_cor_d6t_900->Fit("tsmean900d6t","EM0");
ks_dNdpT_cor_d6t_900->Fit("tsfunc900d6t","EM0");
tsfunc900d6t->SetNpx(2000);
tsfunc900d6t->SetLineColor(1);
tsfunc900d6t->SetLineStyle(2);
TH1F *ks_dNdpT_data_d6t_900 = (TH1F*)ks_dNdpT_cor_d6t_900->Clone("ks_dNdpT_data_d6t_900");
ks_dNdpT_data_d6t_900->Divide(ks_dNdpT_cor_data_900_sys);


TLegend *leg = new TLegend(0.58,0.60,0.90,0.85);
leg->SetTextSize(0.042);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(ks_dNdpT_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
leg->AddEntry(ks_dNdpT_cor_d6t_7000, "     D6T MC" , "l" );
leg->AddEntry(ks_dNdpT_cor_py8_7000, "     Pythia8 MC" , "l" );
leg->AddEntry(ks_dNdpT_cor_p0_7000, "     P0 MC" , "l" );
leg->AddEntry(ks_dNdpT_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );
leg->AddEntry(ks_dNdpT_cor_d6t_900, "     D6T MC" , "l" );


TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

ks_dNdpT_cor_data_7000_sys->Draw();
ks_dNdpT_cor_d6t_7000->Draw("samehist");
//tsfunc7000d6t->Draw("same");
ks_dNdpT_cor_py8_7000->Draw("samehist");
//tsfunc7000py8->Draw("same");
ks_dNdpT_cor_p0_7000->Draw("samehist");
//tsfunc7000p0->Draw("same");
ks_dNdpT_cor_d6t_900->Draw("samehist");
//tsfunc900d6t->Draw("same");
ks_dNdpT_cor_data_7000->Draw("E1same");
//tsfunc7000->Draw("same"); 
ks_dNdpT_cor_data_900_sys->Draw("same");
ks_dNdpT_cor_data_900->Draw("E1same");
//tsfunc900->Draw("same"); 
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

c1->SaveAs(outFile_pt_evt + ".root");
c1->SaveAs(outFile_pt_evt + ".eps");
c1->SaveAs(outFile_pt_evt + ".C");
c1->SaveAs(outFile_pt_evt + ".png");

delete c1;


TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

ks_dNdpT_cor_data_7000_sys->SetMaximum(3.0);
ks_dNdpT_cor_data_7000_sys->SetMinimum(0.000004);
ks_dNdpT_cor_data_7000_sys->Draw();
ks_dNdpT_cor_d6t_7000->Draw("samehist");
//tsfunc7000d6t->Draw("same");
ks_dNdpT_cor_py8_7000->Draw("samehist");
//tsfunc7000py8->Draw("same");
ks_dNdpT_cor_p0_7000->Draw("samehist");
//tsfunc7000p0->Draw("same");
ks_dNdpT_cor_d6t_900->Draw("samehist");
//tsfunc900d6t->Draw("same");
ks_dNdpT_cor_data_7000->Draw("E1same");
//tsfunc7000->Draw("same"); 
ks_dNdpT_cor_data_900_sys->Draw("same");
ks_dNdpT_cor_data_900->Draw("E1same");
//tsfunc900->Draw("same"); 
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

gPad->SetLogy(kTRUE);
c1->SaveAs(outFile_pt_evt + "_log.root");
c1->SaveAs(outFile_pt_evt + "_log.eps");
c1->SaveAs(outFile_pt_evt + "_log.C");
c1->SaveAs(outFile_pt_evt + "_log.png");
gPad->SetLogy(kFALSE);

delete c1;



TLegend *leg = new TLegend(0.58,0.70,0.90,0.85);
leg->SetTextSize(0.042);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(ks_dNdpT_data_d6t_7000, " 7 TeV D6T MC" , "l" );
leg->AddEntry(ks_dNdpT_data_py8_7000, " 7 TeV Pythia8 MC" , "l" );
leg->AddEntry(ks_dNdpT_data_p0_7000, "  7 TeV P0 MC" , "l" );
leg->AddEntry(ks_dNdpT_data_d6t_900, "  900 GeV D6T MC" , "l" );


TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

ks_dNdpT_data_d6t_7000->GetYaxis()->SetTitle("Fully corrected  MC / Data");
ks_dNdpT_data_d6t_7000->SetMaximum(1.6);
ks_dNdpT_data_d6t_7000->SetMinimum(0.2);
ks_dNdpT_data_d6t_7000->Draw();
ks_dNdpT_data_py8_7000->Draw("same");
ks_dNdpT_data_p0_7000->Draw("same");
ks_dNdpT_data_d6t_900->Draw("same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

c1->SaveAs(outFile_pt_evt + "_ratio.root");
c1->SaveAs(outFile_pt_evt + "_ratio.eps");
c1->SaveAs(outFile_pt_evt + "_ratio.C");
c1->SaveAs(outFile_pt_evt + "_ratio.png");

delete c1;





// y for Kshort

inFile = "plots/ksy_cor_data_pythia8_7TeV_weight102_nobinomial2.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_cor_data_7000 = (TH1F*)c1->GetPrimitive("ks_dNdy");
fin->Close();
delete c1;
delete fin;
double yld_7000_y = ks_dNdy_cor_data_7000->Integral();
ks_dNdy_cor_data_7000->Scale(1.0/(trig7000*ybinwid));
double ks_yld_7000_y = calc_integral(ks_dNdy_cor_data_7000, 1.0);
ks_dNdy_cor_data_7000->GetXaxis()->SetTitle("K_{S}^{0}  rapidity");
ks_dNdy_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_cor_data_7000->SetMaximum(0.4);
ks_dNdy_cor_data_7000->SetMinimum(0.);
ks_dNdy_cor_data_7000->SetLineColor(4);
ks_dNdy_cor_data_7000->SetMarkerColor(4);
ks_dNdy_cor_data_7000->SetMarkerStyle(20);
ks_dNdy_cor_data_7000->SetMarkerSize(0.8);
TH1F *ks_dNdy_cor_data_7000_sys = (TH1F*)ks_dNdy_cor_data_7000->Clone("ks_dNdy_cor_data_7000_sys");
//cout << "ks y 7000 integral / event: " << ks_yld_7000_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_cor_data_7000_sys->GetBinContent(i);
  double error = ks_dNdy_cor_data_7000_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i-1] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdy_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i-1] << endl;
}
ks_dNdy_cor_data_7000_sys->SetError(yent);

inFile = "ks_tuneD6T_7000_new2/plots/ks_dNdy_gen_tuneD6T_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_cor_d6t_7000 = (TH1F*)c1->GetPrimitive("ks_dNdy_gen");
fin->Close();
delete c1;
delete fin;
ks_dNdy_cor_d6t_7000->Scale(1.0/(mcevt_d6t_7000*ybinwid));
double ks_yld_d6t_7000_y = calc_integral(ks_dNdy_cor_d6t_7000, 1.0);
ks_dNdy_cor_d6t_7000->GetXaxis()->SetTitle("K_{S}^{0}  y");
ks_dNdy_cor_d6t_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_cor_d6t_7000->SetLineColor(4);
ks_dNdy_cor_d6t_7000->SetMarkerStyle(1);
ks_dNdy_cor_d6t_7000->SetLineStyle(1);
cout << "ks D6T y 7000 integral / event: " << ks_yld_d6t_7000_y << endl;
TH1F *ks_dNdy_data_d6t_7000 = (TH1F*)ks_dNdy_cor_d6t_7000->Clone("ks_dNdy_data_d6t_7000");
ks_dNdy_data_d6t_7000->Divide(ks_dNdy_cor_data_7000_sys);

inFile = "ks_pythia8_7000_new2/plots/ks_dNdy_gen_pythia8_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_cor_py8_7000 = (TH1F*)c1->GetPrimitive("ks_dNdy_gen");
fin->Close();
delete c1;
delete fin;
ks_dNdy_cor_py8_7000->Scale(1.0/(mcevt_py8_7000*ybinwid));
double ks_yld_py8_7000_y = calc_integral(ks_dNdy_cor_py8_7000, 1.0);
ks_dNdy_cor_py8_7000->GetXaxis()->SetTitle("K_{S}^{0}  y");
ks_dNdy_cor_py8_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_cor_py8_7000->SetLineColor(4);
ks_dNdy_cor_py8_7000->SetMarkerStyle(1);
ks_dNdy_cor_py8_7000->SetLineStyle(2);
cout << "ks PY8 y 7000 integral / event: " << ks_yld_py8_7000_y << endl;
TH1F *ks_dNdy_data_py8_7000 = (TH1F*)ks_dNdy_cor_py8_7000->Clone("ks_dNdy_data_py8_7000");
ks_dNdy_data_py8_7000->Divide(ks_dNdy_cor_data_7000_sys);

inFile = "ks_tuneP0_7000_new2/plots/ks_dNdy_gen_tuneP0_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_cor_p0_7000 = (TH1F*)c1->GetPrimitive("ks_dNdy_gen");
fin->Close();
delete c1;
delete fin;
ks_dNdy_cor_p0_7000->Scale(1.0/(mcevt_p0_7000*ybinwid));
double ks_yld_p0_7000_y = calc_integral(ks_dNdy_cor_p0_7000, 1.0);
ks_dNdy_cor_p0_7000->GetXaxis()->SetTitle("K_{S}^{0}  y");
ks_dNdy_cor_p0_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_cor_p0_7000->SetLineColor(4);
ks_dNdy_cor_p0_7000->SetMarkerStyle(1);
ks_dNdy_cor_p0_7000->SetLineStyle(3);
cout << "ks P0 y 7000 integral / event: " << ks_yld_p0_7000_y << endl;
TH1F *ks_dNdy_data_p0_7000 = (TH1F*)ks_dNdy_cor_p0_7000->Clone("ks_dNdy_data_p0_7000");
ks_dNdy_data_p0_7000->Divide(ks_dNdy_cor_data_7000_sys);



inFile = "plots/ksy_cor_data_tuneD6T_900GeV_weight111_nobinomial2.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_cor_data_900 = (TH1F*)c1->GetPrimitive("ks_dNdy");
fin->Close();
delete c1;
delete fin;
double yld_900_y = ks_dNdy_cor_data_900->Integral();
ks_dNdy_cor_data_900->Scale(1.0/(trig900*ybinwid));
double ks_yld_900_y = calc_integral(ks_dNdy_cor_data_900, 1.0);
ks_dNdy_cor_data_900->GetXaxis()->SetTitle("K_{S}^{0}  rapidity");
ks_dNdy_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_cor_data_900->SetMaximum(0.4);
ks_dNdy_cor_data_900->SetMinimum(0.);
ks_dNdy_cor_data_900->SetLineColor(2);
ks_dNdy_cor_data_900->SetMarkerColor(2);
ks_dNdy_cor_data_900->SetMarkerStyle(21);
ks_dNdy_cor_data_900->SetMarkerSize(0.8);
TH1F *ks_dNdy_cor_data_900_sys = (TH1F*)ks_dNdy_cor_data_900->Clone("ks_dNdy_cor_data_7000_sys");
//cout << "ks y 900 integral / event: " << ks_yld_900_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_cor_data_900_sys->GetBinContent(i);
  double error = ks_dNdy_cor_data_900_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i-1] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdy_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i-1] << endl;
}
ks_dNdy_cor_data_900_sys->SetError(yent);


inFile = "ks_tuneD6T_900_new2/plots/ks_dNdy_gen_tuneD6T_900GeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_cor_d6t_900 = (TH1F*)c1->GetPrimitive("ks_dNdy_gen");
fin->Close();
delete c1;
delete fin;
ks_dNdy_cor_d6t_900->Scale(1.0/(mcevt_d6t_900*ybinwid));
double ks_yld_d6t_900_y = calc_integral(ks_dNdy_cor_d6t_900, 1.0);
ks_dNdy_cor_d6t_900->GetXaxis()->SetTitle("K_{S}^{0}  y");
ks_dNdy_cor_d6t_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_cor_d6t_900->SetLineColor(2);
ks_dNdy_cor_d6t_900->SetMarkerStyle(1);
ks_dNdy_cor_d6t_900->SetLineStyle(1);
cout << "ks D6T y 900 integral / event: " << ks_yld_d6t_900_y << endl;
TH1F *ks_dNdy_data_d6t_900 = (TH1F*)ks_dNdy_cor_d6t_900->Clone("ks_dNdy_data_d6t_900");
ks_dNdy_data_d6t_900->Divide(ks_dNdy_cor_data_900_sys);


TLegend *leg1 = new TLegend(0.20,0.23,0.55,0.38);
leg1->SetTextSize(0.042);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetFillStyle(0);
leg1->AddEntry(ks_dNdpT_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
leg1->AddEntry(ks_dNdpT_cor_d6t_7000, "     D6T MC" , "l" );
leg1->AddEntry(ks_dNdpT_cor_py8_7000, "     Pythia8 MC" , "l" );
leg1->AddEntry(ks_dNdpT_cor_p0_7000, "     P0 MC" , "l" );

TLegend *leg2 = new TLegend(0.65,0.25,0.90,0.325);
leg2->SetTextSize(0.042);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetFillStyle(0);
leg2->AddEntry(ks_dNdpT_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );
leg2->AddEntry(ks_dNdpT_cor_d6t_900, "     D6T MC" , "l" );


TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

ks_dNdy_cor_data_7000_sys->Draw();
ks_dNdy_cor_d6t_7000->Draw("samehist");
ks_dNdy_cor_py8_7000->Draw("samehist");
ks_dNdy_cor_p0_7000->Draw("samehist");
ks_dNdy_cor_d6t_900->Draw("samehist");
ks_dNdy_cor_data_7000->Draw("E1same");
ks_dNdy_cor_data_900_sys->Draw("same");
ks_dNdy_cor_data_900->Draw("E1same");
leg1->Draw("same");
leg2->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

c1->SaveAs(outFile_y_evt + ".root");
c1->SaveAs(outFile_y_evt + ".eps");
c1->SaveAs(outFile_y_evt + ".C");
c1->SaveAs(outFile_y_evt + ".png");

delete c1;



TLegend *leg = new TLegend(0.58,0.20,0.90,0.38);
leg->SetTextSize(0.042);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(ks_dNdy_data_d6t_7000, " 7 TeV D6T MC" , "l" );
leg->AddEntry(ks_dNdy_data_py8_7000, " 7 TeV Pythia8 MC" , "l" );
leg->AddEntry(ks_dNdy_data_p0_7000, "  7 TeV P0 MC" , "l" );
leg->AddEntry(ks_dNdy_data_d6t_900, "  900 GeV D6T MC" , "l" );


TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

ks_dNdy_data_d6t_7000->GetYaxis()->SetTitle("Fully corrected  MC / Data");
ks_dNdy_data_d6t_7000->SetMaximum(1.2);
ks_dNdy_data_d6t_7000->SetMinimum(0.4);
ks_dNdy_data_d6t_7000->Draw();
ks_dNdy_data_py8_7000->Draw("same");
ks_dNdy_data_p0_7000->Draw("same");
ks_dNdy_data_d6t_900->Draw("same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

c1->SaveAs(outFile_y_evt + "_ratio.root");
c1->SaveAs(outFile_y_evt + "_ratio.eps");
c1->SaveAs(outFile_y_evt + "_ratio.C");
c1->SaveAs(outFile_y_evt + "_ratio.png");

delete c1;





// Lambdas


  double lamptarray[] = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0,7.0,8.0,10.0};
  int lamptNbins = sizeof(lamptarray)/sizeof(double) -1;
  vector<Double_t> lamptBinMinima;
  for ( int ndx = 0; ndx <= lamptNbins; ndx++ ) {
    lamptBinMinima.push_back( lamptarray[ndx] );
  }
  Double_t lampt_default_binWidth = 0.1;

  Double_t lamptOrigNbins = lamptNbins - 3;

TF1 *lamtsfunc900 = new TF1("lamtsfunc900",lamtsallis,0.0,10.0,3);
lamtsfunc900->SetParNames("Amplitude","T","n");
lamtsfunc900->SetParameters(1.0.,0.2,7.2); 
TF1 *lamtsfunc900d6t = new TF1("lamtsfunc900d6t",lamtsallis,0.0,10.0,3);
lamtsfunc900d6t->SetParNames("Amplitude","T","n");
lamtsfunc900d6t->SetParameters(1.0.,0.2,7.2); 
TF1 *lamtsfunc7000 = new TF1("lamtsfunc7000",lamtsallis,0.0,10.0,3);
lamtsfunc7000->SetParNames("Amplitude","T","n");
lamtsfunc7000->SetParameters(1.0.,0.2,7.2); 
TF1 *lamtsfunc7000d6t = new TF1("lamtsfunc7000d6t",lamtsallis,0.0,10.0,3);
lamtsfunc7000d6t->SetParNames("Amplitude","T","n");
lamtsfunc7000d6t->SetParameters(1.0.,0.2,7.2); 
TF1 *lamtsfunc7000py8 = new TF1("lamtsfunc7000py8",lamtsallis,0.0,10.0,3);
lamtsfunc7000py8->SetParNames("Amplitude","T","n");
lamtsfunc7000py8->SetParameters(1.0.,0.2,7.2); 
TF1 *lamtsfunc7000p0 = new TF1("lamtsfunc7000p0",lamtsallis,0.0,10.0,3);
lamtsfunc7000p0->SetParNames("Amplitude","T","n");
lamtsfunc7000p0->SetParameters(1.0.,0.2,7.2); 

TF1 *lamtsmean900 = new TF1("lamtsmean900",lamtsallismean,0.0,10.0,3);
lamtsmean900->SetParNames("Amplitude","T","n");
lamtsmean900->SetParameters(1.0.,0.9,7.2); 
TF1 *lamtsmean900d6t = new TF1("lamtsmean900d6t",lamtsallismean,0.0,10.0,3);
lamtsmean900d6t->SetParNames("Amplitude","T","n");
lamtsmean900d6t->SetParameters(1.0.,0.9,7.2); 
TF1 *lamtsmean7000 = new TF1("lamtsmean7000",lamtsallismean,0.0,10.0,3);
lamtsmean7000->SetParNames("Amplitude","T","n");
lamtsmean7000->SetParameters(1.0.,0.9,7.2); 
TF1 *lamtsmean7000d6t = new TF1("lamtsmean7000d6t",lamtsallismean,0.0,10.0,3);
lamtsmean7000d6t->SetParNames("Amplitude","T","n");
lamtsmean7000d6t->SetParameters(1.0.,0.9,7.2); 
TF1 *lamtsmean7000py8 = new TF1("lamtsmean7000py8",lamtsallismean,0.0,10.0,3);
lamtsmean7000py8->SetParNames("Amplitude","T","n");
lamtsmean7000py8->SetParameters(1.0.,0.9,7.2); 
TF1 *lamtsmean7000p0 = new TF1("lamtsmean7000p0",lamtsallismean,0.0,10.0,3);
lamtsmean7000p0->SetParNames("Amplitude","T","n");
lamtsmean7000p0->SetParameters(1.0.,0.9,7.2); 

  //TString inFile;
  outFile_pt_evt = "final_plots/lampt_cor_datamc_perevt";

  Double_t ptent[40];

  TString outFile_y_evt = "final_plots/lamy_cor_datamc_perevt";

  Double_t yent[23];


  // Pt for Lambda per event

  inFile = "plots/lampt_cor_data_pythia8_7TeV_weight102_subgeant2_nobinomial2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lam_dNdpT_cor_data_7000 = (TH1F*)c1->GetPrimitive("lam_dNdpT_scale");
  fin->Close();
  delete c1;
  delete fin;
  lam_dNdpT_cor_data_7000->Scale(1.0/(trig7000*ptbinwid));
  double lam_yld_7000_pt = calc_integral(lam_dNdpT_cor_data_7000, 1.0);
  lam_dNdpT_cor_data_7000->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  lam_dNdpT_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  lam_dNdpT_cor_data_7000->SetLineColor(4);
  lam_dNdpT_cor_data_7000->SetMarkerColor(4);
  lam_dNdpT_cor_data_7000->SetMarkerStyle(20);
  lam_dNdpT_cor_data_7000->SetMarkerSize(0.8);
  TH1F *lam_dNdpT_cor_data_7000_sys = (TH1F*)lam_dNdpT_cor_data_7000->Clone("lam_dNdpT_cor_data_7000_sys");
  //cout << "lam pt 7000 integral / event: " << lam_yld_7000_pt << endl;
  for (int i = 1; i <= lam_dNdpT_cor_data_7000_sys->GetXaxis()->GetNbins(); ++i) {
    double value = lam_dNdpT_cor_data_7000_sys->GetBinContent(i);
    double error = lam_dNdpT_cor_data_7000_sys->GetBinError(i);
    double syserr = 0.03*value;
    ptent[i-1] = sqrt(error*error + syserr*syserr);
    //cout << lam_dNdpT_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i-1] << endl;
  }
  lam_dNdpT_cor_data_7000_sys->SetError(ptent);
  lam_dNdpT_cor_data_7000->Fit("lamtsmean7000","EM0");
  lam_dNdpT_cor_data_7000->Fit("lamtsfunc7000","EM0");

  // Change the histogram to have empty bins below 0.3.  Will NOT have the TF1 in it.
  //TH1F* oldHisto;
  //rebinpT( lam_dNdpT_cor_data_7000, oldHisto, 4, 20 );
  TH1F* tmp7 = new TH1F("lam_dNdpT_scale7000", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  TH1F* tmp7sys = new TH1F("lam_dNdpT_scale7000_sys", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  for( int ndx = 1; ndx <= lamptNbins; ndx++ ) {
    if( ndx == 1 || ndx == 2 || ndx == 3) continue;
    int origNdx = ndx - 3;
    tmp7->SetBinContent(ndx, lam_dNdpT_cor_data_7000->GetBinContent(origNdx));
    tmp7->SetBinError(ndx, lam_dNdpT_cor_data_7000->GetBinError(origNdx));
    tmp7sys->SetBinContent(ndx, lam_dNdpT_cor_data_7000_sys->GetBinContent(origNdx));
    tmp7sys->SetBinError(ndx, lam_dNdpT_cor_data_7000_sys->GetBinError(origNdx));
  }
  tmp7->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmp7->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmp7->SetLineColor(4);
  tmp7->SetMarkerColor(4);
  tmp7->SetMarkerStyle(20);
  tmp7->SetMarkerSize(0.8);

  tmp7sys->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmp7sys->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmp7sys->SetLineColor(4);
  tmp7sys->SetMarkerColor(4);
  tmp7sys->SetMarkerStyle(20);
  tmp7sys->SetMarkerSize(0.8);

  TH1F* switch7000 = lam_dNdpT_cor_data_7000;
  lam_dNdpT_cor_data_7000 = tmp7;
  tmp7 = switch7000;
  TH1F* switch7000sys = lam_dNdpT_cor_data_7000_sys;
  lam_dNdpT_cor_data_7000_sys = tmp7sys;
  tmp7sys = switch7000;

  lamtsfunc7000->SetNpx(2000);
  lamtsfunc7000->SetLineColor(4);

  inFile = "plots/lampt_cor_data_tuneD6T_900GeV_weight111_subgeant2_nobinomial2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lam_dNdpT_cor_data_900 = (TH1F*)c1->GetPrimitive("lam_dNdpT_scale");
  fin->Close();
  delete c1;
  delete fin;
  double yld_900_pt = lam_dNdpT_cor_data_900->Integral();
  lam_dNdpT_cor_data_900->Scale(1.0/(trig900*ptbinwid));
  double lam_yld_900_pt = calc_integral(lam_dNdpT_cor_data_900, 1.0);
  lam_dNdpT_cor_data_900->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  lam_dNdpT_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  lam_dNdpT_cor_data_900->SetLineColor(2);
  lam_dNdpT_cor_data_900->SetMarkerColor(2);
  lam_dNdpT_cor_data_900->SetMarkerStyle(21);
  lam_dNdpT_cor_data_900->SetMarkerSize(0.8);
  TH1F *lam_dNdpT_cor_data_900_sys = (TH1F*)lam_dNdpT_cor_data_900->Clone("lam_dNdpT_cor_data_7000_sys");
  //cout << "lam pt 900 integral / event: " << lam_yld_900_pt << endl;
  for (int i = 1; i <= lam_dNdpT_cor_data_900_sys->GetXaxis()->GetNbins(); ++i) {
    double value = lam_dNdpT_cor_data_900_sys->GetBinContent(i);
    double error = lam_dNdpT_cor_data_900_sys->GetBinError(i);
    double syserr = 0.03*value;
    ptent[i-1] = sqrt(error*error + syserr*syserr);
    //cout << lam_dNdpT_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i-1] << endl;
  }
  lam_dNdpT_cor_data_900_sys->SetError(ptent);
  lam_dNdpT_cor_data_900->Fit("lamtsmean900","EM0");
  lam_dNdpT_cor_data_900->Fit("lamtsfunc900","EM0");

  // Change the histogram to have empty bins below 0.3.  Will NOT have the TF1 in it.
  TH1F* tmp9 = new TH1F("lam_dNdpT_scale900", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  TH1F* tmp9sys = new TH1F("lam_dNdpT_scale900sys", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  for( int ndx = 1; ndx <= lamptNbins; ndx++ ) {
    if( ndx == 1 || ndx == 2 || ndx == 3) continue;
    int origNdx = ndx - 3;
    tmp9->SetBinContent(ndx, lam_dNdpT_cor_data_900->GetBinContent(origNdx));
    tmp9->SetBinError(ndx, lam_dNdpT_cor_data_900->GetBinError(origNdx));
    tmp9sys->SetBinContent(ndx, lam_dNdpT_cor_data_900_sys->GetBinContent(origNdx));
    tmp9sys->SetBinError(ndx, lam_dNdpT_cor_data_900_sys->GetBinError(origNdx));
  }
  tmp9->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmp9->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmp9->SetLineColor(2);
  tmp9->SetMarkerColor(2);
  tmp9->SetMarkerStyle(21);
  tmp9->SetMarkerSize(0.8);

  tmp9sys->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmp9sys->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmp9sys->SetLineColor(2);
  tmp9sys->SetMarkerColor(2);
  tmp9sys->SetMarkerStyle(21);
  tmp9sys->SetMarkerSize(0.8);

  TH1F* switch900 = lam_dNdpT_cor_data_900;
  lam_dNdpT_cor_data_900 = tmp9;
  tmp9 = switch900;

  TH1F* switch900sys = lam_dNdpT_cor_data_900_sys;
  lam_dNdpT_cor_data_900_sys = tmp9sys;
  tmp9sys = switch900sys;
  //rebinpT( lam_dNdpT_cor_data_900, oldHisto, 2, 21 );

  lamtsfunc900->SetNpx(2000);
  lamtsfunc900->SetLineColor(2);

inFile = "lam_tuneD6T_7000_new2/plots/lam_dNdpT_gen_tuneD6T_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* lam_dNdpT_cor_d6t_7000 = (TH1F*)c1->GetPrimitive("lam_dNdpT_gen_scale");
fin->Close();
delete c1;
delete fin;
lam_dNdpT_cor_d6t_7000->Scale(1.0/(mcevt_d6t_7000*ptbinwid));
double lam_yld_d6t_7000_pt = calc_integral(lam_dNdpT_cor_d6t_7000, 1.0);
lam_dNdpT_cor_d6t_7000->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
lam_dNdpT_cor_d6t_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
lam_dNdpT_cor_d6t_7000->SetLineColor(4);
lam_dNdpT_cor_d6t_7000->SetMarkerStyle(1);
lam_dNdpT_cor_d6t_7000->SetLineStyle(1);
cout << "lam D6T pt 7000 integral / event: " << lam_yld_d6t_7000_pt << endl;
lam_dNdpT_cor_d6t_7000->Fit("lamtsmean7000d6t","EM0");
lam_dNdpT_cor_d6t_7000->Fit("lamtsfunc7000d6t","EM0");
lamtsfunc7000d6t->SetNpx(2000);
lamtsfunc7000d6t->SetLineColor(4);
lamtsfunc7000d6t->SetLineStyle(1);
TH1F *lam_dNdpT_data_d6t_7000 = (TH1F*)lam_dNdpT_cor_d6t_7000->Clone("lam_dNdpT_data_d6t_7000");
lam_dNdpT_data_d6t_7000->Divide(lam_dNdpT_cor_data_7000_sys);

inFile = "lam_pythia8_7000_new2/plots/lam_dNdpT_gen_pythia8_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* lam_dNdpT_cor_py8_7000 = (TH1F*)c1->GetPrimitive("lam_dNdpT_gen_scale");
fin->Close();
delete c1;
delete fin;
lam_dNdpT_cor_py8_7000->Scale(1.0/(mcevt_py8_7000*ptbinwid));
double lam_yld_py8_7000_pt = calc_integral(lam_dNdpT_cor_py8_7000, 1.0);
lam_dNdpT_cor_py8_7000->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
lam_dNdpT_cor_py8_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
lam_dNdpT_cor_py8_7000->SetLineColor(4);
lam_dNdpT_cor_py8_7000->SetMarkerStyle(1);
lam_dNdpT_cor_py8_7000->SetLineStyle(2);
cout << "lam PY8 pt 7000 integral / event: " << lam_yld_py8_7000_pt << endl;
lam_dNdpT_cor_py8_7000->Fit("lamtsmean7000py8","EM0");
lam_dNdpT_cor_py8_7000->Fit("lamtsfunc7000py8","EM0");
lamtsfunc7000py8->SetNpx(2000);
lamtsfunc7000py8->SetLineColor(4);
lamtsfunc7000py8->SetLineStyle(2);
TH1F *lam_dNdpT_data_py8_7000 = (TH1F*)lam_dNdpT_cor_py8_7000->Clone("lam_dNdpT_data_py8_7000");
lam_dNdpT_data_py8_7000->Divide(lam_dNdpT_cor_data_7000_sys);

inFile = "lam_tuneP0_7000_new2/plots/lam_dNdpT_gen_tuneP0_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* lam_dNdpT_cor_p0_7000 = (TH1F*)c1->GetPrimitive("lam_dNdpT_gen_scale");
fin->Close();
delete c1;
delete fin;
lam_dNdpT_cor_p0_7000->Scale(1.0/(mcevt_p0_7000*ptbinwid));
double lam_yld_p0_7000_pt = calc_integral(lam_dNdpT_cor_p0_7000, 1.0);
lam_dNdpT_cor_p0_7000->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
lam_dNdpT_cor_p0_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
lam_dNdpT_cor_p0_7000->SetLineColor(4);
lam_dNdpT_cor_p0_7000->SetMarkerStyle(1);
lam_dNdpT_cor_p0_7000->SetLineStyle(3);
cout << "lam P0 pt 7000 integral / event: " << lam_yld_p0_7000_pt << endl;
lam_dNdpT_cor_p0_7000->Fit("lamtsmean7000p0","EM0");
lam_dNdpT_cor_p0_7000->Fit("lamtsfunc7000p0","EM0");
lamtsfunc7000p0->SetNpx(2000);
lamtsfunc7000p0->SetLineColor(4);
lamtsfunc7000p0->SetLineStyle(3);
TH1F *lam_dNdpT_data_p0_7000 = (TH1F*)lam_dNdpT_cor_p0_7000->Clone("lam_dNdpT_data_p0_7000");
lam_dNdpT_data_p0_7000->Divide(lam_dNdpT_cor_data_7000_sys);

inFile = "lam_tuneD6T_900_new2/plots/lam_dNdpT_gen_tuneD6T_900GeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* lam_dNdpT_cor_d6t_900 = (TH1F*)c1->GetPrimitive("lam_dNdpT_gen_scale");
fin->Close();
delete c1;
delete fin;
lam_dNdpT_cor_d6t_900->Scale(1.0/(mcevt_d6t_900*ptbinwid));
double lam_yld_d6t_900_pt = calc_integral(lam_dNdpT_cor_d6t_900, 1.0);
lam_dNdpT_cor_d6t_900->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
lam_dNdpT_cor_d6t_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
lam_dNdpT_cor_d6t_900->SetLineColor(2);
lam_dNdpT_cor_d6t_900->SetMarkerStyle(1);
lam_dNdpT_cor_d6t_900->SetLineStyle(1);
cout << "lam D6T pt 900 integral / event: " << lam_yld_d6t_900_pt << endl;
lam_dNdpT_cor_d6t_900->Fit("lamtsmean900d6t","EM0");
lam_dNdpT_cor_d6t_900->Fit("lamtsfunc900d6t","EM0");
lamtsfunc900d6t->SetNpx(2000);
lamtsfunc900d6t->SetLineColor(1);
lamtsfunc900d6t->SetLineStyle(2);
TH1F *lam_dNdpT_data_d6t_900 = (TH1F*)lam_dNdpT_cor_d6t_900->Clone("lam_dNdpT_data_d6t_900");
lam_dNdpT_data_d6t_900->Divide(lam_dNdpT_cor_data_900_sys);




  TLegend *leg = new TLegend(0.58,0.60,0.90,0.85);
  leg->SetTextSize(0.042);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
leg->AddEntry(lam_dNdpT_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
leg->AddEntry(lam_dNdpT_cor_d6t_7000, "     D6T MC" , "l" );
leg->AddEntry(lam_dNdpT_cor_py8_7000, "     Pythia8 MC" , "l" );
leg->AddEntry(lam_dNdpT_cor_p0_7000, "     P0 MC" , "l" );
leg->AddEntry(lam_dNdpT_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );
leg->AddEntry(lam_dNdpT_cor_d6t_900, "     D6T MC" , "l" );


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lam_dNdpT_cor_data_7000->Draw("E1");
  //oldHisto->Draw("E1");

  //c1->SaveAs("test.png");
  lam_dNdpT_cor_data_7000_sys->Draw("same");
lam_dNdpT_cor_d6t_7000->Draw("samehist");
lam_dNdpT_cor_py8_7000->Draw("samehist");
lam_dNdpT_cor_p0_7000->Draw("samehist");
lam_dNdpT_cor_d6t_900->Draw("samehist");
  lam_dNdpT_cor_data_7000->Draw("E1same");
//  lamtsfunc7000->Draw("same"); 
  lam_dNdpT_cor_data_900_sys->Draw("same");
  lam_dNdpT_cor_data_900->Draw("E1same");
//  lamtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

  c1->SaveAs(outFile_pt_evt + ".root");
  c1->SaveAs(outFile_pt_evt + ".eps");
  c1->SaveAs(outFile_pt_evt + ".C");
  c1->SaveAs(outFile_pt_evt + ".png");

  delete c1;


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lam_dNdpT_cor_data_7000->SetMaximum(1.0);
  lam_dNdpT_cor_data_7000->SetMinimum(0.000001);
  lam_dNdpT_cor_data_7000->Draw("E1");
lam_dNdpT_cor_d6t_7000->Draw("samehist");
lam_dNdpT_cor_py8_7000->Draw("samehist");
lam_dNdpT_cor_p0_7000->Draw("samehist");
lam_dNdpT_cor_d6t_900->Draw("samehist");
  lam_dNdpT_cor_data_7000_sys->Draw("same");
  lam_dNdpT_cor_data_7000->Draw("E1same");
//  lamtsfunc7000->Draw("same"); 
  lam_dNdpT_cor_data_900_sys->Draw("same");
  lam_dNdpT_cor_data_900->Draw("E1same");
//  lamtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

  gPad->SetLogy(kTRUE);
  c1->SaveAs(outFile_pt_evt + "_log.root");
  c1->SaveAs(outFile_pt_evt + "_log.eps");
  c1->SaveAs(outFile_pt_evt + "_log.C");
  c1->SaveAs(outFile_pt_evt + "_log.png");
  gPad->SetLogy(kFALSE);

  delete c1;



TLegend *leg = new TLegend(0.18,0.78,0.50,0.93);
leg->SetTextSize(0.042);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(lam_dNdpT_data_d6t_7000, " 7 TeV D6T MC" , "l" );
leg->AddEntry(lam_dNdpT_data_py8_7000, " 7 TeV Pythia8 MC" , "l" );
leg->AddEntry(lam_dNdpT_data_p0_7000, "  7 TeV P0 MC" , "l" );
leg->AddEntry(lam_dNdpT_data_d6t_900, "  900 GeV D6T MC" , "l" );


TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

lam_dNdpT_data_d6t_7000->GetYaxis()->SetTitle("Fully corrected  MC / Data");
lam_dNdpT_data_d6t_7000->SetMaximum(1.9);
lam_dNdpT_data_d6t_7000->SetMinimum(0.2);
lam_dNdpT_data_d6t_7000->Draw();
lam_dNdpT_data_py8_7000->Draw("same");
lam_dNdpT_data_p0_7000->Draw("same");
lam_dNdpT_data_d6t_900->Draw("same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

c1->SaveAs(outFile_pt_evt + "_ratio.root");
c1->SaveAs(outFile_pt_evt + "_ratio.eps");
c1->SaveAs(outFile_pt_evt + "_ratio.C");
c1->SaveAs(outFile_pt_evt + "_ratio.png");

delete c1;


  //return;

  // Switch histograms again
  switch7000 = lam_dNdpT_cor_data_7000;
  lam_dNdpT_cor_data_7000 = tmp7;
  tmp7 = switch7000;

  switch900 = lam_dNdpT_cor_data_900;
  lam_dNdpT_cor_data_900 = tmp9;
  tmp9 = switch900;






  // y for Lambda

  inFile = "plots/lamy_cor_data_pythia8_7TeV_weight102_subgeant2_nobinomial2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lam_dNdy_cor_data_7000 = (TH1F*)c1->GetPrimitive("lam_dNdy");
  fin->Close();
  delete c1;
  delete fin;
  double yld_7000_y = lam_dNdy_cor_data_7000->Integral();
  lam_dNdy_cor_data_7000->Scale(1.0/(trig7000*ybinwid));
  double lam_yld_7000_y = calc_integral(lam_dNdy_cor_data_7000, 1.0);
  lam_dNdy_cor_data_7000->GetXaxis()->SetTitle("#Lambda^{0}  rapidity");
  lam_dNdy_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
  lam_dNdy_cor_data_7000->SetMaximum(0.22);
  lam_dNdy_cor_data_7000->SetMinimum(0.);
  lam_dNdy_cor_data_7000->SetLineColor(4);
  lam_dNdy_cor_data_7000->SetMarkerColor(4);
  lam_dNdy_cor_data_7000->SetMarkerStyle(20);
  lam_dNdy_cor_data_7000->SetMarkerSize(0.8);
  TH1F *lam_dNdy_cor_data_7000_sys = (TH1F*)lam_dNdy_cor_data_7000->Clone("lam_dNdy_cor_data_7000_sys");
  //cout << "lam y 7000 integral / event: " << lam_yld_7000_y << endl;
  for (int i = 1; i <= 23; ++i) {
    double value = lam_dNdy_cor_data_7000_sys->GetBinContent(i);
    double error = lam_dNdy_cor_data_7000_sys->GetBinError(i);
    double syserr = 0.03*value;
    yent[i-1] = sqrt(error*error + syserr*syserr);
    //cout << lam_dNdy_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i-1] << endl;
  }
  lam_dNdy_cor_data_7000_sys->SetError(yent);

inFile = "lam_tuneD6T_7000_new2/plots/lam_dNdy_gen_tuneD6T_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* lam_dNdy_cor_d6t_7000 = (TH1F*)c1->GetPrimitive("lam_dNdy_gen");
fin->Close();
delete c1;
delete fin;
lam_dNdy_cor_d6t_7000->Scale(1.0/(mcevt_d6t_7000*ybinwid));
double lam_yld_d6t_7000_y = calc_integral(lam_dNdy_cor_d6t_7000, 1.0);
lam_dNdy_cor_d6t_7000->GetXaxis()->SetTitle("#Lambda^{0}  y");
lam_dNdy_cor_d6t_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
lam_dNdy_cor_d6t_7000->SetLineColor(4);
lam_dNdy_cor_d6t_7000->SetMarkerStyle(1);
lam_dNdy_cor_d6t_7000->SetLineStyle(1);
cout << "lam D6T y 7000 integral / event: " << lam_yld_d6t_7000_y << endl;
TH1F *lam_dNdy_data_d6t_7000 = (TH1F*)lam_dNdy_cor_d6t_7000->Clone("lam_dNdy_data_d6t_7000");
lam_dNdy_data_d6t_7000->Divide(lam_dNdy_cor_data_7000_sys);

inFile = "lam_pythia8_7000_new2/plots/lam_dNdy_gen_pythia8_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* lam_dNdy_cor_py8_7000 = (TH1F*)c1->GetPrimitive("lam_dNdy_gen");
fin->Close();
delete c1;
delete fin;
lam_dNdy_cor_py8_7000->Scale(1.0/(mcevt_py8_7000*ybinwid));
double lam_yld_py8_7000_y = calc_integral(lam_dNdy_cor_py8_7000, 1.0);
lam_dNdy_cor_py8_7000->GetXaxis()->SetTitle("#Lambda^{0}  y");
lam_dNdy_cor_py8_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
lam_dNdy_cor_py8_7000->SetLineColor(4);
lam_dNdy_cor_py8_7000->SetMarkerStyle(1);
lam_dNdy_cor_py8_7000->SetLineStyle(2);
cout << "lam PY8 y 7000 integral / event: " << lam_yld_py8_7000_y << endl;
TH1F *lam_dNdy_data_py8_7000 = (TH1F*)lam_dNdy_cor_py8_7000->Clone("lam_dNdy_data_py8_7000");
lam_dNdy_data_py8_7000->Divide(lam_dNdy_cor_data_7000_sys);

inFile = "lam_tuneP0_7000_new2/plots/lam_dNdy_gen_tuneP0_7TeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* lam_dNdy_cor_p0_7000 = (TH1F*)c1->GetPrimitive("lam_dNdy_gen");
fin->Close();
delete c1;
delete fin;
lam_dNdy_cor_p0_7000->Scale(1.0/(mcevt_p0_7000*ybinwid));
double lam_yld_p0_7000_y = calc_integral(lam_dNdy_cor_p0_7000, 1.0);
lam_dNdy_cor_p0_7000->GetXaxis()->SetTitle("#Lambda^{0}  y");
lam_dNdy_cor_p0_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
lam_dNdy_cor_p0_7000->SetLineColor(4);
lam_dNdy_cor_p0_7000->SetMarkerStyle(1);
lam_dNdy_cor_p0_7000->SetLineStyle(3);
cout << "lam P0 y 7000 integral / event: " << lam_yld_p0_7000_y << endl;
TH1F *lam_dNdy_data_p0_7000 = (TH1F*)lam_dNdy_cor_p0_7000->Clone("lam_dNdy_data_p0_7000");
lam_dNdy_data_p0_7000->Divide(lam_dNdy_cor_data_7000_sys);



  inFile = "plots/lamy_cor_data_tuneD6T_900GeV_weight111_subgeant2_nobinomial2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lam_dNdy_cor_data_900 = (TH1F*)c1->GetPrimitive("lam_dNdy");
  fin->Close();
  delete c1;
  delete fin;
  double yld_900_y = lam_dNdy_cor_data_900->Integral();
  lam_dNdy_cor_data_900->Scale(1.0/(trig900*ybinwid));
  double lam_yld_900_y = calc_integral(lam_dNdy_cor_data_900, 1.0);
  lam_dNdy_cor_data_900->GetXaxis()->SetTitle("#Lambda^{0}  rapidity");
  lam_dNdy_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
  lam_dNdy_cor_data_900->SetMaximum(0.22);
  lam_dNdy_cor_data_900->SetMinimum(0.);
  lam_dNdy_cor_data_900->SetLineColor(2);
  lam_dNdy_cor_data_900->SetMarkerColor(2);
  lam_dNdy_cor_data_900->SetMarkerStyle(21);
  lam_dNdy_cor_data_900->SetMarkerSize(0.8);
  TH1F *lam_dNdy_cor_data_900_sys = (TH1F*)lam_dNdy_cor_data_900->Clone("lam_dNdy_cor_data_7000_sys");
  //cout << "lam y 900 integral / event: " << lam_yld_900_y << endl;
  for (int i = 1; i <= 23; ++i) {
    double value = lam_dNdy_cor_data_900_sys->GetBinContent(i);
    double error = lam_dNdy_cor_data_900_sys->GetBinError(i);
    double syserr = 0.03*value;
    yent[i-1] = sqrt(error*error + syserr*syserr);
    //cout << lam_dNdy_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i-1] << endl;
  }
  lam_dNdy_cor_data_900_sys->SetError(yent);

inFile = "lam_tuneD6T_900_new2/plots/lam_dNdy_gen_tuneD6T_900GeV_weight0.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* lam_dNdy_cor_d6t_900 = (TH1F*)c1->GetPrimitive("lam_dNdy_gen");
fin->Close();
delete c1;
delete fin;
lam_dNdy_cor_d6t_900->Scale(1.0/(mcevt_d6t_900*ybinwid));
double lam_yld_d6t_900_y = calc_integral(lam_dNdy_cor_d6t_900, 1.0);
lam_dNdy_cor_d6t_900->GetXaxis()->SetTitle("#Lambda^{0}  y");
lam_dNdy_cor_d6t_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
lam_dNdy_cor_d6t_900->SetLineColor(2);
lam_dNdy_cor_d6t_900->SetMarkerStyle(1);
lam_dNdy_cor_d6t_900->SetLineStyle(1);
cout << "lam D6T y 900 integral / event: " << lam_yld_d6t_900_y << endl;
TH1F *lam_dNdy_data_d6t_900 = (TH1F*)lam_dNdy_cor_d6t_900->Clone("lam_dNdy_data_d6t_900");
lam_dNdy_data_d6t_900->Divide(lam_dNdy_cor_data_900_sys);


TLegend *leg1 = new TLegend(0.20,0.23,0.55,0.38);
leg1->SetTextSize(0.042);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetFillStyle(0);
leg1->AddEntry(lam_dNdpT_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
leg1->AddEntry(lam_dNdpT_cor_d6t_7000, "     D6T MC" , "l" );
leg1->AddEntry(lam_dNdpT_cor_py8_7000, "     Pythia8 MC" , "l" );
leg1->AddEntry(lam_dNdpT_cor_p0_7000, "     P0 MC" , "l" );

TLegend *leg2 = new TLegend(0.65,0.25,0.90,0.325);
leg2->SetTextSize(0.042);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetFillStyle(0);
leg2->AddEntry(lam_dNdpT_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );
leg2->AddEntry(lam_dNdpT_cor_d6t_900, "     D6T MC" , "l" );


TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

lam_dNdy_cor_data_7000_sys->Draw();
lam_dNdy_cor_d6t_7000->Draw("samehist");
lam_dNdy_cor_py8_7000->Draw("samehist");
lam_dNdy_cor_p0_7000->Draw("samehist");
lam_dNdy_cor_d6t_900->Draw("samehist");
lam_dNdy_cor_data_7000->Draw("E1same");
lam_dNdy_cor_data_900_sys->Draw("same");
lam_dNdy_cor_data_900->Draw("E1same");
leg1->Draw("same");
leg2->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

c1->SaveAs(outFile_y_evt + ".root");
c1->SaveAs(outFile_y_evt + ".eps");
c1->SaveAs(outFile_y_evt + ".C");
c1->SaveAs(outFile_y_evt + ".png");

delete c1;



TLegend *leg = new TLegend(0.58,0.20,0.90,0.38);
leg->SetTextSize(0.042);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(lam_dNdy_data_d6t_7000, " 7 TeV D6T MC" , "l" );
leg->AddEntry(lam_dNdy_data_py8_7000, " 7 TeV Pythia8 MC" , "l" );
leg->AddEntry(lam_dNdy_data_p0_7000, "  7 TeV P0 MC" , "l" );
leg->AddEntry(lam_dNdy_data_d6t_900, "  900 GeV D6T MC" , "l" );


TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

lam_dNdy_data_d6t_7000->GetYaxis()->SetTitle("Fully corrected  MC / Data");
lam_dNdy_data_d6t_7000->SetMaximum(1.0);
lam_dNdy_data_d6t_7000->SetMinimum(0.0);
lam_dNdy_data_d6t_7000->Draw();
lam_dNdy_data_py8_7000->Draw("same");
lam_dNdy_data_p0_7000->Draw("same");
lam_dNdy_data_d6t_900->Draw("same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

c1->SaveAs(outFile_y_evt + "_ratio.root");
c1->SaveAs(outFile_y_evt + "_ratio.eps");
c1->SaveAs(outFile_y_evt + "_ratio.C");
c1->SaveAs(outFile_y_evt + "_ratio.png");

delete c1;







  // Do ratio plots
  TH1F* lamOverKs_pT_cor_data_7000_sys = (TH1F*)lam_dNdpT_cor_data_7000_sys->Clone("lamOverKs_pT_cor_data_7000_sys");
  TH1F* lamOverKs_pT_cor_data_900_sys = (TH1F*)lam_dNdpT_cor_data_900_sys->Clone("lamOverKs_pT_cor_data_7000_sys");

  lamOverKs_pT_cor_data_7000_sys->Divide(ks_dNdpT_cor_data_7000_sys);
  lamOverKs_pT_cor_data_900_sys->Divide(ks_dNdpT_cor_data_900_sys);

  lamOverKs_pT_cor_data_7000_sys->SetMarkerStyle(20);
  lamOverKs_pT_cor_data_7000_sys->SetMarkerColor(4);
  lamOverKs_pT_cor_data_7000_sys->SetLineColor(4);

  lamOverKs_pT_cor_data_900_sys->SetMarkerStyle(21);
  lamOverKs_pT_cor_data_900_sys->SetMarkerColor(2);
  lamOverKs_pT_cor_data_900_sys->SetLineColor(2);

  lamOverKs_pT_cor_data_7000_sys->SetYTitle("N(#Lambda^{0}) / N(K_{S}^{0})");
  lamOverKs_pT_cor_data_7000_sys->SetXTitle("p_{T} [GeV/c]");


  TH1F* lamOverKs_y_cor_data_7000_sys = (TH1F*)lam_dNdy_cor_data_7000_sys->Clone("lamOverKs_y_cor_data_7000_sys");
  TH1F* lamOverKs_y_cor_data_900_sys = (TH1F*)lam_dNdy_cor_data_900_sys->Clone("lamOverKs_y_cor_data_900_sys");

  lamOverKs_y_cor_data_7000_sys->Divide(ks_dNdy_cor_data_7000_sys);
  lamOverKs_y_cor_data_900_sys->Divide(ks_dNdy_cor_data_900_sys);

  lamOverKs_y_cor_data_7000_sys->SetMarkerStyle(20);
  lamOverKs_y_cor_data_7000_sys->SetMarkerColor(4);
  lamOverKs_y_cor_data_7000_sys->SetLineColor(4);

  lamOverKs_y_cor_data_900_sys->SetMarkerStyle(21);
  lamOverKs_y_cor_data_900_sys->SetMarkerColor(2);
  lamOverKs_y_cor_data_900_sys->SetLineColor(2);

  lamOverKs_y_cor_data_7000_sys->SetYTitle("N(#Lambda^{0}) / N(K_{S}^{0})");
  lamOverKs_y_cor_data_7000_sys->SetXTitle("rapidity");

  double pt7000max = lamOverKs_pT_cor_data_7000_sys->GetMaximum();
  double pt900max = lamOverKs_pT_cor_data_900_sys->GetMaximum();
  double y7000max = lamOverKs_y_cor_data_7000_sys->GetMaximum();
  double y900max = lamOverKs_y_cor_data_900_sys->GetMaximum();

  TLegend *leg_pt = new TLegend(0.2,0.25,0.5,0.35);
  leg_pt->SetTextSize(0.042);
  leg_pt->SetBorderSize(0);
  leg_pt->SetFillColor(0);
  leg_pt->SetFillStyle(0);
  leg_pt->AddEntry(lamOverKs_pT_cor_data_7000_sys, "  #sqrt{s} = 7 TeV" , "pl" );
  leg_pt->AddEntry(lamOverKs_pT_cor_data_900_sys, "  #sqrt{s} = 900 GeV" , "pl" );

  c1 = new TCanvas("c1", "lamOverKs, pT");
  lamOverKs_pT_cor_data_7000_sys->SetMinimum(0.);
  lamOverKs_pT_cor_data_7000_sys->SetMaximum(1.2); //pt7000max > pt900max ? 1.1*pt7000max : 1.1*pt900max);
  lamOverKs_pT_cor_data_7000_sys->Draw("E1");
  lamOverKs_pT_cor_data_900_sys->Draw("E1same");
  leg_pt->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90, "CMS preliminary");

  TString ratioOutpT = "final_plots/lamOverKs_pt_cor_datamc_sys";
  c1->SaveAs( ratioOutpT + ".root" );
  c1->SaveAs( ratioOutpT + ".eps" );
  c1->SaveAs( ratioOutpT + ".C" );
  c1->SaveAs( ratioOutpT + ".png" );

  delete c1;

  TLegend *leg_y = new TLegend(0.20,0.25,0.5,0.35);
  leg_y->SetTextSize(0.042);
  leg_y->SetBorderSize(0);
  leg_y->SetFillColor(0);
  leg_y->SetFillStyle(0);
  leg_y->AddEntry(lamOverKs_y_cor_data_7000_sys, "  #sqrt{s} = 7 TeV" , "pl" );
  leg_y->AddEntry(lamOverKs_y_cor_data_900_sys, "  #sqrt{s} = 900 GeV" , "pl" );

  c1 = new TCanvas("c1", "lamOverKs, y");
  lamOverKs_y_cor_data_7000_sys->SetMinimum(0.);
  lamOverKs_y_cor_data_7000_sys->SetMaximum(1.0); //y7000max > y900max ? 1.1*y7000max : 1.1*y900max);
  lamOverKs_y_cor_data_7000_sys->Draw("E1");
  lamOverKs_y_cor_data_900_sys->Draw("E1same");
  leg_y->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

  TString ratioOuty = "final_plots/lamOverKs_y_cor_datamc_sys";
  c1->SaveAs( ratioOuty + ".root" );
  c1->SaveAs( ratioOuty + ".eps" );
  c1->SaveAs( ratioOuty + ".C" );
  c1->SaveAs( ratioOuty + ".png" );

  delete c1;




ofstream outTextFile;
outTextFile.open ("fitinfo.txt");


outTextFile << "Tsallis info for data\n";
outTextFile << "$K_S^0$ & " << tsfunc900->GetParameter(1) << " $\\pm$ " << tsfunc900->GetParError(1)
  << " & " << tsfunc900->GetParameter(2) << " $\\pm$ " << tsfunc900->GetParError(2) << " & "
<< tsmean900->GetParError(1)/tsmean900->GetParameter(1) << " & "<< tsfunc7000->GetParameter(1) << " $\\pm$ " << tsfunc7000->GetParError(1)
  << " & " << tsfunc7000->GetParameter(2) << " $\\pm$ " << tsfunc7000->GetParError(2) << " & "
<< tsmean7000->GetParError(1)/tsmean7000->GetParameter(1) << " \\\\ " << endl;
outTextFile << "$\\Lambda^0$ & " << lamtsfunc900->GetParameter(1) << " $\\pm$ " << lamtsfunc900->GetParError(1)
  << " & " << lamtsfunc900->GetParameter(2) << " $\\pm$ " << lamtsfunc900->GetParError(2) << " & "
<< lamtsmean900->GetParError(1)/lamtsmean900->GetParameter(1) << " & "<< lamtsfunc7000->GetParameter(1) << " $\\pm$ " << lamtsfunc7000->GetParError(1)
  << " & " << lamtsfunc7000->GetParameter(2) << " $\\pm$ " << lamtsfunc7000->GetParError(2) << " & "
<< lamtsmean7000->GetParError(1)/lamtsmean7000->GetParameter(1) << " \\\\ " << endl;

outTextFile << endl;

outTextFile << "Tsallis info for MC\n";
outTextFile << "D6T $K_S^0$ & " << tsfunc900d6t->GetParameter(1) << " $\\pm$ " << tsfunc900d6t->GetParError(1)
  << " & " << tsfunc900d6t->GetParameter(2) << " $\\pm$ " << tsfunc900d6t->GetParError(2) << " & "
<< tsmean900d6t->GetParError(1)/tsmean900d6t->GetParameter(1) << " & "<< tsfunc7000d6t->GetParameter(1) << " $\\pm$ " << tsfunc7000d6t->GetParError(1)
  << " & " << tsfunc7000d6t->GetParameter(2) << " $\\pm$ " << tsfunc7000d6t->GetParError(2) << " & "
<< tsmean7000d6t->GetParError(1)/tsmean7000d6t->GetParameter(1) << " \\\\ " << endl;
outTextFile << "Pythia8 $K_S^0$ & " << " & " << " & "<< tsfunc7000py8->GetParameter(1) << " $\\pm$ " << tsfunc7000py8->GetParError(1)
  << " & " << tsfunc7000py8->GetParameter(2) << " $\\pm$ " << tsfunc7000py8->GetParError(2) << " & "
<< tsmean7000py8->GetParError(1)/tsmean7000py8->GetParameter(1) << " \\\\ " << endl;
outTextFile << "P0 $K_S^0$ & " << " & " << " & "<< tsfunc7000p0->GetParameter(1) << " $\\pm$ " << tsfunc7000p0->GetParError(1)
  << " & " << tsfunc7000p0->GetParameter(2) << " $\\pm$ " << tsfunc7000p0->GetParError(2) << " & "
<< tsmean7000p0->GetParError(1)/tsmean7000p0->GetParameter(1) << " \\\\ " << endl;
outTextFile << "D6T $\\Lambda^0$ & " << lamtsfunc900d6t->GetParameter(1) << " $\\pm$ " << lamtsfunc900d6t->GetParError(1)
  << " & " << lamtsfunc900d6t->GetParameter(2) << " $\\pm$ " << lamtsfunc900d6t->GetParError(2) << " & "
<< lamtsmean900d6t->GetParError(1)/lamtsmean900d6t->GetParameter(1) << " & "<< lamtsfunc7000d6t->GetParameter(1) << " $\\pm$ " << lamtsfunc7000d6t->GetParError(1)
  << " & " << lamtsfunc7000d6t->GetParameter(2) << " $\\pm$ " << lamtsfunc7000d6t->GetParError(2) << " & "
<< lamtsmean7000d6t->GetParError(1)/lamtsmean7000d6t->GetParameter(1) << " \\\\ " << endl;
outTextFile << "Pythia8 $\Lambda^0$ & " << " & " << " & "<< lamtsfunc7000py8->GetParameter(1) << " $\\pm$ " << lamtsfunc7000py8->GetParError(1)
  << " & " << lamtsfunc7000py8->GetParameter(2) << " $\\pm$ " << lamtsfunc7000py8->GetParError(2) << " & "
<< lamtsmean7000py8->GetParError(1)/lamtsmean7000py8->GetParameter(1) << " \\\\ " << endl;
outTextFile << "P0 $\Lambda^0$ & " << " & " << " & "<< lamtsfunc7000p0->GetParameter(1) << " $\\pm$ " << lamtsfunc7000p0->GetParError(1)
  << " & " << lamtsfunc7000p0->GetParameter(2) << " $\\pm$ " << lamtsfunc7000p0->GetParError(2) << " & "
<< lamtsmean7000p0->GetParError(1)/lamtsmean7000p0->GetParameter(1) << " \\\\ " << endl;

outTextFile << endl;
outTextFile << "Cross section for data\n";

inFile = "plots/kslife_cor_data_tuneD6T_900GeV_weight111_nobinomial2.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TF1* ksct900 = (TF1*)c1->GetPrimitive("ksctFit");
fin->Close();
delete c1;
delete fin;
float ctauint900 = 2.99792458e-2*exp(ksct900->GetParameter(0))*ksct900->GetParameter(1)/(trig900*lifebinwid);
float dndy900 =  ks_dNdy_cor_data_900->GetBinContent(1);
float stat900 = ks_dNdy_cor_data_900->GetBinError(1);
float sys900 = sqrt(ks_dNdy_cor_data_900_sys->GetBinError(1)*ks_dNdy_cor_data_900_sys->GetBinError(1) - stat900*stat900);
float A900 = tsfunc900->GetParameter(0);
float T900 = tsfunc900->GetParameter(1);
float n900 = tsfunc900->GetParameter(2);
float ptint900 = tsallisInt(A900,mks,n900,T900);
inFile = "plots/kslife_cor_data_pythia8_7TeV_weight102_nobinomial2.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TF1* ksct7000 = (TF1*)c1->GetPrimitive("ksctFit");
fin->Close();
delete c1;
delete fin;
float ctauint7000 = 2.99792458e-2*exp(ksct7000->GetParameter(0))*ksct7000->GetParameter(1)/(trig7000*lifebinwid);
float dndy7000 =  ks_dNdy_cor_data_7000->GetBinContent(1);
float stat7000 = ks_dNdy_cor_data_7000->GetBinError(1);
float sys7000 = sqrt(ks_dNdy_cor_data_7000_sys->GetBinError(1)*ks_dNdy_cor_data_7000_sys->GetBinError(1) - stat7000*stat7000);
float A7000 = tsfunc7000->GetParameter(0);
float T7000 = tsfunc7000->GetParameter(1);
float n7000 = tsfunc7000->GetParameter(2);
float ptint7000 = tsallisInt(A7000,mks,n7000,T7000);

outTextFile << "$K_S^0$ & " << dndy900 << " $\\pm$ " << stat900 << " $\\pm$ " << sys900 << " & " << ptint900 << " & " << 2.0*ks_yld_900_y << " & " << ctauint900 
<< " & " << dndy7000 << " $\\pm$ " << stat7000 << " $\\pm$ " << sys7000 << " & " << ptint7000 << " & " << 2.0*ks_yld_7000_y << " & " << ctauint7000 << endl;
outTextFile << "$K_S^0$ & " << dndy900*trig900/lum900 << " $\\pm$ " << stat900*trig900/lum900 << " $\\pm$ " << sys900*trig900/lum900 << " & " 
<< ptint900*trig900/lum900 << " & " << 2.0*ks_yld_900_y*trig900/lum900 << " & " << ctauint900*trig900/lum900 
<< " & " << dndy7000*trig7000/lum7000 << " $\\pm$ " << stat7000*trig7000/lum7000 << " $\\pm$ " << sys7000*trig7000/lum7000 << " & " 
<< ptint7000*trig7000/lum7000 << " & " << 2.0*ks_yld_7000_y*trig7000/lum7000 << " & " << ctauint7000*trig7000/lum7000 << endl;

inFile = "plots/lamlife_cor_data_tuneD6T_900GeV_weight111_subgeant2_nobinomial2.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TF1* lamct900 = (TF1*)c1->GetPrimitive("lamctFit");
fin->Close();
delete c1;
delete fin;
ctauint900 = 2.99792458e-2*exp(lamct900->GetParameter(0))*lamct900->GetParameter(1)/(trig900*lifebinwid);
dndy900 =  lam_dNdy_cor_data_900->GetBinContent(1);
stat900 = lam_dNdy_cor_data_900->GetBinError(1);
sys900 = sqrt(lam_dNdy_cor_data_900_sys->GetBinError(1)*lam_dNdy_cor_data_900_sys->GetBinError(1) - stat900*stat900);
A900 = lamtsfunc900->GetParameter(0);
T900 = lamtsfunc900->GetParameter(1);
n900 = lamtsfunc900->GetParameter(2);
ptint900 = tsallisInt(A900,mlam,n900,T900);

inFile = "plots/lamlife_cor_data_pythia8_7TeV_weight102_subgeant2_nobinomial2.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TF1* lamct7000 = (TF1*)c1->GetPrimitive("lamctFit");
fin->Close();
delete c1;
delete fin;
ctauint7000 = 2.99792458e-2*exp(lamct7000->GetParameter(0))*lamct7000->GetParameter(1)/(trig7000*lifebinwid);
dndy7000 =  lam_dNdy_cor_data_7000->GetBinContent(1);
stat7000 = lam_dNdy_cor_data_7000->GetBinError(1);
sys7000 = sqrt(lam_dNdy_cor_data_7000_sys->GetBinError(1)*lam_dNdy_cor_data_7000_sys->GetBinError(1) - stat7000*stat7000);
A7000 = lamtsfunc7000->GetParameter(0);
T7000 = lamtsfunc7000->GetParameter(1);
n7000 = lamtsfunc7000->GetParameter(2);
ptint7000 = tsallisInt(A7000,mlam,n7000,T7000);

outTextFile << "$\Lambda^0$ & " << dndy900 << " $\\pm$ " << stat900 << " $\\pm$ " << sys900 << " & " << ptint900 << " & " << 2.0*lam_yld_900_y << " & " << ctauint900 
<< " & " << dndy7000 << " $\\pm$ " << stat7000 << " $\\pm$ " << sys7000 << " & " << ptint7000 << " & " << 2.0*lam_yld_7000_y << " & " << ctauint7000 << endl;
outTextFile << "$\Lambda^0$ & " << dndy900*trig900/lum900 << " $\\pm$ " << stat900*trig900/lum900 << " $\\pm$ " << sys900*trig900/lum900 << " & " 
<< ptint900*trig900/lum900 << " & " << 2.0*lam_yld_900_y*trig900/lum900 << " & " << ctauint900*trig900/lum900 
<< " & " << dndy7000*trig7000/lum7000 << " $\\pm$ " << stat7000*trig7000/lum7000 << " $\\pm$ " << sys7000*trig7000/lum7000 << " & " 
<< ptint7000*trig7000/lum7000 << " & " << 2.0*lam_yld_7000_y*trig7000/lum7000 << " & " << ctauint7000*trig7000/lum7000 << endl;


outTextFile << endl;
outTextFile << "Cross section for MC\n";

float ks_dndy900d6t =  ks_dNdy_cor_d6t_900->GetBinContent(1);
float ks_stat900d6t = ks_dNdy_cor_d6t_900->GetBinError(1);
float ks_dndy7000d6t =  ks_dNdy_cor_d6t_7000->GetBinContent(1);
float ks_stat7000d6t = ks_dNdy_cor_d6t_7000->GetBinError(1);
float ks_dndy7000py8 =  ks_dNdy_cor_py8_7000->GetBinContent(1);
float ks_stat7000py8 = ks_dNdy_cor_py8_7000->GetBinError(1);
float ks_dndy7000p0 =  ks_dNdy_cor_p0_7000->GetBinContent(1);
float ks_stat7000p0 = ks_dNdy_cor_p0_7000->GetBinError(1);

outTextFile << "$K_S^0$ & " << ks_dndy900d6t << " &  & " << ks_dndy7000d6t << " &  & " << ks_dndy7000py8 << " &  & " << ks_dndy7000p0 << " &   \\\\ " << endl;

float lam_dndy900d6t =  lam_dNdy_cor_d6t_900->GetBinContent(1);
float lam_stat900d6t = lam_dNdy_cor_d6t_900->GetBinError(1);
float lam_dndy7000d6t =  lam_dNdy_cor_d6t_7000->GetBinContent(1);
float lam_stat7000d6t = lam_dNdy_cor_d6t_7000->GetBinError(1);
float lam_dndy7000py8 =  lam_dNdy_cor_py8_7000->GetBinContent(1);
float lam_stat7000py8 = lam_dNdy_cor_py8_7000->GetBinError(1);
float lam_dndy7000p0 =  lam_dNdy_cor_p0_7000->GetBinContent(1);
float lam_stat7000p0 = lam_dNdy_cor_p0_7000->GetBinError(1);

outTextFile << "$\Lambda^0$ & " << lam_dndy900d6t << " &  & " << lam_dndy7000d6t << " &  & " << lam_dndy7000py8 << " &  & " << lam_dndy7000p0 << " &  \\\\ " << endl;

outTextFile.close();


}
//Pythia 8 7 TeV Ks, Lam: 0.637983 +- 0.000111686 and 0.790338 +- 0.000206491
//D6T 7 TeV Ks, Lam: 0.667729 +- 6.37524e-05 and 0.90363 +- 0.000121645
//P0 7 TeV Ks, Lam: 0.662035 +- 0.000137968 and 0.776658 +- 0.000246882
//D6T 900 GeV Ks, Lam: 0.547424 +- 0.000120983 and 0.721044 +- 0.000220097

//Pythia 8 7 TeV Ks, Lam Yields: 11137377, 1.2800001 and 4226842, 0.48578388
//Pythia 8 7 TeV Ks, Lam averge pt: 0.70895094 +- 0.00019298827 and 0.92374253 +- 0.00038351648
//D6T 7 TeV Ks, Lam Yields: 34736272, 1.134765 and 13662389, 0.44632312
//D6T 7 TeV Ks, Lam average pt: 0.75002491 +- 0.00011501105 and 1.0616008 +- 0.00023929476
//P0 7 TeV Ks, Lam Yields: 9749505, 1.1032481 and 3309597, 0.374512
//P0 7 TeV Ks, Lam average pt: 0.72417134 +- 0.00022099366 and 0.90678853 +- 0.00043408293
//D6T 900 GeV Ks, Lam Yields: 6280275, 0.80099653 and 2382457, 0.30386246
//D6T 900 GeV Ks, Lam average pt: 0.57696122 +- 0.00017924383 and 0.75176877 +- 0.0003590417
