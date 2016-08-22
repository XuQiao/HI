{
#include <tdrstyle.C>
#include "tsallis.C"
#include "calc_integral.C"

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

double trig7000 = 15540204;
double trig900 = 12325647;
double lum7000 = 291850.0; // inverse millibarns
double lum900 = 321970.0; // inverse millibarns
double ptbinwid = 0.1;
double ybinwid = 0.2; // 0.2 since it is a 0.1 bin on each side of 0.

TF1 *tsfunc900 = new TF1("tsfunc900",tsallis,0.0,10.0,3);
tsfunc900->SetParNames("Amplitude","T","n");
tsfunc900->SetParameters(10000000.0.,0.2,7.2); 
TF1 *tsfunc7000 = new TF1("tsfunc7000",tsallis,0.0,10.0,3);
tsfunc7000->SetParNames("Amplitude","T","n");
tsfunc7000->SetParameters(10000000.0.,0.2,7.2); 

TLatex *cmsprelim = new TLatex();
cmsprelim->SetNDC();
cmsprelim->SetTextSize(0.045);
cmsprelim->SetTextAlign(22);

TString inFile;
TString outFile_pt_evt = "final_plots/kspt_cor_data_perevt";
TString outFile_pt_lum = "final_plots/kspt_cor_data_perlum";

Double_t ptent[44];

TString outFile_y_evt = "final_plots/ksy_cor_data_perevt";
TString outFile_y_lum = "final_plots/ksy_cor_data_perlum";
 TString outFile_y_pos_evt = "final_plots/ksy_pos_cor_data_perevt";
 TString outFile_y_pos_lum = "final_plots/ksy_pos_cor_data_perlum";
 TString outFile_y_neg_evt = "final_plots/ksy_neg_cor_data_perevt";
 TString outFile_y_neg_lum = "final_plots/ksy_neg_cor_data_perlum";

Double_t *yent = new Double_t[24];


// Pt for Kshort per event

inFile = "plots/kspt_cor_data_pythia8_7TeV_weight102.root";
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
TH1F *ks_dNdpT_cor_data_7000_sys = (TH1F*)ks_dNdpT_cor_data_7000->Clone("ks_dNdpT_cor_data_700_sys");
//cout << "ks pt 7000 integral / event: " << ks_yld_7000_pt << endl;
for (int i = 1; i <= ks_dNdpT_cor_data_7000_sys->GetXaxis()->GetNbins(); ++i) {
  double value = ks_dNdpT_cor_data_7000_sys->GetBinContent(i);
  double error = ks_dNdpT_cor_data_7000_sys->GetBinError(i);
  double syserr = 0.02*value;
  ptent[i] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdpT_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i] << endl;
}
ks_dNdpT_cor_data_7000_sys->SetError(ptent);
ks_dNdpT_cor_data_7000->Fit("tsfunc7000","EM0");
tsfunc7000->SetNpx(2000);
tsfunc7000->SetLineColor(4);

inFile = "plots/kspt_cor_data_tuneD6T_900GeV_weight111.root";
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
TH1F *ks_dNdpT_cor_data_900_sys = (TH1F*)ks_dNdpT_cor_data_900->Clone("ks_dNdpT_cor_data_700_sys");
//cout << "ks pt 900 integral / event: " << ks_yld_900_pt << endl;
for (int i = 1; i <= ks_dNdpT_cor_data_900_sys->GetXaxis()->GetNbins(); ++i) {
  double value = ks_dNdpT_cor_data_900_sys->GetBinContent(i);
  double error = ks_dNdpT_cor_data_900_sys->GetBinError(i);
  double syserr = 0.02*value;
  ptent[i] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdpT_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i] << endl;
}
ks_dNdpT_cor_data_900_sys->SetError(ptent);
ks_dNdpT_cor_data_900->Fit("tsfunc900","EM0");
tsfunc900->SetNpx(2000);
tsfunc900->SetLineColor(2);

TLegend *leg = new TLegend(0.50,0.70,0.85,0.80);
leg->SetTextSize(0.042);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(ks_dNdpT_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
leg->AddEntry(ks_dNdpT_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

ks_dNdpT_cor_data_7000_sys->Draw();
ks_dNdpT_cor_data_7000->Draw("E1same");
tsfunc7000->Draw("same"); 
ks_dNdpT_cor_data_900_sys->Draw("same");
ks_dNdpT_cor_data_900->Draw("E1same");
tsfunc900->Draw("same"); 
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

c1->SaveAs(outFile_pt_evt + ".root");
c1->SaveAs(outFile_pt_evt + ".eps");
c1->SaveAs(outFile_pt_evt + ".C");
c1->SaveAs(outFile_pt_evt + ".png");

delete c1;


TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

ks_dNdpT_cor_data_7000_sys->SetMaximum(3.0);
ks_dNdpT_cor_data_7000_sys->SetMinimum(0.000004);
ks_dNdpT_cor_data_7000_sys->Draw();
ks_dNdpT_cor_data_7000->Draw("E1same");
tsfunc7000->Draw("same"); 
ks_dNdpT_cor_data_900_sys->Draw("same");
ks_dNdpT_cor_data_900->Draw("E1same");
tsfunc900->Draw("same"); 
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

gPad->SetLogy(kTRUE);
c1->SaveAs(outFile_pt_evt + "_log.root");
c1->SaveAs(outFile_pt_evt + "_log.eps");
c1->SaveAs(outFile_pt_evt + "_log.C");
c1->SaveAs(outFile_pt_evt + "_log.png");
gPad->SetLogy(kFALSE);

delete c1;



// Pt for Kshort per luminosity

ks_dNdpT_cor_data_7000->Scale(trig7000/lum7000);
double ks_xsec_7000_pt = calc_integral(ks_dNdpT_cor_data_7000, 1.0);
//cout << "ks pt 7000 xsec: " << ks_xsec_7000_pt << endl;
ks_dNdpT_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow K_{S}^{0}X) / dp_{T}  [ mb / (GeV/c)]");
ks_dNdpT_cor_data_7000_sys->Scale(trig7000/lum7000);
ks_dNdpT_cor_data_7000->Fit("tsfunc7000","EM0");
tsfunc7000->SetNpx(2000);
tsfunc7000->SetLineColor(4);

ks_dNdpT_cor_data_900->Scale(trig900/lum900);
double ks_xsec_900_pt = calc_integral(ks_dNdpT_cor_data_900, 1.0);
//cout << "ks pt 900 xsec: " << ks_xsec_900_pt << endl;
ks_dNdpT_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow K_{S}^{0}X) / dp_{T}  [ mb / (GeV/c)]");
ks_dNdpT_cor_data_900_sys->Scale(trig900/lum900);
//cout << "pt 900\n";
ks_dNdpT_cor_data_900->Fit("tsfunc900","EM0");
tsfunc900->SetNpx(2000);
tsfunc900->SetLineColor(2);

TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

ks_dNdpT_cor_data_7000_sys->Draw();
ks_dNdpT_cor_data_7000->Draw("E1same");
tsfunc7000->Draw("same"); 
ks_dNdpT_cor_data_900_sys->Draw("same");
ks_dNdpT_cor_data_900->Draw("E1same");
tsfunc900->Draw("same"); 
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

c1->SaveAs(outFile_pt_lum + ".root");
c1->SaveAs(outFile_pt_lum + ".eps");
c1->SaveAs(outFile_pt_lum + ".C");
c1->SaveAs(outFile_pt_lum + ".png");

delete c1;


TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

ks_dNdpT_cor_data_7000_sys->SetMaximum(200.0);
ks_dNdpT_cor_data_7000_sys->SetMinimum(0.0002);
ks_dNdpT_cor_data_7000_sys->Draw();
ks_dNdpT_cor_data_7000->Draw("E1same");
tsfunc7000->Draw("same"); 
ks_dNdpT_cor_data_900_sys->Draw("same");
ks_dNdpT_cor_data_900->Draw("E1same");
tsfunc900->Draw("same"); 
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

gPad->SetLogy(kTRUE);
c1->SaveAs(outFile_pt_lum + "_log.root");
c1->SaveAs(outFile_pt_lum + "_log.eps");
c1->SaveAs(outFile_pt_lum + "_log.C");
c1->SaveAs(outFile_pt_lum + "_log.png");
gPad->SetLogy(kFALSE);

delete c1;







// y for Kshort

inFile = "plots/ksy_cor_data_pythia8_7TeV_weight102.root";
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
TH1F *ks_dNdy_cor_data_7000_sys = (TH1F*)ks_dNdy_cor_data_7000->Clone("ks_dNdy_cor_data_700_sys");
//cout << "ks y 7000 integral / event: " << ks_yld_7000_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_cor_data_7000_sys->GetBinContent(i);
  double error = ks_dNdy_cor_data_7000_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdy_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
}
ks_dNdy_cor_data_7000_sys->SetError(yent);

inFile = "plots/ksy_cor_data_tuneD6T_900GeV_weight111.root";
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
TH1F *ks_dNdy_cor_data_900_sys = (TH1F*)ks_dNdy_cor_data_900->Clone("ks_dNdy_cor_data_700_sys");
//cout << "ks y 900 integral / event: " << ks_yld_900_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_cor_data_900_sys->GetBinContent(i);
  double error = ks_dNdy_cor_data_900_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdy_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
}
ks_dNdy_cor_data_900_sys->SetError(yent);

TLegend *leg = new TLegend(0.50,0.25,0.85,0.35);
leg->SetTextSize(0.042);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(ks_dNdy_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
leg->AddEntry(ks_dNdy_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

ks_dNdy_cor_data_7000_sys->Draw();
ks_dNdy_cor_data_7000->Draw("E1same");
ks_dNdy_cor_data_900_sys->Draw("same");
ks_dNdy_cor_data_900->Draw("E1same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

c1->SaveAs(outFile_y_evt + ".root");
c1->SaveAs(outFile_y_evt + ".eps");
c1->SaveAs(outFile_y_evt + ".C");
c1->SaveAs(outFile_y_evt + ".png");

delete c1;



// y for Kshort per luminosity

ks_dNdy_cor_data_7000->Scale(trig7000/lum7000);
double ks_xsec_7000_y = calc_integral(ks_dNdy_cor_data_7000, 1.0);
//cout << "ks y 7000 xsec: " << ks_xsec_7000_y << endl;
ks_dNdy_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow K_{S}^{0}X) / dy  [ mb ]");
ks_dNdy_cor_data_7000_sys->Scale(trig7000/lum7000);

ks_dNdy_cor_data_900->Scale(trig900/lum900);
double ks_xsec_900_y = calc_integral(ks_dNdy_cor_data_900, 1.0);
//cout << "ks y 900 xsec: " << ks_xsec_900_y << endl;
ks_dNdy_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow K_{S}^{0}X) / dy  [ mb ]");
ks_dNdy_cor_data_900_sys->Scale(trig900/lum900);
//cout << "y 900\n";

TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

ks_dNdy_cor_data_7000_sys->SetMaximum(22.0);
ks_dNdy_cor_data_7000_sys->SetMinimum(0.);
ks_dNdy_cor_data_7000_sys->Draw();
ks_dNdy_cor_data_7000->Draw("E1same");
ks_dNdy_cor_data_900_sys->Draw("same");
ks_dNdy_cor_data_900->Draw("E1same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

c1->SaveAs(outFile_y_lum + ".root");
c1->SaveAs(outFile_y_lum + ".eps");
c1->SaveAs(outFile_y_lum + ".C");
c1->SaveAs(outFile_y_lum + ".png");

delete c1;

//-------------
// y+ for Kshort

inFile = "plots/ksy_pos_cor_data_pythia8_7TeV_weight102.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_pos_cor_data_7000 = (TH1F*)c1->GetPrimitive("ks_dNdy_pos");
fin->Close();
delete c1;
delete fin;
double yld_7000_y_pos = ks_dNdy_pos_cor_data_7000->Integral();
ks_dNdy_pos_cor_data_7000->Scale(1.0/(trig7000*ybinwid));
double ks_yld_7000_y_pos = calc_integral(ks_dNdy_pos_cor_data_7000, 1.0);
ks_dNdy_pos_cor_data_7000->GetXaxis()->SetTitle("K_{S}^{0}  rapidity > 0");
ks_dNdy_pos_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_pos_cor_data_7000->SetMaximum(0.4);
ks_dNdy_pos_cor_data_7000->SetMinimum(0.);
ks_dNdy_pos_cor_data_7000->SetLineColor(4);
ks_dNdy_pos_cor_data_7000->SetMarkerColor(4);
ks_dNdy_pos_cor_data_7000->SetMarkerStyle(20);
ks_dNdy_pos_cor_data_7000->SetMarkerSize(0.8);
TH1F *ks_dNdy_pos_cor_data_7000_sys = (TH1F*)ks_dNdy_pos_cor_data_7000->Clone("ks_dNdy_pos_cor_data_700_sys");
//cout << "ks y 7000 integral / event: " << ks_yld_7000_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_pos_cor_data_7000_sys->GetBinContent(i);
  double error = ks_dNdy_pos_cor_data_7000_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdy_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
}
ks_dNdy_pos_cor_data_7000_sys->SetError(yent);

inFile = "plots/ksy_pos_cor_data_tuneD6T_900GeV_weight111.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_pos_cor_data_900 = (TH1F*)c1->GetPrimitive("ks_dNdy_pos");
fin->Close();
delete c1;
delete fin;
double yld_900_y_pos = ks_dNdy_pos_cor_data_900->Integral();
ks_dNdy_pos_cor_data_900->Scale(1.0/(trig900*ybinwid));
double ks_yld_900_y_pos = calc_integral(ks_dNdy_pos_cor_data_900, 1.0);
ks_dNdy_pos_cor_data_900->GetXaxis()->SetTitle("K_{S}^{0}  rapidity > 0");
ks_dNdy_pos_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_pos_cor_data_900->SetMaximum(0.4);
ks_dNdy_pos_cor_data_900->SetMinimum(0.);
ks_dNdy_pos_cor_data_900->SetLineColor(2);
ks_dNdy_pos_cor_data_900->SetMarkerColor(2);
ks_dNdy_pos_cor_data_900->SetMarkerStyle(21);
ks_dNdy_pos_cor_data_900->SetMarkerSize(0.8);
TH1F *ks_dNdy_pos_cor_data_900_sys = (TH1F*)ks_dNdy_pos_cor_data_900->Clone("ks_dNdy_pos_cor_data_700_sys");
//cout << "ks y 900 integral / event: " << ks_yld_900_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_pos_cor_data_900_sys->GetBinContent(i);
  double error = ks_dNdy_pos_cor_data_900_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdy_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
}
ks_dNdy_pos_cor_data_900_sys->SetError(yent);

TLegend *leg = new TLegend(0.50,0.25,0.85,0.35);
leg->SetTextSize(0.042);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(ks_dNdy_pos_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
leg->AddEntry(ks_dNdy_pos_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

ks_dNdy_pos_cor_data_7000_sys->Draw();
ks_dNdy_pos_cor_data_7000->Draw("E1same");
ks_dNdy_pos_cor_data_900_sys->Draw("same");
ks_dNdy_pos_cor_data_900->Draw("E1same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

c1->SaveAs(outFile_y_pos_evt + ".root");
c1->SaveAs(outFile_y_pos_evt + ".eps");
c1->SaveAs(outFile_y_pos_evt + ".C");
c1->SaveAs(outFile_y_pos_evt + ".png");

delete c1;



// y+ for Kshort per luminosity

ks_dNdy_pos_cor_data_7000->Scale(trig7000/lum7000);
double ks_xsec_7000_y_pos = calc_integral(ks_dNdy_pos_cor_data_7000, 1.0);
//cout << "ks y 7000 xsec: " << ks_xsec_7000_y << endl;
ks_dNdy_pos_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow K_{S}^{0}X) / dy  [ mb ]");
ks_dNdy_pos_cor_data_7000_sys->Scale(trig7000/lum7000);

ks_dNdy_pos_cor_data_900->Scale(trig900/lum900);
double ks_xsec_900_y_pos = calc_integral(ks_dNdy_pos_cor_data_900, 1.0);
//cout << "ks y 900 xsec: " << ks_xsec_900_y << endl;
ks_dNdy_pos_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow K_{S}^{0}X) / dy  [ mb ]");
ks_dNdy_pos_cor_data_900_sys->Scale(trig900/lum900);
//cout << "y 900\n";

TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

ks_dNdy_pos_cor_data_7000_sys->SetMaximum(22.0);
ks_dNdy_pos_cor_data_7000_sys->SetMinimum(0.);
ks_dNdy_pos_cor_data_7000_sys->Draw();
ks_dNdy_pos_cor_data_7000->Draw("E1same");
ks_dNdy_pos_cor_data_900_sys->Draw("same");
ks_dNdy_pos_cor_data_900->Draw("E1same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

c1->SaveAs(outFile_y_pos_lum + ".root");
c1->SaveAs(outFile_y_pos_lum + ".eps");
c1->SaveAs(outFile_y_pos_lum + ".C");
c1->SaveAs(outFile_y_pos_lum + ".png");

delete c1;

//------------------

// y- for Kshort

inFile = "plots/ksy_neg_cor_data_pythia8_7TeV_weight102.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_neg_cor_data_7000 = (TH1F*)c1->GetPrimitive("ks_dNdy_neg");
fin->Close();
delete c1;
delete fin;
double yld_7000_y_neg = ks_dNdy_neg_cor_data_7000->Integral();
ks_dNdy_neg_cor_data_7000->Scale(1.0/(trig7000*ybinwid));
double ks_yld_7000_y_neg = calc_integral(ks_dNdy_neg_cor_data_7000, 1.0);
ks_dNdy_neg_cor_data_7000->GetXaxis()->SetTitle("K_{S}^{0}  rapidity < 0");
ks_dNdy_neg_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_neg_cor_data_7000->SetMaximum(0.4);
ks_dNdy_neg_cor_data_7000->SetMinimum(0.);
ks_dNdy_neg_cor_data_7000->SetLineColor(4);
ks_dNdy_neg_cor_data_7000->SetMarkerColor(4);
ks_dNdy_neg_cor_data_7000->SetMarkerStyle(20);
ks_dNdy_neg_cor_data_7000->SetMarkerSize(0.8);
TH1F *ks_dNdy_neg_cor_data_7000_sys = (TH1F*)ks_dNdy_neg_cor_data_7000->Clone("ks_dNdy_neg_cor_data_700_sys");
//cout << "ks y 7000 integral / event: " << ks_yld_7000_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_neg_cor_data_7000_sys->GetBinContent(i);
  double error = ks_dNdy_neg_cor_data_7000_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdy_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
}
ks_dNdy_neg_cor_data_7000_sys->SetError(yent);

inFile = "plots/ksy_neg_cor_data_tuneD6T_900GeV_weight111.root";
TFile* fin = TFile::Open(inFile, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* ks_dNdy_neg_cor_data_900 = (TH1F*)c1->GetPrimitive("ks_dNdy_neg");
fin->Close();
delete c1;
delete fin;
double yld_900_y_neg = ks_dNdy_neg_cor_data_900->Integral();
ks_dNdy_neg_cor_data_900->Scale(1.0/(trig900*ybinwid));
double ks_yld_900_y_neg = calc_integral(ks_dNdy_neg_cor_data_900, 1.0);
ks_dNdy_neg_cor_data_900->GetXaxis()->SetTitle("K_{S}^{0}  rapidity < 0");
ks_dNdy_neg_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow K_{S}^{0}X) / dy");
ks_dNdy_neg_cor_data_900->SetMaximum(0.4);
ks_dNdy_neg_cor_data_900->SetMinimum(0.);
ks_dNdy_neg_cor_data_900->SetLineColor(2);
ks_dNdy_neg_cor_data_900->SetMarkerColor(2);
ks_dNdy_neg_cor_data_900->SetMarkerStyle(21);
ks_dNdy_neg_cor_data_900->SetMarkerSize(0.8);
TH1F *ks_dNdy_neg_cor_data_900_sys = (TH1F*)ks_dNdy_neg_cor_data_900->Clone("ks_dNdy_neg_cor_data_700_sys");
//cout << "ks y 900 integral / event: " << ks_yld_900_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_neg_cor_data_900_sys->GetBinContent(i);
  double error = ks_dNdy_neg_cor_data_900_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i] = sqrt(error*error + syserr*syserr);
  //cout << ks_dNdy_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
}
ks_dNdy_neg_cor_data_900_sys->SetError(yent);

TLegend *leg = new TLegend(0.50,0.25,0.85,0.35);
leg->SetTextSize(0.042);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(ks_dNdy_neg_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
leg->AddEntry(ks_dNdy_neg_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

ks_dNdy_neg_cor_data_7000_sys->Draw();
ks_dNdy_neg_cor_data_7000->Draw("E1same");
ks_dNdy_neg_cor_data_900_sys->Draw("same");
ks_dNdy_neg_cor_data_900->Draw("E1same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

c1->SaveAs(outFile_y_neg_evt + ".root");
c1->SaveAs(outFile_y_neg_evt + ".eps");
c1->SaveAs(outFile_y_neg_evt + ".C");
c1->SaveAs(outFile_y_neg_evt + ".png");

delete c1;



// y- for Kshort per luminosity

ks_dNdy_neg_cor_data_7000->Scale(trig7000/lum7000);
double ks_xsec_7000_y_neg = calc_integral(ks_dNdy_neg_cor_data_7000, 1.0);
//cout << "ks y 7000 xsec: " << ks_xsec_7000_y << endl;
ks_dNdy_neg_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow K_{S}^{0}X) / dy  [ mb ]");
ks_dNdy_neg_cor_data_7000_sys->Scale(trig7000/lum7000);

ks_dNdy_neg_cor_data_900->Scale(trig900/lum900);
double ks_xsec_900_y_neg = calc_integral(ks_dNdy_neg_cor_data_900, 1.0);
//cout << "ks y 900 xsec: " << ks_xsec_900_y << endl;
ks_dNdy_neg_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow K_{S}^{0}X) / dy  [ mb ]");
ks_dNdy_neg_cor_data_900_sys->Scale(trig900/lum900);
//cout << "y 900\n";

TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

ks_dNdy_neg_cor_data_7000_sys->SetMaximum(22.0);
ks_dNdy_neg_cor_data_7000_sys->SetMinimum(0.);
ks_dNdy_neg_cor_data_7000_sys->Draw();
ks_dNdy_neg_cor_data_7000->Draw("E1same");
ks_dNdy_neg_cor_data_900_sys->Draw("same");
ks_dNdy_neg_cor_data_900->Draw("E1same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

c1->SaveAs(outFile_y_neg_lum + ".root");
c1->SaveAs(outFile_y_neg_lum + ".eps");
c1->SaveAs(outFile_y_neg_lum + ".C");
c1->SaveAs(outFile_y_neg_lum + ".png");

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

  TF1 *lamtsfunc900 = new TF1("lamtsfunc900",tsallis,0.0, 10.0, 3);
  lamtsfunc900->SetParNames("Amplitude","T","n");
  lamtsfunc900->SetParameters(10000000.0., 0.2, 7.2); 
  TF1 *lamtsfunc7000 = new TF1("lamtsfunc7000",tsallis,0.0, 10.0, 3);
  lamtsfunc7000->SetParNames("Amplitude","T","n");
  lamtsfunc7000->SetParameters(10000000.0., 0.2, 7.2); 

  TF1 *lamRegtsfunc900 = new TF1("lamRegtsfunc900",tsallis,0.0, 10.0, 3);
  lamRegtsfunc900->SetParNames("Amplitude","T","n");
  lamRegtsfunc900->SetParameters(10000000.0., 0.2, 7.2); 
  TF1 *lamRegtsfunc7000 = new TF1("lamRegtsfunc7000",tsallis,0.0, 10.0, 3);
  lamRegtsfunc7000->SetParNames("Amplitude","T","n");
  lamRegtsfunc7000->SetParameters(10000000.0., 0.2, 7.2); 

  TF1 *lamBartsfunc900 = new TF1("lamBartsfunc900",tsallis,0.0, 10.0, 3);
  lamBartsfunc900->SetParNames("Amplitude","T","n");
  lamBartsfunc900->SetParameters(10000000.0., 0.2, 7.2); 
  TF1 *lamBartsfunc7000 = new TF1("lamBartsfunc7000",tsallis,0.0, 10.0, 3);
  lamBartsfunc7000->SetParNames("Amplitude","T","n");
  lamBartsfunc7000->SetParameters(10000000.0., 0.2, 7.2); 

  //TString inFile;
  TString outFile_pt_evt = "final_plots/lampt_cor_data_perevt";
  TString outFile_pt_lum = "final_plots/lampt_cor_data_perlum";
  TString outFileReg_pt_evt = "final_plots/lamRegpt_cor_data_perevt";
  TString outFileBar_pt_evt = "final_plots/lamBarpt_cor_data_perevt";
  TString outFileReg_pt_lum = "final_plots/lamRegpt_cor_data_perlum";
  TString outFileBar_pt_lum = "final_plots/lamBarpt_cor_data_perlum";

  Double_t ptent[41];

  TString outFile_y_evt = "final_plots/lamy_cor_data_perevt";
  TString outFile_y_lum = "final_plots/lamy_cor_data_perlum";
  TString outFileReg_y_evt = "final_plots/lamRegy_cor_data_perevt";
  TString outFileBar_y_evt = "final_plots/lamBary_cor_data_perevt";
  TString outFileReg_y_lum = "final_plots/lamRegy_cor_data_perlum";
  TString outFileBar_y_lum = "final_plots/lamBary_cor_data_perlum";

  delete[] yent;
  Double_t *yent = new Double_t[24];


  // Pt for Lambda(inclusive) per event

  inFile = "plots/lampt_cor_data_pythia8_7TeV_weight102_subgeant2.root";
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
  TH1F *lam_dNdpT_cor_data_7000_sys = (TH1F*)lam_dNdpT_cor_data_7000->Clone("lam_dNdpT_cor_data_700_sys");
  //cout << "lam pt 7000 integral / event: " << lam_yld_7000_pt << endl;
  for (int i = 1; i <= lam_dNdpT_cor_data_7000_sys->GetXaxis()->GetNbins(); ++i) {
    double value = lam_dNdpT_cor_data_7000_sys->GetBinContent(i);
    double error = lam_dNdpT_cor_data_7000_sys->GetBinError(i);
    double syserr = 0.03*value;
    ptent[i] = sqrt(error*error + syserr*syserr);
    //cout << lam_dNdpT_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i] << endl;
  }
  lam_dNdpT_cor_data_7000_sys->SetError(ptent);
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

  inFile = "plots/lampt_cor_data_tuneD6T_900GeV_weight111_subgeant2.root";
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
  TH1F *lam_dNdpT_cor_data_900_sys = (TH1F*)lam_dNdpT_cor_data_900->Clone("lam_dNdpT_cor_data_700_sys");
  //cout << "lam pt 900 integral / event: " << lam_yld_900_pt << endl;
  for (int i = 1; i <= lam_dNdpT_cor_data_900_sys->GetXaxis()->GetNbins(); ++i) {
    double value = lam_dNdpT_cor_data_900_sys->GetBinContent(i);
    double error = lam_dNdpT_cor_data_900_sys->GetBinError(i);
    double syserr = 0.03*value;
    ptent[i] = sqrt(error*error + syserr*syserr);
    //cout << lam_dNdpT_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i] << endl;
  }
  lam_dNdpT_cor_data_900_sys->SetError(ptent);
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

  TLegend *leg = new TLegend(0.50,0.70,0.85,0.80);
  leg->SetTextSize(0.042);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(lam_dNdpT_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
  leg->AddEntry(lam_dNdpT_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lam_dNdpT_cor_data_7000->Draw("E1");
  //oldHisto->Draw("E1");

  //c1->SaveAs("test.png");
  lam_dNdpT_cor_data_7000_sys->Draw("same");
  lam_dNdpT_cor_data_7000->Draw("E1same");
  lamtsfunc7000->Draw("same"); 
  lam_dNdpT_cor_data_900_sys->Draw("same");
  lam_dNdpT_cor_data_900->Draw("E1same");
  lamtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFile_pt_evt + ".root");
  c1->SaveAs(outFile_pt_evt + ".eps");
  c1->SaveAs(outFile_pt_evt + ".C");
  c1->SaveAs(outFile_pt_evt + ".png");

  delete c1;


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lam_dNdpT_cor_data_7000->SetMaximum(1.0);
  lam_dNdpT_cor_data_7000->SetMinimum(0.000001);
  lam_dNdpT_cor_data_7000->Draw("E1");
  lam_dNdpT_cor_data_7000_sys->Draw("same");
  lam_dNdpT_cor_data_7000->Draw("E1same");
  lamtsfunc7000->Draw("same"); 
  lam_dNdpT_cor_data_900_sys->Draw("same");
  lam_dNdpT_cor_data_900->Draw("E1same");
  lamtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  gPad->SetLogy(kTRUE);
  c1->SaveAs(outFile_pt_evt + "_log.root");
  c1->SaveAs(outFile_pt_evt + "_log.eps");
  c1->SaveAs(outFile_pt_evt + "_log.C");
  c1->SaveAs(outFile_pt_evt + "_log.png");
  gPad->SetLogy(kFALSE);

  delete c1;

  //return;

  // Switch histograms again
  switch7000 = lam_dNdpT_cor_data_7000;
  lam_dNdpT_cor_data_7000 = tmp7;
  tmp7 = switch7000;

  switch900 = lam_dNdpT_cor_data_900;
  lam_dNdpT_cor_data_900 = tmp9;
  tmp9 = switch900;

  // Pt for Lam per luminosity

  lam_dNdpT_cor_data_7000->Scale(trig7000/lum7000);
  double lam_xsec_7000_pt = calc_integral(lam_dNdpT_cor_data_7000, 1.0);
  //cout << "lam pt 7000 xsec: " << lam_xsec_7000_pt << endl;
  lam_dNdpT_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dp_{T}  [ mb / (GeV/c)]");
  lam_dNdpT_cor_data_7000_sys->Scale(trig7000/lum7000);
  lam_dNdpT_cor_data_7000->Fit("lamtsfunc7000","EM0");
  switch7000 = lam_dNdpT_cor_data_7000;
  lam_dNdpT_cor_data_7000 = tmp7;
  tmp7 = switch7000;
  lam_dNdpT_cor_data_7000->Scale(trig7000/lum7000);
  lamtsfunc7000->SetNpx(2000);
  lamtsfunc7000->SetLineColor(4);

  lam_dNdpT_cor_data_900->Scale(trig900/lum900);
  double lam_xsec_900_pt = calc_integral(lam_dNdpT_cor_data_900, 1.0);
  //cout << "lam pt 900 xsec: " << lam_xsec_900_pt << endl;
  lam_dNdpT_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dp_{T}  [ mb / (GeV/c)]");
  lam_dNdpT_cor_data_900_sys->Scale(trig900/lum900);
  //cout << "pt 900\n";
  lam_dNdpT_cor_data_900->Fit("lamtsfunc900","EM0");
  switch900 = lam_dNdpT_cor_data_900;
  lam_dNdpT_cor_data_900 = tmp9;
  tmp9 = switch900;
  lam_dNdpT_cor_data_900->Scale(trig900/lum900);
  lamtsfunc900->SetNpx(2000);
  lamtsfunc900->SetLineColor(2);

  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lam_dNdpT_cor_data_7000->Draw("E1");
  lam_dNdpT_cor_data_7000_sys->Draw("same");
  lam_dNdpT_cor_data_7000->Draw("E1same");
  lamtsfunc7000->Draw("same"); 
  lam_dNdpT_cor_data_900_sys->Draw("same");
  lam_dNdpT_cor_data_900->Draw("E1same");
  lamtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFile_pt_lum + ".root");
  c1->SaveAs(outFile_pt_lum + ".eps");
  c1->SaveAs(outFile_pt_lum + ".C");
  c1->SaveAs(outFile_pt_lum + ".png");

  delete c1;


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lam_dNdpT_cor_data_7000->SetMaximum(70.);
  lam_dNdpT_cor_data_7000->SetMinimum(0.00007);
  lam_dNdpT_cor_data_7000_sys->SetMaximum(70.0);
  lam_dNdpT_cor_data_7000_sys->SetMinimum(0.00007);
  lam_dNdpT_cor_data_7000->Draw("E1");
  lam_dNdpT_cor_data_7000_sys->Draw("same");
  lam_dNdpT_cor_data_7000->Draw("E1same");
  lamtsfunc7000->Draw("same"); 
  lam_dNdpT_cor_data_900_sys->Draw("same");
  lam_dNdpT_cor_data_900->Draw("E1same");
  lamtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  gPad->SetLogy(kTRUE);
  c1->SaveAs(outFile_pt_lum + "_log.root");
  c1->SaveAs(outFile_pt_lum + "_log.eps");
  c1->SaveAs(outFile_pt_lum + "_log.C");
  c1->SaveAs(outFile_pt_lum + "_log.png");
  gPad->SetLogy(kFALSE);

  delete c1;

  //---------------

  // Pt for Lambda(particle only) per event

  inFile = "plots/lamRegpt_cor_data_pythia8_7TeV_weight102_subgeant2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lamReg_dNdpT_cor_data_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_scale");
  fin->Close();
  delete c1;
  delete fin;
  lamReg_dNdpT_cor_data_7000->Scale(1.0/(trig7000*ptbinwid));
  double lamReg_yld_7000_pt = calc_integral(lamReg_dNdpT_cor_data_7000, 1.0);
  lamReg_dNdpT_cor_data_7000->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  lamReg_dNdpT_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  lamReg_dNdpT_cor_data_7000->SetLineColor(4);
  lamReg_dNdpT_cor_data_7000->SetMarkerColor(4);
  lamReg_dNdpT_cor_data_7000->SetMarkerStyle(20);
  lamReg_dNdpT_cor_data_7000->SetMarkerSize(0.8);
  TH1F *lamReg_dNdpT_cor_data_7000_sys = (TH1F*)lamReg_dNdpT_cor_data_7000->Clone("lamReg_dNdpT_cor_data_700_sys");
  //cout << "lamReg pt 7000 integral / event: " << lamReg_yld_7000_pt << endl;
  for (int i = 1; i <= lamReg_dNdpT_cor_data_7000_sys->GetXaxis()->GetNbins(); ++i) {
    double value = lamReg_dNdpT_cor_data_7000_sys->GetBinContent(i);
    double error = lamReg_dNdpT_cor_data_7000_sys->GetBinError(i);
    double syserr = 0.03*value;
    ptent[i] = sqrt(error*error + syserr*syserr);
    //cout << lamReg_dNdpT_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i] << endl;
  }
  lamReg_dNdpT_cor_data_7000_sys->SetError(ptent);
  lamReg_dNdpT_cor_data_7000->Fit("lamRegtsfunc7000","EM0");

  // Change the histogram to have empty bins below 0.3.  Will NOT have the TF1 in it.
  //TH1F* oldHisto;
  //rebinpT( lamReg_dNdpT_cor_data_7000, oldHisto, 4, 20 );
  TH1F* tmpReg7 = new TH1F("lamReg_dNdpT_scale7000", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  TH1F* tmpReg7sys = new TH1F("lamReg_dNdpT_scale7000_sys", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  for( int ndx = 1; ndx <= lamptNbins; ndx++ ) {
    if( ndx == 1 || ndx == 2 || ndx == 3) continue;
    int origNdx = ndx - 3;
    tmpReg7->SetBinContent(ndx, lamReg_dNdpT_cor_data_7000->GetBinContent(origNdx));
    tmpReg7->SetBinError(ndx, lamReg_dNdpT_cor_data_7000->GetBinError(origNdx));
    tmpReg7sys->SetBinContent(ndx, lamReg_dNdpT_cor_data_7000_sys->GetBinContent(origNdx));
    tmpReg7sys->SetBinError(ndx, lamReg_dNdpT_cor_data_7000_sys->GetBinError(origNdx));
  }
  tmpReg7->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmpReg7->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmpReg7->SetLineColor(4);
  tmpReg7->SetMarkerColor(4);
  tmpReg7->SetMarkerStyle(20);
  tmpReg7->SetMarkerSize(0.8);

  tmpReg7sys->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmpReg7sys->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmpReg7sys->SetLineColor(4);
  tmpReg7sys->SetMarkerColor(4);
  tmpReg7sys->SetMarkerStyle(20);
  tmpReg7sys->SetMarkerSize(0.8);

  TH1F* switchReg7000 = lamReg_dNdpT_cor_data_7000;
  lamReg_dNdpT_cor_data_7000 = tmpReg7;
  tmpReg7 = switchReg7000;
  TH1F* switchReg7000sys = lamReg_dNdpT_cor_data_7000_sys;
  lamReg_dNdpT_cor_data_7000_sys = tmpReg7sys;
  tmpReg7sys = switchReg7000;

  lamRegtsfunc7000->SetNpx(2000);
  lamRegtsfunc7000->SetLineColor(4);

  inFile = "plots/lamRegpt_cor_data_tuneD6T_900GeV_weight111_subgeant2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lamReg_dNdpT_cor_data_900 = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_scale");
  fin->Close();
  delete c1;
  delete fin;
  double yld_900_pt = lamReg_dNdpT_cor_data_900->Integral();
  lamReg_dNdpT_cor_data_900->Scale(1.0/(trig900*ptbinwid));
  double lamReg_yld_900_pt = calc_integral(lamReg_dNdpT_cor_data_900, 1.0);
  lamReg_dNdpT_cor_data_900->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  lamReg_dNdpT_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  lamReg_dNdpT_cor_data_900->SetLineColor(2);
  lamReg_dNdpT_cor_data_900->SetMarkerColor(2);
  lamReg_dNdpT_cor_data_900->SetMarkerStyle(21);
  lamReg_dNdpT_cor_data_900->SetMarkerSize(0.8);
  TH1F *lamReg_dNdpT_cor_data_900_sys = (TH1F*)lamReg_dNdpT_cor_data_900->Clone("lamReg_dNdpT_cor_data_700_sys");
  //cout << "lamReg pt 900 integral / event: " << lamReg_yld_900_pt << endl;
  for (int i = 1; i <= lamReg_dNdpT_cor_data_900_sys->GetXaxis()->GetNbins(); ++i) {
    double value = lamReg_dNdpT_cor_data_900_sys->GetBinContent(i);
    double error = lamReg_dNdpT_cor_data_900_sys->GetBinError(i);
    double syserr = 0.03*value;
    ptent[i] = sqrt(error*error + syserr*syserr);
    //cout << lamReg_dNdpT_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i] << endl;
  }
  lamReg_dNdpT_cor_data_900_sys->SetError(ptent);
  lamReg_dNdpT_cor_data_900->Fit("lamRegtsfunc900","EM0");

  // Change the histogram to have empty bins below 0.3.  Will NOT have the TF1 in it.
  TH1F* tmpReg9 = new TH1F("lamReg_dNdpT_scale900", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  TH1F* tmpReg9sys = new TH1F("lamReg_dNdpT_scale900sys", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  for( int ndx = 1; ndx <= lamptNbins; ndx++ ) {
    if( ndx == 1 || ndx == 2 || ndx == 3) continue;
    int origNdx = ndx - 3;
    tmpReg9->SetBinContent(ndx, lamReg_dNdpT_cor_data_900->GetBinContent(origNdx));
    tmpReg9->SetBinError(ndx, lamReg_dNdpT_cor_data_900->GetBinError(origNdx));
    tmpReg9sys->SetBinContent(ndx, lamReg_dNdpT_cor_data_900_sys->GetBinContent(origNdx));
    tmpReg9sys->SetBinError(ndx, lamReg_dNdpT_cor_data_900_sys->GetBinError(origNdx));
  }
  tmpReg9->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmpReg9->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmpReg9->SetLineColor(2);
  tmpReg9->SetMarkerColor(2);
  tmpReg9->SetMarkerStyle(21);
  tmpReg9->SetMarkerSize(0.8);

  tmpReg9sys->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmpReg9sys->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmpReg9sys->SetLineColor(2);
  tmpReg9sys->SetMarkerColor(2);
  tmpReg9sys->SetMarkerStyle(21);
  tmpReg9sys->SetMarkerSize(0.8);

  TH1F* switchReg900 = lamReg_dNdpT_cor_data_900;
  lamReg_dNdpT_cor_data_900 = tmpReg9;
  tmpReg9 = switchReg900;

  TH1F* switchReg900sys = lamReg_dNdpT_cor_data_900_sys;
  lamReg_dNdpT_cor_data_900_sys = tmpReg9sys;
  tmpReg9sys = switchReg900sys;
  //rebinpT( lamReg_dNdpT_cor_data_900, oldHisto, 2, 21 );



  lamRegtsfunc900->SetNpx(2000);
  lamRegtsfunc900->SetLineColor(2);

  TLegend *leg = new TLegend(0.50,0.70,0.85,0.80);
  leg->SetTextSize(0.042);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(lamReg_dNdpT_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
  leg->AddEntry(lamReg_dNdpT_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lamReg_dNdpT_cor_data_7000->Draw("E1");
  //oldHisto->Draw("E1");

  //c1->SaveAs("test.png");
  lamReg_dNdpT_cor_data_7000_sys->Draw("same");
  lamReg_dNdpT_cor_data_7000->Draw("E1same");
  lamRegtsfunc7000->Draw("same"); 
  lamReg_dNdpT_cor_data_900_sys->Draw("same");
  lamReg_dNdpT_cor_data_900->Draw("E1same");
  lamRegtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFileReg_pt_evt + ".root");
  c1->SaveAs(outFileReg_pt_evt + ".eps");
  c1->SaveAs(outFileReg_pt_evt + ".C");
  c1->SaveAs(outFileReg_pt_evt + ".png");

  delete c1;


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lamReg_dNdpT_cor_data_7000->SetMaximum(1.0);
  lamReg_dNdpT_cor_data_7000->SetMinimum(0.000001);
  lamReg_dNdpT_cor_data_7000->Draw("E1");
  lamReg_dNdpT_cor_data_7000_sys->Draw("same");
  lamReg_dNdpT_cor_data_7000->Draw("E1same");
  lamRegtsfunc7000->Draw("same"); 
  lamReg_dNdpT_cor_data_900_sys->Draw("same");
  lamReg_dNdpT_cor_data_900->Draw("E1same");
  lamRegtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  gPad->SetLogy(kTRUE);
  c1->SaveAs(outFileReg_pt_evt + "_log.root");
  c1->SaveAs(outFileReg_pt_evt + "_log.eps");
  c1->SaveAs(outFileReg_pt_evt + "_log.C");
  c1->SaveAs(outFileReg_pt_evt + "_log.png");
  gPad->SetLogy(kFALSE);

  delete c1;

  //return;

  // Switch histograms again
  switchReg7000 = lamReg_dNdpT_cor_data_7000;
  lamReg_dNdpT_cor_data_7000 = tmpReg7;
  tmpReg7 = switchReg7000;

  switchReg900 = lamReg_dNdpT_cor_data_900;
  lamReg_dNdpT_cor_data_900 = tmpReg9;
  tmpReg9 = switchReg900;

  // Pt for LamReg per luminosity

  lamReg_dNdpT_cor_data_7000->Scale(trig7000/lum7000);
  double lamReg_xsec_7000_pt = calc_integral(lamReg_dNdpT_cor_data_7000, 1.0);
  //cout << "lamReg pt 7000 xsec: " << lamReg_xsec_7000_pt << endl;
  lamReg_dNdpT_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dp_{T}  [ mb / (GeV/c)]");
  lamReg_dNdpT_cor_data_7000_sys->Scale(trig7000/lum7000);
  lamReg_dNdpT_cor_data_7000->Fit("lamRegtsfunc7000","EM0");
  switchReg7000 = lamReg_dNdpT_cor_data_7000;
  lamReg_dNdpT_cor_data_7000 = tmpReg7;
  tmpReg7 = switchReg7000;
  lamReg_dNdpT_cor_data_7000->Scale(trig7000/lum7000);
  lamRegtsfunc7000->SetNpx(2000);
  lamRegtsfunc7000->SetLineColor(4);

  lamReg_dNdpT_cor_data_900->Scale(trig900/lum900);
  double lamReg_xsec_900_pt = calc_integral(lamReg_dNdpT_cor_data_900, 1.0);
  //cout << "lamReg pt 900 xsec: " << lamReg_xsec_900_pt << endl;
  lamReg_dNdpT_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dp_{T}  [ mb / (GeV/c)]");
  lamReg_dNdpT_cor_data_900_sys->Scale(trig900/lum900);
  //cout << "pt 900\n";
  lamReg_dNdpT_cor_data_900->Fit("lamRegtsfunc900","EM0");
  switchReg900 = lamReg_dNdpT_cor_data_900;
  lamReg_dNdpT_cor_data_900 = tmpReg9;
  tmpReg9 = switchReg900;
  lamReg_dNdpT_cor_data_900->Scale(trig900/lum900);
  lamRegtsfunc900->SetNpx(2000);
  lamRegtsfunc900->SetLineColor(2);

  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lamReg_dNdpT_cor_data_7000->Draw("E1");
  lamReg_dNdpT_cor_data_7000_sys->Draw("same");
  lamReg_dNdpT_cor_data_7000->Draw("E1same");
  lamRegtsfunc7000->Draw("same"); 
  lamReg_dNdpT_cor_data_900_sys->Draw("same");
  lamReg_dNdpT_cor_data_900->Draw("E1same");
  lamRegtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFileReg_pt_lum + ".root");
  c1->SaveAs(outFileReg_pt_lum + ".eps");
  c1->SaveAs(outFileReg_pt_lum + ".C");
  c1->SaveAs(outFileReg_pt_lum + ".png");

  delete c1;


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lamReg_dNdpT_cor_data_7000->SetMaximum(70.);
  lamReg_dNdpT_cor_data_7000->SetMinimum(0.00007);
  lamReg_dNdpT_cor_data_7000_sys->SetMaximum(70.0);
  lamReg_dNdpT_cor_data_7000_sys->SetMinimum(0.00007);
  lamReg_dNdpT_cor_data_7000->Draw("E1");
  lamReg_dNdpT_cor_data_7000_sys->Draw("same");
  lamReg_dNdpT_cor_data_7000->Draw("E1same");
  lamRegtsfunc7000->Draw("same"); 
  lamReg_dNdpT_cor_data_900_sys->Draw("same");
  lamReg_dNdpT_cor_data_900->Draw("E1same");
  lamRegtsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  gPad->SetLogy(kTRUE);
  c1->SaveAs(outFileReg_pt_lum + "_log.root");
  c1->SaveAs(outFileReg_pt_lum + "_log.eps");
  c1->SaveAs(outFileReg_pt_lum + "_log.C");
  c1->SaveAs(outFileReg_pt_lum + "_log.png");
  gPad->SetLogy(kFALSE);

  delete c1;


   //---------------

  // Pt for Lambda(antiparticle only) per event

  inFile = "plots/lamBarpt_cor_data_pythia8_7TeV_weight102_subgeant2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lamBar_dNdpT_cor_data_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_scale");
  fin->Close();
  delete c1;
  delete fin;
  lamBar_dNdpT_cor_data_7000->Scale(1.0/(trig7000*ptbinwid));
  double lamBar_yld_7000_pt = calc_integral(lamBar_dNdpT_cor_data_7000, 1.0);
  lamBar_dNdpT_cor_data_7000->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  lamBar_dNdpT_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  lamBar_dNdpT_cor_data_7000->SetLineColor(4);
  lamBar_dNdpT_cor_data_7000->SetMarkerColor(4);
  lamBar_dNdpT_cor_data_7000->SetMarkerStyle(20);
  lamBar_dNdpT_cor_data_7000->SetMarkerSize(0.8);
  TH1F *lamBar_dNdpT_cor_data_7000_sys = (TH1F*)lamBar_dNdpT_cor_data_7000->Clone("lamBar_dNdpT_cor_data_700_sys");
  //cout << "lamBar pt 7000 integral / event: " << lamBar_yld_7000_pt << endl;
  for (int i = 1; i <= lamBar_dNdpT_cor_data_7000_sys->GetXaxis()->GetNbins(); ++i) {
    double value = lamBar_dNdpT_cor_data_7000_sys->GetBinContent(i);
    double error = lamBar_dNdpT_cor_data_7000_sys->GetBinError(i);
    double syserr = 0.03*value;
    ptent[i] = sqrt(error*error + syserr*syserr);
    //cout << lamBar_dNdpT_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i] << endl;
  }
  lamBar_dNdpT_cor_data_7000_sys->SetError(ptent);
  lamBar_dNdpT_cor_data_7000->Fit("lamBartsfunc7000","EM0");

  // Change the histogram to have empty bins below 0.3.  Will NOT have the TF1 in it.
  //TH1F* oldHisto;
  //rebinpT( lamBar_dNdpT_cor_data_7000, oldHisto, 4, 20 );
  TH1F* tmpBar7 = new TH1F("lamBar_dNdpT_scale7000", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  TH1F* tmpBar7sys = new TH1F("lamBar_dNdpT_scale7000_sys", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  for( int ndx = 1; ndx <= lamptNbins; ndx++ ) {
    if( ndx == 1 || ndx == 2 || ndx == 3) continue;
    int origNdx = ndx - 3;
    tmpBar7->SetBinContent(ndx, lamBar_dNdpT_cor_data_7000->GetBinContent(origNdx));
    tmpBar7->SetBinError(ndx, lamBar_dNdpT_cor_data_7000->GetBinError(origNdx));
    tmpBar7sys->SetBinContent(ndx, lamBar_dNdpT_cor_data_7000_sys->GetBinContent(origNdx));
    tmpBar7sys->SetBinError(ndx, lamBar_dNdpT_cor_data_7000_sys->GetBinError(origNdx));
  }
  tmpBar7->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmpBar7->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmpBar7->SetLineColor(4);
  tmpBar7->SetMarkerColor(4);
  tmpBar7->SetMarkerStyle(20);
  tmpBar7->SetMarkerSize(0.8);

  tmpBar7sys->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmpBar7sys->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmpBar7sys->SetLineColor(4);
  tmpBar7sys->SetMarkerColor(4);
  tmpBar7sys->SetMarkerStyle(20);
  tmpBar7sys->SetMarkerSize(0.8);

  TH1F* switchBar7000 = lamBar_dNdpT_cor_data_7000;
  lamBar_dNdpT_cor_data_7000 = tmpBar7;
  tmpBar7 = switchBar7000;
  TH1F* switchBar7000sys = lamBar_dNdpT_cor_data_7000_sys;
  lamBar_dNdpT_cor_data_7000_sys = tmpBar7sys;
  tmpBar7sys = switchBar7000;

  lamBartsfunc7000->SetNpx(2000);
  lamBartsfunc7000->SetLineColor(4);

  inFile = "plots/lamBarpt_cor_data_tuneD6T_900GeV_weight111_subgeant2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lamBar_dNdpT_cor_data_900 = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_scale");
  fin->Close();
  delete c1;
  delete fin;
  double yld_900_pt = lamBar_dNdpT_cor_data_900->Integral();
  lamBar_dNdpT_cor_data_900->Scale(1.0/(trig900*ptbinwid));
  double lamBar_yld_900_pt = calc_integral(lamBar_dNdpT_cor_data_900, 1.0);
  lamBar_dNdpT_cor_data_900->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  lamBar_dNdpT_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  lamBar_dNdpT_cor_data_900->SetLineColor(2);
  lamBar_dNdpT_cor_data_900->SetMarkerColor(2);
  lamBar_dNdpT_cor_data_900->SetMarkerStyle(21);
  lamBar_dNdpT_cor_data_900->SetMarkerSize(0.8);
  TH1F *lamBar_dNdpT_cor_data_900_sys = (TH1F*)lamBar_dNdpT_cor_data_900->Clone("lamBar_dNdpT_cor_data_700_sys");
  //cout << "lamBar pt 900 integral / event: " << lamBar_yld_900_pt << endl;
  for (int i = 1; i <= lamBar_dNdpT_cor_data_900_sys->GetXaxis()->GetNbins(); ++i) {
    double value = lamBar_dNdpT_cor_data_900_sys->GetBinContent(i);
    double error = lamBar_dNdpT_cor_data_900_sys->GetBinError(i);
    double syserr = 0.03*value;
    ptent[i] = sqrt(error*error + syserr*syserr);
    //cout << lamBar_dNdpT_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i] << endl;
  }
  lamBar_dNdpT_cor_data_900_sys->SetError(ptent);
  lamBar_dNdpT_cor_data_900->Fit("lamBartsfunc900","EM0");

  // Change the histogram to have empty bins below 0.3.  Will NOT have the TF1 in it.
  TH1F* tmpBar9 = new TH1F("lamBar_dNdpT_scale900", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  TH1F* tmpBar9sys = new TH1F("lamBar_dNdpT_scale900sys", "dN/dp_{T}", lamptNbins, &lamptBinMinima[0]);
  for( int ndx = 1; ndx <= lamptNbins; ndx++ ) {
    if( ndx == 1 || ndx == 2 || ndx == 3) continue;
    int origNdx = ndx - 3;
    tmpBar9->SetBinContent(ndx, lamBar_dNdpT_cor_data_900->GetBinContent(origNdx));
    tmpBar9->SetBinError(ndx, lamBar_dNdpT_cor_data_900->GetBinError(origNdx));
    tmpBar9sys->SetBinContent(ndx, lamBar_dNdpT_cor_data_900_sys->GetBinContent(origNdx));
    tmpBar9sys->SetBinError(ndx, lamBar_dNdpT_cor_data_900_sys->GetBinError(origNdx));
  }
  tmpBar9->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmpBar9->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmpBar9->SetLineColor(2);
  tmpBar9->SetMarkerColor(2);
  tmpBar9->SetMarkerStyle(21);
  tmpBar9->SetMarkerSize(0.8);

  tmpBar9sys->GetXaxis()->SetTitle("#Lambda^{0}  p_{T}  [GeV/c]");
  tmpBar9sys->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dp_{T}  (GeV/c)^{-1}");
  tmpBar9sys->SetLineColor(2);
  tmpBar9sys->SetMarkerColor(2);
  tmpBar9sys->SetMarkerStyle(21);
  tmpBar9sys->SetMarkerSize(0.8);

  TH1F* switchBar900 = lamBar_dNdpT_cor_data_900;
  lamBar_dNdpT_cor_data_900 = tmpBar9;
  tmpBar9 = switchBar900;

  TH1F* switchBar900sys = lamBar_dNdpT_cor_data_900_sys;
  lamBar_dNdpT_cor_data_900_sys = tmpBar9sys;
  tmpBar9sys = switchBar900sys;
  //rebinpT( lamBar_dNdpT_cor_data_900, oldHisto, 2, 21 );



  lamBartsfunc900->SetNpx(2000);
  lamBartsfunc900->SetLineColor(2);

  TLegend *leg = new TLegend(0.50,0.70,0.85,0.80);
  leg->SetTextSize(0.042);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(lamBar_dNdpT_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
  leg->AddEntry(lamBar_dNdpT_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lamBar_dNdpT_cor_data_7000->Draw("E1");
  //oldHisto->Draw("E1");

  //c1->SaveAs("test.png");
  lamBar_dNdpT_cor_data_7000_sys->Draw("same");
  lamBar_dNdpT_cor_data_7000->Draw("E1same");
  lamBartsfunc7000->Draw("same"); 
  lamBar_dNdpT_cor_data_900_sys->Draw("same");
  lamBar_dNdpT_cor_data_900->Draw("E1same");
  lamBartsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFileBar_pt_evt + ".root");
  c1->SaveAs(outFileBar_pt_evt + ".eps");
  c1->SaveAs(outFileBar_pt_evt + ".C");
  c1->SaveAs(outFileBar_pt_evt + ".png");

  delete c1;


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lamBar_dNdpT_cor_data_7000->SetMaximum(1.0);
  lamBar_dNdpT_cor_data_7000->SetMinimum(0.000001);
  lamBar_dNdpT_cor_data_7000->Draw("E1");
  lamBar_dNdpT_cor_data_7000_sys->Draw("same");
  lamBar_dNdpT_cor_data_7000->Draw("E1same");
  lamBartsfunc7000->Draw("same"); 
  lamBar_dNdpT_cor_data_900_sys->Draw("same");
  lamBar_dNdpT_cor_data_900->Draw("E1same");
  lamBartsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  gPad->SetLogy(kTRUE);
  c1->SaveAs(outFileBar_pt_evt + "_log.root");
  c1->SaveAs(outFileBar_pt_evt + "_log.eps");
  c1->SaveAs(outFileBar_pt_evt + "_log.C");
  c1->SaveAs(outFileBar_pt_evt + "_log.png");
  gPad->SetLogy(kFALSE);

  delete c1;

  //return;

  // Switch histograms again
  switchBar7000 = lamBar_dNdpT_cor_data_7000;
  lamBar_dNdpT_cor_data_7000 = tmpBar7;
  tmpBar7 = switchBar7000;

  switchBar900 = lamBar_dNdpT_cor_data_900;
  lamBar_dNdpT_cor_data_900 = tmpBar9;
  tmpBar9 = switchBar900;

  // Pt for LamBar per luminosity

  lamBar_dNdpT_cor_data_7000->Scale(trig7000/lum7000);
  double lamBar_xsec_7000_pt = calc_integral(lamBar_dNdpT_cor_data_7000, 1.0);
  //cout << "lamBar pt 7000 xsec: " << lamBar_xsec_7000_pt << endl;
  lamBar_dNdpT_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dp_{T}  [ mb / (GeV/c)]");
  lamBar_dNdpT_cor_data_7000_sys->Scale(trig7000/lum7000);
  lamBar_dNdpT_cor_data_7000->Fit("lamBartsfunc7000","EM0");
  switchBar7000 = lamBar_dNdpT_cor_data_7000;
  lamBar_dNdpT_cor_data_7000 = tmpBar7;
  tmpBar7 = switchBar7000;
  lamBar_dNdpT_cor_data_7000->Scale(trig7000/lum7000);
  lamBartsfunc7000->SetNpx(2000);
  lamBartsfunc7000->SetLineColor(4);

  lamBar_dNdpT_cor_data_900->Scale(trig900/lum900);
  double lamBar_xsec_900_pt = calc_integral(lamBar_dNdpT_cor_data_900, 1.0);
  //cout << "lamBar pt 900 xsec: " << lamBar_xsec_900_pt << endl;
  lamBar_dNdpT_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dp_{T}  [ mb / (GeV/c)]");
  lamBar_dNdpT_cor_data_900_sys->Scale(trig900/lum900);
  //cout << "pt 900\n";
  lamBar_dNdpT_cor_data_900->Fit("lamBartsfunc900","EM0");
  switchBar900 = lamBar_dNdpT_cor_data_900;
  lamBar_dNdpT_cor_data_900 = tmpBar9;
  tmpBar9 = switchBar900;
  lamBar_dNdpT_cor_data_900->Scale(trig900/lum900);
  lamBartsfunc900->SetNpx(2000);
  lamBartsfunc900->SetLineColor(2);

  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lamBar_dNdpT_cor_data_7000->Draw("E1");
  lamBar_dNdpT_cor_data_7000_sys->Draw("same");
  lamBar_dNdpT_cor_data_7000->Draw("E1same");
  lamBartsfunc7000->Draw("same"); 
  lamBar_dNdpT_cor_data_900_sys->Draw("same");
  lamBar_dNdpT_cor_data_900->Draw("E1same");
  lamBartsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFileBar_pt_lum + ".root");
  c1->SaveAs(outFileBar_pt_lum + ".eps");
  c1->SaveAs(outFileBar_pt_lum + ".C");
  c1->SaveAs(outFileBar_pt_lum + ".png");

  delete c1;


  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");

  lamBar_dNdpT_cor_data_7000->SetMaximum(70.);
  lamBar_dNdpT_cor_data_7000->SetMinimum(0.00007);
  lamBar_dNdpT_cor_data_7000_sys->SetMaximum(70.0);
  lamBar_dNdpT_cor_data_7000_sys->SetMinimum(0.00007);
  lamBar_dNdpT_cor_data_7000->Draw("E1");
  lamBar_dNdpT_cor_data_7000_sys->Draw("same");
  lamBar_dNdpT_cor_data_7000->Draw("E1same");
  lamBartsfunc7000->Draw("same"); 
  lamBar_dNdpT_cor_data_900_sys->Draw("same");
  lamBar_dNdpT_cor_data_900->Draw("E1same");
  lamBartsfunc900->Draw("same"); 
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  gPad->SetLogy(kTRUE);
  c1->SaveAs(outFileBar_pt_lum + "_log.root");
  c1->SaveAs(outFileBar_pt_lum + "_log.eps");
  c1->SaveAs(outFileBar_pt_lum + "_log.C");
  c1->SaveAs(outFileBar_pt_lum + "_log.png");
  gPad->SetLogy(kFALSE);

  delete c1;







  // y for Lambda (inclusive)

  inFile = "plots/lamy_cor_data_pythia8_7TeV_weight102_subgeant2.root";
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
  TH1F *lam_dNdy_cor_data_7000_sys = (TH1F*)lam_dNdy_cor_data_7000->Clone("lam_dNdy_cor_data_700_sys");
  //cout << "lam y 7000 integral / event: " << lam_yld_7000_y << endl;
  for (int i = 1; i <= 23; ++i) {
    double value = lam_dNdy_cor_data_7000_sys->GetBinContent(i);
    double error = lam_dNdy_cor_data_7000_sys->GetBinError(i);
    double syserr = 0.03*value;
    yent[i] = sqrt(error*error + syserr*syserr);
    cout << lam_dNdy_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
  }
  lam_dNdy_cor_data_7000_sys->SetError(yent);

  inFile = "plots/lamy_cor_data_tuneD6T_900GeV_weight111_subgeant2.root";
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
  TH1F *lam_dNdy_cor_data_900_sys = (TH1F*)lam_dNdy_cor_data_900->Clone("lam_dNdy_cor_data_700_sys");
  //cout << "lam y 900 integral / event: " << lam_yld_900_y << endl;
  for (int i = 1; i <= 23; ++i) {
    double value = lam_dNdy_cor_data_900_sys->GetBinContent(i);
    double error = lam_dNdy_cor_data_900_sys->GetBinError(i);
    double syserr = 0.03*value;
    yent[i] = sqrt(error*error + syserr*syserr);
    //cout << lam_dNdy_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
  }
  lam_dNdy_cor_data_900_sys->SetError(yent);

  TLegend *leg = new TLegend(0.50,0.25,0.85,0.35);
  leg->SetTextSize(0.042);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(lam_dNdy_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
  leg->AddEntry(lam_dNdy_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


  TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

  lam_dNdy_cor_data_7000_sys->Draw();
  lam_dNdy_cor_data_7000->Draw("E1same");
  lam_dNdy_cor_data_900_sys->Draw("same");
  lam_dNdy_cor_data_900->Draw("E1same");
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFile_y_evt + ".root");
  c1->SaveAs(outFile_y_evt + ".eps");
  c1->SaveAs(outFile_y_evt + ".C");
  c1->SaveAs(outFile_y_evt + ".png");

  delete c1;



  // y for Lambda per luminosity

  lam_dNdy_cor_data_7000->Scale(trig7000/lum7000);
  double lam_xsec_7000_y = calc_integral(lam_dNdy_cor_data_7000, 1.0);
  //cout << "lam y 7000 xsec: " << lam_xsec_7000_y << endl;
  lam_dNdy_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dy  [mb]");
  lam_dNdy_cor_data_7000_sys->Scale(trig7000/lum7000);

  lam_dNdy_cor_data_900->Scale(trig900/lum900);
  double lam_xsec_900_y = calc_integral(lam_dNdy_cor_data_900, 1.0);
  //cout << "lam y 900 xsec: " << lam_xsec_900_y << endl;
  lam_dNdy_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dy  [mb]");
  lam_dNdy_cor_data_900_sys->Scale(trig900/lum900);
  //cout << "y 900\n";

  TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

  lam_dNdy_cor_data_7000_sys->SetMaximum(12.0);
  lam_dNdy_cor_data_7000_sys->SetMinimum(0.);
  lam_dNdy_cor_data_7000_sys->Draw();
  lam_dNdy_cor_data_7000->Draw("E1same");
  lam_dNdy_cor_data_900_sys->Draw("same");
  lam_dNdy_cor_data_900->Draw("E1same");
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFile_y_lum + ".root");
  c1->SaveAs(outFile_y_lum + ".eps");
  c1->SaveAs(outFile_y_lum + ".C");
  c1->SaveAs(outFile_y_lum + ".png");

  delete c1;

  //---------------

  // y for Lambda (particle)

  inFile = "plots/lamRegy_cor_data_pythia8_7TeV_weight102_subgeant2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lamReg_dNdy_cor_data_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdy");
  fin->Close();
  delete c1;
  delete fin;
  double yldReg_7000_y = lamReg_dNdy_cor_data_7000->Integral();
  lamReg_dNdy_cor_data_7000->Scale(1.0/(trig7000*ybinwid));
  double lamReg_yld_7000_y = calc_integral(lamReg_dNdy_cor_data_7000, 1.0);
  lamReg_dNdy_cor_data_7000->GetXaxis()->SetTitle("#Lambda^{0}  rapidity");
  lamReg_dNdy_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
  lamReg_dNdy_cor_data_7000->SetMaximum(0.22);
  lamReg_dNdy_cor_data_7000->SetMinimum(0.);
  lamReg_dNdy_cor_data_7000->SetLineColor(4);
  lamReg_dNdy_cor_data_7000->SetMarkerColor(4);
  lamReg_dNdy_cor_data_7000->SetMarkerStyle(20);
  lamReg_dNdy_cor_data_7000->SetMarkerSize(0.8);
  TH1F *lamReg_dNdy_cor_data_7000_sys = (TH1F*)lamReg_dNdy_cor_data_7000->Clone("lamReg_dNdy_cor_data_700_sys");
  //cout << "lamReg y 7000 integral / event: " << lamReg_yld_7000_y << endl;
  for (int i = 1; i <= 23; ++i) {
    double value = lamReg_dNdy_cor_data_7000_sys->GetBinContent(i);
    double error = lamReg_dNdy_cor_data_7000_sys->GetBinError(i);
    double syserr = 0.03*value;
    yent[i] = sqrt(error*error + syserr*syserr);
    //cout << lamReg_dNdy_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
  }
  lamReg_dNdy_cor_data_7000_sys->SetError(yent);

  inFile = "plots/lamRegy_cor_data_tuneD6T_900GeV_weight111_subgeant2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lamReg_dNdy_cor_data_900 = (TH1F*)c1->GetPrimitive("lamReg_dNdy");
  fin->Close();
  delete c1;
  delete fin;
  double yld_900_y = lamReg_dNdy_cor_data_900->Integral();
  lamReg_dNdy_cor_data_900->Scale(1.0/(trig900*ybinwid));
  double lamReg_yld_900_y = calc_integral(lamReg_dNdy_cor_data_900, 1.0);
  lamReg_dNdy_cor_data_900->GetXaxis()->SetTitle("#Lambda^{0}  rapidity");
  lamReg_dNdy_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
  lamReg_dNdy_cor_data_900->SetMaximum(0.22);
  lamReg_dNdy_cor_data_900->SetMinimum(0.);
  lamReg_dNdy_cor_data_900->SetLineColor(2);
  lamReg_dNdy_cor_data_900->SetMarkerColor(2);
  lamReg_dNdy_cor_data_900->SetMarkerStyle(21);
  lamReg_dNdy_cor_data_900->SetMarkerSize(0.8);
  TH1F *lamReg_dNdy_cor_data_900_sys = (TH1F*)lamReg_dNdy_cor_data_900->Clone("lamReg_dNdy_cor_data_700_sys");
  //cout << "lamReg y 900 integral / event: " << lamReg_yld_900_y << endl;
  for (int i = 1; i <= 23; ++i) {
    double value = lamReg_dNdy_cor_data_900_sys->GetBinContent(i);
    double error = lamReg_dNdy_cor_data_900_sys->GetBinError(i);
    double syserr = 0.03*value;
    yent[i] = sqrt(error*error + syserr*syserr);
    //cout << lamReg_dNdy_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
  }
  lamReg_dNdy_cor_data_900_sys->SetError(yent);

  TLegend *leg = new TLegend(0.50,0.25,0.85,0.35);
  leg->SetTextSize(0.042);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(lamReg_dNdy_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
  leg->AddEntry(lamReg_dNdy_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


  TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

  lamReg_dNdy_cor_data_7000_sys->Draw();
  lamReg_dNdy_cor_data_7000->Draw("E1same");
  lamReg_dNdy_cor_data_900_sys->Draw("same");
  lamReg_dNdy_cor_data_900->Draw("E1same");
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFileReg_y_evt + ".root");
  c1->SaveAs(outFileReg_y_evt + ".eps");
  c1->SaveAs(outFileReg_y_evt + ".C");
  c1->SaveAs(outFileReg_y_evt + ".png");

  delete c1;



  // y for Lambda per luminosity

  lamReg_dNdy_cor_data_7000->Scale(trig7000/lum7000);
  double lamReg_xsec_7000_y = calc_integral(lamReg_dNdy_cor_data_7000, 1.0);
  //cout << "lamReg y 7000 xsec: " << lamReg_xsec_7000_y << endl;
  lamReg_dNdy_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dy  [mb]");
  lamReg_dNdy_cor_data_7000_sys->Scale(trig7000/lum7000);

  lamReg_dNdy_cor_data_900->Scale(trig900/lum900);
  double lamReg_xsec_900_y = calc_integral(lamReg_dNdy_cor_data_900, 1.0);
  //cout << "lamReg y 900 xsec: " << lamReg_xsec_900_y << endl;
  lamReg_dNdy_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dy  [mb]");
  lamReg_dNdy_cor_data_900_sys->Scale(trig900/lum900);
  //cout << "y 900\n";

  lamReg_dNdy_cor_data_7000_sys->SetMaximum(12.0);
  lamReg_dNdy_cor_data_7000_sys->SetMinimum(0.);
  c1 = new TCanvas();
  lamReg_dNdy_cor_data_7000_sys->Draw("E1");
  c1->SaveAs("numeratorTest.png");
  lam_dNdy_cor_data_7000_sys->Draw("E1");
  c1->SaveAs("denominatorTest.png");
  delete c1;

  TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

  lamReg_dNdy_cor_data_7000_sys->SetMaximum(12.0);
  lamReg_dNdy_cor_data_7000_sys->SetMinimum(0.);
  lamReg_dNdy_cor_data_7000_sys->Draw();
  lamReg_dNdy_cor_data_7000->Draw("E1same");
  lamReg_dNdy_cor_data_900_sys->Draw("same");
  lamReg_dNdy_cor_data_900->Draw("E1same");
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFileReg_y_lum + ".root");
  c1->SaveAs(outFileReg_y_lum + ".eps");
  c1->SaveAs(outFileReg_y_lum + ".C");
  c1->SaveAs(outFileReg_y_lum + ".png");

  delete c1;



  //----------------

  // y for Lambda (antiparticle)

  inFile = "plots/lamBary_cor_data_pythia8_7TeV_weight102_subgeant2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lamBar_dNdy_cor_data_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdy");
  fin->Close();
  delete c1;
  delete fin;
  double yldBar_7000_y = lamBar_dNdy_cor_data_7000->Integral();
  lamBar_dNdy_cor_data_7000->Scale(1.0/(trig7000*ybinwid));
  double lamBar_yld_7000_y = calc_integral(lamBar_dNdy_cor_data_7000, 1.0);
  lamBar_dNdy_cor_data_7000->GetXaxis()->SetTitle("#Lambda^{0}  rapidity");
  lamBar_dNdy_cor_data_7000->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
  lamBar_dNdy_cor_data_7000->SetMaximum(0.22);
  lamBar_dNdy_cor_data_7000->SetMinimum(0.);
  lamBar_dNdy_cor_data_7000->SetLineColor(4);
  lamBar_dNdy_cor_data_7000->SetMarkerColor(4);
  lamBar_dNdy_cor_data_7000->SetMarkerStyle(20);
  lamBar_dNdy_cor_data_7000->SetMarkerSize(0.8);
  TH1F *lamBar_dNdy_cor_data_7000_sys = (TH1F*)lamBar_dNdy_cor_data_7000->Clone("lamBar_dNdy_cor_data_700_sys");
  //cout << "lamBar y 7000 integral / event: " << lamBar_yld_7000_y << endl;
  for (int i = 1; i <= 23; ++i) {
    double value = lamBar_dNdy_cor_data_7000_sys->GetBinContent(i);
    double error = lamBar_dNdy_cor_data_7000_sys->GetBinError(i);
    double syserr = 0.03*value;
    yent[i] = sqrt(error*error + syserr*syserr);
    //cout << lamBar_dNdy_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
  }
  lamBar_dNdy_cor_data_7000_sys->SetError(yent);

  inFile = "plots/lamBary_cor_data_tuneD6T_900GeV_weight111_subgeant2.root";
  TFile* fin = TFile::Open(inFile, "READ");
  TCanvas *c1 = (TCanvas*)fin->Get("c1");
  TH1F* lamBar_dNdy_cor_data_900 = (TH1F*)c1->GetPrimitive("lamBar_dNdy");
  fin->Close();
  delete c1;
  delete fin;
  double yld_900_y = lamBar_dNdy_cor_data_900->Integral();
  lamBar_dNdy_cor_data_900->Scale(1.0/(trig900*ybinwid));
  double lamBar_yld_900_y = calc_integral(lamBar_dNdy_cor_data_900, 1.0);
  lamBar_dNdy_cor_data_900->GetXaxis()->SetTitle("#Lambda^{0}  rapidity");
  lamBar_dNdy_cor_data_900->GetYaxis()->SetTitle("(1/N_{NSD Event}) dN(pp #rightarrow #Lambda^{0}X) / dy");
  lamBar_dNdy_cor_data_900->SetMaximum(0.22);
  lamBar_dNdy_cor_data_900->SetMinimum(0.);
  lamBar_dNdy_cor_data_900->SetLineColor(2);
  lamBar_dNdy_cor_data_900->SetMarkerColor(2);
  lamBar_dNdy_cor_data_900->SetMarkerStyle(21);
  lamBar_dNdy_cor_data_900->SetMarkerSize(0.8);
  TH1F *lamBar_dNdy_cor_data_900_sys = (TH1F*)lamBar_dNdy_cor_data_900->Clone("lamBar_dNdy_cor_data_700_sys");
  //cout << "lamBar y 900 integral / event: " << lamBar_yld_900_y << endl;
  for (int i = 1; i <= 23; ++i) {
    double value = lamBar_dNdy_cor_data_900_sys->GetBinContent(i);
    double error = lamBar_dNdy_cor_data_900_sys->GetBinError(i);
    double syserr = 0.03*value;
    yent[i] = sqrt(error*error + syserr*syserr);
    //cout << lamBar_dNdy_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i] << endl;
  }
  lamBar_dNdy_cor_data_900_sys->SetError(yent);

  TLegend *leg = new TLegend(0.50,0.25,0.85,0.35);
  leg->SetTextSize(0.042);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(lamBar_dNdy_cor_data_7000, "  #sqrt{s} = 7 TeV" , "pl" );
  leg->AddEntry(lamBar_dNdy_cor_data_900, "  #sqrt{s} = 900 GeV" , "pl" );


  TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

  lamBar_dNdy_cor_data_7000_sys->Draw();
  lamBar_dNdy_cor_data_7000->Draw("E1same");
  lamBar_dNdy_cor_data_900_sys->Draw("same");
  lamBar_dNdy_cor_data_900->Draw("E1same");
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFileBar_y_evt + ".root");
  c1->SaveAs(outFileBar_y_evt + ".eps");
  c1->SaveAs(outFileBar_y_evt + ".C");
  c1->SaveAs(outFileBar_y_evt + ".png");

  delete c1;



  // y for Lambda per luminosity

  lamBar_dNdy_cor_data_7000->Scale(trig7000/lum7000);
  double lamBar_xsec_7000_y = calc_integral(lamBar_dNdy_cor_data_7000, 1.0);
  //cout << "lamBar y 7000 xsec: " << lamBar_xsec_7000_y << endl;
  lamBar_dNdy_cor_data_7000_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dy  [mb]");
  lamBar_dNdy_cor_data_7000_sys->Scale(trig7000/lum7000);

  lamBar_dNdy_cor_data_900->Scale(trig900/lum900);
  double lamBar_xsec_900_y = calc_integral(lamBar_dNdy_cor_data_900, 1.0);
  //cout << "lamBar y 900 xsec: " << lamBar_xsec_900_y << endl;
  lamBar_dNdy_cor_data_900_sys->GetYaxis()->SetTitle("d#sigma(pp #rightarrow #Lambda^{0}X) / dy  [mb]");
  lamBar_dNdy_cor_data_900_sys->Scale(trig900/lum900);
  //cout << "y 900\n";

  TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

  lamBar_dNdy_cor_data_7000_sys->SetMaximum(12.0);
  lamBar_dNdy_cor_data_7000_sys->SetMinimum(0.);
  lamBar_dNdy_cor_data_7000_sys->Draw();
  lamBar_dNdy_cor_data_7000->Draw("E1same");
  lamBar_dNdy_cor_data_900_sys->Draw("same");
  lamBar_dNdy_cor_data_900->Draw("E1same");
  leg->Draw("same");
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  c1->SaveAs(outFileBar_y_lum + ".root");
  c1->SaveAs(outFileBar_y_lum + ".eps");
  c1->SaveAs(outFileBar_y_lum + ".C");
  c1->SaveAs(outFileBar_y_lum + ".png");

  delete c1;

  //----------------

  // Get the generated plots for all MC tunes and uncorrected for pythia8, make asymmetry plots
  // Rapidity first.
  //  GENERATED:
  TFile* fin_y_py8_gen_lall 
    = TFile::Open("lam_pythia8_7000/plots/lam_dNdy_gen_pythia8_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_py8_gen_lall->Get("c1");
  TH1F* lam_dNdy_gen_pythia8_7000 = (TH1F*)c1->GetPrimitive("lam_dNdy_gen");
  fin_y_py8_gen_lall->Close();
  delete c1;

  TFile* fin_y_py8_gen_part
    = TFile::Open("lam_pythia8_7000/plots/lamReg_dNdy_gen_pythia8_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_py8_gen_part->Get("c1");
  TH1F* lamReg_dNdy_gen_pythia8_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdy_gen");
  fin_y_py8_gen_part->Close();
  delete c1;

  TFile* fin_y_py8_gen_anti
    = TFile::Open("lam_pythia8_7000/plots/lamBar_dNdy_gen_pythia8_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_py8_gen_anti->Get("c1");
  TH1F* lamBar_dNdy_gen_pythia8_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdy_gen");
  fin_y_py8_gen_anti->Close();
  delete c1;

  TH1F* lam_dNdy_gen_pythia8_7000_asym = (TH1F*)lamReg_dNdy_gen_pythia8_7000->Clone("lam_dNdy_gen_pythia8_7000_asym");
  double *lam_dNdy_gen_pythia8_7000_asym_err 
    = new double[lam_dNdy_gen_pythia8_7000_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdy_gen_pythia8_7000_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdy_gen_pythia8_7000->GetBinContent(ndx);
    double b = lamBar_dNdy_gen_pythia8_7000->GetBinContent(ndx);
    double aPlusb = lam_dNdy_gen_pythia8_7000->GetBinContent(ndx);
    double aErr = lamReg_dNdy_gen_pythia8_7000->GetBinError(ndx);
    double bErr = lamBar_dNdy_gen_pythia8_7000->GetBinError(ndx);
    double aPlusbErr = lam_dNdy_gen_pythia8_7000->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdy_gen_pythia8_7000_asym_err[ndx] = fError;
  }
  lam_dNdy_gen_pythia8_7000_asym->Add( lamBar_dNdy_gen_pythia8_7000, -1. );
  lam_dNdy_gen_pythia8_7000_asym->Divide( lam_dNdy_gen_pythia8_7000 );
  lam_dNdy_gen_pythia8_7000_asym->SetError(lam_dNdy_gen_pythia8_7000_asym_err);


  //lam_dNdy_gen_pythia8_7000_asym->SetLineColor(41);

  // tuneD6T, 7000 GeV:
  TFile* fin_y_d6t_gen_lall 
    = TFile::Open("lam_tuneD6T_7000/plots/lam_dNdy_gen_tuneD6T_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_d6t_gen_lall->Get("c1");
  TH1F* lam_dNdy_gen_tuneD6T_7000 = (TH1F*)c1->GetPrimitive("lam_dNdy_gen");
  fin_y_d6t_gen_lall->Close();
  delete c1;

  TFile* fin_y_d6t_gen_part
    = TFile::Open("lam_tuneD6T_7000/plots/lamReg_dNdy_gen_tuneD6T_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_d6t_gen_part->Get("c1");
  TH1F* lamReg_dNdy_gen_tuneD6T_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdy_gen");
  fin_y_d6t_gen_part->Close();
  delete c1;

  TFile* fin_y_d6t_gen_anti
    = TFile::Open("lam_tuneD6T_7000/plots/lamBar_dNdy_gen_tuneD6T_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_d6t_gen_anti->Get("c1");
  TH1F* lamBar_dNdy_gen_tuneD6T_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdy_gen");
  fin_y_d6t_gen_anti->Close();
  delete c1;

  TH1F* lam_dNdy_gen_tuneD6T_7000_asym = (TH1F*)lamReg_dNdy_gen_tuneD6T_7000->Clone("lam_dNdy_gen_tuneD6T_7000_asym");
  double *lam_dNdy_gen_tuneD6T_7000_asym_err 
    = new double[lam_dNdy_gen_tuneD6T_7000_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdy_gen_tuneD6T_7000_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdy_gen_tuneD6T_7000->GetBinContent(ndx);
    double b = lamBar_dNdy_gen_tuneD6T_7000->GetBinContent(ndx);
    double aPlusb = lam_dNdy_gen_tuneD6T_7000->GetBinContent(ndx);
    double aErr = lamReg_dNdy_gen_tuneD6T_7000->GetBinError(ndx);
    double bErr = lamBar_dNdy_gen_tuneD6T_7000->GetBinError(ndx);
    double aPlusbErr = lam_dNdy_gen_tuneD6T_7000->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdy_gen_tuneD6T_7000_asym_err[ndx] = fError;
  }
  lam_dNdy_gen_tuneD6T_7000_asym->Add( lamBar_dNdy_gen_tuneD6T_7000, -1. );
  lam_dNdy_gen_tuneD6T_7000_asym->Divide( lam_dNdy_gen_tuneD6T_7000 );
  lam_dNdy_gen_tuneD6T_7000_asym->SetError(lam_dNdy_gen_tuneD6T_7000_asym_err);
  //lam_dNdy_gen_tuneD6T_7000_asym->Add( lamBar_dNdy_gen_tuneD6T_7000, -1. );
  //lam_dNdy_gen_tuneD6T_7000_asym->Divide( lam_dNdy_gen_tuneD6T_7000 );
  //lam_dNdy_gen_tuneD6T_7000_asym->SetLineColor(42);

  // tuneD6T, 7000 GeV:
  TFile* fin_y_d6t900_gen_lall 
    = TFile::Open("lam_tuneD6T_900/plots/lam_dNdy_gen_tuneD6T_900GeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_d6t900_gen_lall->Get("c1");
  TH1F* lam_dNdy_gen_tuneD6T_900 = (TH1F*)c1->GetPrimitive("lam_dNdy_gen");
  fin_y_d6t900_gen_lall->Close();
  delete c1;

  TFile* fin_y_d6t900_gen_part
    = TFile::Open("lam_tuneD6T_900/plots/lamReg_dNdy_gen_tuneD6T_900GeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_d6t900_gen_part->Get("c1");
  TH1F* lamReg_dNdy_gen_tuneD6T_900 = (TH1F*)c1->GetPrimitive("lamReg_dNdy_gen");
  fin_y_d6t900_gen_part->Close();
  delete c1;

  TFile* fin_y_d6t900_gen_anti
    = TFile::Open("lam_tuneD6T_900/plots/lamBar_dNdy_gen_tuneD6T_900GeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_d6t900_gen_anti->Get("c1");
  TH1F* lamBar_dNdy_gen_tuneD6T_900 = (TH1F*)c1->GetPrimitive("lamBar_dNdy_gen");
  fin_y_d6t900_gen_anti->Close();
  delete c1;

  TH1F* lam_dNdy_gen_tuneD6T_900_asym = (TH1F*)lamReg_dNdy_gen_tuneD6T_900->Clone("lam_dNdy_gen_tuneD6T_900_asym");
  double *lam_dNdy_gen_tuneD6T_900_asym_err 
    = new double[lam_dNdy_gen_tuneD6T_900_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdy_gen_tuneD6T_900_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdy_gen_tuneD6T_900->GetBinContent(ndx);
    double b = lamBar_dNdy_gen_tuneD6T_900->GetBinContent(ndx);
    double aPlusb = lam_dNdy_gen_tuneD6T_900->GetBinContent(ndx);
    double aErr = lamReg_dNdy_gen_tuneD6T_900->GetBinError(ndx);
    double bErr = lamBar_dNdy_gen_tuneD6T_900->GetBinError(ndx);
    double aPlusbErr = lam_dNdy_gen_tuneD6T_900->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdy_gen_tuneD6T_900_asym_err[ndx] = fError;
  }
  lam_dNdy_gen_tuneD6T_900_asym->Add( lamBar_dNdy_gen_tuneD6T_900, -1. );
  lam_dNdy_gen_tuneD6T_900_asym->Divide( lam_dNdy_gen_tuneD6T_900 );
  lam_dNdy_gen_tuneD6T_900_asym->SetError(lam_dNdy_gen_tuneD6T_900_asym_err);
  //lam_dNdy_gen_tuneD6T_900_asym->Add( lamBar_dNdy_gen_tuneD6T_900, -1. );
  //lam_dNdy_gen_tuneD6T_900_asym->Divide( lam_dNdy_gen_tuneD6T_900 );
  //lam_dNdy_gen_tuneD6T_900_asym->SetLineColor(42);

  // P0:
  TFile* fin_y_p0_gen_lall 
    = TFile::Open("lam_tuneP0_7000/plots/lam_dNdy_gen_tuneP0_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_p0_gen_lall->Get("c1");
  TH1F* lam_dNdy_gen_tuneP0_7000 = (TH1F*)c1->GetPrimitive("lam_dNdy_gen");
  fin_y_p0_gen_lall->Close();
  delete c1;

  TFile* fin_y_p0_gen_part
    = TFile::Open("lam_tuneP0_7000/plots/lamReg_dNdy_gen_tuneP0_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_p0_gen_part->Get("c1");
  TH1F* lamReg_dNdy_gen_tuneP0_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdy_gen");
  fin_y_p0_gen_part->Close();
  delete c1;

  TFile* fin_y_p0_gen_anti
    = TFile::Open("lam_tuneP0_7000/plots/lamBar_dNdy_gen_tuneP0_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_y_p0_gen_anti->Get("c1");
  TH1F* lamBar_dNdy_gen_tuneP0_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdy_gen");
  fin_y_p0_gen_anti->Close();
  delete c1;

  TH1F* lam_dNdy_gen_tuneP0_7000_asym = (TH1F*)lamReg_dNdy_gen_tuneP0_7000->Clone("lam_dNdy_gen_tuneP0_7000_asym");
  double *lam_dNdy_gen_tuneP0_7000_asym_err 
    = new double[lam_dNdy_gen_tuneP0_7000_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdy_gen_tuneP0_7000_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdy_gen_tuneP0_7000->GetBinContent(ndx);
    double b = lamBar_dNdy_gen_tuneP0_7000->GetBinContent(ndx);
    double aPlusb = lam_dNdy_gen_tuneP0_7000->GetBinContent(ndx);
    double aErr = lamReg_dNdy_gen_tuneP0_7000->GetBinError(ndx);
    double bErr = lamBar_dNdy_gen_tuneP0_7000->GetBinError(ndx);
    double aPlusbErr = lam_dNdy_gen_tuneP0_7000->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdy_gen_tuneP0_7000_asym_err[ndx] = fError;
  }
  lam_dNdy_gen_tuneP0_7000_asym->Add( lamBar_dNdy_gen_tuneP0_7000, -1. );
  lam_dNdy_gen_tuneP0_7000_asym->Divide( lam_dNdy_gen_tuneP0_7000 );
  lam_dNdy_gen_tuneP0_7000_asym->SetError(lam_dNdy_gen_tuneP0_7000_asym_err);
  //lam_dNdy_gen_tuneP0_7000_asym->Add( lamBar_dNdy_gen_tuneP0_7000, -1. );
  //lam_dNdy_gen_tuneP0_7000_asym->Divide( lam_dNdy_gen_tuneP0_7000 );
  //lam_dNdy_gen_tuneP0_7000_asym->SetLineColor(31);
  //----

  cout << "At pT:" << endl;
  // pT:
  TFile* fin_pT_py8_gen_lall 
    = TFile::Open("lam_pythia8_7000/plots/lam_dNdpT_gen_pythia8_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_py8_gen_lall->Get("c1");
  TH1F* lam_dNdpT_gen_pythia8_7000 = (TH1F*)c1->GetPrimitive("lam_dNdpT_gen_scale");
  fin_pT_py8_gen_lall->Close();
  delete c1;

  TFile* fin_pT_py8_gen_part
    = TFile::Open("lam_pythia8_7000/plots/lamReg_dNdpT_gen_pythia8_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_py8_gen_part->Get("c1");
  TH1F* lamReg_dNdpT_gen_pythia8_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_gen_scale");
  fin_pT_py8_gen_part->Close();
  delete c1;

  TFile* fin_pT_py8_gen_anti
    = TFile::Open("lam_pythia8_7000/plots/lamBar_dNdpT_gen_pythia8_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_py8_gen_anti->Get("c1");
  TH1F* lamBar_dNdpT_gen_pythia8_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_gen_scale");
  fin_pT_py8_gen_anti->Close();
  delete c1;

  TH1F* lam_dNdpT_gen_pythia8_7000_asym = (TH1F*)lamReg_dNdpT_gen_pythia8_7000->Clone("lam_dNdpT_gen_pythia8_7000_asym");
  double *lam_dNdpT_gen_pythia8_7000_asym_err 
    = new double[lam_dNdpT_gen_pythia8_7000_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdpT_gen_pythia8_7000_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdpT_gen_pythia8_7000->GetBinContent(ndx);
    double b = lamBar_dNdpT_gen_pythia8_7000->GetBinContent(ndx);
    double aPlusb = lam_dNdpT_gen_pythia8_7000->GetBinContent(ndx);
    double aErr = lamReg_dNdpT_gen_pythia8_7000->GetBinError(ndx);
    double bErr = lamBar_dNdpT_gen_pythia8_7000->GetBinError(ndx);
    double aPlusbErr = lam_dNdpT_gen_pythia8_7000->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdpT_gen_pythia8_7000_asym_err[ndx] = fError;
  }
  lam_dNdpT_gen_pythia8_7000_asym->Add( lamBar_dNdpT_gen_pythia8_7000, -1. );
  lam_dNdpT_gen_pythia8_7000_asym->Divide( lam_dNdpT_gen_pythia8_7000 );
  lam_dNdpT_gen_pythia8_7000_asym->SetError(lam_dNdpT_gen_pythia8_7000_asym_err);
  //lam_dNdpT_gen_pythia8_7000_asym->Add( lamBar_dNdpT_gen_pythia8_7000, -1. );
  //lam_dNdpT_gen_pythia8_7000_asym->Divide( lam_dNdpT_gen_pythia8_7000 );
  //lam_dNdpT_gen_pythia8_7000_asym->SetLineColor(41);

  // tuneD6T, 7000 GeV:
  TFile* fin_pT_d6t_gen_lall 
    = TFile::Open("lam_tuneD6T_7000/plots/lam_dNdpT_gen_tuneD6T_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_d6t_gen_lall->Get("c1");
  TH1F* lam_dNdpT_gen_tuneD6T_7000 = (TH1F*)c1->GetPrimitive("lam_dNdpT_gen_scale");
  fin_pT_d6t_gen_lall->Close();
  delete c1;

  TFile* fin_pT_d6t_gen_part
    = TFile::Open("lam_tuneD6T_7000/plots/lamReg_dNdpT_gen_tuneD6T_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_d6t_gen_part->Get("c1");
  TH1F* lamReg_dNdpT_gen_tuneD6T_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_gen_scale");
  fin_pT_d6t_gen_part->Close();
  delete c1;

  TFile* fin_pT_d6t_gen_anti
    = TFile::Open("lam_tuneD6T_7000/plots/lamBar_dNdpT_gen_tuneD6T_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_d6t_gen_anti->Get("c1");
  TH1F* lamBar_dNdpT_gen_tuneD6T_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_gen_scale");
  fin_pT_d6t_gen_anti->Close();
  delete c1;

  TH1F* lam_dNdpT_gen_tuneD6T_7000_asym = (TH1F*)lamReg_dNdpT_gen_tuneD6T_7000->Clone("lam_dNdpT_gen_tuneD6T_7000_asym");
  double *lam_dNdpT_gen_tuneD6T_7000_asym_err 
    = new double[lam_dNdpT_gen_tuneD6T_7000_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdpT_gen_tuneD6T_7000_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdpT_gen_tuneD6T_7000->GetBinContent(ndx);
    double b = lamBar_dNdpT_gen_tuneD6T_7000->GetBinContent(ndx);
    double aPlusb = lam_dNdpT_gen_tuneD6T_7000->GetBinContent(ndx);
    double aErr = lamReg_dNdpT_gen_tuneD6T_7000->GetBinError(ndx);
    double bErr = lamBar_dNdpT_gen_tuneD6T_7000->GetBinError(ndx);
    double aPlusbErr = lam_dNdpT_gen_tuneD6T_7000->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdpT_gen_tuneD6T_7000_asym_err[ndx] = fError;
  }
  lam_dNdpT_gen_tuneD6T_7000_asym->Add( lamBar_dNdpT_gen_tuneD6T_7000, -1. );
  lam_dNdpT_gen_tuneD6T_7000_asym->Divide( lam_dNdpT_gen_tuneD6T_7000 );
  lam_dNdpT_gen_tuneD6T_7000_asym->SetError(lam_dNdpT_gen_tuneD6T_7000_asym_err);
  //lam_dNdpT_gen_tuneD6T_7000_asym->Add( lamBar_dNdpT_gen_tuneD6T_7000, -1. );
  //lam_dNdpT_gen_tuneD6T_7000_asym->Divide( lam_dNdpT_gen_tuneD6T_7000 );
  //lam_dNdpT_gen_tuneD6T_7000_asym->SetLineColor(42);

  // tuneD6T, 900 GeV:
  TFile* fin_pT_d6t900_gen_lall 
    = TFile::Open("lam_tuneD6T_900/plots/lam_dNdpT_gen_tuneD6T_900GeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_d6t900_gen_lall->Get("c1");
  TH1F* lam_dNdpT_gen_tuneD6T_900 = (TH1F*)c1->GetPrimitive("lam_dNdpT_gen_scale");
  fin_pT_d6t900_gen_lall->Close();
  delete c1;

  TFile* fin_pT_d6t900_gen_part
    = TFile::Open("lam_tuneD6T_900/plots/lamReg_dNdpT_gen_tuneD6T_900GeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_d6t900_gen_part->Get("c1");
  TH1F* lamReg_dNdpT_gen_tuneD6T_900 = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_gen_scale");
  fin_pT_d6t900_gen_part->Close();
  delete c1;

  TFile* fin_pT_d6t900_gen_anti
    = TFile::Open("lam_tuneD6T_900/plots/lamBar_dNdpT_gen_tuneD6T_900GeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_d6t900_gen_anti->Get("c1");
  TH1F* lamBar_dNdpT_gen_tuneD6T_900 = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_gen_scale");
  fin_pT_d6t900_gen_anti->Close();
  delete c1;

  TH1F* lam_dNdpT_gen_tuneD6T_900_asym = (TH1F*)lamReg_dNdpT_gen_tuneD6T_900->Clone("lam_dNdpT_gen_tuneD6T_900_asym");
  double *lam_dNdpT_gen_tuneD6T_900_asym_err 
    = new double[lam_dNdpT_gen_tuneD6T_900_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdpT_gen_tuneD6T_900_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdpT_gen_tuneD6T_900->GetBinContent(ndx);
    double b = lamBar_dNdpT_gen_tuneD6T_900->GetBinContent(ndx);
    double aPlusb = lam_dNdpT_gen_tuneD6T_900->GetBinContent(ndx);
    double aErr = lamReg_dNdpT_gen_tuneD6T_900->GetBinError(ndx);
    double bErr = lamBar_dNdpT_gen_tuneD6T_900->GetBinError(ndx);
    double aPlusbErr = lam_dNdpT_gen_tuneD6T_900->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
    //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdpT_gen_tuneD6T_900_asym_err[ndx] = fError;
  }
  lam_dNdpT_gen_tuneD6T_900_asym->Add( lamBar_dNdpT_gen_tuneD6T_900, -1. );
  lam_dNdpT_gen_tuneD6T_900_asym->Divide( lam_dNdpT_gen_tuneD6T_900 );
  lam_dNdpT_gen_tuneD6T_900_asym->SetError(lam_dNdpT_gen_tuneD6T_900_asym_err);
  //lam_dNdpT_gen_tuneD6T_900_asym->Add( lamBar_dNdpT_gen_tuneD6T_900, -1. );
  //lam_dNdpT_gen_tuneD6T_900_asym->Divide( lam_dNdpT_gen_tuneD6T_900 );
  //lam_dNdpT_gen_tuneD6T_900_asym->SetLineColor(42);

  // P0:
  TFile* fin_pT_p0_gen_lall 
    = TFile::Open("lam_tuneP0_7000/plots/lam_dNdpT_gen_tuneP0_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_p0_gen_lall->Get("c1");
  TH1F* lam_dNdpT_gen_tuneP0_7000 = (TH1F*)c1->GetPrimitive("lam_dNdpT_gen_scale");
  fin_pT_p0_gen_lall->Close();
  delete c1;

  TFile* fin_pT_p0_gen_part
    = TFile::Open("lam_tuneP0_7000/plots/lamReg_dNdpT_gen_tuneP0_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_p0_gen_part->Get("c1");
  TH1F* lamReg_dNdpT_gen_tuneP0_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_gen_scale");
  fin_pT_p0_gen_part->Close();
  delete c1;

  TFile* fin_pT_p0_gen_anti
    = TFile::Open("lam_tuneP0_7000/plots/lamBar_dNdpT_gen_tuneP0_7TeV_weight0_subgeant2.root", "READ");
  c1 = (TCanvas*)fin_pT_p0_gen_anti->Get("c1");
  TH1F* lamBar_dNdpT_gen_tuneP0_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_gen_scale");
  fin_pT_p0_gen_anti->Close();
  delete c1;

  TH1F* lam_dNdpT_gen_tuneP0_7000_asym = (TH1F*)lamReg_dNdpT_gen_tuneP0_7000->Clone("lam_dNdpT_gen_tuneP0_7000_asym");
  double *lam_dNdpT_gen_tuneP0_7000_asym_err 
    = new double[lam_dNdpT_gen_tuneP0_7000_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdpT_gen_tuneP0_7000_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdpT_gen_tuneP0_7000->GetBinContent(ndx);
    double b = lamBar_dNdpT_gen_tuneP0_7000->GetBinContent(ndx);
    double aPlusb = lam_dNdpT_gen_tuneP0_7000->GetBinContent(ndx);
    double aErr = lamReg_dNdpT_gen_tuneP0_7000->GetBinError(ndx);
    double bErr = lamBar_dNdpT_gen_tuneP0_7000->GetBinError(ndx);
    double aPlusbErr = lam_dNdpT_gen_tuneP0_7000->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdpT_gen_tuneP0_7000_asym_err[ndx] = fError;
  }
  lam_dNdpT_gen_tuneP0_7000_asym->Add( lamBar_dNdpT_gen_tuneP0_7000, -1. );
  lam_dNdpT_gen_tuneP0_7000_asym->Divide( lam_dNdpT_gen_tuneP0_7000 );
  lam_dNdpT_gen_tuneP0_7000_asym->SetError(lam_dNdpT_gen_tuneP0_7000_asym_err);
  //lam_dNdpT_gen_tuneP0_7000_asym->Add( lamBar_dNdpT_gen_tuneP0_7000, -1. );
  //lam_dNdpT_gen_tuneP0_7000_asym->Divide( lam_dNdpT_gen_tuneP0_7000 );
  //lam_dNdpT_gen_tuneP0_7000_asym->SetLineColor(31);

  // Now get the UNCORRECTED data lam plots for asymmetries

  TFile* fin_dNdy_uncor_7000
    = TFile::Open("lam_data_7000/plots/lam_dNdy.root", "READ");
  c1 = (TCanvas*)fin_dNdy_uncor_7000->Get("c1");
  TH1F* lam_dNdy_uncor_data_7000 = (TH1F*)c1->GetPrimitive("lam_dNdy");
  fin_dNdy_uncor_7000->Close();
  delete c1;

  TFile* fin_dNdy_uncorReg_7000
    = TFile::Open("lam_data_7000/plots/lamReg_dNdy.root", "READ");
  c1 = (TCanvas*)fin_dNdy_uncorReg_7000->Get("c1");
  TH1F* lamReg_dNdy_uncor_data_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdy");
  fin_dNdy_uncorReg_7000->Close();
  delete c1;

  TFile* fin_dNdy_uncorBar_7000
    = TFile::Open("lam_data_7000/plots/lamBar_dNdy.root", "READ");
  c1 = (TCanvas*)fin_dNdy_uncorBar_7000->Get("c1");
  TH1F* lamBar_dNdy_uncor_data_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdy");
  fin_dNdy_uncorBar_7000->Close();
  delete c1;

  TH1F* lam_dNdy_uncor_data_7000_asym = (TH1F*)lamReg_dNdy_uncor_data_7000->Clone("lam_dNdy_uncor_data_7000_asym");
  double *lam_dNdy_uncor_data_7000_asym_err 
    = new double[lam_dNdy_uncor_data_7000_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdy_uncor_data_7000_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdy_uncor_data_7000->GetBinContent(ndx);
    double b = lamBar_dNdy_uncor_data_7000->GetBinContent(ndx);
    double aPlusb = lam_dNdy_uncor_data_7000->GetBinContent(ndx);
    double aErr = lamReg_dNdy_uncor_data_7000->GetBinError(ndx);
    double bErr = lamBar_dNdy_uncor_data_7000->GetBinError(ndx);
    double aPlusbErr = lam_dNdy_uncor_data_7000->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdy_uncor_data_7000_asym_err[ndx] = fError;
  }
  lam_dNdy_uncor_data_7000_asym->Add( lamBar_dNdy_uncor_data_7000, -1. );
  lam_dNdy_uncor_data_7000_asym->Divide( lam_dNdy_uncor_data_7000 );
  lam_dNdy_uncor_data_7000_asym->SetError(lam_dNdy_uncor_data_7000_asym_err);
  //lam_dNdy_uncor_data_7000_asym->Add(lamBar_dNdy_uncor_data_7000, -1.);
  //lam_dNdy_uncor_data_7000_asym->Divide(lam_dNdy_uncor_data_7000);

  // pT:
  TFile* fin_dNdpT_uncor_7000
    = TFile::Open("lam_data_7000/plots/lam_dNdpT.root", "READ");
  c1 = (TCanvas*)fin_dNdpT_uncor_7000->Get("c1");
  TH1F* lam_dNdpT_uncor_data_7000 = (TH1F*)c1->GetPrimitive("lam_dNdpT_scale");
  fin_dNdpT_uncor_7000->Close();
  delete c1;

  TFile* fin_dNdpT_uncorReg_7000
    = TFile::Open("lam_data_7000/plots/lamReg_dNdpT.root", "READ");
  c1 = (TCanvas*)fin_dNdpT_uncorReg_7000->Get("c1");
  TH1F* lamReg_dNdpT_uncor_data_7000 = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_scale");
  fin_dNdpT_uncorReg_7000->Close();
  delete c1;

  TFile* fin_dNdpT_uncorBar_7000
    = TFile::Open("lam_data_7000/plots/lamBar_dNdpT.root", "READ");
  c1 = (TCanvas*)fin_dNdpT_uncorBar_7000->Get("c1");
  TH1F* lamBar_dNdpT_uncor_data_7000 = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_scale");
  fin_dNdpT_uncorBar_7000->Close();
  delete c1;

  TH1F* lam_dNdpT_uncor_data_7000_asym = (TH1F*)lamReg_dNdpT_uncor_data_7000->Clone("lam_dNdpT_uncor_data_7000_asym");
  double *lam_dNdpT_uncor_data_7000_asym_err 
    = new double[lam_dNdpT_uncor_data_7000_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdpT_uncor_data_7000_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdpT_uncor_data_7000->GetBinContent(ndx);
    double b = lamBar_dNdpT_uncor_data_7000->GetBinContent(ndx);
    double aPlusb = lam_dNdpT_uncor_data_7000->GetBinContent(ndx);
    double aErr = lamReg_dNdpT_uncor_data_7000->GetBinError(ndx);
    double bErr = lamBar_dNdpT_uncor_data_7000->GetBinError(ndx);
    double aPlusbErr = lam_dNdpT_uncor_data_7000->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdpT_uncor_data_7000_asym_err[ndx] = fError;
  }
  lam_dNdpT_uncor_data_7000_asym->Add( lamBar_dNdpT_uncor_data_7000, -1. );
  lam_dNdpT_uncor_data_7000_asym->Divide( lam_dNdpT_uncor_data_7000 );
  lam_dNdpT_uncor_data_7000_asym->SetError(lam_dNdpT_uncor_data_7000_asym_err);
  //lam_dNdpT_uncor_data_7000_asym->Add(lamBar_dNdpT_uncor_data_7000, -1.);
  //lam_dNdpT_uncor_data_7000_asym->Divide(lam_dNdpT_uncor_data_7000);


  // Both for 900 GeV

  TFile* fin_dNdy_uncor_900
    = TFile::Open("lam_data_900/plots/lam_dNdy.root", "READ");
  c1 = (TCanvas*)fin_dNdy_uncor_900->Get("c1");
  TH1F* lam_dNdy_uncor_data_900 = (TH1F*)c1->GetPrimitive("lam_dNdy");
  fin_dNdy_uncor_900->Close();
  delete c1;

  TFile* fin_dNdy_uncorReg_900
    = TFile::Open("lam_data_900/plots/lamReg_dNdy.root", "READ");
  c1 = (TCanvas*)fin_dNdy_uncorReg_900->Get("c1");
  TH1F* lamReg_dNdy_uncor_data_900 = (TH1F*)c1->GetPrimitive("lamReg_dNdy");
  fin_dNdy_uncorReg_900->Close();
  delete c1;

  TFile* fin_dNdy_uncorBar_900
    = TFile::Open("lam_data_900/plots/lamBar_dNdy.root", "READ");
  c1 = (TCanvas*)fin_dNdy_uncorBar_900->Get("c1");
  TH1F* lamBar_dNdy_uncor_data_900 = (TH1F*)c1->GetPrimitive("lamBar_dNdy");
  fin_dNdy_uncorBar_900->Close();
  delete c1;

  TH1F* lam_dNdy_uncor_data_900_asym = (TH1F*)lamReg_dNdy_uncor_data_900->Clone("lam_dNdy_uncor_data_900_asym");
  double *lam_dNdy_uncor_data_900_asym_err 
    = new double[lam_dNdy_uncor_data_900_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdy_uncor_data_900_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdy_uncor_data_900->GetBinContent(ndx);
    double b = lamBar_dNdy_uncor_data_900->GetBinContent(ndx);
    double aPlusb = lam_dNdy_uncor_data_900->GetBinContent(ndx);
    double aErr = lamReg_dNdy_uncor_data_900->GetBinError(ndx);
    double bErr = lamBar_dNdy_uncor_data_900->GetBinError(ndx);
    double aPlusbErr = lam_dNdy_uncor_data_900->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdy_uncor_data_900_asym_err[ndx] = fError;
  }
  lam_dNdy_uncor_data_900_asym->Add( lamBar_dNdy_uncor_data_900, -1. );
  lam_dNdy_uncor_data_900_asym->Divide( lam_dNdy_uncor_data_900 );
  lam_dNdy_uncor_data_900_asym->SetError(lam_dNdy_uncor_data_900_asym_err);
  //lam_dNdy_uncor_data_900_asym->Add(lamBar_dNdy_uncor_data_900, -1.);
  //lam_dNdy_uncor_data_900_asym->Divide(lam_dNdy_uncor_data_900);

  // pT:
  TFile* fin_dNdpT_uncor_900
    = TFile::Open("lam_data_900/plots/lam_dNdpT.root", "READ");
  c1 = (TCanvas*)fin_dNdpT_uncor_900->Get("c1");
  TH1F* lam_dNdpT_uncor_data_900 = (TH1F*)c1->GetPrimitive("lam_dNdpT_scale");
  fin_dNdpT_uncor_900->Close();
  delete c1;

  TFile* fin_dNdpT_uncorReg_900
    = TFile::Open("lam_data_900/plots/lamReg_dNdpT.root", "READ");
  c1 = (TCanvas*)fin_dNdpT_uncorReg_900->Get("c1");
  TH1F* lamReg_dNdpT_uncor_data_900 = (TH1F*)c1->GetPrimitive("lamReg_dNdpT_scale");
  fin_dNdpT_uncorReg_900->Close();
  delete c1;

  TFile* fin_dNdpT_uncorBar_900
    = TFile::Open("lam_data_900/plots/lamBar_dNdpT.root", "READ");
  c1 = (TCanvas*)fin_dNdpT_uncorBar_900->Get("c1");
  TH1F* lamBar_dNdpT_uncor_data_900 = (TH1F*)c1->GetPrimitive("lamBar_dNdpT_scale");
  fin_dNdpT_uncorBar_900->Close();
  delete c1;

  TH1F* lam_dNdpT_uncor_data_900_asym = (TH1F*)lamReg_dNdpT_uncor_data_900->Clone("lam_dNdpT_uncor_data_900_asym");
  double *lam_dNdpT_uncor_data_900_asym_err 
    = new double[lam_dNdpT_uncor_data_900_asym->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdpT_uncor_data_900_asym->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdpT_uncor_data_900->GetBinContent(ndx);
    double b = lamBar_dNdpT_uncor_data_900->GetBinContent(ndx);
    double aPlusb = lam_dNdpT_uncor_data_900->GetBinContent(ndx);
    double aErr = lamReg_dNdpT_uncor_data_900->GetBinError(ndx);
    double bErr = lamBar_dNdpT_uncor_data_900->GetBinError(ndx);
    double aPlusbErr = lam_dNdpT_uncor_data_900->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);// + 
				      //(a - b)*(a - b)/aPlusb/aPlusb*aPlusbErr*aPlusbErr);
    lam_dNdpT_uncor_data_900_asym_err[ndx] = fError;
  }
  lam_dNdpT_uncor_data_900_asym->Add( lamBar_dNdpT_uncor_data_900, -1. );
  lam_dNdpT_uncor_data_900_asym->Divide( lam_dNdpT_uncor_data_900 );
  lam_dNdpT_uncor_data_900_asym->SetError(lam_dNdpT_uncor_data_900_asym_err);
  //lam_dNdpT_uncor_data_900_asym->Add(lamBar_dNdpT_uncor_data_900, -1.);
  //lam_dNdpT_uncor_data_900_asym->Divide(lam_dNdpT_uncor_data_900);


  // That's a lot of crap, so here's the list of what we just made:
  /*
    lam_dNdy_gen_pythia8_7000_asym
    lam_dNdy_gen_tuneD6T_7000_asym
    lam_dNdy_gen_tuneD6T_900_asym
    lam_dNdy_gen_tuneP0_7000_asym

    lam_dNdpT_gen_pythia8_7000_asym
    lam_dNdpT_gen_tuneD6T_7000_asym
    lam_dNdpT_gen_tuneD6T_900_asym
    lam_dNdpT_gen_tuneP0_7000_asym

    lam_dNdy_uncor_data_7000_asym
    lam_dNdpT_uncor_data_7000_asym
    lam_dNdy_uncor_data_900_asym
    lam_dNdpT_uncor_data_900_asym

  */
  lam_dNdy_gen_pythia8_7000_asym->SetLineColor(4);
  lam_dNdy_gen_pythia8_7000_asym->SetMarkerStyle(1);
  lam_dNdy_gen_pythia8_7000_asym->SetLineStyle(2);

  lam_dNdy_gen_tuneD6T_7000_asym->SetLineColor(4);
  lam_dNdy_gen_tuneD6T_7000_asym->SetMarkerStyle(1);
  lam_dNdy_gen_tuneD6T_7000_asym->SetLineStyle(1);

  lam_dNdy_gen_tuneD6T_900_asym->SetLineColor(2);
  lam_dNdy_gen_tuneD6T_900_asym->SetMarkerStyle(1);
  lam_dNdy_gen_tuneD6T_900_asym->SetLineStyle(1);

  lam_dNdy_gen_tuneP0_7000_asym->SetLineColor(4);
  lam_dNdy_gen_tuneP0_7000_asym->SetMarkerStyle(1);
  lam_dNdy_gen_tuneP0_7000_asym->SetLineStyle(3);

  lam_dNdpT_gen_pythia8_7000_asym->SetLineColor(4);
  lam_dNdpT_gen_pythia8_7000_asym->SetMarkerStyle(1);
  lam_dNdpT_gen_pythia8_7000_asym->SetLineStyle(2);

  lam_dNdpT_gen_tuneD6T_7000_asym->SetLineColor(4);
  lam_dNdpT_gen_tuneD6T_7000_asym->SetMarkerStyle(1);
  lam_dNdpT_gen_tuneD6T_7000_asym->SetLineStyle(1);

  lam_dNdpT_gen_tuneD6T_900_asym->SetLineColor(2);
  lam_dNdpT_gen_tuneD6T_900_asym->SetMarkerStyle(1);
  lam_dNdpT_gen_tuneD6T_900_asym->SetLineStyle(1);

  lam_dNdpT_gen_tuneP0_7000_asym->SetLineColor(4);
  lam_dNdpT_gen_tuneP0_7000_asym->SetMarkerStyle(1);
  lam_dNdpT_gen_tuneP0_7000_asym->SetLineStyle(3);


  lam_dNdy_uncor_data_7000_asym->SetLineColor(7);
  lam_dNdy_uncor_data_7000_asym->SetMarkerColor(7);
  lam_dNdy_uncor_data_7000_asym->SetMarkerStyle(23);
  lam_dNdy_uncor_data_7000_asym->SetMarkerSize(0.8);

  lam_dNdpT_uncor_data_7000_asym->SetLineColor(7);
  lam_dNdpT_uncor_data_7000_asym->SetMarkerColor(7);
  lam_dNdpT_uncor_data_7000_asym->SetMarkerStyle(23);
  lam_dNdpT_uncor_data_7000_asym->SetMarkerSize(0.8);

  lam_dNdy_uncor_data_900_asym->SetLineColor(6);
  lam_dNdy_uncor_data_900_asym->SetMarkerColor(6);
  lam_dNdy_uncor_data_900_asym->SetMarkerStyle(22);
  lam_dNdy_uncor_data_900_asym->SetMarkerSize(0.9);

  lam_dNdpT_uncor_data_900_asym->SetLineColor(6);
  lam_dNdpT_uncor_data_900_asym->SetMarkerColor(6);
  lam_dNdpT_uncor_data_900_asym->SetMarkerStyle(22);
  lam_dNdpT_uncor_data_900_asym->SetMarkerSize(0.9);

  TF1* zeroLiney = new TF1("zeroLiney", "[0]", 0., 2.3);
  TF1* zeroLinepT = new TF1("zeroLinepT", "[0]", 0., 10.);
  zeroLiney->SetParameter(0, 0.);
  zeroLinepT->SetParameter(0, 0.);
  zeroLiney->SetLineColor(1);
  zeroLinepT->SetLineColor(1);


  // Do lambda/lambdaBar asymmetry plots
  TH1F* lamAsym_pt_cor_data_7000_sys_num = (TH1F*)lamReg_dNdpT_cor_data_7000_sys->Clone("lamAsym_pt_cor_data_7000_sys_num");
  TH1F* lamAsym_pt_cor_data_7000_sys_den = (TH1F*)lam_dNdpT_cor_data_7000_sys->Clone("lamAsym_pt_cor_data_7000_sys_den");
  double *lam_pt_cor_data_7000_sys_err = new double[lam_dNdpT_cor_data_7000_sys->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdpT_cor_data_7000_sys->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdpT_cor_data_7000_sys->GetBinContent(ndx);
    double b = lamBar_dNdpT_cor_data_7000_sys->GetBinContent(ndx);
    double aPlusb = lam_dNdpT_cor_data_7000_sys->GetBinContent(ndx);
    //cout << "ndx: " << ndx << ", aPlusb: " << aPlusb << endl;
    if( aPlusb <= 0.) continue;
    double aErr = lamReg_dNdpT_cor_data_7000_sys->GetBinError(ndx);
    double bErr = lamBar_dNdpT_cor_data_7000_sys->GetBinError(ndx);
    double aPlusbErr = lam_dNdpT_cor_data_7000_sys->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);
    lam_pt_cor_data_7000_sys_err[ndx] = fError;
  }
  lamAsym_pt_cor_data_7000_sys_num->Add(lamBar_dNdpT_cor_data_7000_sys, -1.);
  lamAsym_pt_cor_data_7000_sys_num->Divide(lamAsym_pt_cor_data_7000_sys_den);
  lamAsym_pt_cor_data_7000_sys_num->SetError(lam_pt_cor_data_7000_sys_err);
  //lamAsym_pt_cor_data_7000_sys_num->Add(lamBar_dNdpT_cor_data_7000_sys, -1.);
  //lamAsym_pt_cor_data_7000_sys_num->Divide(lamAsym_pt_cor_data_7000_sys_den);

  lamAsym_pt_cor_data_7000_sys_num->SetMarkerStyle(20);
  lamAsym_pt_cor_data_7000_sys_num->SetMarkerColor(4);
  lamAsym_pt_cor_data_7000_sys_num->SetLineColor(4);
  lamAsym_pt_cor_data_7000_sys_num->SetYTitle("#sigma(#Lambda^{0}) - #sigma(#bar{#Lambda}^{0}) / #sigma(#Lambda^{0}+#bar{#Lambda}^{0})");
  lamAsym_pt_cor_data_7000_sys_num->SetXTitle("p_{T} [GeV/c]");
  lamAsym_pt_cor_data_7000_sys_num->SetMinimum(-0.4);
  lamAsym_pt_cor_data_7000_sys_num->SetMaximum(0.4);

  TH1F* lamAsym_pt_cor_data_900_sys_num = (TH1F*)lamReg_dNdpT_cor_data_900_sys->Clone("lamAsym_pt_cor_data_900_sys_num");
  TH1F* lamAsym_pt_cor_data_900_sys_den = (TH1F*)lam_dNdpT_cor_data_900_sys->Clone("lamAsym_pt_cor_data_900_sys_den");
  double *lam_pt_cor_data_900_sys_err = new double[lam_dNdpT_cor_data_900_sys->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdpT_cor_data_900_sys->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdpT_cor_data_900_sys->GetBinContent(ndx);
    double b = lamBar_dNdpT_cor_data_900_sys->GetBinContent(ndx);
    double aPlusb = lam_dNdpT_cor_data_900_sys->GetBinContent(ndx);
    if( aPlusb <= 0. ) continue;
    double aErr = lamReg_dNdpT_cor_data_900_sys->GetBinError(ndx);
    double bErr = lamBar_dNdpT_cor_data_900_sys->GetBinError(ndx);
    double aPlusbErr = lam_dNdpT_cor_data_900_sys->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);
    lam_pt_cor_data_900_sys_err[ndx] = fError;
  }
  lamAsym_pt_cor_data_900_sys_num->Add(lamBar_dNdpT_cor_data_900_sys, -1.);
  lamAsym_pt_cor_data_900_sys_num->Divide(lamAsym_pt_cor_data_900_sys_den);
  lamAsym_pt_cor_data_900_sys_num->SetError(lam_pt_cor_data_900_sys_err);
  //lamAsym_pt_cor_data_900_sys_num->Add(lamBar_dNdpT_cor_data_900_sys, -1.);
  //lamAsym_pt_cor_data_900_sys_num->Divide(lamAsym_pt_cor_data_900_sys_den);

  lamAsym_pt_cor_data_900_sys_num->SetMarkerStyle(21);
  lamAsym_pt_cor_data_900_sys_num->SetMarkerColor(2);
  lamAsym_pt_cor_data_900_sys_num->SetLineColor(2);

  TLegend *leg_ptAsym = new TLegend(0.2,0.2,0.55,0.38);
  leg_ptAsym->SetTextSize(0.042);
  leg_ptAsym->SetBorderSize(0);
  leg_ptAsym->SetFillColor(0);
  leg_ptAsym->SetFillStyle(0);
  leg_ptAsym->AddEntry(lamAsym_pt_cor_data_7000_sys_num, "  #sqrt{s} = 7 TeV" , "pl" );
  leg_ptAsym->AddEntry(lam_dNdpT_gen_tuneD6T_7000_asym, "     D6T MC", "l");
  leg_ptAsym->AddEntry(lam_dNdpT_gen_pythia8_7000_asym, "     Pythia8 MC", "l");
  leg_ptAsym->AddEntry(lam_dNdpT_gen_tuneP0_7000_asym, "     P0 MC", "l");
  TLegend *leg_ptAsym900 = new TLegend(0.63, 0.2, 0.97, 0.3);
  leg_ptAsym900->SetTextSize(0.042);
  leg_ptAsym900->SetBorderSize(0);
  leg_ptAsym900->SetFillColor(0);
  leg_ptAsym900->SetFillStyle(0);
  leg_ptAsym900->AddEntry(lamAsym_pt_cor_data_900_sys_num, "  #sqrt{s} = 900 GeV" , "pl" );
  leg_ptAsym900->AddEntry(lam_dNdpT_gen_tuneD6T_900_asym, "     D6T MC", "l");

  c1 = new TCanvas("c1", "asym");
  lamAsym_pt_cor_data_7000_sys_num->Draw("E1");
  zeroLinepT->Draw("same");
  lam_dNdpT_gen_tuneD6T_7000_asym->Draw("histsame");
  lam_dNdpT_gen_pythia8_7000_asym->Draw("histsame");
  lam_dNdpT_gen_tuneP0_7000_asym->Draw("histsame");
  lam_dNdpT_gen_tuneD6T_900_asym->Draw("histsame");
  lamAsym_pt_cor_data_7000_sys_num->Draw("E1same");
  lamAsym_pt_cor_data_900_sys_num->Draw("E1same");
  TString ptAsymName = "final_plots/lamAsym_pt_sys";
  leg_ptAsym->Draw("same");
  leg_ptAsym900->Draw("same");
  cmsprelim->DrawLatex(0.78, 0.90, "CMS Preliminary");

  c1->SaveAs(ptAsymName + ".root");
  c1->SaveAs(ptAsymName + ".eps");
  c1->SaveAs(ptAsymName + ".C");
  c1->SaveAs(ptAsymName + ".png");
  delete c1;



  TH1F* lamAsym_y_cor_data_7000_sys_num = (TH1F*)lamReg_dNdy_cor_data_7000_sys->Clone("lamAsym_y_cor_data_7000_sys_num");
  TH1F* lamAsym_y_cor_data_7000_sys_den = (TH1F*)lam_dNdy_cor_data_7000_sys->Clone("lamAsym_y_cor_data_7000_sys_den");
  double *lam_y_cor_data_7000_sys_err = new double[lam_dNdy_cor_data_7000_sys->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdy_cor_data_7000_sys->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdy_cor_data_7000_sys->GetBinContent(ndx);
    double b = lamBar_dNdy_cor_data_7000_sys->GetBinContent(ndx);
    double aPlusb = lam_dNdy_cor_data_7000_sys->GetBinContent(ndx);
    double aErr = lamReg_dNdy_cor_data_7000_sys->GetBinError(ndx);
    double bErr = lamBar_dNdy_cor_data_7000_sys->GetBinError(ndx);
    double aPlusbErr = lam_dNdy_cor_data_7000_sys->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);
    lam_y_cor_data_7000_sys_err[ndx] = fError;
  }
  lamAsym_y_cor_data_7000_sys_num->Add(lamBar_dNdy_cor_data_7000_sys, -1.);
  lamAsym_y_cor_data_7000_sys_num->Divide(lamAsym_y_cor_data_7000_sys_den);
  lamAsym_y_cor_data_7000_sys_num->SetError(lam_y_cor_data_7000_sys_err);

  lamAsym_y_cor_data_7000_sys_num->SetMarkerStyle(20);
  lamAsym_y_cor_data_7000_sys_num->SetMarkerColor(4);
  lamAsym_y_cor_data_7000_sys_num->SetLineColor(4);
  lamAsym_y_cor_data_7000_sys_num->SetYTitle("#sigma(#Lambda^{0}) - #sigma(#bar{#Lambda}^{0}) / #sigma(#Lambda^{0} + #bar{#Lambda}^{0})");
  lamAsym_y_cor_data_7000_sys_num->SetXTitle("rapidity");
  lamAsym_y_cor_data_7000_sys_num->SetMinimum(-0.2);
  lamAsym_y_cor_data_7000_sys_num->SetMaximum(0.2);

  TH1F* lamAsym_y_cor_data_900_sys_num = (TH1F*)lamReg_dNdy_cor_data_900_sys->Clone("lamAsym_y_cor_data_900_sys_num");
  TH1F* lamAsym_y_cor_data_900_sys_den = (TH1F*)lam_dNdy_cor_data_900_sys->Clone("lamAsym_y_cor_data_900_sys_den");
  double *lam_y_cor_data_900_sys_err = new double[lam_dNdy_cor_data_900_sys->GetNbinsX() + 1];
  for( int ndx = 1; ndx <= lam_dNdy_cor_data_900_sys->GetNbinsX(); ndx++ ) {
    double a = lamReg_dNdy_cor_data_900_sys->GetBinContent(ndx);
    double b = lamBar_dNdy_cor_data_900_sys->GetBinContent(ndx);
    double aPlusb = lam_dNdy_cor_data_900_sys->GetBinContent(ndx);
    double aErr = lamReg_dNdy_cor_data_900_sys->GetBinError(ndx);
    double bErr = lamBar_dNdy_cor_data_900_sys->GetBinError(ndx);
    double aPlusbErr = lam_dNdy_cor_data_900_sys->GetBinError(ndx);
    double fError = 2 / aPlusb / aPlusb * sqrt(b*b*aErr*aErr + a*a*bErr*bErr);
    lam_y_cor_data_900_sys_err[ndx] = fError;
  }
  lamAsym_y_cor_data_900_sys_num->Add(lamBar_dNdy_cor_data_900_sys, -1.);
  lamAsym_y_cor_data_900_sys_num->Divide(lamAsym_y_cor_data_900_sys_den);
  lamAsym_y_cor_data_900_sys_num->SetError(lam_y_cor_data_900_sys_err);
  //lamAsym_y_cor_data_900_sys_num->Add(lamBar_dNdy_cor_data_900_sys, -1.);
  //lamAsym_y_cor_data_900_sys_num->Divide(lamAsym_y_cor_data_900_sys_den);

  lamAsym_y_cor_data_900_sys_num->SetMarkerStyle(21);
  lamAsym_y_cor_data_900_sys_num->SetMarkerColor(2);
  lamAsym_y_cor_data_900_sys_num->SetLineColor(2);

  TLegend *leg_yAsym = new TLegend(0.2,0.2,0.55,0.38);
  leg_yAsym->SetTextSize(0.042);
  leg_yAsym->SetBorderSize(0);
  leg_yAsym->SetFillColor(0);
  leg_yAsym->SetFillStyle(0);
  leg_yAsym->AddEntry(lamAsym_y_cor_data_7000_sys_num, "  #sqrt{s} = 7 TeV" , "pl" );
  leg_yAsym->AddEntry(lam_dNdy_gen_tuneD6T_7000_asym, "     D6T MC", "l");
  leg_yAsym->AddEntry(lam_dNdy_gen_pythia8_7000_asym, "     Pythia8 MC", "l");
  leg_yAsym->AddEntry(lam_dNdy_gen_tuneP0_7000_asym, "     P0 MC", "l");
  TLegend *leg_yAsym900 = new TLegend(0.63, 0.2, 0.97, 0.3);
  leg_yAsym900->SetTextSize(0.042);
  leg_yAsym900->SetBorderSize(0);
  leg_yAsym900->SetFillColor(0);
  leg_yAsym900->SetFillStyle(0);
  leg_yAsym900->AddEntry(lamAsym_y_cor_data_900_sys_num, "  #sqrt{s} = 900 GeV" , "pl" );
  leg_yAsym900->AddEntry(lam_dNdy_gen_tuneD6T_900_asym, "     D6T MC", "l");

  c1 = new TCanvas("c1", "asym");
  lamAsym_y_cor_data_7000_sys_num->Draw("E1");
  zeroLiney->Draw("same");
  lam_dNdy_gen_tuneD6T_7000_asym->Draw("histsame");
  lam_dNdy_gen_pythia8_7000_asym->Draw("histsame");
  lam_dNdy_gen_tuneP0_7000_asym->Draw("histsame");
  lam_dNdy_gen_tuneD6T_900_asym->Draw("histsame");
  lamAsym_y_cor_data_7000_sys_num->Draw("E1same");
  lamAsym_y_cor_data_900_sys_num->Draw("E1same");
  leg_yAsym->Draw("same");
  leg_yAsym900->Draw("same");
  cmsprelim->DrawLatex(0.78, 0.90, "CMS Preliminary");
  TString yAsymName = "final_plots/lamAsym_y_sys";

  c1->SaveAs(yAsymName + ".root");
  c1->SaveAs(yAsymName + ".eps");
  c1->SaveAs(yAsymName + ".C");
  c1->SaveAs(yAsymName + ".png");
  delete c1;

  c1 = new TCanvas("c1", "yAsym with uncorrected");

  TLegend *leg_yAsymUnc = new TLegend(0.2,0.19,0.55,0.37);
  leg_yAsymUnc->SetTextSize(0.04);
  leg_yAsymUnc->SetBorderSize(0);
  leg_yAsymUnc->SetFillColor(0);
  leg_yAsymUnc->SetFillStyle(0);
  leg_yAsymUnc->AddEntry(lamAsym_y_cor_data_7000_sys_num, "  #sqrt{s} = 7 TeV corrected", "pl");
  leg_yAsymUnc->AddEntry(lam_dNdy_uncor_data_7000_asym, "  #sqrt{s} = 7 TeV uncorrected", "pl");
  leg_yAsymUnc->AddEntry(lamAsym_y_cor_data_900_sys_num, "  #sqrt{s} = 900 GeV corrected", "pl");
  leg_yAsymUnc->AddEntry(lam_dNdy_uncor_data_900_asym, "  #sqrt{s} = 900 GeV uncorrected", "pl");


  //lamAsym_y_cor_data_7000_sys_num->SetMaximum(0.2);
  //lamAsym_y_cor_data_7000_sys_num->SetMinimum(-0.2);

  lamAsym_y_cor_data_7000_sys_num->Draw("E1");
  zeroLiney->Draw("same");
  lam_dNdy_uncor_data_7000_asym->Draw("E1same");
  lam_dNdy_uncor_data_900_asym->Draw("E1same");
  lamAsym_y_cor_data_900_sys_num->Draw("E1same");
  lamAsym_y_cor_data_7000_sys_num->Draw("E1same");
  leg_yAsymUnc->Draw("same");
  cmsprelim->DrawLatex(0.78, 0.90, "CMS Preliminary");

  TString yAsymNameUncor = "final_plots/lamAsym_y_withUncor";
  c1->SaveAs(yAsymNameUncor + ".root");
  c1->SaveAs(yAsymNameUncor + ".eps");
  c1->SaveAs(yAsymNameUncor + ".C");
  c1->SaveAs(yAsymNameUncor + ".png");
  delete c1;


  c1 = new TCanvas("c1", "ptAsym with uncorrected");

  TLegend *leg_ptAsymUnc = new TLegend(0.2,0.18,0.55,0.38);
  leg_ptAsymUnc->SetTextSize(0.04);
  leg_ptAsymUnc->SetBorderSize(0);
  leg_ptAsymUnc->SetFillColor(0);
  leg_ptAsymUnc->SetFillStyle(0);
  leg_ptAsymUnc->AddEntry(lamAsym_pt_cor_data_7000_sys_num, "  #sqrt{s} = 7 TeV corrected", "pl");
  leg_ptAsymUnc->AddEntry(lam_dNdpT_uncor_data_7000_asym, "  #sqrt{s} = 7 TeV uncorrected", "pl");
  leg_ptAsymUnc->AddEntry(lamAsym_pt_cor_data_900_sys_num, "  #sqrt{s} = 900 GeV corrected", "pl");
  leg_ptAsymUnc->AddEntry(lam_dNdpT_uncor_data_900_asym, "  #sqrt{s} = 900 GeV uncorrected", "pl");


  //lamAsym_pt_cor_data_7000_sys_num->SetMaximum(0.4);
  //lamAsym_pt_cor_data_7000_sys_num->SetMinimum(-0.4);

  lamAsym_pt_cor_data_7000_sys_num->Draw("E1");
  zeroLinepT->Draw("same");
  lam_dNdpT_uncor_data_7000_asym->Draw("E1same");
  lam_dNdpT_uncor_data_900_asym->Draw("E1same");
  lamAsym_pt_cor_data_900_sys_num->Draw("E1same");
  lamAsym_pt_cor_data_7000_sys_num->Draw("E1same");
  leg_ptAsymUnc->Draw("same");
  cmsprelim->DrawLatex(0.78, 0.90, "CMS Preliminary");

  TString ptAsymNameUncor = "final_plots/lamAsym_pt_withUncor";
  c1->SaveAs(ptAsymNameUncor + ".root");
  c1->SaveAs(ptAsymNameUncor + ".eps");
  c1->SaveAs(ptAsymNameUncor + ".C");
  c1->SaveAs(ptAsymNameUncor + ".png");
  delete c1;

  //---------------------

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
  cmsprelim->DrawLatex(0.78,0.90, "CMS Preliminary");

  TString ratioOutpT = "final_plots/lamOverKs_pt_cor_data_sys";
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
  cmsprelim->DrawLatex(0.78,0.90,"CMS Preliminary");

  TString ratioOuty = "final_plots/lamOverKs_y_cor_data_sys";
  c1->SaveAs( ratioOuty + ".root" );
  c1->SaveAs( ratioOuty + ".eps" );
  c1->SaveAs( ratioOuty + ".C" );
  c1->SaveAs( ratioOuty + ".png" );

  delete c1;

  c1 = new TCanvas("c1", "ks y");
  ks_dNdy_posOverNeg_7000 = (TH1F*)ks_dNdy_pos_cor_data_7000->Clone("ks_dNdy_posOverNeg_7000");
  ks_dNdy_posOverNeg_7000_den = (TH1F*)ks_dNdy_neg_cor_data_7000->Clone("ks_dNdy_posOverNeg_7000_den");
  ks_dNdy_posOverNeg_7000->Divide(ks_dNdy_posOverNeg_7000_den);

  TF1* constOne = new TF1("constOne", "[0]", 0., 2.3);
  constOne->SetParameter(0, 1.);
  //constOne->FixParameter(0, 1.);
  ks_dNdy_posOverNeg_7000->Fit("constOne", "RLE");
  tdrStyle->SetOptFit(1);
  ostringstream fitstrfory;
  fitstrfory.precision(2);
  fitstrfory << "#chi^{2} / NDOF = " << constOne->GetChisquare() << " / " << constOne->GetNDF();
  ostringstream fitstrfory2;
  fitstrfory2.setf(ios::fixed, ios::floatfield);
  fitstrfory2.precision(1);
  fitstrfory2 << "Fit value = " << constOne->GetParameter(0); 
  ostringstream fitstrfory3;
  fitstrfory3.precision(1);
  fitstrfory3 << " #pm " << constOne->GetParError(0);

  ks_dNdy_posOverNeg_900 = (TH1F*)ks_dNdy_pos_cor_data_900->Clone("ks_dNdy_posOverNeg_900");
  ks_dNdy_posOverNeg_900_den = (TH1F*)ks_dNdy_neg_cor_data_900->Clone("ks_dNdy_posOverNeg_900_den");
  ks_dNdy_posOverNeg_900->Divide(ks_dNdy_posOverNeg_900_den);

  ks_dNdy_posOverNeg_7000->SetYTitle("#sigma(K^{0}_{S}, y > 0) / #sigma(K^{0}_{S}, y < 0)");
  ks_dNdy_posOverNeg_900->SetYTitle("#sigma(K^{0}_{S}, y > 0) / #sigma(K^{0}_{S}, y < 0)");
  ks_dNdy_posOverNeg_7000->SetXTitle("rapidity |y|");
  ks_dNdy_posOverNeg_900->SetXTitle("rapidity |y|");

  ks_dNdy_posOverNeg_7000->SetMarkerStyle(20);
  ks_dNdy_posOverNeg_7000->SetMarkerColor(4);
  ks_dNdy_posOverNeg_7000->SetLineColor(4);

  ks_dNdy_posOverNeg_900->SetMarkerStyle(21);
  ks_dNdy_posOverNeg_900->SetMarkerColor(2);
  ks_dNdy_posOverNeg_900->SetLineColor(2);

  ks_dNdy_posOverNeg_7000->SetMinimum(0.75);
  ks_dNdy_posOverNeg_7000->SetMaximum(1.25);

  TLegend* leg_yAsym = new TLegend(0.2, 0.25, 0.5, 0.35);
  leg_yAsym->SetTextSize(0.042);
  leg_yAsym->SetBorderSize(0);
  leg_yAsym->SetFillColor(0);
  leg_yAsym->SetFillStyle(0);
  leg_yAsym->AddEntry(ks_dNdy_posOverNeg_7000, "  #sqrt{s} = 7 TeV" , "pl" );
  leg_yAsym->AddEntry(ks_dNdy_posOverNeg_900, "  #sqrt{s} = 900 GeV" , "pl" );



  ks_dNdy_posOverNeg_7000->Draw("E1");
  //ks_dNdy_posOverNeg_900->Draw("E1same");
  //leg_yAsym->Draw("same");
  cmsprelim->DrawLatex(0.78, 0.90, "CMS Preliminary");
  cmsprelim->DrawLatex(0.7, 0.28, fitstrfory.str().c_str());
  cmsprelim->DrawLatex(0.7, 0.35, (fitstrfory2.str() + fitstrfory3.str()).c_str());

  TString yAsymOut = "final_plots/posYovernegY_ks";
  c1->SaveAs( yAsymOut + ".root" );
  c1->SaveAs( yAsymOut + ".eps" );
  c1->SaveAs( yAsymOut + ".C" );
  c1->SaveAs( yAsymOut + ".png" );

  delete c1;

}
