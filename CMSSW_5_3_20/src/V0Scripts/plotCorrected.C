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
std::cout.precision(3);

double trig7000 = 15540204;
double trig900 = 12325647;
double lum7000 = 291850.0; // inverse millibarns
double lum900 = 321970.0; // inverse millibarns
double ptbinwid = 0.1;
double ybinwid = 0.1;

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

Double_t ptent[43];

TString outFile_y_evt = "final_plots/ksy_cor_data_perevt";
TString outFile_y_lum = "final_plots/ksy_cor_data_perlum";

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
ks_dNdpT_cor_data_7000->GetYaxis()->SetTitle("N(pp #rightarrow K_{S}^{0}X) / Event_{NSD} / (GeV/c)");
ks_dNdpT_cor_data_7000->SetLineColor(4);
ks_dNdpT_cor_data_7000->SetMarkerColor(4);
ks_dNdpT_cor_data_7000->SetMarkerStyle(20);
ks_dNdpT_cor_data_7000->SetMarkerSize(0.8);
TH1F *ks_dNdpT_cor_data_7000_sys = (TH1F*)ks_dNdpT_cor_data_7000->Clone("ks_dNdpT_cor_data_700_sys");
cout << "ks pt 7000 integral / event: " << ks_yld_7000_pt << endl;
for (int i = 1; i <= ks_dNdpT_cor_data_7000_sys->GetXaxis()->GetNbins(); ++i) {
  double value = ks_dNdpT_cor_data_7000_sys->GetBinContent(i);
  double error = ks_dNdpT_cor_data_7000_sys->GetBinError(i);
  double syserr = 0.02*value;
  ptent[i-1] = sqrt(error*error + syserr*syserr);
  cout << ks_dNdpT_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i-1] << endl;
}
ks_dNdpT_cor_data_7000_sys->SetError(ptent);
ks_dNdpT_cor_data_7000->Fit("tsfunc7000","EM0");
tsfunc7000->SetNpx(2000);
tsfunc7000->SetLineColor(4);

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
ks_dNdpT_cor_data_900->GetYaxis()->SetTitle("N(pp #rightarrow K_{S}^{0}X) / Event_{NSD} / (GeV/c)");
ks_dNdpT_cor_data_900->SetLineColor(2);
ks_dNdpT_cor_data_900->SetMarkerColor(2);
ks_dNdpT_cor_data_900->SetMarkerStyle(21);
ks_dNdpT_cor_data_900->SetMarkerSize(0.8);
TH1F *ks_dNdpT_cor_data_900_sys = (TH1F*)ks_dNdpT_cor_data_900->Clone("ks_dNdpT_cor_data_700_sys");
cout << "ks pt 900 integral / event: " << ks_yld_900_pt << endl;
for (int i = 1; i <= ks_dNdpT_cor_data_900_sys->GetXaxis()->GetNbins(); ++i) {
  double value = ks_dNdpT_cor_data_900_sys->GetBinContent(i);
  double error = ks_dNdpT_cor_data_900_sys->GetBinError(i);
  double syserr = 0.02*value;
  ptent[i-1] = sqrt(error*error + syserr*syserr);
  cout << ks_dNdpT_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << ptent[i-1] << endl;
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
ks_dNdpT_cor_data_7000->Draw("E1same");
tsfunc7000->Draw("same"); 
ks_dNdpT_cor_data_900_sys->Draw("same");
ks_dNdpT_cor_data_900->Draw("E1same");
tsfunc900->Draw("same"); 
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

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
cout << "ks pt 7000 xsec: " << ks_xsec_7000_pt << endl;
ks_dNdpT_cor_data_7000_sys->GetYaxis()->SetTitle("#sigma(pp #rightarrow K_{S}^{0}X) [ mb / (GeV/c)]");
ks_dNdpT_cor_data_7000_sys->Scale(trig7000/lum7000);
ks_dNdpT_cor_data_7000->Fit("tsfunc7000","EM0");
tsfunc7000->SetNpx(2000);
tsfunc7000->SetLineColor(4);

ks_dNdpT_cor_data_900->Scale(trig900/lum900);
double ks_xsec_900_pt = calc_integral(ks_dNdpT_cor_data_900, 1.0);
cout << "ks pt 900 xsec: " << ks_xsec_900_pt << endl;
ks_dNdpT_cor_data_900_sys->GetYaxis()->SetTitle("#sigma(pp #rightarrow K_{S}^{0}X) [ mb / (GeV/c)]");
ks_dNdpT_cor_data_900_sys->Scale(trig900/lum900);
cout << "pt 900\n";
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
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

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
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

gPad->SetLogy(kTRUE);
c1->SaveAs(outFile_pt_lum + "_log.root");
c1->SaveAs(outFile_pt_lum + "_log.eps");
c1->SaveAs(outFile_pt_lum + "_log.C");
c1->SaveAs(outFile_pt_lum + "_log.png");
gPad->SetLogy(kFALSE);

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
ks_dNdy_cor_data_7000->GetYaxis()->SetTitle("N(pp #rightarrow K_{S}^{0}X) / Event_{NSD}");
ks_dNdy_cor_data_7000->SetMaximum(0.8);
ks_dNdy_cor_data_7000->SetMinimum(0.);
ks_dNdy_cor_data_7000->SetLineColor(4);
ks_dNdy_cor_data_7000->SetMarkerColor(4);
ks_dNdy_cor_data_7000->SetMarkerStyle(20);
ks_dNdy_cor_data_7000->SetMarkerSize(0.8);
TH1F *ks_dNdy_cor_data_7000_sys = (TH1F*)ks_dNdy_cor_data_7000->Clone("ks_dNdy_cor_data_700_sys");
cout << "ks y 7000 integral / event: " << ks_yld_7000_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_cor_data_7000_sys->GetBinContent(i);
  double error = ks_dNdy_cor_data_7000_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i-1] = sqrt(error*error + syserr*syserr);
  cout << ks_dNdy_cor_data_7000_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i-1] << endl;
}
ks_dNdy_cor_data_7000_sys->SetError(yent);

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
ks_dNdy_cor_data_900->GetYaxis()->SetTitle("N(pp #rightarrow K_{S}^{0}X) / Event_{NSD}");
ks_dNdy_cor_data_900->SetMaximum(0.8);
ks_dNdy_cor_data_900->SetMinimum(0.);
ks_dNdy_cor_data_900->SetLineColor(2);
ks_dNdy_cor_data_900->SetMarkerColor(2);
ks_dNdy_cor_data_900->SetMarkerStyle(21);
ks_dNdy_cor_data_900->SetMarkerSize(0.8);
TH1F *ks_dNdy_cor_data_900_sys = (TH1F*)ks_dNdy_cor_data_900->Clone("ks_dNdy_cor_data_700_sys");
cout << "ks y 900 integral / event: " << ks_yld_900_y << endl;
for (int i = 1; i <= 23; ++i) {
  double value = ks_dNdy_cor_data_900_sys->GetBinContent(i);
  double error = ks_dNdy_cor_data_900_sys->GetBinError(i);
  double syserr = 0.02*value;
  yent[i-1] = sqrt(error*error + syserr*syserr);
  cout << ks_dNdy_cor_data_900_sys->GetBinCenter(i) << " & " << value << " & " << error << " & " << syserr << " & " << yent[i-1] << endl;
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
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

c1->SaveAs(outFile_y_evt + ".root");
c1->SaveAs(outFile_y_evt + ".eps");
c1->SaveAs(outFile_y_evt + ".C");
c1->SaveAs(outFile_y_evt + ".png");

delete c1;



// y for Kshort per luminosity

ks_dNdy_cor_data_7000->Scale(trig7000/lum7000);
double ks_xsec_7000_y = calc_integral(ks_dNdy_cor_data_7000, 1.0);
cout << "ks y 7000 xsec: " << ks_xsec_7000_y << endl;
ks_dNdy_cor_data_7000_sys->GetYaxis()->SetTitle("#sigma(pp #rightarrow K_{S}^{0}X) [ mb ]");
ks_dNdy_cor_data_7000_sys->Scale(trig7000/lum7000);

ks_dNdy_cor_data_900->Scale(trig900/lum900);
double ks_xsec_900_y = calc_integral(ks_dNdy_cor_data_900, 1.0);
cout << "ks y 900 xsec: " << ks_xsec_900_y << endl;
ks_dNdy_cor_data_900_sys->GetYaxis()->SetTitle("#sigma(pp #rightarrow K_{S}^{0}X) [ mb ]");
ks_dNdy_cor_data_900_sys->Scale(trig900/lum900);
cout << "y 900\n";

TCanvas* c1 = new TCanvas("c1", "V0 y distribution");

ks_dNdy_cor_data_7000_sys->SetMaximum(43.0);
ks_dNdy_cor_data_7000_sys->SetMinimum(0.);
ks_dNdy_cor_data_7000_sys->Draw();
ks_dNdy_cor_data_7000->Draw("E1same");
ks_dNdy_cor_data_900_sys->Draw("same");
ks_dNdy_cor_data_900->Draw("E1same");
leg->Draw("same");
cmsprelim->DrawLatex(0.78,0.90,"CMS preliminary");

c1->SaveAs(outFile_y_lum + ".root");
c1->SaveAs(outFile_y_lum + ".eps");
c1->SaveAs(outFile_y_lum + ".C");
c1->SaveAs(outFile_y_lum + ".png");

delete c1;


}
