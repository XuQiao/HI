{
/*#include <tdrstyle.C>
#include <sstream>
//#include <iostream>
#include <string>

gROOT->Reset();
setTDRStyle();
gStyle->SetLabelOffset(0.017, "X");
gStyle->SetPadTopMargin(0.06);
gStyle->SetPadBottomMargin(0.12);
gStyle->SetStatFormat("5.5g");
gStyle->SetStatY(0.94);
gStyle->SetFitFormat("5.5g");
//gStyle->SetStatFormat("NELU",AutoPrecision(3));
*/
TFile* fin = TFile::Open("/scratch/xuq7/strange_particles/v0analysis.root", "READ");
fin->cd("ksAnalyzer");
TCanvas *c1;
TH1F *histo = (TH1F*)ksMass_eta_2->Clone();

histo->Draw();
c1->Print("ksMass_eta_2.gif");
TH1F *histo = (TH1F*)ksMass_eta_-2->Clone();
ksMass_eta_-2->Draw();
c1->Print("ksMass_eta_-2.gif");

}
