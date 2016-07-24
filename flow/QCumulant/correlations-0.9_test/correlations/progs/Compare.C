/**
 * @file   Compare.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr  9 10:19:47 2014
 * 
 * @brief  Compare outputs
 */
/*
 * Multi-particle correlations 
 * Copyright (C) 2013 K.Gulbrandsen, A.Bilandzic, C.H. Christensen.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses.
 */
#ifndef __CINT__
# include <TFile.h>
# include <TH1.h>
# include <TProfile.h>
# include <TCanvas.h>
# include <TError.h>
# include <TLatex.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TClass.h>
#else
class TFile;
class TH1;
class TProfile;
class TCanvas;
#endif

Bool_t
GetHistos(const TString& mode,
          TH1*& h,
          TH1*& real,
          TH1*& imag,
          TProfile*& time)
{
  TString in(mode);
  in.ToLower();
  if (!in.EndsWith(".root")) in.Append(".root");
  TFile* file = TFile::Open(in, "READ");
  if (!file) {
    Error("GetHistos", "Failed to open %s", in.Data());
    return false;
  }
  h    = static_cast<TH1*>(file->Get("harmonics"));
  real = static_cast<TH1*>(file->Get("reals"));
  imag = static_cast<TH1*>(file->Get("imags"));
  time = static_cast<TProfile*>(file->Get("timing"));
  if (!h || !real || !imag || !time) {
    Error("GetHistos", "Failed to find all histograms");
    return false;
  }
  real->SetDirectory(0);
  imag->SetDirectory(0);
  time->SetDirectory(0);
  h   ->SetDirectory(0);

  real->SetYTitle(mode);
  imag->SetYTitle(mode);
  time->SetYTitle(mode);
  h   ->SetYTitle(mode);

  real->Sumw2();
  imag->Sumw2();
  time->Sumw2();
  h   ->Sumw2();

  file->Close();
  return true;
}

void DrawTwoInPad(TVirtualPad* p,
                  Int_t        sub,
                  TH1*         h1,
                  TH1*         h2,
		  Bool_t       ratio,
		  Bool_t       logy=false,
		  Bool_t       legend=false)
{
  TVirtualPad* pp = p->cd(sub);
  pp->SetRightMargin(0.02);
  pp->SetLeftMargin(0.10);
  TVirtualPad* ppp = pp;
  if (ratio) {
    pp->Divide(1,2,0,0);
    ppp = pp->cd(1);				
    ppp->SetRightMargin(0.02);
  }
  if (logy) ppp->SetLogy();
  TH1* hs[] = { h1, h2, 0 };
  if (h1->GetMaximum() < h2->GetMaximum()) {
    hs[0] = h2;
    hs[1] = h1;
  }
  TH1** ph = hs;
  Double_t size = (ratio ? 0.1 : 0.05);
  Double_t off  = (ratio ? 0.6 : 0.5);
  h1->SetFillStyle(3004);
  h2->SetFillStyle(3005);
  while (*ph) {
    TString opt("hist");
    if (ph != hs) opt.Append(" same");

    TH1* copy = (*ph)->DrawCopy(opt);
    copy->GetXaxis()->SetLabelSize(2*size);
    copy->GetYaxis()->SetLabelSize(size);
    copy->GetYaxis()->SetTitleSize(size);
    copy->GetYaxis()->SetTitleOffset(off);
    copy->SetYTitle(copy->GetTitle());
    copy->SetTitle("");
    copy->SetDirectory(0);
    ph++;
  }
  TString s1 = h1->GetYaxis()->GetTitle();
  TString s2 = h2->GetYaxis()->GetTitle();
  
  if (legend) { 
    TLegend* l = new TLegend(0.6, 0.1, 0.9, 0.9);
    l->SetBorderSize(0);
    TLegendEntry* e = l->AddEntry("dummy", s1, "lf");
    l->SetFillColor(kWhite);
    e->SetFillColor(kBlack);
    e->SetFillStyle(h1->GetFillStyle());
    e = l->AddEntry("dummy", s2, "lf");
    e->SetFillColor(kBlack);
    e->SetFillStyle(h2->GetFillStyle());
    l->Draw();
  }
  if (!ratio) return;
  ppp = pp->cd(2);
  ppp->SetRightMargin(0.02);
  TH1* r = static_cast<TH1*>(h1->Clone(Form("ratio%s", h1->GetName())));
  r->SetDirectory(0);
  r->SetTitle("");
  r->GetXaxis()->SetLabelSize(size);
  r->GetYaxis()->SetLabelSize(size);
  r->GetYaxis()->SetTitleSize(0.9*size);
  r->GetYaxis()->SetTitleOffset(0.9*off);
  r->SetMarkerStyle(20);
  r->SetMarkerColor(h1->GetFillColor()+1);
  r->SetFillStyle(3007);
  r->SetYTitle(Form("#frac{%s}{%s}", s1.Data(), s2.Data()));

  // r->Add(h2, -1);
  // r->Divide(h1);
  if (!r->IsA()->InheritsFrom(TProfile::Class())) {
    r->GetSumw2()->Set(0); // r->Sumw2(false);
    h2->GetSumw2()->Set(0); // h2->Sumw2(false);
  }
  r->Divide(h2);
  Printf("%s", r->GetName());
  for (UShort_t bin = 1; bin <= r->GetNbinsX(); bin++) {
    Printf("  bin # %2d: Diff=%g+/-%g", bin, r->GetBinContent(bin),
	   r->GetBinError(bin));
    r->SetBinError(bin, 0);
  }
  r->GetSumw2()->Set(0); //r->Sumw2(false);
  r->SetMarkerSize(4);
  r->SetMaximum(r->GetMaximum()*1.2);
  r->SetMinimum(r->GetMinimum()*0.8);
  r->Draw("hist text30");  
  p->Modified();
  p->Update();
  p->cd();
}

void
Compare(const TString& mode1,
        const TString& mode2)
{
  TH1* h1 = 0;
  TH1* h2 = 0;
  TH1* r1 = 0;
  TH1* i1 = 0;
  TH1* r2 = 0;
  TH1* i2 = 0;
  TProfile* t1 = 0;
  TProfile* t2 = 0;

  if (!GetHistos(mode1, h1, r1, i1, t1)) return;
  if (!GetHistos(mode2, h2, r2, i2, t2)) return;

  TCanvas* can = new TCanvas("comp", "Comparison", 1200, 800);
  can->SetTopMargin(0.15);
  can->SetBottomMargin(0.15);
  can->SetRightMargin(0.03);
  can->SetLeftMargin(0.03);
  can->Divide(2, 2);
  
  DrawTwoInPad(can, 4, t1, t2, true, true);
  DrawTwoInPad(can, 3, i1, i2, true, false, true);
  DrawTwoInPad(can, 2, r1, r2, true);
  DrawTwoInPad(can, 1, h1, h2, false);

  can->Modified();
  can->Update();
  can->cd();

  can->Print(Form("%s_%s.png", mode1.Data(), mode2.Data()));
}

#ifdef AS_PROG
#include <TROOT.h>
#include <TApplication.h>
#include <iostream>
void
usage(const char* progname)
{
  std::cout << "Usage: " << progname << " [OPTIONS]\n\n" << "Options:\n"
            << "  -h          This help\n" << "  -1 NAME     First set\n"
            << "  -2 NAME     Second set\n" << std::flush;
}
int
main(int argc,
     char** argv)
{
  TString a("closed");
  TString b("recursive");
  bool batch = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'h':
        usage(argv[0]);
        return 0;
      case '1':
        a = argv[++i];
        break;
      case '2':
        b = argv[++i];
        break;
      case 'b':
        batch = false;
        break;
      case 'B':
        batch = true;
        break;
      default:
        std::cerr << argv[0] << ": Unknown option " << argv[i] << std::endl;
        return 1;
      }
    }
  }
  TApplication app("app", 0, 0);
  if (batch) gROOT->SetBatch(true);

  Compare(a, b);

  if (!batch) app.Run();
  return 0;
}
#endif
//
// EOF
//
