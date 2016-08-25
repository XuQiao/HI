/*
 * This script should be AcLic compiled - i.e.,
 *
 *   Root> .L test/Test.C++
 *
 * due to the nested loops
 */
/**
 * @file   Analyze.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr  9 10:18:10 2014
 * 
 * @brief  A ROOT script example. 
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
//#ifndef __CINT__
# include <correlations/Types.hh>
# include <correlations/Result.hh>
# include <correlations/QVector.hh>
# include <correlations/recursive/FromQVector.hh>
# include <correlations/recurrence/FromQVector.hh>
# include <correlations/closed/FromQVector.hh>
# include <TH1.h>
# include <TH2.h>
# include <TFile.h>
# include <TTree.h>
# include <TRandom.h>
# include <TProfile.h>
# include <TCanvas.h>
# include <TStopwatch.h>
# include <TRandom.h>
# include <TMath.h>
# include <TString.h>
# include <TLatex.h>
# include <TVectorD.h>
# include "par.h"
//#else
// To force loading of libraries
class TF1;
class TH1F;
class TProfile;
class TFile;
class TCanvas;
//#endif
void DrawInPad(TVirtualPad* p,
               Int_t sub,
               TH1* h,
               Bool_t logy=false)
{
  TVirtualPad* pp = p->cd(sub);
  pp->SetRightMargin(0.02);
  if (logy) pp->SetLogy();
  TH1* copy = h->DrawCopy("hist");
  copy->GetXaxis()->SetLabelSize(0.13);
  copy->GetYaxis()->SetLabelSize(0.08);
  copy->SetDirectory(0);
}
/**
 * A test
 *
 * @param maxH     Max harmonic
 * @param nEvents  Number of events
 * @param minN     Least multiplicity
 * @param maxN     Largest multiplicity
 */

void
Analyze(const TString mode="CLOSED",
        UShort_t       maxH=6,
        Bool_t         /*doLoops*/=false,
        TString filename="",
        int ifile=0   )
{
#ifdef __CINT__
  gROOT->LoadMacro("correlations/Types.hh++");
  gROOT->LoadMacro("correlations/Result.hh++");
  gROOT->LoadMacro("correlations/QVector.hh++");
  gROOT->LoadMacro("correlations/recursive/FromQVector.hh++");
  gROOT->LoadMacro("correlations/recurrence/FromQVector.hh++");
  gROOT->LoadMacro("correlations/closed/FromQVector.hh++");
  gROOT->LoadMacro("correlations/test/ReadData.hh++");
#endif
  // --- Setup of harmonics, etc -------------------------------------
  gRandom->SetSeed(54321);
  UShort_t emode = 0;
  if      (mode.EqualTo("closed",     TString::kIgnoreCase)) emode = 0;
  else if (mode.EqualTo("recurrence", TString::kIgnoreCase)) emode = 1;
  else if (mode.EqualTo("recursive",  TString::kIgnoreCase)) emode = 2;
  else
    Warning("Analyze", "Mode %s unknown, assuming CLOSED", mode.Data());

  correlations::QVector        q[nbin];
  correlations::FromQVector*   c[nbin];
  correlations::HarmonicVector h(maxH);
  for (UShort_t i = 0; i < maxH; i++) {
    // Generate random harmonicx
   // h[i] = -6 + gRandom->Integer(12);
   h[0] = 2;
   h[1] = -2;
   h[2] = 2;
   h[3] = -2;
   h[4] = 2;
   h[5] = -2;
   h[6] = 2;
   h[7] = -2;
   // Printf("h_%d:\t%d", i, h[i]);
  }

  // Resize the Q-vector to fit the harmonics
    for(int ibin=0;ibin<nbin;ibin++){
    q[ibin] = correlations::QVector(0,0,false);
    q[ibin].resize(h);
  switch (emode) {
  case 0: c[ibin] = new correlations::closed::FromQVector(q[ibin]); break;
  case 1: c[ibin] = new correlations::recurrence::FromQVector(q[ibin]); break;
  case 2: c[ibin] = new correlations::recursive::FromQVector(q[ibin]); break;
    }
  }
  //Printf("Correlator: %s", c->name());


  // --- Some histograms ---------------------------------------------
  TH1* sumreals[nbin];
  TH1* sumimags[nbin];
  TH1* weights[nbin];
  TVectorD tottrk;
  TVectorD Nevent;
  tottrk.ResizeTo(nbin);
  tottrk.Zero();
  Nevent.ResizeTo(nbin);
  Nevent.Zero();

  //TH1*      reals  = new TH1D("reals", "Re(C{n})", maxH-2+1, 2+.5, maxH+1+.5);
  //TH1*      imags  = static_cast<TH1*>(reals->Clone("imags"));
    for(int ibin=0;ibin<nbin;ibin++){
  sumreals[ibin]  = new TH1D(Form("sumreals_%d",ibin), "Re(C{n})", maxH-2+1, 2+.5, maxH+1+.5);
  sumimags[ibin]  = static_cast<TH1*>(sumreals[ibin]->Clone(Form("sumimags_%d",ibin)));
  weights[ibin]  = static_cast<TH1*>(sumreals[ibin]->Clone(Form("weights_%d",ibin)));
    }
  TProfile* timing = new TProfile("timing", "Timing", maxH-2+1, 2+.5,maxH+1+.5);
  TH1*      hs     = new TH1I("harmonics", "Harmonics", maxH, 1.5, maxH+1.5);
  /*reals->SetFillColor(kGreen+1);
  reals->SetFillStyle(3001);
  reals->SetStats(0);
  imags->SetTitle("Im(C{n})");
  imags->SetFillColor(kBlue+1);
  imags->SetFillStyle(3001);
  imags->SetStats(0);
  timing->SetFillColor(kRed+1);
  timing->SetFillStyle(3001);
  timing->SetStats(0);
  hs->SetFillColor(kMagenta+1);
  hs->SetFillStyle(3001);
  hs->SetStats(0);*/
  for (UShort_t i = 0; i < maxH-1; i++) {
    TString label = TString::Format("C{%d}", i+2);
//    reals->GetXaxis()->SetBinLabel(i+1, label);
//    imags->GetXaxis()->SetBinLabel(i+1, label);
    timing->GetXaxis()->SetBinLabel(i+1, label);
    hs->GetXaxis()->SetBinLabel(i+1,Form("h_{%d}", i+1));
    hs->SetBinContent(i+1, h[i]);
  }
  hs->GetXaxis()->SetBinLabel(maxH,Form("h_{%d}", maxH));
  hs->SetBinContent(maxH, h[maxH-1]);

  TStopwatch timer;
 
  // --- Setup input ------------------------------------------------
  TFile*   file = TFile::Open(Form("%s",filename.Data()), "READ");
  TFile*   trkEff_pp = TFile::Open("/home/xuq7/HI/flow/LYZ/v2/TrackCorrections_HIJING_538_OFFICIAL_Mar24.root");
  TTree*   tree = static_cast<TTree*>(file->Get("demo/TrackTree"));
  TH2F* rTotalEff3D = (TH2F*)trkEff_pp->Get("rTotalEff3D");
//  TArrayD  phis(0);
  Int_t mult, ntrk;
  double phi[10000],pt[10000],eta[10000];
//  TArrayD  weights(0);
//  Double_t phiR = 0;
//  TArrayD* pPhis = &phis;
//  TArrayD* pWeights = &weights;
  tree->SetBranchAddress("phi", phi);
  tree->SetBranchAddress("pt",pt);
  tree->SetBranchAddress("eta",eta);
  tree->SetBranchAddress("mult", &mult);
  tree->SetBranchAddress("ntrk", &ntrk);
//  tree->SetBranchAddress("weight", &pWeights);
//  tree->SetBranchAddress("event", &phiR);


  // --- The results -------------------------------------------------
  const UShort_t             nQ = maxH - 1;
  correlations::ResultVector qs[nbin];
    for(int ibin=0;ibin<nbin;ibin++)
        qs[ibin] = correlations::ResultVector(nQ);

  // --- Event loop --------------------------------------------------
  Int_t nEvents = tree->GetEntries();
  for (Int_t event = 0; event < nEvents; event++) {
    tree->GetEntry(event);
    int xbin=-1;
    for(int j=0;j<nbin;j++)
        if(ntrk<trkbin[j]&&ntrk>=trkbin[j+1])
            xbin=j;
        if(xbin<0 || xbin==nbin) continue;
    tottrk[xbin]+=ntrk;
    q[xbin].reset();
 //   printf("Event # %4u  %4d particles ", event++, phis.GetSize());
    for (UShort_t pa = 0; pa < mult; pa++){
	if(eta[pa]<-2.40||eta[pa]>2.40) continue;
        if(pt[pa]<ptmin||pt[pa]>ptmax) continue; //event selection
  //  phis.Set(n,pPhis);
       double weight = 1./rTotalEff3D->GetBinContent(rTotalEff3D->FindBin(eta[pa],pt[pa]));
      q[xbin].fill(phi[pa], weight);
    }
    for (UShort_t i = 0; i < nQ; i++) {
      UShort_t n = i + 2;
   //   printf("%s%d", i == 0 ? "" : "..", n);
      timer.Reset();
      timer.Start();
      qs[xbin][i] += c[xbin]->calculate(n, h);
      timer.Stop();
      timing->Fill(n+.5, timer.RealTime());
    }

  //  printf(" done\n");
  Nevent[xbin]++;
  }
  file->Close();

    for(int ibin=0;ibin<nbin;ibin++){
  for (UShort_t i = 0; i < nQ; i++) {
 //   UShort_t iq = i+2;
 //   Double_t              t  = timing->GetBinContent(i+1);
 //   correlations::Complex rc = qs[i].eval();
   // Printf("QC{%2d}: %12g + %12gi  <t>: %10gs",
//	   iq, rc.real(), rc.imag(), t);
   // if(i==0)Printf("v2{%2d}: %3g\n",2,sqrt(qs[0].eval().real()));
   // if(i==2)Printf("v2{%2d}: %3g\n",4,TMath::Power(fabs(qs[2].eval().real()),1./4));
   // if(i==4)Printf("v2{%2d}: %3g\n",6,TMath::Power(fabs(qs[4].eval().real()),1./6));
    sumreals[ibin]->SetBinContent(i+1,qs[ibin][i]._sum.real());
    sumimags[ibin]->SetBinContent(i+1,qs[ibin][i]._sum.imag());
    weights[ibin]->SetBinContent(i+1,qs[ibin][i]._weights);
    //reals->SetBinContent(i+1, rc.real());
    //imags->SetBinContent(i+1, rc.imag());
  }
    }

/*
  TCanvas* can = new TCanvas("C", "C");
  can->SetTopMargin(0.15);
  can->SetBottomMargin(0.15);
  can->SetRightMargin(0.03);
  can->Divide(1,3, 0, 0);

  DrawInPad(can, 3, timing, true);
  DrawInPad(can, 1, reals);
  DrawInPad(can, 2, imags);

  can->cd(0);
  TLatex* ltx = new TLatex(0.5,0.995,c->name());
  ltx->SetNDC(true);
  ltx->SetTextAlign(23);
  ltx->SetTextSize(0.04);
  ltx->Draw();

  can->Modified();
  can->Update();
  can->cd();
*/
  TString out(mode);
  out.ToLower();
  file = TFile::Open(Form("%s/%s_%d.root",outdir.Data(),out.Data(),ifile), "RECREATE");
    for(int ibin=0;ibin<nbin;ibin++){
  sumimags[ibin]->Write();
  sumreals[ibin]->Write();
  weights[ibin]->Write();
    }
  Nevent.Write("Nevent");
  tottrk.Write("tottrk");
  timing->Write();
  hs->Write();
  file->Write();
  file->Close();
    for(int ibin=0;ibin<nbin;ibin++){
  delete sumimags[ibin];
  delete sumreals[ibin];
  delete weights[ibin];
    }
  delete timing;
  delete hs;
}

#ifdef AS_PROG
#include <TApplication.h>
#include <TROOT.h>
#include <string>
#include <algorithm>

void usage(const char* progname)
{
  std::cout << "Usage: " << progname << " [OPTIONS]\n\n"
      << "Options:\n"
      << "  -h          This help\n"
      << "  -t MODE     How to analyze\n"
      << "  -n MAXH     Maximum harmonic to analyze\n"
      << "  -L          Also do nested loop calculations\n"
      << "  -l          Do not do nested loop calculations\n"
      << std::flush;
}

int
main(int argc, char** argv)
{
  bool           loops     = true;
  bool           batch     = false;
  unsigned short maxH      = 6;
  std::string    mode("CLOSED");

  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'h': usage(argv[0]); return 0;
      case 't': mode      = argv[++i];  break;
      case 'l': loops     = false; break;
      case 'L': loops     = true; break;
      case 'b': batch     = false; break;
      case 'B': batch     = true; break;
      case 'n': maxH      = atoi(argv[++i]); break;
      default:
	std::cerr << argv[0] << ": Unknown option " << argv[i] << std::endl;
	return 1;
      }
    }
  }

  TApplication app("app", 0, 0);
  if (batch) gROOT->SetBatch(true);

  Analyze1(mode.c_str(), maxH, loops);

  if (!batch) app.Run();
  return 0;
}
#endif
//
// EOF
//
