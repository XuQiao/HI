/**
 * @file   Write.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr  9 10:20:34 2014
 * 
 * @brief  Write output data to a ROOT tree
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
# include <TH1D.h>
# include <TRandom.h>
# include <TF1.h>
# include <TTree.h>
# include <TFile.h>
# include <TMath.h>
#else
class TH1D;
class TRandom;
class TTree;
class TFile;
#endif

Double_t func(Double_t* xx,
              Double_t* pp)
{
  Double_t phi = xx[0];
  Double_t psi = pp[0];
  Double_t ret = 1;
  for (UShort_t i = 1; i < 7; i++)
    ret += 2 * pp[i] * TMath::Cos(i * (phi-psi));
  return ret;
}

void Write(UShort_t nEvents, UShort_t minN, UShort_t maxN)
{
  TFile*   out        = TFile::Open("data.root", "RECREATE");
  TTree*   tree       = new TTree("T", "T");
  TF1*     phiFunc    = new TF1("phiFunc",&func,0,TMath::TwoPi(),7);
  TH1*     weightDist = new TH1D("weights", "Weights",72,0,TMath::TwoPi());
  TH1*     phiDist    = new TH1D("phis", "#phi-#phi_{R}", 72,0,TMath::TwoPi());
  TArrayD  phis(0);
  TArrayD  weights(0);
  Double_t phiR = 0;

  tree->Branch("phi", &phis);
  tree->Branch("weight", &weights);
  tree->Branch("event", &phiR, "phiR/D");

 // phiFunc->SetParameters(0, 0.01, 0.05, 0.03, 0.02, 0.01, 0.01);
  phiFunc->SetParameters(0, 0., 0.05, 0., 0., 0., 0.);
  for (UShort_t bin = 1; bin <= weightDist->GetNbinsX(); bin++)
    weightDist->SetBinContent(bin, 1 + gRandom->Uniform(-0.144, +0.144));
  weightDist->SetFillColor(kGreen+1);
  weightDist->SetFillStyle(3001);
  weightDist->SetXTitle("#phi");
  weightDist->SetYTitle("w");

  phiDist->SetFillColor(kRed+1);
  phiDist->SetFillStyle(3001);
  phiDist->SetXTitle("#phi - #phi_{R}");
  phiDist->SetYTitle("1/N dN/d(#phi-#phi_{R})");
  phiDist->Sumw2();

  for (UShort_t ev = 0; ev < nEvents; ev++) {
     // printf("%d\n",ev);
    UShort_t mult = UShort_t(gRandom->Uniform(minN,maxN));
    phiR          = gRandom->Uniform(0, TMath::TwoPi());
    phiFunc->SetParameter(0, phiR);
    phis.Set(mult);
    weights.Set(mult);
    phis.Reset(0);
    weights.Reset(0);

    for (UShort_t pa = 0; pa < mult; pa++) {
      Double_t phi    = phiFunc->GetRandom();
      Double_t weight = weightDist->GetBinContent(weightDist->FindBin(phi));
      phis[pa]        = phi;
      weights[pa]     = weight;
      Double_t fill   = phi - phiR;
      if (fill < 0) fill += TMath::TwoPi();
      phiDist->Fill(fill, weight/mult);
    }

    tree->Fill();
  }
  phiFunc->SetParameter(0, 0);
  phiDist->Scale(phiFunc->Integral(0,TMath::TwoPi()) / nEvents, "width");
  phiFunc->Write();
  out->Write();
  out->Close();
}

#ifdef AS_PROG
void usage(const char* progname)
{
  std::cout << "Usage: " << progname << " [OPTIONS]\n\n"
      << "Options:\n"
      << "  -h          This help\n"
      << "  -e NEVENTS  Set number of events\n"
      << "  -m MINIMUM  Set least multiplicity of each event\n"
      << "  -M MAXIMUM  Set largest multiplicity of each event\n"
      << std::flush;
}
int
main(int argc, char** argv)
{
  UShort_t nEvents = 100;
  UShort_t minN    = 80;
  UShort_t maxN    = 100;

  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'h': usage(argv[0]); return 0;
      case 'e': nEvents   = atoi(argv[++i]); break;
      case 'm': minN      = atoi(argv[++i]); break;
      case 'M': maxN      = atoi(argv[++i]); break;
      default:
        std::cerr << argv[0] << ": Unknown option " << argv[i] << std::endl;
        return 1;
      }
    }
  }
  Write(nEvents, minN, maxN);

  return 0;
}
#endif
//
// EOF
//
