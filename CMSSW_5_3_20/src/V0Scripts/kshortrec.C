#ifndef kshortrec_cxx
#define kshortrec_cxx
#include "kshortrec.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
//#include "constants.h"

using namespace std;

vector<double> ksptBinMinima;
vector<double> lamptBinMinima;
vector<double> kspBinMinima;
vector<double> lampBinMinima;
vector<double> kslifeBinMinima;
vector<double> lamlifeBinMinima;
vector<double> ntrkBinMinima;

kshortrec::kshortrec(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  /*if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fin);
      if (!f) {
	f = new TFile(fin);
	if ( kshort ) {
	  f->cd("analyzeKshort");
	} else {
	  f->cd("analyzeLambda");
	}
      }
      tree = (TTree*)gDirectory->Get("ntuple");

      }*/
   Init(tree);

   for ( int ndx = 0; ndx <= ksptNbins; ndx++ ) {
     ksptBinMinima.push_back( ksptarray[ndx] );
   }

   for ( int ndx = 0; ndx <= lamptNbins; ndx++ ) {
     lamptBinMinima.push_back( lamptarray[ndx] );
   }

   for ( int ndx = 0; ndx <= kspNbins; ndx++ ) {
     kspBinMinima.push_back( ksparray[ndx] );
   }

   for ( int ndx = 0; ndx <= lampNbins; ndx++ ) {
     lampBinMinima.push_back( lamparray[ndx] );
   }

   for ( int ndx = 0; ndx <= kslifeNbins; ndx++ ) {
     kslifeBinMinima.push_back( kslifearray[ndx] );
   }

   for ( int ndx = 0; ndx <= lamlifeNbins; ndx++ ) {
     lamlifeBinMinima.push_back( lamlifearray[ndx] );
   }

   for ( int ndx = 0; ndx <= ntrkNbins; ndx++ ) {
     ntrkBinMinima.push_back( ntrkarray[ndx] );
   }
}

kshortrec::~kshortrec()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t kshortrec::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t kshortrec::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
  //cout << "Loading tree" << endl;
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   //cout << "Loaded, centry = " << centry << endl;
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  {
     //cout << "Here it is." << endl;
     return centry;
   }
   //cout << "Checking chain.." << endl;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
     //cout << "Changing current tree.." << endl;
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void kshortrec::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   v0VtxX = 0;
   v0VtxY = 0;
   v0VtxZ = 0;
   v0VtxR = 0;
   v0VtxNormChi2 = 0;
   v0pX = 0;
   v0pY = 0;
   v0pZ = 0;
   v0pT = 0;
   v0Eta = 0;
   v0Phi = 0;
   v0CandMass = 0;
   v0OtherCandMass = 0;
   v0Lifetime = 0;
   v03DLifetime = 0;
   v03DLifetimeError = 0;
   v03DLifetimeScaled = 0;
   v03DLifetimeScaledError = 0;
   v03DLifetimeMPV = 0;
   v03DLifetimeMPVError = 0;
   v0PDG = 0;
   v0VtxSig = 0;
   v0CosThetaPAndLineOfFlight = 0;
   v0PCAPrimary = 0;
   v0PCAPrimaryError = 0;
   v0PCABeamSpot = 0;
   v0PCABeamSpotError = 0;
   v3dIpWrtPrimary = 0;
   v3dIpWrtPrimarySig = 0;
   thePosDau3DIp = 0;
   thePosDau3DIpSig = 0;
   theNegDau3DIp = 0;
   theNegDau3DIpSig = 0;
   posDauNhits = 0;
   negDauNhits = 0;
   priVtxWasRefit = 0;
   priRefitVtxNdof = 0;
   priRefitVtxX = 0;
   priRefitVtxY = 0;
   priRefitVtxZ = 0;
   v0VtxPriDist = 0;
   v0VtxPriDistErr = 0;
   v0VtxBSTransDist = 0;
   v0VtxBSTransDistErr = 0;
   v0PosDaupX = 0;
   v0PosDaupY = 0;
   v0PosDaupZ = 0;
   v0PosDauAtVtxpX = 0;
   v0PosDauAtVtxpY = 0;
   v0PosDauAtVtxpZ = 0;
   v0PosDauEta = 0;
   v0PosDauNValidHits = 0;
   v0PosDauNPixelHits = 0;
   v0PosDauNStripHits = 0;
   v0PosDauChi2 = 0;
   v0PosDauNormChi2 = 0;
   v0PosDauD0 = 0;
   v0DaupT_lowestpT = 0;
   v0NegDaupX = 0;
   v0NegDaupY = 0;
   v0NegDaupZ = 0;
   v0NegDauAtVtxpX = 0;
   v0NegDauAtVtxpY = 0;
   v0NegDauAtVtxpZ = 0;
   v0NegDauEta = 0;
   v0NegDauNValidHits = 0;
   v0NegDauNPixelHits = 0;
   v0NegDauNStripHits = 0;
   v0NegDauChi2 = 0;
   v0NegDauNormChi2 = 0;
   v0NegDauD0 = 0;
   v0DauEta_lowestpT = 0;
   v0DauPhi_lowestpT = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("v0VtxX", &v0VtxX, &b_v0VtxX);
   fChain->SetBranchAddress("v0VtxY", &v0VtxY, &b_v0VtxY);
   fChain->SetBranchAddress("v0VtxZ", &v0VtxZ, &b_v0VtxZ);
   fChain->SetBranchAddress("v0VtxR", &v0VtxR, &b_v0VtxR);
   fChain->SetBranchAddress("v0VtxNormChi2", &v0VtxNormChi2, &b_v0VtxNormChi2);
   fChain->SetBranchAddress("v0pX", &v0pX, &b_v0pX);
   fChain->SetBranchAddress("v0pY", &v0pY, &b_v0pY);
   fChain->SetBranchAddress("v0pZ", &v0pZ, &b_v0pZ);
   fChain->SetBranchAddress("v0pT", &v0pT, &b_v0pT);
   fChain->SetBranchAddress("v0Eta", &v0Eta, &b_v0Eta);
   fChain->SetBranchAddress("v0Phi", &v0Phi, &b_v0Phi);
   fChain->SetBranchAddress("v0CandMass", &v0CandMass, &b_v0CandMass);
   fChain->SetBranchAddress("v0OtherCandMass", &v0OtherCandMass, &b_v0OtherCandMass);
   fChain->SetBranchAddress("v0Lifetime", &v0Lifetime, &b_v0Lifetime);
   fChain->SetBranchAddress("v03DLifetime", &v03DLifetime, &b_v03DLifetime);
   fChain->SetBranchAddress("v03DLifetimeError", &v03DLifetimeError, &b_v03DLifetimeError);
   fChain->SetBranchAddress("v03DLifetimeScaled", &v03DLifetimeScaled, &b_v03DLifetimeScaled);
   fChain->SetBranchAddress("v03DLifetimeScaledError", &v03DLifetimeScaledError, &b_v03DLifetimeScaledError);
   fChain->SetBranchAddress("v03DLifetimeMPV", &v03DLifetimeMPV, &b_v03DLifetimeMPV);
   fChain->SetBranchAddress("v03DLifetimeMPVError", &v03DLifetimeMPVError, &b_v03DLifetimeMPVError);
   fChain->SetBranchAddress("v0PDG", &v0PDG, &b_v0PDG);
   fChain->SetBranchAddress("v0VtxSig", &v0VtxSig, &b_v0VtxSig);
   fChain->SetBranchAddress("v0CosThetaPAndLineOfFlight", &v0CosThetaPAndLineOfFlight, &b_v0CosThetaPAndLineOfFlight);
   fChain->SetBranchAddress("v0PCAPrimary", &v0PCAPrimary, &b_v0PCAPrimary);
   fChain->SetBranchAddress("v0PCAPrimaryError", &v0PCAPrimaryError, &b_v0PCAPrimaryError);
   fChain->SetBranchAddress("v0PCABeamSpot", &v0PCABeamSpot, &b_v0PCABeamSpot);
   fChain->SetBranchAddress("v0PCABeamSpotError", &v0PCABeamSpotError, &b_v0PCABeamSpotError);
   fChain->SetBranchAddress("v3dIpWrtPrimary", &v3dIpWrtPrimary, &b_v3dIpWrtPrimary);
   fChain->SetBranchAddress("v3dIpWrtPrimarySig", &v3dIpWrtPrimarySig, &b_v3dIpWrtPrimarySig);
   fChain->SetBranchAddress("thePosDau3DIp", &thePosDau3DIp, &b_thePosDau3DIp);
   fChain->SetBranchAddress("thePosDau3DIpSig", &thePosDau3DIpSig, &b_thePosDau3DIpSig);
   fChain->SetBranchAddress("theNegDau3DIp", &theNegDau3DIp, &b_theNegDau3DIp);
   fChain->SetBranchAddress("theNegDau3DIpSig", &theNegDau3DIpSig, &b_theNegDau3DIpSig);
   fChain->SetBranchAddress("posDauNhits", &posDauNhits, &b_posDauNhits);
   fChain->SetBranchAddress("negDauNhits", &negDauNhits, &b_negDauNhits);
   fChain->SetBranchAddress("priVtxChi2", &priVtxChi2, &b_priVtxChi2);
   fChain->SetBranchAddress("privtxNdof", &privtxNdof, &b_priVtxNdof);
   fChain->SetBranchAddress("priVtxX", &priVtxX, &b_priVtxX);
   fChain->SetBranchAddress("priVtxY", &priVtxY, &b_priVtxY);
   fChain->SetBranchAddress("priVtxZ", &priVtxZ, &b_priVtxZ);
   fChain->SetBranchAddress("priVtxNtracks", &priVtxNtracks, &b_priVtxNtracks);
   fChain->SetBranchAddress("priVtxIsFake", &priVtxIsFake, &b_priVtxIsFake);
   fChain->SetBranchAddress("priVtxWasRefit", &priVtxWasRefit, &b_priVtxWasRefit);
   fChain->SetBranchAddress("nLooseTracks", &nLooseTracks, &b_nLooseTracks);
   fChain->SetBranchAddress("nHighPurTracks", &nHighPurTracks, &b_nHighPurTracks);
   fChain->SetBranchAddress("nLooseTracks_pt_500", &nLooseTracks_pt_500, &b_nLooseTracks_pt_500);
   fChain->SetBranchAddress("nHighPurTracks_pt_500", &nHighPurTracks_pt_500, &b_nHighPurTracks_pt_500);
   fChain->SetBranchAddress("nV0s", &nV0s, &b_nV0s);
   fChain->SetBranchAddress("trigHF", &trigHF, &b_trigHF);
   fChain->SetBranchAddress("trigTech40_41", &trigTech40_41, &b_trigTech40_41);
   fChain->SetBranchAddress("trigTech34", &trigTech34, &b_trigTech34);
   fChain->SetBranchAddress("trigTech36_39", &trigTech36_39, &b_trigTech36_39);
   fChain->SetBranchAddress("trigScraping", &trigScraping, &b_trigScraping);
   fChain->SetBranchAddress("trigAlgo124", &trigAlgo124, &b_trigAlgo124);
   fChain->SetBranchAddress("trigHLTminBias", &trigHLTminBias, &b_trigHLTminBias);
   fChain->SetBranchAddress("runNum", &runNum, &b_runNum);
   fChain->SetBranchAddress("evtNum", &evtNum, &b_evtNum);
   fChain->SetBranchAddress("processType", &processType, &b_processType);
   fChain->SetBranchAddress("priRefitVtxNdof", &priRefitVtxNdof, &b_priRefitVtxNdof);
   fChain->SetBranchAddress("priRefitVtxX", &priRefitVtxX, &b_priRefitVtxX);
   fChain->SetBranchAddress("priRefitVtxY", &priRefitVtxY, &b_priRefitVtxY);
   fChain->SetBranchAddress("priRefitVtxZ", &priRefitVtxZ, &b_priRefitVtxZ);
   fChain->SetBranchAddress("v0VtxPriDist", &v0VtxPriDist, &b_v0VtxPriDist);
   fChain->SetBranchAddress("v0VtxPriDistErr", &v0VtxPriDistErr, &b_v0VtxPriDistErr);
   fChain->SetBranchAddress("v0VtxBSTransDist", &v0VtxBSTransDist, &b_v0VtxBSTransDist);
   fChain->SetBranchAddress("v0VtxBSTransDistErr", &v0VtxBSTransDistErr, &b_v0VtxBSTransDistErr);
   fChain->SetBranchAddress("v0PosDaupX", &v0PosDaupX, &b_v0PosDaupX);
   fChain->SetBranchAddress("v0PosDaupY", &v0PosDaupY, &b_v0PosDaupY);
   fChain->SetBranchAddress("v0PosDaupZ", &v0PosDaupZ, &b_v0PosDaupZ);
   fChain->SetBranchAddress("v0PosDauAtVtxpX", &v0PosDauAtVtxpX, &b_v0PosDauAtVtxpX);
   fChain->SetBranchAddress("v0PosDauAtVtxpY", &v0PosDauAtVtxpY, &b_v0PosDauAtVtxpY);
   fChain->SetBranchAddress("v0PosDauAtVtxpZ", &v0PosDauAtVtxpZ, &b_v0PosDauAtVtxpZ);
   fChain->SetBranchAddress("v0PosDauEta", &v0PosDauEta, &b_v0PosDauEta);
   fChain->SetBranchAddress("v0PosDauNValidHits", &v0PosDauNValidHits, &b_v0PosDauNValidHits);
   fChain->SetBranchAddress("v0PosDauNPixelHits", &v0PosDauNPixelHits, &b_v0PosDauNPixelHits);
   fChain->SetBranchAddress("v0PosDauNStripHits", &v0PosDauNStripHits, &b_v0PosDauNStripHits);
   fChain->SetBranchAddress("v0PosDauChi2", &v0PosDauChi2, &b_v0PosDauChi2);
   fChain->SetBranchAddress("v0PosDauNormChi2", &v0PosDauNormChi2, &b_v0PosDauNormChi2);
   fChain->SetBranchAddress("v0PosDauD0", &v0PosDauD0, &b_v0PosDauD0);
   fChain->SetBranchAddress("v0DaupT_lowestpT", &v0DaupT_lowestpT, &b_v0DaupT_lowestpT);
   fChain->SetBranchAddress("v0NegDaupX", &v0NegDaupX, &b_v0NegDaupX);
   fChain->SetBranchAddress("v0NegDaupY", &v0NegDaupY, &b_v0NegDaupY);
   fChain->SetBranchAddress("v0NegDaupZ", &v0NegDaupZ, &b_v0NegDaupZ);
   fChain->SetBranchAddress("v0NegDauAtVtxpX", &v0NegDauAtVtxpX, &b_v0NegDauAtVtxpX);
   fChain->SetBranchAddress("v0NegDauAtVtxpY", &v0NegDauAtVtxpY, &b_v0NegDauAtVtxpY);
   fChain->SetBranchAddress("v0NegDauAtVtxpZ", &v0NegDauAtVtxpZ, &b_v0NegDauAtVtxpZ);
   fChain->SetBranchAddress("v0NegDauEta", &v0NegDauEta, &b_v0NegDauEta);
   fChain->SetBranchAddress("v0NegDauNValidHits", &v0NegDauNValidHits, &b_v0NegDauNValidHits);
   fChain->SetBranchAddress("v0NegDauNPixelHits", &v0NegDauNPixelHits, &b_v0NegDauNPixelHits);
   fChain->SetBranchAddress("v0NegDauNStripHits", &v0NegDauNStripHits, &b_v0NegDauNStripHits);
   fChain->SetBranchAddress("v0NegDauChi2", &v0NegDauChi2, &b_v0NegDauChi2);
   fChain->SetBranchAddress("v0NegDauNormChi2", &v0NegDauNormChi2, &b_v0NegDauNormChi2);
   fChain->SetBranchAddress("v0NegDauD0", &v0NegDauD0, &b_v0NegDauD0);
   fChain->SetBranchAddress("v0DauEta_lowestpT", &v0DauEta_lowestpT, &b_v0DauEta_lowestpT);
   fChain->SetBranchAddress("v0DauPhi_lowestpT", &v0DauPhi_lowestpT, &b_v0DauPhi_lowestpT);
   Notify();
}

Bool_t kshortrec::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void kshortrec::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t kshortrec::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  return 1;
}

void kshortrec::Loop(int kshort,
		     //int ksptNbins,
		     //int kspNbins,
		     //int kslifeNbins,
		     //double ksyXmin, double ksyXmax, int ksyNbins,
		     //double ksyXmin_signed, double ksyXmax_signed, int ksyNbins_signed,
		     //int lamptNbins,
		     //int lampNbins,
		     //int lamlifeNbins,
		     //double lamyXmin, double lamyXmax, int lamyNbins,
		     //double lamptXmin, 
		     TH1F* nV0sH,
		     TH1F* numvee,
		     TH1F* ntrk,
		     TH1F* vee_nocuts,
		     TH1F* veeReg_nocuts,
		     TH1F* veeBar_nocuts,
		     TH1F* vee_withcuts,
		     TH1F* veeReg_withcuts,
		     TH1F* veeBar_withcuts,
		     //float* ksptarray,
		     //float* ksparray,
		     //float* kslifearray,
		     TH1F** ksPtrs,
		     TH1F** ksPtrs_signed,
		     TH1F** ksptPtrs,
		     TH1F** kspPtrs,
		     TH1F** kslifePtrs,
		     //float* lamptarray,
		     //float* lamparray,
		     //float* lamlifearray,
		     TH1F** lamPtrs,
		     TH1F** lamRegPtrs,
		     TH1F** lamBarPtrs,
		     TH1F** lamptPtrs,
		     TH1F** lamRegptPtrs,
		     TH1F** lamBarptPtrs,
		     TH1F** lampPtrs,
		     TH1F** lamlifePtrs)
{
//   In a ROOT session, you can do:
//      Root > .L kshortrec.C
//      Root > kshortrec t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  const Double_t ksMassConst = 0.497614;
  const Double_t lamMassConst = 1.115683;
  const Double_t ksMassConst2 = ksMassConst * ksMassConst;
  const Double_t lamMassConst2 = lamMassConst * lamMassConst;
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   int iybin, iybin_signed, iptbin, ipbin, ilifebin;
   double p, pt, //ksy, 
     lifetime, mass;
   //   cout << "Check incoming values: " << ksYXmin << ", " << ksYBinWidth << ", " << ksptNbins << endl;
   cout << "nentries = " << nentries << endl;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     printf("\r%3d%%", (int)(100L * jentry / nentries));
     fflush(stdout);
     //cout << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //      cout << "nb, nbytes = " << nb << ", " << nbytes << endl;
      //      if ( nV0s > 0 ) cout << "number of v0s = " << nV0s << endl;
      nV0sH->Fill(nV0s);
      //      if ( priVtxChi2 <= 0.0 ) continue;
      int ncand = (*v0CandMass).size();
      numvee->Fill(ncand);
      if ( !trigTech34 || !trigScraping || !trigHF || !trigTech36_39 
	   || !trigHLTminBias || priVtxIsFake ) continue;
      ntrk->Fill(nLooseTracks);
      for ( int i = 0; i<ncand; i++ ) {
	int pdgId = (*v0PDG)[i];
	mass = (*v0CandMass)[i];
	vee_nocuts->Fill( mass );
	if( pdgId == 3122 ) veeReg_nocuts->Fill(mass);
	else if( pdgId == -3122 ) veeBar_nocuts->Fill(mass);
	if ( fabs((*v3dIpWrtPrimarySig)[i]) >= 3.0 ) continue;
	if ( kshort && fabs((*v0OtherCandMass)[i]-1.115683)<0.0075 ) continue;
	if ( !kshort && fabs((*v0OtherCandMass)[i]-0.497614)<0.020 ) continue;
	if ( (*v0VtxSig)[i] <= 10.0 ) continue;
	if ( (*v0NegDauNormChi2)[i]<0.0 || (*v0PosDauNormChi2)[i] < 0.0 ) continue;
	lifetime = (*v03DLifetime)[i];
	//	if ( lifetime < 0.02 ) continue;
	//eta = fabs((*v0Eta)[i]);
	//if ( eta >= 2.5 ) continue;
	double p2 = ((*v0pZ)[i]) * ((*v0pZ)[i]) + ((*v0pT)[i]) * ((*v0pT)[i]);
	double pz = (*v0pZ)[i];
	//int pdgId = (*v0PDG)[i];
	//cout << ksy << ", " << lamy << endl;
	pt = (*v0pT)[i];
	p = sqrt( p2 );
	if (kshort) {
	  double ksenergy = sqrt( p2 + ksMassConst2 );
	  double ksy_signed = (0.5*log((ksenergy + pz) / (ksenergy - pz)));
	  double ksy = fabs(ksy_signed);
	  if( ksy > ksyXmax ) continue;
	  vee_withcuts->Fill( mass );
	  iybin = (int) ((ksy-ksyXmin)/ksyBinWidth);
	  iybin_signed = (int) ((ksy_signed - ksyXmin_signed) / ksyBinWidth_signed);
	  for ( iptbin = 0; iptbin < ksptNbins; iptbin++ ) {
	    if ( pt < ksptarray[iptbin+1] ) break;
	  }
	  for ( ipbin = 0; ipbin < kspNbins; ipbin++ ) {
	    if ( p < ksparray[ipbin+1] ) break;
	  }
	  for ( ilifebin = 0; ilifebin < kslifeNbins; ilifebin++ ) {
	    if ( lifetime < kslifearray[ilifebin+1] ) break;
	  }
	  ksPtrs[iybin]->Fill(mass,1.);
	  ksPtrs_signed[iybin_signed]->Fill(mass, 1.);
	  if ( ( pt >= ksptarray[0] ) && ( pt < ksptarray[ksptNbins] ) ) {
	    ksptPtrs[iptbin]->Fill(mass,1.);
	  }
	  if ( ( p >= ksparray[0] ) && ( p < ksparray[kspNbins] ) ) {
	    kspPtrs[ipbin]->Fill(mass,1.);
	  }
	  if ( ( lifetime >= kslifearray[0] ) && ( lifetime < kslifearray[kslifeNbins] ) ) {
	    kslifePtrs[ilifebin]->Fill(mass,1.);
	  }
	} 
	else {//lambda
	  //cout << "In lam part." << endl;
	  double lamenergy = sqrt( p2 + lamMassConst2 );
	  double lamy_signed = (0.5*log((lamenergy + pz) / (lamenergy - pz)));
	  double lamy = fabs(lamy_signed);
	  if( lamy > lamyXmax ) continue;
	  vee_withcuts->Fill( mass );
	  if( pdgId > 0 ) veeReg_withcuts->Fill(mass);
	  else if( pdgId < 0 ) veeBar_withcuts->Fill(mass);
	  iybin = (int) ((lamy-lamyXmin)/lamyBinWidth);
	  //iybin_signed = (int)((lamy_signed - lamyXmin_signed) / lamyBinWidth_signed);
	  for ( iptbin = 0; iptbin < lamptNbins; iptbin++ ) {
	    if ( pt < lamptarray[iptbin+1] ) break;
	  }
	  for ( ipbin = 0; ipbin < lampNbins; ipbin++ ) {
	    if ( p < lamparray[ipbin+1] ) break;
	  }
	  for ( ilifebin = 0; ilifebin < lamlifeNbins; ilifebin++ ) {
	    if ( lifetime < lamlifearray[ilifebin+1] ) break;
	  }
	  lamPtrs[iybin]->Fill(mass,1.);
	  //lamPtrs_signed[iybin_signed]->Fill(mass, 1.);
	  if( pdgId > 0 ) {
	    lamRegPtrs[iybin]->Fill(mass, 1.);
	    //lamRegPtrs_signed[iybin_signed]->Fill(mass, 1.);
	  }
	  else if( pdgId < 0 ) {
	    lamBarPtrs[iybin]->Fill(mass, 1.);
	    //lamBarPtrs_signed[iybin_signed]->Fill(mass, 1.);
	  }
	  if ( ( pt >= lamptarray[0] ) && ( pt < lamptarray[lamptNbins] ) ) {
	    lamptPtrs[iptbin]->Fill(mass,1.);
	    if( pdgId > 0 ) lamRegptPtrs[iptbin]->Fill(mass, 1.);
	    else if( pdgId < 0 ) lamBarptPtrs[iptbin]->Fill(mass, 1.);
	  }
	  if ( ( p >= lamparray[0] ) && ( p < lamparray[lampNbins] ) ) {
	    lampPtrs[ipbin]->Fill(mass,1.);
	  }
	  if ( ( lifetime >= lamlifearray[0] ) && ( lifetime < lamlifearray[lamlifeNbins] ) ) {
	    lamlifePtrs[ilifebin]->Fill(mass,1.);
	  }
	}
      }
   }
}

#endif
