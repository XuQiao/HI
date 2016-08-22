//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 20 15:28:48 2010 by ROOT version 5.27/02
// from TTree ntuple/vee ntuple
// found on file: ../mc_500kEvents.root
//////////////////////////////////////////////////////////

#ifndef kshortrec_masspeak_h
#define kshortrec_masspeak_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class kshortrec_masspeak {
public :
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<float>   *v0VtxX;
   vector<float>   *v0VtxY;
   vector<float>   *v0VtxZ;
   vector<float>   *v0VtxR;
   vector<float>   *v0VtxNormChi2;
   vector<float>   *v0pX;
   vector<float>   *v0pY;
   vector<float>   *v0pZ;
   vector<float>   *v0pT;
   vector<float>   *v0Eta;
   vector<float>   *v0Phi;
   vector<float>   *v0CandMass;
   vector<float>   *v0OtherCandMass;
   vector<float>   *v0Lifetime;
   vector<float>   *v03DLifetime;
   vector<float>   *v03DLifetimeError;
   vector<float>   *v03DLifetimeScaled;
   vector<float>   *v03DLifetimeScaledError;
   vector<float>   *v03DLifetimeMPV;
   vector<float>   *v03DLifetimeMPVError;
   vector<int>     *v0PDG;
   vector<float>   *v0VtxSig;
   vector<float>   *v0CosThetaPAndLineOfFlight;
   vector<float>   *v0PCAPrimary;
   vector<float>   *v0PCAPrimaryError;
   vector<float>   *v0PCABeamSpot;
   vector<float>   *v0PCABeamSpotError;
   vector<float>   *v3dIpWrtPrimary;
   vector<float>   *v3dIpWrtPrimarySig;
   vector<float>   *thePosDau3DIp;
   vector<float>   *thePosDau3DIpSig;
   vector<float>   *theNegDau3DIp;
   vector<float>   *theNegDau3DIpSig;
   vector<int>     *posDauNhits;
   vector<int>     *negDauNhits;
   Float_t         priVtxChi2;
   Float_t         privtxNdof;
   Float_t         priVtxX;
   Float_t         priVtxY;
   Float_t         priVtxZ;
   Int_t           priVtxNtracks;
   Int_t           priVtxIsFake;
   Int_t           priVtxWasRefit;
   Int_t           nLooseTracks;
   Int_t           nHighPurTracks;
   Int_t           nLooseTracks_pt_500;
   Int_t           nHighPurTracks_pt_500;
   Int_t           nV0s;
   Int_t           trigHF;
   Int_t           trigTech40_41;
   Int_t           trigTech34;
   Int_t           trigTech36_39;
   Int_t           trigScraping;
   Int_t           trigAlgo124;
   Int_t           trigHLTminBias;
   Int_t           runNum;
   Int_t           evtNum;
   vector<float>   *v0VtxPriDist;
   vector<float>   *v0VtxPriDistErr;
   vector<float>   *v0VtxBSTransDist;
   vector<float>   *v0VtxBSTransDistErr;
   vector<float>   *v0PosDaupX;
   vector<float>   *v0PosDaupY;
   vector<float>   *v0PosDaupZ;
   vector<float>   *v0PosDauAtVtxpX;
   vector<float>   *v0PosDauAtVtxpY;
   vector<float>   *v0PosDauAtVtxpZ;
   vector<float>   *v0PosDauEta;
   vector<float>   *v0PosDauNValidHits;
   vector<float>   *v0PosDauNPixelHits;
   vector<float>   *v0PosDauNStripHits;
   vector<float>   *v0PosDauChi2;
   vector<float>   *v0PosDauNormChi2;
   vector<float>   *v0PosDauD0;
   vector<float>   *v0DaupT_lowestpT;
   vector<float>   *v0NegDaupX;
   vector<float>   *v0NegDaupY;
   vector<float>   *v0NegDaupZ;
   vector<float>   *v0NegDauAtVtxpX;
   vector<float>   *v0NegDauAtVtxpY;
   vector<float>   *v0NegDauAtVtxpZ;
   vector<float>   *v0NegDauEta;
   vector<float>   *v0NegDauNValidHits;
   vector<float>   *v0NegDauNPixelHits;
   vector<float>   *v0NegDauNStripHits;
   vector<float>   *v0NegDauChi2;
   vector<float>   *v0NegDauNormChi2;
   vector<float>   *v0NegDauD0;
   vector<float>   *v0DauEta_lowestpT;
   vector<float>   *v0DauPhi_lowestpT;

   // List of branches
   TBranch        *b_v0VtxX;   //!
   TBranch        *b_v0VtxY;   //!
   TBranch        *b_v0VtxZ;   //!
   TBranch        *b_v0VtxR;   //!
   TBranch        *b_v0VtxNormChi2;   //!
   TBranch        *b_v0pX;   //!
   TBranch        *b_v0pY;   //!
   TBranch        *b_v0pZ;   //!
   TBranch        *b_v0pT;   //!
   TBranch        *b_v0Eta;   //!
   TBranch        *b_v0Phi;   //!
   TBranch        *b_v0CandMass;   //!
   TBranch        *b_v0OtherCandMass;   //!
   TBranch        *b_v0Lifetime;   //!
   TBranch        *b_v03DLifetime;   //!
   TBranch        *b_v03DLifetimeError;   //!
   TBranch        *b_v03DLifetimeScaled;   //!
   TBranch        *b_v03DLifetimeScaledError;   //!
   TBranch        *b_v03DLifetimeMPV;   //!
   TBranch        *b_v03DLifetimeMPVError;   //!
   TBranch        *b_v0PDG;   //!
   TBranch        *b_v0VtxSig;   //!
   TBranch        *b_v0CosThetaPAndLineOfFlight;   //!
   TBranch        *b_v0PCAPrimary;   //!
   TBranch        *b_v0PCAPrimaryError;   //!
   TBranch        *b_v0PCABeamSpot;   //!
   TBranch        *b_v0PCABeamSpotError;   //!
   TBranch        *b_v3dIpWrtPrimary;   //!
   TBranch        *b_v3dIpWrtPrimarySig;   //!
   TBranch        *b_thePosDau3DIp;   //!
   TBranch        *b_thePosDau3DIpSig;   //!
   TBranch        *b_theNegDau3DIp;   //!
   TBranch        *b_theNegDau3DIpSig;   //!
   TBranch        *b_posDauNhits;   //!
   TBranch        *b_negDauNhits;   //!
   TBranch        *b_priVtxChi2;   //!
   TBranch        *b_priVtxNdof;   //!
   TBranch        *b_priVtxX;   //!
   TBranch        *b_priVtxY;   //!
   TBranch        *b_priVtxZ;   //!
   TBranch        *b_priVtxNtracks;   //!
   TBranch        *b_priVtxIsFake;   //!
   TBranch        *b_priVtxWasRefit;   //!
   TBranch        *b_nLooseTracks;   //!
   TBranch        *b_nHighPurTracks;   //!
   TBranch        *b_nLooseTracks_pt_500;   //!
   TBranch        *b_nHighPurTracks_pt_500;   //!
   TBranch        *b_nV0s;   //!
   TBranch        *b_trigHF;   //!
   TBranch        *b_trigTech40_41;   //!
   TBranch        *b_trigTech34;   //!
   TBranch        *b_trigTech36_39;   //!
   TBranch        *b_trigScraping;   //!
   TBranch        *b_trigAlgo124;   //!
   TBranch        *b_trigHLTminBias;   //!
   TBranch        *b_runNum;   //!
   TBranch        *b_evtNum;   //!
   TBranch        *b_v0VtxPriDist;   //!
   TBranch        *b_v0VtxPriDistErr;   //!
   TBranch        *b_v0VtxBSTransDist;   //!
   TBranch        *b_v0VtxBSTransDistErr;   //!
   TBranch        *b_v0PosDaupX;   //!
   TBranch        *b_v0PosDaupY;   //!
   TBranch        *b_v0PosDaupZ;   //!
   TBranch        *b_v0PosDauAtVtxpX;   //!
   TBranch        *b_v0PosDauAtVtxpY;   //!
   TBranch        *b_v0PosDauAtVtxpZ;   //!
   TBranch        *b_v0PosDauEta;   //!
   TBranch        *b_v0PosDauNValidHits;   //!
   TBranch        *b_v0PosDauNPixelHits;   //!
   TBranch        *b_v0PosDauNStripHits;   //!
   TBranch        *b_v0PosDauChi2;   //!
   TBranch        *b_v0PosDauNormChi2;   //!
   TBranch        *b_v0PosDauD0;   //!
   TBranch        *b_v0DaupT_lowestpT;   //!
   TBranch        *b_v0NegDaupX;   //!
   TBranch        *b_v0NegDaupY;   //!
   TBranch        *b_v0NegDaupZ;   //!
   TBranch        *b_v0NegDauAtVtxpX;   //!
   TBranch        *b_v0NegDauAtVtxpY;   //!
   TBranch        *b_v0NegDauAtVtxpZ;   //!
   TBranch        *b_v0NegDauEta;   //!
   TBranch        *b_v0NegDauNValidHits;   //!
   TBranch        *b_v0NegDauNPixelHits;   //!
   TBranch        *b_v0NegDauNStripHits;   //!
   TBranch        *b_v0NegDauChi2;   //!
   TBranch        *b_v0NegDauNormChi2;   //!
   TBranch        *b_v0NegDauD0;   //!
   TBranch        *b_v0DauEta_lowestpT;   //!
   TBranch        *b_v0DauPhi_lowestpT;   //!

   kshortrec_masspeak(TTree *tree=0);

   virtual ~kshortrec_masspeak();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef kshortrec_masspeak_cxx
kshortrec_masspeak::kshortrec_masspeak(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
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

   }
   Init(tree);
}

kshortrec_masspeak::~kshortrec_masspeak()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t kshortrec_masspeak::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t kshortrec_masspeak::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void kshortrec_masspeak::Init(TTree *tree)
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

Bool_t kshortrec_masspeak::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void kshortrec_masspeak::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t kshortrec_masspeak::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef kshortrec_masspeak_cxx
