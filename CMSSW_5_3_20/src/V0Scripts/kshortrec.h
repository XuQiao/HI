//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 20 15:28:48 2010 by ROOT version 5.27/02
// from TTree ntuple/vee ntuple
// found on file: ../mc_500kEvents.root
//////////////////////////////////////////////////////////

#ifndef kshortrec_h
#define kshortrec_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include "constants.h"
using namespace std;

class kshortrec {
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
   vector<int>     *priVtxWasRefit;
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
   Int_t           processType;
   vector<float>   *priRefitVtxNdof;
   vector<float>   *priRefitVtxX;
   vector<float>   *priRefitVtxY;
   vector<float>   *priRefitVtxZ;
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
   TBranch        *b_processType;   //!
   TBranch        *b_priRefitVtxNdof;   //!
   TBranch        *b_priRefitVtxX;   //!
   TBranch        *b_priRefitVtxY;   //!
   TBranch        *b_priRefitVtxZ;   //!
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

   kshortrec(TTree *tree);
   /*
   virtual ~kshortrec();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   */
    ~kshortrec();
    Int_t    Cut(Long64_t entry);
    Int_t    GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    void     Init(TTree *tree);
    //void     Loop(int, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**);
    void     Loop(int, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**, TH1F**);
    Bool_t   Notify();
    void     Show(Long64_t entry = -1);
};

#endif

