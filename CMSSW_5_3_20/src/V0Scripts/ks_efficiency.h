//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 25 07:00:14 2010 by ROOT version 5.27/02
// from TTree kseff/ks efficiency
// found on file: ../mc_7TeV_pythia8.root
//////////////////////////////////////////////////////////

#ifndef ks_efficiency_h
#define ks_efficiency_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include "constants.h"
using namespace std;

class ks_efficiency {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           kstrigHF;
   Int_t           kstrigTech40_41;
   Int_t           kstrigScraping;
   Int_t           kstrigTech36_39;
   Int_t           kstrigTech34;
   Int_t           kstrigAlgo124;
   Int_t           trigHLTminBias;
   Int_t           nLooseTracks;
   Int_t           nHighPurTracks;
   Int_t           nLooseTracks_pt_500;
   Int_t           nHighPurTracks_pt_500;
   Int_t           runNum;
   Int_t           evtNum;
   Int_t           processType;
   Float_t         priVtxX;
   Float_t         priVtxY;
   Float_t         priVtxZ;
   Float_t         priVtxChi2;
   Float_t         priVtxNdof;
   Int_t           priVtxNtracks;
   Int_t           priVtxIsFake;
   vector<float>   *ksgenEta;
   vector<float>   *ksgenpT;
   vector<float>   *ksgenp;
   vector<float>   *ksgenPhi;
   vector<float>   *ksgenRho;
   vector<float>   *ksgenR;
   vector<float>   *ksgenMass;
   vector<float>   *ksgenCtau;
   vector<float>   *ksgendauMaxEta;
   vector<float>   *ksrecoEta;
   vector<float>   *ksrecopT;
   vector<float>   *ksrecop;
   vector<float>   *ksrecoPhi;
   vector<float>   *ksrecoRho;
   vector<float>   *ksrecoR;
   vector<float>   *ksrecoMass;
   vector<float>   *ksrecoCtau;
   vector<float>   *ksrecoVtxX;
   vector<float>   *ksrecoVtxY;
   vector<float>   *ksrecoVtxZ;
   vector<float>   *ksrecoVtxSig;
   vector<float>   *ksrecoVtxChi2;
   vector<float>   *ksrecoVtxNormChi2;
   vector<float>   *ksrecopX;
   vector<float>   *ksrecopY;
   vector<float>   *ksrecopZ;
   vector<float>   *ksrecoOtherCandMass;
   vector<float>   *ksrecoPosDauChi2;
   vector<float>   *ksrecoPosDauNormChi2;
   vector<float>   *ksrecoPosDauD0;
   vector<float>   *ksrecoNegDauChi2;
   vector<float>   *ksrecoNegDauNormChi2;
   vector<float>   *ksrecoNegDauD0;
   vector<float>   *ks3dIpWrtPrimary;
   vector<float>   *ks3dIpWrtPrimaryError;
   vector<float>   *ks3dIpWrtPrimarySig;
   vector<float>   *ksposDau3DIp;
   vector<float>   *ksposDau3DIpError;
   vector<float>   *ksposDau3DIpSig;
   vector<float>   *ksnegDau3DIp;
   vector<float>   *ksnegDau3DIpError;
   vector<float>   *ksnegDau3DIpSig;
   vector<int>     *ksposDauNhits;
   vector<int>     *ksnegDauNhits;
   vector<float>   *ksposDeltaR;
   vector<float>   *ksposDeltaEta;
   vector<float>   *ksposDeltaPhi;
   vector<float>   *kspospT;
   vector<float>   *ksposp;
   vector<float>   *ksnegDeltaR;
   vector<float>   *ksnegDeltaEta;
   vector<float>   *ksnegDeltaPhi;
   vector<float>   *ksnegpT;
   vector<float>   *ksnegp;
   vector<float>   *ksdeltaL;
   vector<float>   *ksgenProdVtxX;
   vector<float>   *ksgenProdVtxY;
   vector<float>   *ksgenProdVtxZ;
   vector<float>   *ksgenDecVtxX;
   vector<float>   *ksgenDecVtxY;
   vector<float>   *ksgenDecVtxZ;
   vector<int>     *ksv0MatchStatus;
   vector<int>     *ksgenPDG;
   vector<int>     *ksgenPosDauPDG;
   vector<int>     *ksgenNegDauPDG;

   // List of branches
   TBranch        *b_kstrigHF;   //!
   TBranch        *b_kstrigTech40_41;   //!
   TBranch        *b_kstrigScraping;   //!
   TBranch        *b_kstrigTech36_39;   //!
   TBranch        *b_kstrigTech34;   //!
   TBranch        *b_kstrigAlgo124;   //!
   TBranch        *b_trigHLTminBias;   //!
   TBranch        *b_nLooseTracks;   //!
   TBranch        *b_nHighPurTracks;   //!
   TBranch        *b_nLooseTracks_pt_500;   //!
   TBranch        *b_nHighPurTracks_pt_500;   //!
   TBranch        *b_runNum;   //!
   TBranch        *b_evtNum;   //!
   TBranch        *b_processType;   //!
   TBranch        *b_priVtxX;   //!
   TBranch        *b_priVtxY;   //!
   TBranch        *b_priVtxZ;   //!
   TBranch        *b_priVtxChi2;   //!
   TBranch        *b_priVtxNdof;   //!
   TBranch        *b_priVtxNtracks;   //!
   TBranch        *b_priVtxIsFake;   //!
   TBranch        *b_ksgenEta;   //!
   TBranch        *b_ksgenpT;   //!
   TBranch        *b_ksgenp;   //!
   TBranch        *b_ksgenPhi;   //!
   TBranch        *b_ksgenRho;   //!
   TBranch        *b_ksgenR;   //!
   TBranch        *b_ksgenMass;   //!
   TBranch        *b_ksgenCtau;   //!
   TBranch        *b_ksgendauMaxEta;   //!
   TBranch        *b_ksrecoEta;   //!
   TBranch        *b_ksrecopT;   //!
   TBranch        *b_ksrecop;   //!
   TBranch        *b_ksrecoPhi;   //!
   TBranch        *b_ksrecoRho;   //!
   TBranch        *b_ksrecoR;   //!
   TBranch        *b_ksrecoMass;   //!
   TBranch        *b_ksrecoCtau;   //!
   TBranch        *b_ksrecoVtxX;   //!
   TBranch        *b_ksrecoVtxY;   //!
   TBranch        *b_ksrecoVtxZ;   //!
   TBranch        *b_ksrecoVtxSig;   //!
   TBranch        *b_ksrecoVtxChi2;   //!
   TBranch        *b_ksrecoVtxNormChi2;   //!
   TBranch        *b_ksrecopX;   //!
   TBranch        *b_ksrecopY;   //!
   TBranch        *b_ksrecopZ;   //!
   TBranch        *b_ksrecoOtherCandMass;   //!
   TBranch        *b_ksrecoPosDauChi2;   //!
   TBranch        *b_ksrecoPosDauNormChi2;   //!
   TBranch        *b_ksrecoPosDauD0;   //!
   TBranch        *b_ksrecoNegDauChi2;   //!
   TBranch        *b_ksrecoNegDauNormChi2;   //!
   TBranch        *b_ksrecoNegDauD0;   //!
   TBranch        *b_ks3dIpWrtPrimary;   //!
   TBranch        *b_ks3dIpWrtPrimaryError;   //!
   TBranch        *b_ks3dIpWrtPrimarySig;   //!
   TBranch        *b_ksposDau3DIp;   //!
   TBranch        *b_ksposDau3DIpError;   //!
   TBranch        *b_ksposDau3DIpSig;   //!
   TBranch        *b_ksnegDau3DIp;   //!
   TBranch        *b_ksnegDau3DIpError;   //!
   TBranch        *b_ksnegDau3DIpSig;   //!
   TBranch        *b_ksposDauNhits;   //!
   TBranch        *b_ksnegDauNhits;   //!
   TBranch        *b_ksposDeltaR;   //!
   TBranch        *b_ksposDeltaEta;   //!
   TBranch        *b_ksposDeltaPhi;   //!
   TBranch        *b_kspospT;   //!
   TBranch        *b_ksposp;   //!
   TBranch        *b_ksnegDeltaR;   //!
   TBranch        *b_ksnegDeltaEta;   //!
   TBranch        *b_ksnegDeltaPhi;   //!
   TBranch        *b_ksnegpT;   //!
   TBranch        *b_ksnegp;   //!
   TBranch        *b_ksdeltaL;   //!
   TBranch        *b_ksgenProdVtxX;   //!
   TBranch        *b_ksgenProdVtxY;   //!
   TBranch        *b_ksgenProdVtxZ;   //!
   TBranch        *b_ksgenDecVtxX;   //!
   TBranch        *b_ksgenDecVtxY;   //!
   TBranch        *b_ksgenDecVtxZ;   //!
   TBranch        *b_ksv0MatchStatus;   //!
   TBranch        *b_ksgenPDG;   //!
   TBranch        *b_ksgenPosDauPDG;   //!
   TBranch        *b_ksgenNegDauPDG;   //!

   ks_efficiency(TTree *tree);
   virtual ~ks_efficiency();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ks_efficiency_cxx
ks_efficiency::ks_efficiency(TTree *tree)
{
   Init(tree);
}

ks_efficiency::~ks_efficiency()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ks_efficiency::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ks_efficiency::LoadTree(Long64_t entry)
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

void ks_efficiency::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ksgenEta = 0;
   ksgenpT = 0;
   ksgenp = 0;
   ksgenPhi = 0;
   ksgenRho = 0;
   ksgenR = 0;
   ksgenMass = 0;
   ksgenCtau = 0;
   ksgendauMaxEta = 0;
   ksrecoEta = 0;
   ksrecopT = 0;
   ksrecop = 0;
   ksrecoPhi = 0;
   ksrecoRho = 0;
   ksrecoR = 0;
   ksrecoMass = 0;
   ksrecoCtau = 0;
   ksrecoVtxX = 0;
   ksrecoVtxY = 0;
   ksrecoVtxZ = 0;
   ksrecoVtxSig = 0;
   ksrecoVtxChi2 = 0;
   ksrecoVtxNormChi2 = 0;
   ksrecopX = 0;
   ksrecopY = 0;
   ksrecopZ = 0;
   ksrecoOtherCandMass = 0;
   ksrecoPosDauChi2 = 0;
   ksrecoPosDauNormChi2 = 0;
   ksrecoPosDauD0 = 0;
   ksrecoNegDauChi2 = 0;
   ksrecoNegDauNormChi2 = 0;
   ksrecoNegDauD0 = 0;
   ks3dIpWrtPrimary = 0;
   ks3dIpWrtPrimaryError = 0;
   ks3dIpWrtPrimarySig = 0;
   ksposDau3DIp = 0;
   ksposDau3DIpError = 0;
   ksposDau3DIpSig = 0;
   ksnegDau3DIp = 0;
   ksnegDau3DIpError = 0;
   ksnegDau3DIpSig = 0;
   ksposDauNhits = 0;
   ksnegDauNhits = 0;
   ksposDeltaR = 0;
   ksposDeltaEta = 0;
   ksposDeltaPhi = 0;
   kspospT = 0;
   ksposp = 0;
   ksnegDeltaR = 0;
   ksnegDeltaEta = 0;
   ksnegDeltaPhi = 0;
   ksnegpT = 0;
   ksnegp = 0;
   ksdeltaL = 0;
   ksgenProdVtxX = 0;
   ksgenProdVtxY = 0;
   ksgenProdVtxZ = 0;
   ksgenDecVtxX = 0;
   ksgenDecVtxY = 0;
   ksgenDecVtxZ = 0;
   ksv0MatchStatus = 0;
   ksgenPDG = 0;
   ksgenPosDauPDG = 0;
   ksgenNegDauPDG = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("kstrigHF", &kstrigHF, &b_kstrigHF);
   fChain->SetBranchAddress("kstrigTech40_41", &kstrigTech40_41, &b_kstrigTech40_41);
   fChain->SetBranchAddress("kstrigScraping", &kstrigScraping, &b_kstrigScraping);
   fChain->SetBranchAddress("kstrigTech36_39", &kstrigTech36_39, &b_kstrigTech36_39);
   fChain->SetBranchAddress("kstrigTech34", &kstrigTech34, &b_kstrigTech34);
   fChain->SetBranchAddress("kstrigAlgo124", &kstrigAlgo124, &b_kstrigAlgo124);
   fChain->SetBranchAddress("trigHLTminBias", &trigHLTminBias, &b_trigHLTminBias);
   fChain->SetBranchAddress("nLooseTracks", &nLooseTracks, &b_nLooseTracks);
   fChain->SetBranchAddress("nHighPurTracks", &nHighPurTracks, &b_nHighPurTracks);
   fChain->SetBranchAddress("nLooseTracks_pt_500", &nLooseTracks_pt_500, &b_nLooseTracks_pt_500);
   fChain->SetBranchAddress("nHighPurTracks_pt_500", &nHighPurTracks_pt_500, &b_nHighPurTracks_pt_500);
   fChain->SetBranchAddress("runNum", &runNum, &b_runNum);
   fChain->SetBranchAddress("evtNum", &evtNum, &b_evtNum);
   fChain->SetBranchAddress("processType", &processType, &b_processType);
   fChain->SetBranchAddress("priVtxX", &priVtxX, &b_priVtxX);
   fChain->SetBranchAddress("priVtxY", &priVtxY, &b_priVtxY);
   fChain->SetBranchAddress("priVtxZ", &priVtxZ, &b_priVtxZ);
   fChain->SetBranchAddress("priVtxChi2", &priVtxChi2, &b_priVtxChi2);
   fChain->SetBranchAddress("priVtxNdof", &priVtxNdof, &b_priVtxNdof);
   fChain->SetBranchAddress("priVtxNtracks", &priVtxNtracks, &b_priVtxNtracks);
   fChain->SetBranchAddress("priVtxIsFake", &priVtxIsFake, &b_priVtxIsFake);
   fChain->SetBranchAddress("ksgenEta", &ksgenEta, &b_ksgenEta);
   fChain->SetBranchAddress("ksgenpT", &ksgenpT, &b_ksgenpT);
   fChain->SetBranchAddress("ksgenp", &ksgenp, &b_ksgenp);
   fChain->SetBranchAddress("ksgenPhi", &ksgenPhi, &b_ksgenPhi);
   fChain->SetBranchAddress("ksgenRho", &ksgenRho, &b_ksgenRho);
   fChain->SetBranchAddress("ksgenR", &ksgenR, &b_ksgenR);
   fChain->SetBranchAddress("ksgenMass", &ksgenMass, &b_ksgenMass);
   fChain->SetBranchAddress("ksgenCtau", &ksgenCtau, &b_ksgenCtau);
   fChain->SetBranchAddress("ksgendauMaxEta", &ksgendauMaxEta, &b_ksgendauMaxEta);
   fChain->SetBranchAddress("ksrecoEta", &ksrecoEta, &b_ksrecoEta);
   fChain->SetBranchAddress("ksrecopT", &ksrecopT, &b_ksrecopT);
   fChain->SetBranchAddress("ksrecop", &ksrecop, &b_ksrecop);
   fChain->SetBranchAddress("ksrecoPhi", &ksrecoPhi, &b_ksrecoPhi);
   fChain->SetBranchAddress("ksrecoRho", &ksrecoRho, &b_ksrecoRho);
   fChain->SetBranchAddress("ksrecoR", &ksrecoR, &b_ksrecoR);
   fChain->SetBranchAddress("ksrecoMass", &ksrecoMass, &b_ksrecoMass);
   fChain->SetBranchAddress("ksrecoCtau", &ksrecoCtau, &b_ksrecoCtau);
   fChain->SetBranchAddress("ksrecoVtxX", &ksrecoVtxX, &b_ksrecoVtxX);
   fChain->SetBranchAddress("ksrecoVtxY", &ksrecoVtxY, &b_ksrecoVtxY);
   fChain->SetBranchAddress("ksrecoVtxZ", &ksrecoVtxZ, &b_ksrecoVtxZ);
   fChain->SetBranchAddress("ksrecoVtxSig", &ksrecoVtxSig, &b_ksrecoVtxSig);
   fChain->SetBranchAddress("ksrecoVtxChi2", &ksrecoVtxChi2, &b_ksrecoVtxChi2);
   fChain->SetBranchAddress("ksrecoVtxNormChi2", &ksrecoVtxNormChi2, &b_ksrecoVtxNormChi2);
   fChain->SetBranchAddress("ksrecopX", &ksrecopX, &b_ksrecopX);
   fChain->SetBranchAddress("ksrecopY", &ksrecopY, &b_ksrecopY);
   fChain->SetBranchAddress("ksrecopZ", &ksrecopZ, &b_ksrecopZ);
   fChain->SetBranchAddress("ksrecoOtherCandMass", &ksrecoOtherCandMass, &b_ksrecoOtherCandMass);
   fChain->SetBranchAddress("ksrecoPosDauChi2", &ksrecoPosDauChi2, &b_ksrecoPosDauChi2);
   fChain->SetBranchAddress("ksrecoPosDauNormChi2", &ksrecoPosDauNormChi2, &b_ksrecoPosDauNormChi2);
   fChain->SetBranchAddress("ksrecoPosDauD0", &ksrecoPosDauD0, &b_ksrecoPosDauD0);
   fChain->SetBranchAddress("ksrecoNegDauChi2", &ksrecoNegDauChi2, &b_ksrecoNegDauChi2);
   fChain->SetBranchAddress("ksrecoNegDauNormChi2", &ksrecoNegDauNormChi2, &b_ksrecoNegDauNormChi2);
   fChain->SetBranchAddress("ksrecoNegDauD0", &ksrecoNegDauD0, &b_ksrecoNegDauD0);
   fChain->SetBranchAddress("ks3dIpWrtPrimary", &ks3dIpWrtPrimary, &b_ks3dIpWrtPrimary);
   fChain->SetBranchAddress("ks3dIpWrtPrimaryError", &ks3dIpWrtPrimaryError, &b_ks3dIpWrtPrimaryError);
   fChain->SetBranchAddress("ks3dIpWrtPrimarySig", &ks3dIpWrtPrimarySig, &b_ks3dIpWrtPrimarySig);
   fChain->SetBranchAddress("ksposDau3DIp", &ksposDau3DIp, &b_ksposDau3DIp);
   fChain->SetBranchAddress("ksposDau3DIpError", &ksposDau3DIpError, &b_ksposDau3DIpError);
   fChain->SetBranchAddress("ksposDau3DIpSig", &ksposDau3DIpSig, &b_ksposDau3DIpSig);
   fChain->SetBranchAddress("ksnegDau3DIp", &ksnegDau3DIp, &b_ksnegDau3DIp);
   fChain->SetBranchAddress("ksnegDau3DIpError", &ksnegDau3DIpError, &b_ksnegDau3DIpError);
   fChain->SetBranchAddress("ksnegDau3DIpSig", &ksnegDau3DIpSig, &b_ksnegDau3DIpSig);
   fChain->SetBranchAddress("ksposDauNhits", &ksposDauNhits, &b_ksposDauNhits);
   fChain->SetBranchAddress("ksnegDauNhits", &ksnegDauNhits, &b_ksnegDauNhits);
   fChain->SetBranchAddress("ksposDeltaR", &ksposDeltaR, &b_ksposDeltaR);
   fChain->SetBranchAddress("ksposDeltaEta", &ksposDeltaEta, &b_ksposDeltaEta);
   fChain->SetBranchAddress("ksposDeltaPhi", &ksposDeltaPhi, &b_ksposDeltaPhi);
   fChain->SetBranchAddress("kspospT", &kspospT, &b_kspospT);
   fChain->SetBranchAddress("ksposp", &ksposp, &b_ksposp);
   fChain->SetBranchAddress("ksnegDeltaR", &ksnegDeltaR, &b_ksnegDeltaR);
   fChain->SetBranchAddress("ksnegDeltaEta", &ksnegDeltaEta, &b_ksnegDeltaEta);
   fChain->SetBranchAddress("ksnegDeltaPhi", &ksnegDeltaPhi, &b_ksnegDeltaPhi);
   fChain->SetBranchAddress("ksnegpT", &ksnegpT, &b_ksnegpT);
   fChain->SetBranchAddress("ksnegp", &ksnegp, &b_ksnegp);
   fChain->SetBranchAddress("ksdeltaL", &ksdeltaL, &b_ksdeltaL);
   fChain->SetBranchAddress("ksgenProdVtxX", &ksgenProdVtxX, &b_ksgenProdVtxX);
   fChain->SetBranchAddress("ksgenProdVtxY", &ksgenProdVtxY, &b_ksgenProdVtxY);
   fChain->SetBranchAddress("ksgenProdVtxZ", &ksgenProdVtxZ, &b_ksgenProdVtxZ);
   fChain->SetBranchAddress("ksgenDecVtxX", &ksgenDecVtxX, &b_ksgenDecVtxX);
   fChain->SetBranchAddress("ksgenDecVtxY", &ksgenDecVtxY, &b_ksgenDecVtxY);
   fChain->SetBranchAddress("ksgenDecVtxZ", &ksgenDecVtxZ, &b_ksgenDecVtxZ);
   fChain->SetBranchAddress("ksv0MatchStatus", &ksv0MatchStatus, &b_ksv0MatchStatus);
   fChain->SetBranchAddress("ksgenPDG", &ksgenPDG, &b_ksgenPDG);
   fChain->SetBranchAddress("ksgenPosDauPDG", &ksgenPosDauPDG, &b_ksgenPosDauPDG);
   fChain->SetBranchAddress("ksgenNegDauPDG", &ksgenNegDauPDG, &b_ksgenNegDauPDG);
   Notify();
}

Bool_t ks_efficiency::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ks_efficiency::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ks_efficiency::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ks_efficiency_cxx
