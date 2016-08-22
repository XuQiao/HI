//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 25 07:00:14 2010 by ROOT version 5.27/02
// from TTree kseff/ks efficiency
// found on file: ../mc_7TeV_pythia8.root
//////////////////////////////////////////////////////////

#ifndef lam_efficiency_h
#define lam_efficiency_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <TCanvas.h>
#include "constants.h"
using namespace std;

class lam_efficiency {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           lamtrigHF;
   Int_t           lamtrigTech40_41;
   Int_t           lamtrigScraping;
   Int_t           lamtrigTech36_39;
   Int_t           lamtrigTech34;
   Int_t           lamtrigAlgo124;
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
   vector<float>   *lamgenEta;
   vector<float>   *lamgenpT;
   vector<float>   *lamgenp;
   vector<float>   *lamgenPhi;
   vector<float>   *lamgenRho;
   vector<float>   *lamgenR;
   vector<float>   *lamgenMass;
   vector<float>   *lamgenCtau;
   vector<float>   *lamgendauMaxEta;
   vector<float>   *lamrecoEta;
   vector<float>   *lamrecopT;
   vector<float>   *lamrecop;
   vector<float>   *lamrecoPhi;
   vector<float>   *lamrecoRho;
   vector<float>   *lamrecoR;
   vector<float>   *lamrecoMass;
   vector<float>   *lamrecoCtau;
   vector<float>   *lamrecoVtxX;
   vector<float>   *lamrecoVtxY;
   vector<float>   *lamrecoVtxZ;
   vector<float>   *lamrecoVtxSig;
   vector<float>   *lamrecoVtxChi2;
   vector<float>   *lamrecoVtxNormChi2;
   vector<float>   *lamrecopX;
   vector<float>   *lamrecopY;
   vector<float>   *lamrecopZ;
   vector<float>   *lamrecoOtherCandMass;
   vector<float>   *lamrecoPosDauChi2;
   vector<float>   *lamrecoPosDauNormChi2;
   vector<float>   *lamrecoPosDauD0;
   vector<float>   *lamrecoNegDauChi2;
   vector<float>   *lamrecoNegDauNormChi2;
   vector<float>   *lamrecoNegDauD0;
   vector<float>   *lam3dIpWrtPrimary;
   vector<float>   *lam3dIpWrtPrimaryError;
   vector<float>   *lam3dIpWrtPrimarySig;
   vector<float>   *lamposDau3DIp;
   vector<float>   *lamposDau3DIpError;
   vector<float>   *lamposDau3DIpSig;
   vector<float>   *lamnegDau3DIp;
   vector<float>   *lamnegDau3DIpError;
   vector<float>   *lamnegDau3DIpSig;
   vector<int>     *lamposDauNhits;
   vector<int>     *lamnegDauNhits;
   vector<float>   *lamposDeltaR;
   vector<float>   *lamposDeltaEta;
   vector<float>   *lamposDeltaPhi;
   vector<float>   *lampospT;
   vector<float>   *lamposp;
   vector<float>   *lamnegDeltaR;
   vector<float>   *lamnegDeltaEta;
   vector<float>   *lamnegDeltaPhi;
   vector<float>   *lamnegpT;
   vector<float>   *lamnegp;
   vector<float>   *lamdeltaL;
   vector<float>   *lamgenProdVtxX;
   vector<float>   *lamgenProdVtxY;
   vector<float>   *lamgenProdVtxZ;
   vector<float>   *lamgenDecVtxX;
   vector<float>   *lamgenDecVtxY;
   vector<float>   *lamgenDecVtxZ;
   vector<int>     *lamv0MatchStatus;
   vector<int>     *lamgenPDG;
   vector<int>     *lamgenPosDauPDG;
   vector<int>     *lamgenNegDauPDG;

   // List of branches
   TBranch        *b_lamtrigHF;   //!
   TBranch        *b_lamtrigTech40_41;   //!
   TBranch        *b_lamtrigScraping;   //!
   TBranch        *b_lamtrigTech36_39;   //!
   TBranch        *b_lamtrigTech34;   //!
   TBranch        *b_lamtrigAlgo124;   //!
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
   TBranch        *b_lamgenEta;   //!
   TBranch        *b_lamgenpT;   //!
   TBranch        *b_lamgenp;   //!
   TBranch        *b_lamgenPhi;   //!
   TBranch        *b_lamgenRho;   //!
   TBranch        *b_lamgenR;   //!
   TBranch        *b_lamgenMass;   //!
   TBranch        *b_lamgenCtau;   //!
   TBranch        *b_lamgendauMaxEta;   //!
   TBranch        *b_lamrecoEta;   //!
   TBranch        *b_lamrecopT;   //!
   TBranch        *b_lamrecop;   //!
   TBranch        *b_lamrecoPhi;   //!
   TBranch        *b_lamrecoRho;   //!
   TBranch        *b_lamrecoR;   //!
   TBranch        *b_lamrecoMass;   //!
   TBranch        *b_lamrecoCtau;   //!
   TBranch        *b_lamrecoVtxX;   //!
   TBranch        *b_lamrecoVtxY;   //!
   TBranch        *b_lamrecoVtxZ;   //!
   TBranch        *b_lamrecoVtxSig;   //!
   TBranch        *b_lamrecoVtxChi2;   //!
   TBranch        *b_lamrecoVtxNormChi2;   //!
   TBranch        *b_lamrecopX;   //!
   TBranch        *b_lamrecopY;   //!
   TBranch        *b_lamrecopZ;   //!
   TBranch        *b_lamrecoOtherCandMass;   //!
   TBranch        *b_lamrecoPosDauChi2;   //!
   TBranch        *b_lamrecoPosDauNormChi2;   //!
   TBranch        *b_lamrecoPosDauD0;   //!
   TBranch        *b_lamrecoNegDauChi2;   //!
   TBranch        *b_lamrecoNegDauNormChi2;   //!
   TBranch        *b_lamrecoNegDauD0;   //!
   TBranch        *b_lam3dIpWrtPrimary;   //!
   TBranch        *b_lam3dIpWrtPrimaryError;   //!
   TBranch        *b_lam3dIpWrtPrimarySig;   //!
   TBranch        *b_lamposDau3DIp;   //!
   TBranch        *b_lamposDau3DIpError;   //!
   TBranch        *b_lamposDau3DIpSig;   //!
   TBranch        *b_lamnegDau3DIp;   //!
   TBranch        *b_lamnegDau3DIpError;   //!
   TBranch        *b_lamnegDau3DIpSig;   //!
   TBranch        *b_lamposDauNhits;   //!
   TBranch        *b_lamnegDauNhits;   //!
   TBranch        *b_lamposDeltaR;   //!
   TBranch        *b_lamposDeltaEta;   //!
   TBranch        *b_lamposDeltaPhi;   //!
   TBranch        *b_lampospT;   //!
   TBranch        *b_lamposp;   //!
   TBranch        *b_lamnegDeltaR;   //!
   TBranch        *b_lamnegDeltaEta;   //!
   TBranch        *b_lamnegDeltaPhi;   //!
   TBranch        *b_lamnegpT;   //!
   TBranch        *b_lamnegp;   //!
   TBranch        *b_lamdeltaL;   //!
   TBranch        *b_lamgenProdVtxX;   //!
   TBranch        *b_lamgenProdVtxY;   //!
   TBranch        *b_lamgenProdVtxZ;   //!
   TBranch        *b_lamgenDecVtxX;   //!
   TBranch        *b_lamgenDecVtxY;   //!
   TBranch        *b_lamgenDecVtxZ;   //!
   TBranch        *b_lamv0MatchStatus;   //!
   TBranch        *b_lamgenPDG;   //!
   TBranch        *b_lamgenPosDauPDG;   //!
   TBranch        *b_lamgenNegDauPDG;   //!

   lam_efficiency(TTree *tree);
   virtual ~lam_efficiency();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, TH1F*);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef lam_efficiency_cxx
lam_efficiency::lam_efficiency(TTree *tree)
{
   Init(tree);
}

lam_efficiency::~lam_efficiency()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t lam_efficiency::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t lam_efficiency::LoadTree(Long64_t entry)
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

void lam_efficiency::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   lamgenEta = 0;
   lamgenpT = 0;
   lamgenp = 0;
   lamgenPhi = 0;
   lamgenRho = 0;
   lamgenR = 0;
   lamgenMass = 0;
   lamgenCtau = 0;
   lamgendauMaxEta = 0;
   lamrecoEta = 0;
   lamrecopT = 0;
   lamrecop = 0;
   lamrecoPhi = 0;
   lamrecoRho = 0;
   lamrecoR = 0;
   lamrecoMass = 0;
   lamrecoCtau = 0;
   lamrecoVtxX = 0;
   lamrecoVtxY = 0;
   lamrecoVtxZ = 0;
   lamrecoVtxSig = 0;
   lamrecoVtxChi2 = 0;
   lamrecoVtxNormChi2 = 0;
   lamrecopX = 0;
   lamrecopY = 0;
   lamrecopZ = 0;
   lamrecoOtherCandMass = 0;
   lamrecoPosDauChi2 = 0;
   lamrecoPosDauNormChi2 = 0;
   lamrecoPosDauD0 = 0;
   lamrecoNegDauChi2 = 0;
   lamrecoNegDauNormChi2 = 0;
   lamrecoNegDauD0 = 0;
   lam3dIpWrtPrimary = 0;
   lam3dIpWrtPrimaryError = 0;
   lam3dIpWrtPrimarySig = 0;
   lamposDau3DIp = 0;
   lamposDau3DIpError = 0;
   lamposDau3DIpSig = 0;
   lamnegDau3DIp = 0;
   lamnegDau3DIpError = 0;
   lamnegDau3DIpSig = 0;
   lamposDauNhits = 0;
   lamnegDauNhits = 0;
   lamposDeltaR = 0;
   lamposDeltaEta = 0;
   lamposDeltaPhi = 0;
   lampospT = 0;
   lamposp = 0;
   lamnegDeltaR = 0;
   lamnegDeltaEta = 0;
   lamnegDeltaPhi = 0;
   lamnegpT = 0;
   lamnegp = 0;
   lamdeltaL = 0;
   lamgenProdVtxX = 0;
   lamgenProdVtxY = 0;
   lamgenProdVtxZ = 0;
   lamgenDecVtxX = 0;
   lamgenDecVtxY = 0;
   lamgenDecVtxZ = 0;
   lamv0MatchStatus = 0;
   lamgenPDG = 0;
   lamgenPosDauPDG = 0;
   lamgenNegDauPDG = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("lamtrigHF", &lamtrigHF, &b_lamtrigHF);
   fChain->SetBranchAddress("lamtrigTech40_41", &lamtrigTech40_41, &b_lamtrigTech40_41);
   fChain->SetBranchAddress("lamtrigScraping", &lamtrigScraping, &b_lamtrigScraping);
   fChain->SetBranchAddress("lamtrigTech36_39", &lamtrigTech36_39, &b_lamtrigTech36_39);
   fChain->SetBranchAddress("lamtrigTech34", &lamtrigTech34, &b_lamtrigTech34);
   fChain->SetBranchAddress("lamtrigAlgo124", &lamtrigAlgo124, &b_lamtrigAlgo124);
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
   fChain->SetBranchAddress("lamgenEta", &lamgenEta, &b_lamgenEta);
   fChain->SetBranchAddress("lamgenpT", &lamgenpT, &b_lamgenpT);
   fChain->SetBranchAddress("lamgenp", &lamgenp, &b_lamgenp);
   fChain->SetBranchAddress("lamgenPhi", &lamgenPhi, &b_lamgenPhi);
   fChain->SetBranchAddress("lamgenRho", &lamgenRho, &b_lamgenRho);
   fChain->SetBranchAddress("lamgenR", &lamgenR, &b_lamgenR);
   fChain->SetBranchAddress("lamgenMass", &lamgenMass, &b_lamgenMass);
   fChain->SetBranchAddress("lamgenCtau", &lamgenCtau, &b_lamgenCtau);
   fChain->SetBranchAddress("lamgendauMaxEta", &lamgendauMaxEta, &b_lamgendauMaxEta);
   fChain->SetBranchAddress("lamrecoEta", &lamrecoEta, &b_lamrecoEta);
   fChain->SetBranchAddress("lamrecopT", &lamrecopT, &b_lamrecopT);
   fChain->SetBranchAddress("lamrecop", &lamrecop, &b_lamrecop);
   fChain->SetBranchAddress("lamrecoPhi", &lamrecoPhi, &b_lamrecoPhi);
   fChain->SetBranchAddress("lamrecoRho", &lamrecoRho, &b_lamrecoRho);
   fChain->SetBranchAddress("lamrecoR", &lamrecoR, &b_lamrecoR);
   fChain->SetBranchAddress("lamrecoMass", &lamrecoMass, &b_lamrecoMass);
   fChain->SetBranchAddress("lamrecoCtau", &lamrecoCtau, &b_lamrecoCtau);
   fChain->SetBranchAddress("lamrecoVtxX", &lamrecoVtxX, &b_lamrecoVtxX);
   fChain->SetBranchAddress("lamrecoVtxY", &lamrecoVtxY, &b_lamrecoVtxY);
   fChain->SetBranchAddress("lamrecoVtxZ", &lamrecoVtxZ, &b_lamrecoVtxZ);
   fChain->SetBranchAddress("lamrecoVtxSig", &lamrecoVtxSig, &b_lamrecoVtxSig);
   fChain->SetBranchAddress("lamrecoVtxChi2", &lamrecoVtxChi2, &b_lamrecoVtxChi2);
   fChain->SetBranchAddress("lamrecoVtxNormChi2", &lamrecoVtxNormChi2, &b_lamrecoVtxNormChi2);
   fChain->SetBranchAddress("lamrecopX", &lamrecopX, &b_lamrecopX);
   fChain->SetBranchAddress("lamrecopY", &lamrecopY, &b_lamrecopY);
   fChain->SetBranchAddress("lamrecopZ", &lamrecopZ, &b_lamrecopZ);
   fChain->SetBranchAddress("lamrecoOtherCandMass", &lamrecoOtherCandMass, &b_lamrecoOtherCandMass);
   fChain->SetBranchAddress("lamrecoPosDauChi2", &lamrecoPosDauChi2, &b_lamrecoPosDauChi2);
   fChain->SetBranchAddress("lamrecoPosDauNormChi2", &lamrecoPosDauNormChi2, &b_lamrecoPosDauNormChi2);
   fChain->SetBranchAddress("lamrecoPosDauD0", &lamrecoPosDauD0, &b_lamrecoPosDauD0);
   fChain->SetBranchAddress("lamrecoNegDauChi2", &lamrecoNegDauChi2, &b_lamrecoNegDauChi2);
   fChain->SetBranchAddress("lamrecoNegDauNormChi2", &lamrecoNegDauNormChi2, &b_lamrecoNegDauNormChi2);
   fChain->SetBranchAddress("lamrecoNegDauD0", &lamrecoNegDauD0, &b_lamrecoNegDauD0);
   fChain->SetBranchAddress("lam3dIpWrtPrimary", &lam3dIpWrtPrimary, &b_lam3dIpWrtPrimary);
   fChain->SetBranchAddress("lam3dIpWrtPrimaryError", &lam3dIpWrtPrimaryError, &b_lam3dIpWrtPrimaryError);
   fChain->SetBranchAddress("lam3dIpWrtPrimarySig", &lam3dIpWrtPrimarySig, &b_lam3dIpWrtPrimarySig);
   fChain->SetBranchAddress("lamposDau3DIp", &lamposDau3DIp, &b_lamposDau3DIp);
   fChain->SetBranchAddress("lamposDau3DIpError", &lamposDau3DIpError, &b_lamposDau3DIpError);
   fChain->SetBranchAddress("lamposDau3DIpSig", &lamposDau3DIpSig, &b_lamposDau3DIpSig);
   fChain->SetBranchAddress("lamnegDau3DIp", &lamnegDau3DIp, &b_lamnegDau3DIp);
   fChain->SetBranchAddress("lamnegDau3DIpError", &lamnegDau3DIpError, &b_lamnegDau3DIpError);
   fChain->SetBranchAddress("lamnegDau3DIpSig", &lamnegDau3DIpSig, &b_lamnegDau3DIpSig);
   fChain->SetBranchAddress("lamposDauNhits", &lamposDauNhits, &b_lamposDauNhits);
   fChain->SetBranchAddress("lamnegDauNhits", &lamnegDauNhits, &b_lamnegDauNhits);
   fChain->SetBranchAddress("lamposDeltaR", &lamposDeltaR, &b_lamposDeltaR);
   fChain->SetBranchAddress("lamposDeltaEta", &lamposDeltaEta, &b_lamposDeltaEta);
   fChain->SetBranchAddress("lamposDeltaPhi", &lamposDeltaPhi, &b_lamposDeltaPhi);
   fChain->SetBranchAddress("lampospT", &lampospT, &b_lampospT);
   fChain->SetBranchAddress("lamposp", &lamposp, &b_lamposp);
   fChain->SetBranchAddress("lamnegDeltaR", &lamnegDeltaR, &b_lamnegDeltaR);
   fChain->SetBranchAddress("lamnegDeltaEta", &lamnegDeltaEta, &b_lamnegDeltaEta);
   fChain->SetBranchAddress("lamnegDeltaPhi", &lamnegDeltaPhi, &b_lamnegDeltaPhi);
   fChain->SetBranchAddress("lamnegpT", &lamnegpT, &b_lamnegpT);
   fChain->SetBranchAddress("lamnegp", &lamnegp, &b_lamnegp);
   fChain->SetBranchAddress("lamdeltaL", &lamdeltaL, &b_lamdeltaL);
   fChain->SetBranchAddress("lamgenProdVtxX", &lamgenProdVtxX, &b_lamgenProdVtxX);
   fChain->SetBranchAddress("lamgenProdVtxY", &lamgenProdVtxY, &b_lamgenProdVtxY);
   fChain->SetBranchAddress("lamgenProdVtxZ", &lamgenProdVtxZ, &b_lamgenProdVtxZ);
   fChain->SetBranchAddress("lamgenDecVtxX", &lamgenDecVtxX, &b_lamgenDecVtxX);
   fChain->SetBranchAddress("lamgenDecVtxY", &lamgenDecVtxY, &b_lamgenDecVtxY);
   fChain->SetBranchAddress("lamgenDecVtxZ", &lamgenDecVtxZ, &b_lamgenDecVtxZ);
   fChain->SetBranchAddress("lamv0MatchStatus", &lamv0MatchStatus, &b_lamv0MatchStatus);
   fChain->SetBranchAddress("lamgenPDG", &lamgenPDG, &b_lamgenPDG);
   fChain->SetBranchAddress("lamgenPosDauPDG", &lamgenPosDauPDG, &b_lamgenPosDauPDG);
   fChain->SetBranchAddress("lamgenNegDauPDG", &lamgenNegDauPDG, &b_lamgenNegDauPDG);
   Notify();
}

Bool_t lam_efficiency::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void lam_efficiency::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t lam_efficiency::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef lam_efficiency_cxx
