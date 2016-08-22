//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun  8 11:23:43 2010 by ROOT version 5.22/00d
// from TTree ntuple/btosmumu ntuple
// found on file: /nfs/data37/cms/v0ntuple/data_7TeV_june4fixesGood_PARTIAL.root
//////////////////////////////////////////////////////////

#ifndef xiPlotter_h
#define xiPlotter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <vector>
#include <TGraphAsymmErrors.h>

class xiPlotter {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types

   UInt_t          nXI;
   UInt_t          run;
   UInt_t          event;
   Int_t           trigHF;
   Int_t           trigTech40_41;
   Int_t           trigTech36_39;
   Int_t           trigTech34;
   Int_t           trigScraping;
   Int_t           trigAlgo124;
   Int_t           trigHLTminBias;
   Int_t           nLooseTracks;
   Int_t           nHighPurTracks;
   Int_t           nLooseTracks_pt_500;
   Int_t           nHighPurTracks_pt_500;
   Int_t           processType;
   Float_t         priVtxX;
   Float_t         priVtxY;
   Float_t         priVtxZ;
   Float_t         priVtxXE;
   Float_t         priVtxYE;
   Float_t         priVtxZE;
   Float_t         priVtxCL;
   UInt_t          priNTrk;
   Int_t           priNTrk_weightGT_0_5;
   Float_t         bsX;
   Float_t         bsY;
   Float_t         bsZ;
   Float_t         bsXE;
   Float_t         bsYE;
   Float_t         bsZE;
   std::vector<float>   *priRefitVtxX;
   std::vector<float>   *priRefitVtxY;
   std::vector<float>   *priRefitVtxZ;
   std::vector<float>   *priRefitVtxXE;
   std::vector<float>   *priRefitVtxYE;
   std::vector<float>   *priRefitVtxZE;
   std::vector<float>   *priRefitVtxCL;
   std::vector<unsigned int> *priRefitNTrk;
   std::vector<float>   *xiMass;
   std::vector<float>   *omMass;
   std::vector<float>   *xiVtxCL;
   std::vector<float>   *xiPx;
   std::vector<float>   *xiPy;
   std::vector<float>   *xiPz;
   std::vector<double>  *xiPxE;
   std::vector<double>  *xiPyE;
   std::vector<double>  *xiPzE;
   std::vector<float>   *xiDecayVtxX;
   std::vector<float>   *xiDecayVtxY;
   std::vector<float>   *xiDecayVtxZ;
   std::vector<double>  *xiDecayVtxXE;
   std::vector<double>  *xiDecayVtxYE;
   std::vector<double>  *xiDecayVtxZE;
   std::vector<float>   *xictauB;
   std::vector<float>   *xictauBE;
   std::vector<float>   *xictau3D;
   std::vector<float>   *xictau3DE;
   std::vector<float>   *xictauMPV;
   std::vector<float>   *xictauMPVE;
   std::vector<float>   *VMass;
   std::vector<float>   *VVtxCL;
   std::vector<float>   *VPx;
   std::vector<float>   *VPy;
   std::vector<float>   *VPz;
   std::vector<float>   *VDecayVtxX;
   std::vector<float>   *VDecayVtxY;
   std::vector<float>   *VDecayVtxZ;
   std::vector<float>   *VDecayVtxXE;
   std::vector<float>   *VDecayVtxYE;
   std::vector<float>   *VDecayVtxZE;
   std::vector<float>   *batPiPx;
   std::vector<float>   *batPiPy;
   std::vector<float>   *batPiPz;
   std::vector<float>   *batPiEta;
   std::vector<float>   *batPiPhi;
   std::vector<float>   *batPiD0;
   std::vector<float>   *batPiD0E;
   std::vector<float>   *batPiDz;
   std::vector<float>   *batPiDzE;
   std::vector<float>   *batPiPVweight;
   std::vector<int>     *batPiQ;
   std::vector<int>     *batPiNValidHits;
   std::vector<int>     *batPiTrkQual;
   std::vector<int>     *nTracks;
   std::vector<float>   *VTrkPTrk1Mass;
   std::vector<float>   *VTrkPPx;
   std::vector<float>   *VTrkPPy;
   std::vector<float>   *VTrkPPz;
   std::vector<float>   *VTrkPEta;
   std::vector<float>   *VTrkPPhi;
   std::vector<float>   *VTrkPD0;
   std::vector<float>   *VTrkPD0E;
   std::vector<float>   *VTrkPDz;
   std::vector<float>   *VTrkPDzE;
   std::vector<int>     *VTrkPQ;
   std::vector<float>   *VTrkPPVweight;
   std::vector<int>     *VTrkPNValidHits;
   std::vector<float>   *VTrkPiMass;
   std::vector<float>   *VTrkPiPx;
   std::vector<float>   *VTrkPiPy;
   std::vector<float>   *VTrkPiPz;
   std::vector<float>   *VTrkPiEta;
   std::vector<float>   *VTrkPiPhi;
   std::vector<float>   *VTrkPiD0;
   std::vector<float>   *VTrkPiD0E;
   std::vector<float>   *VTrkPiDz;
   std::vector<float>   *VTrkPiDzE;
   std::vector<int>     *VTrkPiQ;
   std::vector<float>   *VTrkPiPVweight;
   std::vector<int>     *VTrkPiNValidHits;
   std::vector<float>   *VTransversePCAPrimary;
   std::vector<float>   *VTransversePCAPrimaryError;
   std::vector<float>   *VLongitudinalPCAPrimary;
   std::vector<float>   *VLongitudinalPCAPrimaryError;
   std::vector<float>   *VTransversePCABeamSpot;
   std::vector<float>   *VTransversePCABeamSpotError;
   std::vector<float>   *VFLSig2D;
   std::vector<float>   *VFLSig3D;
   std::vector<float>   *VCosThetaPAndveeVertexToPrimaryVector;
   std::vector<float>   *XiTransversePCAPrimary;
   std::vector<float>   *XiTransversePCAPrimaryError;
   std::vector<float>   *XiLongitudinalPCAPrimary;
   std::vector<float>   *XiLongitudinalPCAPrimaryError;
   std::vector<float>   *XiTransversePCABeamSpot;
   std::vector<float>   *XiTransversePCABeamSpotError;
   std::vector<float>   *XiCosThetaPAndLineOfFlight;
   std::vector<double>  *V3dIpWrtPrimary;
   std::vector<double>  *V3dIpWrtPrimarySig;
   std::vector<double>  *VTrkPi3DIp;
   std::vector<double>  *VTrkPi3DIpSig;
   std::vector<double>  *VTrkP3DIp;
   std::vector<double>  *VTrkP3DIpSig;
   std::vector<double>  *batPi3DIp;
   std::vector<double>  *batPi3DIpSig;
   std::vector<bool>    *pionInV0;
   UInt_t          nGenXi;
   std::vector<float>   *genBatPiEta;
   std::vector<float>   *genBatPiPhi;
   std::vector<float>   *genLamPiEta;
   std::vector<float>   *genLamPiPhi;
   std::vector<float>   *genLamPEta;
   std::vector<float>   *genLamPPhi;
   std::vector<float>   *genXiP;
   std::vector<float>   *genXiEta;
   std::vector<float>   *genXiL;
   std::vector<float>   *genXiR;
   std::vector<float>   *genXiDecayVX;
   std::vector<float>   *genXiDecayVY;
   std::vector<float>   *genXiDecayVZ;
   std::vector<float>   *genXiMomL;
   std::vector<int>     *genXiQ;

   // List of branches
   TBranch        *b_nXI;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_trigHF;   //!
   TBranch        *b_trigTech40_41;   //!
   TBranch        *b_trigTech36_39;   //!
   TBranch        *b_trigTech34;   //!
   TBranch        *b_trigScraping;   //!
   TBranch        *b_trigAlgo124;   //!
   TBranch        *b_trigHLTminBias;   //!
   TBranch        *b_nLooseTracks;   //!
   TBranch        *b_nHighPurTracks;   //!
   TBranch        *b_nLooseTracks_pt_500;   //!
   TBranch        *b_nHighPurTracks_pt_500;   //!
   TBranch        *b_processType;   //!
   TBranch        *b_priVtxX;   //!
   TBranch        *b_priVtxY;   //!
   TBranch        *b_priVtxZ;   //!
   TBranch        *b_priVtxXE;   //!
   TBranch        *b_priVtxYE;   //!
   TBranch        *b_priVtxZE;   //!
   TBranch        *b_priVtxCL;   //!
   TBranch        *b_priNTrk;   //!
   TBranch        *b_priNTrk_weightGT_0_5;   //!
   TBranch        *b_bsX;   //!
   TBranch        *b_bsY;   //!
   TBranch        *b_bsZ;   //!
   TBranch        *b_bsXE;   //!
   TBranch        *b_bsYE;   //!
   TBranch        *b_bsZE;   //!
   TBranch        *b_priRefitVtxX;   //!
   TBranch        *b_priRefitVtxY;   //!
   TBranch        *b_priRefitVtxZ;   //!
   TBranch        *b_priRefitVtxXE;   //!
   TBranch        *b_priRefitVtxYE;   //!
   TBranch        *b_priRefitVtxZE;   //!
   TBranch        *b_priRefitVtxCL;   //!
   TBranch        *b_priRefitNTrk;   //!
   TBranch        *b_xiMass;   //!
   TBranch        *b_omMass;   //!
   TBranch        *b_xiVtxCL;   //!
   TBranch        *b_xiPx;   //!
   TBranch        *b_xiPy;   //!
   TBranch        *b_xiPz;   //!
   TBranch        *b_xiPxE;   //!
   TBranch        *b_xiPyE;   //!
   TBranch        *b_xiPzE;   //!
   TBranch        *b_xiDecayVtxX;   //!
   TBranch        *b_xiDecayVtxY;   //!
   TBranch        *b_xiDecayVtxZ;   //!
   TBranch        *b_xiDecayVtxXE;   //!
   TBranch        *b_xiDecayVtxYE;   //!
   TBranch        *b_xiDecayVtxZE;   //!
   TBranch        *b_xictauB;   //!
   TBranch        *b_xictauBE;   //!
   TBranch        *b_xictau3D;   //!
   TBranch        *b_xictau3DE;   //!
   TBranch        *b_xictauMPV;   //!
   TBranch        *b_xictauMPVE;   //!
   TBranch        *b_VMass;   //!
   TBranch        *b_VVtxCL;   //!
   TBranch        *b_VPx;   //!
   TBranch        *b_VPy;   //!
   TBranch        *b_VPz;   //!
   TBranch        *b_VDecayVtxX;   //!
   TBranch        *b_VDecayVtxY;   //!
   TBranch        *b_VDecayVtxZ;   //!
   TBranch        *b_VDecayVtxXE;   //!
   TBranch        *b_VDecayVtxYE;   //!
   TBranch        *b_VDecayVtxZE;   //!
   TBranch        *b_batPiPx;   //!
   TBranch        *b_batPiPy;   //!
   TBranch        *b_batPiPz;   //!
   TBranch        *b_batPiEta;   //!
   TBranch        *b_batPiPhi;   //!
   TBranch        *b_batPiD0;   //!
   TBranch        *b_batPiD0E;   //!
   TBranch        *b_batPiDz;   //!
   TBranch        *b_batPiDzE;   //!
   TBranch        *b_batPiPVweight;   //!
   TBranch        *b_batPiQ;   //!
   TBranch        *b_batPiNValidHits;   //!
   TBranch        *b_batPiTrkQual;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_VTrkPTrk1Mass;   //!
   TBranch        *b_VTrkPPx;   //!
   TBranch        *b_VTrkPPy;   //!
   TBranch        *b_VTrkPPz;   //!
   TBranch        *b_VTrkPEta;   //!
   TBranch        *b_VTrkPPhi;   //!
   TBranch        *b_VTrkPD0;   //!
   TBranch        *b_VTrkPD0E;   //!
   TBranch        *b_VTrkPDz;   //!
   TBranch        *b_VTrkPDzE;   //!
   TBranch        *b_VTrkPQ;   //!
   TBranch        *b_VTrkPPVweight;   //!
   TBranch        *b_VTrkPNValidHits;   //!
   TBranch        *b_VTrkPiMass;   //!
   TBranch        *b_VTrkPiPx;   //!
   TBranch        *b_VTrkPiPy;   //!
   TBranch        *b_VTrkPiPz;   //!
   TBranch        *b_VTrkPiEta;   //!
   TBranch        *b_VTrkPiPhi;   //!
   TBranch        *b_VTrkPiD0;   //!
   TBranch        *b_VTrkPiD0E;   //!
   TBranch        *b_VTrkPiDz;   //!
   TBranch        *b_VTrkPiDzE;   //!
   TBranch        *b_VTrkPiQ;   //!
   TBranch        *b_VTrkPiPVweight;   //!
   TBranch        *b_VTrkPiNValidHits;   //!
   TBranch        *b_VTransversePCAPrimary;   //!
   TBranch        *b_VTransversePCAPrimaryError;   //!
   TBranch        *b_VLongitudinalPCAPrimary;   //!
   TBranch        *b_VLongitudinalPCAPrimaryError;   //!
   TBranch        *b_VTransversePCABeamSpot;   //!
   TBranch        *b_VTransversePCABeamSpotError;   //!
   TBranch        *b_VFLSig2D;   //!
   TBranch        *b_VFLSig3D;   //!
   TBranch        *b_VCosThetaPAndveeVertexToPrimaryVector;   //!
   TBranch        *b_XiTransversePCAPrimary;   //!
   TBranch        *b_XiTransversePCAPrimaryError;   //!
   TBranch        *b_XiLongitudinalPCAPrimary;   //!
   TBranch        *b_XiLongitudinalPCAPrimaryError;   //!
   TBranch        *b_XiTransversePCABeamSpot;   //!
   TBranch        *b_XiTransversePCABeamSpotError;   //!
   TBranch        *b_XiCosThetaPAndLineOfFlight;   //!
   TBranch        *b_V3dIpWrtPrimary;   //!
   TBranch        *b_V3dIpWrtPrimarySig;   //!
   TBranch        *b_VTrkPi3DIp;   //!
   TBranch        *b_VTrkPi3DIpSig;   //!
   TBranch        *b_VTrkP3DIp;   //!
   TBranch        *b_VTrkP3DIpSig;   //!
   TBranch        *b_batPi3DIp;   //!
   TBranch        *b_batPi3DIpSig;   //!
   TBranch        *b_pionInV0;   //!
   TBranch        *b_nGenXi;   //!
   TBranch        *b_genBatPiEta;   //!
   TBranch        *b_genBatPiPhi;   //!
   TBranch        *b_genLamPiEta;   //!
   TBranch        *b_genLamPiPhi;   //!
   TBranch        *b_genLamPEta;   //!
   TBranch        *b_genLamPPhi;   //!
   TBranch        *b_genXiP;   //!
   TBranch        *b_genXiEta;   //!
   TBranch        *b_genXiL;   //!
   TBranch        *b_genXiR;   //!
   TBranch        *b_genXiDecayVX;   //!
   TBranch        *b_genXiDecayVY;   //!
   TBranch        *b_genXiDecayVZ;   //!
   TBranch        *b_genXiMomL;   //!
   TBranch        *b_genXiQ;   //!

   xiPlotter(TTree *tree=0);
   virtual ~xiPlotter();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void SetNEvents(int);
   void bookHistos();
   void SetBins();
   void SetStyle();
   void SetOutputDir(std::string theOutputDir);
   void SetCutValues(float, float, float, float, float, float);
   void SetEffFile(std::string);
   void WriteHistos();
   void CleanVectors();
   void FitXi(TH1F*);
   void GetNumEff(TH1F*, TH1F*, TH1F*, bool=false);
   void DoNumFit(std::vector<TH1F*>, TH1F*, std::string);
   void GetAsymm(TH1F*, TH1F*, TGraphAsymmErrors*);
   
private:
   std::vector<TH1F*> thePtHistoVector;
   std::vector<TH1F*> thePosPtHistoVector;
   std::vector<TH1F*> theNegPtHistoVector;
   
   
   std::vector<TH1F*> theYHistoVector;
   std::vector<TH1F*> thePosYHistoVector;
   std::vector<TH1F*> theNegYHistoVector;
  
   
   std::vector<TH1F*> theCTauHistoVector;
   std::vector<TF1*> theFitVector;
   std::vector<float> thePtBins;
   std::vector<float> theYBins;
   std::vector<float> theCTauBins;
   int nEvents;
   
   float Xi3dIpCut;
   float batPi3dIpCut;
   float v0Pi3dIpCut; 
   float v0Proton3dIpCut;
   float XiSepCut;
   float vtxCLCut;
   
   
   double theFitYield;
   double theFitYieldError;
   
   std::string theOutDir;
   std::string theEffFileName;
   TFile* theOutRootFile;
   TFile* theEffFile;
   
   TH1F* data;
   TH1F* dataOm;
   TH1F* dataBad;
   
   TH1F* numEffpTScaled;
   TH1F* numEffpTScaledUncor;
   TH1F* numFitPt;
   
   TH1F* numEffPosXipTScaled;
   TH1F* numEffPosXipTScaledUncor;
   TH1F* numFitPosXiPt;
   
   TH1F* numEffNegXipTScaled;
   TH1F* numEffNegXipTScaledUncor;
   TH1F* numFitNegXiPt;
   
   TH1F* posNumFitPt;
   TH1F* negNumFitPt;
   TH1F* asymmFitPt;
   
   TH1F* numEffYScaled;
   TH1F* numEffYScaledUncor;
   TH1F* numFitY;
   
   TH1F* posNumFitY;
   TH1F* negNumFitY;
   TH1F* asymmFitY;
   
   TH1F* numEffPosXiYScaled;
   TH1F* numEffPosXiYScaledUncor;
   TH1F* numFitPosXiY;
   
   TH1F* numEffNegXiYScaled;
   TH1F* numEffNegXiYScaledUncor;
   TH1F* numFitNegXiY;
   
   
   TH1F* numEffCtauScaled;
   TH1F* numEffCtauScaledUncor;
   TH1F* numFitCTau;
   
   TH1F* effCtau;
   TH1F* effCtauUncor;
   
   TH1F* effY;
   TH1F* effYUncor;
   
   TH1F* effPt;
   TH1F* effPtUncor;
   
   TH1F* effPosXiPt;
   TH1F* effPosXiPtUncor;
   
   TH1F* effNegXiPt;
   TH1F* effNegXiPtUncor;
   
   TH1F* effPosXiY;
   TH1F* effPosXiYUncor;
   
   TH1F* effNegXiY;
   TH1F* effNegXiYUncor;
  
};

#endif

#ifdef xiPlotter_cxx
xiPlotter::xiPlotter(TTree *tree){
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/nfs/data37/cms/v0ntuple/data_7TeV_june4fixesGood_PARTIAL.root");
      if (!f) {
         f = new TFile("/nfs/data37/cms/v0ntuple/data_7TeV_june4fixesGood_PARTIAL.root");
         f->cd("/nfs/data37/cms/v0ntuple/data_7TeV_june4fixesGood_PARTIAL.root:/ntupleXiOm");
      }
      tree = (TTree*)gDirectory->Get("ntuple");

   }
   
   Init(tree);
   
   nEvents = -1;
   theOutDir = "plots"; 
   
   Xi3dIpCut = 3.;
   batPi3dIpCut = 3.;
   v0Pi3dIpCut = 3.; 
   v0Proton3dIpCut = 3.;
   XiSepCut = 4.;
   vtxCLCut = 0.01;
   
}

xiPlotter::~xiPlotter(){
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t xiPlotter::GetEntry(Long64_t entry){
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t xiPlotter::LoadTree(Long64_t entry){
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

void xiPlotter::Init(TTree *tree){
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer

   priRefitVtxX = 0;
   priRefitVtxY = 0;
   priRefitVtxZ = 0;
   priRefitVtxXE = 0;
   priRefitVtxYE = 0;
   priRefitVtxZE = 0;
   priRefitVtxCL = 0;
   priRefitNTrk = 0;
   xiMass = 0;
   omMass = 0;
   xiVtxCL = 0;
   xiPx = 0;
   xiPy = 0;
   xiPz = 0;
   xiPxE = 0;
   xiPyE = 0;
   xiPzE = 0;
   xiDecayVtxX = 0;
   xiDecayVtxY = 0;
   xiDecayVtxZ = 0;
   xiDecayVtxXE = 0;
   xiDecayVtxYE = 0;
   xiDecayVtxZE = 0;
   xictauB = 0;
   xictauBE = 0;
   xictau3D = 0;
   xictau3DE = 0;
   xictauMPV = 0;
   xictauMPVE = 0;
   VMass = 0;
   VVtxCL = 0;
   VPx = 0;
   VPy = 0;
   VPz = 0;
   VDecayVtxX = 0;
   VDecayVtxY = 0;
   VDecayVtxZ = 0;
   VDecayVtxXE = 0;
   VDecayVtxYE = 0;
   VDecayVtxZE = 0;
   batPiPx = 0;
   batPiPy = 0;
   batPiPz = 0;
   batPiEta = 0;
   batPiPhi = 0;
   batPiD0 = 0;
   batPiD0E = 0;
   batPiDz = 0;
   batPiDzE = 0;
   batPiPVweight = 0;
   batPiQ = 0;
   batPiNValidHits = 0;
   batPiTrkQual = 0;
   nTracks = 0;
   VTrkPTrk1Mass = 0;
   VTrkPPx = 0;
   VTrkPPy = 0;
   VTrkPPz = 0;
   VTrkPEta = 0;
   VTrkPPhi = 0;
   VTrkPD0 = 0;
   VTrkPD0E = 0;
   VTrkPDz = 0;
   VTrkPDzE = 0;
   VTrkPQ = 0;
   VTrkPPVweight = 0;
   VTrkPNValidHits = 0;
   VTrkPiMass = 0;
   VTrkPiPx = 0;
   VTrkPiPy = 0;
   VTrkPiPz = 0;
   VTrkPiEta = 0;
   VTrkPiPhi = 0;
   VTrkPiD0 = 0;
   VTrkPiD0E = 0;
   VTrkPiDz = 0;
   VTrkPiDzE = 0;
   VTrkPiQ = 0;
   VTrkPiPVweight = 0;
   VTrkPiNValidHits = 0;
   VTransversePCAPrimary = 0;
   VTransversePCAPrimaryError = 0;
   VLongitudinalPCAPrimary = 0;
   VLongitudinalPCAPrimaryError = 0;
   VTransversePCABeamSpot = 0;
   VTransversePCABeamSpotError = 0;
   VFLSig2D = 0;
   VFLSig3D = 0;
   VCosThetaPAndveeVertexToPrimaryVector = 0;
   XiTransversePCAPrimary = 0;
   XiTransversePCAPrimaryError = 0;
   XiLongitudinalPCAPrimary = 0;
   XiLongitudinalPCAPrimaryError = 0;
   XiTransversePCABeamSpot = 0;
   XiTransversePCABeamSpotError = 0;
   XiCosThetaPAndLineOfFlight = 0;
   V3dIpWrtPrimary = 0;
   V3dIpWrtPrimarySig = 0;
   VTrkPi3DIp = 0;
   VTrkPi3DIpSig = 0;
   VTrkP3DIp = 0;
   VTrkP3DIpSig = 0;
   batPi3DIp = 0;
   batPi3DIpSig = 0;
   pionInV0 = 0;
   genBatPiEta = 0;
   genBatPiPhi = 0;
   genLamPiEta = 0;
   genLamPiPhi = 0;
   genLamPEta = 0;
   genLamPPhi = 0;
   genXiP = 0;
   genXiEta = 0;
   genXiL = 0;
   genXiR = 0;
   genXiDecayVX = 0;
   genXiDecayVY = 0;
   genXiDecayVZ = 0;
   genXiMomL = 0;
   genXiQ = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nXI", &nXI, &b_nXI);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("trigHF", &trigHF, &b_trigHF);
   fChain->SetBranchAddress("trigTech40_41", &trigTech40_41, &b_trigTech40_41);
   fChain->SetBranchAddress("trigTech36_39", &trigTech36_39, &b_trigTech36_39);
   fChain->SetBranchAddress("trigTech34", &trigTech34, &b_trigTech34);
   fChain->SetBranchAddress("trigScraping", &trigScraping, &b_trigScraping);
   fChain->SetBranchAddress("trigAlgo124", &trigAlgo124, &b_trigAlgo124);
   fChain->SetBranchAddress("trigHLTminBias", &trigHLTminBias, &b_trigHLTminBias);
   fChain->SetBranchAddress("nLooseTracks", &nLooseTracks, &b_nLooseTracks);
   fChain->SetBranchAddress("nHighPurTracks", &nHighPurTracks, &b_nHighPurTracks);
   fChain->SetBranchAddress("nLooseTracks_pt_500", &nLooseTracks_pt_500, &b_nLooseTracks_pt_500);
   fChain->SetBranchAddress("nHighPurTracks_pt_500", &nHighPurTracks_pt_500, &b_nHighPurTracks_pt_500);
   fChain->SetBranchAddress("processType", &processType, &b_processType);
   fChain->SetBranchAddress("priVtxX", &priVtxX, &b_priVtxX);
   fChain->SetBranchAddress("priVtxY", &priVtxY, &b_priVtxY);
   fChain->SetBranchAddress("priVtxZ", &priVtxZ, &b_priVtxZ);
   fChain->SetBranchAddress("priVtxXE", &priVtxXE, &b_priVtxXE);
   fChain->SetBranchAddress("priVtxYE", &priVtxYE, &b_priVtxYE);
   fChain->SetBranchAddress("priVtxZE", &priVtxZE, &b_priVtxZE);
   fChain->SetBranchAddress("priVtxCL", &priVtxCL, &b_priVtxCL);
   fChain->SetBranchAddress("priNTrk", &priNTrk, &b_priNTrk);
   fChain->SetBranchAddress("priNTrk_weightGT_0_5", &priNTrk_weightGT_0_5, &b_priNTrk_weightGT_0_5);
   fChain->SetBranchAddress("bsX", &bsX, &b_bsX);
   fChain->SetBranchAddress("bsY", &bsY, &b_bsY);
   fChain->SetBranchAddress("bsZ", &bsZ, &b_bsZ);
   fChain->SetBranchAddress("bsXE", &bsXE, &b_bsXE);
   fChain->SetBranchAddress("bsYE", &bsYE, &b_bsYE);
   fChain->SetBranchAddress("bsZE", &bsZE, &b_bsZE);
   fChain->SetBranchAddress("priRefitVtxX", &priRefitVtxX, &b_priRefitVtxX);
   fChain->SetBranchAddress("priRefitVtxY", &priRefitVtxY, &b_priRefitVtxY);
   fChain->SetBranchAddress("priRefitVtxZ", &priRefitVtxZ, &b_priRefitVtxZ);
   fChain->SetBranchAddress("priRefitVtxXE", &priRefitVtxXE, &b_priRefitVtxXE);
   fChain->SetBranchAddress("priRefitVtxYE", &priRefitVtxYE, &b_priRefitVtxYE);
   fChain->SetBranchAddress("priRefitVtxZE", &priRefitVtxZE, &b_priRefitVtxZE);
   fChain->SetBranchAddress("priRefitVtxCL", &priRefitVtxCL, &b_priRefitVtxCL);
   fChain->SetBranchAddress("priRefitNTrk", &priRefitNTrk, &b_priRefitNTrk);
   fChain->SetBranchAddress("xiMass", &xiMass, &b_xiMass);
   fChain->SetBranchAddress("omMass", &omMass, &b_omMass);
   fChain->SetBranchAddress("xiVtxCL", &xiVtxCL, &b_xiVtxCL);
   fChain->SetBranchAddress("xiPx", &xiPx, &b_xiPx);
   fChain->SetBranchAddress("xiPy", &xiPy, &b_xiPy);
   fChain->SetBranchAddress("xiPz", &xiPz, &b_xiPz);
   fChain->SetBranchAddress("xiPxE", &xiPxE, &b_xiPxE);
   fChain->SetBranchAddress("xiPyE", &xiPyE, &b_xiPyE);
   fChain->SetBranchAddress("xiPzE", &xiPzE, &b_xiPzE);
   fChain->SetBranchAddress("xiDecayVtxX", &xiDecayVtxX, &b_xiDecayVtxX);
   fChain->SetBranchAddress("xiDecayVtxY", &xiDecayVtxY, &b_xiDecayVtxY);
   fChain->SetBranchAddress("xiDecayVtxZ", &xiDecayVtxZ, &b_xiDecayVtxZ);
   fChain->SetBranchAddress("xiDecayVtxXE", &xiDecayVtxXE, &b_xiDecayVtxXE);
   fChain->SetBranchAddress("xiDecayVtxYE", &xiDecayVtxYE, &b_xiDecayVtxYE);
   fChain->SetBranchAddress("xiDecayVtxZE", &xiDecayVtxZE, &b_xiDecayVtxZE);
   fChain->SetBranchAddress("xictauB", &xictauB, &b_xictauB);
   fChain->SetBranchAddress("xictauBE", &xictauBE, &b_xictauBE);
   fChain->SetBranchAddress("xictau3D", &xictau3D, &b_xictau3D);
   fChain->SetBranchAddress("xictau3DE", &xictau3DE, &b_xictau3DE);
   fChain->SetBranchAddress("xictauMPV", &xictauMPV, &b_xictauMPV);
   fChain->SetBranchAddress("xictauMPVE", &xictauMPVE, &b_xictauMPVE);
   fChain->SetBranchAddress("VMass", &VMass, &b_VMass);
   fChain->SetBranchAddress("VVtxCL", &VVtxCL, &b_VVtxCL);
   fChain->SetBranchAddress("VPx", &VPx, &b_VPx);
   fChain->SetBranchAddress("VPy", &VPy, &b_VPy);
   fChain->SetBranchAddress("VPz", &VPz, &b_VPz);
   fChain->SetBranchAddress("VDecayVtxX", &VDecayVtxX, &b_VDecayVtxX);
   fChain->SetBranchAddress("VDecayVtxY", &VDecayVtxY, &b_VDecayVtxY);
   fChain->SetBranchAddress("VDecayVtxZ", &VDecayVtxZ, &b_VDecayVtxZ);
   fChain->SetBranchAddress("VDecayVtxXE", &VDecayVtxXE, &b_VDecayVtxXE);
   fChain->SetBranchAddress("VDecayVtxYE", &VDecayVtxYE, &b_VDecayVtxYE);
   fChain->SetBranchAddress("VDecayVtxZE", &VDecayVtxZE, &b_VDecayVtxZE);
   fChain->SetBranchAddress("batPiPx", &batPiPx, &b_batPiPx);
   fChain->SetBranchAddress("batPiPy", &batPiPy, &b_batPiPy);
   fChain->SetBranchAddress("batPiPz", &batPiPz, &b_batPiPz);
   fChain->SetBranchAddress("batPiEta", &batPiEta, &b_batPiEta);
   fChain->SetBranchAddress("batPiPhi", &batPiPhi, &b_batPiPhi);
   fChain->SetBranchAddress("batPiD0", &batPiD0, &b_batPiD0);
   fChain->SetBranchAddress("batPiD0E", &batPiD0E, &b_batPiD0E);
   fChain->SetBranchAddress("batPiDz", &batPiDz, &b_batPiDz);
   fChain->SetBranchAddress("batPiDzE", &batPiDzE, &b_batPiDzE);
   fChain->SetBranchAddress("batPiPVweight", &batPiPVweight, &b_batPiPVweight);
   fChain->SetBranchAddress("batPiQ", &batPiQ, &b_batPiQ);
   fChain->SetBranchAddress("batPiNValidHits", &batPiNValidHits, &b_batPiNValidHits);
   fChain->SetBranchAddress("batPiTrkQual", &batPiTrkQual, &b_batPiTrkQual);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("VTrkPTrk1Mass", &VTrkPTrk1Mass, &b_VTrkPTrk1Mass);
   fChain->SetBranchAddress("VTrkPPx", &VTrkPPx, &b_VTrkPPx);
   fChain->SetBranchAddress("VTrkPPy", &VTrkPPy, &b_VTrkPPy);
   fChain->SetBranchAddress("VTrkPPz", &VTrkPPz, &b_VTrkPPz);
   fChain->SetBranchAddress("VTrkPEta", &VTrkPEta, &b_VTrkPEta);
   fChain->SetBranchAddress("VTrkPPhi", &VTrkPPhi, &b_VTrkPPhi);
   fChain->SetBranchAddress("VTrkPD0", &VTrkPD0, &b_VTrkPD0);
   fChain->SetBranchAddress("VTrkPD0E", &VTrkPD0E, &b_VTrkPD0E);
   fChain->SetBranchAddress("VTrkPDz", &VTrkPDz, &b_VTrkPDz);
   fChain->SetBranchAddress("VTrkPDzE", &VTrkPDzE, &b_VTrkPDzE);
   fChain->SetBranchAddress("VTrkPQ", &VTrkPQ, &b_VTrkPQ);
   fChain->SetBranchAddress("VTrkPPVweight", &VTrkPPVweight, &b_VTrkPPVweight);
   fChain->SetBranchAddress("VTrkPNValidHits", &VTrkPNValidHits, &b_VTrkPNValidHits);
   fChain->SetBranchAddress("VTrkPiMass", &VTrkPiMass, &b_VTrkPiMass);
   fChain->SetBranchAddress("VTrkPiPx", &VTrkPiPx, &b_VTrkPiPx);
   fChain->SetBranchAddress("VTrkPiPy", &VTrkPiPy, &b_VTrkPiPy);
   fChain->SetBranchAddress("VTrkPiPz", &VTrkPiPz, &b_VTrkPiPz);
   fChain->SetBranchAddress("VTrkPiEta", &VTrkPiEta, &b_VTrkPiEta);
   fChain->SetBranchAddress("VTrkPiPhi", &VTrkPiPhi, &b_VTrkPiPhi);
   fChain->SetBranchAddress("VTrkPiD0", &VTrkPiD0, &b_VTrkPiD0);
   fChain->SetBranchAddress("VTrkPiD0E", &VTrkPiD0E, &b_VTrkPiD0E);
   fChain->SetBranchAddress("VTrkPiDz", &VTrkPiDz, &b_VTrkPiDz);
   fChain->SetBranchAddress("VTrkPiDzE", &VTrkPiDzE, &b_VTrkPiDzE);
   fChain->SetBranchAddress("VTrkPiQ", &VTrkPiQ, &b_VTrkPiQ);
   fChain->SetBranchAddress("VTrkPiPVweight", &VTrkPiPVweight, &b_VTrkPiPVweight);
   fChain->SetBranchAddress("VTrkPiNValidHits", &VTrkPiNValidHits, &b_VTrkPiNValidHits);
   fChain->SetBranchAddress("VTransversePCAPrimary", &VTransversePCAPrimary, &b_VTransversePCAPrimary);
   fChain->SetBranchAddress("VTransversePCAPrimaryError", &VTransversePCAPrimaryError, &b_VTransversePCAPrimaryError);
   fChain->SetBranchAddress("VLongitudinalPCAPrimary", &VLongitudinalPCAPrimary, &b_VLongitudinalPCAPrimary);
   fChain->SetBranchAddress("VLongitudinalPCAPrimaryError", &VLongitudinalPCAPrimaryError, &b_VLongitudinalPCAPrimaryError);
   fChain->SetBranchAddress("VTransversePCABeamSpot", &VTransversePCABeamSpot, &b_VTransversePCABeamSpot);
   fChain->SetBranchAddress("VTransversePCABeamSpotError", &VTransversePCABeamSpotError, &b_VTransversePCABeamSpotError);
   fChain->SetBranchAddress("VFLSig2D", &VFLSig2D, &b_VFLSig2D);
   fChain->SetBranchAddress("VFLSig3D", &VFLSig3D, &b_VFLSig3D);
   fChain->SetBranchAddress("VCosThetaPAndveeVertexToPrimaryVector", &VCosThetaPAndveeVertexToPrimaryVector, &b_VCosThetaPAndveeVertexToPrimaryVector);
   fChain->SetBranchAddress("XiTransversePCAPrimary", &XiTransversePCAPrimary, &b_XiTransversePCAPrimary);
   fChain->SetBranchAddress("XiTransversePCAPrimaryError", &XiTransversePCAPrimaryError, &b_XiTransversePCAPrimaryError);
   fChain->SetBranchAddress("XiLongitudinalPCAPrimary", &XiLongitudinalPCAPrimary, &b_XiLongitudinalPCAPrimary);
   fChain->SetBranchAddress("XiLongitudinalPCAPrimaryError", &XiLongitudinalPCAPrimaryError, &b_XiLongitudinalPCAPrimaryError);
   fChain->SetBranchAddress("XiTransversePCABeamSpot", &XiTransversePCABeamSpot, &b_XiTransversePCABeamSpot);
   fChain->SetBranchAddress("XiTransversePCABeamSpotError", &XiTransversePCABeamSpotError, &b_XiTransversePCABeamSpotError);
   fChain->SetBranchAddress("XiCosThetaPAndLineOfFlight", &XiCosThetaPAndLineOfFlight, &b_XiCosThetaPAndLineOfFlight);
   fChain->SetBranchAddress("V3dIpWrtPrimary", &V3dIpWrtPrimary, &b_V3dIpWrtPrimary);
   fChain->SetBranchAddress("V3dIpWrtPrimarySig", &V3dIpWrtPrimarySig, &b_V3dIpWrtPrimarySig);
   fChain->SetBranchAddress("VTrkPi3DIp", &VTrkPi3DIp, &b_VTrkPi3DIp);
   fChain->SetBranchAddress("VTrkPi3DIpSig", &VTrkPi3DIpSig, &b_VTrkPi3DIpSig);
   fChain->SetBranchAddress("VTrkP3DIp", &VTrkP3DIp, &b_VTrkP3DIp);
   fChain->SetBranchAddress("VTrkP3DIpSig", &VTrkP3DIpSig, &b_VTrkP3DIpSig);
   fChain->SetBranchAddress("batPi3DIp", &batPi3DIp, &b_batPi3DIp);
   fChain->SetBranchAddress("batPi3DIpSig", &batPi3DIpSig, &b_batPi3DIpSig);
   fChain->SetBranchAddress("pionInV0", &pionInV0, &b_pionInV0);
   fChain->SetBranchAddress("nGenXi", &nGenXi, &b_nGenXi);
   fChain->SetBranchAddress("genBatPiEta", &genBatPiEta, &b_genBatPiEta);
   fChain->SetBranchAddress("genBatPiPhi", &genBatPiPhi, &b_genBatPiPhi);
   fChain->SetBranchAddress("genLamPiEta", &genLamPiEta, &b_genLamPiEta);
   fChain->SetBranchAddress("genLamPiPhi", &genLamPiPhi, &b_genLamPiPhi);
   fChain->SetBranchAddress("genLamPEta", &genLamPEta, &b_genLamPEta);
   fChain->SetBranchAddress("genLamPPhi", &genLamPPhi, &b_genLamPPhi);
   fChain->SetBranchAddress("genXiP", &genXiP, &b_genXiP);
   fChain->SetBranchAddress("genXiEta", &genXiEta, &b_genXiEta);
   fChain->SetBranchAddress("genXiL", &genXiL, &b_genXiL);
   fChain->SetBranchAddress("genXiR", &genXiR, &b_genXiR);
   fChain->SetBranchAddress("genXiDecayVX", &genXiDecayVX, &b_genXiDecayVX);
   fChain->SetBranchAddress("genXiDecayVY", &genXiDecayVY, &b_genXiDecayVY);
   fChain->SetBranchAddress("genXiDecayVZ", &genXiDecayVZ, &b_genXiDecayVZ);
   fChain->SetBranchAddress("genXiMomL", &genXiMomL, &b_genXiMomL);
   fChain->SetBranchAddress("genXiQ", &genXiQ, &b_genXiQ);
   Notify();
}

Bool_t xiPlotter::Notify(){
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void xiPlotter::Show(Long64_t entry){
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

void xiPlotter::SetNEvents(int theNumberOfEvents){
  nEvents = theNumberOfEvents;
}

void xiPlotter::WriteHistos(){
  theOutRootFile->Write();
  theOutRootFile->Close();
}

void xiPlotter::SetCutValues(float xi3dCut, float batpi3dCut, float v0pi3dCut, float v0proton3dCut, float sepCut, float vxtclCut){
   
   Xi3dIpCut = xi3dCut;
   batPi3dIpCut = batpi3dCut;
   v0Pi3dIpCut = v0pi3dCut; 
   v0Proton3dIpCut = v0proton3dCut;
   XiSepCut = sepCut;
   vtxCLCut = vxtclCut;
  
}

void xiPlotter::SetStyle(){
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetFillColor(0);
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}

Int_t xiPlotter::Cut(Long64_t entry){
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef xiPlotter_cxx
