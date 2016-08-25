//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul  8 19:18:37 2012 by ROOT version 5.27/06b
// from TTree HltTree/
// found on file: /d102/yjlee/hiForest2MC/Pythia80_HydjetDrum_mix01_HiForest2_v20.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Hlts {
public :
   Hlts(){};
   ~Hlts(){};

   // Declaration of leaf types
   Int_t           NL1IsolEm;
   Float_t         L1IsolEmEt[58];   //[NL1IsolEm]
   Float_t         L1IsolEmE[58];   //[NL1IsolEm]
   Float_t         L1IsolEmEta[58];   //[NL1IsolEm]
   Float_t         L1IsolEmPhi[58];   //[NL1IsolEm]
   Int_t           NL1NIsolEm;
   Float_t         L1NIsolEmEt[58];   //[NL1NIsolEm]
   Float_t         L1NIsolEmE[58];   //[NL1NIsolEm]
   Float_t         L1NIsolEmEta[58];   //[NL1NIsolEm]
   Float_t         L1NIsolEmPhi[58];   //[NL1NIsolEm]
   Int_t           NL1Mu;
   Float_t         L1MuPt[58];   //[NL1Mu]
   Float_t         L1MuE[58];   //[NL1Mu]
   Float_t         L1MuEta[58];   //[NL1Mu]
   Float_t         L1MuPhi[58];   //[NL1Mu]
   Int_t           L1MuIsol[58];   //[NL1Mu]
   Int_t           L1MuMip[58];   //[NL1Mu]
   Int_t           L1MuFor[58];   //[NL1Mu]
   Int_t           L1MuRPC[58];   //[NL1Mu]
   Int_t           L1MuQal[58];   //[NL1Mu]
   Int_t           L1MuChg[58];   //[NL1Mu]
   Int_t           NL1CenJet;
   Float_t         L1CenJetEt[58];   //[NL1CenJet]
   Float_t         L1CenJetE[58];   //[NL1CenJet]
   Float_t         L1CenJetEta[58];   //[NL1CenJet]
   Float_t         L1CenJetPhi[58];   //[NL1CenJet]
   Int_t           NL1ForJet;
   Float_t         L1ForJetEt[58];   //[NL1ForJet]
   Float_t         L1ForJetE[58];   //[NL1ForJet]
   Float_t         L1ForJetEta[58];   //[NL1ForJet]
   Float_t         L1ForJetPhi[58];   //[NL1ForJet]
   Int_t           NL1Tau;
   Float_t         L1TauEt[58];   //[NL1Tau]
   Float_t         L1TauE[58];   //[NL1Tau]
   Float_t         L1TauEta[58];   //[NL1Tau]
   Float_t         L1TauPhi[58];   //[NL1Tau]
   Float_t         L1Met;
   Float_t         L1MetPhi;
   Float_t         L1EtTot;
   Float_t         L1Mht;
   Float_t         L1MhtPhi;
   Float_t         L1EtHad;
   Int_t           L1HfRing1EtSumPositiveEta;
   Int_t           L1HfRing2EtSumPositiveEta;
   Int_t           L1HfRing1EtSumNegativeEta;
   Int_t           L1HfRing2EtSumNegativeEta;
   Int_t           L1HfTowerCountPositiveEtaRing1;
   Int_t           L1HfTowerCountNegativeEtaRing1;
   Int_t           L1HfTowerCountPositiveEtaRing2;
   Int_t           L1HfTowerCountNegativeEtaRing2;
   Int_t           Run;
   Int_t           Event;
   Int_t           LumiBlock;
   Int_t           Bx;
   Int_t           Orbit;
   Double_t        AvgInstDelLumi;
   Int_t           raw2digi_step;
   Int_t           raw2digi_step_Prescl;
   Int_t           reconstruction_step;
   Int_t           reconstruction_step_Prescl;
   Int_t           HLTriggerFirstPath;
   Int_t           HLTriggerFirstPath_Prescl;
   Int_t           HLT_HIDTCalibration_v1;
   Int_t           HLT_HIDTCalibration_v1_Prescl;
   Int_t           HLT_HIEcalCalibration_v1;
   Int_t           HLT_HIEcalCalibration_v1_Prescl;
   Int_t           HLT_HIHcalCalibration_v1;
   Int_t           HLT_HIHcalCalibration_v1_Prescl;
   Int_t           HLT_HIZeroBias_v1;
   Int_t           HLT_HIZeroBias_v1_Prescl;
   Int_t           HLT_HIZeroBiasXOR_v1;
   Int_t           HLT_HIZeroBiasXOR_v1_Prescl;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v1;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIMinBiasBSC_v1;
   Int_t           HLT_HIMinBiasBSC_v1_Prescl;
   Int_t           HLT_HIMinBiasBSC_OR_v1;
   Int_t           HLT_HIMinBiasBSC_OR_v1_Prescl;
   Int_t           HLT_HIMinBiasHF_v1;
   Int_t           HLT_HIMinBiasHF_v1_Prescl;
   Int_t           HLT_HIMinBiasHf_OR_v1;
   Int_t           HLT_HIMinBiasHf_OR_v1_Prescl;
   Int_t           HLT_HIMinBiasHfOrBSC_v1;
   Int_t           HLT_HIMinBiasHfOrBSC_v1_Prescl;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v1;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_v1;
   Int_t           HLT_HIMinBiasZDC_Calo_v1_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v1;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIBptxXOR_v1;
   Int_t           HLT_HIBptxXOR_v1_Prescl;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v1;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl;
   Int_t           HLT_HIL1DoubleMuOpen_v1;
   Int_t           HLT_HIL1DoubleMuOpen_v1_Prescl;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v1;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v1_Prescl;
   Int_t           HLT_HIL2Mu3_v1;
   Int_t           HLT_HIL2Mu3_v1_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ_v1;
   Int_t           HLT_HIL2Mu3_NHitQ_v1_Prescl;
   Int_t           HLT_HIL2Mu7_v1;
   Int_t           HLT_HIL2Mu7_v1_Prescl;
   Int_t           HLT_HIL2Mu15_v1;
   Int_t           HLT_HIL2Mu15_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu0_v1;
   Int_t           HLT_HIL2DoubleMu0_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v1;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu3_v1;
   Int_t           HLT_HIL2DoubleMu3_v1_Prescl;
   Int_t           HLT_HIL3Mu3_v1;
   Int_t           HLT_HIL3Mu3_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_v1;
   Int_t           HLT_HIL3DoubleMuOpen_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v1;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v1;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v1;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl;
   Int_t           HLT_HISinglePhoton15_v1;
   Int_t           HLT_HISinglePhoton15_v1_Prescl;
   Int_t           HLT_HISinglePhoton20_v1;
   Int_t           HLT_HISinglePhoton20_v1_Prescl;
   Int_t           HLT_HISinglePhoton30_v1;
   Int_t           HLT_HISinglePhoton30_v1_Prescl;
   Int_t           HLT_HISinglePhoton40_v1;
   Int_t           HLT_HISinglePhoton40_v1_Prescl;
   Int_t           HLT_HIPhoton10_Photon15_v1;
   Int_t           HLT_HIPhoton10_Photon15_v1_Prescl;
   Int_t           HLT_HIPhoton15_Photon20_v1;
   Int_t           HLT_HIPhoton15_Photon20_v1_Prescl;
   Int_t           HLT_HIDoublePhoton10_v1;
   Int_t           HLT_HIDoublePhoton10_v1_Prescl;
   Int_t           HLT_HIDoublePhoton15_v1;
   Int_t           HLT_HIDoublePhoton15_v1_Prescl;
   Int_t           HLT_HIDoublePhoton20_v1;
   Int_t           HLT_HIDoublePhoton20_v1_Prescl;
   Int_t           HLT_HIJet55_v1;
   Int_t           HLT_HIJet55_v1_Prescl;
   Int_t           HLT_HIJet65_v1;
   Int_t           HLT_HIJet65_v1_Prescl;
   Int_t           HLT_Jet60_v1;
   Int_t           HLT_PAJet20_NoJetID_v1;
   Int_t           HLT_PAJet40_NoJetID_v1;
   Int_t           HLT_PAJet60_NoJetID_v1;
   Int_t           HLT_PAJet80_NoJetID_v1;
   Int_t           HLT_PAJet100_NoJetID_v1;
   Int_t           HLT_PAZeroBias_v1;
   Int_t           HLT_PAZeroBias_v1_Prescl;
   Int_t           HLT_PAZeroBiasPixel_SingleTrack_v1;
   Int_t           HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_PAFullTrack12_v2;
   Int_t           HLT_PAFullTrack20_v2;
   Int_t           HLT_PAFullTrack30_v2;
   Int_t           HLT_HIJet80_v1;
   Int_t           HLT_HIJet80_v1_Prescl;
   Int_t           HLT_HIJet95_v1;
   Int_t           HLT_HIJet95_v1_Prescl;
   Int_t           HLT_HIDiJet55_v1;
   Int_t           HLT_HIDiJet55_v1_Prescl;
   Int_t           HLT_HIJet65_Jet55_v1;
   Int_t           HLT_HIJet65_Jet55_v1_Prescl;
   Int_t           HLT_HIJetE30_NoBPTX_v1;
   Int_t           HLT_HIJetE30_NoBPTX_v1_Prescl;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v1;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl;
   Int_t           HLT_HIActivityHF_Coincidence3_v1;
   Int_t           HLT_HIActivityHF_Coincidence3_v1_Prescl;
   Int_t           HLT_HIActivityHF_Single3_v1;
   Int_t           HLT_HIActivityHF_Single3_v1_Prescl;
   Int_t           HLT_HIClusterVertexCompatibility_v1;
   Int_t           HLT_HIClusterVertexCompatibility_v1_Prescl;
   Int_t           HLT_HICentralityVeto_v1;
   Int_t           HLT_HICentralityVeto_v1_Prescl;
   Int_t           HLT_HIFullTrack12_L1Central_v1;
   Int_t           HLT_HIFullTrack12_L1Central_v1_Prescl;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v1;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v1_Prescl;
   Int_t           HLT_HIFullTrack14_L1Central_v1;
   Int_t           HLT_HIFullTrack14_L1Central_v1_Prescl;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v1;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v1_Prescl;
   Int_t           HLT_HIFullTrack20_L1Central_v1;
   Int_t           HLT_HIFullTrack20_L1Central_v1_Prescl;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v1;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v1_Prescl;
   Int_t           HLT_HIFullTrack25_L1Central_v1;
   Int_t           HLT_HIFullTrack25_L1Central_v1_Prescl;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v1;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v1_Prescl;
   Int_t           HLT_HIRandom_v1;
   Int_t           HLT_HIRandom_v1_Prescl;
   Int_t           HLT_HIUCC010_v1;
   Int_t           HLT_HIUCC010_v1_Prescl;
   Int_t           HLT_HIUCC015_v1;
   Int_t           HLT_HIUCC015_v1_Prescl;
   Int_t           HLT_HICentral10_v1;
   Int_t           HLT_HICentral10_v1_Prescl;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl;
   Int_t           HLTriggerFinalPath;
   Int_t           HLTriggerFinalPath_Prescl;
   Int_t           HLT_HIMET120_v1;
   Int_t           HLT_HIMET120_v1_Prescl;
   Int_t           HLT_HIMET200_v1;
   Int_t           HLT_HIMET200_v1_Prescl;
   Int_t           HLT_HIMET220_v1;
   Int_t           HLT_HIMET220_v1_Prescl;
   Int_t           HLT_HIPhysics_v1;
   Int_t           HLT_HIPhysics_v1_Prescl;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v2;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v2_Prescl;
   Int_t           HLT_HIL1DoubleMuOpen_v2;
   Int_t           HLT_HIL1DoubleMuOpen_v2_Prescl;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v2;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v2_Prescl;
   Int_t           HLT_HIL2DoubleMu0_v2;
   Int_t           HLT_HIL2DoubleMu0_v2_Prescl;
   Int_t           HLT_HIL2DoubleMu3_v2;
   Int_t           HLT_HIL2DoubleMu3_v2_Prescl;
   Int_t           HLT_HIL2Mu15_v2;
   Int_t           HLT_HIL2Mu15_v2_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ_v2;
   Int_t           HLT_HIL2Mu3_NHitQ_v2_Prescl;
   Int_t           HLT_HIL2Mu3_v2;
   Int_t           HLT_HIL2Mu3_v2_Prescl;
   Int_t           HLT_HIL2Mu7_v2;
   Int_t           HLT_HIL2Mu7_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v2;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v2;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v2;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_v2;
   Int_t           HLT_HIL3DoubleMuOpen_v2_Prescl;
   Int_t           HLT_HIL3Mu3_v2;
   Int_t           HLT_HIL3Mu3_v2_Prescl;
   Int_t           HLT_HIDiJet55_v2;
   Int_t           HLT_HIDiJet55_v2_Prescl;
   Int_t           HLT_HIDoublePhoton10_v2;
   Int_t           HLT_HIDoublePhoton10_v2_Prescl;
   Int_t           HLT_HIDoublePhoton15_v2;
   Int_t           HLT_HIDoublePhoton15_v2_Prescl;
   Int_t           HLT_HIDoublePhoton20_v2;
   Int_t           HLT_HIDoublePhoton20_v2_Prescl;
   Int_t           HLT_HIFullTrack12_L1Central_v2;
   Int_t           HLT_HIFullTrack12_L1Central_v2_Prescl;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v2;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v2_Prescl;
   Int_t           HLT_HIFullTrack14_L1Central_v2;
   Int_t           HLT_HIFullTrack14_L1Central_v2_Prescl;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v2;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v2_Prescl;
   Int_t           HLT_HIFullTrack20_L1Central_v2;
   Int_t           HLT_HIFullTrack20_L1Central_v2_Prescl;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v2;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v2_Prescl;
   Int_t           HLT_HIFullTrack25_L1Central_v2;
   Int_t           HLT_HIFullTrack25_L1Central_v2_Prescl;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v2;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v2_Prescl;
   Int_t           HLT_HIJet55_v2;
   Int_t           HLT_HIJet55_v2_Prescl;
   Int_t           HLT_HIJet65_Jet55_v2;
   Int_t           HLT_HIJet65_Jet55_v2_Prescl;
   Int_t           HLT_HIJet65_v2;
   Int_t           HLT_HIJet65_v2_Prescl;
   Int_t           HLT_HIJet80_v2;
   Int_t           HLT_HIJet80_v2_Prescl;
   Int_t           HLT_HIJet95_v2;
   Int_t           HLT_HIJet95_v2_Prescl;
   Int_t           HLT_HIJetE30_NoBPTX_v2;
   Int_t           HLT_HIJetE30_NoBPTX_v2_Prescl;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v2;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl;
   Int_t           HLT_HIMET120_v2;
   Int_t           HLT_HIMET120_v2_Prescl;
   Int_t           HLT_HIMET200_v2;
   Int_t           HLT_HIMET200_v2_Prescl;
   Int_t           HLT_HIMET220_v2;
   Int_t           HLT_HIMET220_v2_Prescl;
   Int_t           HLT_HIPhoton10_Photon15_v2;
   Int_t           HLT_HIPhoton10_Photon15_v2_Prescl;
   Int_t           HLT_HIPhoton15_Photon20_v2;
   Int_t           HLT_HIPhoton15_Photon20_v2_Prescl;
   Int_t           HLT_HISinglePhoton15_v2;
   Int_t           HLT_HISinglePhoton15_v2_Prescl;
   Int_t           HLT_HISinglePhoton20_v2;
   Int_t           HLT_HISinglePhoton20_v2_Prescl;
   Int_t           HLT_HISinglePhoton30_v2;
   Int_t           HLT_HISinglePhoton30_v2_Prescl;
   Int_t           HLT_HISinglePhoton40_v2;
   Int_t           HLT_HISinglePhoton40_v2_Prescl;
   Int_t           HLT_HIBptxXOR_v2;
   Int_t           HLT_HIBptxXOR_v2_Prescl;
   Int_t           HLT_HICentral10_v2;
   Int_t           HLT_HICentral10_v2_Prescl;
   Int_t           HLT_HICentralityVeto_v2;
   Int_t           HLT_HICentralityVeto_v2_Prescl;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v2;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl;
   Int_t           HLT_HIMinBiasBSC_OR_v2;
   Int_t           HLT_HIMinBiasBSC_OR_v2_Prescl;
   Int_t           HLT_HIMinBiasBSC_v2;
   Int_t           HLT_HIMinBiasBSC_v2_Prescl;
   Int_t           HLT_HIMinBiasHF_v2;
   Int_t           HLT_HIMinBiasHF_v2_Prescl;
   Int_t           HLT_HIMinBiasHfOrBSC_v2;
   Int_t           HLT_HIMinBiasHfOrBSC_v2_Prescl;
   Int_t           HLT_HIMinBiasHf_OR_v2;
   Int_t           HLT_HIMinBiasHf_OR_v2_Prescl;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v2;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v2;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_v2;
   Int_t           HLT_HIMinBiasZDC_Calo_v2_Prescl;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIPhysics_v2;
   Int_t           HLT_HIPhysics_v2_Prescl;
   Int_t           HLT_HIRandom_v2;
   Int_t           HLT_HIRandom_v2_Prescl;
   Int_t           HLT_HIUCC010_v2;
   Int_t           HLT_HIUCC010_v2_Prescl;
   Int_t           HLT_HIUCC015_v2;
   Int_t           HLT_HIUCC015_v2_Prescl;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v2;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIZeroBiasXOR_v2;
   Int_t           HLT_HIZeroBiasXOR_v2_Prescl;
   Int_t           HLT_HIZeroBias_v2;
   Int_t           HLT_HIZeroBias_v2_Prescl;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v3;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v3_Prescl;
   Int_t           HLT_HIL1DoubleMuOpen_v3;
   Int_t           HLT_HIL1DoubleMuOpen_v3_Prescl;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v3;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v3_Prescl;
   Int_t           HLT_HIL2DoubleMu0_v3;
   Int_t           HLT_HIL2DoubleMu0_v3_Prescl;
   Int_t           HLT_HIL2DoubleMu3_v3;
   Int_t           HLT_HIL2DoubleMu3_v3_Prescl;
   Int_t           HLT_HIL2Mu15_v3;
   Int_t           HLT_HIL2Mu15_v3_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ_v3;
   Int_t           HLT_HIL2Mu3_NHitQ_v3_Prescl;
   Int_t           HLT_HIL2Mu3_v3;
   Int_t           HLT_HIL2Mu3_v3_Prescl;
   Int_t           HLT_HIL2Mu7_v3;
   Int_t           HLT_HIL2Mu7_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v3;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v3;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v3;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_v3;
   Int_t           HLT_HIL3DoubleMuOpen_v3_Prescl;
   Int_t           HLT_HIL3Mu3_v3;
   Int_t           HLT_HIL3Mu3_v3_Prescl;
   Int_t           HLT_HIDiJet55_v3;
   Int_t           HLT_HIDiJet55_v3_Prescl;
   Int_t           HLT_HIDoublePhoton10_v3;
   Int_t           HLT_HIDoublePhoton10_v3_Prescl;
   Int_t           HLT_HIDoublePhoton15_v3;
   Int_t           HLT_HIDoublePhoton15_v3_Prescl;
   Int_t           HLT_HIDoublePhoton20_v3;
   Int_t           HLT_HIDoublePhoton20_v3_Prescl;
   Int_t           HLT_HIFullTrack12_L1Central_v3;
   Int_t           HLT_HIFullTrack12_L1Central_v3_Prescl;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v3;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v3_Prescl;
   Int_t           HLT_HIFullTrack14_L1Central_v3;
   Int_t           HLT_HIFullTrack14_L1Central_v3_Prescl;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v3;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v3_Prescl;
   Int_t           HLT_HIFullTrack20_L1Central_v3;
   Int_t           HLT_HIFullTrack20_L1Central_v3_Prescl;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v3;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v3_Prescl;
   Int_t           HLT_HIFullTrack25_L1Central_v3;
   Int_t           HLT_HIFullTrack25_L1Central_v3_Prescl;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v3;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v3_Prescl;
   Int_t           HLT_HIJet55_v3;
   Int_t           HLT_HIJet55_v3_Prescl;
   Int_t           HLT_HIJet65_Jet55_v3;
   Int_t           HLT_HIJet65_Jet55_v3_Prescl;
   Int_t           HLT_HIJet65_v3;
   Int_t           HLT_HIJet65_v3_Prescl;
   Int_t           HLT_HIJet80_v3;
   Int_t           HLT_HIJet80_v3_Prescl;
   Int_t           HLT_HIJet95_v3;
   Int_t           HLT_HIJet95_v3_Prescl;
   Int_t           HLT_HIJetE30_NoBPTX_v3;
   Int_t           HLT_HIJetE30_NoBPTX_v3_Prescl;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v3;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl;
   Int_t           HLT_HIMET120_v3;
   Int_t           HLT_HIMET120_v3_Prescl;
   Int_t           HLT_HIMET200_v3;
   Int_t           HLT_HIMET200_v3_Prescl;
   Int_t           HLT_HIMET220_v3;
   Int_t           HLT_HIMET220_v3_Prescl;
   Int_t           HLT_HIPhoton10_Photon15_v3;
   Int_t           HLT_HIPhoton10_Photon15_v3_Prescl;
   Int_t           HLT_HIPhoton15_Photon20_v3;
   Int_t           HLT_HIPhoton15_Photon20_v3_Prescl;
   Int_t           HLT_HISinglePhoton15_v3;
   Int_t           HLT_HISinglePhoton15_v3_Prescl;
   Int_t           HLT_HISinglePhoton20_v3;
   Int_t           HLT_HISinglePhoton20_v3_Prescl;
   Int_t           HLT_HISinglePhoton30_v3;
   Int_t           HLT_HISinglePhoton30_v3_Prescl;
   Int_t           HLT_HISinglePhoton40_v3;
   Int_t           HLT_HISinglePhoton40_v3_Prescl;
   Int_t           HLT_HIBptxXOR_v3;
   Int_t           HLT_HIBptxXOR_v3_Prescl;
   Int_t           HLT_HICentral10_v3;
   Int_t           HLT_HICentral10_v3_Prescl;
   Int_t           HLT_HICentralityVeto_v3;
   Int_t           HLT_HICentralityVeto_v3_Prescl;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v3;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl;
   Int_t           HLT_HIMinBiasBSC_OR_v3;
   Int_t           HLT_HIMinBiasBSC_OR_v3_Prescl;
   Int_t           HLT_HIMinBiasBSC_v3;
   Int_t           HLT_HIMinBiasBSC_v3_Prescl;
   Int_t           HLT_HIMinBiasHF_v3;
   Int_t           HLT_HIMinBiasHF_v3_Prescl;
   Int_t           HLT_HIMinBiasHfOrBSC_v3;
   Int_t           HLT_HIMinBiasHfOrBSC_v3_Prescl;
   Int_t           HLT_HIMinBiasHf_OR_v3;
   Int_t           HLT_HIMinBiasHf_OR_v3_Prescl;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v3;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v3;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_v3;
   Int_t           HLT_HIMinBiasZDC_Calo_v3_Prescl;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIPhysics_v3;
   Int_t           HLT_HIPhysics_v3_Prescl;
   Int_t           HLT_HIRandom_v3;
   Int_t           HLT_HIRandom_v3_Prescl;
   Int_t           HLT_HIUCC010_v3;
   Int_t           HLT_HIUCC010_v3_Prescl;
   Int_t           HLT_HIUCC015_v3;
   Int_t           HLT_HIUCC015_v3_Prescl;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v3;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIZeroBiasXOR_v3;
   Int_t           HLT_HIZeroBiasXOR_v3_Prescl;
   Int_t           HLT_HIZeroBias_v3;
   Int_t           HLT_HIZeroBias_v3_Prescl;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v4;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v4_Prescl;
   Int_t           HLT_HIL1DoubleMuOpen_v4;
   Int_t           HLT_HIL1DoubleMuOpen_v4_Prescl;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v4;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v4_Prescl;
   Int_t           HLT_HIL2DoubleMu0_v4;
   Int_t           HLT_HIL2DoubleMu0_v4_Prescl;
   Int_t           HLT_HIL2DoubleMu3_v4;
   Int_t           HLT_HIL2DoubleMu3_v4_Prescl;
   Int_t           HLT_HIL2Mu15_v4;
   Int_t           HLT_HIL2Mu15_v4_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ_v4;
   Int_t           HLT_HIL2Mu3_NHitQ_v4_Prescl;
   Int_t           HLT_HIL2Mu3_v4;
   Int_t           HLT_HIL2Mu3_v4_Prescl;
   Int_t           HLT_HIL2Mu7_v4;
   Int_t           HLT_HIL2Mu7_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v4;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v4;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v4;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_v4;
   Int_t           HLT_HIL3DoubleMuOpen_v4_Prescl;
   Int_t           HLT_HIL3Mu3_v4;
   Int_t           HLT_HIL3Mu3_v4_Prescl;
   Int_t           HLT_HIDiJet55_v4;
   Int_t           HLT_HIDiJet55_v4_Prescl;
   Int_t           HLT_HIDoublePhoton10_v4;
   Int_t           HLT_HIDoublePhoton10_v4_Prescl;
   Int_t           HLT_HIDoublePhoton15_v4;
   Int_t           HLT_HIDoublePhoton15_v4_Prescl;
   Int_t           HLT_HIDoublePhoton20_v4;
   Int_t           HLT_HIDoublePhoton20_v4_Prescl;
   Int_t           HLT_HIFullTrack12_L1Central_v4;
   Int_t           HLT_HIFullTrack12_L1Central_v4_Prescl;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v4;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v4_Prescl;
   Int_t           HLT_HIFullTrack14_L1Central_v4;
   Int_t           HLT_HIFullTrack14_L1Central_v4_Prescl;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v4;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v4_Prescl;
   Int_t           HLT_HIFullTrack20_L1Central_v4;
   Int_t           HLT_HIFullTrack20_L1Central_v4_Prescl;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v4;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v4_Prescl;
   Int_t           HLT_HIFullTrack25_L1Central_v4;
   Int_t           HLT_HIFullTrack25_L1Central_v4_Prescl;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v4;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v4_Prescl;
   Int_t           HLT_HIJet55_v4;
   Int_t           HLT_HIJet55_v4_Prescl;
   Int_t           HLT_HIJet65_Jet55_v4;
   Int_t           HLT_HIJet65_Jet55_v4_Prescl;
   Int_t           HLT_HIJet65_v4;
   Int_t           HLT_HIJet65_v4_Prescl;
   Int_t           HLT_HIJet80_v4;
   Int_t           HLT_HIJet80_v4_Prescl;
   Int_t           HLT_HIJet95_v4;
   Int_t           HLT_HIJet95_v4_Prescl;
   Int_t           HLT_HIJetE30_NoBPTX_v4;
   Int_t           HLT_HIJetE30_NoBPTX_v4_Prescl;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v4;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl;
   Int_t           HLT_HIMET120_v4;
   Int_t           HLT_HIMET120_v4_Prescl;
   Int_t           HLT_HIMET200_v4;
   Int_t           HLT_HIMET200_v4_Prescl;
   Int_t           HLT_HIMET220_v4;
   Int_t           HLT_HIMET220_v4_Prescl;
   Int_t           HLT_HIPhoton10_Photon15_v4;
   Int_t           HLT_HIPhoton10_Photon15_v4_Prescl;
   Int_t           HLT_HIPhoton15_Photon20_v4;
   Int_t           HLT_HIPhoton15_Photon20_v4_Prescl;
   Int_t           HLT_HISinglePhoton15_v4;
   Int_t           HLT_HISinglePhoton15_v4_Prescl;
   Int_t           HLT_HISinglePhoton20_v4;
   Int_t           HLT_HISinglePhoton20_v4_Prescl;
   Int_t           HLT_HISinglePhoton30_v4;
   Int_t           HLT_HISinglePhoton30_v4_Prescl;
   Int_t           HLT_HISinglePhoton40_v4;
   Int_t           HLT_HISinglePhoton40_v4_Prescl;
   Int_t           HLT_HIBptxXOR_v4;
   Int_t           HLT_HIBptxXOR_v4_Prescl;
   Int_t           HLT_HICentral10_v4;
   Int_t           HLT_HICentral10_v4_Prescl;
   Int_t           HLT_HICentralityVeto_v4;
   Int_t           HLT_HICentralityVeto_v4_Prescl;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v4;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl;
   Int_t           HLT_HIMinBiasBSC_OR_v4;
   Int_t           HLT_HIMinBiasBSC_OR_v4_Prescl;
   Int_t           HLT_HIMinBiasBSC_v4;
   Int_t           HLT_HIMinBiasBSC_v4_Prescl;
   Int_t           HLT_HIMinBiasHF_v4;
   Int_t           HLT_HIMinBiasHF_v4_Prescl;
   Int_t           HLT_HIMinBiasHfOrBSC_v4;
   Int_t           HLT_HIMinBiasHfOrBSC_v4_Prescl;
   Int_t           HLT_HIMinBiasHf_OR_v4;
   Int_t           HLT_HIMinBiasHf_OR_v4_Prescl;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v4;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v4;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_v4;
   Int_t           HLT_HIMinBiasZDC_Calo_v4_Prescl;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIPhysics_v4;
   Int_t           HLT_HIPhysics_v4_Prescl;
   Int_t           HLT_HIRandom_v4;
   Int_t           HLT_HIRandom_v4_Prescl;
   Int_t           HLT_HIUCC010_v4;
   Int_t           HLT_HIUCC010_v4_Prescl;
   Int_t           HLT_HIUCC015_v4;
   Int_t           HLT_HIUCC015_v4_Prescl;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v4;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIZeroBiasXOR_v4;
   Int_t           HLT_HIZeroBiasXOR_v4_Prescl;
   Int_t           HLT_HIZeroBias_v4;
   Int_t           HLT_HIZeroBias_v4_Prescl;
   Int_t           L1_BptxMinus;
   Int_t           L1_BptxMinus_Prescl;
   Int_t           L1_BptxMinus_NotBptxPlus;
   Int_t           L1_BptxMinus_NotBptxPlus_Prescl;
   Int_t           L1_BptxPlus;
   Int_t           L1_BptxPlus_Prescl;
   Int_t           L1_BptxPlusANDMinus;
   Int_t           L1_BptxPlusANDMinus_Prescl;
   Int_t           L1_BptxPlusORMinus;
   Int_t           L1_BptxPlusORMinus_Prescl;
   Int_t           L1_BptxPlusORMinus_instance1;
   Int_t           L1_BptxPlusORMinus_instance1_Prescl;
   Int_t           L1_BptxPlus_NotBptxMinus;
   Int_t           L1_BptxPlus_NotBptxMinus_Prescl;
   Int_t           L1_BptxXOR;
   Int_t           L1_BptxXOR_Prescl;
   Int_t           L1_BptxXOR_BscMinBiasOR;
   Int_t           L1_BptxXOR_BscMinBiasOR_Prescl;
   Int_t           L1_Bsc2Minus_BptxMinus;
   Int_t           L1_Bsc2Minus_BptxMinus_Prescl;
   Int_t           L1_Bsc2Plus_BptxPlus;
   Int_t           L1_Bsc2Plus_BptxPlus_Prescl;
   Int_t           L1_BscMinBiasInnerThreshold1;
   Int_t           L1_BscMinBiasInnerThreshold1_Prescl;
   Int_t           L1_BscMinBiasInnerThreshold1_BptxAND;
   Int_t           L1_BscMinBiasInnerThreshold1_BptxAND_Prescl;
   Int_t           L1_BscMinBiasInnerThreshold2_BptxAND;
   Int_t           L1_BscMinBiasInnerThreshold2_BptxAND_Prescl;
   Int_t           L1_BscMinBiasOR_BptxAND;
   Int_t           L1_BscMinBiasOR_BptxAND_Prescl;
   Int_t           L1_BscMinBiasOR_BptxPlusORMinus;
   Int_t           L1_BscMinBiasOR_BptxPlusORMinus_Prescl;
   Int_t           L1_BscMinBiasThreshold1;
   Int_t           L1_BscMinBiasThreshold1_Prescl;
   Int_t           L1_BscMinBiasThreshold1_BptxAND;
   Int_t           L1_BscMinBiasThreshold1_BptxAND_Prescl;
   Int_t           L1_BscMinBiasThreshold2_BptxAND;
   Int_t           L1_BscMinBiasThreshold2_BptxAND_Prescl;
   Int_t           L1_DoubleEG10;
   Int_t           L1_DoubleEG10_Prescl;
   Int_t           L1_DoubleEG3;
   Int_t           L1_DoubleEG3_Prescl;
   Int_t           L1_DoubleEG5_BptxAND;
   Int_t           L1_DoubleEG5_BptxAND_Prescl;
   Int_t           L1_DoubleEG8_BptxAND;
   Int_t           L1_DoubleEG8_BptxAND_Prescl;
   Int_t           L1_DoubleForJet32_EtaOpp;
   Int_t           L1_DoubleForJet32_EtaOpp_Prescl;
   Int_t           L1_DoubleForJet44_EtaOpp;
   Int_t           L1_DoubleForJet44_EtaOpp_Prescl;
   Int_t           L1_DoubleHfBitCountsRing1_P1N1;
   Int_t           L1_DoubleHfBitCountsRing1_P1N1_Prescl;
   Int_t           L1_DoubleHfBitCountsRing2_P1N1;
   Int_t           L1_DoubleHfBitCountsRing2_P1N1_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing1_P200N200;
   Int_t           L1_DoubleHfRingEtSumsRing1_P200N200_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing1_P4N4;
   Int_t           L1_DoubleHfRingEtSumsRing1_P4N4_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing2_P200N200;
   Int_t           L1_DoubleHfRingEtSumsRing2_P200N200_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing2_P4N4;
   Int_t           L1_DoubleHfRingEtSumsRing2_P4N4_Prescl;
   Int_t           L1_DoubleMu3_BptxAND;
   Int_t           L1_DoubleMu3_BptxAND_Prescl;
   Int_t           L1_DoubleMuOpen_BptxAND;
   Int_t           L1_DoubleMuOpen_BptxAND_Prescl;
   Int_t           L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND;
   Int_t           L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;
   Int_t           L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;
   Int_t           L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;
   Int_t           L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND;
   Int_t           L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;
   Int_t           L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;
   Int_t           L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;
   Int_t           L1_ETM30_BptxAND;
   Int_t           L1_ETM30_BptxAND_Prescl;
   Int_t           L1_ETM50_BptxAND;
   Int_t           L1_ETM50_BptxAND_Prescl;
   Int_t           L1_ETT100_BptxAND;
   Int_t           L1_ETT100_BptxAND_Prescl;
   Int_t           L1_ETT140_BptxAND;
   Int_t           L1_ETT140_BptxAND_Prescl;
   Int_t           L1_ETT2000;
   Int_t           L1_ETT2000_Prescl;
   Int_t           L1_ETT220_BptxAND;
   Int_t           L1_ETT220_BptxAND_Prescl;
   Int_t           L1_ETT70_BptxAND;
   Int_t           L1_ETT70_BptxAND_Prescl;
   Int_t           L1_ETT800;
   Int_t           L1_ETT800_Prescl;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl;
   Int_t           L1_HcalHfCoincidencePm;
   Int_t           L1_HcalHfCoincidencePm_Prescl;
   Int_t           L1_HcalHfCoincidencePm_BptxAND;
   Int_t           L1_HcalHfCoincidencePm_BptxAND_Prescl;
   Int_t           L1_HcalHfMmOrPpOrPm;
   Int_t           L1_HcalHfMmOrPpOrPm_Prescl;
   Int_t           L1_HcalHfMmOrPpOrPm_BptxAND;
   Int_t           L1_HcalHfMmOrPpOrPm_BptxAND_Prescl;
   Int_t           L1_HcalHfMmpOrMpp;
   Int_t           L1_HcalHfMmpOrMpp_Prescl;
   Int_t           L1_HcalHfMmpOrMpp_BptxAND;
   Int_t           L1_HcalHfMmpOrMpp_BptxAND_Prescl;
   Int_t           L1_HcalHoTotalOR;
   Int_t           L1_HcalHoTotalOR_Prescl;
   Int_t           L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND;
   Int_t           L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;
   Int_t           L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;
   Int_t           L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasOR;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasThreshold1;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasThreshold2;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl;
   Int_t           L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1;
   Int_t           L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl;
   Int_t           L1_NotBsc2_BptxAND_HcalHfCoincidencePm;
   Int_t           L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl;
   Int_t           L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm;
   Int_t           L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl;
   Int_t           L1_NotBsc2_BptxAND_HcalHfMmpOrMpp;
   Int_t           L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasInnerThreshold1;
   Int_t           L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasInnerThreshold2;
   Int_t           L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasOR;
   Int_t           L1_NotBsc2_BscMinBiasOR_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasThreshold1;
   Int_t           L1_NotBsc2_BscMinBiasThreshold1_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasThreshold2;
   Int_t           L1_NotBsc2_BscMinBiasThreshold2_Prescl;
   Int_t           L1_NotBsc2_HcalHfCoincidencePm;
   Int_t           L1_NotBsc2_HcalHfCoincidencePm_Prescl;
   Int_t           L1_NotBsc2_HcalHfMmOrPpOrPm;
   Int_t           L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl;
   Int_t           L1_NotBsc2_HcalHfMmpOrMpp;
   Int_t           L1_NotBsc2_HcalHfMmpOrMpp_Prescl;
   Int_t           L1_SingleEG12;
   Int_t           L1_SingleEG12_Prescl;
   Int_t           L1_SingleEG15;
   Int_t           L1_SingleEG15_Prescl;
   Int_t           L1_SingleEG18;
   Int_t           L1_SingleEG18_Prescl;
   Int_t           L1_SingleEG20;
   Int_t           L1_SingleEG20_Prescl;
   Int_t           L1_SingleEG22;
   Int_t           L1_SingleEG22_Prescl;
   Int_t           L1_SingleEG2_BptxAND;
   Int_t           L1_SingleEG2_BptxAND_Prescl;
   Int_t           L1_SingleEG30;
   Int_t           L1_SingleEG30_Prescl;
   Int_t           L1_SingleEG5_BptxAND;
   Int_t           L1_SingleEG5_BptxAND_Prescl;
   Int_t           L1_SingleEG8_BptxAND;
   Int_t           L1_SingleEG8_BptxAND_Prescl;
   Int_t           L1_SingleJet128;
   Int_t           L1_SingleJet128_Prescl;
   Int_t           L1_SingleJet128_BptxAND;
   Int_t           L1_SingleJet128_BptxAND_Prescl;
   Int_t           L1_SingleJet128_Central_BptxAND;
   Int_t           L1_SingleJet128_Central_BptxAND_Prescl;
   Int_t           L1_SingleJet16;
   Int_t           L1_SingleJet16_Prescl;
   Int_t           L1_SingleJet16_Central_NotETT140_BptxAND;
   Int_t           L1_SingleJet16_Central_NotETT140_BptxAND_Prescl;
   Int_t           L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo;
   Int_t           L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl;
   Int_t           L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo;
   Int_t           L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl;
   Int_t           L1_SingleJet36_BptxAND;
   Int_t           L1_SingleJet36_BptxAND_Prescl;
   Int_t           L1_SingleJet52_BptxAND;
   Int_t           L1_SingleJet52_BptxAND_Prescl;
   Int_t           L1_SingleJet52_Central_BptxAND;
   Int_t           L1_SingleJet52_Central_BptxAND_Prescl;
   Int_t           L1_SingleJet68_BptxAND;
   Int_t           L1_SingleJet68_BptxAND_Prescl;
   Int_t           L1_SingleJet80_BptxAND;
   Int_t           L1_SingleJet80_BptxAND_Prescl;
   Int_t           L1_SingleJet80_Central;
   Int_t           L1_SingleJet80_Central_Prescl;
   Int_t           L1_SingleJet80_Central_BptxAND;
   Int_t           L1_SingleJet80_Central_BptxAND_Prescl;
   Int_t           L1_SingleJet92;
   Int_t           L1_SingleJet92_Prescl;
   Int_t           L1_SingleJet92_BptxAND;
   Int_t           L1_SingleJet92_BptxAND_Prescl;
   Int_t           L1_SingleJet92_Central;
   Int_t           L1_SingleJet92_Central_Prescl;
   Int_t           L1_SingleJet92_Central_BptxAND;
   Int_t           L1_SingleJet92_Central_BptxAND_Prescl;
   Int_t           L1_SingleMu0_BptxAND;
   Int_t           L1_SingleMu0_BptxAND_Prescl;
   Int_t           L1_SingleMu10;
   Int_t           L1_SingleMu10_Prescl;
   Int_t           L1_SingleMu14;
   Int_t           L1_SingleMu14_Prescl;
   Int_t           L1_SingleMu20;
   Int_t           L1_SingleMu20_Prescl;
   Int_t           L1_SingleMu3_BptxAND;
   Int_t           L1_SingleMu3_BptxAND_Prescl;
   Int_t           L1_SingleMu5_BptxAND;
   Int_t           L1_SingleMu5_BptxAND_Prescl;
   Int_t           L1_SingleMu7;
   Int_t           L1_SingleMu7_Prescl;
   Int_t           L1_SingleMuBeamHalo;
   Int_t           L1_SingleMuBeamHalo_Prescl;
   Int_t           L1_SingleMuOpen_BptxAND;
   Int_t           L1_SingleMuOpen_BptxAND_Prescl;
   Int_t           L1_ZdcCaloMinus;
   Int_t           L1_ZdcCaloMinus_Prescl;
   Int_t           L1_ZdcCaloMinus_BptxAND;
   Int_t           L1_ZdcCaloMinus_BptxAND_Prescl;
   Int_t           L1_ZdcCaloPlus;
   Int_t           L1_ZdcCaloPlus_Prescl;
   Int_t           L1_ZdcCaloPlus_BptxAND;
   Int_t           L1_ZdcCaloPlus_BptxAND_Prescl;
   Int_t           L1_ZdcCaloPlus_ZdcCaloMinus;
   Int_t           L1_ZdcCaloPlus_ZdcCaloMinus_Prescl;
   Int_t           L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND;
   Int_t           L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl;
   Int_t           L1_ZeroBias;
   Int_t           L1_ZeroBias_Prescl;
   Int_t           L1Tech_BPTX_PreBPTX_v0;
   Int_t           L1Tech_BPTX_PreBPTX_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_v0;
   Int_t           L1Tech_BPTX_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_v0;
   Int_t           L1Tech_BPTX_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_quiet_v0;
   Int_t           L1Tech_BPTX_quiet_v0_Prescl;
   Int_t           L1Tech_BSC_HighMultiplicity_v0;
   Int_t           L1Tech_BSC_HighMultiplicity_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam1_inner_v0;
   Int_t           L1Tech_BSC_halo_beam1_inner_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam1_outer_v0;
   Int_t           L1Tech_BSC_halo_beam1_outer_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam2_inner_v0;
   Int_t           L1Tech_BSC_halo_beam2_inner_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam2_outer_v0;
   Int_t           L1Tech_BSC_halo_beam2_outer_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_OR_v0;
   Int_t           L1Tech_BSC_minBias_OR_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_inner_threshold1_v0;
   Int_t           L1Tech_BSC_minBias_inner_threshold1_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_inner_threshold2_v0;
   Int_t           L1Tech_BSC_minBias_inner_threshold2_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_threshold1_v0;
   Int_t           L1Tech_BSC_minBias_threshold1_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_threshold2_v0;
   Int_t           L1Tech_BSC_minBias_threshold2_v0_Prescl;
   Int_t           L1Tech_BSC_splash_beam1_v0;
   Int_t           L1Tech_BSC_splash_beam1_v0_Prescl;
   Int_t           L1Tech_BSC_splash_beam2_v0;
   Int_t           L1Tech_BSC_splash_beam2_v0_Prescl;
   Int_t           L1Tech_CASTOR_HaloMuon_v0;
   Int_t           L1Tech_CASTOR_HaloMuon_v0_Prescl;
   Int_t           L1Tech_DT_GlobalOR_v0;
   Int_t           L1Tech_DT_GlobalOR_v0_Prescl;
   Int_t           L1Tech_FSC_St1Sect45_down_v0;
   Int_t           L1Tech_FSC_St1Sect45_down_v0_Prescl;
   Int_t           L1Tech_FSC_St1Sect45_upp_v0;
   Int_t           L1Tech_FSC_St1Sect45_upp_v0_Prescl;
   Int_t           L1Tech_FSC_St1Sect56_down_v0;
   Int_t           L1Tech_FSC_St1Sect56_down_v0_Prescl;
   Int_t           L1Tech_FSC_St1Sect56_upp_v0;
   Int_t           L1Tech_FSC_St1Sect56_upp_v0_Prescl;
   Int_t           L1Tech_FSC_St2Sect45_down_v0;
   Int_t           L1Tech_FSC_St2Sect45_down_v0_Prescl;
   Int_t           L1Tech_FSC_St2Sect45_upp_v0;
   Int_t           L1Tech_FSC_St2Sect45_upp_v0_Prescl;
   Int_t           L1Tech_FSC_St2Sect56_down_v0;
   Int_t           L1Tech_FSC_St2Sect56_down_v0_Prescl;
   Int_t           L1Tech_FSC_St2Sect56_upp_v0;
   Int_t           L1Tech_FSC_St2Sect56_upp_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_downLeft_v0;
   Int_t           L1Tech_FSC_St3Sect45_downLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_downRight_v0;
   Int_t           L1Tech_FSC_St3Sect45_downRight_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_uppLeft_v0;
   Int_t           L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_uppRight_v0;
   Int_t           L1Tech_FSC_St3Sect45_uppRight_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect56_downLeft_v0;
   Int_t           L1Tech_FSC_St3Sect56_downLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect56_downRight_v0;
   Int_t           L1Tech_FSC_St3Sect56_downRight_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect56_uppLeft_v0;
   Int_t           L1Tech_FSC_St3Sect56_uppLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect56_uppRight_v0;
   Int_t           L1Tech_FSC_St3Sect56_uppRight_v0_Prescl;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0_Prescl;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;
   Int_t           L1Tech_HCAL_HF_single_channel_v0;
   Int_t           L1Tech_HCAL_HF_single_channel_v0_Prescl;
   Int_t           L1Tech_HCAL_HO_totalOR_v0;
   Int_t           L1Tech_HCAL_HO_totalOR_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RB0_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBminus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBminus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBst1_collisions_v0;
   Int_t           L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;

   // List of branches
   TBranch        *b_NL1IsolEm;   //!
   TBranch        *b_L1IsolEmEt;   //!
   TBranch        *b_L1IsolEmE;   //!
   TBranch        *b_L1IsolEmEta;   //!
   TBranch        *b_L1IsolEmPhi;   //!
   TBranch        *b_NL1NIsolEm;   //!
   TBranch        *b_L1NIsolEmEt;   //!
   TBranch        *b_L1NIsolEmE;   //!
   TBranch        *b_L1NIsolEmEta;   //!
   TBranch        *b_L1NIsolEmPhi;   //!
   TBranch        *b_NL1Mu;   //!
   TBranch        *b_L1MuPt;   //!
   TBranch        *b_L1MuE;   //!
   TBranch        *b_L1MuEta;   //!
   TBranch        *b_L1MuPhi;   //!
   TBranch        *b_L1MuIsol;   //!
   TBranch        *b_L1MuMip;   //!
   TBranch        *b_L1MuFor;   //!
   TBranch        *b_L1MuRPC;   //!
   TBranch        *b_L1MuQal;   //!
   TBranch        *b_L1MuChg;   //!
   TBranch        *b_NL1CenJet;   //!
   TBranch        *b_L1CenJetEt;   //!
   TBranch        *b_L1CenJetE;   //!
   TBranch        *b_L1CenJetEta;   //!
   TBranch        *b_L1CenJetPhi;   //!
   TBranch        *b_NL1ForJet;   //!
   TBranch        *b_L1ForJetEt;   //!
   TBranch        *b_L1ForJetE;   //!
   TBranch        *b_L1ForJetEta;   //!
   TBranch        *b_L1ForJetPhi;   //!
   TBranch        *b_NL1Tau;   //!
   TBranch        *b_L1TauEt;   //!
   TBranch        *b_L1TauE;   //!
   TBranch        *b_L1TauEta;   //!
   TBranch        *b_L1TauPhi;   //!
   TBranch        *b_L1Met;   //!
   TBranch        *b_L1MetPhi;   //!
   TBranch        *b_L1EtTot;   //!
   TBranch        *b_L1Mht;   //!
   TBranch        *b_L1MhtPhi;   //!
   TBranch        *b_L1EtHad;   //!
   TBranch        *b_L1HfRing1EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing2EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing1EtSumNegativeEta;   //!
   TBranch        *b_L1HfRing2EtSumNegativeEta;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing1;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing1;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing2;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing2;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Bx;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_AvgInstDelLumi;   //!
   TBranch        *b_raw2digi_step;   //!
   TBranch        *b_raw2digi_step_Prescl;   //!
   TBranch        *b_reconstruction_step;   //!
   TBranch        *b_reconstruction_step_Prescl;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLTriggerFirstPath_Prescl;   //!
   TBranch        *b_HLT_HIDTCalibration_v1;   //!
   TBranch        *b_HLT_HIDTCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_HIEcalCalibration_v1;   //!
   TBranch        *b_HLT_HIEcalCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_HIHcalCalibration_v1;   //!
   TBranch        *b_HLT_HIHcalCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_HIZeroBias_v1;   //!
   TBranch        *b_HLT_HIZeroBias_v1_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v1;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v1_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v1;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v1;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHF_v1;   //!
   TBranch        *b_HLT_HIMinBiasHF_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v1;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v1;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v1;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIBptxXOR_v1;   //!
   TBranch        *b_HLT_HIBptxXOR_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v1;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v1;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v1;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_v1;   //!
   TBranch        *b_HLT_HIL2Mu3_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v1;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_v1;   //!
   TBranch        *b_HLT_HIL2Mu7_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v1;   //!
   TBranch        *b_HLT_HIL2Mu15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_v1;   //!
   TBranch        *b_HLT_HIL3Mu3_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_v1;   //!
   TBranch        *b_HLT_HISinglePhoton15_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_v1;   //!
   TBranch        *b_HLT_HISinglePhoton20_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_v1;   //!
   TBranch        *b_HLT_HISinglePhoton30_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_v1;   //!
   TBranch        *b_HLT_HISinglePhoton40_v1_Prescl;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v1;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v1_Prescl;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v1;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v1_Prescl;   //!
   TBranch        *b_HLT_HIJet55_v1;   //!
   TBranch        *b_HLT_HIJet55_v1_Prescl;   //!
   TBranch        *b_HLT_HIJet65_v1;   //!
   TBranch        *b_HLT_HIJet65_v1_Prescl;   //!
   TBranch        *b_HLT_Jet60_v1;   //!
   TBranch        *b_HLT_PAJet20_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet40_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet60_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet80_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet100_NoJetID_v1;   //!
   TBranch        *b_HLT_PAZeroBias_v1;   //!
   TBranch        *b_HLT_PAZeroBias_v1_Prescl;   //!
   TBranch        *b_HLT_PAZeroBiasPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_PAFullTrack12_v2;   //!
   TBranch        *b_HLT_PAFullTrack20_v2;   //!
   TBranch        *b_HLT_PAFullTrack30_v2;   //!
   TBranch        *b_HLT_HIJet80_v1;   //!
   TBranch        *b_HLT_HIJet80_v1_Prescl;   //!
   TBranch        *b_HLT_HIJet95_v1;   //!
   TBranch        *b_HLT_HIJet95_v1_Prescl;   //!
   TBranch        *b_HLT_HIDiJet55_v1;   //!
   TBranch        *b_HLT_HIDiJet55_v1_Prescl;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v1;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v1_Prescl;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v1;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v1_Prescl;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v1;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl;   //!
   TBranch        *b_HLT_HIActivityHF_Coincidence3_v1;   //!
   TBranch        *b_HLT_HIActivityHF_Coincidence3_v1_Prescl;   //!
   TBranch        *b_HLT_HIActivityHF_Single3_v1;   //!
   TBranch        *b_HLT_HIActivityHF_Single3_v1_Prescl;   //!
   TBranch        *b_HLT_HIClusterVertexCompatibility_v1;   //!
   TBranch        *b_HLT_HIClusterVertexCompatibility_v1_Prescl;   //!
   TBranch        *b_HLT_HICentralityVeto_v1;   //!
   TBranch        *b_HLT_HICentralityVeto_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v1;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v1;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v1;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v1;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v1;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v1;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v1;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v1;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v1_Prescl;   //!
   TBranch        *b_HLT_HIRandom_v1;   //!
   TBranch        *b_HLT_HIRandom_v1_Prescl;   //!
   TBranch        *b_HLT_HIUCC010_v1;   //!
   TBranch        *b_HLT_HIUCC010_v1_Prescl;   //!
   TBranch        *b_HLT_HIUCC015_v1;   //!
   TBranch        *b_HLT_HIUCC015_v1_Prescl;   //!
   TBranch        *b_HLT_HICentral10_v1;   //!
   TBranch        *b_HLT_HICentral10_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFinalPath_Prescl;   //!
   TBranch        *b_HLT_HIMET120_v1;   //!
   TBranch        *b_HLT_HIMET120_v1_Prescl;   //!
   TBranch        *b_HLT_HIMET200_v1;   //!
   TBranch        *b_HLT_HIMET200_v1_Prescl;   //!
   TBranch        *b_HLT_HIMET220_v1;   //!
   TBranch        *b_HLT_HIMET220_v1_Prescl;   //!
   TBranch        *b_HLT_HIPhysics_v1;   //!
   TBranch        *b_HLT_HIPhysics_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v2;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v2_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v2;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v2;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v2;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v2;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v2;   //!
   TBranch        *b_HLT_HIL2Mu15_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v2;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_v2;   //!
   TBranch        *b_HLT_HIL2Mu3_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_v2;   //!
   TBranch        *b_HLT_HIL2Mu7_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_v2;   //!
   TBranch        *b_HLT_HIL3Mu3_v2_Prescl;   //!
   TBranch        *b_HLT_HIDiJet55_v2;   //!
   TBranch        *b_HLT_HIDiJet55_v2_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v2;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v2_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v2;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v2_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v2;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v2;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v2;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v2;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v2;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v2;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v2;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v2;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v2;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet55_v2;   //!
   TBranch        *b_HLT_HIJet55_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v2;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet65_v2;   //!
   TBranch        *b_HLT_HIJet65_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet80_v2;   //!
   TBranch        *b_HLT_HIJet80_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet95_v2;   //!
   TBranch        *b_HLT_HIJet95_v2_Prescl;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v2;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v2_Prescl;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v2;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl;   //!
   TBranch        *b_HLT_HIMET120_v2;   //!
   TBranch        *b_HLT_HIMET120_v2_Prescl;   //!
   TBranch        *b_HLT_HIMET200_v2;   //!
   TBranch        *b_HLT_HIMET200_v2_Prescl;   //!
   TBranch        *b_HLT_HIMET220_v2;   //!
   TBranch        *b_HLT_HIMET220_v2_Prescl;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v2;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v2_Prescl;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v2;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v2_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_v2;   //!
   TBranch        *b_HLT_HISinglePhoton15_v2_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_v2;   //!
   TBranch        *b_HLT_HISinglePhoton20_v2_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_v2;   //!
   TBranch        *b_HLT_HISinglePhoton30_v2_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_v2;   //!
   TBranch        *b_HLT_HISinglePhoton40_v2_Prescl;   //!
   TBranch        *b_HLT_HIBptxXOR_v2;   //!
   TBranch        *b_HLT_HIBptxXOR_v2_Prescl;   //!
   TBranch        *b_HLT_HICentral10_v2;   //!
   TBranch        *b_HLT_HICentral10_v2_Prescl;   //!
   TBranch        *b_HLT_HICentralityVeto_v2;   //!
   TBranch        *b_HLT_HICentralityVeto_v2_Prescl;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v2;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v2;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v2;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHF_v2;   //!
   TBranch        *b_HLT_HIMinBiasHF_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v2;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v2;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v2;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIPhysics_v2;   //!
   TBranch        *b_HLT_HIPhysics_v2_Prescl;   //!
   TBranch        *b_HLT_HIRandom_v2;   //!
   TBranch        *b_HLT_HIRandom_v2_Prescl;   //!
   TBranch        *b_HLT_HIUCC010_v2;   //!
   TBranch        *b_HLT_HIUCC010_v2_Prescl;   //!
   TBranch        *b_HLT_HIUCC015_v2;   //!
   TBranch        *b_HLT_HIUCC015_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v2;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v2_Prescl;   //!
   TBranch        *b_HLT_HIZeroBias_v2;   //!
   TBranch        *b_HLT_HIZeroBias_v2_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v3;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v3_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v3;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v3;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v3;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v3;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v3;   //!
   TBranch        *b_HLT_HIL2Mu15_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v3;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_v3;   //!
   TBranch        *b_HLT_HIL2Mu3_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_v3;   //!
   TBranch        *b_HLT_HIL2Mu7_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_v3;   //!
   TBranch        *b_HLT_HIL3Mu3_v3_Prescl;   //!
   TBranch        *b_HLT_HIDiJet55_v3;   //!
   TBranch        *b_HLT_HIDiJet55_v3_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v3;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v3_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v3;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v3_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v3;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v3;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v3;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v3;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v3;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v3;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v3;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v3;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v3;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet55_v3;   //!
   TBranch        *b_HLT_HIJet55_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v3;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet65_v3;   //!
   TBranch        *b_HLT_HIJet65_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet80_v3;   //!
   TBranch        *b_HLT_HIJet80_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet95_v3;   //!
   TBranch        *b_HLT_HIJet95_v3_Prescl;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v3;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v3_Prescl;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v3;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl;   //!
   TBranch        *b_HLT_HIMET120_v3;   //!
   TBranch        *b_HLT_HIMET120_v3_Prescl;   //!
   TBranch        *b_HLT_HIMET200_v3;   //!
   TBranch        *b_HLT_HIMET200_v3_Prescl;   //!
   TBranch        *b_HLT_HIMET220_v3;   //!
   TBranch        *b_HLT_HIMET220_v3_Prescl;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v3;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v3_Prescl;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v3;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v3_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_v3;   //!
   TBranch        *b_HLT_HISinglePhoton15_v3_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_v3;   //!
   TBranch        *b_HLT_HISinglePhoton20_v3_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_v3;   //!
   TBranch        *b_HLT_HISinglePhoton30_v3_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_v3;   //!
   TBranch        *b_HLT_HISinglePhoton40_v3_Prescl;   //!
   TBranch        *b_HLT_HIBptxXOR_v3;   //!
   TBranch        *b_HLT_HIBptxXOR_v3_Prescl;   //!
   TBranch        *b_HLT_HICentral10_v3;   //!
   TBranch        *b_HLT_HICentral10_v3_Prescl;   //!
   TBranch        *b_HLT_HICentralityVeto_v3;   //!
   TBranch        *b_HLT_HICentralityVeto_v3_Prescl;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v3;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v3;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v3;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHF_v3;   //!
   TBranch        *b_HLT_HIMinBiasHF_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v3;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v3;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v3;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIPhysics_v3;   //!
   TBranch        *b_HLT_HIPhysics_v3_Prescl;   //!
   TBranch        *b_HLT_HIRandom_v3;   //!
   TBranch        *b_HLT_HIRandom_v3_Prescl;   //!
   TBranch        *b_HLT_HIUCC010_v3;   //!
   TBranch        *b_HLT_HIUCC010_v3_Prescl;   //!
   TBranch        *b_HLT_HIUCC015_v3;   //!
   TBranch        *b_HLT_HIUCC015_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v3;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v3_Prescl;   //!
   TBranch        *b_HLT_HIZeroBias_v3;   //!
   TBranch        *b_HLT_HIZeroBias_v3_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v4;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v4_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v4;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v4;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v4;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v4;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v4;   //!
   TBranch        *b_HLT_HIL2Mu15_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v4;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_v4;   //!
   TBranch        *b_HLT_HIL2Mu3_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_v4;   //!
   TBranch        *b_HLT_HIL2Mu7_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_v4;   //!
   TBranch        *b_HLT_HIL3Mu3_v4_Prescl;   //!
   TBranch        *b_HLT_HIDiJet55_v4;   //!
   TBranch        *b_HLT_HIDiJet55_v4_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v4;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v4_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v4;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v4_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v4;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v4;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v4;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v4;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v4;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v4;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v4;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v4;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v4;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet55_v4;   //!
   TBranch        *b_HLT_HIJet55_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v4;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet65_v4;   //!
   TBranch        *b_HLT_HIJet65_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet80_v4;   //!
   TBranch        *b_HLT_HIJet80_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet95_v4;   //!
   TBranch        *b_HLT_HIJet95_v4_Prescl;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v4;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v4_Prescl;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v4;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl;   //!
   TBranch        *b_HLT_HIMET120_v4;   //!
   TBranch        *b_HLT_HIMET120_v4_Prescl;   //!
   TBranch        *b_HLT_HIMET200_v4;   //!
   TBranch        *b_HLT_HIMET200_v4_Prescl;   //!
   TBranch        *b_HLT_HIMET220_v4;   //!
   TBranch        *b_HLT_HIMET220_v4_Prescl;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v4;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v4_Prescl;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v4;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v4_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_v4;   //!
   TBranch        *b_HLT_HISinglePhoton15_v4_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_v4;   //!
   TBranch        *b_HLT_HISinglePhoton20_v4_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_v4;   //!
   TBranch        *b_HLT_HISinglePhoton30_v4_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_v4;   //!
   TBranch        *b_HLT_HISinglePhoton40_v4_Prescl;   //!
   TBranch        *b_HLT_HIBptxXOR_v4;   //!
   TBranch        *b_HLT_HIBptxXOR_v4_Prescl;   //!
   TBranch        *b_HLT_HICentral10_v4;   //!
   TBranch        *b_HLT_HICentral10_v4_Prescl;   //!
   TBranch        *b_HLT_HICentralityVeto_v4;   //!
   TBranch        *b_HLT_HICentralityVeto_v4_Prescl;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v4;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v4;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v4;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHF_v4;   //!
   TBranch        *b_HLT_HIMinBiasHF_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v4;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v4;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v4;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIPhysics_v4;   //!
   TBranch        *b_HLT_HIPhysics_v4_Prescl;   //!
   TBranch        *b_HLT_HIRandom_v4;   //!
   TBranch        *b_HLT_HIRandom_v4_Prescl;   //!
   TBranch        *b_HLT_HIUCC010_v4;   //!
   TBranch        *b_HLT_HIUCC010_v4_Prescl;   //!
   TBranch        *b_HLT_HIUCC015_v4;   //!
   TBranch        *b_HLT_HIUCC015_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v4;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v4_Prescl;   //!
   TBranch        *b_HLT_HIZeroBias_v4;   //!
   TBranch        *b_HLT_HIZeroBias_v4_Prescl;   //!
   TBranch        *b_L1_BptxMinus;   //!
   TBranch        *b_L1_BptxMinus_Prescl;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus_Prescl;   //!
   TBranch        *b_L1_BptxPlus;   //!
   TBranch        *b_L1_BptxPlus_Prescl;   //!
   TBranch        *b_L1_BptxPlusANDMinus;   //!
   TBranch        *b_L1_BptxPlusANDMinus_Prescl;   //!
   TBranch        *b_L1_BptxPlusORMinus;   //!
   TBranch        *b_L1_BptxPlusORMinus_Prescl;   //!
   TBranch        *b_L1_BptxPlusORMinus_instance1;   //!
   TBranch        *b_L1_BptxPlusORMinus_instance1_Prescl;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus_Prescl;   //!
   TBranch        *b_L1_BptxXOR;   //!
   TBranch        *b_L1_BptxXOR_Prescl;   //!
   TBranch        *b_L1_BptxXOR_BscMinBiasOR;   //!
   TBranch        *b_L1_BptxXOR_BscMinBiasOR_Prescl;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus_Prescl;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus_Prescl;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1_Prescl;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold2_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold2_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusORMinus;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusORMinus_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold1;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold2_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasThreshold2_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleEG10;   //!
   TBranch        *b_L1_DoubleEG10_Prescl;   //!
   TBranch        *b_L1_DoubleEG3;   //!
   TBranch        *b_L1_DoubleEG3_Prescl;   //!
   TBranch        *b_L1_DoubleEG5_BptxAND;   //!
   TBranch        *b_L1_DoubleEG5_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleEG8_BptxAND;   //!
   TBranch        *b_L1_DoubleEG8_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleForJet32_EtaOpp;   //!
   TBranch        *b_L1_DoubleForJet32_EtaOpp_Prescl;   //!
   TBranch        *b_L1_DoubleForJet44_EtaOpp;   //!
   TBranch        *b_L1_DoubleForJet44_EtaOpp_Prescl;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing1_P1N1;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing1_P1N1_Prescl;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing2_P1N1;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing2_P1N1_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P200N200;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P200N200_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P4N4;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P4N4_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P200N200;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P200N200_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P4N4;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P4N4_Prescl;   //!
   TBranch        *b_L1_DoubleMu3_BptxAND;   //!
   TBranch        *b_L1_DoubleMu3_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND_Prescl;   //!
   TBranch        *b_L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND;   //!
   TBranch        *b_L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;   //!
   TBranch        *b_L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;   //!
   TBranch        *b_L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;   //!
   TBranch        *b_L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND;   //!
   TBranch        *b_L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;   //!
   TBranch        *b_L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;   //!
   TBranch        *b_L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETM30_BptxAND;   //!
   TBranch        *b_L1_ETM30_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETM50_BptxAND;   //!
   TBranch        *b_L1_ETM50_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT100_BptxAND;   //!
   TBranch        *b_L1_ETT100_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT140_BptxAND;   //!
   TBranch        *b_L1_ETT140_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT2000;   //!
   TBranch        *b_L1_ETT2000_Prescl;   //!
   TBranch        *b_L1_ETT220_BptxAND;   //!
   TBranch        *b_L1_ETT220_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT70_BptxAND;   //!
   TBranch        *b_L1_ETT70_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT800;   //!
   TBranch        *b_L1_ETT800_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincidencePm;   //!
   TBranch        *b_L1_HcalHfCoincidencePm_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincidencePm_BptxAND;   //!
   TBranch        *b_L1_HcalHfCoincidencePm_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHfMmOrPpOrPm;   //!
   TBranch        *b_L1_HcalHfMmOrPpOrPm_Prescl;   //!
   TBranch        *b_L1_HcalHfMmOrPpOrPm_BptxAND;   //!
   TBranch        *b_L1_HcalHfMmOrPpOrPm_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHfMmpOrMpp;   //!
   TBranch        *b_L1_HcalHfMmpOrMpp_Prescl;   //!
   TBranch        *b_L1_HcalHfMmpOrMpp_BptxAND;   //!
   TBranch        *b_L1_HcalHfMmpOrMpp_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHoTotalOR;   //!
   TBranch        *b_L1_HcalHoTotalOR_Prescl;   //!
   TBranch        *b_L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND;   //!
   TBranch        *b_L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;   //!
   TBranch        *b_L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;   //!
   TBranch        *b_L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasOR;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasThreshold1;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasThreshold2;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfCoincidencePm;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfMmpOrMpp;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasInnerThreshold1;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasInnerThreshold2;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasOR;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasOR_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasThreshold1;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasThreshold1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasThreshold2;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasThreshold2_Prescl;   //!
   TBranch        *b_L1_NotBsc2_HcalHfCoincidencePm;   //!
   TBranch        *b_L1_NotBsc2_HcalHfCoincidencePm_Prescl;   //!
   TBranch        *b_L1_NotBsc2_HcalHfMmOrPpOrPm;   //!
   TBranch        *b_L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl;   //!
   TBranch        *b_L1_NotBsc2_HcalHfMmpOrMpp;   //!
   TBranch        *b_L1_NotBsc2_HcalHfMmpOrMpp_Prescl;   //!
   TBranch        *b_L1_SingleEG12;   //!
   TBranch        *b_L1_SingleEG12_Prescl;   //!
   TBranch        *b_L1_SingleEG15;   //!
   TBranch        *b_L1_SingleEG15_Prescl;   //!
   TBranch        *b_L1_SingleEG18;   //!
   TBranch        *b_L1_SingleEG18_Prescl;   //!
   TBranch        *b_L1_SingleEG20;   //!
   TBranch        *b_L1_SingleEG20_Prescl;   //!
   TBranch        *b_L1_SingleEG22;   //!
   TBranch        *b_L1_SingleEG22_Prescl;   //!
   TBranch        *b_L1_SingleEG2_BptxAND;   //!
   TBranch        *b_L1_SingleEG2_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG30;   //!
   TBranch        *b_L1_SingleEG30_Prescl;   //!
   TBranch        *b_L1_SingleEG5_BptxAND;   //!
   TBranch        *b_L1_SingleEG5_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG8_BptxAND;   //!
   TBranch        *b_L1_SingleEG8_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet128;   //!
   TBranch        *b_L1_SingleJet128_Prescl;   //!
   TBranch        *b_L1_SingleJet128_BptxAND;   //!
   TBranch        *b_L1_SingleJet128_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet128_Central_BptxAND;   //!
   TBranch        *b_L1_SingleJet128_Central_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet16;   //!
   TBranch        *b_L1_SingleJet16_Prescl;   //!
   TBranch        *b_L1_SingleJet16_Central_NotETT140_BptxAND;   //!
   TBranch        *b_L1_SingleJet16_Central_NotETT140_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo;   //!
   TBranch        *b_L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo;   //!
   TBranch        *b_L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleJet36_BptxAND;   //!
   TBranch        *b_L1_SingleJet36_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet52_BptxAND;   //!
   TBranch        *b_L1_SingleJet52_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet52_Central_BptxAND;   //!
   TBranch        *b_L1_SingleJet52_Central_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet68_BptxAND;   //!
   TBranch        *b_L1_SingleJet68_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet80_BptxAND;   //!
   TBranch        *b_L1_SingleJet80_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet80_Central;   //!
   TBranch        *b_L1_SingleJet80_Central_Prescl;   //!
   TBranch        *b_L1_SingleJet80_Central_BptxAND;   //!
   TBranch        *b_L1_SingleJet80_Central_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet92;   //!
   TBranch        *b_L1_SingleJet92_Prescl;   //!
   TBranch        *b_L1_SingleJet92_BptxAND;   //!
   TBranch        *b_L1_SingleJet92_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet92_Central;   //!
   TBranch        *b_L1_SingleJet92_Central_Prescl;   //!
   TBranch        *b_L1_SingleJet92_Central_BptxAND;   //!
   TBranch        *b_L1_SingleJet92_Central_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu0_BptxAND;   //!
   TBranch        *b_L1_SingleMu0_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu10;   //!
   TBranch        *b_L1_SingleMu10_Prescl;   //!
   TBranch        *b_L1_SingleMu14;   //!
   TBranch        *b_L1_SingleMu14_Prescl;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu20_Prescl;   //!
   TBranch        *b_L1_SingleMu3_BptxAND;   //!
   TBranch        *b_L1_SingleMu3_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu5_BptxAND;   //!
   TBranch        *b_L1_SingleMu5_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMu7_Prescl;   //!
   TBranch        *b_L1_SingleMuBeamHalo;   //!
   TBranch        *b_L1_SingleMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen_BptxAND;   //!
   TBranch        *b_L1_SingleMuOpen_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZdcCaloMinus;   //!
   TBranch        *b_L1_ZdcCaloMinus_Prescl;   //!
   TBranch        *b_L1_ZdcCaloMinus_BptxAND;   //!
   TBranch        *b_L1_ZdcCaloMinus_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZdcCaloPlus;   //!
   TBranch        *b_L1_ZdcCaloPlus_Prescl;   //!
   TBranch        *b_L1_ZdcCaloPlus_BptxAND;   //!
   TBranch        *b_L1_ZdcCaloPlus_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZdcCaloPlus_ZdcCaloMinus;   //!
   TBranch        *b_L1_ZdcCaloPlus_ZdcCaloMinus_Prescl;   //!
   TBranch        *b_L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND;   //!
   TBranch        *b_L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_HighMultiplicity_v0;   //!
   TBranch        *b_L1Tech_BSC_HighMultiplicity_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_inner_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_inner_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_outer_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_outer_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_inner_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_inner_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_outer_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_outer_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_OR_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_OR_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold1_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold2_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold2_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold1_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold2_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold2_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_splash_beam1_v0;   //!
   TBranch        *b_L1Tech_BSC_splash_beam1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_splash_beam2_v0;   //!
   TBranch        *b_L1Tech_BSC_splash_beam2_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0_Prescl;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St1Sect45_down_v0;   //!
   TBranch        *b_L1Tech_FSC_St1Sect45_down_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St1Sect45_upp_v0;   //!
   TBranch        *b_L1Tech_FSC_St1Sect45_upp_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St1Sect56_down_v0;   //!
   TBranch        *b_L1Tech_FSC_St1Sect56_down_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St1Sect56_upp_v0;   //!
   TBranch        *b_L1Tech_FSC_St1Sect56_upp_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St2Sect45_down_v0;   //!
   TBranch        *b_L1Tech_FSC_St2Sect45_down_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St2Sect45_upp_v0;   //!
   TBranch        *b_L1Tech_FSC_St2Sect45_upp_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St2Sect56_down_v0;   //!
   TBranch        *b_L1Tech_FSC_St2Sect56_down_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St2Sect56_upp_v0;   //!
   TBranch        *b_L1Tech_FSC_St2Sect56_upp_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_downLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_downLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_downRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_downRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_uppLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_uppLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_uppRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_uppRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RB0_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBst1_collisions_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;   //!

};


void setupHltTree(TTree *t,Hlts &tHlts,bool doCheck = 1)
{
   // Set branch addresses and branch pointers
   if (t->GetBranch("NL1IsolEm")) t->SetBranchAddress("NL1IsolEm", &tHlts.NL1IsolEm, &tHlts.b_NL1IsolEm);
   if (t->GetBranch("L1IsolEmEt")) t->SetBranchAddress("L1IsolEmEt", tHlts.L1IsolEmEt, &tHlts.b_L1IsolEmEt);
   if (t->GetBranch("L1IsolEmE")) t->SetBranchAddress("L1IsolEmE", tHlts.L1IsolEmE, &tHlts.b_L1IsolEmE);
   if (t->GetBranch("L1IsolEmEta")) t->SetBranchAddress("L1IsolEmEta", tHlts.L1IsolEmEta, &tHlts.b_L1IsolEmEta);
   if (t->GetBranch("L1IsolEmPhi")) t->SetBranchAddress("L1IsolEmPhi", tHlts.L1IsolEmPhi, &tHlts.b_L1IsolEmPhi);
   if (t->GetBranch("NL1NIsolEm")) t->SetBranchAddress("NL1NIsolEm", &tHlts.NL1NIsolEm, &tHlts.b_NL1NIsolEm);
   if (t->GetBranch("L1NIsolEmEt")) t->SetBranchAddress("L1NIsolEmEt", tHlts.L1NIsolEmEt, &tHlts.b_L1NIsolEmEt);
   if (t->GetBranch("L1NIsolEmE")) t->SetBranchAddress("L1NIsolEmE", tHlts.L1NIsolEmE, &tHlts.b_L1NIsolEmE);
   if (t->GetBranch("L1NIsolEmEta")) t->SetBranchAddress("L1NIsolEmEta", tHlts.L1NIsolEmEta, &tHlts.b_L1NIsolEmEta);
   if (t->GetBranch("L1NIsolEmPhi")) t->SetBranchAddress("L1NIsolEmPhi", tHlts.L1NIsolEmPhi, &tHlts.b_L1NIsolEmPhi);
   if (t->GetBranch("NL1Mu")) t->SetBranchAddress("NL1Mu", &tHlts.NL1Mu, &tHlts.b_NL1Mu);
   if (t->GetBranch("L1MuPt")) t->SetBranchAddress("L1MuPt", tHlts.L1MuPt, &tHlts.b_L1MuPt);
   if (t->GetBranch("L1MuE")) t->SetBranchAddress("L1MuE", tHlts.L1MuE, &tHlts.b_L1MuE);
   if (t->GetBranch("L1MuEta")) t->SetBranchAddress("L1MuEta", tHlts.L1MuEta, &tHlts.b_L1MuEta);
   if (t->GetBranch("L1MuPhi")) t->SetBranchAddress("L1MuPhi", tHlts.L1MuPhi, &tHlts.b_L1MuPhi);
   if (t->GetBranch("L1MuIsol")) t->SetBranchAddress("L1MuIsol", tHlts.L1MuIsol, &tHlts.b_L1MuIsol);
   if (t->GetBranch("L1MuMip")) t->SetBranchAddress("L1MuMip", tHlts.L1MuMip, &tHlts.b_L1MuMip);
   if (t->GetBranch("L1MuFor")) t->SetBranchAddress("L1MuFor", tHlts.L1MuFor, &tHlts.b_L1MuFor);
   if (t->GetBranch("L1MuRPC")) t->SetBranchAddress("L1MuRPC", tHlts.L1MuRPC, &tHlts.b_L1MuRPC);
   if (t->GetBranch("L1MuQal")) t->SetBranchAddress("L1MuQal", tHlts.L1MuQal, &tHlts.b_L1MuQal);
   if (t->GetBranch("L1MuChg")) t->SetBranchAddress("L1MuChg", tHlts.L1MuChg, &tHlts.b_L1MuChg);
   if (t->GetBranch("NL1CenJet")) t->SetBranchAddress("NL1CenJet", &tHlts.NL1CenJet, &tHlts.b_NL1CenJet);
   if (t->GetBranch("L1CenJetEt")) t->SetBranchAddress("L1CenJetEt", tHlts.L1CenJetEt, &tHlts.b_L1CenJetEt);
   if (t->GetBranch("L1CenJetE")) t->SetBranchAddress("L1CenJetE", tHlts.L1CenJetE, &tHlts.b_L1CenJetE);
   if (t->GetBranch("L1CenJetEta")) t->SetBranchAddress("L1CenJetEta", tHlts.L1CenJetEta, &tHlts.b_L1CenJetEta);
   if (t->GetBranch("L1CenJetPhi")) t->SetBranchAddress("L1CenJetPhi", tHlts.L1CenJetPhi, &tHlts.b_L1CenJetPhi);
   if (t->GetBranch("NL1ForJet")) t->SetBranchAddress("NL1ForJet", &tHlts.NL1ForJet, &tHlts.b_NL1ForJet);
   if (t->GetBranch("L1ForJetEt")) t->SetBranchAddress("L1ForJetEt", tHlts.L1ForJetEt, &tHlts.b_L1ForJetEt);
   if (t->GetBranch("L1ForJetE")) t->SetBranchAddress("L1ForJetE", tHlts.L1ForJetE, &tHlts.b_L1ForJetE);
   if (t->GetBranch("L1ForJetEta")) t->SetBranchAddress("L1ForJetEta", tHlts.L1ForJetEta, &tHlts.b_L1ForJetEta);
   if (t->GetBranch("L1ForJetPhi")) t->SetBranchAddress("L1ForJetPhi", tHlts.L1ForJetPhi, &tHlts.b_L1ForJetPhi);
   if (t->GetBranch("NL1Tau")) t->SetBranchAddress("NL1Tau", &tHlts.NL1Tau, &tHlts.b_NL1Tau);
   if (t->GetBranch("L1TauEt")) t->SetBranchAddress("L1TauEt", tHlts.L1TauEt, &tHlts.b_L1TauEt);
   if (t->GetBranch("L1TauE")) t->SetBranchAddress("L1TauE", tHlts.L1TauE, &tHlts.b_L1TauE);
   if (t->GetBranch("L1TauEta")) t->SetBranchAddress("L1TauEta", tHlts.L1TauEta, &tHlts.b_L1TauEta);
   if (t->GetBranch("L1TauPhi")) t->SetBranchAddress("L1TauPhi", tHlts.L1TauPhi, &tHlts.b_L1TauPhi);
   if (t->GetBranch("L1Met")) t->SetBranchAddress("L1Met", &tHlts.L1Met, &tHlts.b_L1Met);
   if (t->GetBranch("L1MetPhi")) t->SetBranchAddress("L1MetPhi", &tHlts.L1MetPhi, &tHlts.b_L1MetPhi);
   if (t->GetBranch("L1EtTot")) t->SetBranchAddress("L1EtTot", &tHlts.L1EtTot, &tHlts.b_L1EtTot);
   if (t->GetBranch("L1Mht")) t->SetBranchAddress("L1Mht", &tHlts.L1Mht, &tHlts.b_L1Mht);
   if (t->GetBranch("L1MhtPhi")) t->SetBranchAddress("L1MhtPhi", &tHlts.L1MhtPhi, &tHlts.b_L1MhtPhi);
   if (t->GetBranch("L1EtHad")) t->SetBranchAddress("L1EtHad", &tHlts.L1EtHad, &tHlts.b_L1EtHad);
   if (t->GetBranch("L1HfRing1EtSumPositiveEta")) t->SetBranchAddress("L1HfRing1EtSumPositiveEta", &tHlts.L1HfRing1EtSumPositiveEta, &tHlts.b_L1HfRing1EtSumPositiveEta);
   if (t->GetBranch("L1HfRing2EtSumPositiveEta")) t->SetBranchAddress("L1HfRing2EtSumPositiveEta", &tHlts.L1HfRing2EtSumPositiveEta, &tHlts.b_L1HfRing2EtSumPositiveEta);
   if (t->GetBranch("L1HfRing1EtSumNegativeEta")) t->SetBranchAddress("L1HfRing1EtSumNegativeEta", &tHlts.L1HfRing1EtSumNegativeEta, &tHlts.b_L1HfRing1EtSumNegativeEta);
   if (t->GetBranch("L1HfRing2EtSumNegativeEta")) t->SetBranchAddress("L1HfRing2EtSumNegativeEta", &tHlts.L1HfRing2EtSumNegativeEta, &tHlts.b_L1HfRing2EtSumNegativeEta);
   if (t->GetBranch("L1HfTowerCountPositiveEtaRing1")) t->SetBranchAddress("L1HfTowerCountPositiveEtaRing1", &tHlts.L1HfTowerCountPositiveEtaRing1, &tHlts.b_L1HfTowerCountPositiveEtaRing1);
   if (t->GetBranch("L1HfTowerCountNegativeEtaRing1")) t->SetBranchAddress("L1HfTowerCountNegativeEtaRing1", &tHlts.L1HfTowerCountNegativeEtaRing1, &tHlts.b_L1HfTowerCountNegativeEtaRing1);
   if (t->GetBranch("L1HfTowerCountPositiveEtaRing2")) t->SetBranchAddress("L1HfTowerCountPositiveEtaRing2", &tHlts.L1HfTowerCountPositiveEtaRing2, &tHlts.b_L1HfTowerCountPositiveEtaRing2);
   if (t->GetBranch("L1HfTowerCountNegativeEtaRing2")) t->SetBranchAddress("L1HfTowerCountNegativeEtaRing2", &tHlts.L1HfTowerCountNegativeEtaRing2, &tHlts.b_L1HfTowerCountNegativeEtaRing2);
   if (t->GetBranch("Run")) t->SetBranchAddress("Run", &tHlts.Run, &tHlts.b_Run);
   if (t->GetBranch("Event")) t->SetBranchAddress("Event", &tHlts.Event, &tHlts.b_Event);
   if (t->GetBranch("LumiBlock")) t->SetBranchAddress("LumiBlock", &tHlts.LumiBlock, &tHlts.b_LumiBlock);
   if (t->GetBranch("Bx")) t->SetBranchAddress("Bx", &tHlts.Bx, &tHlts.b_Bx);
   if (t->GetBranch("Orbit")) t->SetBranchAddress("Orbit", &tHlts.Orbit, &tHlts.b_Orbit);
   if (t->GetBranch("AvgInstDelLumi")) t->SetBranchAddress("AvgInstDelLumi", &tHlts.AvgInstDelLumi, &tHlts.b_AvgInstDelLumi);
   if (t->GetBranch("raw2digi_step")) t->SetBranchAddress("raw2digi_step", &tHlts.raw2digi_step, &tHlts.b_raw2digi_step);
   if (t->GetBranch("raw2digi_step_Prescl")) t->SetBranchAddress("raw2digi_step_Prescl", &tHlts.raw2digi_step_Prescl, &tHlts.b_raw2digi_step_Prescl);
   if (t->GetBranch("reconstruction_step")) t->SetBranchAddress("reconstruction_step", &tHlts.reconstruction_step, &tHlts.b_reconstruction_step);
   if (t->GetBranch("reconstruction_step_Prescl")) t->SetBranchAddress("reconstruction_step_Prescl", &tHlts.reconstruction_step_Prescl, &tHlts.b_reconstruction_step_Prescl);
   if (t->GetBranch("HLTriggerFirstPath")) t->SetBranchAddress("HLTriggerFirstPath", &tHlts.HLTriggerFirstPath, &tHlts.b_HLTriggerFirstPath);
   if (t->GetBranch("HLTriggerFirstPath_Prescl")) t->SetBranchAddress("HLTriggerFirstPath_Prescl", &tHlts.HLTriggerFirstPath_Prescl, &tHlts.b_HLTriggerFirstPath_Prescl);
   if (t->GetBranch("HLT_HIDTCalibration_v1")) t->SetBranchAddress("HLT_HIDTCalibration_v1", &tHlts.HLT_HIDTCalibration_v1, &tHlts.b_HLT_HIDTCalibration_v1);
   if (t->GetBranch("HLT_HIDTCalibration_v1_Prescl")) t->SetBranchAddress("HLT_HIDTCalibration_v1_Prescl", &tHlts.HLT_HIDTCalibration_v1_Prescl, &tHlts.b_HLT_HIDTCalibration_v1_Prescl);
   if (t->GetBranch("HLT_HIEcalCalibration_v1")) t->SetBranchAddress("HLT_HIEcalCalibration_v1", &tHlts.HLT_HIEcalCalibration_v1, &tHlts.b_HLT_HIEcalCalibration_v1);
   if (t->GetBranch("HLT_HIEcalCalibration_v1_Prescl")) t->SetBranchAddress("HLT_HIEcalCalibration_v1_Prescl", &tHlts.HLT_HIEcalCalibration_v1_Prescl, &tHlts.b_HLT_HIEcalCalibration_v1_Prescl);
   if (t->GetBranch("HLT_HIHcalCalibration_v1")) t->SetBranchAddress("HLT_HIHcalCalibration_v1", &tHlts.HLT_HIHcalCalibration_v1, &tHlts.b_HLT_HIHcalCalibration_v1);
   if (t->GetBranch("HLT_HIHcalCalibration_v1_Prescl")) t->SetBranchAddress("HLT_HIHcalCalibration_v1_Prescl", &tHlts.HLT_HIHcalCalibration_v1_Prescl, &tHlts.b_HLT_HIHcalCalibration_v1_Prescl);
   if (t->GetBranch("HLT_HIZeroBias_v1")) t->SetBranchAddress("HLT_HIZeroBias_v1", &tHlts.HLT_HIZeroBias_v1, &tHlts.b_HLT_HIZeroBias_v1);
   if (t->GetBranch("HLT_HIZeroBias_v1_Prescl")) t->SetBranchAddress("HLT_HIZeroBias_v1_Prescl", &tHlts.HLT_HIZeroBias_v1_Prescl, &tHlts.b_HLT_HIZeroBias_v1_Prescl);
   if (t->GetBranch("HLT_HIZeroBiasXOR_v1")) t->SetBranchAddress("HLT_HIZeroBiasXOR_v1", &tHlts.HLT_HIZeroBiasXOR_v1, &tHlts.b_HLT_HIZeroBiasXOR_v1);
   if (t->GetBranch("HLT_HIZeroBiasXOR_v1_Prescl")) t->SetBranchAddress("HLT_HIZeroBiasXOR_v1_Prescl", &tHlts.HLT_HIZeroBiasXOR_v1_Prescl, &tHlts.b_HLT_HIZeroBiasXOR_v1_Prescl);
   if (t->GetBranch("HLT_HIZeroBiasPixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v1", &tHlts.HLT_HIZeroBiasPixel_SingleTrack_v1, &tHlts.b_HLT_HIZeroBiasPixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasBSC_v1")) t->SetBranchAddress("HLT_HIMinBiasBSC_v1", &tHlts.HLT_HIMinBiasBSC_v1, &tHlts.b_HLT_HIMinBiasBSC_v1);
   if (t->GetBranch("HLT_HIMinBiasBSC_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasBSC_v1_Prescl", &tHlts.HLT_HIMinBiasBSC_v1_Prescl, &tHlts.b_HLT_HIMinBiasBSC_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasBSC_OR_v1")) t->SetBranchAddress("HLT_HIMinBiasBSC_OR_v1", &tHlts.HLT_HIMinBiasBSC_OR_v1, &tHlts.b_HLT_HIMinBiasBSC_OR_v1);
   if (t->GetBranch("HLT_HIMinBiasBSC_OR_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasBSC_OR_v1_Prescl", &tHlts.HLT_HIMinBiasBSC_OR_v1_Prescl, &tHlts.b_HLT_HIMinBiasBSC_OR_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHF_v1")) t->SetBranchAddress("HLT_HIMinBiasHF_v1", &tHlts.HLT_HIMinBiasHF_v1, &tHlts.b_HLT_HIMinBiasHF_v1);
   if (t->GetBranch("HLT_HIMinBiasHF_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHF_v1_Prescl", &tHlts.HLT_HIMinBiasHF_v1_Prescl, &tHlts.b_HLT_HIMinBiasHF_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHf_OR_v1")) t->SetBranchAddress("HLT_HIMinBiasHf_OR_v1", &tHlts.HLT_HIMinBiasHf_OR_v1, &tHlts.b_HLT_HIMinBiasHf_OR_v1);
   if (t->GetBranch("HLT_HIMinBiasHf_OR_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHf_OR_v1_Prescl", &tHlts.HLT_HIMinBiasHf_OR_v1_Prescl, &tHlts.b_HLT_HIMinBiasHf_OR_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHfOrBSC_v1")) t->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1", &tHlts.HLT_HIMinBiasHfOrBSC_v1, &tHlts.b_HLT_HIMinBiasHfOrBSC_v1);
   if (t->GetBranch("HLT_HIMinBiasHfOrBSC_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1_Prescl", &tHlts.HLT_HIMinBiasHfOrBSC_v1_Prescl, &tHlts.b_HLT_HIMinBiasHfOrBSC_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasPixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v1", &tHlts.HLT_HIMinBiasPixel_SingleTrack_v1, &tHlts.b_HLT_HIMinBiasPixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIMinBiasPixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIMinBiasPixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIMinBiasPixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_v1")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v1", &tHlts.HLT_HIMinBiasZDC_Calo_v1, &tHlts.b_HLT_HIMinBiasZDC_Calo_v1);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v1_Prescl", &tHlts.HLT_HIMinBiasZDC_Calo_v1_Prescl, &tHlts.b_HLT_HIMinBiasZDC_Calo_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1", &tHlts.HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1, &tHlts.b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl", &tHlts.HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl, &tHlts.b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDCPixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v1", &tHlts.HLT_HIMinBiasZDCPixel_SingleTrack_v1, &tHlts.b_HLT_HIMinBiasZDCPixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1", &tHlts.HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1, &tHlts.b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLT_HIBptxXOR_v1")) t->SetBranchAddress("HLT_HIBptxXOR_v1", &tHlts.HLT_HIBptxXOR_v1, &tHlts.b_HLT_HIBptxXOR_v1);
   if (t->GetBranch("HLT_HIBptxXOR_v1_Prescl")) t->SetBranchAddress("HLT_HIBptxXOR_v1_Prescl", &tHlts.HLT_HIBptxXOR_v1_Prescl, &tHlts.b_HLT_HIBptxXOR_v1_Prescl);
   if (t->GetBranch("HLT_HIL1Algo_BptxXOR_BSC_OR_v1")) t->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v1", &tHlts.HLT_HIL1Algo_BptxXOR_BSC_OR_v1, &tHlts.b_HLT_HIL1Algo_BptxXOR_BSC_OR_v1);
   if (t->GetBranch("HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl")) t->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl", &tHlts.HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl, &tHlts.b_HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl);
   if (t->GetBranch("HLT_HIL1DoubleMuOpen_v1")) t->SetBranchAddress("HLT_HIL1DoubleMuOpen_v1", &tHlts.HLT_HIL1DoubleMuOpen_v1, &tHlts.b_HLT_HIL1DoubleMuOpen_v1);
   if (t->GetBranch("HLT_HIL1DoubleMuOpen_v1_Prescl")) t->SetBranchAddress("HLT_HIL1DoubleMuOpen_v1_Prescl", &tHlts.HLT_HIL1DoubleMuOpen_v1_Prescl, &tHlts.b_HLT_HIL1DoubleMuOpen_v1_Prescl);
   if (t->GetBranch("HLT_HIL1DoubleMu0_HighQ_v1")) t->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v1", &tHlts.HLT_HIL1DoubleMu0_HighQ_v1, &tHlts.b_HLT_HIL1DoubleMu0_HighQ_v1);
   if (t->GetBranch("HLT_HIL1DoubleMu0_HighQ_v1_Prescl")) t->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v1_Prescl", &tHlts.HLT_HIL1DoubleMu0_HighQ_v1_Prescl, &tHlts.b_HLT_HIL1DoubleMu0_HighQ_v1_Prescl);
   if (t->GetBranch("HLT_HIL2Mu3_v1")) t->SetBranchAddress("HLT_HIL2Mu3_v1", &tHlts.HLT_HIL2Mu3_v1, &tHlts.b_HLT_HIL2Mu3_v1);
   if (t->GetBranch("HLT_HIL2Mu3_v1_Prescl")) t->SetBranchAddress("HLT_HIL2Mu3_v1_Prescl", &tHlts.HLT_HIL2Mu3_v1_Prescl, &tHlts.b_HLT_HIL2Mu3_v1_Prescl);
   if (t->GetBranch("HLT_HIL2Mu3_NHitQ_v1")) t->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v1", &tHlts.HLT_HIL2Mu3_NHitQ_v1, &tHlts.b_HLT_HIL2Mu3_NHitQ_v1);
   if (t->GetBranch("HLT_HIL2Mu3_NHitQ_v1_Prescl")) t->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v1_Prescl", &tHlts.HLT_HIL2Mu3_NHitQ_v1_Prescl, &tHlts.b_HLT_HIL2Mu3_NHitQ_v1_Prescl);
   if (t->GetBranch("HLT_HIL2Mu7_v1")) t->SetBranchAddress("HLT_HIL2Mu7_v1", &tHlts.HLT_HIL2Mu7_v1, &tHlts.b_HLT_HIL2Mu7_v1);
   if (t->GetBranch("HLT_HIL2Mu7_v1_Prescl")) t->SetBranchAddress("HLT_HIL2Mu7_v1_Prescl", &tHlts.HLT_HIL2Mu7_v1_Prescl, &tHlts.b_HLT_HIL2Mu7_v1_Prescl);
   if (t->GetBranch("HLT_HIL2Mu15_v1")) t->SetBranchAddress("HLT_HIL2Mu15_v1", &tHlts.HLT_HIL2Mu15_v1, &tHlts.b_HLT_HIL2Mu15_v1);
   if (t->GetBranch("HLT_HIL2Mu15_v1_Prescl")) t->SetBranchAddress("HLT_HIL2Mu15_v1_Prescl", &tHlts.HLT_HIL2Mu15_v1_Prescl, &tHlts.b_HLT_HIL2Mu15_v1_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_v1")) t->SetBranchAddress("HLT_HIL2DoubleMu0_v1", &tHlts.HLT_HIL2DoubleMu0_v1, &tHlts.b_HLT_HIL2DoubleMu0_v1);
   if (t->GetBranch("HLT_HIL2DoubleMu0_v1_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_v1_Prescl", &tHlts.HLT_HIL2DoubleMu0_v1_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_v1_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_NHitQ_v1")) t->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v1", &tHlts.HLT_HIL2DoubleMu0_NHitQ_v1, &tHlts.b_HLT_HIL2DoubleMu0_NHitQ_v1);
   if (t->GetBranch("HLT_HIL2DoubleMu0_NHitQ_v1_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v1_Prescl", &tHlts.HLT_HIL2DoubleMu0_NHitQ_v1_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_NHitQ_v1_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1")) t->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1", &tHlts.HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1, &tHlts.b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1);
   if (t->GetBranch("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl", &tHlts.HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu3_v1")) t->SetBranchAddress("HLT_HIL2DoubleMu3_v1", &tHlts.HLT_HIL2DoubleMu3_v1, &tHlts.b_HLT_HIL2DoubleMu3_v1);
   if (t->GetBranch("HLT_HIL2DoubleMu3_v1_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu3_v1_Prescl", &tHlts.HLT_HIL2DoubleMu3_v1_Prescl, &tHlts.b_HLT_HIL2DoubleMu3_v1_Prescl);
   if (t->GetBranch("HLT_HIL3Mu3_v1")) t->SetBranchAddress("HLT_HIL3Mu3_v1", &tHlts.HLT_HIL3Mu3_v1, &tHlts.b_HLT_HIL3Mu3_v1);
   if (t->GetBranch("HLT_HIL3Mu3_v1_Prescl")) t->SetBranchAddress("HLT_HIL3Mu3_v1_Prescl", &tHlts.HLT_HIL3Mu3_v1_Prescl, &tHlts.b_HLT_HIL3Mu3_v1_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_v1")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_v1", &tHlts.HLT_HIL3DoubleMuOpen_v1, &tHlts.b_HLT_HIL3DoubleMuOpen_v1);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_v1_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_v1_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_v1_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_v1_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_v1")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v1", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_v1, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_v1);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_SS_v1")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v1", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_SS_v1, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v1);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_v1")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v1", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_v1, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v1);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton15_v1")) t->SetBranchAddress("HLT_HISinglePhoton15_v1", &tHlts.HLT_HISinglePhoton15_v1, &tHlts.b_HLT_HISinglePhoton15_v1);
   if (t->GetBranch("HLT_HISinglePhoton15_v1_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton15_v1_Prescl", &tHlts.HLT_HISinglePhoton15_v1_Prescl, &tHlts.b_HLT_HISinglePhoton15_v1_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton20_v1")) t->SetBranchAddress("HLT_HISinglePhoton20_v1", &tHlts.HLT_HISinglePhoton20_v1, &tHlts.b_HLT_HISinglePhoton20_v1);
   if (t->GetBranch("HLT_HISinglePhoton20_v1_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton20_v1_Prescl", &tHlts.HLT_HISinglePhoton20_v1_Prescl, &tHlts.b_HLT_HISinglePhoton20_v1_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton30_v1")) t->SetBranchAddress("HLT_HISinglePhoton30_v1", &tHlts.HLT_HISinglePhoton30_v1, &tHlts.b_HLT_HISinglePhoton30_v1);
   if (t->GetBranch("HLT_HISinglePhoton30_v1_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton30_v1_Prescl", &tHlts.HLT_HISinglePhoton30_v1_Prescl, &tHlts.b_HLT_HISinglePhoton30_v1_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton40_v1")) t->SetBranchAddress("HLT_HISinglePhoton40_v1", &tHlts.HLT_HISinglePhoton40_v1, &tHlts.b_HLT_HISinglePhoton40_v1);
   if (t->GetBranch("HLT_HISinglePhoton40_v1_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton40_v1_Prescl", &tHlts.HLT_HISinglePhoton40_v1_Prescl, &tHlts.b_HLT_HISinglePhoton40_v1_Prescl);
   if (t->GetBranch("HLT_HIPhoton10_Photon15_v1")) t->SetBranchAddress("HLT_HIPhoton10_Photon15_v1", &tHlts.HLT_HIPhoton10_Photon15_v1, &tHlts.b_HLT_HIPhoton10_Photon15_v1);
   if (t->GetBranch("HLT_HIPhoton10_Photon15_v1_Prescl")) t->SetBranchAddress("HLT_HIPhoton10_Photon15_v1_Prescl", &tHlts.HLT_HIPhoton10_Photon15_v1_Prescl, &tHlts.b_HLT_HIPhoton10_Photon15_v1_Prescl);
   if (t->GetBranch("HLT_HIPhoton15_Photon20_v1")) t->SetBranchAddress("HLT_HIPhoton15_Photon20_v1", &tHlts.HLT_HIPhoton15_Photon20_v1, &tHlts.b_HLT_HIPhoton15_Photon20_v1);
   if (t->GetBranch("HLT_HIPhoton15_Photon20_v1_Prescl")) t->SetBranchAddress("HLT_HIPhoton15_Photon20_v1_Prescl", &tHlts.HLT_HIPhoton15_Photon20_v1_Prescl, &tHlts.b_HLT_HIPhoton15_Photon20_v1_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton10_v1")) t->SetBranchAddress("HLT_HIDoublePhoton10_v1", &tHlts.HLT_HIDoublePhoton10_v1, &tHlts.b_HLT_HIDoublePhoton10_v1);
   if (t->GetBranch("HLT_HIDoublePhoton10_v1_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton10_v1_Prescl", &tHlts.HLT_HIDoublePhoton10_v1_Prescl, &tHlts.b_HLT_HIDoublePhoton10_v1_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton15_v1")) t->SetBranchAddress("HLT_HIDoublePhoton15_v1", &tHlts.HLT_HIDoublePhoton15_v1, &tHlts.b_HLT_HIDoublePhoton15_v1);
   if (t->GetBranch("HLT_HIDoublePhoton15_v1_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton15_v1_Prescl", &tHlts.HLT_HIDoublePhoton15_v1_Prescl, &tHlts.b_HLT_HIDoublePhoton15_v1_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton20_v1")) t->SetBranchAddress("HLT_HIDoublePhoton20_v1", &tHlts.HLT_HIDoublePhoton20_v1, &tHlts.b_HLT_HIDoublePhoton20_v1);
   if (t->GetBranch("HLT_HIDoublePhoton20_v1_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton20_v1_Prescl", &tHlts.HLT_HIDoublePhoton20_v1_Prescl, &tHlts.b_HLT_HIDoublePhoton20_v1_Prescl);
   if (t->GetBranch("HLT_HIJet55_v1")) t->SetBranchAddress("HLT_HIJet55_v1", &tHlts.HLT_HIJet55_v1, &tHlts.b_HLT_HIJet55_v1);
   if (t->GetBranch("HLT_HIJet55_v1_Prescl")) t->SetBranchAddress("HLT_HIJet55_v1_Prescl", &tHlts.HLT_HIJet55_v1_Prescl, &tHlts.b_HLT_HIJet55_v1_Prescl);
   if (t->GetBranch("HLT_HIJet65_v1")) t->SetBranchAddress("HLT_HIJet65_v1", &tHlts.HLT_HIJet65_v1, &tHlts.b_HLT_HIJet65_v1);
   if (t->GetBranch("HLT_HIJet65_v1_Prescl")) t->SetBranchAddress("HLT_HIJet65_v1_Prescl", &tHlts.HLT_HIJet65_v1_Prescl, &tHlts.b_HLT_HIJet65_v1_Prescl);
   if (t->GetBranch("HLT_Jet60_v1")) t->SetBranchAddress("HLT_Jet60_v1", &tHlts.HLT_Jet60_v1, &tHlts.b_HLT_Jet60_v1);
   if (t->GetBranch("HLT_PAJet20_NoJetID_v1")) t->SetBranchAddress("HLT_PAJet20_NoJetID_v1", &tHlts.HLT_PAJet20_NoJetID_v1, &tHlts.b_HLT_PAJet20_NoJetID_v1);
   if (t->GetBranch("HLT_PAJet40_NoJetID_v1")) t->SetBranchAddress("HLT_PAJet40_NoJetID_v1", &tHlts.HLT_PAJet40_NoJetID_v1, &tHlts.b_HLT_PAJet40_NoJetID_v1);
   if (t->GetBranch("HLT_PAJet60_NoJetID_v1")) t->SetBranchAddress("HLT_PAJet60_NoJetID_v1", &tHlts.HLT_PAJet60_NoJetID_v1, &tHlts.b_HLT_PAJet60_NoJetID_v1);
   if (t->GetBranch("HLT_PAJet80_NoJetID_v1")) t->SetBranchAddress("HLT_PAJet80_NoJetID_v1", &tHlts.HLT_PAJet80_NoJetID_v1, &tHlts.b_HLT_PAJet80_NoJetID_v1);
   if (t->GetBranch("HLT_PAJet100_NoJetID_v1")) t->SetBranchAddress("HLT_PAJet100_NoJetID_v1", &tHlts.HLT_PAJet100_NoJetID_v1, &tHlts.b_HLT_PAJet100_NoJetID_v1);
  if (t->GetBranch("HLT_PAZeroBias_v1"))  t->SetBranchAddress("HLT_PAZeroBias_v1", &tHlts.HLT_PAZeroBias_v1, &tHlts.b_HLT_PAZeroBias_v1);
  if (t->GetBranch("HLT_PAZeroBias_v1_Prescl"))  t->SetBranchAddress("HLT_PAZeroBias_v1_Prescl", &tHlts.HLT_PAZeroBias_v1_Prescl, &tHlts.b_HLT_PAZeroBias_v1_Prescl);
  if (t->GetBranch("HLT_PAZeroBiasPixel_SingleTrack_v1"))  t->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1", &tHlts.HLT_PAZeroBiasPixel_SingleTrack_v1, &tHlts.b_HLT_PAZeroBiasPixel_SingleTrack_v1);
   if (t->GetBranch("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl", &tHlts.HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl);
  if (t->GetBranch("HLT_PAFullTrack12_v2"))  t->SetBranchAddress("HLT_PAFullTrack12_v2", &tHlts.HLT_PAFullTrack12_v2, &tHlts.b_HLT_PAFullTrack12_v2);
  if (t->GetBranch("HLT_PAFullTrack20_v2"))  t->SetBranchAddress("HLT_PAFullTrack20_v2", &tHlts.HLT_PAFullTrack20_v2, &tHlts.b_HLT_PAFullTrack20_v2);
  if (t->GetBranch("HLT_PAFullTrack30_v2"))  t->SetBranchAddress("HLT_PAFullTrack30_v2", &tHlts.HLT_PAFullTrack30_v2, &tHlts.b_HLT_PAFullTrack30_v2);
   if (t->GetBranch("HLT_HIJet80_v1")) t->SetBranchAddress("HLT_HIJet80_v1", &tHlts.HLT_HIJet80_v1, &tHlts.b_HLT_HIJet80_v1);
   if (t->GetBranch("HLT_HIJet80_v1_Prescl")) t->SetBranchAddress("HLT_HIJet80_v1_Prescl", &tHlts.HLT_HIJet80_v1_Prescl, &tHlts.b_HLT_HIJet80_v1_Prescl);
   if (t->GetBranch("HLT_HIJet95_v1")) t->SetBranchAddress("HLT_HIJet95_v1", &tHlts.HLT_HIJet95_v1, &tHlts.b_HLT_HIJet95_v1);
   if (t->GetBranch("HLT_HIJet95_v1_Prescl")) t->SetBranchAddress("HLT_HIJet95_v1_Prescl", &tHlts.HLT_HIJet95_v1_Prescl, &tHlts.b_HLT_HIJet95_v1_Prescl);
   if (t->GetBranch("HLT_HIDiJet55_v1")) t->SetBranchAddress("HLT_HIDiJet55_v1", &tHlts.HLT_HIDiJet55_v1, &tHlts.b_HLT_HIDiJet55_v1);
   if (t->GetBranch("HLT_HIDiJet55_v1_Prescl")) t->SetBranchAddress("HLT_HIDiJet55_v1_Prescl", &tHlts.HLT_HIDiJet55_v1_Prescl, &tHlts.b_HLT_HIDiJet55_v1_Prescl);
   if (t->GetBranch("HLT_HIJet65_Jet55_v1")) t->SetBranchAddress("HLT_HIJet65_Jet55_v1", &tHlts.HLT_HIJet65_Jet55_v1, &tHlts.b_HLT_HIJet65_Jet55_v1);
   if (t->GetBranch("HLT_HIJet65_Jet55_v1_Prescl")) t->SetBranchAddress("HLT_HIJet65_Jet55_v1_Prescl", &tHlts.HLT_HIJet65_Jet55_v1_Prescl, &tHlts.b_HLT_HIJet65_Jet55_v1_Prescl);
   if (t->GetBranch("HLT_HIJetE30_NoBPTX_v1")) t->SetBranchAddress("HLT_HIJetE30_NoBPTX_v1", &tHlts.HLT_HIJetE30_NoBPTX_v1, &tHlts.b_HLT_HIJetE30_NoBPTX_v1);
   if (t->GetBranch("HLT_HIJetE30_NoBPTX_v1_Prescl")) t->SetBranchAddress("HLT_HIJetE30_NoBPTX_v1_Prescl", &tHlts.HLT_HIJetE30_NoBPTX_v1_Prescl, &tHlts.b_HLT_HIJetE30_NoBPTX_v1_Prescl);
   if (t->GetBranch("HLT_HIJetE50_NoBPTX3BX_NoHalo_v1")) t->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v1", &tHlts.HLT_HIJetE50_NoBPTX3BX_NoHalo_v1, &tHlts.b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v1);
   if (t->GetBranch("HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl")) t->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl", &tHlts.HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl, &tHlts.b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl);
   if (t->GetBranch("HLT_HIActivityHF_Coincidence3_v1")) t->SetBranchAddress("HLT_HIActivityHF_Coincidence3_v1", &tHlts.HLT_HIActivityHF_Coincidence3_v1, &tHlts.b_HLT_HIActivityHF_Coincidence3_v1);
   if (t->GetBranch("HLT_HIActivityHF_Coincidence3_v1_Prescl")) t->SetBranchAddress("HLT_HIActivityHF_Coincidence3_v1_Prescl", &tHlts.HLT_HIActivityHF_Coincidence3_v1_Prescl, &tHlts.b_HLT_HIActivityHF_Coincidence3_v1_Prescl);
   if (t->GetBranch("HLT_HIActivityHF_Single3_v1")) t->SetBranchAddress("HLT_HIActivityHF_Single3_v1", &tHlts.HLT_HIActivityHF_Single3_v1, &tHlts.b_HLT_HIActivityHF_Single3_v1);
   if (t->GetBranch("HLT_HIActivityHF_Single3_v1_Prescl")) t->SetBranchAddress("HLT_HIActivityHF_Single3_v1_Prescl", &tHlts.HLT_HIActivityHF_Single3_v1_Prescl, &tHlts.b_HLT_HIActivityHF_Single3_v1_Prescl);
   if (t->GetBranch("HLT_HIClusterVertexCompatibility_v1")) t->SetBranchAddress("HLT_HIClusterVertexCompatibility_v1", &tHlts.HLT_HIClusterVertexCompatibility_v1, &tHlts.b_HLT_HIClusterVertexCompatibility_v1);
   if (t->GetBranch("HLT_HIClusterVertexCompatibility_v1_Prescl")) t->SetBranchAddress("HLT_HIClusterVertexCompatibility_v1_Prescl", &tHlts.HLT_HIClusterVertexCompatibility_v1_Prescl, &tHlts.b_HLT_HIClusterVertexCompatibility_v1_Prescl);
   if (t->GetBranch("HLT_HICentralityVeto_v1")) t->SetBranchAddress("HLT_HICentralityVeto_v1", &tHlts.HLT_HICentralityVeto_v1, &tHlts.b_HLT_HICentralityVeto_v1);
   if (t->GetBranch("HLT_HICentralityVeto_v1_Prescl")) t->SetBranchAddress("HLT_HICentralityVeto_v1_Prescl", &tHlts.HLT_HICentralityVeto_v1_Prescl, &tHlts.b_HLT_HICentralityVeto_v1_Prescl);
   if (t->GetBranch("HLT_HIFullTrack12_L1Central_v1")) t->SetBranchAddress("HLT_HIFullTrack12_L1Central_v1", &tHlts.HLT_HIFullTrack12_L1Central_v1, &tHlts.b_HLT_HIFullTrack12_L1Central_v1);
   if (t->GetBranch("HLT_HIFullTrack12_L1Central_v1_Prescl")) t->SetBranchAddress("HLT_HIFullTrack12_L1Central_v1_Prescl", &tHlts.HLT_HIFullTrack12_L1Central_v1_Prescl, &tHlts.b_HLT_HIFullTrack12_L1Central_v1_Prescl);
   if (t->GetBranch("HLT_HIFullTrack12_L1Peripheral_v1")) t->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v1", &tHlts.HLT_HIFullTrack12_L1Peripheral_v1, &tHlts.b_HLT_HIFullTrack12_L1Peripheral_v1);
   if (t->GetBranch("HLT_HIFullTrack12_L1Peripheral_v1_Prescl")) t->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v1_Prescl", &tHlts.HLT_HIFullTrack12_L1Peripheral_v1_Prescl, &tHlts.b_HLT_HIFullTrack12_L1Peripheral_v1_Prescl);
   if (t->GetBranch("HLT_HIFullTrack14_L1Central_v1")) t->SetBranchAddress("HLT_HIFullTrack14_L1Central_v1", &tHlts.HLT_HIFullTrack14_L1Central_v1, &tHlts.b_HLT_HIFullTrack14_L1Central_v1);
   if (t->GetBranch("HLT_HIFullTrack14_L1Central_v1_Prescl")) t->SetBranchAddress("HLT_HIFullTrack14_L1Central_v1_Prescl", &tHlts.HLT_HIFullTrack14_L1Central_v1_Prescl, &tHlts.b_HLT_HIFullTrack14_L1Central_v1_Prescl);
   if (t->GetBranch("HLT_HIFullTrack14_L1Peripheral_v1")) t->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v1", &tHlts.HLT_HIFullTrack14_L1Peripheral_v1, &tHlts.b_HLT_HIFullTrack14_L1Peripheral_v1);
   if (t->GetBranch("HLT_HIFullTrack14_L1Peripheral_v1_Prescl")) t->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v1_Prescl", &tHlts.HLT_HIFullTrack14_L1Peripheral_v1_Prescl, &tHlts.b_HLT_HIFullTrack14_L1Peripheral_v1_Prescl);
   if (t->GetBranch("HLT_HIFullTrack20_L1Central_v1")) t->SetBranchAddress("HLT_HIFullTrack20_L1Central_v1", &tHlts.HLT_HIFullTrack20_L1Central_v1, &tHlts.b_HLT_HIFullTrack20_L1Central_v1);
   if (t->GetBranch("HLT_HIFullTrack20_L1Central_v1_Prescl")) t->SetBranchAddress("HLT_HIFullTrack20_L1Central_v1_Prescl", &tHlts.HLT_HIFullTrack20_L1Central_v1_Prescl, &tHlts.b_HLT_HIFullTrack20_L1Central_v1_Prescl);
   if (t->GetBranch("HLT_HIFullTrack20_L1Peripheral_v1")) t->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v1", &tHlts.HLT_HIFullTrack20_L1Peripheral_v1, &tHlts.b_HLT_HIFullTrack20_L1Peripheral_v1);
   if (t->GetBranch("HLT_HIFullTrack20_L1Peripheral_v1_Prescl")) t->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v1_Prescl", &tHlts.HLT_HIFullTrack20_L1Peripheral_v1_Prescl, &tHlts.b_HLT_HIFullTrack20_L1Peripheral_v1_Prescl);
   if (t->GetBranch("HLT_HIFullTrack25_L1Central_v1")) t->SetBranchAddress("HLT_HIFullTrack25_L1Central_v1", &tHlts.HLT_HIFullTrack25_L1Central_v1, &tHlts.b_HLT_HIFullTrack25_L1Central_v1);
   if (t->GetBranch("HLT_HIFullTrack25_L1Central_v1_Prescl")) t->SetBranchAddress("HLT_HIFullTrack25_L1Central_v1_Prescl", &tHlts.HLT_HIFullTrack25_L1Central_v1_Prescl, &tHlts.b_HLT_HIFullTrack25_L1Central_v1_Prescl);
   if (t->GetBranch("HLT_HIFullTrack25_L1Peripheral_v1")) t->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v1", &tHlts.HLT_HIFullTrack25_L1Peripheral_v1, &tHlts.b_HLT_HIFullTrack25_L1Peripheral_v1);
   if (t->GetBranch("HLT_HIFullTrack25_L1Peripheral_v1_Prescl")) t->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v1_Prescl", &tHlts.HLT_HIFullTrack25_L1Peripheral_v1_Prescl, &tHlts.b_HLT_HIFullTrack25_L1Peripheral_v1_Prescl);
   if (t->GetBranch("HLT_HIRandom_v1")) t->SetBranchAddress("HLT_HIRandom_v1", &tHlts.HLT_HIRandom_v1, &tHlts.b_HLT_HIRandom_v1);
   if (t->GetBranch("HLT_HIRandom_v1_Prescl")) t->SetBranchAddress("HLT_HIRandom_v1_Prescl", &tHlts.HLT_HIRandom_v1_Prescl, &tHlts.b_HLT_HIRandom_v1_Prescl);
   if (t->GetBranch("HLT_HIUCC010_v1")) t->SetBranchAddress("HLT_HIUCC010_v1", &tHlts.HLT_HIUCC010_v1, &tHlts.b_HLT_HIUCC010_v1);
   if (t->GetBranch("HLT_HIUCC010_v1_Prescl")) t->SetBranchAddress("HLT_HIUCC010_v1_Prescl", &tHlts.HLT_HIUCC010_v1_Prescl, &tHlts.b_HLT_HIUCC010_v1_Prescl);
   if (t->GetBranch("HLT_HIUCC015_v1")) t->SetBranchAddress("HLT_HIUCC015_v1", &tHlts.HLT_HIUCC015_v1, &tHlts.b_HLT_HIUCC015_v1);
   if (t->GetBranch("HLT_HIUCC015_v1_Prescl")) t->SetBranchAddress("HLT_HIUCC015_v1_Prescl", &tHlts.HLT_HIUCC015_v1_Prescl, &tHlts.b_HLT_HIUCC015_v1_Prescl);
   if (t->GetBranch("HLT_HICentral10_v1")) t->SetBranchAddress("HLT_HICentral10_v1", &tHlts.HLT_HICentral10_v1, &tHlts.b_HLT_HICentral10_v1);
   if (t->GetBranch("HLT_HICentral10_v1_Prescl")) t->SetBranchAddress("HLT_HICentral10_v1_Prescl", &tHlts.HLT_HICentral10_v1_Prescl, &tHlts.b_HLT_HICentral10_v1_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuMuPixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v1", &tHlts.HLT_HIUPCNeuMuPixel_SingleTrack_v1, &tHlts.b_HLT_HIUPCNeuMuPixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuEG2Pixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v1", &tHlts.HLT_HIUPCNeuEG2Pixel_SingleTrack_v1, &tHlts.b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuEG5Pixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v1", &tHlts.HLT_HIUPCNeuEG5Pixel_SingleTrack_v1, &tHlts.b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1", &tHlts.HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1, &tHlts.b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1", &tHlts.HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1, &tHlts.b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1", &tHlts.HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1, &tHlts.b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl", &tHlts.HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl);
   if (t->GetBranch("HLTriggerFinalPath")) t->SetBranchAddress("HLTriggerFinalPath", &tHlts.HLTriggerFinalPath, &tHlts.b_HLTriggerFinalPath);
   if (t->GetBranch("HLTriggerFinalPath_Prescl")) t->SetBranchAddress("HLTriggerFinalPath_Prescl", &tHlts.HLTriggerFinalPath_Prescl, &tHlts.b_HLTriggerFinalPath_Prescl);
   if (t->GetBranch("HLT_HIMET120_v1")) t->SetBranchAddress("HLT_HIMET120_v1", &tHlts.HLT_HIMET120_v1, &tHlts.b_HLT_HIMET120_v1);
   if (t->GetBranch("HLT_HIMET120_v1_Prescl")) t->SetBranchAddress("HLT_HIMET120_v1_Prescl", &tHlts.HLT_HIMET120_v1_Prescl, &tHlts.b_HLT_HIMET120_v1_Prescl);
   if (t->GetBranch("HLT_HIMET200_v1")) t->SetBranchAddress("HLT_HIMET200_v1", &tHlts.HLT_HIMET200_v1, &tHlts.b_HLT_HIMET200_v1);
   if (t->GetBranch("HLT_HIMET200_v1_Prescl")) t->SetBranchAddress("HLT_HIMET200_v1_Prescl", &tHlts.HLT_HIMET200_v1_Prescl, &tHlts.b_HLT_HIMET200_v1_Prescl);
   if (t->GetBranch("HLT_HIMET220_v1")) t->SetBranchAddress("HLT_HIMET220_v1", &tHlts.HLT_HIMET220_v1, &tHlts.b_HLT_HIMET220_v1);
   if (t->GetBranch("HLT_HIMET220_v1_Prescl")) t->SetBranchAddress("HLT_HIMET220_v1_Prescl", &tHlts.HLT_HIMET220_v1_Prescl, &tHlts.b_HLT_HIMET220_v1_Prescl);
   if (t->GetBranch("HLT_HIPhysics_v1")) t->SetBranchAddress("HLT_HIPhysics_v1", &tHlts.HLT_HIPhysics_v1, &tHlts.b_HLT_HIPhysics_v1);
   if (t->GetBranch("HLT_HIPhysics_v1_Prescl")) t->SetBranchAddress("HLT_HIPhysics_v1_Prescl", &tHlts.HLT_HIPhysics_v1_Prescl, &tHlts.b_HLT_HIPhysics_v1_Prescl);
   if (t->GetBranch("HLT_HIL1DoubleMu0_HighQ_v2")) t->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v2", &tHlts.HLT_HIL1DoubleMu0_HighQ_v2, &tHlts.b_HLT_HIL1DoubleMu0_HighQ_v2);
   if (t->GetBranch("HLT_HIL1DoubleMu0_HighQ_v2_Prescl")) t->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v2_Prescl", &tHlts.HLT_HIL1DoubleMu0_HighQ_v2_Prescl, &tHlts.b_HLT_HIL1DoubleMu0_HighQ_v2_Prescl);
   if (t->GetBranch("HLT_HIL1DoubleMuOpen_v2")) t->SetBranchAddress("HLT_HIL1DoubleMuOpen_v2", &tHlts.HLT_HIL1DoubleMuOpen_v2, &tHlts.b_HLT_HIL1DoubleMuOpen_v2);
   if (t->GetBranch("HLT_HIL1DoubleMuOpen_v2_Prescl")) t->SetBranchAddress("HLT_HIL1DoubleMuOpen_v2_Prescl", &tHlts.HLT_HIL1DoubleMuOpen_v2_Prescl, &tHlts.b_HLT_HIL1DoubleMuOpen_v2_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2")) t->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2", &tHlts.HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2, &tHlts.b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2);
   if (t->GetBranch("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl", &tHlts.HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_NHitQ_v2")) t->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v2", &tHlts.HLT_HIL2DoubleMu0_NHitQ_v2, &tHlts.b_HLT_HIL2DoubleMu0_NHitQ_v2);
   if (t->GetBranch("HLT_HIL2DoubleMu0_NHitQ_v2_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v2_Prescl", &tHlts.HLT_HIL2DoubleMu0_NHitQ_v2_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_NHitQ_v2_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_v2")) t->SetBranchAddress("HLT_HIL2DoubleMu0_v2", &tHlts.HLT_HIL2DoubleMu0_v2, &tHlts.b_HLT_HIL2DoubleMu0_v2);
   if (t->GetBranch("HLT_HIL2DoubleMu0_v2_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_v2_Prescl", &tHlts.HLT_HIL2DoubleMu0_v2_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_v2_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu3_v2")) t->SetBranchAddress("HLT_HIL2DoubleMu3_v2", &tHlts.HLT_HIL2DoubleMu3_v2, &tHlts.b_HLT_HIL2DoubleMu3_v2);
   if (t->GetBranch("HLT_HIL2DoubleMu3_v2_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu3_v2_Prescl", &tHlts.HLT_HIL2DoubleMu3_v2_Prescl, &tHlts.b_HLT_HIL2DoubleMu3_v2_Prescl);
   if (t->GetBranch("HLT_HIL2Mu15_v2")) t->SetBranchAddress("HLT_HIL2Mu15_v2", &tHlts.HLT_HIL2Mu15_v2, &tHlts.b_HLT_HIL2Mu15_v2);
   if (t->GetBranch("HLT_HIL2Mu15_v2_Prescl")) t->SetBranchAddress("HLT_HIL2Mu15_v2_Prescl", &tHlts.HLT_HIL2Mu15_v2_Prescl, &tHlts.b_HLT_HIL2Mu15_v2_Prescl);
   if (t->GetBranch("HLT_HIL2Mu3_NHitQ_v2")) t->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v2", &tHlts.HLT_HIL2Mu3_NHitQ_v2, &tHlts.b_HLT_HIL2Mu3_NHitQ_v2);
   if (t->GetBranch("HLT_HIL2Mu3_NHitQ_v2_Prescl")) t->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v2_Prescl", &tHlts.HLT_HIL2Mu3_NHitQ_v2_Prescl, &tHlts.b_HLT_HIL2Mu3_NHitQ_v2_Prescl);
   if (t->GetBranch("HLT_HIL2Mu3_v2")) t->SetBranchAddress("HLT_HIL2Mu3_v2", &tHlts.HLT_HIL2Mu3_v2, &tHlts.b_HLT_HIL2Mu3_v2);
   if (t->GetBranch("HLT_HIL2Mu3_v2_Prescl")) t->SetBranchAddress("HLT_HIL2Mu3_v2_Prescl", &tHlts.HLT_HIL2Mu3_v2_Prescl, &tHlts.b_HLT_HIL2Mu3_v2_Prescl);
   if (t->GetBranch("HLT_HIL2Mu7_v2")) t->SetBranchAddress("HLT_HIL2Mu7_v2", &tHlts.HLT_HIL2Mu7_v2, &tHlts.b_HLT_HIL2Mu7_v2);
   if (t->GetBranch("HLT_HIL2Mu7_v2_Prescl")) t->SetBranchAddress("HLT_HIL2Mu7_v2_Prescl", &tHlts.HLT_HIL2Mu7_v2_Prescl, &tHlts.b_HLT_HIL2Mu7_v2_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_v2")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v2", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_v2, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v2);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_SS_v2")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v2", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_SS_v2, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v2);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_v2")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v2", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_v2, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_v2);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_v2")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_v2", &tHlts.HLT_HIL3DoubleMuOpen_v2, &tHlts.b_HLT_HIL3DoubleMuOpen_v2);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_v2_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_v2_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_v2_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_v2_Prescl);
   if (t->GetBranch("HLT_HIL3Mu3_v2")) t->SetBranchAddress("HLT_HIL3Mu3_v2", &tHlts.HLT_HIL3Mu3_v2, &tHlts.b_HLT_HIL3Mu3_v2);
   if (t->GetBranch("HLT_HIL3Mu3_v2_Prescl")) t->SetBranchAddress("HLT_HIL3Mu3_v2_Prescl", &tHlts.HLT_HIL3Mu3_v2_Prescl, &tHlts.b_HLT_HIL3Mu3_v2_Prescl);
   if (t->GetBranch("HLT_HIDiJet55_v2")) t->SetBranchAddress("HLT_HIDiJet55_v2", &tHlts.HLT_HIDiJet55_v2, &tHlts.b_HLT_HIDiJet55_v2);
   if (t->GetBranch("HLT_HIDiJet55_v2_Prescl")) t->SetBranchAddress("HLT_HIDiJet55_v2_Prescl", &tHlts.HLT_HIDiJet55_v2_Prescl, &tHlts.b_HLT_HIDiJet55_v2_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton10_v2")) t->SetBranchAddress("HLT_HIDoublePhoton10_v2", &tHlts.HLT_HIDoublePhoton10_v2, &tHlts.b_HLT_HIDoublePhoton10_v2);
   if (t->GetBranch("HLT_HIDoublePhoton10_v2_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton10_v2_Prescl", &tHlts.HLT_HIDoublePhoton10_v2_Prescl, &tHlts.b_HLT_HIDoublePhoton10_v2_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton15_v2")) t->SetBranchAddress("HLT_HIDoublePhoton15_v2", &tHlts.HLT_HIDoublePhoton15_v2, &tHlts.b_HLT_HIDoublePhoton15_v2);
   if (t->GetBranch("HLT_HIDoublePhoton15_v2_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton15_v2_Prescl", &tHlts.HLT_HIDoublePhoton15_v2_Prescl, &tHlts.b_HLT_HIDoublePhoton15_v2_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton20_v2")) t->SetBranchAddress("HLT_HIDoublePhoton20_v2", &tHlts.HLT_HIDoublePhoton20_v2, &tHlts.b_HLT_HIDoublePhoton20_v2);
   if (t->GetBranch("HLT_HIDoublePhoton20_v2_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton20_v2_Prescl", &tHlts.HLT_HIDoublePhoton20_v2_Prescl, &tHlts.b_HLT_HIDoublePhoton20_v2_Prescl);
   if (t->GetBranch("HLT_HIFullTrack12_L1Central_v2")) t->SetBranchAddress("HLT_HIFullTrack12_L1Central_v2", &tHlts.HLT_HIFullTrack12_L1Central_v2, &tHlts.b_HLT_HIFullTrack12_L1Central_v2);
   if (t->GetBranch("HLT_HIFullTrack12_L1Central_v2_Prescl")) t->SetBranchAddress("HLT_HIFullTrack12_L1Central_v2_Prescl", &tHlts.HLT_HIFullTrack12_L1Central_v2_Prescl, &tHlts.b_HLT_HIFullTrack12_L1Central_v2_Prescl);
   if (t->GetBranch("HLT_HIFullTrack12_L1Peripheral_v2")) t->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v2", &tHlts.HLT_HIFullTrack12_L1Peripheral_v2, &tHlts.b_HLT_HIFullTrack12_L1Peripheral_v2);
   if (t->GetBranch("HLT_HIFullTrack12_L1Peripheral_v2_Prescl")) t->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v2_Prescl", &tHlts.HLT_HIFullTrack12_L1Peripheral_v2_Prescl, &tHlts.b_HLT_HIFullTrack12_L1Peripheral_v2_Prescl);
   if (t->GetBranch("HLT_HIFullTrack14_L1Central_v2")) t->SetBranchAddress("HLT_HIFullTrack14_L1Central_v2", &tHlts.HLT_HIFullTrack14_L1Central_v2, &tHlts.b_HLT_HIFullTrack14_L1Central_v2);
   if (t->GetBranch("HLT_HIFullTrack14_L1Central_v2_Prescl")) t->SetBranchAddress("HLT_HIFullTrack14_L1Central_v2_Prescl", &tHlts.HLT_HIFullTrack14_L1Central_v2_Prescl, &tHlts.b_HLT_HIFullTrack14_L1Central_v2_Prescl);
   if (t->GetBranch("HLT_HIFullTrack14_L1Peripheral_v2")) t->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v2", &tHlts.HLT_HIFullTrack14_L1Peripheral_v2, &tHlts.b_HLT_HIFullTrack14_L1Peripheral_v2);
   if (t->GetBranch("HLT_HIFullTrack14_L1Peripheral_v2_Prescl")) t->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v2_Prescl", &tHlts.HLT_HIFullTrack14_L1Peripheral_v2_Prescl, &tHlts.b_HLT_HIFullTrack14_L1Peripheral_v2_Prescl);
   if (t->GetBranch("HLT_HIFullTrack20_L1Central_v2")) t->SetBranchAddress("HLT_HIFullTrack20_L1Central_v2", &tHlts.HLT_HIFullTrack20_L1Central_v2, &tHlts.b_HLT_HIFullTrack20_L1Central_v2);
   if (t->GetBranch("HLT_HIFullTrack20_L1Central_v2_Prescl")) t->SetBranchAddress("HLT_HIFullTrack20_L1Central_v2_Prescl", &tHlts.HLT_HIFullTrack20_L1Central_v2_Prescl, &tHlts.b_HLT_HIFullTrack20_L1Central_v2_Prescl);
   if (t->GetBranch("HLT_HIFullTrack20_L1Peripheral_v2")) t->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v2", &tHlts.HLT_HIFullTrack20_L1Peripheral_v2, &tHlts.b_HLT_HIFullTrack20_L1Peripheral_v2);
   if (t->GetBranch("HLT_HIFullTrack20_L1Peripheral_v2_Prescl")) t->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v2_Prescl", &tHlts.HLT_HIFullTrack20_L1Peripheral_v2_Prescl, &tHlts.b_HLT_HIFullTrack20_L1Peripheral_v2_Prescl);
   if (t->GetBranch("HLT_HIFullTrack25_L1Central_v2")) t->SetBranchAddress("HLT_HIFullTrack25_L1Central_v2", &tHlts.HLT_HIFullTrack25_L1Central_v2, &tHlts.b_HLT_HIFullTrack25_L1Central_v2);
   if (t->GetBranch("HLT_HIFullTrack25_L1Central_v2_Prescl")) t->SetBranchAddress("HLT_HIFullTrack25_L1Central_v2_Prescl", &tHlts.HLT_HIFullTrack25_L1Central_v2_Prescl, &tHlts.b_HLT_HIFullTrack25_L1Central_v2_Prescl);
   if (t->GetBranch("HLT_HIFullTrack25_L1Peripheral_v2")) t->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v2", &tHlts.HLT_HIFullTrack25_L1Peripheral_v2, &tHlts.b_HLT_HIFullTrack25_L1Peripheral_v2);
   if (t->GetBranch("HLT_HIFullTrack25_L1Peripheral_v2_Prescl")) t->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v2_Prescl", &tHlts.HLT_HIFullTrack25_L1Peripheral_v2_Prescl, &tHlts.b_HLT_HIFullTrack25_L1Peripheral_v2_Prescl);
   if (t->GetBranch("HLT_HIJet55_v2")) t->SetBranchAddress("HLT_HIJet55_v2", &tHlts.HLT_HIJet55_v2, &tHlts.b_HLT_HIJet55_v2);
   if (t->GetBranch("HLT_HIJet55_v2_Prescl")) t->SetBranchAddress("HLT_HIJet55_v2_Prescl", &tHlts.HLT_HIJet55_v2_Prescl, &tHlts.b_HLT_HIJet55_v2_Prescl);
   if (t->GetBranch("HLT_HIJet65_Jet55_v2")) t->SetBranchAddress("HLT_HIJet65_Jet55_v2", &tHlts.HLT_HIJet65_Jet55_v2, &tHlts.b_HLT_HIJet65_Jet55_v2);
   if (t->GetBranch("HLT_HIJet65_Jet55_v2_Prescl")) t->SetBranchAddress("HLT_HIJet65_Jet55_v2_Prescl", &tHlts.HLT_HIJet65_Jet55_v2_Prescl, &tHlts.b_HLT_HIJet65_Jet55_v2_Prescl);
   if (t->GetBranch("HLT_HIJet65_v2")) t->SetBranchAddress("HLT_HIJet65_v2", &tHlts.HLT_HIJet65_v2, &tHlts.b_HLT_HIJet65_v2);
   if (t->GetBranch("HLT_HIJet65_v2_Prescl")) t->SetBranchAddress("HLT_HIJet65_v2_Prescl", &tHlts.HLT_HIJet65_v2_Prescl, &tHlts.b_HLT_HIJet65_v2_Prescl);
   if (t->GetBranch("HLT_HIJet80_v2")) t->SetBranchAddress("HLT_HIJet80_v2", &tHlts.HLT_HIJet80_v2, &tHlts.b_HLT_HIJet80_v2);
   if (t->GetBranch("HLT_HIJet80_v2_Prescl")) t->SetBranchAddress("HLT_HIJet80_v2_Prescl", &tHlts.HLT_HIJet80_v2_Prescl, &tHlts.b_HLT_HIJet80_v2_Prescl);
   if (t->GetBranch("HLT_HIJet95_v2")) t->SetBranchAddress("HLT_HIJet95_v2", &tHlts.HLT_HIJet95_v2, &tHlts.b_HLT_HIJet95_v2);
   if (t->GetBranch("HLT_HIJet95_v2_Prescl")) t->SetBranchAddress("HLT_HIJet95_v2_Prescl", &tHlts.HLT_HIJet95_v2_Prescl, &tHlts.b_HLT_HIJet95_v2_Prescl);
   if (t->GetBranch("HLT_HIJetE30_NoBPTX_v2")) t->SetBranchAddress("HLT_HIJetE30_NoBPTX_v2", &tHlts.HLT_HIJetE30_NoBPTX_v2, &tHlts.b_HLT_HIJetE30_NoBPTX_v2);
   if (t->GetBranch("HLT_HIJetE30_NoBPTX_v2_Prescl")) t->SetBranchAddress("HLT_HIJetE30_NoBPTX_v2_Prescl", &tHlts.HLT_HIJetE30_NoBPTX_v2_Prescl, &tHlts.b_HLT_HIJetE30_NoBPTX_v2_Prescl);
   if (t->GetBranch("HLT_HIJetE50_NoBPTX3BX_NoHalo_v2")) t->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v2", &tHlts.HLT_HIJetE50_NoBPTX3BX_NoHalo_v2, &tHlts.b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v2);
   if (t->GetBranch("HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl")) t->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl", &tHlts.HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl, &tHlts.b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl);
   if (t->GetBranch("HLT_HIMET120_v2")) t->SetBranchAddress("HLT_HIMET120_v2", &tHlts.HLT_HIMET120_v2, &tHlts.b_HLT_HIMET120_v2);
   if (t->GetBranch("HLT_HIMET120_v2_Prescl")) t->SetBranchAddress("HLT_HIMET120_v2_Prescl", &tHlts.HLT_HIMET120_v2_Prescl, &tHlts.b_HLT_HIMET120_v2_Prescl);
   if (t->GetBranch("HLT_HIMET200_v2")) t->SetBranchAddress("HLT_HIMET200_v2", &tHlts.HLT_HIMET200_v2, &tHlts.b_HLT_HIMET200_v2);
   if (t->GetBranch("HLT_HIMET200_v2_Prescl")) t->SetBranchAddress("HLT_HIMET200_v2_Prescl", &tHlts.HLT_HIMET200_v2_Prescl, &tHlts.b_HLT_HIMET200_v2_Prescl);
   if (t->GetBranch("HLT_HIMET220_v2")) t->SetBranchAddress("HLT_HIMET220_v2", &tHlts.HLT_HIMET220_v2, &tHlts.b_HLT_HIMET220_v2);
   if (t->GetBranch("HLT_HIMET220_v2_Prescl")) t->SetBranchAddress("HLT_HIMET220_v2_Prescl", &tHlts.HLT_HIMET220_v2_Prescl, &tHlts.b_HLT_HIMET220_v2_Prescl);
   if (t->GetBranch("HLT_HIPhoton10_Photon15_v2")) t->SetBranchAddress("HLT_HIPhoton10_Photon15_v2", &tHlts.HLT_HIPhoton10_Photon15_v2, &tHlts.b_HLT_HIPhoton10_Photon15_v2);
   if (t->GetBranch("HLT_HIPhoton10_Photon15_v2_Prescl")) t->SetBranchAddress("HLT_HIPhoton10_Photon15_v2_Prescl", &tHlts.HLT_HIPhoton10_Photon15_v2_Prescl, &tHlts.b_HLT_HIPhoton10_Photon15_v2_Prescl);
   if (t->GetBranch("HLT_HIPhoton15_Photon20_v2")) t->SetBranchAddress("HLT_HIPhoton15_Photon20_v2", &tHlts.HLT_HIPhoton15_Photon20_v2, &tHlts.b_HLT_HIPhoton15_Photon20_v2);
   if (t->GetBranch("HLT_HIPhoton15_Photon20_v2_Prescl")) t->SetBranchAddress("HLT_HIPhoton15_Photon20_v2_Prescl", &tHlts.HLT_HIPhoton15_Photon20_v2_Prescl, &tHlts.b_HLT_HIPhoton15_Photon20_v2_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton15_v2")) t->SetBranchAddress("HLT_HISinglePhoton15_v2", &tHlts.HLT_HISinglePhoton15_v2, &tHlts.b_HLT_HISinglePhoton15_v2);
   if (t->GetBranch("HLT_HISinglePhoton15_v2_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton15_v2_Prescl", &tHlts.HLT_HISinglePhoton15_v2_Prescl, &tHlts.b_HLT_HISinglePhoton15_v2_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton20_v2")) t->SetBranchAddress("HLT_HISinglePhoton20_v2", &tHlts.HLT_HISinglePhoton20_v2, &tHlts.b_HLT_HISinglePhoton20_v2);
   if (t->GetBranch("HLT_HISinglePhoton20_v2_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton20_v2_Prescl", &tHlts.HLT_HISinglePhoton20_v2_Prescl, &tHlts.b_HLT_HISinglePhoton20_v2_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton30_v2")) t->SetBranchAddress("HLT_HISinglePhoton30_v2", &tHlts.HLT_HISinglePhoton30_v2, &tHlts.b_HLT_HISinglePhoton30_v2);
   if (t->GetBranch("HLT_HISinglePhoton30_v2_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton30_v2_Prescl", &tHlts.HLT_HISinglePhoton30_v2_Prescl, &tHlts.b_HLT_HISinglePhoton30_v2_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton40_v2")) t->SetBranchAddress("HLT_HISinglePhoton40_v2", &tHlts.HLT_HISinglePhoton40_v2, &tHlts.b_HLT_HISinglePhoton40_v2);
   if (t->GetBranch("HLT_HISinglePhoton40_v2_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton40_v2_Prescl", &tHlts.HLT_HISinglePhoton40_v2_Prescl, &tHlts.b_HLT_HISinglePhoton40_v2_Prescl);
   if (t->GetBranch("HLT_HIBptxXOR_v2")) t->SetBranchAddress("HLT_HIBptxXOR_v2", &tHlts.HLT_HIBptxXOR_v2, &tHlts.b_HLT_HIBptxXOR_v2);
   if (t->GetBranch("HLT_HIBptxXOR_v2_Prescl")) t->SetBranchAddress("HLT_HIBptxXOR_v2_Prescl", &tHlts.HLT_HIBptxXOR_v2_Prescl, &tHlts.b_HLT_HIBptxXOR_v2_Prescl);
   if (t->GetBranch("HLT_HICentral10_v2")) t->SetBranchAddress("HLT_HICentral10_v2", &tHlts.HLT_HICentral10_v2, &tHlts.b_HLT_HICentral10_v2);
   if (t->GetBranch("HLT_HICentral10_v2_Prescl")) t->SetBranchAddress("HLT_HICentral10_v2_Prescl", &tHlts.HLT_HICentral10_v2_Prescl, &tHlts.b_HLT_HICentral10_v2_Prescl);
   if (t->GetBranch("HLT_HICentralityVeto_v2")) t->SetBranchAddress("HLT_HICentralityVeto_v2", &tHlts.HLT_HICentralityVeto_v2, &tHlts.b_HLT_HICentralityVeto_v2);
   if (t->GetBranch("HLT_HICentralityVeto_v2_Prescl")) t->SetBranchAddress("HLT_HICentralityVeto_v2_Prescl", &tHlts.HLT_HICentralityVeto_v2_Prescl, &tHlts.b_HLT_HICentralityVeto_v2_Prescl);
   if (t->GetBranch("HLT_HIL1Algo_BptxXOR_BSC_OR_v2")) t->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v2", &tHlts.HLT_HIL1Algo_BptxXOR_BSC_OR_v2, &tHlts.b_HLT_HIL1Algo_BptxXOR_BSC_OR_v2);
   if (t->GetBranch("HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl")) t->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl", &tHlts.HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl, &tHlts.b_HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasBSC_OR_v2")) t->SetBranchAddress("HLT_HIMinBiasBSC_OR_v2", &tHlts.HLT_HIMinBiasBSC_OR_v2, &tHlts.b_HLT_HIMinBiasBSC_OR_v2);
   if (t->GetBranch("HLT_HIMinBiasBSC_OR_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasBSC_OR_v2_Prescl", &tHlts.HLT_HIMinBiasBSC_OR_v2_Prescl, &tHlts.b_HLT_HIMinBiasBSC_OR_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasBSC_v2")) t->SetBranchAddress("HLT_HIMinBiasBSC_v2", &tHlts.HLT_HIMinBiasBSC_v2, &tHlts.b_HLT_HIMinBiasBSC_v2);
   if (t->GetBranch("HLT_HIMinBiasBSC_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasBSC_v2_Prescl", &tHlts.HLT_HIMinBiasBSC_v2_Prescl, &tHlts.b_HLT_HIMinBiasBSC_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHF_v2")) t->SetBranchAddress("HLT_HIMinBiasHF_v2", &tHlts.HLT_HIMinBiasHF_v2, &tHlts.b_HLT_HIMinBiasHF_v2);
   if (t->GetBranch("HLT_HIMinBiasHF_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHF_v2_Prescl", &tHlts.HLT_HIMinBiasHF_v2_Prescl, &tHlts.b_HLT_HIMinBiasHF_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHfOrBSC_v2")) t->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v2", &tHlts.HLT_HIMinBiasHfOrBSC_v2, &tHlts.b_HLT_HIMinBiasHfOrBSC_v2);
   if (t->GetBranch("HLT_HIMinBiasHfOrBSC_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v2_Prescl", &tHlts.HLT_HIMinBiasHfOrBSC_v2_Prescl, &tHlts.b_HLT_HIMinBiasHfOrBSC_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHf_OR_v2")) t->SetBranchAddress("HLT_HIMinBiasHf_OR_v2", &tHlts.HLT_HIMinBiasHf_OR_v2, &tHlts.b_HLT_HIMinBiasHf_OR_v2);
   if (t->GetBranch("HLT_HIMinBiasHf_OR_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHf_OR_v2_Prescl", &tHlts.HLT_HIMinBiasHf_OR_v2_Prescl, &tHlts.b_HLT_HIMinBiasHf_OR_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasPixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v2", &tHlts.HLT_HIMinBiasPixel_SingleTrack_v2, &tHlts.b_HLT_HIMinBiasPixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIMinBiasPixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIMinBiasPixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIMinBiasPixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDCPixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v2", &tHlts.HLT_HIMinBiasZDCPixel_SingleTrack_v2, &tHlts.b_HLT_HIMinBiasZDCPixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2", &tHlts.HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2, &tHlts.b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl", &tHlts.HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl, &tHlts.b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_v2")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v2", &tHlts.HLT_HIMinBiasZDC_Calo_v2, &tHlts.b_HLT_HIMinBiasZDC_Calo_v2);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v2_Prescl", &tHlts.HLT_HIMinBiasZDC_Calo_v2_Prescl, &tHlts.b_HLT_HIMinBiasZDC_Calo_v2_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2", &tHlts.HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2, &tHlts.b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIPhysics_v2")) t->SetBranchAddress("HLT_HIPhysics_v2", &tHlts.HLT_HIPhysics_v2, &tHlts.b_HLT_HIPhysics_v2);
   if (t->GetBranch("HLT_HIPhysics_v2_Prescl")) t->SetBranchAddress("HLT_HIPhysics_v2_Prescl", &tHlts.HLT_HIPhysics_v2_Prescl, &tHlts.b_HLT_HIPhysics_v2_Prescl);
   if (t->GetBranch("HLT_HIRandom_v2")) t->SetBranchAddress("HLT_HIRandom_v2", &tHlts.HLT_HIRandom_v2, &tHlts.b_HLT_HIRandom_v2);
   if (t->GetBranch("HLT_HIRandom_v2_Prescl")) t->SetBranchAddress("HLT_HIRandom_v2_Prescl", &tHlts.HLT_HIRandom_v2_Prescl, &tHlts.b_HLT_HIRandom_v2_Prescl);
   if (t->GetBranch("HLT_HIUCC010_v2")) t->SetBranchAddress("HLT_HIUCC010_v2", &tHlts.HLT_HIUCC010_v2, &tHlts.b_HLT_HIUCC010_v2);
   if (t->GetBranch("HLT_HIUCC010_v2_Prescl")) t->SetBranchAddress("HLT_HIUCC010_v2_Prescl", &tHlts.HLT_HIUCC010_v2_Prescl, &tHlts.b_HLT_HIUCC010_v2_Prescl);
   if (t->GetBranch("HLT_HIUCC015_v2")) t->SetBranchAddress("HLT_HIUCC015_v2", &tHlts.HLT_HIUCC015_v2, &tHlts.b_HLT_HIUCC015_v2);
   if (t->GetBranch("HLT_HIUCC015_v2_Prescl")) t->SetBranchAddress("HLT_HIUCC015_v2_Prescl", &tHlts.HLT_HIUCC015_v2_Prescl, &tHlts.b_HLT_HIUCC015_v2_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuEG2Pixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v2", &tHlts.HLT_HIUPCNeuEG2Pixel_SingleTrack_v2, &tHlts.b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuEG5Pixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v2", &tHlts.HLT_HIUPCNeuEG5Pixel_SingleTrack_v2, &tHlts.b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2", &tHlts.HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2, &tHlts.b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2", &tHlts.HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2, &tHlts.b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2", &tHlts.HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2, &tHlts.b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuMuPixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v2", &tHlts.HLT_HIUPCNeuMuPixel_SingleTrack_v2, &tHlts.b_HLT_HIUPCNeuMuPixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIZeroBiasPixel_SingleTrack_v2")) t->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v2", &tHlts.HLT_HIZeroBiasPixel_SingleTrack_v2, &tHlts.b_HLT_HIZeroBiasPixel_SingleTrack_v2);
   if (t->GetBranch("HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl")) t->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl", &tHlts.HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl, &tHlts.b_HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl);
   if (t->GetBranch("HLT_HIZeroBiasXOR_v2")) t->SetBranchAddress("HLT_HIZeroBiasXOR_v2", &tHlts.HLT_HIZeroBiasXOR_v2, &tHlts.b_HLT_HIZeroBiasXOR_v2);
   if (t->GetBranch("HLT_HIZeroBiasXOR_v2_Prescl")) t->SetBranchAddress("HLT_HIZeroBiasXOR_v2_Prescl", &tHlts.HLT_HIZeroBiasXOR_v2_Prescl, &tHlts.b_HLT_HIZeroBiasXOR_v2_Prescl);
   if (t->GetBranch("HLT_HIZeroBias_v2")) t->SetBranchAddress("HLT_HIZeroBias_v2", &tHlts.HLT_HIZeroBias_v2, &tHlts.b_HLT_HIZeroBias_v2);
   if (t->GetBranch("HLT_HIZeroBias_v2_Prescl")) t->SetBranchAddress("HLT_HIZeroBias_v2_Prescl", &tHlts.HLT_HIZeroBias_v2_Prescl, &tHlts.b_HLT_HIZeroBias_v2_Prescl);
   if (t->GetBranch("HLT_HIL1DoubleMu0_HighQ_v3")) t->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v3", &tHlts.HLT_HIL1DoubleMu0_HighQ_v3, &tHlts.b_HLT_HIL1DoubleMu0_HighQ_v3);
   if (t->GetBranch("HLT_HIL1DoubleMu0_HighQ_v3_Prescl")) t->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v3_Prescl", &tHlts.HLT_HIL1DoubleMu0_HighQ_v3_Prescl, &tHlts.b_HLT_HIL1DoubleMu0_HighQ_v3_Prescl);
   if (t->GetBranch("HLT_HIL1DoubleMuOpen_v3")) t->SetBranchAddress("HLT_HIL1DoubleMuOpen_v3", &tHlts.HLT_HIL1DoubleMuOpen_v3, &tHlts.b_HLT_HIL1DoubleMuOpen_v3);
   if (t->GetBranch("HLT_HIL1DoubleMuOpen_v3_Prescl")) t->SetBranchAddress("HLT_HIL1DoubleMuOpen_v3_Prescl", &tHlts.HLT_HIL1DoubleMuOpen_v3_Prescl, &tHlts.b_HLT_HIL1DoubleMuOpen_v3_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3")) t->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3", &tHlts.HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3, &tHlts.b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3);
   if (t->GetBranch("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl", &tHlts.HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_NHitQ_v3")) t->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v3", &tHlts.HLT_HIL2DoubleMu0_NHitQ_v3, &tHlts.b_HLT_HIL2DoubleMu0_NHitQ_v3);
   if (t->GetBranch("HLT_HIL2DoubleMu0_NHitQ_v3_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v3_Prescl", &tHlts.HLT_HIL2DoubleMu0_NHitQ_v3_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_NHitQ_v3_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_v3")) t->SetBranchAddress("HLT_HIL2DoubleMu0_v3", &tHlts.HLT_HIL2DoubleMu0_v3, &tHlts.b_HLT_HIL2DoubleMu0_v3);
   if (t->GetBranch("HLT_HIL2DoubleMu0_v3_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_v3_Prescl", &tHlts.HLT_HIL2DoubleMu0_v3_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_v3_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu3_v3")) t->SetBranchAddress("HLT_HIL2DoubleMu3_v3", &tHlts.HLT_HIL2DoubleMu3_v3, &tHlts.b_HLT_HIL2DoubleMu3_v3);
   if (t->GetBranch("HLT_HIL2DoubleMu3_v3_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu3_v3_Prescl", &tHlts.HLT_HIL2DoubleMu3_v3_Prescl, &tHlts.b_HLT_HIL2DoubleMu3_v3_Prescl);
   if (t->GetBranch("HLT_HIL2Mu15_v3")) t->SetBranchAddress("HLT_HIL2Mu15_v3", &tHlts.HLT_HIL2Mu15_v3, &tHlts.b_HLT_HIL2Mu15_v3);
   if (t->GetBranch("HLT_HIL2Mu15_v3_Prescl")) t->SetBranchAddress("HLT_HIL2Mu15_v3_Prescl", &tHlts.HLT_HIL2Mu15_v3_Prescl, &tHlts.b_HLT_HIL2Mu15_v3_Prescl);
   if (t->GetBranch("HLT_HIL2Mu3_NHitQ_v3")) t->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v3", &tHlts.HLT_HIL2Mu3_NHitQ_v3, &tHlts.b_HLT_HIL2Mu3_NHitQ_v3);
   if (t->GetBranch("HLT_HIL2Mu3_NHitQ_v3_Prescl")) t->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v3_Prescl", &tHlts.HLT_HIL2Mu3_NHitQ_v3_Prescl, &tHlts.b_HLT_HIL2Mu3_NHitQ_v3_Prescl);
   if (t->GetBranch("HLT_HIL2Mu3_v3")) t->SetBranchAddress("HLT_HIL2Mu3_v3", &tHlts.HLT_HIL2Mu3_v3, &tHlts.b_HLT_HIL2Mu3_v3);
   if (t->GetBranch("HLT_HIL2Mu3_v3_Prescl")) t->SetBranchAddress("HLT_HIL2Mu3_v3_Prescl", &tHlts.HLT_HIL2Mu3_v3_Prescl, &tHlts.b_HLT_HIL2Mu3_v3_Prescl);
   if (t->GetBranch("HLT_HIL2Mu7_v3")) t->SetBranchAddress("HLT_HIL2Mu7_v3", &tHlts.HLT_HIL2Mu7_v3, &tHlts.b_HLT_HIL2Mu7_v3);
   if (t->GetBranch("HLT_HIL2Mu7_v3_Prescl")) t->SetBranchAddress("HLT_HIL2Mu7_v3_Prescl", &tHlts.HLT_HIL2Mu7_v3_Prescl, &tHlts.b_HLT_HIL2Mu7_v3_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_v3")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v3", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_v3, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v3);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_SS_v3")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v3", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_SS_v3, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v3);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_v3")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v3", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_v3, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_v3);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_v3")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_v3", &tHlts.HLT_HIL3DoubleMuOpen_v3, &tHlts.b_HLT_HIL3DoubleMuOpen_v3);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_v3_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_v3_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_v3_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_v3_Prescl);
   if (t->GetBranch("HLT_HIL3Mu3_v3")) t->SetBranchAddress("HLT_HIL3Mu3_v3", &tHlts.HLT_HIL3Mu3_v3, &tHlts.b_HLT_HIL3Mu3_v3);
   if (t->GetBranch("HLT_HIL3Mu3_v3_Prescl")) t->SetBranchAddress("HLT_HIL3Mu3_v3_Prescl", &tHlts.HLT_HIL3Mu3_v3_Prescl, &tHlts.b_HLT_HIL3Mu3_v3_Prescl);
   if (t->GetBranch("HLT_HIDiJet55_v3")) t->SetBranchAddress("HLT_HIDiJet55_v3", &tHlts.HLT_HIDiJet55_v3, &tHlts.b_HLT_HIDiJet55_v3);
   if (t->GetBranch("HLT_HIDiJet55_v3_Prescl")) t->SetBranchAddress("HLT_HIDiJet55_v3_Prescl", &tHlts.HLT_HIDiJet55_v3_Prescl, &tHlts.b_HLT_HIDiJet55_v3_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton10_v3")) t->SetBranchAddress("HLT_HIDoublePhoton10_v3", &tHlts.HLT_HIDoublePhoton10_v3, &tHlts.b_HLT_HIDoublePhoton10_v3);
   if (t->GetBranch("HLT_HIDoublePhoton10_v3_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton10_v3_Prescl", &tHlts.HLT_HIDoublePhoton10_v3_Prescl, &tHlts.b_HLT_HIDoublePhoton10_v3_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton15_v3")) t->SetBranchAddress("HLT_HIDoublePhoton15_v3", &tHlts.HLT_HIDoublePhoton15_v3, &tHlts.b_HLT_HIDoublePhoton15_v3);
   if (t->GetBranch("HLT_HIDoublePhoton15_v3_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton15_v3_Prescl", &tHlts.HLT_HIDoublePhoton15_v3_Prescl, &tHlts.b_HLT_HIDoublePhoton15_v3_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton20_v3")) t->SetBranchAddress("HLT_HIDoublePhoton20_v3", &tHlts.HLT_HIDoublePhoton20_v3, &tHlts.b_HLT_HIDoublePhoton20_v3);
   if (t->GetBranch("HLT_HIDoublePhoton20_v3_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton20_v3_Prescl", &tHlts.HLT_HIDoublePhoton20_v3_Prescl, &tHlts.b_HLT_HIDoublePhoton20_v3_Prescl);
   if (t->GetBranch("HLT_HIFullTrack12_L1Central_v3")) t->SetBranchAddress("HLT_HIFullTrack12_L1Central_v3", &tHlts.HLT_HIFullTrack12_L1Central_v3, &tHlts.b_HLT_HIFullTrack12_L1Central_v3);
   if (t->GetBranch("HLT_HIFullTrack12_L1Central_v3_Prescl")) t->SetBranchAddress("HLT_HIFullTrack12_L1Central_v3_Prescl", &tHlts.HLT_HIFullTrack12_L1Central_v3_Prescl, &tHlts.b_HLT_HIFullTrack12_L1Central_v3_Prescl);
   if (t->GetBranch("HLT_HIFullTrack12_L1Peripheral_v3")) t->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v3", &tHlts.HLT_HIFullTrack12_L1Peripheral_v3, &tHlts.b_HLT_HIFullTrack12_L1Peripheral_v3);
   if (t->GetBranch("HLT_HIFullTrack12_L1Peripheral_v3_Prescl")) t->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v3_Prescl", &tHlts.HLT_HIFullTrack12_L1Peripheral_v3_Prescl, &tHlts.b_HLT_HIFullTrack12_L1Peripheral_v3_Prescl);
   if (t->GetBranch("HLT_HIFullTrack14_L1Central_v3")) t->SetBranchAddress("HLT_HIFullTrack14_L1Central_v3", &tHlts.HLT_HIFullTrack14_L1Central_v3, &tHlts.b_HLT_HIFullTrack14_L1Central_v3);
   if (t->GetBranch("HLT_HIFullTrack14_L1Central_v3_Prescl")) t->SetBranchAddress("HLT_HIFullTrack14_L1Central_v3_Prescl", &tHlts.HLT_HIFullTrack14_L1Central_v3_Prescl, &tHlts.b_HLT_HIFullTrack14_L1Central_v3_Prescl);
   if (t->GetBranch("HLT_HIFullTrack14_L1Peripheral_v3")) t->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v3", &tHlts.HLT_HIFullTrack14_L1Peripheral_v3, &tHlts.b_HLT_HIFullTrack14_L1Peripheral_v3);
   if (t->GetBranch("HLT_HIFullTrack14_L1Peripheral_v3_Prescl")) t->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v3_Prescl", &tHlts.HLT_HIFullTrack14_L1Peripheral_v3_Prescl, &tHlts.b_HLT_HIFullTrack14_L1Peripheral_v3_Prescl);
   if (t->GetBranch("HLT_HIFullTrack20_L1Central_v3")) t->SetBranchAddress("HLT_HIFullTrack20_L1Central_v3", &tHlts.HLT_HIFullTrack20_L1Central_v3, &tHlts.b_HLT_HIFullTrack20_L1Central_v3);
   if (t->GetBranch("HLT_HIFullTrack20_L1Central_v3_Prescl")) t->SetBranchAddress("HLT_HIFullTrack20_L1Central_v3_Prescl", &tHlts.HLT_HIFullTrack20_L1Central_v3_Prescl, &tHlts.b_HLT_HIFullTrack20_L1Central_v3_Prescl);
   if (t->GetBranch("HLT_HIFullTrack20_L1Peripheral_v3")) t->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v3", &tHlts.HLT_HIFullTrack20_L1Peripheral_v3, &tHlts.b_HLT_HIFullTrack20_L1Peripheral_v3);
   if (t->GetBranch("HLT_HIFullTrack20_L1Peripheral_v3_Prescl")) t->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v3_Prescl", &tHlts.HLT_HIFullTrack20_L1Peripheral_v3_Prescl, &tHlts.b_HLT_HIFullTrack20_L1Peripheral_v3_Prescl);
   if (t->GetBranch("HLT_HIFullTrack25_L1Central_v3")) t->SetBranchAddress("HLT_HIFullTrack25_L1Central_v3", &tHlts.HLT_HIFullTrack25_L1Central_v3, &tHlts.b_HLT_HIFullTrack25_L1Central_v3);
   if (t->GetBranch("HLT_HIFullTrack25_L1Central_v3_Prescl")) t->SetBranchAddress("HLT_HIFullTrack25_L1Central_v3_Prescl", &tHlts.HLT_HIFullTrack25_L1Central_v3_Prescl, &tHlts.b_HLT_HIFullTrack25_L1Central_v3_Prescl);
   if (t->GetBranch("HLT_HIFullTrack25_L1Peripheral_v3")) t->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v3", &tHlts.HLT_HIFullTrack25_L1Peripheral_v3, &tHlts.b_HLT_HIFullTrack25_L1Peripheral_v3);
   if (t->GetBranch("HLT_HIFullTrack25_L1Peripheral_v3_Prescl")) t->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v3_Prescl", &tHlts.HLT_HIFullTrack25_L1Peripheral_v3_Prescl, &tHlts.b_HLT_HIFullTrack25_L1Peripheral_v3_Prescl);
   if (t->GetBranch("HLT_HIJet55_v3")) t->SetBranchAddress("HLT_HIJet55_v3", &tHlts.HLT_HIJet55_v3, &tHlts.b_HLT_HIJet55_v3);
   if (t->GetBranch("HLT_HIJet55_v3_Prescl")) t->SetBranchAddress("HLT_HIJet55_v3_Prescl", &tHlts.HLT_HIJet55_v3_Prescl, &tHlts.b_HLT_HIJet55_v3_Prescl);
   if (t->GetBranch("HLT_HIJet65_Jet55_v3")) t->SetBranchAddress("HLT_HIJet65_Jet55_v3", &tHlts.HLT_HIJet65_Jet55_v3, &tHlts.b_HLT_HIJet65_Jet55_v3);
   if (t->GetBranch("HLT_HIJet65_Jet55_v3_Prescl")) t->SetBranchAddress("HLT_HIJet65_Jet55_v3_Prescl", &tHlts.HLT_HIJet65_Jet55_v3_Prescl, &tHlts.b_HLT_HIJet65_Jet55_v3_Prescl);
   if (t->GetBranch("HLT_HIJet65_v3")) t->SetBranchAddress("HLT_HIJet65_v3", &tHlts.HLT_HIJet65_v3, &tHlts.b_HLT_HIJet65_v3);
   if (t->GetBranch("HLT_HIJet65_v3_Prescl")) t->SetBranchAddress("HLT_HIJet65_v3_Prescl", &tHlts.HLT_HIJet65_v3_Prescl, &tHlts.b_HLT_HIJet65_v3_Prescl);
   if (t->GetBranch("HLT_HIJet80_v3")) t->SetBranchAddress("HLT_HIJet80_v3", &tHlts.HLT_HIJet80_v3, &tHlts.b_HLT_HIJet80_v3);
   if (t->GetBranch("HLT_HIJet80_v3_Prescl")) t->SetBranchAddress("HLT_HIJet80_v3_Prescl", &tHlts.HLT_HIJet80_v3_Prescl, &tHlts.b_HLT_HIJet80_v3_Prescl);
   if (t->GetBranch("HLT_HIJet95_v3")) t->SetBranchAddress("HLT_HIJet95_v3", &tHlts.HLT_HIJet95_v3, &tHlts.b_HLT_HIJet95_v3);
   if (t->GetBranch("HLT_HIJet95_v3_Prescl")) t->SetBranchAddress("HLT_HIJet95_v3_Prescl", &tHlts.HLT_HIJet95_v3_Prescl, &tHlts.b_HLT_HIJet95_v3_Prescl);
   if (t->GetBranch("HLT_HIJetE30_NoBPTX_v3")) t->SetBranchAddress("HLT_HIJetE30_NoBPTX_v3", &tHlts.HLT_HIJetE30_NoBPTX_v3, &tHlts.b_HLT_HIJetE30_NoBPTX_v3);
   if (t->GetBranch("HLT_HIJetE30_NoBPTX_v3_Prescl")) t->SetBranchAddress("HLT_HIJetE30_NoBPTX_v3_Prescl", &tHlts.HLT_HIJetE30_NoBPTX_v3_Prescl, &tHlts.b_HLT_HIJetE30_NoBPTX_v3_Prescl);
   if (t->GetBranch("HLT_HIJetE50_NoBPTX3BX_NoHalo_v3")) t->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v3", &tHlts.HLT_HIJetE50_NoBPTX3BX_NoHalo_v3, &tHlts.b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v3);
   if (t->GetBranch("HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl")) t->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl", &tHlts.HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl, &tHlts.b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl);
   if (t->GetBranch("HLT_HIMET120_v3")) t->SetBranchAddress("HLT_HIMET120_v3", &tHlts.HLT_HIMET120_v3, &tHlts.b_HLT_HIMET120_v3);
   if (t->GetBranch("HLT_HIMET120_v3_Prescl")) t->SetBranchAddress("HLT_HIMET120_v3_Prescl", &tHlts.HLT_HIMET120_v3_Prescl, &tHlts.b_HLT_HIMET120_v3_Prescl);
   if (t->GetBranch("HLT_HIMET200_v3")) t->SetBranchAddress("HLT_HIMET200_v3", &tHlts.HLT_HIMET200_v3, &tHlts.b_HLT_HIMET200_v3);
   if (t->GetBranch("HLT_HIMET200_v3_Prescl")) t->SetBranchAddress("HLT_HIMET200_v3_Prescl", &tHlts.HLT_HIMET200_v3_Prescl, &tHlts.b_HLT_HIMET200_v3_Prescl);
   if (t->GetBranch("HLT_HIMET220_v3")) t->SetBranchAddress("HLT_HIMET220_v3", &tHlts.HLT_HIMET220_v3, &tHlts.b_HLT_HIMET220_v3);
   if (t->GetBranch("HLT_HIMET220_v3_Prescl")) t->SetBranchAddress("HLT_HIMET220_v3_Prescl", &tHlts.HLT_HIMET220_v3_Prescl, &tHlts.b_HLT_HIMET220_v3_Prescl);
   if (t->GetBranch("HLT_HIPhoton10_Photon15_v3")) t->SetBranchAddress("HLT_HIPhoton10_Photon15_v3", &tHlts.HLT_HIPhoton10_Photon15_v3, &tHlts.b_HLT_HIPhoton10_Photon15_v3);
   if (t->GetBranch("HLT_HIPhoton10_Photon15_v3_Prescl")) t->SetBranchAddress("HLT_HIPhoton10_Photon15_v3_Prescl", &tHlts.HLT_HIPhoton10_Photon15_v3_Prescl, &tHlts.b_HLT_HIPhoton10_Photon15_v3_Prescl);
   if (t->GetBranch("HLT_HIPhoton15_Photon20_v3")) t->SetBranchAddress("HLT_HIPhoton15_Photon20_v3", &tHlts.HLT_HIPhoton15_Photon20_v3, &tHlts.b_HLT_HIPhoton15_Photon20_v3);
   if (t->GetBranch("HLT_HIPhoton15_Photon20_v3_Prescl")) t->SetBranchAddress("HLT_HIPhoton15_Photon20_v3_Prescl", &tHlts.HLT_HIPhoton15_Photon20_v3_Prescl, &tHlts.b_HLT_HIPhoton15_Photon20_v3_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton15_v3")) t->SetBranchAddress("HLT_HISinglePhoton15_v3", &tHlts.HLT_HISinglePhoton15_v3, &tHlts.b_HLT_HISinglePhoton15_v3);
   if (t->GetBranch("HLT_HISinglePhoton15_v3_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton15_v3_Prescl", &tHlts.HLT_HISinglePhoton15_v3_Prescl, &tHlts.b_HLT_HISinglePhoton15_v3_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton20_v3")) t->SetBranchAddress("HLT_HISinglePhoton20_v3", &tHlts.HLT_HISinglePhoton20_v3, &tHlts.b_HLT_HISinglePhoton20_v3);
   if (t->GetBranch("HLT_HISinglePhoton20_v3_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton20_v3_Prescl", &tHlts.HLT_HISinglePhoton20_v3_Prescl, &tHlts.b_HLT_HISinglePhoton20_v3_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton30_v3")) t->SetBranchAddress("HLT_HISinglePhoton30_v3", &tHlts.HLT_HISinglePhoton30_v3, &tHlts.b_HLT_HISinglePhoton30_v3);
   if (t->GetBranch("HLT_HISinglePhoton30_v3_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton30_v3_Prescl", &tHlts.HLT_HISinglePhoton30_v3_Prescl, &tHlts.b_HLT_HISinglePhoton30_v3_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton40_v3")) t->SetBranchAddress("HLT_HISinglePhoton40_v3", &tHlts.HLT_HISinglePhoton40_v3, &tHlts.b_HLT_HISinglePhoton40_v3);
   if (t->GetBranch("HLT_HISinglePhoton40_v3_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton40_v3_Prescl", &tHlts.HLT_HISinglePhoton40_v3_Prescl, &tHlts.b_HLT_HISinglePhoton40_v3_Prescl);
   if (t->GetBranch("HLT_HIBptxXOR_v3")) t->SetBranchAddress("HLT_HIBptxXOR_v3", &tHlts.HLT_HIBptxXOR_v3, &tHlts.b_HLT_HIBptxXOR_v3);
   if (t->GetBranch("HLT_HIBptxXOR_v3_Prescl")) t->SetBranchAddress("HLT_HIBptxXOR_v3_Prescl", &tHlts.HLT_HIBptxXOR_v3_Prescl, &tHlts.b_HLT_HIBptxXOR_v3_Prescl);
   if (t->GetBranch("HLT_HICentral10_v3")) t->SetBranchAddress("HLT_HICentral10_v3", &tHlts.HLT_HICentral10_v3, &tHlts.b_HLT_HICentral10_v3);
   if (t->GetBranch("HLT_HICentral10_v3_Prescl")) t->SetBranchAddress("HLT_HICentral10_v3_Prescl", &tHlts.HLT_HICentral10_v3_Prescl, &tHlts.b_HLT_HICentral10_v3_Prescl);
   if (t->GetBranch("HLT_HICentralityVeto_v3")) t->SetBranchAddress("HLT_HICentralityVeto_v3", &tHlts.HLT_HICentralityVeto_v3, &tHlts.b_HLT_HICentralityVeto_v3);
   if (t->GetBranch("HLT_HICentralityVeto_v3_Prescl")) t->SetBranchAddress("HLT_HICentralityVeto_v3_Prescl", &tHlts.HLT_HICentralityVeto_v3_Prescl, &tHlts.b_HLT_HICentralityVeto_v3_Prescl);
   if (t->GetBranch("HLT_HIL1Algo_BptxXOR_BSC_OR_v3")) t->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v3", &tHlts.HLT_HIL1Algo_BptxXOR_BSC_OR_v3, &tHlts.b_HLT_HIL1Algo_BptxXOR_BSC_OR_v3);
   if (t->GetBranch("HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl")) t->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl", &tHlts.HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl, &tHlts.b_HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasBSC_OR_v3")) t->SetBranchAddress("HLT_HIMinBiasBSC_OR_v3", &tHlts.HLT_HIMinBiasBSC_OR_v3, &tHlts.b_HLT_HIMinBiasBSC_OR_v3);
   if (t->GetBranch("HLT_HIMinBiasBSC_OR_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasBSC_OR_v3_Prescl", &tHlts.HLT_HIMinBiasBSC_OR_v3_Prescl, &tHlts.b_HLT_HIMinBiasBSC_OR_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasBSC_v3")) t->SetBranchAddress("HLT_HIMinBiasBSC_v3", &tHlts.HLT_HIMinBiasBSC_v3, &tHlts.b_HLT_HIMinBiasBSC_v3);
   if (t->GetBranch("HLT_HIMinBiasBSC_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasBSC_v3_Prescl", &tHlts.HLT_HIMinBiasBSC_v3_Prescl, &tHlts.b_HLT_HIMinBiasBSC_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHF_v3")) t->SetBranchAddress("HLT_HIMinBiasHF_v3", &tHlts.HLT_HIMinBiasHF_v3, &tHlts.b_HLT_HIMinBiasHF_v3);
   if (t->GetBranch("HLT_HIMinBiasHF_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHF_v3_Prescl", &tHlts.HLT_HIMinBiasHF_v3_Prescl, &tHlts.b_HLT_HIMinBiasHF_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHfOrBSC_v3")) t->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v3", &tHlts.HLT_HIMinBiasHfOrBSC_v3, &tHlts.b_HLT_HIMinBiasHfOrBSC_v3);
   if (t->GetBranch("HLT_HIMinBiasHfOrBSC_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v3_Prescl", &tHlts.HLT_HIMinBiasHfOrBSC_v3_Prescl, &tHlts.b_HLT_HIMinBiasHfOrBSC_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHf_OR_v3")) t->SetBranchAddress("HLT_HIMinBiasHf_OR_v3", &tHlts.HLT_HIMinBiasHf_OR_v3, &tHlts.b_HLT_HIMinBiasHf_OR_v3);
   if (t->GetBranch("HLT_HIMinBiasHf_OR_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHf_OR_v3_Prescl", &tHlts.HLT_HIMinBiasHf_OR_v3_Prescl, &tHlts.b_HLT_HIMinBiasHf_OR_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasPixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v3", &tHlts.HLT_HIMinBiasPixel_SingleTrack_v3, &tHlts.b_HLT_HIMinBiasPixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIMinBiasPixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIMinBiasPixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIMinBiasPixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDCPixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v3", &tHlts.HLT_HIMinBiasZDCPixel_SingleTrack_v3, &tHlts.b_HLT_HIMinBiasZDCPixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3", &tHlts.HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3, &tHlts.b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl", &tHlts.HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl, &tHlts.b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_v3")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v3", &tHlts.HLT_HIMinBiasZDC_Calo_v3, &tHlts.b_HLT_HIMinBiasZDC_Calo_v3);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v3_Prescl", &tHlts.HLT_HIMinBiasZDC_Calo_v3_Prescl, &tHlts.b_HLT_HIMinBiasZDC_Calo_v3_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3", &tHlts.HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3, &tHlts.b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIPhysics_v3")) t->SetBranchAddress("HLT_HIPhysics_v3", &tHlts.HLT_HIPhysics_v3, &tHlts.b_HLT_HIPhysics_v3);
   if (t->GetBranch("HLT_HIPhysics_v3_Prescl")) t->SetBranchAddress("HLT_HIPhysics_v3_Prescl", &tHlts.HLT_HIPhysics_v3_Prescl, &tHlts.b_HLT_HIPhysics_v3_Prescl);
   if (t->GetBranch("HLT_HIRandom_v3")) t->SetBranchAddress("HLT_HIRandom_v3", &tHlts.HLT_HIRandom_v3, &tHlts.b_HLT_HIRandom_v3);
   if (t->GetBranch("HLT_HIRandom_v3_Prescl")) t->SetBranchAddress("HLT_HIRandom_v3_Prescl", &tHlts.HLT_HIRandom_v3_Prescl, &tHlts.b_HLT_HIRandom_v3_Prescl);
   if (t->GetBranch("HLT_HIUCC010_v3")) t->SetBranchAddress("HLT_HIUCC010_v3", &tHlts.HLT_HIUCC010_v3, &tHlts.b_HLT_HIUCC010_v3);
   if (t->GetBranch("HLT_HIUCC010_v3_Prescl")) t->SetBranchAddress("HLT_HIUCC010_v3_Prescl", &tHlts.HLT_HIUCC010_v3_Prescl, &tHlts.b_HLT_HIUCC010_v3_Prescl);
   if (t->GetBranch("HLT_HIUCC015_v3")) t->SetBranchAddress("HLT_HIUCC015_v3", &tHlts.HLT_HIUCC015_v3, &tHlts.b_HLT_HIUCC015_v3);
   if (t->GetBranch("HLT_HIUCC015_v3_Prescl")) t->SetBranchAddress("HLT_HIUCC015_v3_Prescl", &tHlts.HLT_HIUCC015_v3_Prescl, &tHlts.b_HLT_HIUCC015_v3_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuEG2Pixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v3", &tHlts.HLT_HIUPCNeuEG2Pixel_SingleTrack_v3, &tHlts.b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuEG5Pixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v3", &tHlts.HLT_HIUPCNeuEG5Pixel_SingleTrack_v3, &tHlts.b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3", &tHlts.HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3, &tHlts.b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3", &tHlts.HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3, &tHlts.b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3", &tHlts.HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3, &tHlts.b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuMuPixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v3", &tHlts.HLT_HIUPCNeuMuPixel_SingleTrack_v3, &tHlts.b_HLT_HIUPCNeuMuPixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIZeroBiasPixel_SingleTrack_v3")) t->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v3", &tHlts.HLT_HIZeroBiasPixel_SingleTrack_v3, &tHlts.b_HLT_HIZeroBiasPixel_SingleTrack_v3);
   if (t->GetBranch("HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl")) t->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl", &tHlts.HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl, &tHlts.b_HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl);
   if (t->GetBranch("HLT_HIZeroBiasXOR_v3")) t->SetBranchAddress("HLT_HIZeroBiasXOR_v3", &tHlts.HLT_HIZeroBiasXOR_v3, &tHlts.b_HLT_HIZeroBiasXOR_v3);
   if (t->GetBranch("HLT_HIZeroBiasXOR_v3_Prescl")) t->SetBranchAddress("HLT_HIZeroBiasXOR_v3_Prescl", &tHlts.HLT_HIZeroBiasXOR_v3_Prescl, &tHlts.b_HLT_HIZeroBiasXOR_v3_Prescl);
   if (t->GetBranch("HLT_HIZeroBias_v3")) t->SetBranchAddress("HLT_HIZeroBias_v3", &tHlts.HLT_HIZeroBias_v3, &tHlts.b_HLT_HIZeroBias_v3);
   if (t->GetBranch("HLT_HIZeroBias_v3_Prescl")) t->SetBranchAddress("HLT_HIZeroBias_v3_Prescl", &tHlts.HLT_HIZeroBias_v3_Prescl, &tHlts.b_HLT_HIZeroBias_v3_Prescl);
   if (t->GetBranch("HLT_HIL1DoubleMu0_HighQ_v4")) t->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v4", &tHlts.HLT_HIL1DoubleMu0_HighQ_v4, &tHlts.b_HLT_HIL1DoubleMu0_HighQ_v4);
   if (t->GetBranch("HLT_HIL1DoubleMu0_HighQ_v4_Prescl")) t->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v4_Prescl", &tHlts.HLT_HIL1DoubleMu0_HighQ_v4_Prescl, &tHlts.b_HLT_HIL1DoubleMu0_HighQ_v4_Prescl);
   if (t->GetBranch("HLT_HIL1DoubleMuOpen_v4")) t->SetBranchAddress("HLT_HIL1DoubleMuOpen_v4", &tHlts.HLT_HIL1DoubleMuOpen_v4, &tHlts.b_HLT_HIL1DoubleMuOpen_v4);
   if (t->GetBranch("HLT_HIL1DoubleMuOpen_v4_Prescl")) t->SetBranchAddress("HLT_HIL1DoubleMuOpen_v4_Prescl", &tHlts.HLT_HIL1DoubleMuOpen_v4_Prescl, &tHlts.b_HLT_HIL1DoubleMuOpen_v4_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4")) t->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4", &tHlts.HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4, &tHlts.b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4);
   if (t->GetBranch("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl", &tHlts.HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_NHitQ_v4")) t->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v4", &tHlts.HLT_HIL2DoubleMu0_NHitQ_v4, &tHlts.b_HLT_HIL2DoubleMu0_NHitQ_v4);
   if (t->GetBranch("HLT_HIL2DoubleMu0_NHitQ_v4_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v4_Prescl", &tHlts.HLT_HIL2DoubleMu0_NHitQ_v4_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_NHitQ_v4_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu0_v4")) t->SetBranchAddress("HLT_HIL2DoubleMu0_v4", &tHlts.HLT_HIL2DoubleMu0_v4, &tHlts.b_HLT_HIL2DoubleMu0_v4);
   if (t->GetBranch("HLT_HIL2DoubleMu0_v4_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu0_v4_Prescl", &tHlts.HLT_HIL2DoubleMu0_v4_Prescl, &tHlts.b_HLT_HIL2DoubleMu0_v4_Prescl);
   if (t->GetBranch("HLT_HIL2DoubleMu3_v4")) t->SetBranchAddress("HLT_HIL2DoubleMu3_v4", &tHlts.HLT_HIL2DoubleMu3_v4, &tHlts.b_HLT_HIL2DoubleMu3_v4);
   if (t->GetBranch("HLT_HIL2DoubleMu3_v4_Prescl")) t->SetBranchAddress("HLT_HIL2DoubleMu3_v4_Prescl", &tHlts.HLT_HIL2DoubleMu3_v4_Prescl, &tHlts.b_HLT_HIL2DoubleMu3_v4_Prescl);
   if (t->GetBranch("HLT_HIL2Mu15_v4")) t->SetBranchAddress("HLT_HIL2Mu15_v4", &tHlts.HLT_HIL2Mu15_v4, &tHlts.b_HLT_HIL2Mu15_v4);
   if (t->GetBranch("HLT_HIL2Mu15_v4_Prescl")) t->SetBranchAddress("HLT_HIL2Mu15_v4_Prescl", &tHlts.HLT_HIL2Mu15_v4_Prescl, &tHlts.b_HLT_HIL2Mu15_v4_Prescl);
   if (t->GetBranch("HLT_HIL2Mu3_NHitQ_v4")) t->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v4", &tHlts.HLT_HIL2Mu3_NHitQ_v4, &tHlts.b_HLT_HIL2Mu3_NHitQ_v4);
   if (t->GetBranch("HLT_HIL2Mu3_NHitQ_v4_Prescl")) t->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v4_Prescl", &tHlts.HLT_HIL2Mu3_NHitQ_v4_Prescl, &tHlts.b_HLT_HIL2Mu3_NHitQ_v4_Prescl);
   if (t->GetBranch("HLT_HIL2Mu3_v4")) t->SetBranchAddress("HLT_HIL2Mu3_v4", &tHlts.HLT_HIL2Mu3_v4, &tHlts.b_HLT_HIL2Mu3_v4);
   if (t->GetBranch("HLT_HIL2Mu3_v4_Prescl")) t->SetBranchAddress("HLT_HIL2Mu3_v4_Prescl", &tHlts.HLT_HIL2Mu3_v4_Prescl, &tHlts.b_HLT_HIL2Mu3_v4_Prescl);
   if (t->GetBranch("HLT_HIL2Mu7_v4")) t->SetBranchAddress("HLT_HIL2Mu7_v4", &tHlts.HLT_HIL2Mu7_v4, &tHlts.b_HLT_HIL2Mu7_v4);
   if (t->GetBranch("HLT_HIL2Mu7_v4_Prescl")) t->SetBranchAddress("HLT_HIL2Mu7_v4_Prescl", &tHlts.HLT_HIL2Mu7_v4_Prescl, &tHlts.b_HLT_HIL2Mu7_v4_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_v4")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v4", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_v4, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v4);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_SS_v4")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v4", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_SS_v4, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v4);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_v4")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v4", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_v4, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_v4);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_v4")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_v4", &tHlts.HLT_HIL3DoubleMuOpen_v4, &tHlts.b_HLT_HIL3DoubleMuOpen_v4);
   if (t->GetBranch("HLT_HIL3DoubleMuOpen_v4_Prescl")) t->SetBranchAddress("HLT_HIL3DoubleMuOpen_v4_Prescl", &tHlts.HLT_HIL3DoubleMuOpen_v4_Prescl, &tHlts.b_HLT_HIL3DoubleMuOpen_v4_Prescl);
   if (t->GetBranch("HLT_HIL3Mu3_v4")) t->SetBranchAddress("HLT_HIL3Mu3_v4", &tHlts.HLT_HIL3Mu3_v4, &tHlts.b_HLT_HIL3Mu3_v4);
   if (t->GetBranch("HLT_HIL3Mu3_v4_Prescl")) t->SetBranchAddress("HLT_HIL3Mu3_v4_Prescl", &tHlts.HLT_HIL3Mu3_v4_Prescl, &tHlts.b_HLT_HIL3Mu3_v4_Prescl);
   if (t->GetBranch("HLT_HIDiJet55_v4")) t->SetBranchAddress("HLT_HIDiJet55_v4", &tHlts.HLT_HIDiJet55_v4, &tHlts.b_HLT_HIDiJet55_v4);
   if (t->GetBranch("HLT_HIDiJet55_v4_Prescl")) t->SetBranchAddress("HLT_HIDiJet55_v4_Prescl", &tHlts.HLT_HIDiJet55_v4_Prescl, &tHlts.b_HLT_HIDiJet55_v4_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton10_v4")) t->SetBranchAddress("HLT_HIDoublePhoton10_v4", &tHlts.HLT_HIDoublePhoton10_v4, &tHlts.b_HLT_HIDoublePhoton10_v4);
   if (t->GetBranch("HLT_HIDoublePhoton10_v4_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton10_v4_Prescl", &tHlts.HLT_HIDoublePhoton10_v4_Prescl, &tHlts.b_HLT_HIDoublePhoton10_v4_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton15_v4")) t->SetBranchAddress("HLT_HIDoublePhoton15_v4", &tHlts.HLT_HIDoublePhoton15_v4, &tHlts.b_HLT_HIDoublePhoton15_v4);
   if (t->GetBranch("HLT_HIDoublePhoton15_v4_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton15_v4_Prescl", &tHlts.HLT_HIDoublePhoton15_v4_Prescl, &tHlts.b_HLT_HIDoublePhoton15_v4_Prescl);
   if (t->GetBranch("HLT_HIDoublePhoton20_v4")) t->SetBranchAddress("HLT_HIDoublePhoton20_v4", &tHlts.HLT_HIDoublePhoton20_v4, &tHlts.b_HLT_HIDoublePhoton20_v4);
   if (t->GetBranch("HLT_HIDoublePhoton20_v4_Prescl")) t->SetBranchAddress("HLT_HIDoublePhoton20_v4_Prescl", &tHlts.HLT_HIDoublePhoton20_v4_Prescl, &tHlts.b_HLT_HIDoublePhoton20_v4_Prescl);
   if (t->GetBranch("HLT_HIFullTrack12_L1Central_v4")) t->SetBranchAddress("HLT_HIFullTrack12_L1Central_v4", &tHlts.HLT_HIFullTrack12_L1Central_v4, &tHlts.b_HLT_HIFullTrack12_L1Central_v4);
   if (t->GetBranch("HLT_HIFullTrack12_L1Central_v4_Prescl")) t->SetBranchAddress("HLT_HIFullTrack12_L1Central_v4_Prescl", &tHlts.HLT_HIFullTrack12_L1Central_v4_Prescl, &tHlts.b_HLT_HIFullTrack12_L1Central_v4_Prescl);
   if (t->GetBranch("HLT_HIFullTrack12_L1Peripheral_v4")) t->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v4", &tHlts.HLT_HIFullTrack12_L1Peripheral_v4, &tHlts.b_HLT_HIFullTrack12_L1Peripheral_v4);
   if (t->GetBranch("HLT_HIFullTrack12_L1Peripheral_v4_Prescl")) t->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v4_Prescl", &tHlts.HLT_HIFullTrack12_L1Peripheral_v4_Prescl, &tHlts.b_HLT_HIFullTrack12_L1Peripheral_v4_Prescl);
   if (t->GetBranch("HLT_HIFullTrack14_L1Central_v4")) t->SetBranchAddress("HLT_HIFullTrack14_L1Central_v4", &tHlts.HLT_HIFullTrack14_L1Central_v4, &tHlts.b_HLT_HIFullTrack14_L1Central_v4);
   if (t->GetBranch("HLT_HIFullTrack14_L1Central_v4_Prescl")) t->SetBranchAddress("HLT_HIFullTrack14_L1Central_v4_Prescl", &tHlts.HLT_HIFullTrack14_L1Central_v4_Prescl, &tHlts.b_HLT_HIFullTrack14_L1Central_v4_Prescl);
   if (t->GetBranch("HLT_HIFullTrack14_L1Peripheral_v4")) t->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v4", &tHlts.HLT_HIFullTrack14_L1Peripheral_v4, &tHlts.b_HLT_HIFullTrack14_L1Peripheral_v4);
   if (t->GetBranch("HLT_HIFullTrack14_L1Peripheral_v4_Prescl")) t->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v4_Prescl", &tHlts.HLT_HIFullTrack14_L1Peripheral_v4_Prescl, &tHlts.b_HLT_HIFullTrack14_L1Peripheral_v4_Prescl);
   if (t->GetBranch("HLT_HIFullTrack20_L1Central_v4")) t->SetBranchAddress("HLT_HIFullTrack20_L1Central_v4", &tHlts.HLT_HIFullTrack20_L1Central_v4, &tHlts.b_HLT_HIFullTrack20_L1Central_v4);
   if (t->GetBranch("HLT_HIFullTrack20_L1Central_v4_Prescl")) t->SetBranchAddress("HLT_HIFullTrack20_L1Central_v4_Prescl", &tHlts.HLT_HIFullTrack20_L1Central_v4_Prescl, &tHlts.b_HLT_HIFullTrack20_L1Central_v4_Prescl);
   if (t->GetBranch("HLT_HIFullTrack20_L1Peripheral_v4")) t->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v4", &tHlts.HLT_HIFullTrack20_L1Peripheral_v4, &tHlts.b_HLT_HIFullTrack20_L1Peripheral_v4);
   if (t->GetBranch("HLT_HIFullTrack20_L1Peripheral_v4_Prescl")) t->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v4_Prescl", &tHlts.HLT_HIFullTrack20_L1Peripheral_v4_Prescl, &tHlts.b_HLT_HIFullTrack20_L1Peripheral_v4_Prescl);
   if (t->GetBranch("HLT_HIFullTrack25_L1Central_v4")) t->SetBranchAddress("HLT_HIFullTrack25_L1Central_v4", &tHlts.HLT_HIFullTrack25_L1Central_v4, &tHlts.b_HLT_HIFullTrack25_L1Central_v4);
   if (t->GetBranch("HLT_HIFullTrack25_L1Central_v4_Prescl")) t->SetBranchAddress("HLT_HIFullTrack25_L1Central_v4_Prescl", &tHlts.HLT_HIFullTrack25_L1Central_v4_Prescl, &tHlts.b_HLT_HIFullTrack25_L1Central_v4_Prescl);
   if (t->GetBranch("HLT_HIFullTrack25_L1Peripheral_v4")) t->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v4", &tHlts.HLT_HIFullTrack25_L1Peripheral_v4, &tHlts.b_HLT_HIFullTrack25_L1Peripheral_v4);
   if (t->GetBranch("HLT_HIFullTrack25_L1Peripheral_v4_Prescl")) t->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v4_Prescl", &tHlts.HLT_HIFullTrack25_L1Peripheral_v4_Prescl, &tHlts.b_HLT_HIFullTrack25_L1Peripheral_v4_Prescl);
   if (t->GetBranch("HLT_HIJet55_v4")) t->SetBranchAddress("HLT_HIJet55_v4", &tHlts.HLT_HIJet55_v4, &tHlts.b_HLT_HIJet55_v4);
   if (t->GetBranch("HLT_HIJet55_v4_Prescl")) t->SetBranchAddress("HLT_HIJet55_v4_Prescl", &tHlts.HLT_HIJet55_v4_Prescl, &tHlts.b_HLT_HIJet55_v4_Prescl);
   if (t->GetBranch("HLT_HIJet65_Jet55_v4")) t->SetBranchAddress("HLT_HIJet65_Jet55_v4", &tHlts.HLT_HIJet65_Jet55_v4, &tHlts.b_HLT_HIJet65_Jet55_v4);
   if (t->GetBranch("HLT_HIJet65_Jet55_v4_Prescl")) t->SetBranchAddress("HLT_HIJet65_Jet55_v4_Prescl", &tHlts.HLT_HIJet65_Jet55_v4_Prescl, &tHlts.b_HLT_HIJet65_Jet55_v4_Prescl);
   if (t->GetBranch("HLT_HIJet65_v4")) t->SetBranchAddress("HLT_HIJet65_v4", &tHlts.HLT_HIJet65_v4, &tHlts.b_HLT_HIJet65_v4);
   if (t->GetBranch("HLT_HIJet65_v4_Prescl")) t->SetBranchAddress("HLT_HIJet65_v4_Prescl", &tHlts.HLT_HIJet65_v4_Prescl, &tHlts.b_HLT_HIJet65_v4_Prescl);
   if (t->GetBranch("HLT_HIJet80_v4")) t->SetBranchAddress("HLT_HIJet80_v4", &tHlts.HLT_HIJet80_v4, &tHlts.b_HLT_HIJet80_v4);
   if (t->GetBranch("HLT_HIJet80_v4_Prescl")) t->SetBranchAddress("HLT_HIJet80_v4_Prescl", &tHlts.HLT_HIJet80_v4_Prescl, &tHlts.b_HLT_HIJet80_v4_Prescl);
   if (t->GetBranch("HLT_HIJet95_v4")) t->SetBranchAddress("HLT_HIJet95_v4", &tHlts.HLT_HIJet95_v4, &tHlts.b_HLT_HIJet95_v4);
   if (t->GetBranch("HLT_HIJet95_v4_Prescl")) t->SetBranchAddress("HLT_HIJet95_v4_Prescl", &tHlts.HLT_HIJet95_v4_Prescl, &tHlts.b_HLT_HIJet95_v4_Prescl);
   if (t->GetBranch("HLT_HIJetE30_NoBPTX_v4")) t->SetBranchAddress("HLT_HIJetE30_NoBPTX_v4", &tHlts.HLT_HIJetE30_NoBPTX_v4, &tHlts.b_HLT_HIJetE30_NoBPTX_v4);
   if (t->GetBranch("HLT_HIJetE30_NoBPTX_v4_Prescl")) t->SetBranchAddress("HLT_HIJetE30_NoBPTX_v4_Prescl", &tHlts.HLT_HIJetE30_NoBPTX_v4_Prescl, &tHlts.b_HLT_HIJetE30_NoBPTX_v4_Prescl);
   if (t->GetBranch("HLT_HIJetE50_NoBPTX3BX_NoHalo_v4")) t->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v4", &tHlts.HLT_HIJetE50_NoBPTX3BX_NoHalo_v4, &tHlts.b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v4);
   if (t->GetBranch("HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl")) t->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl", &tHlts.HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl, &tHlts.b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl);
   if (t->GetBranch("HLT_HIMET120_v4")) t->SetBranchAddress("HLT_HIMET120_v4", &tHlts.HLT_HIMET120_v4, &tHlts.b_HLT_HIMET120_v4);
   if (t->GetBranch("HLT_HIMET120_v4_Prescl")) t->SetBranchAddress("HLT_HIMET120_v4_Prescl", &tHlts.HLT_HIMET120_v4_Prescl, &tHlts.b_HLT_HIMET120_v4_Prescl);
   if (t->GetBranch("HLT_HIMET200_v4")) t->SetBranchAddress("HLT_HIMET200_v4", &tHlts.HLT_HIMET200_v4, &tHlts.b_HLT_HIMET200_v4);
   if (t->GetBranch("HLT_HIMET200_v4_Prescl")) t->SetBranchAddress("HLT_HIMET200_v4_Prescl", &tHlts.HLT_HIMET200_v4_Prescl, &tHlts.b_HLT_HIMET200_v4_Prescl);
   if (t->GetBranch("HLT_HIMET220_v4")) t->SetBranchAddress("HLT_HIMET220_v4", &tHlts.HLT_HIMET220_v4, &tHlts.b_HLT_HIMET220_v4);
   if (t->GetBranch("HLT_HIMET220_v4_Prescl")) t->SetBranchAddress("HLT_HIMET220_v4_Prescl", &tHlts.HLT_HIMET220_v4_Prescl, &tHlts.b_HLT_HIMET220_v4_Prescl);
   if (t->GetBranch("HLT_HIPhoton10_Photon15_v4")) t->SetBranchAddress("HLT_HIPhoton10_Photon15_v4", &tHlts.HLT_HIPhoton10_Photon15_v4, &tHlts.b_HLT_HIPhoton10_Photon15_v4);
   if (t->GetBranch("HLT_HIPhoton10_Photon15_v4_Prescl")) t->SetBranchAddress("HLT_HIPhoton10_Photon15_v4_Prescl", &tHlts.HLT_HIPhoton10_Photon15_v4_Prescl, &tHlts.b_HLT_HIPhoton10_Photon15_v4_Prescl);
   if (t->GetBranch("HLT_HIPhoton15_Photon20_v4")) t->SetBranchAddress("HLT_HIPhoton15_Photon20_v4", &tHlts.HLT_HIPhoton15_Photon20_v4, &tHlts.b_HLT_HIPhoton15_Photon20_v4);
   if (t->GetBranch("HLT_HIPhoton15_Photon20_v4_Prescl")) t->SetBranchAddress("HLT_HIPhoton15_Photon20_v4_Prescl", &tHlts.HLT_HIPhoton15_Photon20_v4_Prescl, &tHlts.b_HLT_HIPhoton15_Photon20_v4_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton15_v4")) t->SetBranchAddress("HLT_HISinglePhoton15_v4", &tHlts.HLT_HISinglePhoton15_v4, &tHlts.b_HLT_HISinglePhoton15_v4);
   if (t->GetBranch("HLT_HISinglePhoton15_v4_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton15_v4_Prescl", &tHlts.HLT_HISinglePhoton15_v4_Prescl, &tHlts.b_HLT_HISinglePhoton15_v4_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton20_v4")) t->SetBranchAddress("HLT_HISinglePhoton20_v4", &tHlts.HLT_HISinglePhoton20_v4, &tHlts.b_HLT_HISinglePhoton20_v4);
   if (t->GetBranch("HLT_HISinglePhoton20_v4_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton20_v4_Prescl", &tHlts.HLT_HISinglePhoton20_v4_Prescl, &tHlts.b_HLT_HISinglePhoton20_v4_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton30_v4")) t->SetBranchAddress("HLT_HISinglePhoton30_v4", &tHlts.HLT_HISinglePhoton30_v4, &tHlts.b_HLT_HISinglePhoton30_v4);
   if (t->GetBranch("HLT_HISinglePhoton30_v4_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton30_v4_Prescl", &tHlts.HLT_HISinglePhoton30_v4_Prescl, &tHlts.b_HLT_HISinglePhoton30_v4_Prescl);
   if (t->GetBranch("HLT_HISinglePhoton40_v4")) t->SetBranchAddress("HLT_HISinglePhoton40_v4", &tHlts.HLT_HISinglePhoton40_v4, &tHlts.b_HLT_HISinglePhoton40_v4);
   if (t->GetBranch("HLT_HISinglePhoton40_v4_Prescl")) t->SetBranchAddress("HLT_HISinglePhoton40_v4_Prescl", &tHlts.HLT_HISinglePhoton40_v4_Prescl, &tHlts.b_HLT_HISinglePhoton40_v4_Prescl);
   if (t->GetBranch("HLT_HIBptxXOR_v4")) t->SetBranchAddress("HLT_HIBptxXOR_v4", &tHlts.HLT_HIBptxXOR_v4, &tHlts.b_HLT_HIBptxXOR_v4);
   if (t->GetBranch("HLT_HIBptxXOR_v4_Prescl")) t->SetBranchAddress("HLT_HIBptxXOR_v4_Prescl", &tHlts.HLT_HIBptxXOR_v4_Prescl, &tHlts.b_HLT_HIBptxXOR_v4_Prescl);
   if (t->GetBranch("HLT_HICentral10_v4")) t->SetBranchAddress("HLT_HICentral10_v4", &tHlts.HLT_HICentral10_v4, &tHlts.b_HLT_HICentral10_v4);
   if (t->GetBranch("HLT_HICentral10_v4_Prescl")) t->SetBranchAddress("HLT_HICentral10_v4_Prescl", &tHlts.HLT_HICentral10_v4_Prescl, &tHlts.b_HLT_HICentral10_v4_Prescl);
   if (t->GetBranch("HLT_HICentralityVeto_v4")) t->SetBranchAddress("HLT_HICentralityVeto_v4", &tHlts.HLT_HICentralityVeto_v4, &tHlts.b_HLT_HICentralityVeto_v4);
   if (t->GetBranch("HLT_HICentralityVeto_v4_Prescl")) t->SetBranchAddress("HLT_HICentralityVeto_v4_Prescl", &tHlts.HLT_HICentralityVeto_v4_Prescl, &tHlts.b_HLT_HICentralityVeto_v4_Prescl);
   if (t->GetBranch("HLT_HIL1Algo_BptxXOR_BSC_OR_v4")) t->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v4", &tHlts.HLT_HIL1Algo_BptxXOR_BSC_OR_v4, &tHlts.b_HLT_HIL1Algo_BptxXOR_BSC_OR_v4);
   if (t->GetBranch("HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl")) t->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl", &tHlts.HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl, &tHlts.b_HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasBSC_OR_v4")) t->SetBranchAddress("HLT_HIMinBiasBSC_OR_v4", &tHlts.HLT_HIMinBiasBSC_OR_v4, &tHlts.b_HLT_HIMinBiasBSC_OR_v4);
   if (t->GetBranch("HLT_HIMinBiasBSC_OR_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasBSC_OR_v4_Prescl", &tHlts.HLT_HIMinBiasBSC_OR_v4_Prescl, &tHlts.b_HLT_HIMinBiasBSC_OR_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasBSC_v4")) t->SetBranchAddress("HLT_HIMinBiasBSC_v4", &tHlts.HLT_HIMinBiasBSC_v4, &tHlts.b_HLT_HIMinBiasBSC_v4);
   if (t->GetBranch("HLT_HIMinBiasBSC_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasBSC_v4_Prescl", &tHlts.HLT_HIMinBiasBSC_v4_Prescl, &tHlts.b_HLT_HIMinBiasBSC_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHF_v4")) t->SetBranchAddress("HLT_HIMinBiasHF_v4", &tHlts.HLT_HIMinBiasHF_v4, &tHlts.b_HLT_HIMinBiasHF_v4);
   if (t->GetBranch("HLT_HIMinBiasHF_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHF_v4_Prescl", &tHlts.HLT_HIMinBiasHF_v4_Prescl, &tHlts.b_HLT_HIMinBiasHF_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHfOrBSC_v4")) t->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v4", &tHlts.HLT_HIMinBiasHfOrBSC_v4, &tHlts.b_HLT_HIMinBiasHfOrBSC_v4);
   if (t->GetBranch("HLT_HIMinBiasHfOrBSC_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v4_Prescl", &tHlts.HLT_HIMinBiasHfOrBSC_v4_Prescl, &tHlts.b_HLT_HIMinBiasHfOrBSC_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasHf_OR_v4")) t->SetBranchAddress("HLT_HIMinBiasHf_OR_v4", &tHlts.HLT_HIMinBiasHf_OR_v4, &tHlts.b_HLT_HIMinBiasHf_OR_v4);
   if (t->GetBranch("HLT_HIMinBiasHf_OR_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasHf_OR_v4_Prescl", &tHlts.HLT_HIMinBiasHf_OR_v4_Prescl, &tHlts.b_HLT_HIMinBiasHf_OR_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasPixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v4", &tHlts.HLT_HIMinBiasPixel_SingleTrack_v4, &tHlts.b_HLT_HIMinBiasPixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIMinBiasPixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIMinBiasPixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIMinBiasPixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDCPixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v4", &tHlts.HLT_HIMinBiasZDCPixel_SingleTrack_v4, &tHlts.b_HLT_HIMinBiasZDCPixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4", &tHlts.HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4, &tHlts.b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl", &tHlts.HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl, &tHlts.b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_v4")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v4", &tHlts.HLT_HIMinBiasZDC_Calo_v4, &tHlts.b_HLT_HIMinBiasZDC_Calo_v4);
   if (t->GetBranch("HLT_HIMinBiasZDC_Calo_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v4_Prescl", &tHlts.HLT_HIMinBiasZDC_Calo_v4_Prescl, &tHlts.b_HLT_HIMinBiasZDC_Calo_v4_Prescl);
   if (t->GetBranch("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4", &tHlts.HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4, &tHlts.b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIPhysics_v4")) t->SetBranchAddress("HLT_HIPhysics_v4", &tHlts.HLT_HIPhysics_v4, &tHlts.b_HLT_HIPhysics_v4);
   if (t->GetBranch("HLT_HIPhysics_v4_Prescl")) t->SetBranchAddress("HLT_HIPhysics_v4_Prescl", &tHlts.HLT_HIPhysics_v4_Prescl, &tHlts.b_HLT_HIPhysics_v4_Prescl);
   if (t->GetBranch("HLT_HIRandom_v4")) t->SetBranchAddress("HLT_HIRandom_v4", &tHlts.HLT_HIRandom_v4, &tHlts.b_HLT_HIRandom_v4);
   if (t->GetBranch("HLT_HIRandom_v4_Prescl")) t->SetBranchAddress("HLT_HIRandom_v4_Prescl", &tHlts.HLT_HIRandom_v4_Prescl, &tHlts.b_HLT_HIRandom_v4_Prescl);
   if (t->GetBranch("HLT_HIUCC010_v4")) t->SetBranchAddress("HLT_HIUCC010_v4", &tHlts.HLT_HIUCC010_v4, &tHlts.b_HLT_HIUCC010_v4);
   if (t->GetBranch("HLT_HIUCC010_v4_Prescl")) t->SetBranchAddress("HLT_HIUCC010_v4_Prescl", &tHlts.HLT_HIUCC010_v4_Prescl, &tHlts.b_HLT_HIUCC010_v4_Prescl);
   if (t->GetBranch("HLT_HIUCC015_v4")) t->SetBranchAddress("HLT_HIUCC015_v4", &tHlts.HLT_HIUCC015_v4, &tHlts.b_HLT_HIUCC015_v4);
   if (t->GetBranch("HLT_HIUCC015_v4_Prescl")) t->SetBranchAddress("HLT_HIUCC015_v4_Prescl", &tHlts.HLT_HIUCC015_v4_Prescl, &tHlts.b_HLT_HIUCC015_v4_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuEG2Pixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v4", &tHlts.HLT_HIUPCNeuEG2Pixel_SingleTrack_v4, &tHlts.b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuEG5Pixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v4", &tHlts.HLT_HIUPCNeuEG5Pixel_SingleTrack_v4, &tHlts.b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4", &tHlts.HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4, &tHlts.b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4", &tHlts.HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4, &tHlts.b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4", &tHlts.HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4, &tHlts.b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIUPCNeuMuPixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v4", &tHlts.HLT_HIUPCNeuMuPixel_SingleTrack_v4, &tHlts.b_HLT_HIUPCNeuMuPixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIZeroBiasPixel_SingleTrack_v4")) t->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v4", &tHlts.HLT_HIZeroBiasPixel_SingleTrack_v4, &tHlts.b_HLT_HIZeroBiasPixel_SingleTrack_v4);
   if (t->GetBranch("HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl")) t->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl", &tHlts.HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl, &tHlts.b_HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl);
   if (t->GetBranch("HLT_HIZeroBiasXOR_v4")) t->SetBranchAddress("HLT_HIZeroBiasXOR_v4", &tHlts.HLT_HIZeroBiasXOR_v4, &tHlts.b_HLT_HIZeroBiasXOR_v4);
   if (t->GetBranch("HLT_HIZeroBiasXOR_v4_Prescl")) t->SetBranchAddress("HLT_HIZeroBiasXOR_v4_Prescl", &tHlts.HLT_HIZeroBiasXOR_v4_Prescl, &tHlts.b_HLT_HIZeroBiasXOR_v4_Prescl);
   if (t->GetBranch("HLT_HIZeroBias_v4")) t->SetBranchAddress("HLT_HIZeroBias_v4", &tHlts.HLT_HIZeroBias_v4, &tHlts.b_HLT_HIZeroBias_v4);
   if (t->GetBranch("HLT_HIZeroBias_v4_Prescl")) t->SetBranchAddress("HLT_HIZeroBias_v4_Prescl", &tHlts.HLT_HIZeroBias_v4_Prescl, &tHlts.b_HLT_HIZeroBias_v4_Prescl);
   if (t->GetBranch("L1_BptxMinus")) t->SetBranchAddress("L1_BptxMinus", &tHlts.L1_BptxMinus, &tHlts.b_L1_BptxMinus);
   if (t->GetBranch("L1_BptxMinus_Prescl")) t->SetBranchAddress("L1_BptxMinus_Prescl", &tHlts.L1_BptxMinus_Prescl, &tHlts.b_L1_BptxMinus_Prescl);
   if (t->GetBranch("L1_BptxMinus_NotBptxPlus")) t->SetBranchAddress("L1_BptxMinus_NotBptxPlus", &tHlts.L1_BptxMinus_NotBptxPlus, &tHlts.b_L1_BptxMinus_NotBptxPlus);
   if (t->GetBranch("L1_BptxMinus_NotBptxPlus_Prescl")) t->SetBranchAddress("L1_BptxMinus_NotBptxPlus_Prescl", &tHlts.L1_BptxMinus_NotBptxPlus_Prescl, &tHlts.b_L1_BptxMinus_NotBptxPlus_Prescl);
   if (t->GetBranch("L1_BptxPlus")) t->SetBranchAddress("L1_BptxPlus", &tHlts.L1_BptxPlus, &tHlts.b_L1_BptxPlus);
   if (t->GetBranch("L1_BptxPlus_Prescl")) t->SetBranchAddress("L1_BptxPlus_Prescl", &tHlts.L1_BptxPlus_Prescl, &tHlts.b_L1_BptxPlus_Prescl);
   if (t->GetBranch("L1_BptxPlusANDMinus")) t->SetBranchAddress("L1_BptxPlusANDMinus", &tHlts.L1_BptxPlusANDMinus, &tHlts.b_L1_BptxPlusANDMinus);
   if (t->GetBranch("L1_BptxPlusANDMinus_Prescl")) t->SetBranchAddress("L1_BptxPlusANDMinus_Prescl", &tHlts.L1_BptxPlusANDMinus_Prescl, &tHlts.b_L1_BptxPlusANDMinus_Prescl);
   if (t->GetBranch("L1_BptxPlusORMinus")) t->SetBranchAddress("L1_BptxPlusORMinus", &tHlts.L1_BptxPlusORMinus, &tHlts.b_L1_BptxPlusORMinus);
   if (t->GetBranch("L1_BptxPlusORMinus_Prescl")) t->SetBranchAddress("L1_BptxPlusORMinus_Prescl", &tHlts.L1_BptxPlusORMinus_Prescl, &tHlts.b_L1_BptxPlusORMinus_Prescl);
   if (t->GetBranch("L1_BptxPlusORMinus_instance1")) t->SetBranchAddress("L1_BptxPlusORMinus_instance1", &tHlts.L1_BptxPlusORMinus_instance1, &tHlts.b_L1_BptxPlusORMinus_instance1);
   if (t->GetBranch("L1_BptxPlusORMinus_instance1_Prescl")) t->SetBranchAddress("L1_BptxPlusORMinus_instance1_Prescl", &tHlts.L1_BptxPlusORMinus_instance1_Prescl, &tHlts.b_L1_BptxPlusORMinus_instance1_Prescl);
   if (t->GetBranch("L1_BptxPlus_NotBptxMinus")) t->SetBranchAddress("L1_BptxPlus_NotBptxMinus", &tHlts.L1_BptxPlus_NotBptxMinus, &tHlts.b_L1_BptxPlus_NotBptxMinus);
   if (t->GetBranch("L1_BptxPlus_NotBptxMinus_Prescl")) t->SetBranchAddress("L1_BptxPlus_NotBptxMinus_Prescl", &tHlts.L1_BptxPlus_NotBptxMinus_Prescl, &tHlts.b_L1_BptxPlus_NotBptxMinus_Prescl);
   if (t->GetBranch("L1_BptxXOR")) t->SetBranchAddress("L1_BptxXOR", &tHlts.L1_BptxXOR, &tHlts.b_L1_BptxXOR);
   if (t->GetBranch("L1_BptxXOR_Prescl")) t->SetBranchAddress("L1_BptxXOR_Prescl", &tHlts.L1_BptxXOR_Prescl, &tHlts.b_L1_BptxXOR_Prescl);
   if (t->GetBranch("L1_BptxXOR_BscMinBiasOR")) t->SetBranchAddress("L1_BptxXOR_BscMinBiasOR", &tHlts.L1_BptxXOR_BscMinBiasOR, &tHlts.b_L1_BptxXOR_BscMinBiasOR);
   if (t->GetBranch("L1_BptxXOR_BscMinBiasOR_Prescl")) t->SetBranchAddress("L1_BptxXOR_BscMinBiasOR_Prescl", &tHlts.L1_BptxXOR_BscMinBiasOR_Prescl, &tHlts.b_L1_BptxXOR_BscMinBiasOR_Prescl);
   if (t->GetBranch("L1_Bsc2Minus_BptxMinus")) t->SetBranchAddress("L1_Bsc2Minus_BptxMinus", &tHlts.L1_Bsc2Minus_BptxMinus, &tHlts.b_L1_Bsc2Minus_BptxMinus);
   if (t->GetBranch("L1_Bsc2Minus_BptxMinus_Prescl")) t->SetBranchAddress("L1_Bsc2Minus_BptxMinus_Prescl", &tHlts.L1_Bsc2Minus_BptxMinus_Prescl, &tHlts.b_L1_Bsc2Minus_BptxMinus_Prescl);
   if (t->GetBranch("L1_Bsc2Plus_BptxPlus")) t->SetBranchAddress("L1_Bsc2Plus_BptxPlus", &tHlts.L1_Bsc2Plus_BptxPlus, &tHlts.b_L1_Bsc2Plus_BptxPlus);
   if (t->GetBranch("L1_Bsc2Plus_BptxPlus_Prescl")) t->SetBranchAddress("L1_Bsc2Plus_BptxPlus_Prescl", &tHlts.L1_Bsc2Plus_BptxPlus_Prescl, &tHlts.b_L1_Bsc2Plus_BptxPlus_Prescl);
   if (t->GetBranch("L1_BscMinBiasInnerThreshold1")) t->SetBranchAddress("L1_BscMinBiasInnerThreshold1", &tHlts.L1_BscMinBiasInnerThreshold1, &tHlts.b_L1_BscMinBiasInnerThreshold1);
   if (t->GetBranch("L1_BscMinBiasInnerThreshold1_Prescl")) t->SetBranchAddress("L1_BscMinBiasInnerThreshold1_Prescl", &tHlts.L1_BscMinBiasInnerThreshold1_Prescl, &tHlts.b_L1_BscMinBiasInnerThreshold1_Prescl);
   if (t->GetBranch("L1_BscMinBiasInnerThreshold1_BptxAND")) t->SetBranchAddress("L1_BscMinBiasInnerThreshold1_BptxAND", &tHlts.L1_BscMinBiasInnerThreshold1_BptxAND, &tHlts.b_L1_BscMinBiasInnerThreshold1_BptxAND);
   if (t->GetBranch("L1_BscMinBiasInnerThreshold1_BptxAND_Prescl")) t->SetBranchAddress("L1_BscMinBiasInnerThreshold1_BptxAND_Prescl", &tHlts.L1_BscMinBiasInnerThreshold1_BptxAND_Prescl, &tHlts.b_L1_BscMinBiasInnerThreshold1_BptxAND_Prescl);
   if (t->GetBranch("L1_BscMinBiasInnerThreshold2_BptxAND")) t->SetBranchAddress("L1_BscMinBiasInnerThreshold2_BptxAND", &tHlts.L1_BscMinBiasInnerThreshold2_BptxAND, &tHlts.b_L1_BscMinBiasInnerThreshold2_BptxAND);
   if (t->GetBranch("L1_BscMinBiasInnerThreshold2_BptxAND_Prescl")) t->SetBranchAddress("L1_BscMinBiasInnerThreshold2_BptxAND_Prescl", &tHlts.L1_BscMinBiasInnerThreshold2_BptxAND_Prescl, &tHlts.b_L1_BscMinBiasInnerThreshold2_BptxAND_Prescl);
   if (t->GetBranch("L1_BscMinBiasOR_BptxAND")) t->SetBranchAddress("L1_BscMinBiasOR_BptxAND", &tHlts.L1_BscMinBiasOR_BptxAND, &tHlts.b_L1_BscMinBiasOR_BptxAND);
   if (t->GetBranch("L1_BscMinBiasOR_BptxAND_Prescl")) t->SetBranchAddress("L1_BscMinBiasOR_BptxAND_Prescl", &tHlts.L1_BscMinBiasOR_BptxAND_Prescl, &tHlts.b_L1_BscMinBiasOR_BptxAND_Prescl);
   if (t->GetBranch("L1_BscMinBiasOR_BptxPlusORMinus")) t->SetBranchAddress("L1_BscMinBiasOR_BptxPlusORMinus", &tHlts.L1_BscMinBiasOR_BptxPlusORMinus, &tHlts.b_L1_BscMinBiasOR_BptxPlusORMinus);
   if (t->GetBranch("L1_BscMinBiasOR_BptxPlusORMinus_Prescl")) t->SetBranchAddress("L1_BscMinBiasOR_BptxPlusORMinus_Prescl", &tHlts.L1_BscMinBiasOR_BptxPlusORMinus_Prescl, &tHlts.b_L1_BscMinBiasOR_BptxPlusORMinus_Prescl);
   if (t->GetBranch("L1_BscMinBiasThreshold1")) t->SetBranchAddress("L1_BscMinBiasThreshold1", &tHlts.L1_BscMinBiasThreshold1, &tHlts.b_L1_BscMinBiasThreshold1);
   if (t->GetBranch("L1_BscMinBiasThreshold1_Prescl")) t->SetBranchAddress("L1_BscMinBiasThreshold1_Prescl", &tHlts.L1_BscMinBiasThreshold1_Prescl, &tHlts.b_L1_BscMinBiasThreshold1_Prescl);
   if (t->GetBranch("L1_BscMinBiasThreshold1_BptxAND")) t->SetBranchAddress("L1_BscMinBiasThreshold1_BptxAND", &tHlts.L1_BscMinBiasThreshold1_BptxAND, &tHlts.b_L1_BscMinBiasThreshold1_BptxAND);
   if (t->GetBranch("L1_BscMinBiasThreshold1_BptxAND_Prescl")) t->SetBranchAddress("L1_BscMinBiasThreshold1_BptxAND_Prescl", &tHlts.L1_BscMinBiasThreshold1_BptxAND_Prescl, &tHlts.b_L1_BscMinBiasThreshold1_BptxAND_Prescl);
   if (t->GetBranch("L1_BscMinBiasThreshold2_BptxAND")) t->SetBranchAddress("L1_BscMinBiasThreshold2_BptxAND", &tHlts.L1_BscMinBiasThreshold2_BptxAND, &tHlts.b_L1_BscMinBiasThreshold2_BptxAND);
   if (t->GetBranch("L1_BscMinBiasThreshold2_BptxAND_Prescl")) t->SetBranchAddress("L1_BscMinBiasThreshold2_BptxAND_Prescl", &tHlts.L1_BscMinBiasThreshold2_BptxAND_Prescl, &tHlts.b_L1_BscMinBiasThreshold2_BptxAND_Prescl);
   if (t->GetBranch("L1_DoubleEG10")) t->SetBranchAddress("L1_DoubleEG10", &tHlts.L1_DoubleEG10, &tHlts.b_L1_DoubleEG10);
   if (t->GetBranch("L1_DoubleEG10_Prescl")) t->SetBranchAddress("L1_DoubleEG10_Prescl", &tHlts.L1_DoubleEG10_Prescl, &tHlts.b_L1_DoubleEG10_Prescl);
   if (t->GetBranch("L1_DoubleEG3")) t->SetBranchAddress("L1_DoubleEG3", &tHlts.L1_DoubleEG3, &tHlts.b_L1_DoubleEG3);
   if (t->GetBranch("L1_DoubleEG3_Prescl")) t->SetBranchAddress("L1_DoubleEG3_Prescl", &tHlts.L1_DoubleEG3_Prescl, &tHlts.b_L1_DoubleEG3_Prescl);
   if (t->GetBranch("L1_DoubleEG5_BptxAND")) t->SetBranchAddress("L1_DoubleEG5_BptxAND", &tHlts.L1_DoubleEG5_BptxAND, &tHlts.b_L1_DoubleEG5_BptxAND);
   if (t->GetBranch("L1_DoubleEG5_BptxAND_Prescl")) t->SetBranchAddress("L1_DoubleEG5_BptxAND_Prescl", &tHlts.L1_DoubleEG5_BptxAND_Prescl, &tHlts.b_L1_DoubleEG5_BptxAND_Prescl);
   if (t->GetBranch("L1_DoubleEG8_BptxAND")) t->SetBranchAddress("L1_DoubleEG8_BptxAND", &tHlts.L1_DoubleEG8_BptxAND, &tHlts.b_L1_DoubleEG8_BptxAND);
   if (t->GetBranch("L1_DoubleEG8_BptxAND_Prescl")) t->SetBranchAddress("L1_DoubleEG8_BptxAND_Prescl", &tHlts.L1_DoubleEG8_BptxAND_Prescl, &tHlts.b_L1_DoubleEG8_BptxAND_Prescl);
   if (t->GetBranch("L1_DoubleForJet32_EtaOpp")) t->SetBranchAddress("L1_DoubleForJet32_EtaOpp", &tHlts.L1_DoubleForJet32_EtaOpp, &tHlts.b_L1_DoubleForJet32_EtaOpp);
   if (t->GetBranch("L1_DoubleForJet32_EtaOpp_Prescl")) t->SetBranchAddress("L1_DoubleForJet32_EtaOpp_Prescl", &tHlts.L1_DoubleForJet32_EtaOpp_Prescl, &tHlts.b_L1_DoubleForJet32_EtaOpp_Prescl);
   if (t->GetBranch("L1_DoubleForJet44_EtaOpp")) t->SetBranchAddress("L1_DoubleForJet44_EtaOpp", &tHlts.L1_DoubleForJet44_EtaOpp, &tHlts.b_L1_DoubleForJet44_EtaOpp);
   if (t->GetBranch("L1_DoubleForJet44_EtaOpp_Prescl")) t->SetBranchAddress("L1_DoubleForJet44_EtaOpp_Prescl", &tHlts.L1_DoubleForJet44_EtaOpp_Prescl, &tHlts.b_L1_DoubleForJet44_EtaOpp_Prescl);
   if (t->GetBranch("L1_DoubleHfBitCountsRing1_P1N1")) t->SetBranchAddress("L1_DoubleHfBitCountsRing1_P1N1", &tHlts.L1_DoubleHfBitCountsRing1_P1N1, &tHlts.b_L1_DoubleHfBitCountsRing1_P1N1);
   if (t->GetBranch("L1_DoubleHfBitCountsRing1_P1N1_Prescl")) t->SetBranchAddress("L1_DoubleHfBitCountsRing1_P1N1_Prescl", &tHlts.L1_DoubleHfBitCountsRing1_P1N1_Prescl, &tHlts.b_L1_DoubleHfBitCountsRing1_P1N1_Prescl);
   if (t->GetBranch("L1_DoubleHfBitCountsRing2_P1N1")) t->SetBranchAddress("L1_DoubleHfBitCountsRing2_P1N1", &tHlts.L1_DoubleHfBitCountsRing2_P1N1, &tHlts.b_L1_DoubleHfBitCountsRing2_P1N1);
   if (t->GetBranch("L1_DoubleHfBitCountsRing2_P1N1_Prescl")) t->SetBranchAddress("L1_DoubleHfBitCountsRing2_P1N1_Prescl", &tHlts.L1_DoubleHfBitCountsRing2_P1N1_Prescl, &tHlts.b_L1_DoubleHfBitCountsRing2_P1N1_Prescl);
   if (t->GetBranch("L1_DoubleHfRingEtSumsRing1_P200N200")) t->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P200N200", &tHlts.L1_DoubleHfRingEtSumsRing1_P200N200, &tHlts.b_L1_DoubleHfRingEtSumsRing1_P200N200);
   if (t->GetBranch("L1_DoubleHfRingEtSumsRing1_P200N200_Prescl")) t->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P200N200_Prescl", &tHlts.L1_DoubleHfRingEtSumsRing1_P200N200_Prescl, &tHlts.b_L1_DoubleHfRingEtSumsRing1_P200N200_Prescl);
   if (t->GetBranch("L1_DoubleHfRingEtSumsRing1_P4N4")) t->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P4N4", &tHlts.L1_DoubleHfRingEtSumsRing1_P4N4, &tHlts.b_L1_DoubleHfRingEtSumsRing1_P4N4);
   if (t->GetBranch("L1_DoubleHfRingEtSumsRing1_P4N4_Prescl")) t->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P4N4_Prescl", &tHlts.L1_DoubleHfRingEtSumsRing1_P4N4_Prescl, &tHlts.b_L1_DoubleHfRingEtSumsRing1_P4N4_Prescl);
   if (t->GetBranch("L1_DoubleHfRingEtSumsRing2_P200N200")) t->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P200N200", &tHlts.L1_DoubleHfRingEtSumsRing2_P200N200, &tHlts.b_L1_DoubleHfRingEtSumsRing2_P200N200);
   if (t->GetBranch("L1_DoubleHfRingEtSumsRing2_P200N200_Prescl")) t->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P200N200_Prescl", &tHlts.L1_DoubleHfRingEtSumsRing2_P200N200_Prescl, &tHlts.b_L1_DoubleHfRingEtSumsRing2_P200N200_Prescl);
   if (t->GetBranch("L1_DoubleHfRingEtSumsRing2_P4N4")) t->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P4N4", &tHlts.L1_DoubleHfRingEtSumsRing2_P4N4, &tHlts.b_L1_DoubleHfRingEtSumsRing2_P4N4);
   if (t->GetBranch("L1_DoubleHfRingEtSumsRing2_P4N4_Prescl")) t->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P4N4_Prescl", &tHlts.L1_DoubleHfRingEtSumsRing2_P4N4_Prescl, &tHlts.b_L1_DoubleHfRingEtSumsRing2_P4N4_Prescl);
   if (t->GetBranch("L1_DoubleMu3_BptxAND")) t->SetBranchAddress("L1_DoubleMu3_BptxAND", &tHlts.L1_DoubleMu3_BptxAND, &tHlts.b_L1_DoubleMu3_BptxAND);
   if (t->GetBranch("L1_DoubleMu3_BptxAND_Prescl")) t->SetBranchAddress("L1_DoubleMu3_BptxAND_Prescl", &tHlts.L1_DoubleMu3_BptxAND_Prescl, &tHlts.b_L1_DoubleMu3_BptxAND_Prescl);
   if (t->GetBranch("L1_DoubleMuOpen_BptxAND")) t->SetBranchAddress("L1_DoubleMuOpen_BptxAND", &tHlts.L1_DoubleMuOpen_BptxAND, &tHlts.b_L1_DoubleMuOpen_BptxAND);
   if (t->GetBranch("L1_DoubleMuOpen_BptxAND_Prescl")) t->SetBranchAddress("L1_DoubleMuOpen_BptxAND_Prescl", &tHlts.L1_DoubleMuOpen_BptxAND_Prescl, &tHlts.b_L1_DoubleMuOpen_BptxAND_Prescl);
   if (t->GetBranch("L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND")) t->SetBranchAddress("L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND", &tHlts.L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND, &tHlts.b_L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND);
   if (t->GetBranch("L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl")) t->SetBranchAddress("L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl", &tHlts.L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl, &tHlts.b_L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl);
   if (t->GetBranch("L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND")) t->SetBranchAddress("L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND", &tHlts.L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND, &tHlts.b_L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND);
   if (t->GetBranch("L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl")) t->SetBranchAddress("L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl", &tHlts.L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl, &tHlts.b_L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl);
   if (t->GetBranch("L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND")) t->SetBranchAddress("L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND", &tHlts.L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND, &tHlts.b_L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND);
   if (t->GetBranch("L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl")) t->SetBranchAddress("L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl", &tHlts.L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl, &tHlts.b_L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl);
   if (t->GetBranch("L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND")) t->SetBranchAddress("L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND", &tHlts.L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND, &tHlts.b_L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND);
   if (t->GetBranch("L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl")) t->SetBranchAddress("L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl", &tHlts.L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl, &tHlts.b_L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl);
   if (t->GetBranch("L1_ETM30_BptxAND")) t->SetBranchAddress("L1_ETM30_BptxAND", &tHlts.L1_ETM30_BptxAND, &tHlts.b_L1_ETM30_BptxAND);
   if (t->GetBranch("L1_ETM30_BptxAND_Prescl")) t->SetBranchAddress("L1_ETM30_BptxAND_Prescl", &tHlts.L1_ETM30_BptxAND_Prescl, &tHlts.b_L1_ETM30_BptxAND_Prescl);
   if (t->GetBranch("L1_ETM50_BptxAND")) t->SetBranchAddress("L1_ETM50_BptxAND", &tHlts.L1_ETM50_BptxAND, &tHlts.b_L1_ETM50_BptxAND);
   if (t->GetBranch("L1_ETM50_BptxAND_Prescl")) t->SetBranchAddress("L1_ETM50_BptxAND_Prescl", &tHlts.L1_ETM50_BptxAND_Prescl, &tHlts.b_L1_ETM50_BptxAND_Prescl);
   if (t->GetBranch("L1_ETT100_BptxAND")) t->SetBranchAddress("L1_ETT100_BptxAND", &tHlts.L1_ETT100_BptxAND, &tHlts.b_L1_ETT100_BptxAND);
   if (t->GetBranch("L1_ETT100_BptxAND_Prescl")) t->SetBranchAddress("L1_ETT100_BptxAND_Prescl", &tHlts.L1_ETT100_BptxAND_Prescl, &tHlts.b_L1_ETT100_BptxAND_Prescl);
   if (t->GetBranch("L1_ETT140_BptxAND")) t->SetBranchAddress("L1_ETT140_BptxAND", &tHlts.L1_ETT140_BptxAND, &tHlts.b_L1_ETT140_BptxAND);
   if (t->GetBranch("L1_ETT140_BptxAND_Prescl")) t->SetBranchAddress("L1_ETT140_BptxAND_Prescl", &tHlts.L1_ETT140_BptxAND_Prescl, &tHlts.b_L1_ETT140_BptxAND_Prescl);
   if (t->GetBranch("L1_ETT2000")) t->SetBranchAddress("L1_ETT2000", &tHlts.L1_ETT2000, &tHlts.b_L1_ETT2000);
   if (t->GetBranch("L1_ETT2000_Prescl")) t->SetBranchAddress("L1_ETT2000_Prescl", &tHlts.L1_ETT2000_Prescl, &tHlts.b_L1_ETT2000_Prescl);
   if (t->GetBranch("L1_ETT220_BptxAND")) t->SetBranchAddress("L1_ETT220_BptxAND", &tHlts.L1_ETT220_BptxAND, &tHlts.b_L1_ETT220_BptxAND);
   if (t->GetBranch("L1_ETT220_BptxAND_Prescl")) t->SetBranchAddress("L1_ETT220_BptxAND_Prescl", &tHlts.L1_ETT220_BptxAND_Prescl, &tHlts.b_L1_ETT220_BptxAND_Prescl);
   if (t->GetBranch("L1_ETT70_BptxAND")) t->SetBranchAddress("L1_ETT70_BptxAND", &tHlts.L1_ETT70_BptxAND, &tHlts.b_L1_ETT70_BptxAND);
   if (t->GetBranch("L1_ETT70_BptxAND_Prescl")) t->SetBranchAddress("L1_ETT70_BptxAND_Prescl", &tHlts.L1_ETT70_BptxAND_Prescl, &tHlts.b_L1_ETT70_BptxAND_Prescl);
   if (t->GetBranch("L1_ETT800")) t->SetBranchAddress("L1_ETT800", &tHlts.L1_ETT800, &tHlts.b_L1_ETT800);
   if (t->GetBranch("L1_ETT800_Prescl")) t->SetBranchAddress("L1_ETT800_Prescl", &tHlts.L1_ETT800_Prescl, &tHlts.b_L1_ETT800_Prescl);
   if (t->GetBranch("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND")) t->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND", &tHlts.L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND, &tHlts.b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND);
   if (t->GetBranch("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl")) t->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl", &tHlts.L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl, &tHlts.b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl);
   if (t->GetBranch("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1")) t->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1", &tHlts.L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1, &tHlts.b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1);
   if (t->GetBranch("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl")) t->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl", &tHlts.L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl, &tHlts.b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl);
   if (t->GetBranch("L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND")) t->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND", &tHlts.L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND, &tHlts.b_L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND);
   if (t->GetBranch("L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl")) t->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl", &tHlts.L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl, &tHlts.b_L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl);
   if (t->GetBranch("L1_HcalHfCoincidencePm")) t->SetBranchAddress("L1_HcalHfCoincidencePm", &tHlts.L1_HcalHfCoincidencePm, &tHlts.b_L1_HcalHfCoincidencePm);
   if (t->GetBranch("L1_HcalHfCoincidencePm_Prescl")) t->SetBranchAddress("L1_HcalHfCoincidencePm_Prescl", &tHlts.L1_HcalHfCoincidencePm_Prescl, &tHlts.b_L1_HcalHfCoincidencePm_Prescl);
   if (t->GetBranch("L1_HcalHfCoincidencePm_BptxAND")) t->SetBranchAddress("L1_HcalHfCoincidencePm_BptxAND", &tHlts.L1_HcalHfCoincidencePm_BptxAND, &tHlts.b_L1_HcalHfCoincidencePm_BptxAND);
   if (t->GetBranch("L1_HcalHfCoincidencePm_BptxAND_Prescl")) t->SetBranchAddress("L1_HcalHfCoincidencePm_BptxAND_Prescl", &tHlts.L1_HcalHfCoincidencePm_BptxAND_Prescl, &tHlts.b_L1_HcalHfCoincidencePm_BptxAND_Prescl);
   if (t->GetBranch("L1_HcalHfMmOrPpOrPm")) t->SetBranchAddress("L1_HcalHfMmOrPpOrPm", &tHlts.L1_HcalHfMmOrPpOrPm, &tHlts.b_L1_HcalHfMmOrPpOrPm);
   if (t->GetBranch("L1_HcalHfMmOrPpOrPm_Prescl")) t->SetBranchAddress("L1_HcalHfMmOrPpOrPm_Prescl", &tHlts.L1_HcalHfMmOrPpOrPm_Prescl, &tHlts.b_L1_HcalHfMmOrPpOrPm_Prescl);
   if (t->GetBranch("L1_HcalHfMmOrPpOrPm_BptxAND")) t->SetBranchAddress("L1_HcalHfMmOrPpOrPm_BptxAND", &tHlts.L1_HcalHfMmOrPpOrPm_BptxAND, &tHlts.b_L1_HcalHfMmOrPpOrPm_BptxAND);
   if (t->GetBranch("L1_HcalHfMmOrPpOrPm_BptxAND_Prescl")) t->SetBranchAddress("L1_HcalHfMmOrPpOrPm_BptxAND_Prescl", &tHlts.L1_HcalHfMmOrPpOrPm_BptxAND_Prescl, &tHlts.b_L1_HcalHfMmOrPpOrPm_BptxAND_Prescl);
   if (t->GetBranch("L1_HcalHfMmpOrMpp")) t->SetBranchAddress("L1_HcalHfMmpOrMpp", &tHlts.L1_HcalHfMmpOrMpp, &tHlts.b_L1_HcalHfMmpOrMpp);
   if (t->GetBranch("L1_HcalHfMmpOrMpp_Prescl")) t->SetBranchAddress("L1_HcalHfMmpOrMpp_Prescl", &tHlts.L1_HcalHfMmpOrMpp_Prescl, &tHlts.b_L1_HcalHfMmpOrMpp_Prescl);
   if (t->GetBranch("L1_HcalHfMmpOrMpp_BptxAND")) t->SetBranchAddress("L1_HcalHfMmpOrMpp_BptxAND", &tHlts.L1_HcalHfMmpOrMpp_BptxAND, &tHlts.b_L1_HcalHfMmpOrMpp_BptxAND);
   if (t->GetBranch("L1_HcalHfMmpOrMpp_BptxAND_Prescl")) t->SetBranchAddress("L1_HcalHfMmpOrMpp_BptxAND_Prescl", &tHlts.L1_HcalHfMmpOrMpp_BptxAND_Prescl, &tHlts.b_L1_HcalHfMmpOrMpp_BptxAND_Prescl);
   if (t->GetBranch("L1_HcalHoTotalOR")) t->SetBranchAddress("L1_HcalHoTotalOR", &tHlts.L1_HcalHoTotalOR, &tHlts.b_L1_HcalHoTotalOR);
   if (t->GetBranch("L1_HcalHoTotalOR_Prescl")) t->SetBranchAddress("L1_HcalHoTotalOR_Prescl", &tHlts.L1_HcalHoTotalOR_Prescl, &tHlts.b_L1_HcalHoTotalOR_Prescl);
   if (t->GetBranch("L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND")) t->SetBranchAddress("L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND", &tHlts.L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND, &tHlts.b_L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND);
   if (t->GetBranch("L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl")) t->SetBranchAddress("L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl", &tHlts.L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl, &tHlts.b_L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl);
   if (t->GetBranch("L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND")) t->SetBranchAddress("L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND", &tHlts.L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND, &tHlts.b_L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND);
   if (t->GetBranch("L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl")) t->SetBranchAddress("L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl", &tHlts.L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl, &tHlts.b_L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasOR")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasOR", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasOR, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasOR);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasThreshold1")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasThreshold1", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasThreshold1, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasThreshold1);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasThreshold2")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasThreshold2", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasThreshold2, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasThreshold2);
   if (t->GetBranch("L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl")) t->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl", &tHlts.L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl, &tHlts.b_L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl);
   if (t->GetBranch("L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1")) t->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1", &tHlts.L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1, &tHlts.b_L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1);
   if (t->GetBranch("L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl")) t->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl", &tHlts.L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl, &tHlts.b_L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl);
   if (t->GetBranch("L1_NotBsc2_BptxAND_HcalHfCoincidencePm")) t->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfCoincidencePm", &tHlts.L1_NotBsc2_BptxAND_HcalHfCoincidencePm, &tHlts.b_L1_NotBsc2_BptxAND_HcalHfCoincidencePm);
   if (t->GetBranch("L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl")) t->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl", &tHlts.L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl, &tHlts.b_L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl);
   if (t->GetBranch("L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm")) t->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm", &tHlts.L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm, &tHlts.b_L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm);
   if (t->GetBranch("L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl")) t->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl", &tHlts.L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl, &tHlts.b_L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl);
   if (t->GetBranch("L1_NotBsc2_BptxAND_HcalHfMmpOrMpp")) t->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfMmpOrMpp", &tHlts.L1_NotBsc2_BptxAND_HcalHfMmpOrMpp, &tHlts.b_L1_NotBsc2_BptxAND_HcalHfMmpOrMpp);
   if (t->GetBranch("L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl")) t->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl", &tHlts.L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl, &tHlts.b_L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasInnerThreshold1")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasInnerThreshold1", &tHlts.L1_NotBsc2_BscMinBiasInnerThreshold1, &tHlts.b_L1_NotBsc2_BscMinBiasInnerThreshold1);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl", &tHlts.L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl, &tHlts.b_L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasInnerThreshold2")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasInnerThreshold2", &tHlts.L1_NotBsc2_BscMinBiasInnerThreshold2, &tHlts.b_L1_NotBsc2_BscMinBiasInnerThreshold2);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl", &tHlts.L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl, &tHlts.b_L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasOR")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasOR", &tHlts.L1_NotBsc2_BscMinBiasOR, &tHlts.b_L1_NotBsc2_BscMinBiasOR);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasOR_Prescl")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasOR_Prescl", &tHlts.L1_NotBsc2_BscMinBiasOR_Prescl, &tHlts.b_L1_NotBsc2_BscMinBiasOR_Prescl);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasThreshold1")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasThreshold1", &tHlts.L1_NotBsc2_BscMinBiasThreshold1, &tHlts.b_L1_NotBsc2_BscMinBiasThreshold1);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasThreshold1_Prescl")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasThreshold1_Prescl", &tHlts.L1_NotBsc2_BscMinBiasThreshold1_Prescl, &tHlts.b_L1_NotBsc2_BscMinBiasThreshold1_Prescl);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasThreshold2")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasThreshold2", &tHlts.L1_NotBsc2_BscMinBiasThreshold2, &tHlts.b_L1_NotBsc2_BscMinBiasThreshold2);
   if (t->GetBranch("L1_NotBsc2_BscMinBiasThreshold2_Prescl")) t->SetBranchAddress("L1_NotBsc2_BscMinBiasThreshold2_Prescl", &tHlts.L1_NotBsc2_BscMinBiasThreshold2_Prescl, &tHlts.b_L1_NotBsc2_BscMinBiasThreshold2_Prescl);
   if (t->GetBranch("L1_NotBsc2_HcalHfCoincidencePm")) t->SetBranchAddress("L1_NotBsc2_HcalHfCoincidencePm", &tHlts.L1_NotBsc2_HcalHfCoincidencePm, &tHlts.b_L1_NotBsc2_HcalHfCoincidencePm);
   if (t->GetBranch("L1_NotBsc2_HcalHfCoincidencePm_Prescl")) t->SetBranchAddress("L1_NotBsc2_HcalHfCoincidencePm_Prescl", &tHlts.L1_NotBsc2_HcalHfCoincidencePm_Prescl, &tHlts.b_L1_NotBsc2_HcalHfCoincidencePm_Prescl);
   if (t->GetBranch("L1_NotBsc2_HcalHfMmOrPpOrPm")) t->SetBranchAddress("L1_NotBsc2_HcalHfMmOrPpOrPm", &tHlts.L1_NotBsc2_HcalHfMmOrPpOrPm, &tHlts.b_L1_NotBsc2_HcalHfMmOrPpOrPm);
   if (t->GetBranch("L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl")) t->SetBranchAddress("L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl", &tHlts.L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl, &tHlts.b_L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl);
   if (t->GetBranch("L1_NotBsc2_HcalHfMmpOrMpp")) t->SetBranchAddress("L1_NotBsc2_HcalHfMmpOrMpp", &tHlts.L1_NotBsc2_HcalHfMmpOrMpp, &tHlts.b_L1_NotBsc2_HcalHfMmpOrMpp);
   if (t->GetBranch("L1_NotBsc2_HcalHfMmpOrMpp_Prescl")) t->SetBranchAddress("L1_NotBsc2_HcalHfMmpOrMpp_Prescl", &tHlts.L1_NotBsc2_HcalHfMmpOrMpp_Prescl, &tHlts.b_L1_NotBsc2_HcalHfMmpOrMpp_Prescl);
   if (t->GetBranch("L1_SingleEG12")) t->SetBranchAddress("L1_SingleEG12", &tHlts.L1_SingleEG12, &tHlts.b_L1_SingleEG12);
   if (t->GetBranch("L1_SingleEG12_Prescl")) t->SetBranchAddress("L1_SingleEG12_Prescl", &tHlts.L1_SingleEG12_Prescl, &tHlts.b_L1_SingleEG12_Prescl);
   if (t->GetBranch("L1_SingleEG15")) t->SetBranchAddress("L1_SingleEG15", &tHlts.L1_SingleEG15, &tHlts.b_L1_SingleEG15);
   if (t->GetBranch("L1_SingleEG15_Prescl")) t->SetBranchAddress("L1_SingleEG15_Prescl", &tHlts.L1_SingleEG15_Prescl, &tHlts.b_L1_SingleEG15_Prescl);
   if (t->GetBranch("L1_SingleEG18")) t->SetBranchAddress("L1_SingleEG18", &tHlts.L1_SingleEG18, &tHlts.b_L1_SingleEG18);
   if (t->GetBranch("L1_SingleEG18_Prescl")) t->SetBranchAddress("L1_SingleEG18_Prescl", &tHlts.L1_SingleEG18_Prescl, &tHlts.b_L1_SingleEG18_Prescl);
   if (t->GetBranch("L1_SingleEG20")) t->SetBranchAddress("L1_SingleEG20", &tHlts.L1_SingleEG20, &tHlts.b_L1_SingleEG20);
   if (t->GetBranch("L1_SingleEG20_Prescl")) t->SetBranchAddress("L1_SingleEG20_Prescl", &tHlts.L1_SingleEG20_Prescl, &tHlts.b_L1_SingleEG20_Prescl);
   if (t->GetBranch("L1_SingleEG22")) t->SetBranchAddress("L1_SingleEG22", &tHlts.L1_SingleEG22, &tHlts.b_L1_SingleEG22);
   if (t->GetBranch("L1_SingleEG22_Prescl")) t->SetBranchAddress("L1_SingleEG22_Prescl", &tHlts.L1_SingleEG22_Prescl, &tHlts.b_L1_SingleEG22_Prescl);
   if (t->GetBranch("L1_SingleEG2_BptxAND")) t->SetBranchAddress("L1_SingleEG2_BptxAND", &tHlts.L1_SingleEG2_BptxAND, &tHlts.b_L1_SingleEG2_BptxAND);
   if (t->GetBranch("L1_SingleEG2_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleEG2_BptxAND_Prescl", &tHlts.L1_SingleEG2_BptxAND_Prescl, &tHlts.b_L1_SingleEG2_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleEG30")) t->SetBranchAddress("L1_SingleEG30", &tHlts.L1_SingleEG30, &tHlts.b_L1_SingleEG30);
   if (t->GetBranch("L1_SingleEG30_Prescl")) t->SetBranchAddress("L1_SingleEG30_Prescl", &tHlts.L1_SingleEG30_Prescl, &tHlts.b_L1_SingleEG30_Prescl);
   if (t->GetBranch("L1_SingleEG5_BptxAND")) t->SetBranchAddress("L1_SingleEG5_BptxAND", &tHlts.L1_SingleEG5_BptxAND, &tHlts.b_L1_SingleEG5_BptxAND);
   if (t->GetBranch("L1_SingleEG5_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleEG5_BptxAND_Prescl", &tHlts.L1_SingleEG5_BptxAND_Prescl, &tHlts.b_L1_SingleEG5_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleEG8_BptxAND")) t->SetBranchAddress("L1_SingleEG8_BptxAND", &tHlts.L1_SingleEG8_BptxAND, &tHlts.b_L1_SingleEG8_BptxAND);
   if (t->GetBranch("L1_SingleEG8_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleEG8_BptxAND_Prescl", &tHlts.L1_SingleEG8_BptxAND_Prescl, &tHlts.b_L1_SingleEG8_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet128")) t->SetBranchAddress("L1_SingleJet128", &tHlts.L1_SingleJet128, &tHlts.b_L1_SingleJet128);
   if (t->GetBranch("L1_SingleJet128_Prescl")) t->SetBranchAddress("L1_SingleJet128_Prescl", &tHlts.L1_SingleJet128_Prescl, &tHlts.b_L1_SingleJet128_Prescl);
   if (t->GetBranch("L1_SingleJet128_BptxAND")) t->SetBranchAddress("L1_SingleJet128_BptxAND", &tHlts.L1_SingleJet128_BptxAND, &tHlts.b_L1_SingleJet128_BptxAND);
   if (t->GetBranch("L1_SingleJet128_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet128_BptxAND_Prescl", &tHlts.L1_SingleJet128_BptxAND_Prescl, &tHlts.b_L1_SingleJet128_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet128_Central_BptxAND")) t->SetBranchAddress("L1_SingleJet128_Central_BptxAND", &tHlts.L1_SingleJet128_Central_BptxAND, &tHlts.b_L1_SingleJet128_Central_BptxAND);
   if (t->GetBranch("L1_SingleJet128_Central_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet128_Central_BptxAND_Prescl", &tHlts.L1_SingleJet128_Central_BptxAND_Prescl, &tHlts.b_L1_SingleJet128_Central_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet16")) t->SetBranchAddress("L1_SingleJet16", &tHlts.L1_SingleJet16, &tHlts.b_L1_SingleJet16);
   if (t->GetBranch("L1_SingleJet16_Prescl")) t->SetBranchAddress("L1_SingleJet16_Prescl", &tHlts.L1_SingleJet16_Prescl, &tHlts.b_L1_SingleJet16_Prescl);
   if (t->GetBranch("L1_SingleJet16_Central_NotETT140_BptxAND")) t->SetBranchAddress("L1_SingleJet16_Central_NotETT140_BptxAND", &tHlts.L1_SingleJet16_Central_NotETT140_BptxAND, &tHlts.b_L1_SingleJet16_Central_NotETT140_BptxAND);
   if (t->GetBranch("L1_SingleJet16_Central_NotETT140_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet16_Central_NotETT140_BptxAND_Prescl", &tHlts.L1_SingleJet16_Central_NotETT140_BptxAND_Prescl, &tHlts.b_L1_SingleJet16_Central_NotETT140_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo")) t->SetBranchAddress("L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo", &tHlts.L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo, &tHlts.b_L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo);
   if (t->GetBranch("L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl")) t->SetBranchAddress("L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl", &tHlts.L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl, &tHlts.b_L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl);
   if (t->GetBranch("L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo")) t->SetBranchAddress("L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo", &tHlts.L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo, &tHlts.b_L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo);
   if (t->GetBranch("L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl")) t->SetBranchAddress("L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl", &tHlts.L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl, &tHlts.b_L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl);
   if (t->GetBranch("L1_SingleJet36_BptxAND")) t->SetBranchAddress("L1_SingleJet36_BptxAND", &tHlts.L1_SingleJet36_BptxAND, &tHlts.b_L1_SingleJet36_BptxAND);
   if (t->GetBranch("L1_SingleJet36_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl", &tHlts.L1_SingleJet36_BptxAND_Prescl, &tHlts.b_L1_SingleJet36_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet52_BptxAND")) t->SetBranchAddress("L1_SingleJet52_BptxAND", &tHlts.L1_SingleJet52_BptxAND, &tHlts.b_L1_SingleJet52_BptxAND);
   if (t->GetBranch("L1_SingleJet52_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl", &tHlts.L1_SingleJet52_BptxAND_Prescl, &tHlts.b_L1_SingleJet52_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet52_Central_BptxAND")) t->SetBranchAddress("L1_SingleJet52_Central_BptxAND", &tHlts.L1_SingleJet52_Central_BptxAND, &tHlts.b_L1_SingleJet52_Central_BptxAND);
   if (t->GetBranch("L1_SingleJet52_Central_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet52_Central_BptxAND_Prescl", &tHlts.L1_SingleJet52_Central_BptxAND_Prescl, &tHlts.b_L1_SingleJet52_Central_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet68_BptxAND")) t->SetBranchAddress("L1_SingleJet68_BptxAND", &tHlts.L1_SingleJet68_BptxAND, &tHlts.b_L1_SingleJet68_BptxAND);
   if (t->GetBranch("L1_SingleJet68_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet68_BptxAND_Prescl", &tHlts.L1_SingleJet68_BptxAND_Prescl, &tHlts.b_L1_SingleJet68_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet80_BptxAND")) t->SetBranchAddress("L1_SingleJet80_BptxAND", &tHlts.L1_SingleJet80_BptxAND, &tHlts.b_L1_SingleJet80_BptxAND);
   if (t->GetBranch("L1_SingleJet80_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet80_BptxAND_Prescl", &tHlts.L1_SingleJet80_BptxAND_Prescl, &tHlts.b_L1_SingleJet80_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet80_Central")) t->SetBranchAddress("L1_SingleJet80_Central", &tHlts.L1_SingleJet80_Central, &tHlts.b_L1_SingleJet80_Central);
   if (t->GetBranch("L1_SingleJet80_Central_Prescl")) t->SetBranchAddress("L1_SingleJet80_Central_Prescl", &tHlts.L1_SingleJet80_Central_Prescl, &tHlts.b_L1_SingleJet80_Central_Prescl);
   if (t->GetBranch("L1_SingleJet80_Central_BptxAND")) t->SetBranchAddress("L1_SingleJet80_Central_BptxAND", &tHlts.L1_SingleJet80_Central_BptxAND, &tHlts.b_L1_SingleJet80_Central_BptxAND);
   if (t->GetBranch("L1_SingleJet80_Central_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet80_Central_BptxAND_Prescl", &tHlts.L1_SingleJet80_Central_BptxAND_Prescl, &tHlts.b_L1_SingleJet80_Central_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet92")) t->SetBranchAddress("L1_SingleJet92", &tHlts.L1_SingleJet92, &tHlts.b_L1_SingleJet92);
   if (t->GetBranch("L1_SingleJet92_Prescl")) t->SetBranchAddress("L1_SingleJet92_Prescl", &tHlts.L1_SingleJet92_Prescl, &tHlts.b_L1_SingleJet92_Prescl);
   if (t->GetBranch("L1_SingleJet92_BptxAND")) t->SetBranchAddress("L1_SingleJet92_BptxAND", &tHlts.L1_SingleJet92_BptxAND, &tHlts.b_L1_SingleJet92_BptxAND);
   if (t->GetBranch("L1_SingleJet92_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet92_BptxAND_Prescl", &tHlts.L1_SingleJet92_BptxAND_Prescl, &tHlts.b_L1_SingleJet92_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleJet92_Central")) t->SetBranchAddress("L1_SingleJet92_Central", &tHlts.L1_SingleJet92_Central, &tHlts.b_L1_SingleJet92_Central);
   if (t->GetBranch("L1_SingleJet92_Central_Prescl")) t->SetBranchAddress("L1_SingleJet92_Central_Prescl", &tHlts.L1_SingleJet92_Central_Prescl, &tHlts.b_L1_SingleJet92_Central_Prescl);
   if (t->GetBranch("L1_SingleJet92_Central_BptxAND")) t->SetBranchAddress("L1_SingleJet92_Central_BptxAND", &tHlts.L1_SingleJet92_Central_BptxAND, &tHlts.b_L1_SingleJet92_Central_BptxAND);
   if (t->GetBranch("L1_SingleJet92_Central_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleJet92_Central_BptxAND_Prescl", &tHlts.L1_SingleJet92_Central_BptxAND_Prescl, &tHlts.b_L1_SingleJet92_Central_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleMu0_BptxAND")) t->SetBranchAddress("L1_SingleMu0_BptxAND", &tHlts.L1_SingleMu0_BptxAND, &tHlts.b_L1_SingleMu0_BptxAND);
   if (t->GetBranch("L1_SingleMu0_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleMu0_BptxAND_Prescl", &tHlts.L1_SingleMu0_BptxAND_Prescl, &tHlts.b_L1_SingleMu0_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleMu10")) t->SetBranchAddress("L1_SingleMu10", &tHlts.L1_SingleMu10, &tHlts.b_L1_SingleMu10);
   if (t->GetBranch("L1_SingleMu10_Prescl")) t->SetBranchAddress("L1_SingleMu10_Prescl", &tHlts.L1_SingleMu10_Prescl, &tHlts.b_L1_SingleMu10_Prescl);
   if (t->GetBranch("L1_SingleMu14")) t->SetBranchAddress("L1_SingleMu14", &tHlts.L1_SingleMu14, &tHlts.b_L1_SingleMu14);
   if (t->GetBranch("L1_SingleMu14_Prescl")) t->SetBranchAddress("L1_SingleMu14_Prescl", &tHlts.L1_SingleMu14_Prescl, &tHlts.b_L1_SingleMu14_Prescl);
   if (t->GetBranch("L1_SingleMu20")) t->SetBranchAddress("L1_SingleMu20", &tHlts.L1_SingleMu20, &tHlts.b_L1_SingleMu20);
   if (t->GetBranch("L1_SingleMu20_Prescl")) t->SetBranchAddress("L1_SingleMu20_Prescl", &tHlts.L1_SingleMu20_Prescl, &tHlts.b_L1_SingleMu20_Prescl);
   if (t->GetBranch("L1_SingleMu3_BptxAND")) t->SetBranchAddress("L1_SingleMu3_BptxAND", &tHlts.L1_SingleMu3_BptxAND, &tHlts.b_L1_SingleMu3_BptxAND);
   if (t->GetBranch("L1_SingleMu3_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleMu3_BptxAND_Prescl", &tHlts.L1_SingleMu3_BptxAND_Prescl, &tHlts.b_L1_SingleMu3_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleMu5_BptxAND")) t->SetBranchAddress("L1_SingleMu5_BptxAND", &tHlts.L1_SingleMu5_BptxAND, &tHlts.b_L1_SingleMu5_BptxAND);
   if (t->GetBranch("L1_SingleMu5_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleMu5_BptxAND_Prescl", &tHlts.L1_SingleMu5_BptxAND_Prescl, &tHlts.b_L1_SingleMu5_BptxAND_Prescl);
   if (t->GetBranch("L1_SingleMu7")) t->SetBranchAddress("L1_SingleMu7", &tHlts.L1_SingleMu7, &tHlts.b_L1_SingleMu7);
   if (t->GetBranch("L1_SingleMu7_Prescl")) t->SetBranchAddress("L1_SingleMu7_Prescl", &tHlts.L1_SingleMu7_Prescl, &tHlts.b_L1_SingleMu7_Prescl);
   if (t->GetBranch("L1_SingleMuBeamHalo")) t->SetBranchAddress("L1_SingleMuBeamHalo", &tHlts.L1_SingleMuBeamHalo, &tHlts.b_L1_SingleMuBeamHalo);
   if (t->GetBranch("L1_SingleMuBeamHalo_Prescl")) t->SetBranchAddress("L1_SingleMuBeamHalo_Prescl", &tHlts.L1_SingleMuBeamHalo_Prescl, &tHlts.b_L1_SingleMuBeamHalo_Prescl);
   if (t->GetBranch("L1_SingleMuOpen_BptxAND")) t->SetBranchAddress("L1_SingleMuOpen_BptxAND", &tHlts.L1_SingleMuOpen_BptxAND, &tHlts.b_L1_SingleMuOpen_BptxAND);
   if (t->GetBranch("L1_SingleMuOpen_BptxAND_Prescl")) t->SetBranchAddress("L1_SingleMuOpen_BptxAND_Prescl", &tHlts.L1_SingleMuOpen_BptxAND_Prescl, &tHlts.b_L1_SingleMuOpen_BptxAND_Prescl);
   if (t->GetBranch("L1_ZdcCaloMinus")) t->SetBranchAddress("L1_ZdcCaloMinus", &tHlts.L1_ZdcCaloMinus, &tHlts.b_L1_ZdcCaloMinus);
   if (t->GetBranch("L1_ZdcCaloMinus_Prescl")) t->SetBranchAddress("L1_ZdcCaloMinus_Prescl", &tHlts.L1_ZdcCaloMinus_Prescl, &tHlts.b_L1_ZdcCaloMinus_Prescl);
   if (t->GetBranch("L1_ZdcCaloMinus_BptxAND")) t->SetBranchAddress("L1_ZdcCaloMinus_BptxAND", &tHlts.L1_ZdcCaloMinus_BptxAND, &tHlts.b_L1_ZdcCaloMinus_BptxAND);
   if (t->GetBranch("L1_ZdcCaloMinus_BptxAND_Prescl")) t->SetBranchAddress("L1_ZdcCaloMinus_BptxAND_Prescl", &tHlts.L1_ZdcCaloMinus_BptxAND_Prescl, &tHlts.b_L1_ZdcCaloMinus_BptxAND_Prescl);
   if (t->GetBranch("L1_ZdcCaloPlus")) t->SetBranchAddress("L1_ZdcCaloPlus", &tHlts.L1_ZdcCaloPlus, &tHlts.b_L1_ZdcCaloPlus);
   if (t->GetBranch("L1_ZdcCaloPlus_Prescl")) t->SetBranchAddress("L1_ZdcCaloPlus_Prescl", &tHlts.L1_ZdcCaloPlus_Prescl, &tHlts.b_L1_ZdcCaloPlus_Prescl);
   if (t->GetBranch("L1_ZdcCaloPlus_BptxAND")) t->SetBranchAddress("L1_ZdcCaloPlus_BptxAND", &tHlts.L1_ZdcCaloPlus_BptxAND, &tHlts.b_L1_ZdcCaloPlus_BptxAND);
   if (t->GetBranch("L1_ZdcCaloPlus_BptxAND_Prescl")) t->SetBranchAddress("L1_ZdcCaloPlus_BptxAND_Prescl", &tHlts.L1_ZdcCaloPlus_BptxAND_Prescl, &tHlts.b_L1_ZdcCaloPlus_BptxAND_Prescl);
   if (t->GetBranch("L1_ZdcCaloPlus_ZdcCaloMinus")) t->SetBranchAddress("L1_ZdcCaloPlus_ZdcCaloMinus", &tHlts.L1_ZdcCaloPlus_ZdcCaloMinus, &tHlts.b_L1_ZdcCaloPlus_ZdcCaloMinus);
   if (t->GetBranch("L1_ZdcCaloPlus_ZdcCaloMinus_Prescl")) t->SetBranchAddress("L1_ZdcCaloPlus_ZdcCaloMinus_Prescl", &tHlts.L1_ZdcCaloPlus_ZdcCaloMinus_Prescl, &tHlts.b_L1_ZdcCaloPlus_ZdcCaloMinus_Prescl);
   if (t->GetBranch("L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND")) t->SetBranchAddress("L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND", &tHlts.L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND, &tHlts.b_L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND);
   if (t->GetBranch("L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl")) t->SetBranchAddress("L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl", &tHlts.L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl, &tHlts.b_L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl);
   if (t->GetBranch("L1_ZeroBias")) t->SetBranchAddress("L1_ZeroBias", &tHlts.L1_ZeroBias, &tHlts.b_L1_ZeroBias);
   if (t->GetBranch("L1_ZeroBias_Prescl")) t->SetBranchAddress("L1_ZeroBias_Prescl", &tHlts.L1_ZeroBias_Prescl, &tHlts.b_L1_ZeroBias_Prescl);
   if (t->GetBranch("L1Tech_BPTX_PreBPTX.v0")) t->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0", &tHlts.L1Tech_BPTX_PreBPTX_v0, &tHlts.b_L1Tech_BPTX_PreBPTX_v0);
   if (t->GetBranch("L1Tech_BPTX_PreBPTX.v0_Prescl")) t->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0_Prescl", &tHlts.L1Tech_BPTX_PreBPTX_v0_Prescl, &tHlts.b_L1Tech_BPTX_PreBPTX_v0_Prescl);
   if (t->GetBranch("L1Tech_BPTX_minus.v0")) t->SetBranchAddress("L1Tech_BPTX_minus.v0", &tHlts.L1Tech_BPTX_minus_v0, &tHlts.b_L1Tech_BPTX_minus_v0);
   if (t->GetBranch("L1Tech_BPTX_minus.v0_Prescl")) t->SetBranchAddress("L1Tech_BPTX_minus.v0_Prescl", &tHlts.L1Tech_BPTX_minus_v0_Prescl, &tHlts.b_L1Tech_BPTX_minus_v0_Prescl);
   if (t->GetBranch("L1Tech_BPTX_minus_AND_not_plus.v0")) t->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0", &tHlts.L1Tech_BPTX_minus_AND_not_plus_v0, &tHlts.b_L1Tech_BPTX_minus_AND_not_plus_v0);
   if (t->GetBranch("L1Tech_BPTX_minus_AND_not_plus.v0_Prescl")) t->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_Prescl", &tHlts.L1Tech_BPTX_minus_AND_not_plus_v0_Prescl, &tHlts.b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl);
   if (t->GetBranch("L1Tech_BPTX_plus.v0")) t->SetBranchAddress("L1Tech_BPTX_plus.v0", &tHlts.L1Tech_BPTX_plus_v0, &tHlts.b_L1Tech_BPTX_plus_v0);
   if (t->GetBranch("L1Tech_BPTX_plus.v0_Prescl")) t->SetBranchAddress("L1Tech_BPTX_plus.v0_Prescl", &tHlts.L1Tech_BPTX_plus_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_v0_Prescl);
   if (t->GetBranch("L1Tech_BPTX_plus_AND_NOT_minus.v0")) t->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0", &tHlts.L1Tech_BPTX_plus_AND_NOT_minus_v0, &tHlts.b_L1Tech_BPTX_plus_AND_NOT_minus_v0);
   if (t->GetBranch("L1Tech_BPTX_plus_AND_NOT_minus.v0_Prescl")) t->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_Prescl", &tHlts.L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl);
   if (t->GetBranch("L1Tech_BPTX_plus_AND_minus.v0")) t->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0", &tHlts.L1Tech_BPTX_plus_AND_minus_v0, &tHlts.b_L1Tech_BPTX_plus_AND_minus_v0);
   if (t->GetBranch("L1Tech_BPTX_plus_AND_minus.v0_Prescl")) t->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_Prescl", &tHlts.L1Tech_BPTX_plus_AND_minus_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_AND_minus_v0_Prescl);
   if (t->GetBranch("L1Tech_BPTX_plus_AND_minus_instance1.v0")) t->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0", &tHlts.L1Tech_BPTX_plus_AND_minus_instance1_v0, &tHlts.b_L1Tech_BPTX_plus_AND_minus_instance1_v0);
   if (t->GetBranch("L1Tech_BPTX_plus_AND_minus_instance1.v0_Prescl")) t->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0_Prescl", &tHlts.L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl);
   if (t->GetBranch("L1Tech_BPTX_plus_OR_minus.v0")) t->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0", &tHlts.L1Tech_BPTX_plus_OR_minus_v0, &tHlts.b_L1Tech_BPTX_plus_OR_minus_v0);
   if (t->GetBranch("L1Tech_BPTX_plus_OR_minus.v0_Prescl")) t->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0_Prescl", &tHlts.L1Tech_BPTX_plus_OR_minus_v0_Prescl, &tHlts.b_L1Tech_BPTX_plus_OR_minus_v0_Prescl);
   if (t->GetBranch("L1Tech_BPTX_quiet.v0")) t->SetBranchAddress("L1Tech_BPTX_quiet.v0", &tHlts.L1Tech_BPTX_quiet_v0, &tHlts.b_L1Tech_BPTX_quiet_v0);
   if (t->GetBranch("L1Tech_BPTX_quiet.v0_Prescl")) t->SetBranchAddress("L1Tech_BPTX_quiet.v0_Prescl", &tHlts.L1Tech_BPTX_quiet_v0_Prescl, &tHlts.b_L1Tech_BPTX_quiet_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_HighMultiplicity.v0")) t->SetBranchAddress("L1Tech_BSC_HighMultiplicity.v0", &tHlts.L1Tech_BSC_HighMultiplicity_v0, &tHlts.b_L1Tech_BSC_HighMultiplicity_v0);
   if (t->GetBranch("L1Tech_BSC_HighMultiplicity.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_HighMultiplicity.v0_Prescl", &tHlts.L1Tech_BSC_HighMultiplicity_v0_Prescl, &tHlts.b_L1Tech_BSC_HighMultiplicity_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_halo_beam1_inner.v0")) t->SetBranchAddress("L1Tech_BSC_halo_beam1_inner.v0", &tHlts.L1Tech_BSC_halo_beam1_inner_v0, &tHlts.b_L1Tech_BSC_halo_beam1_inner_v0);
   if (t->GetBranch("L1Tech_BSC_halo_beam1_inner.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_halo_beam1_inner.v0_Prescl", &tHlts.L1Tech_BSC_halo_beam1_inner_v0_Prescl, &tHlts.b_L1Tech_BSC_halo_beam1_inner_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_halo_beam1_outer.v0")) t->SetBranchAddress("L1Tech_BSC_halo_beam1_outer.v0", &tHlts.L1Tech_BSC_halo_beam1_outer_v0, &tHlts.b_L1Tech_BSC_halo_beam1_outer_v0);
   if (t->GetBranch("L1Tech_BSC_halo_beam1_outer.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_halo_beam1_outer.v0_Prescl", &tHlts.L1Tech_BSC_halo_beam1_outer_v0_Prescl, &tHlts.b_L1Tech_BSC_halo_beam1_outer_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_halo_beam2_inner.v0")) t->SetBranchAddress("L1Tech_BSC_halo_beam2_inner.v0", &tHlts.L1Tech_BSC_halo_beam2_inner_v0, &tHlts.b_L1Tech_BSC_halo_beam2_inner_v0);
   if (t->GetBranch("L1Tech_BSC_halo_beam2_inner.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_halo_beam2_inner.v0_Prescl", &tHlts.L1Tech_BSC_halo_beam2_inner_v0_Prescl, &tHlts.b_L1Tech_BSC_halo_beam2_inner_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_halo_beam2_outer.v0")) t->SetBranchAddress("L1Tech_BSC_halo_beam2_outer.v0", &tHlts.L1Tech_BSC_halo_beam2_outer_v0, &tHlts.b_L1Tech_BSC_halo_beam2_outer_v0);
   if (t->GetBranch("L1Tech_BSC_halo_beam2_outer.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_halo_beam2_outer.v0_Prescl", &tHlts.L1Tech_BSC_halo_beam2_outer_v0_Prescl, &tHlts.b_L1Tech_BSC_halo_beam2_outer_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_minBias_OR.v0")) t->SetBranchAddress("L1Tech_BSC_minBias_OR.v0", &tHlts.L1Tech_BSC_minBias_OR_v0, &tHlts.b_L1Tech_BSC_minBias_OR_v0);
   if (t->GetBranch("L1Tech_BSC_minBias_OR.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_minBias_OR.v0_Prescl", &tHlts.L1Tech_BSC_minBias_OR_v0_Prescl, &tHlts.b_L1Tech_BSC_minBias_OR_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_minBias_inner_threshold1.v0")) t->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold1.v0", &tHlts.L1Tech_BSC_minBias_inner_threshold1_v0, &tHlts.b_L1Tech_BSC_minBias_inner_threshold1_v0);
   if (t->GetBranch("L1Tech_BSC_minBias_inner_threshold1.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold1.v0_Prescl", &tHlts.L1Tech_BSC_minBias_inner_threshold1_v0_Prescl, &tHlts.b_L1Tech_BSC_minBias_inner_threshold1_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_minBias_inner_threshold2.v0")) t->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold2.v0", &tHlts.L1Tech_BSC_minBias_inner_threshold2_v0, &tHlts.b_L1Tech_BSC_minBias_inner_threshold2_v0);
   if (t->GetBranch("L1Tech_BSC_minBias_inner_threshold2.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold2.v0_Prescl", &tHlts.L1Tech_BSC_minBias_inner_threshold2_v0_Prescl, &tHlts.b_L1Tech_BSC_minBias_inner_threshold2_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_minBias_threshold1.v0")) t->SetBranchAddress("L1Tech_BSC_minBias_threshold1.v0", &tHlts.L1Tech_BSC_minBias_threshold1_v0, &tHlts.b_L1Tech_BSC_minBias_threshold1_v0);
   if (t->GetBranch("L1Tech_BSC_minBias_threshold1.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_minBias_threshold1.v0_Prescl", &tHlts.L1Tech_BSC_minBias_threshold1_v0_Prescl, &tHlts.b_L1Tech_BSC_minBias_threshold1_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_minBias_threshold2.v0")) t->SetBranchAddress("L1Tech_BSC_minBias_threshold2.v0", &tHlts.L1Tech_BSC_minBias_threshold2_v0, &tHlts.b_L1Tech_BSC_minBias_threshold2_v0);
   if (t->GetBranch("L1Tech_BSC_minBias_threshold2.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_minBias_threshold2.v0_Prescl", &tHlts.L1Tech_BSC_minBias_threshold2_v0_Prescl, &tHlts.b_L1Tech_BSC_minBias_threshold2_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_splash_beam1.v0")) t->SetBranchAddress("L1Tech_BSC_splash_beam1.v0", &tHlts.L1Tech_BSC_splash_beam1_v0, &tHlts.b_L1Tech_BSC_splash_beam1_v0);
   if (t->GetBranch("L1Tech_BSC_splash_beam1.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_splash_beam1.v0_Prescl", &tHlts.L1Tech_BSC_splash_beam1_v0_Prescl, &tHlts.b_L1Tech_BSC_splash_beam1_v0_Prescl);
   if (t->GetBranch("L1Tech_BSC_splash_beam2.v0")) t->SetBranchAddress("L1Tech_BSC_splash_beam2.v0", &tHlts.L1Tech_BSC_splash_beam2_v0, &tHlts.b_L1Tech_BSC_splash_beam2_v0);
   if (t->GetBranch("L1Tech_BSC_splash_beam2.v0_Prescl")) t->SetBranchAddress("L1Tech_BSC_splash_beam2.v0_Prescl", &tHlts.L1Tech_BSC_splash_beam2_v0_Prescl, &tHlts.b_L1Tech_BSC_splash_beam2_v0_Prescl);
   if (t->GetBranch("L1Tech_CASTOR_HaloMuon.v0")) t->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0", &tHlts.L1Tech_CASTOR_HaloMuon_v0, &tHlts.b_L1Tech_CASTOR_HaloMuon_v0);
   if (t->GetBranch("L1Tech_CASTOR_HaloMuon.v0_Prescl")) t->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0_Prescl", &tHlts.L1Tech_CASTOR_HaloMuon_v0_Prescl, &tHlts.b_L1Tech_CASTOR_HaloMuon_v0_Prescl);
   if (t->GetBranch("L1Tech_DT_GlobalOR.v0")) t->SetBranchAddress("L1Tech_DT_GlobalOR.v0", &tHlts.L1Tech_DT_GlobalOR_v0, &tHlts.b_L1Tech_DT_GlobalOR_v0);
   if (t->GetBranch("L1Tech_DT_GlobalOR.v0_Prescl")) t->SetBranchAddress("L1Tech_DT_GlobalOR.v0_Prescl", &tHlts.L1Tech_DT_GlobalOR_v0_Prescl, &tHlts.b_L1Tech_DT_GlobalOR_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St1Sect45_down.v0")) t->SetBranchAddress("L1Tech_FSC_St1Sect45_down.v0", &tHlts.L1Tech_FSC_St1Sect45_down_v0, &tHlts.b_L1Tech_FSC_St1Sect45_down_v0);
   if (t->GetBranch("L1Tech_FSC_St1Sect45_down.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St1Sect45_down.v0_Prescl", &tHlts.L1Tech_FSC_St1Sect45_down_v0_Prescl, &tHlts.b_L1Tech_FSC_St1Sect45_down_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St1Sect45_upp.v0")) t->SetBranchAddress("L1Tech_FSC_St1Sect45_upp.v0", &tHlts.L1Tech_FSC_St1Sect45_upp_v0, &tHlts.b_L1Tech_FSC_St1Sect45_upp_v0);
   if (t->GetBranch("L1Tech_FSC_St1Sect45_upp.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St1Sect45_upp.v0_Prescl", &tHlts.L1Tech_FSC_St1Sect45_upp_v0_Prescl, &tHlts.b_L1Tech_FSC_St1Sect45_upp_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St1Sect56_down.v0")) t->SetBranchAddress("L1Tech_FSC_St1Sect56_down.v0", &tHlts.L1Tech_FSC_St1Sect56_down_v0, &tHlts.b_L1Tech_FSC_St1Sect56_down_v0);
   if (t->GetBranch("L1Tech_FSC_St1Sect56_down.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St1Sect56_down.v0_Prescl", &tHlts.L1Tech_FSC_St1Sect56_down_v0_Prescl, &tHlts.b_L1Tech_FSC_St1Sect56_down_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St1Sect56_upp.v0")) t->SetBranchAddress("L1Tech_FSC_St1Sect56_upp.v0", &tHlts.L1Tech_FSC_St1Sect56_upp_v0, &tHlts.b_L1Tech_FSC_St1Sect56_upp_v0);
   if (t->GetBranch("L1Tech_FSC_St1Sect56_upp.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St1Sect56_upp.v0_Prescl", &tHlts.L1Tech_FSC_St1Sect56_upp_v0_Prescl, &tHlts.b_L1Tech_FSC_St1Sect56_upp_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St2Sect45_down.v0")) t->SetBranchAddress("L1Tech_FSC_St2Sect45_down.v0", &tHlts.L1Tech_FSC_St2Sect45_down_v0, &tHlts.b_L1Tech_FSC_St2Sect45_down_v0);
   if (t->GetBranch("L1Tech_FSC_St2Sect45_down.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St2Sect45_down.v0_Prescl", &tHlts.L1Tech_FSC_St2Sect45_down_v0_Prescl, &tHlts.b_L1Tech_FSC_St2Sect45_down_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St2Sect45_upp.v0")) t->SetBranchAddress("L1Tech_FSC_St2Sect45_upp.v0", &tHlts.L1Tech_FSC_St2Sect45_upp_v0, &tHlts.b_L1Tech_FSC_St2Sect45_upp_v0);
   if (t->GetBranch("L1Tech_FSC_St2Sect45_upp.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St2Sect45_upp.v0_Prescl", &tHlts.L1Tech_FSC_St2Sect45_upp_v0_Prescl, &tHlts.b_L1Tech_FSC_St2Sect45_upp_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St2Sect56_down.v0")) t->SetBranchAddress("L1Tech_FSC_St2Sect56_down.v0", &tHlts.L1Tech_FSC_St2Sect56_down_v0, &tHlts.b_L1Tech_FSC_St2Sect56_down_v0);
   if (t->GetBranch("L1Tech_FSC_St2Sect56_down.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St2Sect56_down.v0_Prescl", &tHlts.L1Tech_FSC_St2Sect56_down_v0_Prescl, &tHlts.b_L1Tech_FSC_St2Sect56_down_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St2Sect56_upp.v0")) t->SetBranchAddress("L1Tech_FSC_St2Sect56_upp.v0", &tHlts.L1Tech_FSC_St2Sect56_upp_v0, &tHlts.b_L1Tech_FSC_St2Sect56_upp_v0);
   if (t->GetBranch("L1Tech_FSC_St2Sect56_upp.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St2Sect56_upp.v0_Prescl", &tHlts.L1Tech_FSC_St2Sect56_upp_v0_Prescl, &tHlts.b_L1Tech_FSC_St2Sect56_upp_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St3Sect45_downLeft.v0")) t->SetBranchAddress("L1Tech_FSC_St3Sect45_downLeft.v0", &tHlts.L1Tech_FSC_St3Sect45_downLeft_v0, &tHlts.b_L1Tech_FSC_St3Sect45_downLeft_v0);
   if (t->GetBranch("L1Tech_FSC_St3Sect45_downLeft.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St3Sect45_downLeft.v0_Prescl", &tHlts.L1Tech_FSC_St3Sect45_downLeft_v0_Prescl, &tHlts.b_L1Tech_FSC_St3Sect45_downLeft_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St3Sect45_downRight.v0")) t->SetBranchAddress("L1Tech_FSC_St3Sect45_downRight.v0", &tHlts.L1Tech_FSC_St3Sect45_downRight_v0, &tHlts.b_L1Tech_FSC_St3Sect45_downRight_v0);
   if (t->GetBranch("L1Tech_FSC_St3Sect45_downRight.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St3Sect45_downRight.v0_Prescl", &tHlts.L1Tech_FSC_St3Sect45_downRight_v0_Prescl, &tHlts.b_L1Tech_FSC_St3Sect45_downRight_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St3Sect45_uppLeft.v0")) t->SetBranchAddress("L1Tech_FSC_St3Sect45_uppLeft.v0", &tHlts.L1Tech_FSC_St3Sect45_uppLeft_v0, &tHlts.b_L1Tech_FSC_St3Sect45_uppLeft_v0);
   if (t->GetBranch("L1Tech_FSC_St3Sect45_uppLeft.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St3Sect45_uppLeft.v0_Prescl", &tHlts.L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl, &tHlts.b_L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St3Sect45_uppRight.v0")) t->SetBranchAddress("L1Tech_FSC_St3Sect45_uppRight.v0", &tHlts.L1Tech_FSC_St3Sect45_uppRight_v0, &tHlts.b_L1Tech_FSC_St3Sect45_uppRight_v0);
   if (t->GetBranch("L1Tech_FSC_St3Sect45_uppRight.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St3Sect45_uppRight.v0_Prescl", &tHlts.L1Tech_FSC_St3Sect45_uppRight_v0_Prescl, &tHlts.b_L1Tech_FSC_St3Sect45_uppRight_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St3Sect56_downLeft.v0")) t->SetBranchAddress("L1Tech_FSC_St3Sect56_downLeft.v0", &tHlts.L1Tech_FSC_St3Sect56_downLeft_v0, &tHlts.b_L1Tech_FSC_St3Sect56_downLeft_v0);
   if (t->GetBranch("L1Tech_FSC_St3Sect56_downLeft.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St3Sect56_downLeft.v0_Prescl", &tHlts.L1Tech_FSC_St3Sect56_downLeft_v0_Prescl, &tHlts.b_L1Tech_FSC_St3Sect56_downLeft_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St3Sect56_downRight.v0")) t->SetBranchAddress("L1Tech_FSC_St3Sect56_downRight.v0", &tHlts.L1Tech_FSC_St3Sect56_downRight_v0, &tHlts.b_L1Tech_FSC_St3Sect56_downRight_v0);
   if (t->GetBranch("L1Tech_FSC_St3Sect56_downRight.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St3Sect56_downRight.v0_Prescl", &tHlts.L1Tech_FSC_St3Sect56_downRight_v0_Prescl, &tHlts.b_L1Tech_FSC_St3Sect56_downRight_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St3Sect56_uppLeft.v0")) t->SetBranchAddress("L1Tech_FSC_St3Sect56_uppLeft.v0", &tHlts.L1Tech_FSC_St3Sect56_uppLeft_v0, &tHlts.b_L1Tech_FSC_St3Sect56_uppLeft_v0);
   if (t->GetBranch("L1Tech_FSC_St3Sect56_uppLeft.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St3Sect56_uppLeft.v0_Prescl", &tHlts.L1Tech_FSC_St3Sect56_uppLeft_v0_Prescl, &tHlts.b_L1Tech_FSC_St3Sect56_uppLeft_v0_Prescl);
   if (t->GetBranch("L1Tech_FSC_St3Sect56_uppRight.v0")) t->SetBranchAddress("L1Tech_FSC_St3Sect56_uppRight.v0", &tHlts.L1Tech_FSC_St3Sect56_uppRight_v0, &tHlts.b_L1Tech_FSC_St3Sect56_uppRight_v0);
   if (t->GetBranch("L1Tech_FSC_St3Sect56_uppRight.v0_Prescl")) t->SetBranchAddress("L1Tech_FSC_St3Sect56_uppRight.v0_Prescl", &tHlts.L1Tech_FSC_St3Sect56_uppRight_v0_Prescl, &tHlts.b_L1Tech_FSC_St3Sect56_uppRight_v0_Prescl);
   if (t->GetBranch("L1Tech_HCAL_HBHE_totalOR.v0")) t->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0", &tHlts.L1Tech_HCAL_HBHE_totalOR_v0, &tHlts.b_L1Tech_HCAL_HBHE_totalOR_v0);
   if (t->GetBranch("L1Tech_HCAL_HBHE_totalOR.v0_Prescl")) t->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0_Prescl", &tHlts.L1Tech_HCAL_HBHE_totalOR_v0_Prescl, &tHlts.b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl);
   if (t->GetBranch("L1Tech_HCAL_HF_MMP_or_MPP.v1")) t->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1", &tHlts.L1Tech_HCAL_HF_MMP_or_MPP_v1, &tHlts.b_L1Tech_HCAL_HF_MMP_or_MPP_v1);
   if (t->GetBranch("L1Tech_HCAL_HF_MMP_or_MPP.v1_Prescl")) t->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1_Prescl", &tHlts.L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl, &tHlts.b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl);
   if (t->GetBranch("L1Tech_HCAL_HF_coincidence_PM.v2")) t->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2", &tHlts.L1Tech_HCAL_HF_coincidence_PM_v2, &tHlts.b_L1Tech_HCAL_HF_coincidence_PM_v2);
   if (t->GetBranch("L1Tech_HCAL_HF_coincidence_PM.v2_Prescl")) t->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2_Prescl", &tHlts.L1Tech_HCAL_HF_coincidence_PM_v2_Prescl, &tHlts.b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl);
   if (t->GetBranch("L1Tech_HCAL_HF_single_channel.v0")) t->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0", &tHlts.L1Tech_HCAL_HF_single_channel_v0, &tHlts.b_L1Tech_HCAL_HF_single_channel_v0);
   if (t->GetBranch("L1Tech_HCAL_HF_single_channel.v0_Prescl")) t->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0_Prescl", &tHlts.L1Tech_HCAL_HF_single_channel_v0_Prescl, &tHlts.b_L1Tech_HCAL_HF_single_channel_v0_Prescl);
   if (t->GetBranch("L1Tech_HCAL_HO_totalOR.v0")) t->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0", &tHlts.L1Tech_HCAL_HO_totalOR_v0, &tHlts.b_L1Tech_HCAL_HO_totalOR_v0);
   if (t->GetBranch("L1Tech_HCAL_HO_totalOR.v0_Prescl")) t->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0_Prescl", &tHlts.L1Tech_HCAL_HO_totalOR_v0_Prescl, &tHlts.b_L1Tech_HCAL_HO_totalOR_v0_Prescl);
   if (t->GetBranch("L1Tech_RPC_TTU_RB0_Cosmics.v0")) t->SetBranchAddress("L1Tech_RPC_TTU_RB0_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_RB0_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_RB0_Cosmics_v0);
   if (t->GetBranch("L1Tech_RPC_TTU_RB0_Cosmics.v0_Prescl")) t->SetBranchAddress("L1Tech_RPC_TTU_RB0_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl);
   if (t->GetBranch("L1Tech_RPC_TTU_RBminus1_Cosmics.v0")) t->SetBranchAddress("L1Tech_RPC_TTU_RBminus1_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_RBminus1_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0);
   if (t->GetBranch("L1Tech_RPC_TTU_RBminus1_Cosmics.v0_Prescl")) t->SetBranchAddress("L1Tech_RPC_TTU_RBminus1_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl);
   if (t->GetBranch("L1Tech_RPC_TTU_RBminus2_Cosmics.v0")) t->SetBranchAddress("L1Tech_RPC_TTU_RBminus2_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_RBminus2_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0);
   if (t->GetBranch("L1Tech_RPC_TTU_RBminus2_Cosmics.v0_Prescl")) t->SetBranchAddress("L1Tech_RPC_TTU_RBminus2_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl);
   if (t->GetBranch("L1Tech_RPC_TTU_RBplus1_Cosmics.v0")) t->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_RBplus1_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0);
   if (t->GetBranch("L1Tech_RPC_TTU_RBplus1_Cosmics.v0_Prescl")) t->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl);
   if (t->GetBranch("L1Tech_RPC_TTU_RBplus2_Cosmics.v0")) t->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_RBplus2_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0);
   if (t->GetBranch("L1Tech_RPC_TTU_RBplus2_Cosmics.v0_Prescl")) t->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl);
   if (t->GetBranch("L1Tech_RPC_TTU_RBst1_collisions.v0")) t->SetBranchAddress("L1Tech_RPC_TTU_RBst1_collisions.v0", &tHlts.L1Tech_RPC_TTU_RBst1_collisions_v0, &tHlts.b_L1Tech_RPC_TTU_RBst1_collisions_v0);
   if (t->GetBranch("L1Tech_RPC_TTU_RBst1_collisions.v0_Prescl")) t->SetBranchAddress("L1Tech_RPC_TTU_RBst1_collisions.v0_Prescl", &tHlts.L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl);
   if (t->GetBranch("L1Tech_RPC_TTU_barrel_Cosmics.v0")) t->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_barrel_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_barrel_Cosmics_v0);
   if (t->GetBranch("L1Tech_RPC_TTU_barrel_Cosmics.v0_Prescl")) t->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl);
   if (t->GetBranch("L1Tech_RPC_TTU_pointing_Cosmics.v0")) t->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0", &tHlts.L1Tech_RPC_TTU_pointing_Cosmics_v0, &tHlts.b_L1Tech_RPC_TTU_pointing_Cosmics_v0);
   if (t->GetBranch("L1Tech_RPC_TTU_pointing_Cosmics.v0_Prescl")) t->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0_Prescl", &tHlts.L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl, &tHlts.b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl);
   if (doCheck) {
      if (t->GetMaximum("NL1IsolEm")>58) { cout <<"FATAL ERROR: Arrary size of NL1IsolEm too small!!!  "<<t->GetMaximum("NL1IsolEm")<<endl; exit(0);
 }      if (t->GetMaximum("NL1NIsolEm")>58) { cout <<"FATAL ERROR: Arrary size of NL1NIsolEm too small!!!  "<<t->GetMaximum("NL1NIsolEm")<<endl; exit(0);
 }      if (t->GetMaximum("NL1Mu")>58) { cout <<"FATAL ERROR: Arrary size of NL1Mu too small!!!  "<<t->GetMaximum("NL1Mu")<<endl; exit(0);
 }      if (t->GetMaximum("NL1CenJet")>58) { cout <<"FATAL ERROR: Arrary size of NL1CenJet too small!!!  "<<t->GetMaximum("NL1CenJet")<<endl; exit(0);
 }      if (t->GetMaximum("NL1ForJet")>58) { cout <<"FATAL ERROR: Arrary size of NL1ForJet too small!!!  "<<t->GetMaximum("NL1ForJet")<<endl; exit(0);
 }      if (t->GetMaximum("NL1Tau")>58) { cout <<"FATAL ERROR: Arrary size of NL1Tau too small!!!  "<<t->GetMaximum("NL1Tau")<<endl; exit(0);
 }   }
}

