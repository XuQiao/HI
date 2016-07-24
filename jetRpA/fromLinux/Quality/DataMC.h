//--------------------------------------------------------for Data/MC comparison----------------------------------------------------------

TString infilename[] = {
"JetTrig/DataMC/FromYX/DATAPPb_MBAllAkPu3PFJet2013HistIterTrkCorrtest1EtaCut10_PA2013_HiForest_PromptReco_KrisztianMB_JSonPPb_forestv84.root",
"TrackTrig/DATAPPb_MBAllAllAkPu3PFJet2013HistIterTrkCorrtest1EtaCut10_PA2013_HiForest_PromptReco_KrisztianMB_JSonPPb_forestv84.root",

"JetTrig/DataMC/Own/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest1EtaCut10_pt15_HP04_prod25_v85_merged_forest_0.root",
"JetTrig/DataMC/Own/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest1EtaCut10_pt15_HP04_hiforest77_hiSignal.root",

"JetTrig/DataMC/Own/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt15_HP04_prod25_v85_merged_forest_0.root",
"JetTrig/DataMC/Own/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt15_HP04_hiforest77_hiSignal.root",

"JetTrig/DataMC/Own/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt30_HP04_prod25_v85_merged_forest_0.root",
"JetTrig/DataMC/Own/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt30_HP04_hiforest77_hiSignal.root",

"JetTrig/DataMC/Own/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt30_HP04_prod25_v85_merged_forest_0.root",
"JetTrig/DataMC/Own/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt50_HP04_hiforest77_hiSignal.root",

"JetTrig/DataMC/Own/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt80_HP04_prod25_v85_merged_forest_0.root",
"JetTrig/DataMC/Own/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt80_HP04_hiforest77_hiSignal.root",

"JetTrig/DataMC/Own/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt120_HP04_prod25_v85_merged_forest_0.root",
"JetTrig/DataMC/Own/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt120_HP04_hiforest77_hiSignal.root",

"JetTrig/DataMC/Own/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt170_HP04_prod25_v85_merged_forest_0.root",
"JetTrig/DataMC/Own/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt170_HP04_hiforest77_hiSignal.root",

"JetTrig/DataMC/Own/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt220_HP04_prod25_v85_merged_forest_0.root",
"JetTrig/DataMC/Own/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest0EtaCut10_pt220_HP04_hiforest77_hiSignal.root",

"TrackTrig/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest1EtaCut10_pt15_HP04_prod25_v85_merged_forest_0.root",
"TrackTrig/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest1EtaCut10_pt15_HP04_hiforest77_hiSignal.root",

"TrackTrig/MCPPb_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest1EtaCut10_pt30_HP04_prod25_v85_merged_forest_0.root",
"TrackTrig/MCPP_AllAkPu3PFJetGenCharge2013HistIterTrkCorrtest1EtaCut10_pt30_HP04_hiforest77_hiSignal.root",

};

TString infileabb[]={
"Jet_DataPPb_trkcut1_MB",
"Trk_DataPPb_trkcut1_MBAll",

"Jet_MCPPb_trkcut1_pthat15",
"Jet_MCPP_trkcut1_pthat15",

"Jet_MCPPb_trkcut0_pthat15",
"Jet_MCPP_trkcut0_pthat15",

"Jet_MCPPb_trkcut0_pthat30",
"Jet_MCPP_trkcut0_pthat30",

"Jet_MCPPb_trkcut0_pthat50",
"Jet_MCPP_trkcut0_pthat50",

"Jet_MCPPb_trkcut0_pthat80",
"Jet_MCPP_trkcut0_pthat80",

"Jet_MCPPb_trkcut0_pthat120",
"Jet_MCPP_trkcut0_pthat120",

"Jet_MCPPb_trkcut0_pthat170",
"Jet_MCPP_trkcut0_pthat170",

"Jet_MCPPb_trkcut0_pthat220",
"Jet_MCPP_trkcut0_pthat220",

"Trk_MCPPb_trkcut1_pthat15",
"Trk_MCPP_trkcut1_pthat15",

"Trk_MCPPb_trkcut1_pthat30",
"Trk_MCPP_trkcut1_pthat30",

};

double pro[][2]={{100,300},{-1.465,0.535}};

TString histoname_[] = {"jetptEta_0-100%","jetptEta_0-100%","jetpt_0-100%"};
TString xtitle_[] = {"#eta","p_{T} [GeV/c]","p_{T} [GeV/c]"};
double xrange_[][2] = {{-1.465,0.535},{100,300},{100,300}};
double dimension_[]={2.5,1.5,1};
int Nbin_[]={13,17,17};
double binbound_[][100]={{-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3},{30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300},{30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300}};
bool islogy_[]={kFALSE,kTRUE,kTRUE};

int overlayfile[][4]={{2,3,1,0},{4,5,1,0},{6,7,1,0},{8,9,1,0},{10,11,1,0},{12,13,1,0},{14,15,1,0},{16,17,1,0},{18,19,1,0},{20,21,1,0},{22,23,0,1},{22,24,0,1},{23,24,1,1}};
int overlayfile3[][4]={{22,23,24,1}};
int overlayhisto[]={0,1,2};

int nfile=(int)(sizeof(infilename)/sizeof(TString));
int nhisto=(int)(sizeof(histoname_)/sizeof(TString));
int noverlayfile=(int)(sizeof(overlayfile)/sizeof(int))/4;
int noverlayfile3=(int)(sizeof(overlayfile3)/sizeof(int))/4;
int noverlayhisto=(int)(sizeof(overlayhisto)/sizeof(int));

//----------------------------------------------------------pt eta-------------------------------------------------------

TString datafile ="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/PPbCombinedJetTrigAkPu3PFJetHFsumEta4Bin1AllHist.root";
//TString datafile = "~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/CombinedJetTrigAkPu3PFJetAllHist.root";

//TString PPbfile = "~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/mergedCSdiff_MCPPb_AkPu3PFIncJet_2013HistIterTrkCorrtest0EtaCut10_HP04_prod16_v77_merged_forest_0.root";
TString PPbfile = "~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/mergedCSdiff_MCPPb_AllAkPu3PFIncJetNoResidual2013HistIterTrkCorrtest0EtaCut10_pthatcutnewbin_HP04_prod16_v77_merged_forest_0.root";


//TString PPfile = "~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/mergedCSdiff_MCPP_AkPu3PFIncJetNoResidual_2013HistIterTrkCorrtest0EtaBin7_HFsumEta4Bin1_HP04_hiforest77_hiSignal.root";
//TString PPfile = "~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/mergedCSdiff_MCPP_AkPu3PFIncJet_2013HistIterTrkCorrtest0EtaCut10_HP04_hiforest77_hiSignal.root";
TString PPfile = "~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/mergedCSdiff_MCPP_AllAkPu3PFIncJetNoResidual2013HistIterTrkCorrtest0EtaCut10_pthatcutnewbin_HP04_hiforest77_hiSignal.root";
