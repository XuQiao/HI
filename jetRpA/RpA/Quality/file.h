TString algo, datafile, PPfile, PPbfile;

algo="ak3PF";


//if(algo=="ak3PF"){
//Ttring datafile=Form("~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/akPu3PF/old/mergedCSdiff_MCPP_AllAkPu3PFIncJetNoResidual2013HistIterTrkCorrtest0EtaCut10_newbin_HP04_hiforest77_hiSignal.root";

datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/CombinedJetTrigak3PFJetAllHist.root";

PPfile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/mergedCSdiff_MCPP_ak3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_hiforest77_hiSignal.root";
//TString PPfile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/akPu3PF/old/mergedCSdiff_MCPP_AllAkPu3PFIncJetNoResidual2013HistIterTrkCorrtest0EtaCut10_newbin_HP04_hiforest77_hiSignal.root";
//TString PPfile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/akPu3PF/NoResidual/mergedCSdiff_MCPP_akPu3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_hiforest77_hiSignal.root";

PPbfile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/mergedCSdiff_MCPPb_ak3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";
//TString PPbfile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/akPu3PF/old/mergedCSdiff_MCPPb_AllAkPu3PFIncJetNoResidual2013HistIterTrkCorrtest0EtaCut10_newbin_HP04_prod16_v77_merged_forest_0.root";
//TString PPbfile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/akPu3PF/NoResidual/mergedCSdiff_MCPPb_akPu3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";

//}

TString JetIDNameList[18]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN"};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
double binbound_JetID[];
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
double xrange_pt[2]={84+1e-4,97-1e-4};

TFile *fdata=TFile::Open(datafile);
TFile *fPPb=TFile::Open(PPbfile);
TFile *fPP=TFile::Open(PPfile);
int ilist=atoi(getenv("LIST"));
//int ilist=7;
TString JetID,Unit;     TString JetIDName=JetIDNameList[ilist];
TString histonameID=Form("jetpt%s_0-100%%",JetIDName.Data());
//TString histonameID=Form("jetptEta_0-100%%");
TH2F* hdata2F=(TH2F*)fdata->Get(histonameID);
TH2F* hPPb2F=(TH2F*)fPPb->Get(histonameID);
TH2F* hPP2F=(TH2F*)fPP->Get(histonameID);
