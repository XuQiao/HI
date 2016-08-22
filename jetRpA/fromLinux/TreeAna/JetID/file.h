TString algo="akPu3PF";
TString datafile="/home/xuq/Documents/HI/RpA/output/TreeAna/DATAPPbakPu3PFskimJetID.root";
TString PPbfile="/home/xuq/Documents/HI/RpA/output/TreeAna/MCPPbakPu3PFskimJetID.root";
TString PPfile="/home/xuq/Documents/HI/RpA/output/TreeAna/MCPPakPu3PFskimJetID.root";
//TString PPbfile="/home/xuq/Documents/HI/RpA/output/akPu3PF/MCPPbakPu3PFjetana_merged.root";
//TString PPfile="/home/xuq/Documents/HI/RpA/output/akPu3PF/MCPPakPu3PFjetana_merged.root";
//TString PPbfile="/home/xuq/Documents/HI/RpA/output/akPu3PF/MCPPbNoGplus_akPu3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root";
//TString PPfile="/home/xuq/Documents/HI/RpA/output/akPu3PF/MCPPNoGplus_akPu3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root";

TString JetIDNameList[19]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN","SumSum"};
//double binbound_pt[1000];
//double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
double xrange_pt[2]={30+1e-4,600-1e-4};
double xrange_JetIDcut[2];
double JetIDcut[2];

TFile *fdata=TFile::Open(datafile);
TFile *fPP=TFile::Open(PPfile);
TFile *fPPb=TFile::Open(PPbfile);
TString JetID,Unit;
