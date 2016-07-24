TString algo="akPu3PF";
TString datafile="/cms/store/user/qixu/jetRpA/RpA/NewMC/DATAPPbakPu3PFskimUniv.root";
TString dataJetIDfile="/cms/store/user/qixu/jetRpA/RpA/NewMC/DATAPPbakPu3PFskimJetID.root";
TString PPbfile="/cms/store/user/qixu/jetRpA/RpA/NewMC/MCPPbakPu3PF_useskim.root";
//TString PPfile="/scratch/xuq7/RpA/TreeAna/MCPPbakPu3PF.root";

TString JetIDNameList[]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN","PPcut","PPcutTight","neuMaxr1"};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
double xrange_pt[2]={30+1e-4,600-1e-4};
double xrange_JetIDcut[2];
double JetIDcut[2];

const int Neta=8;
const TString etabinnameswap[Neta]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const TString etabinname[Neta]={"-20_-15","-15_-10","-10_-5","-5_5","5_10","10_15","15_20","-10_10"};
const double etabin[Neta]={0.5,0.5,0.5,1.0,0.5,0.5,0.5,2.0};
const TString etastring[Neta]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};

TFile *fdata=TFile::Open(datafile);
TFile *fdataJetID=TFile::Open(dataJetIDfile);
//TFile *fPP=TFile::Open(PPfile);
TFile *fPPb=TFile::Open(PPbfile);
TString JetID,Unit;
