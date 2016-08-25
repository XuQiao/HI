TString algo, datafile, PPfile, PPbfile;

algo="akPu3PF";

TString DataPPbfile="/scratch/xuq7/RpA/fromLinux/TreeAna/DATAPPbakPu3PFskimUniv.root";
TString DataPPbJetIDfile="/scratch/xuq7/RpA/fromLinux/TreeAna/DATAPPbakPu3PFskimJetID.root";
TString DataPbPJetIDfile="/scratch/xuq7/RpA/fromLinux/TreeAna/DATAPbPakPu3PFskimJetID.root";
TString DataPbPfile="/scratch/xuq7/RpA/fromLinux/TreeAna/DATAPbPakPu3PFskimUniv.root";
TString MCPPbfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPbakPu3PFskimUniv.root";
TString MCPbPfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPbPakPu3PFskimUniv.root";
TString MCPPbxSecfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPbakPu3PFxSecskimUniv.root";
//TString MCPbPfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPbPakPu3PFskimUniv.root";

TString MCPPxSecfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPakPu3PFxSecskimUniv.root";
TString MCPPxSecJetIDfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPakPu3PFxSecskimJetID.root";

TString MCPPbJetIDfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPbakPu3PFskimJetID.root";
TString MCPPbxSecJetIDfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPbakPu3PFxSecskimJetID.root";
TString MCPbPJetIDfile="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPbPakPu3PFskimJetID.root";

TString DataPPbfileOld="/scratch/xuq7/RpA/fromLinux/TreeAna/DATAPPbakPu3PFUniv.root";
TString DataPbPfileOld="/scratch/xuq7/RpA/fromLinux/TreeAna/DATAPbPakPu3PFUniv.root";
TString DataPbPJetIDfileOld="/scratch/xuq7/RpA/fromLinux/TreeAna/DATAPbPakPu3PFJetID.root";
TString MCPPbfileOld="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPbakPu3PFUniv.root";
TString MCPPbJetIDfileOld="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPbakPu3PFJetID.root";
TString MCPbPfileOld="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPbPakPu3PFUniv.root";
TString MCPbPJetIDfileOld="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPbPakPu3PFJetID.root";

//TString 
//TString MCPPfileOld="/scratch/xuq7/RpA/fromLinux/TreeAna/MCPPakPu3PFUniv.root";

const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
TString JetIDNameList[19]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN","SumSum"};
double xrange_pt[2]={50+1e-4,600-1e-4};
double xrange_pt_pre[2]={47.5+1e-4,480-1e-4};
TFile *fDataPPb=new TFile(DataPPbfile);
TFile *fDataPbP=new TFile(DataPbPfile);
TFile *fDataPPbJetID=new TFile(DataPPbJetIDfile);
//TFile *fDataPbPJetID=new TFile(DataPbPJetIDfile);
TFile *fMCPPxSec=new TFile(MCPPxSecfile);
TFile *fMCPPxSecJetID=new TFile(MCPPxSecJetIDfile);

TFile *fMCPPb=new TFile(MCPPbfile);
TFile *fMCPPbxSec=new TFile(MCPPbxSecfile);
TFile *fMCPPbJetID=new TFile(MCPPbJetIDfile);
TFile *fMCPPbxSecJetID=new TFile(MCPPbxSecJetIDfile);
TFile *fMCPbP=new TFile(MCPbPfile);
TFile *fMCPbPJetID=new TFile(MCPbPJetIDfile);
//TFile *fDataPbP=new TFile(DataPbPfile);

TFile *fDataPPbOld=new TFile(DataPPbfileOld);
TFile *fDataPbPOld=new TFile(DataPbPfileOld);
TFile *fDataPbPJetIDOld=new TFile(DataPbPJetIDfileOld);
TFile *fMCPPbOld=new TFile(MCPPbfileOld);
TFile *fMCPPbJetIDOld=new TFile(MCPPbJetIDfileOld);
TFile *fMCPbPOld=new TFile(MCPbPfileOld);
TFile *fMCPbPJetIDOld=new TFile(MCPbPJetIDfileOld);

const int Neta=8;
const TString etabinnameswap[Neta]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const TString etabinname[Neta]={"-20_-15","-15_-10","-10_-5","-5_5","5_10","10_15","15_20","-10_10"};
const double etabin[Neta]={0.5,0.5,0.5,1.0,0.5,0.5,0.5,2.0};
const TString etastring[Neta]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};
const int canvas[Neta]={0,1,2,3,6,5,4,7};
static const int nColor = 9;
static const int colorCode[nColor] = {
    1,2,4,46,6,7,8,1,2
//    2, 4, 6, 7, 8, 9, 46,1,2
//    1,2,4,1,4,2,1,1    //for spectra and Asy
 //   2, 4, 6, 7, 8, 9, 46,1
//    2, 4, 6, 1, 8, 9, 4,1
//    4, 9, 8, 1, 6, 4, 2,1  //for RpA color
//    2, 7, 9, 1, 6, 4, 46,1  //for kEtaDep
};
static const int markerCode[nColor] = {
//   33, 34, 29, 21, 30, 28,27,20
    20,34,33,25,27,28,24,21,33
//   34,29,20,25,27,28,30,20   //for asymmetry kAsy
//    29,34,33,25,27,28,30,20   //for spectrum marker
//    27,28,24,21,20,34,33,20   //for
 //   20,34,33,25,27,28,30,20
 //   33, 34, 30, 20, 29, 28,27,20   //for Forward-backward Asymmetry, kEtaDep,
 //   27, 28, 30, 20, 29, 34,33,20  //for RpA
};
