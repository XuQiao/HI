TString algo, datafile, PPfile, PPbfile;

algo="akPu3PF";


TString datafile="/home/xuq/Documents/HI/RpA/output/TreeAna/JetTrigvzandcent.root";
TString PPfile="/home/xuq/Documents/HI/RpA/output/akPu3PF/MCPPakPu3PFjetana_merged.root";
TString PPbfile="/home/xuq/Documents/HI/RpA/output/akPu3PF/MCPPbakPu3PFjetana_merged.root";

TString JetIDNameList[18]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN"};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
double xrange_pt[2]={20+1e-4,600-1e-4};
double xrange_pt_pre[2]={15,692};

TFile *fdata=TFile::Open(datafile);
TFile *fPPb=TFile::Open(PPbfile);
TFile *fPP=TFile::Open(PPfile);
