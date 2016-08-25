TString algo, datafile, PPfile, PPbfile;

//algo="akPu3PF";
algo="akPu3PF";


//TString PPbfile="/scratch/xuq7/RpA/JetID/akPu3PF/NoResidual/NoGplusNewpthat/newSample/MCPPbNoGpluscW_akPu3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root";
TString PPbfile="/scratch/xuq7/RpA/akPu3PF/NoResidual/NoGplusNewpthat/newSample/MCPPbNoGpluscW_akPu3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root";
TString PPbfile="/scratch/xuq7/RpA/akPu3PF/NoResidual/NoGplusNewpthat/newSample/MCPPbakPu3PFjetana_merged.root";
//TString PPbfileOld="/scratch/xuq7/RpA/akPu3PF/NoResidual/NoGplusNewpthat/MCPPbNoGplus_akPu3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root"; //Old sample
TString datafile="/scratch/xuq7/RpA/TreeAna/Datacombined.root";

//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
double xrange_pt[2]={50+1e-4,600-1e-4};
double xrange_pt_pre[2]={47,692};
TFile *fPPb=new TFile(PPbfile);
//TFile *fPPbOld=new TFile(PPbfileOld);
TFile *fdata=new TFile(datafile);
