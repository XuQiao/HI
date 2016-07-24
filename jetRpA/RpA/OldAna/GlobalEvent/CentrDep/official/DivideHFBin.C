void DivideHFBin(){
double goal;
bool IsMC=kTRUE, Isweight=kFALSE;
if(IsMC){
if(Isweight){
TString dir="/afs/cern.ch/user/q/qixu/CMSSW_5_3_8_HI/src/jetRpA/RpA/CentrDep/output/weight";
TString inputfile="Closure_Hijing_DS_bak.root";}
else{
TString dir="/afs/cern.ch/user/q/qixu/CMSSW_5_3_8_HI/src/jetRpA/RpA/CentrDep/output/noweight";
TString inputfile="Closure_Hijing_DS_bak.root";}
}
else{
TString dir="/afs/cern.ch/user/q/qixu/CMSSW_5_3_8_HI/src/jetRpA/RpA/output/JetTrig/CentrDep/";
if(!Isweight)
TString inputfile="MBCen_noweight.root";
else
TString inputfile="MBCen_weight_bak_1.root";
}
ofstream fstr;
//fstr<<setprecision(2)<<fixed<<endl;
TFile *f=TFile::Open(Form("%s/%s",dir.Data(),inputfile.Data()));
cout<<dir<<"/"<<inputfile<<endl;
if(IsMC){
TString Var[4]={"HFEnergy4","HFEnergyPlus","HFEnergyPlus4","Ntrk"};
TString corr[5]={"_DS","_corr_Eric","_corr_QiaoNtrk","_corr_QiaoHFPlus4","_corr_QiaoReHFPlus4"};
TH1F* hHF[4][5];
}
else{
TString Var[4]={"HF4","HFPlus","HFPlus4","Ntrk"};
TString corr[4]={"Wt","_QiaoNtrkWt","_QiaoHFPlus4Wt","_QiaoReHFPlus4Wt"};
TH1F* hHF[4][4];
}
if(IsMC){
if(Isweight)
fstr.open(Form("HFbin_Hijing_weight_bak.txt"));
else
fstr.open(Form("HFbin_Hijing_noweight_bak.txt"));
}
else{
if(Isweight)
fstr.open(Form("HFbin_Data_weight_bak.txt"));
else
fstr.open(Form("HFbin_Data_noweight.txt"));
}
const int nVar=4;	TH1F* hHF_tr[nVar];
int ncorr=sizeof(corr)/sizeof(TString)-1;
const int N=5;
double cent[N]={0.1,0.2,0.3,0.6,0.9};
for(int k=0;k<nVar;k++){
for(int j=0;j<=ncorr;j++){
if(IsMC){
hHF[k][j]=(TH1F*)f->Get(Form("h%s%s",Var[k].Data(),corr[j].Data()));
hHF_tr[k]=(TH1F*)f->Get(Form("h%s_tr",Var[k].Data()));}
else{
hHF[k][j]=(TH1F*)f->Get(Form("pfjet_%s%s",Var[k].Data(),corr[j].Data()));
hHF_tr[k]=(TH1F*)f->Get(Form("pfjet_%s",Var[k].Data()));
}
fstr<<hHF[k][j]->GetName()<<" Integral="<<hHF_tr[k]->Integral()/hHF[k][j]->Integral()<<endl;
for(int icent=0;icent<N;icent++){
//fstr<<"start from 90-100%:"<<endl;
for(int ibin=1;ibin<hHF[k][j]->GetNbinsX();ibin++){
goal=cent[icent]*hHF[k][j]->Integral();
if(hHF[k][j]->Integral(ibin,hHF[k][j]->GetNbinsX())<=goal){
fstr<<icent<<"th: "<<hHF[k][j]->GetXaxis()->GetBinLowEdge(ibin)<<'\t'<<cent[icent]<<endl;
break;
}
}
}
}
}

}
