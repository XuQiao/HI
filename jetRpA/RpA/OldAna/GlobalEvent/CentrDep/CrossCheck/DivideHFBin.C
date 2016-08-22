
double partition(double array[],int low, int high){

double pivot=array[low];
while(low<high){
while(low<high&&array[high]>=pivot)
high--;
array[low]=array[high];
while(low<high&&array[low]<=pivot)
low++;
array[high]=array[low];
}
array[low]=pivot;
return low;
}

double findkth(double array[],int low,int high, int k){

double kth;
if(low == high) kth = array[low];
else{
 double pivot = array[low];
 int splitpoint = partition(array,low,high);
if(k-1 == splitpoint)   kth = pivot;
else if(k-1<splitpoint) kth = findkth(array,low,splitpoint-1,k);
else if(k-1>splitpoint) kth = findkth(array,splitpoint+1,high,k);
}


return kth;
}

void DivideHFBin(){
double goal;
TString dir=".";
TString inputfile="pPbHist_Hijing_DS";
TFile *f=TFile::Open(Form("%s/%s.root",dir.Data(),inputfile.Data()));
TH1F* hHF[3];
if(inputfile.Contains("Pbp")){
TString Var[3]={"HFEnergy4","HFEnergyMinus","HFEnergyMinus4"};
if(inputfile.Contains("withGplus"))
ofstream fstr(Form("HFbin_Pbp_withGplus.txt"));
else
ofstream fstr(Form("HFbin_Pbp.txt"));
}
else{
TString Var[3]={"HFEnergy4","HFEnergyPlus","HFEnergyPlus4"};
if(inputfile.Contains("withGplus"))
ofstream fstr(Form("HFbin_pPb_withGplus.txt"));
else
ofstream fstr(Form("HFbin_pPb.txt"));
}

fstr<<setprecision(2)<<fixed;
const int N=5;
double cent[N]={0.1,0.2,0.3,0.6,0.9};
for(int k=0;k<3;k++){
hHF[k]=(TH1F*)f->Get(Form("h%s_tr",Var[k].Data()));
cout<<hHF[k]->GetName()<<" Integral="<<hHF[k]->Integral()<<endl;
fstr<<hHF[k]->GetName()<<endl;
for(int j=0;j<N;j++){
//fstr<<"start from 90-100%:"<<endl;
for(int i=1;i<hHF[k]->GetNbinsX();i++){
goal=cent[j]*hHF[k]->Integral();
if(hHF[k]->Integral(i,hHF[k]->GetNbinsX())<=goal){
fstr<<j<<"th: "<<hHF[k]->GetXaxis()->GetBinLowEdge(i)<<'\t'<<goal<<endl;
break;
}
}
/*
fstr<<"start from 0-10%:"<<endl;
for(int i=1;i<hHF[k]->GetNbinsX();i++){
goal=(1-cent[N-j-1])*hHF[k]->Integral();
if(hHF[k]->Integral(1,i)>=goal){
fstr<<j<<"th: "<<hHF[k]->GetXaxis()->GetBinCenter(i)<<'\t'<<goal<<endl;
break;
}
}
*/
}


}
/*
TString Epos_file = "/store/user/tuos/pPb_MC_MinBiasTree_v4_Epos_189k.root";
TString Hijing_file="/store/user/tuos/pPb_MC_MinBiasTree_v4_Hijing_180k.root";
TString Var="PlusEta4";
ofstream fstr(Form("HFbin_DS_%s.txt",Var.Data()));
TFile *fin=TFile::Open(Hijing_file);
fin->cd();
TTree *tree;
tree=(TTree*)hiEvtAnalyzer->Get("HiTree");
tree->AddFriend("skimTree=skimanalysis/HltTree");
tree->AddFriend("hltTree=hltanalysis/HltTree");
tree->AddFriend("genpartTree=genparticles/hi");

const int Maxmult=10000;
Float_t hiHF, vz, hiHFplus,eta[Maxmult];
Float_t hiHFminusEta4, hiHFplusEta4,HFEnergy;
Long_t Ev,Nevent;
Int_t n, HLT_PAZeroBiasPixel_SingleTrack_v1, pprimaryvertexFilter, pVertexFilterCutGplus, pPAcollisionEventSelectionPA, pHBHENoiseFilter,phfPosFilter1, phfNegFilter1,pBeamScrapingFilter;
Int_t Nskim=0;  double weight;  int mult, nDS=0, nPlus35, nMinus35;     int sta[Maxmult];
Nevent=tree->GetEntries();

tree->SetBranchAddress("hiHF",&hiHF);//HF energy
tree->SetBranchAddress("hiHFplus",&hiHFplus);//HF energy positive
tree->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);//HF energy positive 4 to 5
tree->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);//HF energy negative -5 to -4
tree->SetBranchAddress("vz",&vz);
tree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
tree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
tree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
//tree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
tree->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
tree->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
tree->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
tree->SetBranchAddress("eta",&eta);
tree->SetBranchAddress("sta",&sta);
tree->SetBranchAddress("mult",&mult);
float Minx=0., Maxx=100.;	double centbin[]={0.1,0.2,0.3,0.6,0.9};
const int N=5;
Int_t i, j;
double array[1000000];
Long_t count_Ev=0;
for(Ev=0;Ev<Nevent;Ev++){

tree->GetEntry(Ev);
weight=1.;
//weight*=fVz->Eval(vz);
if(Ev%5000==0)
cout<<"Have run "<<Ev<<" Events"<<endl;
//fstr<<"weight="<<weight<<endl;
tree->GetEntry(Ev);
if(Var=="PlusEta4")
HFEnergy=hiHFplusEta4;
else if(Var=="Plus")
HFEnergy=hiHFplus;
else if(Var=="Eta4")
HFEnergy=hiHFplusEta4+hiHFminusEta4;
else{exit();}
nPlus35=0; nMinus35=0;
for(int j=0; j<mult;j++){
if(eta[j]>3&&eta[j]<5&&sta[j]==1) nPlus35++;
if(eta[j]>-5&&eta[j]<-3&&sta[j]==1) nMinus35++;
}
if(nPlus35>=1 && nMinus35>=1 && TMath::Abs(vz)<15){
if(HFEnergy>=Minx&&HFEnergy<=Maxx){
count_Ev++;
array[count_Ev-1]=HFEnergy;
//if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<15))
//if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && TMath::Abs(vz)<15))
}
}
}

fstr<<count_Ev<<endl;
int low=0;
int high=count_Ev-1;
int k0[N+1];
double kpoint[N+1];
kpoint[0]=Minx;
kpoint[N]=Maxx;
fstr<<"i=0"<<'\t'<<"kpoint="<<kpoint[0]<<'\t'<<endl;
for(i=1;i<N;i++){
k0[i]=count_Ev*(1-centbin[N-i]);
kpoint[i]=findkth(array,low,high,k0[i]);
fstr<<"i="<<i<<'\t'<<"kpoint="<<kpoint[i]<<'\t'<<endl;
count_Ev=0;
for(Ev=0;Ev<Nevent;Ev++){
tree->GetEntry(Ev);
if(Var=="PlusEta4")
HFEnergy=hiHFplusEta4;
else if(Var=="Plus")
HFEnergy=hiHFplus;
else if(Var=="Eta4")
HFEnergy=hiHFplusEta4+hiHFminusEta4;
else{exit();}
nPlus35=0; nMinus35=0;
for(int j=0; j<mult;j++){
if(eta[j]>3&&eta[j]<5&&sta[j]==1) nPlus35++;
if(eta[j]>-5&&eta[j]<-3&&sta[j]==1) nMinus35++;
}
if(nPlus35>=1 && nMinus35>=1 && TMath::Abs(vz)<15){
if(HFEnergy>=Minx&&HFEnergy<=Maxx){
count_Ev++;
array[count_Ev-1]=HFEnergy;
}
}
}
}
fstr<<"i="<<N<<'\t'<<"kpoint="<<kpoint[i]<<'\t'<<endl;

int PartEvent[N],TotalEvent;

for (Ev=0; Ev<Nevent; Ev++){

                      if(Ev%10000==0) cout<<"Have run "<<Ev/1000<<" thousand events"<<endl;

                        tree->GetEntry(Ev);
if(Var=="PlusEta4")
HFEnergy=hiHFplusEta4;
else if(Var=="Plus")
HFEnergy=hiHFplus;
else if(Var=="Eta4")
HFEnergy=hiHFplusEta4+hiHFminusEta4;
else{exit();}
nPlus35=0; nMinus35=0;
for(int j=0; j<mult;j++){
if(eta[j]>3&&eta[j]<5&&sta[j]==1) nPlus35++;
if(eta[j]>-5&&eta[j]<-3&&sta[j]==1) nMinus35++;
}
if(nPlus35>=1 && nMinus35>=1 && TMath::Abs(vz)<15){
                        for(i=0; i<N; i++)
                        if(HFEnergy>kpoint[i]&&HFEnergy<=kpoint[i+1]){
                                PartEvent[i]++;
                                TotalEvent++;
                        }
}
}
for(i=0;i<N;i++)
fstr<<"SumEvent Percent="<<(double)PartEvent[i]/TotalEvent<<endl;

fin->Close();
*/

}
