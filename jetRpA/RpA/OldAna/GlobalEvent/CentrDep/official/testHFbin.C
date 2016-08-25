void testHFbin(){
TFile *f=TFile::Open("output/weight/pPbHist_Hijing_DS.root");
TH1F* hHFEnergyPlus4_DS=(TH1F*)f->Get("hHFEnergyPlus4_DS");
const int N=50;
double cent[N+1],HFbin[N+1];
for(int icent=N;icent>=0;icent--){
cent[icent]=(double)icent/N;
cout<<icent<<":"<<cent[icent]<<endl;
}
cout<<"----------------"<<endl;
for(int icent=0;icent<=N;icent++){
//fstr<<"start from 90-100%:"<<endl;
for(int ibin=1;ibin<hHFEnergyPlus4_DS->GetNbinsX();ibin++){
double goal=cent[N-icent]*hHFEnergyPlus4_DS->Integral();
if(hHFEnergyPlus4_DS->Integral(ibin,hHFEnergyPlus4_DS->GetNbinsX())<=goal){
cout<<icent<<"th: "<<hHFEnergyPlus4->GetXaxis()->GetBinLowEdge(ibin)<<'\t'<<cent[icent]<<endl;
break;
}
}
HFbin[icent]=hHFEnergyPlus4->GetXaxis()->GetBinLowEdge(ibin);
}
cout<<"----------------"<<endl;
for(int icent=0;icent<=N;icent++){
cout<<HFbin[icent]<<",";
}
cout<<"----------------"<<endl;
int nHFbin=sizeof(HFbin)/sizeof(double)-1;
}
