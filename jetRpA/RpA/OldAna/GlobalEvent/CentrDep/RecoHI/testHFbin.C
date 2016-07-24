void testHFbin(){
TFile *f=TFile::Open("pPbHijing_EffCorr_VRw.root");
TH1F* hHFp4DS=(TH1F*)f->Get("hHFp4DS");
const int N=50;
double cent[N+1],HFbin[N+1];
for(int icent=N;icent>=0;icent--){
cent[icent]=(double)icent/N;
cout<<icent<<":"<<cent[icent]<<endl;
}
cout<<"----------------"<<endl;
for(int icent=0;icent<=N;icent++){
//fstr<<"start from 90-100%:"<<endl;
for(int ibin=1;ibin<hHFp4DS->GetNbinsX();ibin++){
double goal=cent[N-icent]*hHFp4DS->Integral();
if(hHFp4DS->Integral(ibin,hHFp4DS->GetNbinsX())<=goal){
cout<<icent<<"th: "<<hHFp4DS->GetXaxis()->GetBinLowEdge(ibin)<<'\t'<<cent[icent]<<endl;
break;
}
}
HFbin[icent]=hHFp4DS->GetXaxis()->GetBinLowEdge(ibin);
}
cout<<"----------------"<<endl;
for(int icent=0;icent<=N;icent++){
cout<<HFbin[icent]<<",";
}
cout<<"----------------"<<endl;
int nHFbin=sizeof(HFbin)/sizeof(double)-1;
}
