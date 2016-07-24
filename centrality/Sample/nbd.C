{

TF1 *fit1 = new TF1("fit1","[0]*(TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1])))*(TMath::Power(([2]/[1]),x))*(TMath::Power((1+([2]/[1])),-x-[1]))",0,160);
fit1->FixParameter(0,1.);// normalization constant
fit1->FixParameter(1,0.464); // k parameter
fit1->FixParameter(2,3.038);  // mean multiplicity
  Float_t Ncoll;
  Float_t Npart;
  Float_t B;
  double BBCcharge=0;
  TH1D *BBC = new TH1D("BBC","",140,0,140);

     char *file1 = "dAu/glauber_dAu.root";
     TFile *f1 = new TFile(file1);
     TTree *t1 = (TTree*)f1->Get("nt_d_Au");
     t1->SetBranchAddress("Ncoll", &Ncoll);
     t1->SetBranchAddress("Npart", &Npart);
     t1->SetBranchAddress("B", &B);
     long neve=t1->GetEntries();
     cout<<"Nevts = "<<neve<<endl;
     //neve= 10000;
     for(long ne=0; ne<neve; ne++){
       if(ne%10000==0)  cout<<"Have run "<<ne<<" events "<<endl;
       t1->GetEntry(ne);
       //cout<<"Ncoll="<<Ncoll<<endl;
       BBCcharge=0;
       for(int i=0;i<Ncoll;i++){
        BBCcharge+=fit1->GetRandom();
       }
       //cout<<"Ncoll="<<Ncoll<<",  BBCcharge="<<BBCcharge<<" ,  BBCcharge/Ncoll="<<BBCcharge/Ncoll<<endl;
       BBC->Fill(BBCcharge);
      }

const int N=4;
double binv[N+1]; binv[0]=0;
int binh[N+1]; binh[0]=2;
int color[N]={5,4,3,2};
double sum=0;
for(int i=0;i<N;i++)  binv[i+1]=(i+1)*1.0*neve/4;
for(int j=0;j<BBC->GetNbinsX();j++){
sum+=BBC->GetBinContent(j);
for(int i=0;i<N;i++) if((sum>binv[i+1])&&(sum-BBC->GetBinContent(j))<binv[i+1]){ binh[i+1]=j;cout<<"j="<<j<<endl; }
}

//cout<<"BBC->GetNbinsX()="<<BBC->GetNbinsX()<<endl;
//cout<<BBC->GetBinContent(3)<<endl;
TCanvas = new TCanvas();
c1->SetLogy();
BBC->Draw();
for(int i=0;i<N;i++){
   TH1D *h1 = (TH1D*)BBC->Clone();
   h1->SetFillColor(color[i]);
   h1->GetXaxis()->SetRange(binh[i],binh[i+1]);
   h1->Draw("same");
}
BBC->GetXaxis()->SetTitle("BBC Charge South");
c1->Print("bbcCharge.gif");
c1->Print("bbcCharge.pdf");

}


