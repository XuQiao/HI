void calccorr(int i=1,string t="ptfiner"){
float const PI = acos(-1.0);
if(i==0){
const int ncent = 6;
const int npt = 25;
TString dire="north";
}
else if(i==1){
const int ncent = 6;
const int npt = 25;
TString dire="south";
}
else if(i==2){
const int ncent = 6;
const int npt = 1;
TString dire="sn";
}

TFile *f=TFile::Open("../../../work/39GeV/output_3corr.root");
TH1F* kforebbcw[ncent][npt];
TH1F* hforebbcw[ncent][npt];
TH1F* kbackbbcw2[ncent][npt];
TH1F* kforebbcw_In;
TH1F* hforebbcw_In;
TH1F* kbackbbcw2_In;
double ptbin[26] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
double centbin[7] = {0,0.05,0.1,0.2,0.4,0.6,1.0};
//double centbin[ncent+1] = {0,1,2};
TString type = t.c_str();
if(type=="ptIn25_4"){
double selptbin[] = {2.5,4.0};
double selcentbin[ncent+1] = {0,0.01,0.05,0.1,0.2,0.3,0.4,0.6,1.0};
}
else if(type=="ptIn"){
    if(i==2){
double selptbin[] = {0.0,0.2};
    }
    else{
double selptbin[] = {0.4,3.0};
    }
double selcentbin[] = {0,0.05,0.1,0.2,0.4,0.6,1.0};
}
else if(type=="ptcoarser"){
double selptbin[] = {0.2,1.0,2.0,3.0,5.0};
double selcentbin[] = {0,0.01,0.05,0.1,0.2,0.3,0.4,0.6,1.0};
}
else if(type=="ptfiner"){
double selptbin[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
double selcentbin[] = {0,0.05,0.1,0.2,0.4,0.6,1.0};
//double selcentbin[ncent+1] = {0,1,2};
}
else if(type=="centIn"){
double selptbin[] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
double selcentbin[] = {0,1.0};
}
else if(type=="ptccentc"){
double selptbin[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
double selcentbin[] = {0,0.05};
}
else{exit(0);};

/*double ptmin = 1.0; double ptmax = 3.0;//has to be boundary
for(int ipt=0; ipt<npt; ipt++){
if(ptmin >= ptbin[ipt] && ptmin < ptbin[ipt+1]){int xptmin = ipt; continue;}
if(ptmax >= ptbin[ipt] && ptmax < ptbin[ipt+1]){int xptmax = ipt; continue;}
}*/
  ofstream fout(Form("c1_c2_central_%s_%s.dat",type.Data(),dire.Data()));
for(int icent=0; icent<ncent; icent++){
for(int ipt=0; ipt<npt; ipt++){
kforebbcw[icent][ipt] = (TH1F*)f->Get(Form("kfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
//htemp = (TH1F*)f->Get(Form("kforenorthbbcw_%d_%d",icent,ipt));
//kforebbcw[icent][ipt]->Add(htemp,-1);
hforebbcw[icent][ipt] = (TH1F*)f->Get(Form("hfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
//htemp = (TH1F*)f->Get(Form("hforenorthbbcw_%d_%d",icent,ipt));
//hforebbcw[icent][ipt]->Add(htemp,-1);
kbackbbcw2[icent][ipt] = (TH1F*)f->Get(Form("kback%sbbcw2_%d_%d",dire.Data(),icent,ipt));
//htemp = (TH1F*)f->Get(Form("kbacknorthbbcw2_%d_%d",icent,ipt));
//kbackbbcw2[icent][ipt]->Add(htemp,-1);
}
}
int ncent_a = sizeof(selcentbin)/sizeof(double)-1;
int npt_a = sizeof(selptbin)/sizeof(double)-1;
for(int icent_a=0;icent_a<ncent+1;icent_a++){
for(int icent_b=0; icent_b<ncent+1; icent_b++)
if(selcentbin[icent_a] == centbin[icent_b]) break;
int xcentmin = icent_b;
for(int icent_b=0; icent_b<ncent+1; icent_b++)
if(selcentbin[icent_a+1] == centbin[icent_b]) break;
int xcentmax = icent_b;
if((xcentmin == ncent+1) || (xcentmax == ncent+1)) exit(0);
for(int ipt_a=0;ipt_a<npt_a;ipt_a++){
for(int ipt_b=0; ipt_b<npt+1; ipt_b++)
if(selptbin[ipt_a] == ptbin[ipt_b]) break;
int xptmin = ipt_b;
for(int ipt_b=0; ipt_b<npt+1; ipt_b++)
if(selptbin[ipt_a+1] == ptbin[ipt_b]) break;
int xptmax = ipt_b;
if((xptmin == npt+1) || (xptmax == npt+1)) exit(0);
cout<<xcentmin<<"\t"<<xcentmax<<"\t"<<endl;
cout<<xptmin<<"\t"<<xptmax<<"\t"<<endl;
kforebbcw_In = (TH1F*)kforebbcw[xcentmin][xptmin]->Clone();
hforebbcw_In = (TH1F*)hforebbcw[xcentmin][xptmin]->Clone();
kbackbbcw2_In = (TH1F*)kbackbbcw2[xcentmin][xptmin]->Clone();
kforebbcw_In->Reset();
hforebbcw_In->Reset();
kbackbbcw2_In->Reset();
for(int icent=xcentmin; icent<xcentmax; icent++){
for(int ipt=xptmin; ipt<xptmax; ipt++){
kforebbcw_In->Add(kforebbcw[icent][ipt]);
hforebbcw_In->Add(hforebbcw[icent][ipt]);
kbackbbcw2_In->Add(kbackbbcw2[icent][ipt]);
}
}
kforebbcw_In->Rebin(2);
hforebbcw_In->Rebin(2);
kbackbbcw2_In->Rebin(2);


TH1F* hpp;
TH1F* hbackpp;
  hpp = (TH1F*)kforebbcw_In->Clone();
  hbackpp = (TH1F*)kbackbbcw2_In->Clone();
  float nbackpp = hbackpp->Integral()/2.0/PI;
  float nforepp = hpp->Integral()/2.0/PI;
  //float ntrig0 = ptforedis_0->Integral(11,30);
   for(int i=0; i<20; i++){
     float pp_cont = 1.0*hpp->GetBinContent(i+1);
     //float pp0_err = 1.0*hpp->GetBinError(i+1);
     float weight2 = sqrt(1.0*hforebbcw_In->GetBinContent(i+1));

     float backpp_cont = 1.0*hbackpp->GetBinContent(i+1);
    
     float con = pp_cont/backpp_cont*nbackpp/nforepp;
     float err = weight2/backpp_cont*nbackpp/nforepp;

     hpp->SetBinContent(i+1, con);
     hpp->SetBinError(i+1, err);

   }
  TF1 *fun0 = new TF1("fun0","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)+2*[4]*cos(4*x))", -0.5*PI, 1.5*PI);

  fun0->SetLineColor(1);
  hpp->Fit("fun0","NORQ");

  TF1 *fun1 = new TF1("fun1","[0]*(1+2*[1]*cos(x))",   -0.5*PI, 1.5*PI);
  TF1 *fun2 = new TF1("fun2","[0]*(1+2*[1]*cos(2*x))", -0.5*PI, 1.5*PI);
  TF1 *fun3 = new TF1("fun3","[0]*(1+2*[1]*cos(3*x))", -0.5*PI, 1.5*PI);
  TF1 *fun4 = new TF1("fun4","[0]*(1+2*[1]*cos(4*x))", -0.5*PI, 1.5*PI);

//const double ptmean[npt] = {0.360943, 0.691833, 1.1911, 1.69654, 2.20117, 2.70571, 3.2097, 3.71372, 4.21814, 4.72014};
  fout<<fun0->GetParameter(1)<<" "<<fun0->GetParError(1)<<" "
      <<fun0->GetParameter(2)<<" "<<fun0->GetParError(2)<<" "//<<0.05*fun0->GetParameter(2)<<" "<<endl;
      <<fun0->GetParameter(3)<<" "<<fun0->GetParError(3)<<endl;
}
}
  fout.close();

/*  
  cout<<"*************** v2 ***********"<<endl;
  float v2_0 = funvn0->GetParameter(1)/(zym_pp0 + funvn0->GetParameter(0));
  float v2_1 = funvn1->GetParameter(1)/(zym_pp1 + funvn1->GetParameter(0));
  float v2_2 = funvn2->GetParameter(1)/(zym_pp2 + funvn2->GetParameter(0));
  float v2_3 = funvn3->GetParameter(1)/(zym_pp3 + funvn3->GetParameter(0));
  

  cout<<v2_0<<" "<<v2_1<<" "<<v2_2<<" "<<v2_3<<endl;
 
  
  cout<<funvn0->GetParameter(0)*funvn0->GetParameter(1)<<" "
      <<funvn1->GetParameter(0)*funvn1->GetParameter(1)<<" "
      <<funvn2->GetParameter(0)*funvn2->GetParameter(1)<<" "
      <<funvn3->GetParameter(0)*funvn3->GetParameter(1)<<endl;
  

  cout<<"*************** v2 ***********"<<endl;
  cout<<funvn0->GetParameter(2)<<" "<<funvn1->GetParameter(2)<<" "<<funvn2->GetParameter(2)<<" "<<funvn3->GetParameter(2)<<endl;

  cout<<"*************** v3 ***********"<<endl;
  cout<<funvn0->GetParameter(3)<<" "<<funvn1->GetParameter(3)<<" "<<funvn2->GetParameter(3)<<" "<<funvn3->GetParameter(3)<<endl;
*/
}
