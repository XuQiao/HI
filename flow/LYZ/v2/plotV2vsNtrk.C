ofstream fout("V2value.txt");
void plotV2vsNtrk(){

const int nbin24 = 12;
const double avgtrkbin[nbin24]={44.36,54.33,68.63,88.39,108.2,131.3,162.1,196.6,227.5,247.2,269.2,301.2};
const double V24[nbin24]={0.02965,0.03913,0.04832,0.04941,0.04822,0.04955,0.049,0.04805,0.04709,0.04665,0.04772,0.04797};
const double V24err[nbin24]={0.005967,0.004262,0.001496,0.001351,0.001599,0.0003922,0.0003065,0.0002939,0.0004568,0.0008684,0.001471,0.004329};
const double avgtrkbin006[5]={131.3,162.1,196.6,231.5,270.5};
const double V2006[5]={0.049,0.0474,0.0461,0.0451,0.0469};
const double V2006err[5]={0.00147,0.0013,0.0008684,0.004,0.0052};
TCanvas *c1 = new TCanvas;
TGraphErrors *grv24=new TGraphErrors(nbin24,avgtrkbin,V24,0,V24err);
TGraphErrors *grv2006=new TGraphErrors(5,avgtrkbin006,V2006,0,V2006err);
grv24->SetTitle("V_{2} vs Ntrkoffline");
grv24->GetXaxis()->SetTitle("Ntrkoffline");
grv24->GetYaxis()->SetTitle("V_{2}");
grv24->SetMaximum(0.07);
grv24->SetMinimum(0.0);
grv24->SetMarkerSize(1.2);
grv24->SetMarkerColor(1);
grv24->SetLineColor(1);
grv24->SetMarkerStyle(20);
grv2006->SetMarkerColor(6);
grv2006->SetMarkerSize(1.2);
grv2006->SetLineColor(6);
grv2006->SetMarkerStyle(23);
TGraphErrors *grtr03to6 = plot("tracknormcpt03to6",0,2,20,0);
TGraphErrors *grtr03to6fin4 = plot("tracknormcpt03to6/finalbins5",0,4,27,4);
TGraphErrors *grtr03to6fin = plot("tracknormcpt03to6/finalbins5",0,2,30,0);
//TGraphErrors *grtr03to6d3 = plot("tracknormcpt03to6",0,2,30,3);
//TGraphErrors *grtheta10tr03to6 = plot("theta10/tracknormcpt03to6",0,2,30,0);
//TGraphErrors *grtheta15tr03to6 = plot("theta15/tracknormcpt03to6",0,7,27,0);
//TGraphErrors *grtr03to3 = plot("tracknormcpt03to3",0,2,30,0);
//TGraphErrors *grtr01to10 = plot("tracklcpt01to10",0,2,24,0);
//TGraphErrors *gr1Sum = plot("tracknormcpt03to6",1,1,24,0);
//TGraphErrors *gr11 = plot("tracknormcpt03to6",0,3,29,1);
//TGraphErrors *gr12 = plot("tracknormcpt03to6",0,4,30,2);
//TGraphErrors *gr13 = plot("tracknormcpt03to6",0,3,21,3);
//TGraphErrors *gr1f = plot("tracknormcpt03to6/finalbins5",2,24,1);
//TGraphErrors *gr2 = plot("tracknormcpt03to3",4,29,1);
//TGraphErrors *gr2f = plot("tracknormcpt03to3/finalbins5",4,30,1);
//TGraphErrors *gr3 = plot("tracklcpt01to10",2,34,0);
//TGraphErrors *grpf = plot("PFcandpt01to10",4,27,0);
//TGraphErrors *grpf03to6 = plot("PFcandpt03to6",0,4,20,0);
//TGraphErrors *grpf03to3 = plot("PFcandpt03to3tracknormcpt03to6",0,4,30,0);
//TGraphErrors *grpf01to10 = plot("PFcandpt01to10",0,4,24,0);
//TGraphErrors *grpfon = plot("PFcandpt01to10/OnlyCharged",4,33);
TLegend *leg = new TLegend(0.2,0.1,0.4,0.45);
leg->SetTextSize(0.045);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(grtr03to6,Form("track, %.1f < p_{T} < %.1f (GeV/c)",0.3,6.0),"P");
leg->AddEntry(grtr03to6fin4,Form("track final bins, %.1f < p_{T} < %.1f (GeV/c)",0.3,6.0),"P");
leg->AddEntry(grtr03to6fin,Form("track final bins avg, %.1f < p_{T} < %.1f (GeV/c)",0.3,6.0),"P");
//leg->AddEntry(grtheta10tr03to6,Form("track 10 #theta s, %.1f < p_{T} < %.1f (GeV/c)",0.3,6.0),"P");
//leg->AddEntry(grtheta15tr03to6,Form("track 15 #theta s, %.1f < p_{T} < %.1f (GeV/c)",0.3,6.0),"P");
//leg->AddEntry(grtr03to6d3,Form("track with sp and eff corr, %.1f < p_{T} < %.1f (GeV/c)",0.3,6.0),"P");
//leg->AddEntry(grtr01to10,Form("track, %.1f < p_{T} < %.1f (GeV/c)",0.1,10.0),"P");
//leg->AddEntry(gr1Sum,Form("track normal cut Sum reference"),"P");
//leg->AddEntry(gr10,Form("%.1f < p_{T} < %.1f (GeV/c)",0.3,6.0),"");
//leg->AddEntry(gr11,Form("track normal cut integral with v2"),"P");
//leg->AddEntry(gr12,Form("track normal cut with eff corr"),"P");
//leg->AddEntry(gr13,Form("track normal cut with sp and eff corr"),"P");
//leg->AddEntry(gr13,Form("%.1f < p_{T} < %.1f (GeV/c)",0.3,3.0),"");
//leg->AddEntry(gr1f,Form("track normal cut final,%.1f < p_{T} < %.1f (GeV/c)",0.3,6.0),"P");
//leg->AddEntry(gr2,Form("track normal cut,%.1f < p_{T} < %.1f (GeV/c)",0.3,3.0),"P");
//leg->AddEntry(gr2f,Form("track normal cut final,%.1f < p_{T} < %.1f (GeV/c)",0.3,3.0),"P");
//leg->AddEntry(gr3,Form("track loose cut,%.1f < p_{T} < %.1f (GeV/c)",0.1,10.0),"P");
//leg->AddEntry(grpf,Form("PF full,%.1f < p_{T} < %.1f (GeV/c)",0.1,10.0),"P");
//leg->AddEntry(grpf03to6,Form("PF candidate,%.1f < p_{T} < %.1f (GeV/c)",0.3,6.0),"P");
//leg->AddEntry(grpf03to3,Form("PF candidate,%.1f < p_{T} < %.1f (GeV/c)",0.3,3.0),"P");
//leg->AddEntry(grpf01to10,Form("PF candidate,%.1f < p_{T} < %.1f (GeV/c)",0.1,10.0),"P");
//leg->AddEntry(grpfon,Form("PF only charged,%.1f < p_{T} < %.1f (GeV/c)",0.1,10.0),"P");
leg->AddEntry(grv24,"4-particle cumulant","P");
//leg->AddEntry(grv2006,"HIN-14-006","P");
grv24->Draw("AP");
//grv2006->Draw("Psame");
grtr03to6->Draw("Psame");
grtr03to6fin->Draw("Psame");
grtr03to6fin4->Draw("Psame");
//grtheta10tr03to6->Draw("Psame");
//grtheta15tr03to6->Draw("Psame");
//grtr03to6d3->Draw("Psame");
//grtr03to3->Draw("Psame");
//grtr01to10->Draw("Psame");
//gr1Sum->Draw("Psame");
//gr11->Draw("Psame");
//gr12->Draw("Psame");
//gr13->Draw("Psame");
//gr1f->Draw("Psame");
//gr2->Draw("Psame");
//gr2f->Draw("Psame");
//gr3->Draw("Psame");
//grpfon->Draw("Psame");
//grpf->Draw("Psame");
//grpf03to6->Draw("Psame");
//grpf03to3->Draw("Psame");
//grpf01to10->Draw("Psame");
leg->Draw("same");
c1->Print("V2vsNtrk.png");
}

TGraphErrors* plot(TString dir,int isSum, int color, int marker,int debug){
const int ntheta=5;
const int ntotbin=5;
const int nbin[ntotbin] = {6,7,7,8,8};
const int nbinadd[ntotbin] = {0,6,13,20,28};
const double spcorr[ntotbin]={0.938962,0.997107,0.996962,0.995796,1.00314};
TString tbin[ntotbin] = {"M150120","M185150","M220185","M260220","M300260"};
int xbin=0;
int xpt=0;
double Nevent[ntotbin], Ntrk[ntotbin], Mult[ntotbin];
double Neventfin[36], Ntrkfin[36], V2mean_0fin[36], V2meanerr_0fin[36];
double  V2mean_0[ntotbin], V2meanerr_0[ntotbin];
double  V2_1[ntotbin], V2err_1[ntotbin];
double  V2_2[ntotbin], V2err_2[ntotbin];
double  V2_3[ntotbin], V2err_3[ntotbin];
double  V2_0[ntotbin][ntheta];
TVectorD *vecV2_0;
for(int i=0;i<ntotbin;i++){
if(isSum)TFile *mergedV = TFile::Open(Form("%s/%s/mergedV_Sum.root",dir.Data(),tbin[i].Data()));
else TFile *mergedV = TFile::Open(Form("%s/%s/mergedV_Prod.root",dir.Data(),tbin[i].Data()));
if(isSum)TFile *mergedv = TFile::Open(Form("%s/%s/mergedv_Prod.root",dir.Data(),tbin[i].Data()));
else TFile *mergedv = TFile::Open(Form("%s/%s/mergedv_Prod2.root",dir.Data(),tbin[i].Data()));
Nevent[i]=0;	Ntrk[i]=0; Mult[i]=0; 
for(int itheta=0;itheta<ntheta;itheta++)
V2_0[i][itheta]=0;        
V2mean_0[i]=0;
V2meanerr_0[i]=0;
if(debug!=0){
TVectorD *vecV2_2=(TVectorD*)mergedv->Get(Form("V_intcorr"));
TVectorD *vecV2err_2=(TVectorD*)mergedv->Get(Form("V_intcorrerr"));
TVectorD *vecV2_1=(TVectorD*)mergedv->Get(Form("V_int"));
TVectorD *vecV2err_1=(TVectorD*)mergedv->Get(Form("V_interr"));
}
if(dir.Contains("final")){
for(int ibin=0;ibin<nbin[i];ibin++){
TVectorD *vecNtrk = (TVectorD*)mergedV->Get("tottrk");
TVectorD *vecMult = (TVectorD*)mergedV->Get("totmultall");
TVectorD *vecNevent = (TVectorD*)mergedV->Get("Nevent");
if(!dir.Contains("tracknormcpt03to6"))
TVectorD *vecV2mean_0=(TVectorD*)mergedV->Get(Form("D_%d/Vmean",ibin));
else
TVectorD *vecV2mean_0=(TVectorD*)mergedV->Get(Form("Vmean",ibin));
if(!dir.Contains("tracknormcpt03to6"))
TVectorD *vecV2meanerr_0=(TVectorD*)mergedV->Get(Form("D_%d/deltaVmean",ibin));
else
TVectorD *vecV2meanerr_0=(TVectorD*)mergedV->Get(Form("deltaVmean"));
Ntrk[i]+=(*vecNtrk)[ibin];
Mult[i]+=(*vecMult)[ibin];
Nevent[i]+=(*vecNevent)[ibin];
V2mean_0[i]+=(*vecV2mean_0)[ibin]*(*vecNevent)[ibin];
V2meanerr_0[i]+=(*vecV2meanerr_0)[ibin]*(*vecNevent)[ibin];
Ntrkfin[nbinadd[i]+ibin]=(*vecNtrk)[ibin];
Ntrkfin[nbinadd[i]+ibin]/=(*vecNevent)[ibin];
V2mean_0fin[nbinadd[i]+ibin]=(*vecV2mean_0)[ibin];
V2meanerr_0fin[nbinadd[i]+ibin]=(*vecV2meanerr_0)[ibin];
}
Ntrk[i]/=Nevent[i]; 	Mult[i]/=Nevent[i];    
V2mean_0[i]/=Nevent[i];	V2meanerr_0[i]/=Nevent[i];
}
else{
TVectorD *vecNtrk = (TVectorD*)mergedV->Get("avgtrk");
TVectorD *vecavgmult = (TVectorD*)mergedV->Get("avgmultall");
TVectorD *vecNevent = (TVectorD*)mergedV->Get("Nevent");
if(!dir.Contains("tracknormcpt03to6"))
TVectorD *vecV2mean_0=(TVectorD*)mergedV->Get(Form("D_%d/Vmean",xbin));
else
TVectorD *vecV2mean_0=(TVectorD*)mergedV->Get(Form("Vmean",xbin));
//vecV2_0=(TVectorD*)mergedV->Get(Form("D_%d/D_%d/V",ibin,xpt));
for(int itheta=0;itheta<ntheta;itheta++)
//V2_0[i][itheta]=(*vecV2_0)[itheta];
if(!dir.Contains("tracknormcpt03to6"))
TVectorD *vecV2meanerr_0=(TVectorD*)mergedV->Get(Form("D_%d/deltaVmean",xbin));
else
TVectorD *vecV2meanerr_0=(TVectorD*)mergedV->Get(Form("deltaVmean",xbin));
Ntrk[i]=(*vecNtrk)[xbin];
Mult[i]=(*vecavgmult)[xbin];
Nevent[i]=(*vecNevent)[ibin];
V2mean_0[i]=(*vecV2mean_0)[xpt];
V2meanerr_0[i]=(*vecV2meanerr_0)[xpt];
}
if(debug!=0){
V2_1[i]=(*vecV2_1)[0];
V2err_1[i]=(*vecV2err_1)[0];
V2_2[i]=(*vecV2_2)[0];
V2err_2[i]=(*vecV2err_2)[0];
V2_3[i]=(*vecV2_2)[0]*spcorr[i];
V2err_3[i]=(*vecV2err_2)[0]*spcorr[i];
//fout<<Ntrk[i]<<"\t"<<V2_0[i]<<"\t"<<V2_1[i]<<"\t"<<V2_2[i]<<"\t"<<V2_3[i]<<endl;
//fout<<"\t"<<V2meanerr_0[i]<<"\t"<<V2err_1[i]<<"\t"<<V2err_2[i]<<"\t"<<V2err_3[i]<<endl;
}
//cout<<Nevent[i]<<"\t"<<tbin[i]<<"\t"<<dir<<"\t"<<Mult[i]<<endl;
for(int itheta=0;itheta<ntheta;itheta++){
fout<<itheta<<"\t"<<V2_0[i][itheta]<<"\t";
}
fout<<endl;
}
fout<<endl;
//const double V2[nbin]={0.0465,0.0498,0.0430,0.0447,0.0462};
//const double V2[nbin]={0.0536,0.0514,0.0491,0.0482,0.0462};
//const double V2err[nbin]={0.0013,0.0015,0.0015,0.0070,0.0027};
if(debug==0)TGraphErrors *gr=new TGraphErrors(ntotbin,Ntrk,V2mean_0,0,V2meanerr_0);
if(debug==1)TGraphErrors *gr=new TGraphErrors(ntotbin,Ntrk,V2_1,0,V2err_1);
if(debug==2)TGraphErrors *gr=new TGraphErrors(ntotbin,Ntrk,V2_2,0,V2err_2);
if(debug==3) TGraphErrors *gr=new TGraphErrors(ntotbin,Ntrk,V2_3,0,0);//modify
if(debug==4)TGraphErrors *gr=new TGraphErrors(36,Ntrkfin,V2mean_0fin,0,V2meanerr_0fin);
gr->SetMarkerColor(color);
gr->SetMarkerSize(1.3);
gr->SetLineColor(color);
gr->SetMarkerStyle(marker);
return gr;

}

