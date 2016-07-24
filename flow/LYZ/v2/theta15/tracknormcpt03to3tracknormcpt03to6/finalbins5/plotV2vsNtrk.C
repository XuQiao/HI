void plotV2vsNtrk(){

const int ntotbin=5;
const int nbin[ntotbin] = {6,7,7,8,8};
const int Nbin = 36;
const int nbin24 = 12;
const double avgtrkbin[nbin24]={44.36,54.33,68.63,88.39,108.2,131.3,162.1,196.6,227.5,247.2,269.2,301.2};
const double V24[nbin24]={0.02965,0.03913,0.04832,0.04941,0.04822,0.04955,0.049,0.04805,0.04709,0.04665,0.04772,0.04797};
const double V24err[nbin24]={0.005967,0.004262,0.001496,0.001351,0.001599,0.0003922,0.0003065,0.0002939,0.0004568,0.0008684,0.001471,0.004329};
TString tbin[ntotbin] = {"M150120","M185150","M220185","M260220","M300260"};
int xpt=0;
int count=0;
double Ntrk[Nbin], V2_Sum[Nbin], V2err_Sum[Nbin],  V2_Prod[Nbin], V2err_Prod[Nbin];
double Ntrkavg[ntotbin], Nevent[ntotbin], Nmult[ntotbin], V2avg_Sum[ntotbin], V2erravg_Sum[ntotbin],  V2avg_Prod[ntotbin], V2erravg_Prod[ntotbin];
for(int i=0;i<ntotbin;i++){
TFile *mergedV_Sum = TFile::Open(Form("%s/mergedV_Sum.root",tbin[i].Data()));
TFile *mergedV_Prod = TFile::Open(Form("%s/mergedV_Prod.root",tbin[i].Data()));
TVectorD *vecNtrk = (TVectorD*)mergedV_Sum->Get("tottrk");
TVectorD *vecNtrk_back = (TVectorD*)mergedV_Prod->Get("tottrk");
TVectorD *vectotmult = (TVectorD*)mergedV_Prod->Get("totmultall");
TVectorD *vecNevent = (TVectorD*)mergedV_Sum->Get("Nevent");
Ntrkavg[i]=0; Nevent[i]=0;	Nmult[i]=0; V2avg_Sum[i]=0;	V2erravg_Sum[i]=0;	V2avg_Prod[i]=0;	V2erravg_Prod[i]=0;
for(int ibin=0;ibin<nbin[i];ibin++){
count++;
TVectorD *vecV2_Sum=(TVectorD*)mergedV_Sum->Get(Form("D_%d/Vmean",ibin));
TVectorD *vecV2err_Sum=(TVectorD*)mergedV_Sum->Get(Form("D_%d/deltaVmean",ibin));
TVectorD *vecV2_Prod=(TVectorD*)mergedV_Prod->Get(Form("D_%d/Vmean",ibin));
TVectorD *vecV2err_Prod=(TVectorD*)mergedV_Prod->Get(Form("D_%d/deltaVmean",ibin));
Ntrk[count-1]=(*vecNtrk)[ibin]/(*vecNevent)[ibin];
Ntrkavg[i]+=(*vecNtrk)[ibin];
Nevent[i]+=(*vecNevent)[ibin];
Nmult[i]+=(*vectotmult)[ibin];
V2_Sum[count-1]=(*vecV2_Sum)[xpt];
cout<<(*vectotmult)[ibin]<<endl;
V2avg_Sum[i]+=(*vecV2_Sum)[xpt]*(*vectotmult)[ibin];
V2err_Sum[count-1]=(*vecV2err_Sum)[xpt];
V2erravg_Sum[i]+=(*vecV2err_Sum)[xpt]*(*vectotmult)[ibin];
V2_Prod[count-1]=(*vecV2_Prod)[xpt];
V2avg_Prod[i]+=(*vecV2_Prod)[xpt]*(*vectotmult)[ibin];
V2err_Prod[count-1]=(*vecV2err_Prod)[xpt];
V2erravg_Prod[i]+=(*vecV2err_Prod)[xpt]*(*vectotmult)[ibin];
}
Ntrkavg[i]/=Nevent[i];	V2avg_Sum[i]/=Nmult[i];	V2erravg_Sum[i]/=Nmult[i];	V2avg_Prod[i]/=Nmult[i];	V2erravg_Prod[i]/=Nmult[i];
}
//const double V2[nbin]={0.0465,0.0498,0.0430,0.0447,0.0462};
//const double V2Prod[nbin]={0.0536,0.0514,0.0491,0.0482,0.0462};
//const double V2Proderr[nbin]={0.0013,0.0015,0.0015,0.0070,0.0027};
TCanvas *c1 = new TCanvas;
TGraphErrors *grSum=new TGraphErrors(Nbin,Ntrk,V2_Sum,0,V2err_Sum);
TGraphErrors *gravgSum=new TGraphErrors(ntotbin,Ntrkavg,V2avg_Sum,0,V2erravg_Sum);
TGraphErrors *grProd=new TGraphErrors(Nbin,Ntrk,V2_Prod,0,V2err_Prod);
TGraphErrors *gravgProd=new TGraphErrors(ntotbin,Ntrkavg,V2avg_Prod,0,V2erravg_Prod);
TGraphErrors *grv24=new TGraphErrors(nbin24,avgtrkbin,V24,0,V24err);
TLatex *tl1 = new TLatex(0.2,0.3,Form("track normal cut"));
TLatex *tl2 = new TLatex(0.2,0.2,Form("%.1f < p_{T} < %.1f (GeV/c)",0.3,6.0));
tl1->SetNDC();
tl2->SetNDC();
grv24->SetTitle("V_{2} vs Ntrkoffline");
grv24->GetXaxis()->SetTitle("Ntrkoffline");
grv24->GetYaxis()->SetTitle("V_{2}");
grv24->SetMaximum(0.08);
grv24->SetMinimum(0.0);
grSum->SetMarkerSize(1);
gravgSum->SetMarkerSize(1);
grv24->SetMarkerSize(1);
grSum->SetMarkerColor(1);
gravgSum->SetMarkerColor(1);
grProd->SetMarkerColor(3);
grProd->SetLineColor(3);
gravgProd->SetMarkerColor(4);
gravgProd->SetLineColor(4);
grSum->SetMarkerStyle(20);
gravgSum->SetMarkerStyle(24);
grProd->SetMarkerStyle(29);
gravgProd->SetMarkerStyle(30);
grv24->SetMarkerStyle(20);
grv24->SetMarkerColor(2);
grv24->SetLineColor(2);
TLegend *leg = new TLegend(0.2,0.7,0.4,0.85);
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
//leg->AddEntry(grSum,"LYZ Sum method final bin","P");
//leg->AddEntry(gravgSum,"LYZ Sum method average","P");
leg->AddEntry(grProd,"LYZ Prod method final bin","P");
leg->AddEntry(gravgProd,"LYZ Prod method average","P");
leg->AddEntry(grv24,"4-particle cumulant","P");
grv24->Draw("AP");
//grSum->Draw("Psame");
//gravgSum->Draw("Psame");
grProd->Draw("Psame");
gravgProd->Draw("Psame");
tl1->Draw("same");
tl2->Draw("same");
leg->Draw("same");
c1->Print("V2vsNtrk.png");

}

