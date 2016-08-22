#include "/afs/cern.ch/user/q/qixu/CMSSW_5_3_8_HI/src/jetRpA/RpA/Quality/root_setting.h"

void ProduceEff(){
gStyle->SetOptStat(kFALSE);
TString Isweight="weight";
TString Gen[2]={"Hijing","Epos"};
int k=0;
TString Var[4]={"HFEnergy4","HFEnergyPlus","HFEnergyPlus4","Ntrk"};
TString etarange[4]={" HF E_{T} |#eta|>4"," HF E_{T} #eta>0"," HF E_{T} #eta>4",""};
ofstream fstr(Form("Efftable_%s_%s.txt",Gen[k].Data(),Isweight.Data()));

TFile *f=TFile::Open(Form("output/%s/pPbHist_%s_DS.root",Isweight.Data(),Gen[k].Data()));
TFile *fEff=TFile::Open(Form("output/%s/pPb%s_DS_EffCorr.root",Isweight.Data(),Gen[k].Data()),"RECREATE");
fEff->cd();
c1 = new TCanvas("c1"," ",600,600);
makeMultiPanelCanvas(c1,1,1,0,0,0.05,0.05,0.03);
c1->cd();
//c1->SetLogx();
TH1D * hFrame = new TH1D("","",10000,0,100.);

hFrame->GetYaxis()->SetTitle("Efficiency");
hFrame->SetTitle("");
hFrame->SetMinimum(0);
hFrame->SetMaximum(1.1);
hFrame->GetXaxis()->SetRangeUser(0,50);
fixedFontHist(hFrame,1,1);
T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.05);
T->SetTextColor(1);
T->SetTextFont(42);
//double HFbin[]={0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,10,11,12,14,16,20,25,30,40,50};
//double HFbin[]={0,0.5,1,2,3,4,6,8,10,15,20,30,40,100};
double HFbin[]={0,0.24,0.47,0.67,0.86,1.05,1.24,1.43,1.62,1.83,2.03,2.24,2.47,2.7,2.95,3.21,3.48,3.76,4.06,4.38,4.71,5.05,5.42,5.79,6.18,6.58,6.99,7.41,7.84,8.3,8.76,9.24,9.73,10.24,10.77,11.31,11.9,12.49,13.12,13.79,14.5,15.26,16.08,16.97,17.96,19.11,20.44,22.12,24.34,27.91,74.73};
int nHFbin=sizeof(HFbin)/sizeof(double)-1;
//hFrame->GetXaxis()->SetLimits(HFbin[0],hFrame->GetXaxis()->FindBin(20));
TLine *l=new TLine(HFbin[0],1,20,1);
l->SetLineStyle(2);
TH1F* hHF_tr[4]; TH1F* hHF_ds[4]; TH1F* hHF_eff[4]; TH1F* rehHF_tr[4]; TH1F* rehHF_ds[4];TH1F* rehHF_eff[4];//	TGraphAsymmErrors *gEff[4];
TGraph* gEff[4];
TF1 *fitEff[4];
double chi2,ndf,par[5],err[5];
int i0=3;
if(i0!=3)
hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
else
hFrame->GetXaxis()->SetTitle("track Mult");
hFrame->DrawCopy();
TLegend *leg = new TLegend(0.55,0.2,0.75,0.4);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.04);
for(int i=0;i<4;i++){
hHF_tr[i]=(TH1F*)f->Get(Form("h%s_tr",Var[i].Data()));
hHF_ds[i]=(TH1F*)f->Get(Form("h%s_DS",Var[i].Data()));
hHF_eff[i]=(TH1F*)hHF_tr[i]->Clone(Form("hEff_%s",Var[i].Data()));
hHF_eff[i]->Divide(hHF_ds[i]);
fixedFontHist(hHF_eff[i],1,1);
if(i!=3){
rehHF_tr[i]=(TH1F*)hHF_tr[i]->Rebin(nHFbin,hHF_tr[i]->GetName(),HFbin);
rehHF_ds[i]=(TH1F*)hHF_ds[i]->Rebin(nHFbin,hHF_ds[i]->GetName(),HFbin);
rehHF_eff[i]=(TH1F*)rehHF_tr[i]->Clone(Form("rehEff_%s",Var[i].Data()));
rehHF_eff[i]->Divide(rehHF_ds[i]);
gEff[i]=new TGraph(rehHF_eff[i]);
}
else{
gEff[i]=new TGraph(hHF_eff[i]);
//gff[i]->SetPoint(0,0.5,1);
}
gEff[i]->SetName(Form("gEff_%s",Var[i].Data()));
//gEff[i]->BayesDivide(rehHF_tr[i],rehHF_ds[i]);
gEff[i]->SetMarkerSize(1.2);
gEff[i]->SetMarkerStyle(20);
gEff[i]->SetMarkerColor(1);
//gEff[i]->GetXaxis()->SetLimits(0,100);
if(i!=3){
fitEff[i]=new TF1(Form("fitEff_%s",Var[i].Data()),"[0]-[1]*exp([2]*x)-[3]*exp(-[4]*x)",HFbin[0],HFbin[nHFbin]);
if(k==0){
fitEff[i]->SetParameters(1.,0.3,-2,0.3,0.6);
}
else if(k==1){
if(i==1)
fitEff[i]->FixParameter(1,0.5);
fitEff[i]->SetParameters(1.,0.5,-2,0.5,0.5);
}
//gEff[i]->Fit(Form("fitEff_%s",Var[i].Data()),"R");
chi2 = fitEff[i]->GetChisquare();
ndf = fitEff[i]->GetNDF();
fstr<<"chi2="<<chi2<<'\t'<<"NDF="<<ndf<<endl;
for(int j=0;j<5;j++){
par[j]=fitEff[i]->GetParameter(j);
err[j]=fitEff[i]->GetParError(j);
fstr<<"par"<<j<<"="<<par[j]<<'\t'<<"err"<<j<<"="<<err[j]<<endl;
}
fitEff[i]->SetLineColor(2+i);
fitEff[i]->SetLineStyle(2+i);
fitEff[i]->Write();
rehHF_eff[i]->Write();
}
if(i==i0){
if(i0!=3)
leg->AddEntry(gEff[i],etarange[i].Data(),"lp");
gEff[i]->Draw("Psame");
}
//fitEff[i]->Draw("same");
gEff[i]->Write();
fstr<<Var[i]<<endl;
hHF_eff[i]->Write();
}
l->Draw("same");
//T->DrawLatex(0.2,0.85,Var[i]);
leg->Draw("same");
//c1->Print(Form("pic/Efficiency_%s.png",Var[i].Data()));
c1->Print(Form("pic/Efficiency_%s_%s.png",Gen[k].Data(),Var[i0].Data()));

}
