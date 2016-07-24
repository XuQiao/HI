#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void ProduceEff(){
gStyle->SetOptStat(kFALSE);
TString Var[3]={"HFEnergy4","HFEnergyPlus","HFEnergyPlus4"};
ofstream fstr("Eff_table.txt");
TFile *f=TFile::Open("output/weight/pPbHist_Hijing_DS.root");
TFile *fEff=TFile::Open("output/weight/pPbHijing_DS_EffCorr.root","RECREATE");
fEff->cd();
c1 = new TCanvas("c1"," ",600,600);
makeMultiPanelCanvas(c1,1,1,0,0,0.05,0.05,0.03);
c1->cd();
c1->SetLogx();
TH1D * hFrame = new TH1D("","",100,0,100.);
hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Efficiency");
hFrame->SetTitle("");
hFrame->SetMinimum(0);
hFrame->SetMaximum(1.1);
fixedFontHist(hFrame,1,1);
T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.05);
T->SetTextColor(1);
T->SetTextFont(42);
double HFbin[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,10,11,12,14,16,20,25,30,40,50};
int nHFbin=sizeof(HFbin)/sizeof(double)-1;
       hFrame->GetXaxis()->SetLimits(HFbin[0],HFbin[nHFbin]);
TLine *l=new TLine(HFbin[0],1,HFbin[nHFbin],1);
l->SetLineStyle(2);
TH1F* hHF_tr[3]; TH1F* hHF_ds[3]; TH1F* hHF_eff[3]; TH1F* rehHF_tr[3]; TH1F* rehHF_ds[3];	TGraphAsymmErrors *gEff[3];
for(int i=0;i<3;i++){
hHF_tr[i]=(TH1F*)f->Get(Form("h%s_tr",Var[i].Data()));
hHF_ds[i]=(TH1F*)f->Get(Form("h%s_DS",Var[i].Data()));
rehHF_tr[i]=(TH1F*)hHF_tr[i]->Clone(hHF_tr[i]->GetName()); //Rebin(nHFbin,hHF_tr[i]->GetName(),HFbin);
rehHF_ds[i]=(TH1F*)hHF_ds[i]->Clone(hHF_ds[i]->GetName()); //Rebin(nHFbin,hHF_ds[i]->GetName(),HFbin);
gEff[i]=new TGraphAsymmErrors();	gEff[i]->SetName(Form("gEff_%s",Var[i].Data()));
gEff[i]->BayesDivide(rehHF_tr[i],rehHF_ds[i]);
gEff[i]->SetMarkerSize(1.2);
gEff[i]->SetMarkerStyle(20);
hHF_eff[i]=(TH1F*)rehHF_tr[i]->Clone(Form("rehHF_eff_%s",Var[i].Data()));
hHF_eff[i]->Divide(rehHF_ds[i]);
hHF_eff[i]->Write();
gEff[i]->Write();
fixedFontHist(hHF_eff[i],1,1);
hFrame->Draw();
gEff[i]->Draw("Psame");
l->Draw("same");
T->DrawLatex(0.2,0.85,Var[i]);
c1->Print(Form("pic/Efficiency_%s.png",Var[i].Data()));
fstr<<Var[i]<<endl;
fstr<<"bin"<<'\t'<<"HF range"<<'\t'<<"# before"<<'\t'<<"# after"<<'\t'<<"Effciency"<<'\t'<<"Error"<<endl;
/*for(int nbin=1;nbin<hHF_eff[i]->GetNbinsX();nbin++){
if(rehHF_tr[i]->GetBinContent(nbin))
//fstr<<nbin<<'\t'<<HFbin[nbin-1]<<'-'<<HFbin[nbin]<<'\t'<<rehHF_ds[i]->GetBinContent(nbin)<<'\t'<<rehHF_tr[i]->GetBinContent(nbin)<<'\t'<<hHF_eff[i]->GetBinContent(nbin)<<'\t'<<hHF_eff[i]->GetBinError(nbin)<<endl;
}*/
}
}
