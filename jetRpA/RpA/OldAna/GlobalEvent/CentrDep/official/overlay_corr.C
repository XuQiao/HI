#include "../Quality/root_setting.h"

void overlay_corr(){
gStyle->SetOptStat(0);
gStyle->SetErrorX(0);
double HFbin[]={0,0.1,0.2,0.3,0.5,0.7,1,1.5,2,2.5,3,4,5,6,7,8,10,12,15,20,30};
/*bool isQiao=kTRUE;
if(isQiao){
double HFbin[]={0,0.82,4.62,11.12,14.52,19.80,100};}
else{
double HFbin[]={0,0.92,4.72,11.21,14.60,19.88,100};}
*/
int nHFbin=sizeof(HFbin)/sizeof(double)-1;


TFile *f=TFile::Open("output/weight/checkfun_Hijing_DS_bakpointAll.root");
TString histoname[3]={"HFEnergy4","HFEnergyPlus","HFEnergyPlus4"};
TH1F* fDS[3]; TH1F* fcorr[3]; TH1F* fDS_Cl[3];
c1 = new TCanvas("c1"," ",500,700);
makeMultiPanelCanvas(c1,1,2,-0.05,0,0.17,0.15,0.03);
for(int i=0;i<3;i++){
fDS[i]=(TH1F*)f->Get(Form("h%s_DS",histoname[i].Data()));
fcorr[i]=(TH1F*)f->Get(Form("h%s_corr",histoname[i].Data()));
fDS[i]=(TH1F*)fDS[i]->Rebin(nHFbin,Form("DS_%s",fDS[i]->GetName()),HFbin);
fcorr[i]=(TH1F*)fcorr[i]->Rebin(nHFbin,Form("corr_%s",fcorr[i]->GetName()),HFbin);
fixedFontHist(fDS[i],2.1,2.5);
fixedFontHist(fcorr[i],2.1,2.5);
normalizeByBinWidth(fDS[i]);
normalizeByBinWidth(fcorr[i]);
//fDS[i]->Scale(1./fDS[i]->Integral());
//fcorr[i]->Scale(1./fcorr[i]->Integral());
fDS_Cl[i]=(TH1F*)fDS[i]->Clone(Form("ratio_%s",fDS[i]->GetName()));
fDS_Cl[i]->Divide(fcorr[i]);
fDS[i]->SetMarkerSize(1.2);
fDS[i]->SetMarkerColor(1);
fDS[i]->SetLineColor(1);
fDS[i]->SetMarkerStyle(20);

fDS_Cl[i]->SetMarkerSize(1.2);
fDS_Cl[i]->SetMarkerColor(1);
fDS_Cl[i]->SetLineColor(1);
fDS_Cl[i]->SetMarkerStyle(20);
fDS_Cl[i]->GetYaxis()->SetRangeUser(0.8,1.19);

fcorr[i]->GetXaxis()->SetTitle("");
fcorr[i]->GetYaxis()->SetTitle("# of Events");
fDS_Cl[i]->GetXaxis()->SetTitle("HF #Sigma E_{T} (GeV/c) #eta>4");
fDS_Cl[i]->GetYaxis()->SetTitle("DS/corr");

fcorr[i]->SetMarkerSize(1.2);
fcorr[i]->SetMarkerColor(2);
fcorr[i]->SetLineColor(2);
fcorr[i]->SetMarkerStyle(24);
}

i=2;
TLegend *leg=new TLegend(0.2,0.2,0.45,0.4);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.05);
leg->AddEntry(fDS[i],"DS","lp");
leg->AddEntry(fcorr[i],"corr","lp");
c1->cd(1)->SetLogy();
//c1->cd(1)->SetLogx();
fcorr[i]->Draw();
cout<<fcorr[i]->Integral()<<endl;
fDS[i]->Draw("same");
cout<<fDS[i]->Integral()<<endl;
leg->Draw("same");
/*if(isQiao)
TLatex *T=new TLatex(0.4,0.9,"Use Qiao's Centrality HFbinning");
else*/
TLatex *T=new TLatex(0.4,0.9,"");
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.05);
T->SetTextColor(1);
T->SetTextFont(42);
T->Draw("same");

c1->cd(2);
//->SetLogx();
fDS_Cl[i]->Draw();
TLine *l=new TLine(HFbin[0],1,HFbin[nHFbin],1);
l->SetLineStyle(2);
l->Draw("same");
/*if(isQiao)
c1->Print("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/CentrDep/pic/CompareWt_Qiao.pdf");
else*/
//c1->Print("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/CentrDep/pic/CompareWt.pdf");
}

