#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void DrawVzcut(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
const int NColor=7;
static const int Color[NColor]={
   1, 2, 4, 46,6, 7,8
};
static const int Style[NColor] = {
    20, 34, 33, 25,27, 28//,24
};

const double pPbLumi = 20.7 ;

const int Ntrig = 5;const int NVz=4;
TString TrigName[Ntrig]={"Jet20","Jet40","Jet60","Jet80","Jet100"};
TFile *f[Ntrig];

int Vzname[NVz]={2,3,4,7};	
TString Vzstring[NVz]={"HLT","Beam","pprim","EvtSel"};
TH1F* hVz[Ntrig][NVz];	TH1F* hRatio1[NVz];
int trig0=4;
TH1F* hFrame=new TH1F("","",1000,-50,50);
hFrame->GetXaxis()->SetRangeUser(-30,30);
hFrame->GetYaxis()->SetTitle("Number of Events");
fixedFontHist(hFrame,1.4,1.6);
TLegend *leg=new TLegend(0.60,0.75,0.8,0.95);
TLegend *leg_=new TLegend(0.60,0.75,0.8,0.95);
//TLine *l=new TLine(-30,1,30,1);
TLine *l=new TLine(-15,1,15,1);
l->SetLineStyle(2);
l->SetLineColor(1);

TLine *l1=new TLine(-15,1.01e-6,-15,2e-1);
TLine *l3=new TLine(-15,0,-15,2);
TLine *l2=new TLine(15,1.01e-6,15,2e-1);
TLine *l4=new TLine(15,0,15,2);
l1->SetLineStyle(2);
l2->SetLineStyle(2);
l3->SetLineStyle(2);
l4->SetLineStyle(2);
l1->SetLineColor(2);
l2->SetLineColor(2);
l3->SetLineColor(2);
l4->SetLineColor(2);
l1->SetLineWidth(4);
l2->SetLineWidth(4);
l3->SetLineWidth(4);
l4->SetLineWidth(4);

leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextFont(42);
leg->SetTextSize(0.035);
/*if(k0==0){
leg->SetHeader("Before Event Selection");
leg_->SetHeader("Before Event Selection");
}
lse{
leg->SetHeader("After Event Selection");
leg_->SetHeader("After Event Selection");
}*/
//leg->SetTextAlign(32);
leg_->SetBorderSize(0);
leg_->SetFillColor(0);
leg_->SetTextFont(42);
leg_->SetTextSize(0.035);
for(int itrig=4;itrig<5;itrig++){
for(int iVz=0;iVz<NVz;iVz++){
f[itrig]=TFile::Open(Form("/scratch/xuq7/RpA/GlobalEvent/EventSel/%sEvSel.root",TrigName[itrig].Data()));
hVz[itrig][iVz]=(TH1F*)f[itrig]->Get(Form("hVz_cut%d",Vzname[iVz]));
hVz[itrig][iVz]->Rebin(10);
//normalizeByBinWidth(hVz[trig0][iVz]);
//hVz[trig0][iVz]->Scale((double)1.0/hVz[trig0][iVz]->Integral());
}
}
c1 = new TCanvas("c1"," ",600,600);
c2 = new TCanvas("c2"," ",600,600);
makeMultiPanelCanvas(c1,1,1,-0.05,0,0.12,0.12,0.03);
makeMultiPanelCanvas(c2,1,1,-0.05,0,0.12,0.12,0.03);

c1->cd(1)->SetLogy();
hFrame->GetXaxis()->SetTitle("Vz (cm)");
hFrame->GetYaxis()->SetRangeUser(1,1e8);
hFrame->DrawCopy();
for(iVz=0;iVz<NVz;iVz++){
hVz[trig0][iVz]->SetMarkerColor(Color[iVz]);
hVz[trig0][iVz]->SetLineColor(Color[iVz]);
hVz[trig0][iVz]->SetMarkerStyle(Style[iVz]);
hVz[trig0][iVz]->SetMarkerSize(1.2);
leg->AddEntry(hVz[trig0][iVz],Vzstring[iVz],"lp");
hVz[trig0][iVz]->Draw("same");
}
leg->Draw("same");
drawCMS(0.2,0.9,pPbLumi);
c1->Print("Vzcut.pdf");
//l1->Draw("same");
//l2->Draw("same");
c2->cd(1);

hFrame->GetXaxis()->SetTitle("Vz (cm)");
//hFrame->GetYaxis()->SetRangeUser(0.9,1.25);
hFrame->GetYaxis()->SetRangeUser(0.7,1.2);
hFrame->GetYaxis()->SetTitle("Ratio");
hFrame->DrawCopy();
int Selbase=0;
for(iVz=0;iVz<NVz;iVz++){
if(iVz!=Selbase){
hRatio1[iVz]=(TH1F*)hVz[trig0][iVz]->Clone(Form("hRatio1_%d",iVz));
hRatio1[iVz]->SetMarkerColor(Color[iVz]);
hRatio1[iVz]->SetLineColor(Color[iVz]);
hRatio1[iVz]->SetMarkerStyle(Style[iVz]);
hRatio1[iVz]->SetMarkerSize(1.2);
leg_->AddEntry(hRatio1[iVz],Form("%s/%s",Vzstring[iVz].Data(),Vzstring[Selbase].Data()),"lp");
hRatio1[iVz]->Divide(hVz[trig0][Selbase]);
hRatio1[iVz]->Draw("same");
}
}
//l3->Draw("same");
//l4->Draw("same");
leg_->Draw("same");
l->Draw("same");
drawCMS(0.2,0.9,pPbLumi);
c2->Print("Vzcut_ratio.pdf");

/*
int i0=5, k1=0, k2=1;
c2 = new TCanvas("c2"," ",600,800);
makeMultiPanelCanvas(c2,1,2,-0.05,0,0.12,0.12,0.03);
c2->cd(1)->SetLogy();
//hFrame->GetYaxis()->SetRangeUser(1e-5,2e-1);
hFrame->GetYaxis()->SetRangeUser(0.5,1e7);
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetTitle("Number of Event");
hVz[i0][k1]->SetMarkerColor(1);
hVz[i0][k1]->SetLineColor(1);
hVz[i0][k1]->SetMarkerStyle(20);
hVz[i0][k1]->SetMarkerSize(1.2);
hVz[i0][k2]->SetMarkerColor(2);
hVz[i0][k2]->SetLineColor(2);
hVz[i0][k2]->SetMarkerStyle(24);
hVz[i0][k2]->SetMarkerSize(1.2);
leg->Clear();
leg->SetHeader(TrigName[i0]);
leg->AddEntry(hVz[i0][k1],"Before ES","lp");
leg->AddEntry(hVz[i0][k2],"After ES","lp");
hFrame->DrawCopy();
hVz[i0][k1]->Draw("same");
hVz[i0][k2]->Draw("same");
leg->Draw("same");
c2->cd(2);
hFrame->GetYaxis()->SetRangeUser(0,1.99);
hFrame->GetXaxis()->SetTitle("Vz (cm)");
hFrame->GetYaxis()->SetRangeUser(0.1,1.89);
hFrame->GetYaxis()->SetTitle(Form("After/Before"));
hFrame->DrawCopy();
TH1F* hRatio2=(TH1F*)hVz[i0][k2]->Clone("hRatio2");
hRatio2->Divide(hVz[i0][k1]);
hRatio2->SetMarkerColor(1);
hRatio2->SetLineColor(1);
hRatio2->SetMarkerStyle(20);
hRatio2->Draw("same");
l->Draw("same");
c2->Print(Form("%sEvSel.png",TrigName[i0].Data()));
*/
}
