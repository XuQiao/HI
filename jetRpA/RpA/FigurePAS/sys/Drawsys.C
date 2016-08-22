#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void Drawsys(){
	static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};
    TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.11);
    c1->SetFrameFillStyle(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
const int N=8;
const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12","-10_10"};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
const int Nbin_pt=30;
const int Nbin_ptpre=16;
ifstream fstr1[N];
ifstream fstr2[N];
ifstream fstr3[N];
TGraph *g[3][N];
for(int i=0;i<N;i++){
fstr1[i].open(Form("jetsysIter%s.txt",etabinname[i].Data()));
fstr2[i].open(Form("jetaxissys%s.txt",etabinname[i].Data()));
fstr3[i].open(Form("jetUnfsys%s.txt",etabinname[i].Data()));
}

double sys1[Nbin_ptpre],sys2[Nbin_ptpre],sys3[Nbin_ptpre];
double pt[Nbin_ptpre];
int i=7;
for(int j=0;j<Nbin_ptpre;j++){
fstr1[i]>>pt[j];
fstr2[i]>>pt[j];
fstr3[i]>>pt[j];
fstr1[i]>>sys1[j];
fstr2[i]>>sys2[j];
fstr3[i]>>sys3[j];
}
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(50,600);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetRangeUser(0,30);
hFrame->GetYaxis()->SetTitle("systematics(%)");
fixedFontHist(hFrame,1.1,1.4);
hFrame->DrawCopy();
TLegend* leg1=new TLegend(0.6,0.75,0.78,0.92);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.03);
g[0][i]=new TGraph(Nbin_ptpre,pt,sys1);	
g[0][i]->SetMarkerColor(colorCode[0]);	g[0][i]->SetLineColor(colorCode[0]);	
g[0][i]->SetMarkerStyle(markerCode[0]); g[0][i]->SetMarkerSize(1.2);
g[0][i]->Draw("Psame");
g[1][i]=new TGraph(Nbin_ptpre,pt,sys2);
g[1][i]->SetMarkerColor(colorCode[1]); 	g[1][i]->SetLineColor(colorCode[1]);   
g[1][i]->SetMarkerStyle(markerCode[1]); g[1][i]->SetMarkerSize(1.2);
g[1][i]->Draw("Psame");
g[2][i]=new TGraph(Nbin_ptpre,pt,sys3);	
g[2][i]->SetMarkerColor(colorCode[2]); g[2][i]->SetLineColor(colorCode[2]); 
g[2][i]->SetMarkerStyle(markerCode[2]);	g[2][i]->SetMarkerSize(1.2);
g[2][i]->Draw("Psame");
leg1->AddEntry(g[0][i],"Iteration","lp");
leg1->AddEntry(g[1][i],"Jet pointing","lp");
leg1->AddEntry(g[2][i],"Unfolding Closure","lp");
leg1->Draw("same");
TLatex *com3 = new TLatex(0.45,0.24,etastring[i]);
com3->SetTextFont(63);
com3->SetTextSize(17);
com3->SetNDC();
com3->Draw();

ofstream fstrout[N];
//for(int i=0;i<N;i++)
fstrout[i].open(Form("jetsys_%s.txt",etabinname[i].Data()));
for(int j=0;j<Nbin_ptpre;j++){
fstrout[i]<<pt[j]<<'\t';
if(pt[j]==51)
fstrout[i]<<sys1[j];
else if(pt[j]<200)
fstrout[i]<<"4";
else if(pt[j]<300)
fstrout[i]<<"3";
else
fstrout[i]<<"0";
fstrout[i]<<endl;
}
}
