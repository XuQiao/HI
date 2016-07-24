#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

ReverseXAxis (TGraphAsymmErrors *g, double  size)
{
// Remove the current axis
   g->GetXaxis()->SetLabelOffset(999);
   g->GetXaxis()->SetTickLength(0);

//  Redraw the new axis 
   gPad->Update();

	g->GetYaxis()->SetTitle("Event Selection Efficiency");
   TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
                                gPad->GetUymin(),
                                gPad->GetUxmin(),
                                gPad->GetUymin(),
                                g->GetXaxis()->GetXmin(),
                                g->GetXaxis()->GetXmax(),
                                510,"-");
   newaxis->SetLabelOffset(-0.04);
   newaxis->SetTitle("Centrality");
   newaxis->SetTitleFont(42);
   newaxis->SetTitleSize(size);
   newaxis->SetLabelSize(size);
   newaxis->CenterTitle();
   newaxis->SetLabelFont(42);
   newaxis->Draw();


}

ReverseXGraph (TGraphAsymmErrors *g)
{
 //  Create a new graph
   Int_t n = g->GetN();
double x[100],xerrhigh[100],xerrlow[100],y[100],yerrhigh[100],yerrlow[100],xr[100];
const double cent[100]={5,15,25,45,75,95};
const double cent_bin[100]={0,10,20,30,60,90,100};
  double dx=g->GetXaxis()->GetXmin()+g->GetXaxis()->GetXmax();
        for(int i=0;i<n;i++){
   xr[i] = dx-cent[n-i-1];
   xerrhigh[i] = cent_bin[n-i]-cent[n-i-1];
   xerrlow[i] = cent[n-i-1]-cent_bin[n-i-1];
   y[i] = g->GetY()[i];
   yerrhigh[i] = g->GetErrorYhigh(i);
   yerrlow[i] = g->GetErrorYlow(i);
   }
  TGraphAsymmErrors*  gr = new TGraphAsymmErrors(n,xr,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
   gr->SetMarkerSize(g->GetMarkerSize());
   gr->SetMarkerStyle(g->GetMarkerStyle());
   gr->SetMarkerColor(g->GetMarkerColor());
   gr->Draw("P");
}

void Draw(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
bool isHF=kTRUE;
//bool isHF=kFALSE;
TString dir[2]={"weight","noweight"};
const int Nhist=3;	const int NGen=2; const int Ncent=6;
const double HFbin[]={0,0.1,0.2,0.4,0.6,0.8,1,2,4,6,9,12,15,20,30,50,75,100};
const double xerr[Ncent]={0,0,0,0,0,0};
const double HFbin_cent[Nhist][Ncent+1]={
{0,2.92,10.1,19.3,23.7,30.2,100},
{0,1.03,5.00,11.6,15.0,20.3,100},
{0,2.7,12.6,28.6,36.9,49.6,100}};
const double cent[Ncent]={95,75,45,25,15,5};
const int nHFbin=sizeof(HFbin)/sizeof(double)-1;
double color[NGen][Nhist]={{1,2,4},{1,2,4}};
double Style[NGen][Nhist]={{20,24,28},{20,24,28}};
TString Gen[NGen]={"Hijing","Epos"};
TString histoname[Nhist]={"hHFEnergy4","hHFEnergyPlus4","hHFEnergyPlus"};
TString etarange[Nhist]={" HF E_{T} |#eta|>4"," HF E_{T} #eta>4"," HF E_{T} #eta>0"};
TString filename[2][NGen];
TFile* f[2][NGen];
TH1F* h_or[2][NGen][Nhist];	TH1F* reh_or[2][NGen][Nhist];
TH1F* h_tr[2][NGen][Nhist];	TH1F* reh_tr[2][NGen][Nhist];
TH1F* h_ds[2][NGen][Nhist];	TH1F* reh_ds[2][NGen][Nhist];
TLatex* T1[Nhist];	TLatex* T2;	TLatex* T3[NGen];
TGraphAsymmErrors *gEff[2][NGen][Nhist];
c1 = new TCanvas("c1"," ",800,500);
makeMultiPanelCanvas(c1,2,1,0.03,0.03,0.14,0.14,0.03);

for(int k=0;k<NGen;k++){
T3[k]=new TLatex(0.45,0.92,Gen[k]);
T3[k]->SetNDC();
T3[k]->SetTextAlign(12);
T3[k]->SetTextSize(0.05);
T3[k]->SetTextColor(1);
T3[k]->SetTextFont(42);
}
for(int j=0;j<Nhist;j++){
T1[j]=new TLatex(0.4,0.85,etarange[j]);
T1[j]->SetNDC();
T1[j]->SetTextAlign(12);
T1[j]->SetTextSize(0.05);
T1[j]->SetTextColor(1);
T1[j]->SetTextFont(42);
}
T2=new TLatex(0.2,0.25,"Double Side");
T2->SetNDC();
T2->SetTextAlign(12);
T2->SetTextSize(0.06);
T2->SetTextColor(1);
T2->SetTextFont(42);


TH1D * hFrame = new TH1D("","",100,0,100.);
       hFrame->GetXaxis()->SetLimits(0.1-1e-4,99.9+1e-4);
if(!isHF){
   hFrame->GetXaxis()->SetLabelOffset(999);
   hFrame->GetXaxis()->SetTickLength(0);}
   hFrame->GetYaxis()->SetTitle("Event Selection Efficiency");
	hFrame->SetMaximum(1.09);
	hFrame->SetMinimum(0.25);
if(isHF)
hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
   hFrame->GetYaxis()->SetTitleSize(0.05);
        fixedFontHist(hFrame,1,1);
  hFrame->GetYaxis()->SetTitleOffset(1.5);
hFrame->GetYaxis()->SetTitle("Event Selection Efficiency");
hFrame->GetYaxis()->SetLabelSize(18);
hFrame->GetYaxis()->SetTitleOffset(1);

TLegend *leg=new TLegend(0.7,0.2,0.85,0.4);
if(isHF)
c1->cd()->SetLogx();
else
c1->cd();
hFrame->DrawCopy();
for(int i=0;i<2;i++){
j=0;k=1;
filename[i][k]=Form("output/%s/pPbHist_%s_DS",dir[i].Data(),Gen[k].Data());
f[i][k]=TFile::Open(Form("%s.root",filename[i][k].Data()));
h_or[i][k][j]=(TH1F*)f[i][k]->Get(histoname[j]);
h_tr[i][k][j]=(TH1F*)f[i][k]->Get(Form("%s_tr",histoname[j].Data()));
h_ds[i][k][j]=(TH1F*)f[i][k]->Get(Form("%s_DS",histoname[j].Data()));
//reh_or[i][k][j]=(TH1F*)h_or[i][k][j]->Rebin(nHFbin,Form("re%s%s%s",Gen[k].Data(),histoname[j].Data(),dir[i].Data()),HFbin);
//reh_tr[i][k][j]=(TH1F*)h_tr[i][k][j]->Rebin(nHFbin,Form("re%s%s%s_tr",Gen[k].Data(),histoname[j].Data(),dir[i].Data()),HFbin);
//reh_ds[i][k][j]=(TH1F*)h_ds[i][k][j]->Rebin(nHFbin,Form("re%s%s%s_ds",Gen[k].Data(),histoname[j].Data(),dir[i].Data()),HFbin);
reh_or[i][k][j]=(TH1F*)h_or[i][k][j]->Rebin(Ncent,Form("re%s%s%s",Gen[k].Data(),histoname[j].Data(),dir[i].Data()),HFbin_cent[j]);
reh_tr[i][k][j]=(TH1F*)h_tr[i][k][j]->Rebin(Ncent,Form("re%s%s%s_tr",Gen[k].Data(),histoname[j].Data(),dir[i].Data()),HFbin_cent[j]);
reh_ds[i][k][j]=(TH1F*)h_ds[i][k][j]->Rebin(Ncent,Form("re%s%s%s_ds",Gen[k].Data(),histoname[j].Data(),dir[i].Data()),HFbin_cent[j]);
normalizeByBinWidth(reh_or[i][k][j]);
normalizeByBinWidth(reh_tr[i][k][j]);
normalizeByBinWidth(reh_ds[i][k][j]);
gEff[i][k][j]=new TGraphAsymmErrors();     gEff[i][k][j]->SetName(Form("gEff_%s%s%s",Gen[k].Data(),histoname[j].Data(),dir[i].Data()));
gEff[i][k][j]->BayesDivide(reh_tr[i][k][j],reh_ds[i][k][j]);
gEff[i][k][j]->SetMarkerColor(1+i);
gEff[i][k][j]->SetLineColor(1+i);
gEff[i][k][j]->SetMarkerStyle(20+i*4);
gEff[i][k][j]->SetMarkerSize(1.2);
gEff[i][k][j]->SetTitle("");
gEff[i][k][j]->GetXaxis()->SetLimits(0.1-1e-4,99.9+1e-4);
if(isHF)
gEff[i][k][j]->Draw("Psame");
else{
ReverseXAxis(gEff[i][k][j],0.05);
ReverseXGraph(gEff[i][k][j]);
}
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(gEff[i][k][j],dir[i],"p");
TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);
}
T1[j]->Draw("same");
leg->Draw("same");
l->Draw("same");

if(isHF)
c1->Print(Form("pic/pPbHist_%s_DS_weightvsnoweight_%s_HF.png",Gen[k].Data(),histoname[j].Data()));
else
c1->Print(Form("pic/pPbHist_%s_DS_weightvsnoweight_%s_cent.png",Gen[k].Data(),histoname[j].Data()));

}
