#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

ReverseXAxis (TGraphAsymmErrors *g, double  size)
{
 //   Remove the current axis
   g->GetXaxis()->SetLabelOffset(999);
   g->GetXaxis()->SetTickLength(0);

 //   Redraw the new axis 
   gPad->Update();

//g->GetYaxis()->SetTitle("Event Selection Efficiency");
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
double x[100],xerrhigh[100],xerrlow[i],y[100],yerrhigh[100],yerrlow,xr[100];
const double cent[Ncent]={5,15,25,45,75,95};
  double dx=g->GetXaxis()->GetXmin()+g->GetXaxis()->GetXmax();
	for(int i=0;i<n;i++){
   xr[i] = dx-cent[n-i-1];
   xerrhigh[i] = g->GetErrorXhigh(i);
   xerrlow[i] = g->GetErrorXlow(i);
   y[i] = g->GetY()[i];
   yerrhigh[i] = g->GetErrorYhigh(i);
   yerrlow[i] = g->GetErrorYlow(i);
   }
   gr = new TGraphAsymmErrors(n,xr,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
   gr->SetMarkerSize(g->GetMarkerSize());
   gr->SetMarkerStyle(g->GetMarkerStyle());
   gr->SetMarkerColor(g->GetMarkerColor());
   gr->Draw("P");
}

void DrawEff(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
bool isweight=kFALSE;
//bool isweight=kTRUE;
bool isEvo=kTRUE;
//bool isEvo=kFALSE;
bool isHF=kTRUE;
//bool isHF=kFALSE;
if(isweight)
TString dir="weight";
else
TString dir="noweight";
const int Ncut=6;	const int Nhist=3;	const int NGen=2; const int Ncent=6;
const double HFbin[]={0,0.1,0.2,0.4,0.6,0.8,1,2,4,6,9,12,15,20,30,50,75,100};
const double xerr[Ncent]={0,0,0,0,0,0};
const double HFbin_cent[Nhist][Ncent+1]={
{0,2.92,10.1,19.3,23.7,30.2,100},
{0,1.03,5.00,11.6,15.0,20.3,100},
{0,2.7,12.6,28.6,36.9,49.6,100}};
const double cent[Ncent]={95,75,45,25,15,5};
const int nHFbin=sizeof(HFbin)/sizeof(double)-1;
double color[NGen][Nhist]={{2,1,4},{2,1,4}};
double Style[NGen][Nhist]={{20,24,28},{20,24,28}};
if(!isEvo){
TString cut[Ncut]={"HLT","CutGplus","pPA","pprim","vz",""};
TString cutName[Ncut]={"HLT","CutG","pPA","pprim","vz","All cuts"};
}
else{
TString cut[Ncut]={"HLT","HLT_CutGplus","HLT_CutGplus_pPA","HLT_CutGplus_pPA_pprim","",""};
TString cutName[Ncut]={"HLT","HLT+CutGplus","HLT+CutGplus+pPA","#splitline{HLT+CutG}{+pPA+pprim}","#splitline{HLT+CutG+pPA}{+pprim+vz}",""};
}
TString Gen[NGen]={"Hijing","Epos"};
TString histoname[Nhist]={"hHFEnergy4","hHFEnergyPlus4","hHFEnergyPlus"};
TString etarange[Nhist]={" HF E_{T} |#eta|>4"," HF E_{T} #eta>4"," #splitline{HF E_{T}#eta>0}{(Pb going side)}"};
TString filename[Ncut][NGen];
TFile* f[Ncut][NGen];
TH1F* h_or[Ncut][NGen][Nhist];	TH1F* reh_or[Ncut][NGen][Nhist];
TH1F* h_tr[Ncut][NGen][Nhist];	TH1F* reh_tr[Ncut][NGen][Nhist];
TGraphAsymmErrors *gEff[Ncut][NGen][Nhist];
TH1F* rehEff[Ncut][NGen][Nhist];	TLatex* T1[Nhist];	TLatex* T2[Ncut];	TLatex* T3[NGen];
c1 = new TCanvas("c1"," ",800,500);
c2 = new TCanvas("c2"," ",800,700);
c3 = new TCanvas("c3"," ",1000,800);
makeMultiPanelCanvas(c1,2,1,0.02,0.02,0.12,0.12,0.03);
makeMultiPanelCanvas(c2,1,1,0.03,0.03,0.03,0.03,0.03);
makeMultiPanelCanvas(c3,3,2,0.02,0.02,0.15,0.15,0.03);

for(int k=0;k<NGen;k++){
T3[k]=new TLatex(0.45,0.9,Gen[k]);
T3[k]->SetNDC();
T3[k]->SetTextAlign(12);
T3[k]->SetTextSize(0.06);
T3[k]->SetTextColor(1);
T3[k]->SetTextFont(42);
}
for(int j=0;j<Nhist;j++){
T1[j]=new TLatex(0.4,0.85,etarange[j]);
T1[j]->SetNDC();
T1[j]->SetTextAlign(12);
T1[j]->SetTextSize(0.06);
T1[j]->SetTextColor(1);
T1[j]->SetTextFont(42);
}
for(int i=0;i<Ncut;i++){
T2[i]=new TLatex(0.2,0.25,cutName[i]);
T2[i]->SetNDC();
T2[i]->SetTextAlign(12);
T2[i]->SetTextSize(0.08);
T2[i]->SetTextColor(1);
T2[i]->SetTextFont(42);
}
for(int k=0;k<NGen;k++){
for(int j=0;j<Nhist;j++){
for(int i=0;i<Ncut;i++){
if(cut[i]=="")
filename[i][k]=Form("output/%s/pPbHist_%s%s",dir.Data(),Gen[k].Data(),cut[i].Data());
else
filename[i][k]=Form("output/%s/pPbHist_%s_%s",dir.Data(),Gen[k].Data(),cut[i].Data());
f[i][k]=TFile::Open(Form("%s.root",filename[i][k].Data()));
h_or[i][k][j]=(TH1F*)f[i][k]->Get(histoname[j]);
h_tr[i][k][j]=(TH1F*)f[i][k]->Get(Form("%s_tr",histoname[j].Data()));
fixedFontHist(h_or[i][k][j],1.8,2);
fixedFontHist(h_tr[i][k][j],1.8,2);
reh_or[i][k][j]=(TH1F*)h_or[i][k][j]->Rebin(Ncent,Form("re%s%s_%d",Gen[k].Data(),histoname[j].Data(),i),HFbin_cent[j]);
reh_tr[i][k][j]=(TH1F*)h_tr[i][k][j]->Rebin(Ncent,Form("re%s%s_tr_%d",Gen[k].Data(),histoname[j].Data(),i),HFbin_cent[j]);
//reh_or[i][k][j]=(TH1F*)h_or[i][k][j]->Rebin(nHFbin,Form("re%s%s_%d",Gen[k].Data(),histoname[j].Data(),i),HFbin);
//reh_tr[i][k][j]=(TH1F*)h_tr[i][k][j]->Rebin(nHFbin,Form("re%s%s_tr_%d",Gen[k].Data(),histoname[j].Data(),i),HFbin);
normalizeByBinWidth(reh_or[i][k][j]);
normalizeByBinWidth(reh_tr[i][k][j]);
gEff[i][k][j]=new TGraphAsymmErrors();     gEff[i][k][j]->SetName(Form("gEff_%s%s%s",cut[i].Data(),Gen[k].Data(),histoname[j].Data()));
gEff[i][k][j]->BayesDivide(reh_tr[i][k][j],reh_or[i][k][j]);
gEff[i][k][j]->SetMarkerColor(color[k][j]);
gEff[i][k][j]->SetMarkerStyle(Style[k][j]);
gEff[i][k][j]->SetMarkerSize(1.2);
gEff[i][k][j]->GetXaxis()->SetLimits(0.1,99);
gEff[i][k][j]->SetMaximum(1.19);
gEff[i][k][j]->SetMinimum(0.4);
gEff[i][k][j]->SetTitle("");
}
}
}

i=Ncut-1;	//All cut
TLegend *leg=new TLegend(0.40,0.2,0.55,0.4);
 TH1D * hFrame = new TH1D("","",100,0,100.);
      // hFrame-> SetAxisRange(0,100.,"X");
       hFrame->GetXaxis()->SetLimits(0.1,99.9);
  //  hFrame->SetAxisRange(0.35, 1.2, "Y") ;
  	hFrame->SetMinimum(0.35);
  	hFrame->SetMaximum(1.09);
if(!isHF){
   hFrame->GetXaxis()->SetLabelOffset(999);
   hFrame->GetXaxis()->SetTickLength(0);
hFrame->GetXaxis()->SetTitle("Centrality");
}
else	{hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");}
   hFrame->GetYaxis()->SetTitle("Event Selection Efficiency");
   hFrame->GetYaxis()->SetTitleSize(0.05);
	fixedFontHist(hFrame);
  hFrame->GetYaxis()->SetTitleOffset(1.2);
for(k=0;k<NGen;k++){
if(isHF)
c1->cd(k+1)->SetLogx();
else
c1->cd(k+1);
if(k!=0){
hFrame->GetYaxis()->SetTitle("");
hFrame->GetYaxis()->SetLabelSize(0);
}
hFrame->DrawCopy();
for(j=0;j<Nhist;j++){
if(isHF)
gEff[i][k][j]->Draw("Psame");
else{
ReverseXAxis(gEff[i][k][j],0.05);
ReverseXGraph(gEff[i][k][j]);
}
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(gEff[i][k][j],etarange[j],"p");
TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);
}
l->Draw("same");
leg->Draw("same");
T3[k]->Draw("same");
}
if(isHF)
c1->Print(Form("pic/pPbHist_%s_%s_HF.png",dir.Data(),cut[i].Data()));
else
c1->Print(Form("pic/pPbHist_%s_%s_cent.png",dir.Data(),cut[i].Data()));

//------------------------------------------------------------------------------------------

i=Ncut-1;j=1;
if(isHF)
c2->cd()->SetLogx();
else
c2->cd();
TLegend *leg=new TLegend(0.5,0.2,0.65,0.3);
hFrame->GetYaxis()->SetLabelSize(18);
hFrame->GetYaxis()->SetTitle("Event Selection Efficiency");
hFrame->DrawCopy();
for(k=0;k<NGen;k++){
if(k==0){
gEff[i][k][j]->SetMarkerColor(1);
gEff[i][k][j]->SetLineColor(1);
gEff[i][k][j]->SetMarkerStyle(24);
gEff[i][k][j]->SetMarkerSize(1.6);
}
else{
gEff[i][k][j]->SetMarkerColor(2);
gEff[i][k][j]->SetLineColor(2);
gEff[i][k][j]->SetMarkerStyle(20);
}
if(isHF)
gEff[i][k][j]->Draw("Psame");
else{
ReverseXAxis(gEff[i][k][j],0.04);
ReverseXGraph(gEff[i][k][j]);
}
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(gEff[i][k][j],Gen[k],"p");
TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);
}
T1[j]->Draw("same");
leg->Draw("same");
l->Draw("same");

if(isHF)
c2->Print(Form("pic/pPbHist_%s_%s_%s_HF.png",dir.Data(),cut[i].Data(),histoname[j].Data()));
else
c2->Print(Form("pic/pPbHist_%s_%s_%s_cent.png",dir.Data(),cut[i].Data(),histoname[j].Data()));

//----------------------------------------------------------------------------------

j=1;k=0;
for(i=0;i<Ncut;i++){
gEff[i][k][j]->SetMarkerStyle(20);
gEff[i][k][j]->SetMarkerColor(1);
if(isHF)
c3->cd(i+1)->SetLogx();
else
c3->cd(i+1);
if(i==0){
hFrame->GetYaxis()->SetTitle("Event Selection Efficiency");
hFrame->GetYaxis()->SetLabelSize(18);
hFrame->GetXaxis()->SetTitle("");
hFrame->GetXaxis()->SetLabelSize(0);
hFrame->DrawCopy();
}
else if(i==1 || i==2){
hFrame->GetYaxis()->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetLabelSize(0);
hFrame->GetXaxis()->SetLabelSize(0);
hFrame->DrawCopy();
}
else if(i==3){
hFrame->GetYaxis()->SetTitle("Event Selection Efficiency");
hFrame->GetYaxis()->SetLabelSize(1);
if(isHF)
hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
else
hFrame->GetXaxis()->SetTitle("Centrality");
hFrame->GetXaxis()->SetLabelSize(18);
hFrame->DrawCopy();
}
else{
if(isHF)
hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
else
hFrame->GetXaxis()->SetTitle("Centrality");
hFrame->GetXaxis()->SetLabelSize(18);
hFrame->GetYaxis()->SetLabelSize(0);
hFrame->GetYaxis()->SetTitle("");
hFrame->DrawCopy();
}
if(isHF)
gEff[i][k][j]->Draw("Psame");
else{
ReverseXAxis(gEff[i][k][j],0.06);
ReverseXGraph(gEff[i][k][j]);
}

TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);
l->Draw("same");
T2[i]->Draw("same");
}
T3[k]->Draw("same");

if(isEvo)
if(isHF)
c3->Print(Form("pic/pPbHist_%s_%s_%s_Evo_HF.png",dir.Data(),Gen[k].Data(),histoname[j].Data()));
else
c3->Print(Form("pic/pPbHist_%s_%s_%s_Evo_cent.png",dir.Data(),Gen[k].Data(),histoname[j].Data()));
else
if(isHF)
c3->Print(Form("pic/pPbHist_%s_%s_%s_HF.png",dir.Data(),Gen[k].Data(),histoname[j].Data()));
else
c3->Print(Form("pic/pPbHist_%s_%s_%s_cent.png",dir.Data(),Gen[k].Data(),histoname[j].Data()));
}
