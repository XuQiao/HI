#include "/afs/cern.ch/user/q/qixu/CMSSW_5_3_8_HI/src/jetRpA/RpA/Quality/root_setting.h"

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

void DrawEff_DS(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
//bool isweight=kFALSE;
bool isweight=kTRUE;
bool isHF=kTRUE;
//bool isHF=kFALSE;
if(isweight)
TString dir="weight";
else
TString dir="noweight";
const int Nhist=4;	const int NGen=1; const int Ncent=6;
//const double HFbin[]={0,0.1,0.2,0.4,0.6,0.8,1,2,4,6,9,12,15,20,30,50,75,100};
const double xerr[Ncent]={0,0,0,0,0,0};
/*
const double HFbin_cent[Nhist][Ncent+1]={
{0,2.92,10.1,19.3,23.7,30.2,100},
{0,1.03,5.00,11.6,15.0,20.3,100},
{0,2.7,12.6,28.6,36.9,49.6,100}};
*/
const double HFbin_cent[Nhist][Ncent+1]={
{0,2.725,10.065,19.565,23.975,30.455,100},
{0,0.975,5.005,11.775,15.285,20.725,100},
{0,2.515,12.465,28.855,37.235,50.085,100}};
const double cent[Ncent]={95,75,45,25,15,5};
const int nHFbin=sizeof(HFbin)/sizeof(double)-1;
//double color[NGen][Nhist]={{2,1,4},{2,1,4}};
//double Style[NGen][Nhist]={{20,24,28},{20,24,28}};
double color[NGen][Nhist]={{2,1,4}};
double Style[NGen][Nhist]={{20,24,28}};
//TString Gen[NGen]={"Hijing","Epos"};
TString Gen[NGen]={"Hijing"};
TString histoname[Nhist]={"HFEnergy4","HFEnergyPlus4","HFEnergyPlus","Ntrk"};
TString etarange[Nhist]={" HF E_{T} |#eta|>4"," HF E_{T} #eta>4"," HF E_{T} #eta>0",""};
TString filename[NGen];
TFile* f[NGen];
TH1F* h_or[NGen][Nhist];	TH1F* reh_or[NGen][Nhist];
TH1F* h_tr[NGen][Nhist];	TH1F* reh_tr[NGen][Nhist];
TH1F* h_ds[NGen][Nhist];	TH1F* reh_ds[NGen][Nhist];
TH1F* h_eff[NGen][Nhist];	TH1F* reh_eff[NGen][Nhist];

TLatex* T1[Nhist];	TLatex* T2;	TLatex* T3[NGen];
TF1 *fitEff[NGen][Nhist];
//TGraphAsymmErrors *gEff[NGen][Nhist];
TGraph *gEff[NGen][Nhist];
c1 = new TCanvas("c1"," ",800,500);
c2 = new TCanvas("c2"," ",800,700);
makeMultiPanelCanvas(c1,NGen,1,0.03,0.03,0.14,0.14,0.03);
makeMultiPanelCanvas(c2,1,1,0.03,0.03,0.03,0.03,0.03);

//for |eta|>4 HFEnergy:
/*
for(int k=0;k<NGen;k++)
for(int j=0;j<Nhist;j++)
fitEff[k][j]= new TF1(Form("fitEff_%s%s",Gen[k].Data(),histoname[j]),"[0]-[1]*exp([2]*x)-[3]*exp(-[4]*x)",0,300);

fitEff[0][0]->SetParameters(0.999876,0.44648,-1.78911,0.242036,0.549932);
fitEff[0][1]->SetParameters(1.,0.367894,-3.36629,0.104605,0.690272);
fitEff[0][2]->SetParameters(0.99981,0.697501,-2.59569,0.15323,0.411757);
fitEff[1][0]->SetParameters(0.999918,0.691843,-0.785042,0.0352768,0.449521);
fitEff[1][1]->SetParameters(0.999929,0.284,-2.25636,0.264114,0.951221);
fitEff[1][2]->SetParameters(0.999887,0.5,-2.39395,0.373215,0.550851);
*/

for(int k=0;k<NGen;k++){
T3[k]=new TLatex(0.45,0.93,Gen[k]);
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
j=3;
for(int k=0;k<NGen;k++){
//for(int j=0;j<Nhist;j++){
filename[k]=Form("output/%s/pPbHist_%s_DS",dir.Data(),Gen[k].Data());
f[k]=TFile::Open(Form("%s.root",filename[k].Data()));
h_or[k][j]=(TH1F*)f[k]->Get(Form("h%s",histoname[j].Data()));
h_tr[k][j]=(TH1F*)f[k]->Get(Form("h%s_tr",histoname[j].Data()));
h_ds[k][j]=(TH1F*)f[k]->Get(Form("h%s_DS",histoname[j].Data()));
h_eff[k][j]=(TH1F*)f[k]->Get(Form("hEff_%s",histoname[j].Data()));
fixedFontHist(h_or[k][j],1.8,2);
fixedFontHist(h_tr[k][j],1.8,2);
fixedFontHist(h_ds[k][j],1.8,2);
//reh_or[k][j]=(TH1F*)h_or[k][j]->Clone(Form("re%s%s",Gen[k].Data(),histoname[j].Data()));
//reh_tr[k][j]=(TH1F*)h_tr[k][j]->Clone(Form("re%s%s_tr",Gen[k].Data(),histoname[j].Data()));
//reh_ds[k][j]=(TH1F*)h_ds[k][j]->Clone(Form("re%s%s_ds",Gen[k].Data(),histoname[j].Data()));

//reh_or[k][j]=(TH1F*)h_or[k][j]->Rebin(Ncent,Form("re%s%s",Gen[k].Data(),histoname[j].Data()),HFbin_cent[j]);
//reh_tr[k][j]=(TH1F*)h_tr[k][j]->Rebin(Ncent,Form("re%s%s_tr",Gen[k].Data(),histoname[j].Data()),HFbin_cent[j]);
//reh_ds[k][j]=(TH1F*)h_ds[k][j]->Rebin(Ncent,Form("re%s%s_ds",Gen[k].Data(),histoname[j].Data()),HFbin_cent[j]);
reh_or[k][j]=(TH1F*)h_or[k][j]->Rebin(nHFbin,Form("re%s%s",Gen[k].Data(),histoname[j].Data()),HFbin);
reh_tr[k][j]=(TH1F*)h_tr[k][j]->Rebin(nHFbin,Form("re%s%s_tr",Gen[k].Data(),histoname[j].Data()),HFbin);
reh_ds[k][j]=(TH1F*)h_ds[k][j]->Rebin(nHFbin,Form("re%s%s_ds",Gen[k].Data(),histoname[j].Data()),HFbin);
//normalizeByBinWidth(reh_or[k][j]);
//normalizeByBinWidth(reh_tr[k][j]);
//normalizeByBinWidth(reh_ds[k][j]);

//gEff[k][j]=new TGraphAsymmErrors();
gEff[k][j]=new TGraph(h_eff[k][j]);
gEff[k][j]->SetPoint(0,0,1);
gEff[k][j]->SetName(Form("gEff_%s%s",Gen[k].Data(),histoname[j].Data()));
//gEff[k][j]->BayesDivide(reh_tr[k][j],reh_ds[k][j]);
gEff[k][j]->SetMarkerColor(color[k][j]);
gEff[k][j]->SetLineColor(color[k][j]);
gEff[k][j]->SetMarkerStyle(Style[k][j]);
gEff[k][j]->SetMarkerSize(1.2);
gEff[k][j]->GetXaxis()->SetLimits(0.01,99.9);
gEff[k][j]->SetTitle("");
//}
}

TH1D * hFrame = new TH1D("","",100,0,100.);
if(!isHF){
   hFrame->GetXaxis()->SetLabelOffset(999);
   hFrame->GetXaxis()->SetTickLength(0);}
   hFrame->GetYaxis()->SetTitle("Event Selection Efficiency");
	hFrame->SetMaximum(1.09);
	hFrame->SetMinimum(0);
if(isHF)
hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
if(j==3)
hFrame->GetXaxis()->SetTitle("track Mult");
   hFrame->GetYaxis()->SetTitleSize(0.05);
        fixedFontHist(hFrame);
  hFrame->GetYaxis()->SetTitleOffset(1.5);

for(k=0;k<NGen;k++){
if(isHF)
c1->cd(k+1)->SetLogx();
else
c1->cd(k+1);
TLegend *leg=new TLegend(0.60,0.2,0.75,0.4);
if(k!=0){
hFrame->GetYaxis()->SetTitle("");
hFrame->GetYaxis()->SetLabelSize(0);
}
hFrame->DrawCopy();
//for(j=0;j<Nhist;j++){
if(isHF)
gEff[k][j]->Draw("Psame");
else{
ReverseXAxis(gEff[k][j],0.05);
ReverseXGraph(gEff[k][j]);
}
//fitEff[k][j]->SetLineStyle(color[k][j]);
//fitEff[k][j]->SetLineColor(color[k][j]);
//fitEff[k][j]->Draw("same");
leg->AddEntry(gEff[k][j],etarange[j],"lp");
//}
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);

l->Draw("same");
leg->Draw("same");
T3[k]->Draw("same");
}
if(isHF)
c1->Print(Form("pic/OwnBinning/pPbHist_%s_DS_HF.png",dir.Data()));
else
c1->Print(Form("pic/OwnBinning/pPbHist_%s_DS_cent.png",dir.Data()));


//------------------------------------------------------------------------------------------

if(isHF)
//c2->cd()->SetLogx();
//else
c2->cd();
TLegend *leg=new TLegend(0.7,0.2,0.85,0.4);
hFrame->GetYaxis()->SetTitle("Event Selection Efficiency");
hFrame->GetYaxis()->SetLabelSize(18);
hFrame->GetYaxis()->SetTitleOffset(1);
hFrame->DrawCopy();
for(k=0;k<NGen;k++){
if(k==0){
gEff[k][j]->SetMarkerColor(1);
gEff[k][j]->SetLineColor(1);
gEff[k][j]->SetMarkerStyle(20);
gEff[k][j]->SetMarkerSize(1.6);
}
else{
gEff[k][j]->SetMarkerColor(2);
gEff[k][j]->SetLineColor(2);
gEff[k][j]->SetMarkerStyle(24);
gEff[k][j]->SetMarkerSize(1.6);
}
if(isHF)
gEff[k][j]->Draw("Psame");
else{
ReverseXAxis(gEff[k][j],0.05);
ReverseXGraph(gEff[k][j]);
}
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(gEff[k][j],Gen[k],"p");
TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);
}
T1[j]->Draw("same");
leg->Draw("same");
l->Draw("same");

if(isHF)
c2->Print(Form("pic/OwnBinning/pPbHist_%s_DS_%s_HF.png",dir.Data(),histoname[j].Data()));
else
c2->Print(Form("pic/OwnBinning/pPbHist_%s_DS_%s_cent.png",dir.Data(),histoname[j].Data()));

}
