#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
#include "/home/xuq7/HI/jetRpA/RpA/NewMC/produceandcheck/file.h"

TLatex *T=new TLatex();

int ieta=7;
TGraphAsymmErrors* makegraph(TH1* histo1, TH1* histo2){
const int Npoint=1000;
double x[Npoint],y[Npoint],xerrlow[Npoint],xerrhigh[Npoint],yerrlow[Npoint],yerrhigh[Npoint];
TGraphAsymmErrors *g=new TGraphAsymmErrors();
g->BayesDivide(histo1,histo2);
for(int ipoint=0;ipoint<g->GetN();ipoint++){
x[ipoint]=g->GetX()[ipoint];
y[ipoint]=g->GetY()[ipoint];
yerrhigh[ipoint]=g->GetErrorYhigh(ipoint);
yerrlow[ipoint]=g->GetErrorYlow(ipoint);
xerrhigh[ipoint]=0;
xerrlow[ipoint]=0;
}
g0=new TGraphAsymmErrors(g->GetN(),x,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
return g0;
}

TGraphAsymmErrors* makeJetIDcut(int ilist, int dataType){
int ilistw;
if(ilist== 6)       {   JetID = "Max p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0]=0.05; JetIDcut[1]=1; ilistw=1;}
else if(ilist== 7)  {    JetID = "#Sum p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "Max p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "#Sum p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10) {   JetID = "Max p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "#Sum p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{jet}";JetIDcut[0]=0; JetIDcut[1]=1.01; ilistw=3;}
else if(ilist== 13) {    JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{raw}"; JetIDcut[0]=0;JetIDcut[1]=1.13;}
else if(ilist== 14) {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})";JetIDcut[0]=0;JetIDcut[1]=0.975;ilistw=2;}

if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};

int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;

TString JetIDName=JetIDNameList[ilist];
TString histonameIDData=Form("jetpt%s",JetIDName.Data());
TH2F* hdata2F=(TH2F*)fDataPPbJetID->Get(histonameIDData);
TH2F *hPPb2D= (TH2F*)fMCPPb->Get(Form("jetpt%s",JetIDName.Data()));
TH2F *hPPb2D_real= (TH2F*)fMCPPb->Get(Form("jetpt%s_real",JetIDName.Data()));
TH2F *hPPb2D_fake= (TH2F*)fMCPPb->Get(Form("jetpt%s_fake",JetIDName.Data()));

const double binbound_pt_coarse[]={0,30,70,150,200,600};
const int Nbin_pt_coarse=sizeof(binbound_pt_coarse)/sizeof(double)-1;
const double xrange_JetID[2]={binbound_JetID[0]+1e-4,binbound_JetID[Nbin_JetID]-1e-4};
xrange_JetIDcut[0]=JetIDcut[0]+1e-4;
xrange_JetIDcut[1]=JetIDcut[1]-1e-4;

//TH1D* hPPb_pt=hPPb2D->ProjectionX("hPPb_pt",hPPb2D->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPPb_pt=hPPb2D->ProjectionX("hPPb_pt");
TH1D* hPPb_pt_real=hPPb2D_real->ProjectionX("hPPb_pt_real");
TH1D* hPPb_pt_fake=hPPb2D_fake->ProjectionX("hPPb_pt_fake");
TH1D* hPPb_JetIDcut_pt=hPPb2D->ProjectionX("hPPb_JetIDcut_pt",hPPb2D->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* hPPb_JetIDcut_pt_real=hPPb2D_real->ProjectionX("hPPb_JetIDcut_pt_real",hPPb2D_real->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_real->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* hPPb_JetIDcut_pt_fake=hPPb2D_fake->ProjectionX("hPPb_JetIDcut_pt_fake",hPPb2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[1]));

hPPb_JetIDcut_pt_real = (TH1D*)hPPb_JetIDcut_pt_real->Rebin(Nbin_pt,hPPb_JetIDcut_pt_real->GetName(),binbound_pt);
hPPb_JetIDcut_pt = (TH1D*)hPPb_JetIDcut_pt->Rebin(Nbin_pt,hPPb_JetIDcut_pt->GetName(),binbound_pt);
hPPb_pt_real = (TH1D*)hPPb_pt_real->Rebin(Nbin_pt,hPPb_pt_real->GetName(),binbound_pt);
hPPb_pt = (TH1D*)hPPb_pt->Rebin(Nbin_pt,hPPb_pt->GetName(),binbound_pt);

TGraphAsymmErrors *g = new TGraphAsymmErrors();
TH1D* histo1=(TH1D*)hdata2F->ProjectionX("histo1",hdata2F->GetYaxis()->FindBin(xrange_JetIDcut[0]),hdata2F->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* histo2=(TH1D*)hdata2F->ProjectionX("histo2");
histo1=(TH1D*)histo1->Rebin(Nbin_pt,"histo1",binbound_pt);
histo2=(TH1D*)histo2->Rebin(Nbin_pt,"histo2",binbound_pt);
normalizeByBinWidth(histo1);
normalizeByBinWidth(histo2);
g=makegraph(histo1,histo2);
g->SetMarkerSize(1.5);

TGraphAsymmErrors *gfake = new TGraphAsymmErrors();
TGraphAsymmErrors *greal = new TGraphAsymmErrors();
TGraphAsymmErrors *ginc = new TGraphAsymmErrors();
TH1D* ratio_hPPb_JetIDcutvsnocut_pt=(TH1D*)hPPb_JetIDcut_pt->Clone("ratio_hPPb_JetIDcut_ptcutvsnocut_pt");
ginc=makegraph(ratio_hPPb_JetIDcutvsnocut_pt,hPPb_pt);
ginc->SetMarkerSize(1.5);

/*TH1D* ratio_hPPb_JetIDcutvsnocut_pt_fake=(TH1D*)hPPb_JetIDcut_pt_fake->Clone("ratio_hPPb_JetIDcutvsnocut_pt_fake");
TH1D* cl_hPPb_pt_fake=(TH1D*)hPPb_pt_fake->Clone("cl_hPPb_pt_fake");
cl_hPPb_pt_fake=(TH1D*)cl_hPPb_pt_fake->Rebin(Nbin_pt_coarse,"cl_hPPb_pt_fake",binbound_pt_coarse);
ratio_hPPb_JetIDcutvsnocut_pt_fake=(TH1D*)ratio_hPPb_JetIDcutvsnocut_pt_fake->Rebin(Nbin_pt_coarse,"ratio_hPPb_JetIDcutvsnocut_pt_fake",binbound_pt_coarse);
normalizeByBinWidth(ratio_hPPb_JetIDcutvsnocut_pt_fake);
normalizeByBinWidth(cl_hPPb_pt_fake);
gfake=makegraph(ratio_hPPb_JetIDcutvsnocut_pt_fake,cl_hPPb_pt_fake);
gfake->SetMarkerSize(1.5);
gfake->SetMarkerStyle(20);
gfake->SetMarkerColor(2);
gfake->SetLineColor(2);

TH1D* ratio_hPPb_JetIDcutvsnocut_pt_real=(TH1D*)hPPb_JetIDcut_pt_real->Clone("ratio_hPPb_JetIDcutvsnocut_pt_real");
greal=makegraph(ratio_hPPb_JetIDcutvsnocut_pt_real,hPPb_pt_real);
greal->SetMarkerSize(1.5);
greal->SetMarkerStyle(29);
greal->SetMarkerColor(4);
greal->SetLineColor(4);
*/
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.03);
T->SetTextColor(1);
T->SetTextFont(42);
if(ilistw==2){
if(dataType==1) T->DrawLatex(0.22,0.98-ilistw*0.05,Form("Cut %d: %.2f<%s<%.3f",ilistw,JetID.Data(),xrange_JetIDcut[0],xrange_JetIDcut[1]));}
else{
if(dataType==1) T->DrawLatex(0.22,0.98-ilistw*0.05,Form("Cut %d: %.2f<%s<%.2f",ilistw,JetID.Data(),xrange_JetIDcut[0],xrange_JetIDcut[1]));}

if(dataType==1) return g;
else return ginc;
}


void overlayJetIDcut(){
gStyle->SetOptStat(kFALSE);
c1 = new TCanvas("c1","",600,600);
makeMultiPanelCanvas(c1,1,1,-0.12,0,0.13,0.1,0.03);
TH1F* hFrame=new TH1F("","",1000,0,1000);
c1->cd(1)->SetGridx();
hFrame->SetLineColor(0);
hFrame->GetXaxis()->SetLimits(28,599.9);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Cut Efficiency");
hFrame->GetYaxis()->SetRangeUser(0.92,1.02);
fixedFontHist(hFrame,1.2,1.8);
hFrame->DrawCopy();
TGraphAsymmErrors *g6Data = makeJetIDcut(6,1);
g6Data->SetMarkerStyle(20);
g6Data->SetMarkerColor(1);
g6Data->SetLineColor(1);
g6Data->Draw("Psame");
TGraphAsymmErrors *g6MC = makeJetIDcut(6,0);
g6MC->SetMarkerStyle(24);
g6MC->SetMarkerColor(1);
g6MC->SetLineColor(1);
TGraphAsymmErrors *g12Data = makeJetIDcut(12,1);
g12Data->SetMarkerStyle(20);
g12Data->SetMarkerColor(1);
g12Data->SetLineColor(1);
g12Data->Draw("Psame");
TGraphAsymmErrors *g12MC = makeJetIDcut(12,0);
g12MC->SetMarkerStyle(24);
g12MC->SetMarkerColor(1);
g12MC->SetLineColor(1);
g12MC->Draw("Psame");
TGraphAsymmErrors *g14Data = makeJetIDcut(14,1);
g14Data->SetMarkerStyle(21);
g14Data->SetMarkerColor(4);
g14Data->SetLineColor(4);
g14Data->Draw("Psame");
TGraphAsymmErrors *g14MC = makeJetIDcut(14,0);
g14MC->SetMarkerStyle(25);
g14MC->SetMarkerColor(4);
g14MC->SetLineColor(4);
g14MC->Draw("Psame");
g6MC->Draw("Psame");
TLegend *leg1=new TLegend(0.52,0.15,0.8,0.45);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg1->AddEntry(g6Data,"cut 1 Data","lp");
leg1->AddEntry(g6MC,"cut 1 PYTHIA+HIJING","lp");
leg1->AddEntry(g14Data,"cut 2 Data","lp");
leg1->AddEntry(g14MC,"cut 2 PYTHIA+HIJING","lp");
//leg1->AddEntry(g12Data,"cut 3 Data","lp");
//leg1->AddEntry(g12MC,"cut 3 PYTHIA+HIJING","lp");
leg1->Draw("same");
T->SetTextSize(0.04);
T->DrawLatex(0.30,0.20,etastring[ieta]);
TLine *l =new TLine(25,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");

c1->Print(Form("pic/overlay_JetIDcut.png"));
c1->Print(Form("pic/overlay_JetIDcut.pdf"));
}
