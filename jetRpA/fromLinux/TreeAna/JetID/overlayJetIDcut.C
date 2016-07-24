#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

TLatex *T=new TLatex();
TString coll="PPb";
int ieta=0;
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
TString JetID, Unit;
int ilistw;
double JetIDcut[2];
if(ilist== 0)       {   JetID = "Max p_{T}^{charged}"; JetIDcut[0]=4.00; JetIDcut[1]=300.00; ilistw=1;}
if(ilist== 6)       {   JetID = "Max p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0]=0.05; JetIDcut[1]=1; ilistw=2;}
else if(ilist== 7)  {    JetID = "#Sum p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "Max p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "#Sum p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10) {   JetID = "Max p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "#Sum p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{jet}";JetIDcut[0]=0; JetIDcut[1]=1.01; ilistw=3;}
else if(ilist== 13) {    JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{raw}"; JetIDcut[0]=0;JetIDcut[1]=1.13;}
else if(ilist== 14) {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})";JetIDcut[0]=0;JetIDcut[1]=0.975;ilistw=3;}
if(ilist==0){
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,300};}

else if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};

int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;

TString JetIDName=JetIDNameList[ilist];
TString histonameID=Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[ieta].Data());
if(coll=="PPb"){
TH2F* hdata2F=(TH2F*)fDataPPbJetID->Get(histonameID);
TH2F *hPPb2D= (TH2F*)fMCPPbJetID->Get(Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[ieta].Data()));
}
else{
TH2F* hdata2F=(TH2F*)fDataPbPJetIDOld->Get(histonameID);
TH2F *hPPb2D= (TH2F*)fMCPbPJetIDOld->Get(Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinname[ieta].Data()));
}
//TH2F *hPPb2D_real= (TH2F*)fMCPPbJetID->Get(Form("jetpt%sEtaBin_real%s",JetIDName.Data(),etabinnameswap[ieta].Data()));
//TH2F *hPPb2D_fake= (TH2F*)fMCPPbJetID->Get(Form("jetpt%sEtaBin_fake%s",JetIDName.Data(),etabinnameswap[ieta].Data()));

const double binbound_pt_coarse[]={0,30,70,150,200,500};
const int Nbin_pt_coarse=sizeof(binbound_pt_coarse)/sizeof(double)-1;
const double xrange_JetID[2]={binbound_JetID[0]+1e-4,binbound_JetID[Nbin_JetID]-1e-4};
double xrange_JetIDcut[2];
xrange_JetIDcut[0]=JetIDcut[0]+1e-4;
xrange_JetIDcut[1]=JetIDcut[1]-1e-4;

TH1D* hPPb_pt=hPPb2D->ProjectionX("hPPb_pt",hPPb2D->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D->GetYaxis()->FindBin(xrange_JetID[1]));
//TH1D* hPPb_pt_real=hPPb2D_real->ProjectionX("hPPb_pt_real",hPPb2D_real->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D_real->GetYaxis()->FindBin(xrange_JetID[1]));
//TH1D* hPPb_pt_fake=hPPb2D_fake->ProjectionX("hPPb_pt_fake",hPPb2D_fake->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D_fake->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPPb_JetIDcut_pt=hPPb2D->ProjectionX("hPPb_JetIDcut_pt",hPPb2D->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D->GetYaxis()->FindBin(xrange_JetIDcut[1]));
//TH1D* hPPb_JetIDcut_pt_real=hPPb2D_real->ProjectionX("hPPb_JetIDcut_pt_real",hPPb2D_real->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_real->GetYaxis()->FindBin(xrange_JetIDcut[1]));
//TH1D* hPPb_JetIDcut_pt_fake=hPPb2D_fake->ProjectionX("hPPb_JetIDcut_pt_fake",hPPb2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[1]));

//hPPb_JetIDcut_pt_real = (TH1D*)hPPb_JetIDcut_pt_real->Rebin(Nbin_pt,hPPb_JetIDcut_pt_real->GetName(),binbound_pt);
hPPb_JetIDcut_pt = (TH1D*)hPPb_JetIDcut_pt->Rebin(Nbin_pt,hPPb_JetIDcut_pt->GetName(),binbound_pt);
//hPPb_pt_real = (TH1D*)hPPb_pt_real->Rebin(Nbin_pt,hPPb_pt_real->GetName(),binbound_pt);
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
/*
TH1D* ratio_hPPb_JetIDcutvsnocut_pt_fake=(TH1D*)hPPb_JetIDcut_pt_fake->Clone("ratio_hPPb_JetIDcutvsnocut_pt_fake");
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
if(dataType==1){
if(ilist==14)
T->DrawLatex(0.22,0.98-ilistw*0.05,Form("Cut %d: %.2f<%s<%.3f",ilistw,JetID.Data(),xrange_JetIDcut[0],xrange_JetIDcut[1]));
else 
T->DrawLatex(0.22,0.98-ilistw*0.05,Form("Cut %d: %.2f<%s<%.2f",ilistw,JetID.Data(),xrange_JetIDcut[0],xrange_JetIDcut[1]));
return g;
}
else return ginc;
}


void overlayJetIDcut(){
gStyle->SetOptStat(kFALSE);
c1 = new TCanvas("c1","",600,600);
makeMultiPanelCanvas(c1,1,1,-0.12,0,0.13,0.1,0.03);
TH1F* hFrame=new TH1F("","",1000,0,1000);
c1->cd(1);
hFrame->SetLineColor(0);
if(ieta!=0 && ieta!=6)
hFrame->GetXaxis()->SetLimits(45,480);
else
hFrame->GetXaxis()->SetLimits(45,370);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Cut Efficiency");
hFrame->GetYaxis()->SetRangeUser(0.78,1.06);
fixedFontHist(hFrame,1.2,1.8);
hFrame->DrawCopy();
TGraphAsymmErrors *g6Data = makeJetIDcut(6,1);
g6Data->SetMarkerStyle(21);
g6Data->SetMarkerColor(4);
g6Data->SetLineColor(4);
TGraphAsymmErrors *g6MC = makeJetIDcut(6,0);
g6MC->SetMarkerStyle(25);
g6MC->SetMarkerColor(4);
g6MC->SetLineColor(4);
TGraphAsymmErrors *g12Data = makeJetIDcut(0,1);
g12Data->SetMarkerStyle(20);
g12Data->SetMarkerColor(1);
g12Data->SetLineColor(1);
TGraphAsymmErrors *g12MC = makeJetIDcut(0,0);
g12MC->SetMarkerStyle(24);
g12MC->SetMarkerColor(1);
g12MC->SetLineColor(1);
TGraphAsymmErrors *g14Data = makeJetIDcut(14,1);
g14Data->SetMarkerStyle(33);
g14Data->SetMarkerColor(2);
g14Data->SetLineColor(2);
TGraphAsymmErrors *g14MC = makeJetIDcut(14,0);
g14MC->SetMarkerStyle(27);
g14MC->SetMarkerColor(2);
g14MC->SetLineColor(2);
g6Data->Draw("Psame");
g12Data->Draw("Psame");
g14Data->Draw("Psame");
g6MC->Draw("Psame");
g12MC->Draw("Psame");
g14MC->Draw("Psame");
TLegend *leg1=new TLegend(0.42,0.15,0.67,0.45);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
if(coll=="PPb")
TString mc = "PYTHIA+HIJING";
else
TString mc = "HIJING+PYTHIA";
leg1->AddEntry(g12Data,Form("cut 1 %s Data",coll.Data()),"lp");
leg1->AddEntry(g12MC,Form("cut 1 %s",mc.Data()),"lp");
leg1->AddEntry(g6Data,Form("cut 2 %s Data",coll.Data()),"lp");
leg1->AddEntry(g6MC,Form("cut 2 %s",mc.Data()),"lp");
leg1->AddEntry(g14Data,Form("cut 3 %s Data",coll.Data()),"lp");
leg1->AddEntry(g14MC,Form("cut 3 %s",mc.Data()),"lp");
leg1->Draw("same");
T->SetTextSize(0.04);
T->DrawLatex(0.18,0.20,etastring[ieta]);
TLine *l =new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");

c1->Print(Form("%s_overlay_JetIDcut_%s.png",coll.Data(),etabinname[ieta].Data()));
c1->Print(Form("%s_overlay_JetIDcut_%s.pdf",coll.Data(),etabinname[ieta].Data()));
}
