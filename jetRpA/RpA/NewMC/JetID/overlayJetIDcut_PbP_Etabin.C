#include "../../Quality/root_setting.h"
#include "../produceandcheck/file.h"
#include "extract.C"


//int ieta=7;
TGraphAsymmErrors* makegraph(TH1* histo1, TH1* histo2){
const int Npoint=1000;
double x[Npoint],y[Npoint],xerrlow[Npoint],xerrhigh[Npoint],yerrlow[Npoint],yerrhigh[Npoint];
TGraphAsymmErrors *g=new TGraphAsymmErrors();
g->BayesDivide(histo1,histo2);
for(int ipoint=0;ipoint<g->GetN();ipoint++){
x[ipoint]=g->GetX()[ipoint];
y[ipoint]=g->GetY()[ipoint];
yerrhigh[ipoint]=0;
yerrlow[ipoint]=0;
//yerrhigh[ipoint]=g->GetErrorYhigh(ipoint);
//yerrlow[ipoint]=g->GetErrorYlow(ipoint);
xerrhigh[ipoint]=0;
xerrlow[ipoint]=0;
}
g0=new TGraphAsymmErrors(g->GetN(),x,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
return g0;
}
double JetIDcut[2];

TH1D* makeJetIDcut(int ilist, int dataType, int ieta){
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
else if(ilist== 21)   {  JetID = "PP cut Tighter True or False";JetIDcut[0]=2;JetIDcut[1]=5;}

if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};

int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;

TString JetIDName=JetIDNameList[ilist];
TString histoname, histonameID;
if(ieta==7){
    histonameID = Form("jetpt%s",JetIDName.Data());
    histoname = Form("jetpt");
}
else{
    histonameID = Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinname[ieta].Data());
    histoname = Form("jetptEtaBin%s",etabinname[ieta].Data());
}

TH2F* hdata2F=(TH2F*)fDataPbPJetID->Get(histonameID);
TH2F *hPbP2D= (TH2F*)fMCPbP->Get(histonameID);
TH2F *hPbP2D_real= (TH2F*)fMCPbP->Get(Form("%s_real1",histonameID.Data()));
TH2F *hPbP2D_fake= (TH2F*)fMCPbP->Get(Form("%s_fake1",histonameID.Data()));
TH1D* hdata_pt=(TH1D*)fDataPbPJetID->Get(histoname);
TH1D *hPbP_pt= (TH1D*)fMCPbP->Get(histoname);
TH1D *hPbP_pt_real= (TH1D*)fMCPbP->Get(Form("%s_real1",histoname.Data()));
TH1D *hPbP_pt_fake= (TH1D*)fMCPbP->Get(Form("%s_fake1",histoname.Data()));
if(dataType==2 || dataType == 4){
TH1D* hPbP_JetIDcut_pt_real = new TH1D(Form("hPbP_JetIDcut_pt_real_%d_%d",ieta,dataType),Form("hPbP_JetIDcut_pt_real_%d",ieta),1000,0.,1000.);
hPbP_JetIDcut_pt_real->Sumw2();
}
if(dataType==3){
TH1D* hPbP_JetIDcut_pt_fake = new TH1D(Form("hPbP_JetIDcut_pt_fake_%d_%d",ieta,dataType),Form("hPbP_JetIDcut_pt_fake_%d",ieta),1000,0.,1000.);
hPbP_JetIDcut_pt_fake->Sumw2();
}
if(dataType==1 || dataType == 4){
TH1D* hPbP_JetIDcut_pt = new TH1D(Form("hPbP_JetIDcut_pt_%d_%d",ieta,dataType),Form("hPbP_JetIDcut_pt_%d",ieta),1000,0.,1000.);
hPbP_JetIDcut_pt->Sumw2();
}
if(dataType==0){
TH1D* hdata_JetIDcut_pt = new TH1D(Form("hdata_JetIDcut_pt_%d_%d",ieta,dataType),Form("hdata_JetIDcut_pt_%d",ieta),1000,0.,1000.);
hdata_JetIDcut_pt->Sumw2();
}
if(ilist==21){
for(int ibin=1;ibin<hPbP2D_real->GetYaxis()->GetNbins();ibin++){
	double binxcenter = hPbP2D_real->GetYaxis()->GetBinCenter(ibin);
		TVectorD* vecx2 = (TVectorD*)extract(binxcenter);
                TH1D* histo_1=hPbP2D_real->ProjectionX(Form("histo1_%d",ibin),ibin,ibin);
                TH1D* histo_2=hPbP2D_fake->ProjectionX(Form("histo2_%d",ibin),ibin,ibin);
                TH1D* histo_0=hPbP2D->ProjectionX(Form("histo0_%d",ibin),ibin,ibin);
                TH1D* datahisto=hdata2F->ProjectionX(Form("datahisto0_%d",ibin),ibin,ibin);
                int a=(int)(JetIDcut[0]);
		if((*vecx2)[a]!=0){
                if(dataType==2 || dataType ==4)
		hPbP_JetIDcut_pt_real->Add(histo_1);
                if(dataType==3)
		hPbP_JetIDcut_pt_fake->Add(histo_2);
                if(dataType==1 || dataType ==4 )
		hPbP_JetIDcut_pt->Add(histo_0);
                if(dataType==0)
		hdata_JetIDcut_pt->Add(datahisto);
                }
	}
}/*
else{
    double xrange_JetID[2];
    xrange_JetIDcut[0]=JetIDcut[0]+1e-4;
    xrange_JetIDcut[1]=JetIDcut[1]-1e-4;
    TH1D* hPbP_JetIDcut_pt_real=hPbP2D_real->ProjectionX(Form("hPbP_JetIDcut_pt_real_%d",ieta),hPbP2D_real->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPbP2D_real->GetYaxis()->FindBin(xrange_JetIDcut[1]));
    TH1D* hPbP_JetIDcut_pt_fake=hPbP2D_fake->ProjectionX(Form("hPbP_JetIDcut_pt_fake_%d",ieta),hPbP2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPbP2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[1]));
    TH1D* hPbP_JetIDcut_pt=hPbP2D->ProjectionX(Form("hPbP_JetIDcut_pt_%d",ieta),hPbP2D->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPbP2D->GetYaxis()->FindBin(xrange_JetIDcut[1]));
    TH1D* hdata_JetIDcut_pt=hdata2F->ProjectionX(Form("hdata_JetIDcut_pt_%d",ieta),hdata2F->GetYaxis()->FindBin(xrange_JetIDcut[0]),hdata2F->GetYaxis()->FindBin(xrange_JetIDcut[1]));
}
*/
if(dataType==0){
hdata_JetIDcut_pt = (TH1D*)hdata_JetIDcut_pt->Rebin(Nbin_pt,hdata_JetIDcut_pt->GetName(),binbound_pt);
hdata_JetIDcutvsnocut_pt = (TH1D*)hdata_JetIDcut_pt->Clone("hdata_JetIDcutvsnocut_pt");
hdata_pt = (TH1D*)hdata_pt->Rebin(Nbin_pt,hdata_pt->GetName(),binbound_pt);
hdata_JetIDcutvsnocut_pt->Divide(hdata_pt);
getRidYError(hdata_JetIDcutvsnocut_pt);
}
if(dataType==1 || dataType ==4){
hPbP_JetIDcut_pt = (TH1D*)hPbP_JetIDcut_pt->Rebin(Nbin_pt,hPbP_JetIDcut_pt->GetName(),binbound_pt);
hPbP_JetIDcutvsnocut_pt = (TH1D*)hPbP_JetIDcut_pt->Clone("hPbP_JetIDcutvsnocut_pt");
hPbP_pt = (TH1D*)hPbP_pt->Rebin(Nbin_pt,hPbP_pt->GetName(),binbound_pt);
hPbP_JetIDcutvsnocut_pt->Divide(hPbP_pt);
getRidYError(hPbP_JetIDcutvsnocut_pt);
}
if(dataType==2 || dataType==4){
hPbP_JetIDcut_pt_real = (TH1D*)hPbP_JetIDcut_pt_real->Rebin(Nbin_pt,hPbP_JetIDcut_pt_real->GetName(),binbound_pt);
hPbP_JetIDcutvsnocut_pt_real = (TH1D*)hPbP_JetIDcut_pt_real->Clone("hPbP_JetIDcutvsnocut_pt_real");
hPbP_pt_real = (TH1D*)hPbP_pt_real->Rebin(Nbin_pt,hPbP_pt_real->GetName(),binbound_pt);
hPbP_JetIDcutvsnocut_pt_real->Divide(hPbP_pt_real);
getRidYError(hPbP_JetIDcutvsnocut_pt_real);
}

if(dataType==3){
hPbP_JetIDcut_pt_fake = (TH1D*)hPbP_JetIDcut_pt_fake->Rebin(Nbin_pt,hPbP_JetIDcut_pt_fake->GetName(),binbound_pt);
hPbP_JetIDcutvsnocut_pt_fake = (TH1D*)hPbP_JetIDcut_pt_fake->Clone("hPbP_JetIDcutvsnocut_pt_fake");
hPbP_pt_fake = (TH1D*)hPbP_pt_fake->Rebin(Nbin_pt,hPbP_pt_fake->GetName(),binbound_pt);
hPbP_JetIDcutvsnocut_pt_fake->Divide(hPbP_pt_fake);
getRidYError(hPbP_JetIDcutvsnocut_pt_fake);
}
if(dataType==4){
TH1D* ratio_hPbP_JetIDcut_pt_realvsInc = (TH1D*)hPbP_JetIDcut_pt_real->Clone("ratio_hPbP_JetIDcut_pt_realvsInc");
ratio_hPbP_JetIDcut_pt_realvsInc->Divide(hPbP_JetIDcut_pt);
getRidYError(ratio_hPbP_JetIDcut_pt_realvsInc);

}

//TGraphAsymmErrors *g = new TGraphAsymmErrors();

if(dataType==0) return hdata_JetIDcutvsnocut_pt;
    //g=makegraph(hdata_JetIDcut_pt,hdata_pt);   
else if(dataType==1) return hPbP_JetIDcutvsnocut_pt;
    //g=makegraph(hPbP_JetIDcut_pt,hPbP_pt);
else if(dataType==2) return hPbP_JetIDcutvsnocut_pt_real;
    //g=makegraph(hPbP_JetIDcut_pt_real,hPbP_pt_real);
else if(dataType==3) return hPbP_JetIDcutvsnocut_pt_fake;
    //g=makegraph(hPbP_JetIDcut_pt_fake,hPbP_pt_fake);
else if(dataType==4) return ratio_hPbP_JetIDcut_pt_realvsInc;
else exit();
}


void overlayJetIDcut_PbP_Etabin(){
gStyle->SetOptStat(kFALSE);
c1 = new TCanvas("c1"," ",1200,600);
c2 = new TCanvas("c2"," ",1200,600);
makeMultiPanelCanvas(c1,4,2,0,0,0.25,0.2,0.03);
makeMultiPanelCanvas(c2,4,2,0,0,0.25,0.2,0.03);
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->SetLineColor(0);
hFrame->GetXaxis()->SetLimits(56,599);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Cut Efficiency");
fixedFontHist(hFrame,2.0,3.0);
TLegend *leg1=new TLegend(0.32,0.35,0.72,0.75);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.05);
TLatex *T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextColor(1);
T->SetTextFont(42);
T->SetTextSize(0.065);
for(int i=0;i<Neta;i++){
    c1->cd(canvas[i]+1)->SetGridx();
    if(i==0  || i==6){
        hFrame->GetYaxis()->SetTitle("Ratio");
   //     hFrame->GetXaxis()->SetNdivisions(510);
    }
    else{
 //       hFrame->GetXaxis()->SetNdivisions(505);
    }
    if(i>=4)
        hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
    else
        hFrame->GetXaxis()->SetTitle("");

    hFrame->SetMaximum(1.005);
    hFrame->SetMinimum(0.905);
    hFrame->DrawCopy();
/*TGraphAsymmErrors *g6Data = makeJetIDcut(6,1,i);
g6Data->SetMarkerStyle(20);
g6Data->SetMarkerColor(1);
g6Data->SetLineColor(1);
TGraphAsymmErrors *g6MC = makeJetIDcut(6,0,i);
g6MC->SetMarkerStyle(24);
g6MC->SetMarkerColor(1);
g6MC->SetLineColor(1);
//g6MC->Draw("Psame");
TGraphAsymmErrors *g12Data = makeJetIDcut(12,1,i);
g12Data->SetMarkerStyle(20);
g12Data->SetMarkerColor(1);
g12Data->SetLineColor(1);
TGraphAsymmErrors *g12MC = makeJetIDcut(12,0,i);
g12MC->SetMarkerStyle(24);
g12MC->SetMarkerColor(1);
g12MC->SetLineColor(1);
TGraphAsymmErrors *g14Data = makeJetIDcut(14,1,i);
g14Data->SetMarkerStyle(28);
g14Data->SetMarkerColor(4);
g14Data->SetLineColor(4);
TGraphAsymmErrors *g14MC = makeJetIDcut(14,0,i);
g14MC->SetMarkerStyle(25);
g14MC->SetMarkerColor(4);
g14MC->SetLineColor(4);
*/
TH1D *g21Data = makeJetIDcut(21,0,i);
g21Data->SetMarkerStyle(20);
g21Data->SetMarkerColor(1);
g21Data->SetLineColor(1);
TH1D *g21MC = makeJetIDcut(21,1,i);
g21MC->SetMarkerStyle(24);
g21MC->SetMarkerColor(1);
g21MC->SetLineColor(1);
TH1D *g21MC_real = makeJetIDcut(21,2,i);
g21MC_real->SetMarkerStyle(29);
g21MC_real->SetMarkerColor(4);
g21MC_real->SetLineColor(4);
TH1D *g21MC_fake = makeJetIDcut(21,3,i);
g21MC_fake->SetMarkerStyle(20);
g21MC_fake->SetMarkerColor(2);
g21MC_fake->SetLineColor(2);
TH1D *g21MC_realvsInc = makeJetIDcut(21,4,i);
g21MC_realvsInc->SetMarkerStyle(20);
g21MC_realvsInc->SetMarkerColor(2);
g21MC_realvsInc->SetLineColor(2);
TLine *l =new TLine(30,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
if(canvas[i]!=4){
//g14Data->Draw("Psame");
//g12Data->Draw("Psame");
//g6Data->Draw("Psame");
//g12MC->Draw("Psame");
//g14MC->Draw("Psame");
g21Data->Draw("Psame");
g21MC->Draw("Psame");
g21MC_real->Draw("Psame");
g21MC_realvsInc->Draw("Psame");
l->Draw("same");
T->DrawLatex(0.40,0.30,etastring[i]);
}
if(canvas[i]==4){
//leg1->AddEntry(g6Data,"cut 1 Data","lp");
//leg1->AddEntry(g6MC,"cut 1 PYTHIA+HIJING","lp");
//leg1->AddEntry(g14Data,"cut 2 Data","lp");
//leg1->AddEntry(g14MC,"cut 2 PYTHIA+HIJING","lp");
//leg1->AddEntry(g12Data,"cut 3 Data","lp");
//leg1->AddEntry(g12MC,"cut 3 PYTHIA+HIJING","lp");
leg1->AddEntry(g21Data,"Data","lp");
leg1->AddEntry(g21MC,"MC incluvise cut efficiency","lp");
leg1->AddEntry(g21MC_real,"MC real cut efficiency","lp");
leg1->AddEntry(g21MC_realvsInc,"MC real contribution after cut","lp");
leg1->Draw("same");
}
    c2->cd(canvas[i]+1)->SetGridx();
    if(i==0  || i==6){
        hFrame->GetYaxis()->SetTitle("Data/MC");
    }
    if(i>=4)
        hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
    else
        hFrame->GetXaxis()->SetTitle("");

    hFrame->SetMaximum(1.012);
    hFrame->SetMinimum(0.972);
    hFrame->DrawCopy();
    TH1D* ratio = (TH1D*)g21Data->Clone();
    ratio->Divide(g21MC);
if(canvas[i]==4){
//TLegend* leg2 = new TLegend(*leg1);
//leg2->Clear();
//leg2->AddEntry(ratio,"double ratio Data/MC","lp");
//leg2->Draw("same");
}
if(canvas[i]!=4){
    ratio->Draw("Psame");
    l->Draw("same");
    T->DrawLatex(0.40,0.30,etastring[i]);
}
}
c1->Print(Form("pic/overlay_JetIDcut_PbP_Etabin.png"));
}

