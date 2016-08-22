#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};
const double pPbLumi = 35 ;

void DrawRatio(){
    TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
    TCanvas *c2 = new TCanvas("c2a", "c2",0,0,600,600);
    TCanvas *c3 = new TCanvas("c3a", "c3",0,0,600,600);
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
    c1->SetBottomMargin(0.12);
    c2->SetLeftMargin(0.13);
    c2->SetRightMargin(0.06);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.12);
    c2->SetFrameFillStyle(0);
    c2->SetFrameBorderMode(0);
    c3->SetLeftMargin(0.13);
    c3->SetRightMargin(0.06);
    c3->SetTopMargin(0.05);
    c3->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);


const int N=8;
//const TString etabinname[N]={"-22_-12","-12_-7","-7_-3","-3_3","3_7","7_12","12_22","-10_10"};
const TString etabinname[N]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const double etabin[N]={0.5,0.5,0.5,1.0,0.5,0.5,0.5,2};
const TString etastring[N]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};
TFile *f_Reco[N];
TFile *f_Gen[N];
TH1D *hBay_Reco[N];
TH1D *hBay_Reco_Cl[N];
TH1D *hBay_Reco_Cl1[N];
TH1D *hBay_Reco_Cl2[N];
TH1D *hBay_Gen[N];
ofstream fstr[N];
for(int i=0;i<N;i++){
fstr[i].open(Form("jetaxissys%s.txt",etabinname[i].Data()));
fstr[i]<<setprecision(4)<<fixed;
f_Gen[i] = TFile::Open(Form("/home/xuq/Documents/HI/RpA/output/UnfoldedOutputs/PPb_UnfoPriorGen_akPu3PFKurtForestJetPtEtaWeightNoIDCutGenEtaMatrix_jtpt20_EtaBin%s_Inc_v16.root",etabinname[i].Data()));
f_Reco[i] = TFile::Open(Form("/home/xuq/Documents/HI/RpA/output/UnfoldedOutputs/PPb_UnfoPriorGen_akPu3PFKurtForestEtaWeightNoIDCut_jtpt20_EtaBin%s_Inc_v13.root",etabinname[i].Data()));
hBay_Gen[i] = (TH1D*)f_Gen[i]->Get("hReco0");hBay_Gen[i]->Scale(1./etabin[i]);
hBay_Reco[i] = (TH1D*)f_Reco[i]->Get("hReco0");hBay_Reco[i]->Scale(1./etabin[i]);
}

//-------------------------
c1->cd();
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(47,500);
//hFrame->GetXaxis()->SetLimits(50,600);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetRangeUser(0.9,1.1);
hFrame->GetYaxis()->SetTitle("jet axis uncertainty #eta_{reco}/#eta_{gen}");
fixedFontHist(hFrame,1.2,1.6);
hFrame->DrawCopy();

TLegend* leg1=new TLegend(0.6,0.65,0.78,0.92);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.03);

TLegend* leg2=new TLegend(0.5,0.8,0.68,0.92);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetTextSize(0.03);

TLegend* leg3=new TLegend(0.5,0.68,0.68,0.92);
leg3->SetBorderSize(0);
leg3->SetFillColor(0);
leg3->SetTextSize(0.03);

for(int i=0;i<N;i++){
hBay_Reco_Cl[i]=(TH1D*)hBay_Reco[i]->Clone(Form("RecoRatio_%d",i));
hBay_Reco_Cl[i]->Divide(hBay_Gen[i]);
//fstr[i]<<etastring[i]<<endl<<"Data jet axis"<<endl;
for(int ibin=1;ibin<=hBay_Reco_Cl[i]->GetNbinsX();ibin++){
//hBay_Reco_Cl[i]->SetBinContent(ibin,hBay_Reco_Cl[i]->GetBinContent(ibin));
//fstr[i]<<hBay_Reco_Cl[i]->GetXaxis()->GetBinLowEdge(ibin)<<"to"<<hBay_Reco_Cl[i]->GetXaxis()->GetBinUpEdge(ibin)<<'\t';
if(hBay_Reco_Cl[i]->GetBinCenter(ibin)>50 && hBay_Reco_Cl[i]->GetBinCenter(ibin)<600){
fstr[i]<<hBay_Reco_Cl[i]->GetBinCenter(ibin)<<'\t';
double value=TMath::Abs((hBay_Reco_Cl[i]->GetBinContent(ibin)-1))*100;
fstr[i]<<value<<endl;
}
hBay_Reco_Cl[i]->SetBinError(ibin,1e-10);
}
hBay_Reco_Cl[i]->SetMarkerColor(colorCode[i]);
hBay_Reco_Cl[i]->SetLineColor(colorCode[i]);
hBay_Reco_Cl[i]->SetLineWidth(2);
//hBay_Reco_Cl[i]->SetMarkerStyle(markerCode[i]);
hBay_Reco_Cl[i]->SetMarkerStyle(0);
hBay_Reco_Cl[i]->SetMarkerSize(1.2);
hBay_Reco_Cl1[i]=(TH1D*)hBay_Reco_Cl[i]->Clone(Form("RecodoubleRatio_%d",i));
hBay_Reco_Cl2[i]=(TH1D*)hBay_Reco_Cl[i]->Clone(Form("RecoFMRatio_%d",i));
hBay_Reco_Cl[i]->DrawCopy("Histsame");
leg1->AddEntry(hBay_Reco_Cl[i],etastring[i],"lp");
}
for(int i=0;i<N-1;i++)	hBay_Reco_Cl[i]->DrawCopy("same");
leg1->Draw("same");
drawCMS(0.2,0.85,pPbLumi);

/*TLatex *com0 = new TLatex(0.57,0.90,"PYTHIA+HIJING");
com0->SetTextFont(42);
com0->SetTextSize(0.04);
com0->SetNDC();
com0->Draw();

TLatex *com1 = new TLatex(0.65,0.4,"CMS simulation");
com1->SetTextFont(63);
com1->SetTextSize(17);
com1->SetNDC();
com1->Draw();*/
TLatex *com3 = new TLatex(0.45,0.24,"Anti-k_{T} Particle Flow PU Jets: R=0.3");
com3->SetTextFont(63);
com3->SetTextSize(17);
com3->SetNDC();
com3->Draw();

TLine *l=new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
l->Draw("same");

c1->Print("RecoGensys.pdf");

//----------------------------------------Double Ratio
c2->cd();
const int fweta[3]={2,1,0};
const int bweta[3]={4,5,6};
const int meta[2]={3,7};
hFrame->GetYaxis()->SetTitle("jet axis uncertainty Forward/Backward");
hFrame->DrawCopy();

for(int k=0;k<3;k++){
hBay_Reco_Cl1[fweta[k]]->Divide(hBay_Reco_Cl1[bweta[k]]);
hBay_Reco_Cl1[fweta[k]]->SetMarkerColor(colorCode[fweta[2-k]]);
hBay_Reco_Cl1[fweta[k]]->SetLineColor(colorCode[fweta[2-k]]);
hBay_Reco_Cl1[fweta[k]]->SetLineWidth(2);
//hBay_Reco_Cl1[fweta[k]]->SetMarkerStyle(markerCode[fweta[2-k]]);
hBay_Reco_Cl1[fweta[k]]->SetMarkerStyle(0);
hBay_Reco_Cl1[fweta[k]]->DrawCopy("sameHist");
leg2->AddEntry(hBay_Reco_Cl1[fweta[k]],Form("%s/%s",etastring[fweta[k]].Data(),etastring[bweta[k]].Data()),"lp");
}

leg2->Draw("same");
drawCMS(0.2,0.85,pPbLumi);
com3->Draw();
l->Draw("same");
c2->Print("RecoGendoubleratioFBsys.pdf");

c3->cd();
hFrame->GetYaxis()->SetTitle("jet axis uncertainty Double side/Middle");
hFrame->DrawCopy();

for(int m=0;m<N;m++){
if(m!=meta[1]){
hBay_Reco_Cl2[m]->Divide(hBay_Reco_Cl2[meta[1]]);
hBay_Reco_Cl2[m]->SetMarkerColor(colorCode[m]);
hBay_Reco_Cl2[m]->SetLineColor(colorCode[m]);
hBay_Reco_Cl2[m]->SetLineWidth(2);
//hBay_Reco_Cl2[m]->SetMarkerStyle(markerCode[m]);
hBay_Reco_Cl2[m]->SetMarkerStyle(0);
hBay_Reco_Cl2[m]->DrawCopy("sameHist");
leg3->AddEntry(hBay_Reco_Cl2[m],Form("%s/%s",etastring[m].Data(),etastring[meta[1]].Data()),"lp");
}}

leg3->Draw("same");
drawCMS(0.2,0.85,pPbLumi);
com3->Draw();
l->Draw("same");
c3->Print("RecoGendoubleratioFMsys.pdf");
}

