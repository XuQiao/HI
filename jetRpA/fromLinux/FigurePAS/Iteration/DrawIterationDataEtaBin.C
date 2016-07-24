#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

const double pPbLumi = 14.0 ; 
//const double pPbLumi = 20.7 ;
void DrawIterationDataEtaBin(){

    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);   
 /*   TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
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
    c1->SetFrameFillStyle(0);
    c1->SetFrameBorderMode(0);
*/
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
TCanvas* c5 = new TCanvas("c5"," ",1200,600);
makeMultiPanelCanvas(c5,4,2,0,0,0.25,0.20,0.03);

const int N=9;
for(int ieta=0;ieta<Neta;ieta++){
//f = TFile::Open(Form("/home/xuq/Documents/HI/RpA/output/UnfoldedOutputs/PPb_UnfoPriorGen_akPu3PFKurtForestEtaWeightNoIDCut_jtpt20_EtaBin%s_Inc_v13.root",etabinnameswap[ieta].Data()));
f = TFile::Open(Form("/home/xuq/Documents/HI/RpA/output/UnfoldedOutputs/PbP_UnfoPriorGen_akPu3PFKurtForestJetPtEtaWeightNoIDCutJECv17_jtpt20_EtaBin%s_Inc_v17.root",etabinname[ieta].Data()));
TH1D* hReco_cent[N];
TH1D* hReco_cent_Cl[N];
TH1D* hReco0=(TH1D*)f->Get("hReco0");
TH1D* hReco_cent0=(TH1D*)f->Get("hReco_cent0");
TH1D* hUnfoldedJeCsys_cent0=(TH1D*)f->Get("UnfoldedJeCSys_cent0");
for(int i=0;i<N;i++){
//int i=0;
hReco_cent[i]=(TH1D*)f->Get(Form("hRecoRAA_IterSys%d_cent0",i+2));
hReco_cent_Cl[i]= (TH1D*)hReco_cent[i]->Clone(Form("hRecoRAA_IterSys%d_cent0_Cl",i+2));
hReco_cent_Cl[i]->Divide(hReco_cent0);
hReco_cent_Cl[i]->SetMarkerStyle(markerCode[i]);
hReco_cent_Cl[i]->SetMarkerColor(colorCode[i]);
hReco_cent_Cl[i]->SetLineColor(colorCode[i]);
hReco_cent_Cl[i]->SetMarkerSize(1.2);
}
TH1D* hUnfoldedJeCsys_cent0_Cl=(TH1D*)hUnfoldedJeCsys_cent0->Clone("UnfoldedJeCsys_cent0_Cl");
hUnfoldedJeCsys_cent0_Cl->Divide(hReco_cent0);
hUnfoldedJeCsys_cent0_Cl->SetMarkerStyle(20);
hUnfoldedJeCsys_cent0_Cl->SetMarkerColor(1);
hUnfoldedJeCsys_cent0_Cl->SetMarkerSize(1);
TH1F* hFrame=new TH1F("","",1000,0,1000);
if(ieta!=0 && ieta!=6) 
hFrame->GetXaxis()->SetLimits(47.5,475);
else
hFrame->GetXaxis()->SetLimits(47.5,370);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetRangeUser(0.65,1.45);
hFrame->GetYaxis()->SetTitle("Ratio (Unfolded/Nominal)");
fixedFontHist(hFrame,1.4,1.7);
hFrame->GetYaxis()->SetTitleSize(0.06);
hFrame->GetYaxis()->SetLabelSize(0.08);
hFrame->GetXaxis()->SetTitleSize(0.06);
hFrame->GetXaxis()->SetLabelSize(0.08);
hFrame->GetXaxis()->SetNdivisions(505);

c5->cd(canvas[ieta]+1);
hFrame->DrawCopy();

//TLegend* leg1=new TLegend(0.45,0.74,0.62,0.88);
//TLegend* leg2=new TLegend(0.65,0.74,0.82,0.88);
TLegend* leg2=new TLegend(0.40,0.65,0.70,0.95);
TLegend *leg1=new TLegend(0.05,0.65,0.30,0.95);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.06);
leg2->SetTextSize(0.06);
if(ieta==6) continue;
for(int i=0;i<N;i++){
//i=0;
hReco_cent_Cl[i]->DrawCopy("same");
if(i<5)
leg1->AddEntry(hReco_cent_Cl[i],Form("Iteration %d",i+2),"lp");
else
leg2->AddEntry(hReco_cent_Cl[i],Form("Iteration %d",i+2),"lp");
}
for(int i=0;i<N-1;i++) if(i!=7)
hReco_cent_Cl[i]->DrawCopy("same");
//hUnfoldedJeCsys_cent0_Cl->DrawCopy("same");
//leg2->AddEntry(hUnfoldedJeCsys_cent0_Cl,"UnfoldedJeCsys","lp");
if(ieta==0){
drawCMS(0.30,0.86,pPbLumi);
}

if(ieta==3){
//TLatex *com0 = new TLatex(0.17,0.78,"Pb going positive data");
TLatex *com0 = new TLatex(0.10,0.82,"proton going positive data");
com0->SetTextFont(43);
com0->SetTextSize(20);
com0->SetNDC();
com0->Draw();
}

if(ieta==7){
leg1->Draw("same");
leg2->Draw("same");

TLatex *com3 = new TLatex(0.05,0.26,"Anti-k_{T} Particle Flow Jets: R=0.3");
com3->SetTextFont(43);
com3->SetTextSize(17);
com3->SetNDC();
com3->Draw();
}

if(ieta<4)
TLatex *com2 = new TLatex(0.33,0.18,etastring[ieta]);
else
TLatex *com2 = new TLatex(0.33,0.35,etastring[ieta]);
com2->SetTextFont(43);
com2->SetTextSize(20);
com2->SetNDC();

com2->Draw();

TLine *l=new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
l->Draw("same");
}

//c5->Print(Form("PPbIteration_DataEtaBin.pdf"));
c5->Print(Form("PbPIteration_DataEtaBin.pdf"));

}
