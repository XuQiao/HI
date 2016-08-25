#include "/home/xuq7/HI/jetRpA/fromLinux/Quality/root_setting.h"
#include "/home/xuq7/HI/jetRpA/fromLinux/TreeAna/produceandcheck/file.h"

void DrawUnfoldingClosure_Etabin(){


c1 = new TCanvas("c1"," ",1200,600);
makeMultiPanelCanvas(c1,4,2,0,0,0.25,0.2,0.03);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);   
   /* TCanvas *c1 = new TCanvas("c1a", "c1",0,0,1000,600);
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

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);*/
//const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
//TFile *fout = new TFile("Ratio_Unfolded_Gen.root","recreate");
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetYaxis()->SetRangeUser(0.51,1.49);
fixedFontHist(hFrame,1.3,1.8);
hFrame->GetYaxis()->SetTitleSize(0.06);
hFrame->GetYaxis()->SetLabelSize(0.07);
hFrame->GetXaxis()->SetTitleSize(0.06);
hFrame->GetXaxis()->SetLabelSize(0.07);
hFrame->GetXaxis()->SetNdivisions(505);

TLatex *com0 = new TLatex(0.1,0.92,"PYTHIA+HIJING");
//TLatex *com0 = new TLatex(0.1,0.92,"HIJING+PYTHIA");
com0->SetTextFont(63);
com0->SetTextSize(15);
com0->SetNDC();
TLatex *com2 = new TLatex(0.5,0.32,"CMS simulation");
com2->SetTextFont(63);
com2->SetTextSize(15);
com2->SetNDC();
TLatex *com3 = new TLatex(0.1,0.24,"Anti-k_{T} Particle Flow Jets: R=0.3");
com3->SetTextFont(63);
com3->SetTextSize(15);
com3->SetNDC();

TLatex *T = new TLatex();
T->SetTextFont(42);
T->SetTextSize(0.06);
T->SetNDC();

TLine *l=new TLine(47.5,1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);

TFile *f[Neta];
TFile *f1[Neta];
TH1D *hBay[Neta]; TH1F *hGen_PPb[Neta]; TH1F* hBay_Cl[Neta];TH1F* hGen[Neta];
TLegend* leg=new TLegend(0.3,0.05,0.48,0.35);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.06);
//ofstream fstr[N];
TH1D* hBay[Neta];	TH1F* hMeas[Neta];	TH1F* hBin[Neta];
for(int i=0;i<Neta;i++){
c1->cd(canvas[i]+1);
if(i==0  || i==6){
hFrame->GetYaxis()->SetTitle("Reco / Truth");
hFrame->GetXaxis()->SetNdivisions(510);
hFrame->GetXaxis()->SetLimits(47.5,299);
}
else{
hFrame->GetXaxis()->SetNdivisions(505);
hFrame->GetYaxis()->SetTitle("");
if(i==1 || i==5)
hFrame->GetXaxis()->SetLimits(47.5,360);
else if(i==2 || i==4)
hFrame->GetXaxis()->SetLimits(47.5,420);
else
hFrame->GetXaxis()->SetLimits(47.5,455);
}
if(i>=4)
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
else
hFrame->GetXaxis()->SetTitle("");

hFrame->DrawCopy();
//fstr[i]<<setprecision(4)<<fixed;
//fstr[i].open(Form("jetUnfsys%s.txt",etabinname[i].Data()));

f[i] = TFile::Open(Form("/scratch/xuq7/RpA/UnfoldedOutputs/PPb_UnfoPriorGen_akPu3PFKurtMCNoIDCut_MC_jtpt20_EtaBin%s_Inc_v8.root",etabinnameswap[i].Data()));
//f1[i] = TFile::Open(Form("/home/xuq7/HI/RpA/output/UnfoldedOutputs/PbP_UnfoPriorGen_akPu3PFKurtMCNoIDCut_MCJECv17_jtpt20_GenEtaBin%s_Inc_v17.root",etabinname[i].Data()));
//f[i] = TFile::Open(Form("/scratch/xuq7/RpA/UnfoldedOutputs/PbP_UnfoPriorGen_akPu3PFKurtMCNoIDCut_MCJECv17_jtpt20_EtaBin%s_Inc_v17.root",etabinname[i].Data()));
hBay[i] = (TH1D*)f[i]->Get("hRecoScaled_Cent0");
//hBay[i] = (TH1D*)f[i]->Get("hReco0");
hMeas[i] = (TH1F*)f[i]->Get("hMeasScaled_Cent0");
hBin[i] = (TH1F*)f[i]->Get("hRecoBinByBinScaled_Cent0");
//hGen[i] = (TH1F*)f1[i]->Get("hGen0");

hMeas[i]->SetMarkerStyle(20);
hMeas[i]->SetLineColor(2);
hMeas[i]->SetMarkerColor(2);
hMeas[i]->SetMarkerSize(1.2);

hBay[i]->SetMarkerStyle(25);
hBay[i]->SetLineColor(1);
hBay[i]->SetMarkerColor(1);
hBay[i]->SetMarkerSize(1.2);

hBin[i]->SetMarkerStyle(28);
hBin[i]->SetLineColor(4);
hBin[i]->SetMarkerColor(4);
hBin[i]->SetMarkerSize(1.2);
/*
fout->cd();
hBay[i]->Divide(hGen[i]);
hBay[i]->SetName(Form("Ratio_UnfoldedGen_%s",etabinname[i].Data()));
hBay[i]->Write();
*/
hMeas[i]->DrawCopy("same");
hBay[i]->DrawCopy("same");
hBin[i]->DrawCopy("same");
/*for(int ibin=1;ibin<=hBay_Cl[i]->GetNbinsX();ibin++){
if(hBay_Cl[i]->GetBinCenter(ibin)>50 && hBay_Cl[i]->GetBinCenter(ibin)<600){
fstr[i]<<hBay_Cl[i]->GetBinCenter(ibin)<<'\t';
//fstr[i]<<hBay_Cl[i]->GetXaxis()->GetBinLowEdge(ibin)<<"to"<<hBay_Cl[i]->GetXaxis()->GetBinUpEdge(ibin)<<'\t';
fstr[i]<<100*TMath::Abs((hBay_Cl[i]->GetBinContent(ibin)-1))<<endl;
}
}
*/

if(i==0){
leg->AddEntry(hBay[i],"Bayesian","lp");
leg->AddEntry(hBin[i],"Bin-by-bin","lp");
leg->AddEntry(hMeas[i],"no unfolding","lp");
leg->Draw("same");
}
T->DrawLatex(0.5,0.82,etastring[i]);

l->Draw("same");

}
//for(int i=0;i<N-1;i++) 	hBay_Cl[i]->DrawCopy("same");
com0->Draw();
com2->Draw();
com3->Draw();
//c1->Print(Form("PbPUnfoldingClosure_Etabin.pdf"));
c1->Print(Form("PPbUnfoldingClosure_Etabin.pdf"));
}

