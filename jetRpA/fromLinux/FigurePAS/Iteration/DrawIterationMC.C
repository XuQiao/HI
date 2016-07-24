#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

double getMax(double arr[],int N){
double maxv=0;
for(int i=0;i<N;i++){
maxv=TMath::Max(maxv,arr[i]);
}
return maxv;
}

void DrawIterationMC(){

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
const int Nbin_pt=40;
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(47.5,470);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetRangeUser(0.98,1.13);
hFrame->GetYaxis()->SetTitle("Ratio (Unfolded/Nominal)");
fixedFontHist(hFrame,1.2,1.6);
hFrame->DrawCopy();
TFile *f[Neta];
TH1D *hBay[Neta]; TH1F *hGen_PPb[Neta]; TH1F* hBay_Cl[Neta];
TLegend* leg=new TLegend(0.2,0.68,0.38,0.95);
const int N=9;
TH1D* hReco_cent[Neta][N];
TH1D* hReco_cent0[Neta];
TH1D* IterMax[Neta];
TH1D* hReco_cent_Cl[Neta][N];
double IterMaxvalue[Neta][Nbin_pt][N];
//ofstream fstr[Neta];
for(int i=0;i<Neta;i++){
//fstr[i].open(Form("jetsysIter%s.txt",etabinname[i].Data()));
f[i] = TFile::Open(Form("/home/xuq/Documents/HI/RpA/output/UnfoldedOutputs/PPb_UnfoPriorGen_akPu3PFKurtMCNoIDCut_MC_jtpt20_EtaBin%s_Inc_v8.root",etabinname[i].Data()));
//f[i] = TFile::Open(Form("/home/xuq/Documents/HI/RpA/output/UnfoldedOutputs/PbP_UnfoPriorGen_akPu3PFKurtMCNoIDCut_MCJECv17_jtpt20_EtaBin%s_Inc_v17.root",etabinname[i].Data()));
//TH1D* hReco0=(TH1D*)f->Get("hReco0");
hReco_cent0[i]=(TH1D*)f[i]->Get("hReco_cent0");
//TH1D* hUnfoldedJeCsys_cent0=(TH1D*)f->Get("UnfoldedJeCSys_cent0");
for(int j=0;j<N;j++){
//int i=0;
hReco_cent[i][j]=(TH1D*)f[i]->Get(Form("hRecoRAA_IterSys%d_cent0",j+2));
hReco_cent_Cl[i][j]= (TH1D*)hReco_cent[i][j]->Clone(Form("hRecoRAA_IterSys%d_cent0_Cl_%d",j+2,i));
hReco_cent_Cl[i][j]->Divide(hReco_cent0[i]);
}
	IterMax[i]=(TH1D*)hReco_cent_Cl[i][0]->Clone(Form("IterMax_%d",i));
	for(int ibin=1;ibin<=IterMax[i]->GetNbinsX();ibin++){
		for(int j=0;j<N;j++)
		IterMaxvalue[i][ibin][j]=hReco_cent_Cl[i][j]->GetBinContent(ibin);
		IterMax[i]->SetBinContent(ibin,getMax(IterMaxvalue[i][ibin],N));
		IterMax[i]->SetBinError(ibin,1e-10);
	}
/*
	for(int ibin=1;ibin<=IterMax[i]->GetNbinsX();ibin++){
	fstr[i]<<IterMax[i]->GetXaxis()->GetBinLowEdge(ibin)<<"to"<<IterMax[i]->GetXaxis()->GetBinUpEdge(ibin)<<'\t';
		fstr[i]<<IterMax[i]->GetBinContent(ibin)<<endl;
}
*/
}	

TLegend* leg1=new TLegend(0.2,0.64,0.38,0.92);
//TLegend* leg2=new TLegend(0.65,0.70,0.82,0.85);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.035);
for(int i=0;i<Neta;i++){
IterMax[i]->SetMarkerStyle(markerCode[i]);
IterMax[i]->SetMarkerColor(colorCode[i]);
IterMax[i]->SetLineColor(colorCode[i]);
IterMax[i]->SetMarkerSize(1.2);
IterMax[i]->DrawCopy("same");
leg1->AddEntry(IterMax[i],etastring[i],"lp");
}
for(int i=0;i<Neta-1;i++) IterMax[i]->DrawCopy("same");
//hUnfoldedJeCsys_cent0_Cl->DrawCopy("same");
//leg2->AddEntry(hUnfoldedJeCsys_cent0_Cl,"UnfoldedJeCsys","lp");
leg1->Draw("same");

//drawCMS(0.2,0.88,pPbLumi);
TLatex *com2 = new TLatex(0.60,0.88,"CMS simulation");
com2->SetTextFont(42);
com2->SetTextSize(0.04);
com2->SetNDC();
com2->Draw();
TLatex *com0 = new TLatex(0.60,0.83,"PYTHIA+HIJING");
//TLatex *com0 = new TLatex(0.60,0.83,"HIJING+PYTHIA");
com0->SetTextFont(42);
com0->SetTextSize(0.04);
com0->SetNDC();
com0->Draw();
TLatex *com3 = new TLatex(0.30,0.15,"Anti-k_{T} Particle Flow Jets: R=0.3");
com3->SetTextFont(63);
com3->SetTextSize(17);
com3->SetNDC();
com3->Draw();


TLine *l=new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
l->Draw("same");
c1->Print("Iteration_PPb_Etabin.pdf");
//c1->Print("Iteration_PbP_Etabin.pdf");

}
