#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};

void FitTruthNtuple(){

    TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
    gStyle->SetOptFit(1111);
  //  gStyle->SetOptStat(1111);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);   
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLeftMargin(0.16);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.1);
    c1->SetFrameFillStyle(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
	c1->SetLogy();
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

const int N=8;
//const TString etabinname[N]={"-22_-12","-12_-7","-7_-3","-3_3","3_7","7_12","12_22","-10_10"};
const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12","-10_10"};
const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
//const TString etastring[N]={"-2.2<#eta<-1.2","-1.2<#eta<-0.7","-0.7<#eta<-0.3","-0.3<#eta<0.3","0.3<#eta<0.7","0.7<#eta<1.2","1.2<#eta<2.2","-1.0<#eta<1.0"};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
const double etalow[N]={-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,-1.0};
const double etaup[N]={-1.2,-0.7,-0.3,0.3,0.7,1.2,-1.0,1.0};
TH1F* hFrame=new TH1F("","",1000,0,1000);
//hFrame->GetXaxis()->SetRangeUser(26,600);
hFrame->GetXaxis()->SetTitle("p_{T}^{ref} (GeV/c)");
hFrame->GetYaxis()->SetRangeUser(1e-11,1e-2);
hFrame->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{T}d#eta}");
fixedFontHist(hFrame,1.2,1.5);
hFrame->DrawCopy();

TFile *f; TF1 *fnlo[N];
TH1F *hGen_PPb[N];
TLegend* leg=new TLegend(0.6,0.45,0.78,0.62);
f = TFile::Open(Form("/scratch/xuq7/RpA/FromYX/PPbMCKurtForestNoVzCut_ppReco_akPu3PF_QCDjetTrig_JetPt0noIPupperCut.root"));
TTree *t = (TTree*)f->Get("nt");
for(int i=0;i<1;i++){
hGen_PPb[i]=new TH1F(Form("hGen_PPb_%d",i),"",1000,0,1000);
t->Draw(Form("refpt>>hGen_PPb_%d",i),Form("weight*((jteta+0.465>=%f) && (jteta+0.465<%f) && rawpt>22 && abs(vz)<15 && refpt >20)",etalow[i],etaup[i]));
//t->Draw(Form("refpt>>hGen_PPb_%d",i),"weight* ((jteta+0.465)>=-2.2 && (jteta+0.465)<3 && rawpt>22 && abs(vz)<15 && refpt >20)");

fnlo[i] = new TF1(Form("fus%d",i),"[0]*exp([1]/x)*pow(x,[2])*pow(1-x*cosh([4])/3500.,[3])", 20., 1000.);
//fnlo[i] = new TF1(Form("fus%d",i),"[0]*exp([1]/x)*pow(x,[2])", 20., 1000.);
fnlo[i]->SetParameters(4.657e6,-41.96,-5.566,2.922,-1.908);
/*fnlo[i]->FixParameter(0,4.657e6);
fnlo[i]->FixParameter(1,-41.96);
fnlo[i]->FixParameter(2,-5.566);
fnlo[i]->FixParameter(3,2.922);
fnlo[i]->FixParameter(4,-1.908);*/
fnlo[i]->SetLineStyle(1);
fnlo[i]->SetLineColor(colorCode[i]);
//hGen_PPb[i] = (TH1F*)f[i]->Get("hGen0");
//hGen_PPb[i]->Scale(1./etabin[i]);
hGen_PPb[i]->GetXaxis()->SetRangeUser(0,600);
hGen_PPb[i]->GetYaxis()->SetRangeUser(1e-11,1e-2);
hGen_PPb[i]->SetMarkerStyle(markerCode[i]);
hGen_PPb[i]->SetMarkerStyle(markerCode[i]);
hGen_PPb[i]->SetLineColor(1);
hGen_PPb[i]->SetMarkerColor(1);
hGen_PPb[i]->SetMarkerSize(1.2);
hGen_PPb[i]->DrawCopy("APsame");
hGen_PPb[i]->Fit(fnlo[i],"","same",30,800);
cout<<"NDF="<<fnlo[i]->GetNDF()<<'\t'<<"Chi2="<<fnlo[i]->GetChisquare()<<endl;
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.035);
leg->AddEntry(hGen_PPb[i],etastring[i],"lp");
//fnlo[i]->Draw("same");
TLatex *com0 = new TLatex(0.2,0.88,"PYTHIA+HIJING");
com0->SetTextFont(42);
com0->SetTextSize(0.04);
com0->SetNDC();
TLatex *com2 = new TLatex(0.65,0.32,"CMS simulation");
com2->SetTextFont(63);
com2->SetTextSize(17);
com2->SetNDC();
TLatex *com3 = new TLatex(0.45,0.24,"Anti-k_{T} Particle Flow PU Jets R=0.3");
com3->SetTextFont(63);
com3->SetTextSize(17);
com3->SetNDC();


TLine *l=new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
l->Draw("same");

}
//for(int i=0;i<N-1;i++) 	hGen_PPb[i]->DrawCopy("same");
com0->Draw();
com2->Draw();
//drawCMS(0.2,0.85);
leg->Draw("same");
com3->Draw();
c1->Print(Form("TruthFit.pdf"));
}

