#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include "/home/xuq7/HI/jetRpA/RpA/NewMC/produceandcheck/file.h"

void CompareetaSpectrapPbPbp()
{
  bool SavePlots=kTRUE;
   double etacut[]={-3, -2.5, -2.172, -1.93, -1.74, -1.566, -1.392, -1.218,-1.044, -0.879, -0.696, -0.522, -0.348, -0.174, 0, 0.174, 0.348, 0.522, 0.696, 0.879, 1.044, 1.218, 1.392, 1.566, 1.74, 1.93, 2.172, 2.5, 3.0};
   int netabin = sizeof(etacut)/sizeof(double)-1;
 int ptmin=50, ptmax=600;
 TTree *nt1 = (TTree*)fMCPPbjetloop->Get("nt");
 TTree *nt2 = (TTree*)fMCPbPjetloop->Get("nt");
 TH2F *h1pteta = new TH2F("h1pteta","h1pteta",1000,0,1000,netabin,etacut);
 TH2F *h2pteta = new TH2F("h2pteta","h2pteta",1000,0,1000,netabin,etacut);
 h1pteta->Sumw2();
 h2pteta->Sumw2();
 nt1->Draw("refeta:refpt>>h1pteta","xSecWeight*(TMath::Abs(jteta)<3 && rawpt>20 && refpt > 20)");
 nt2->Draw("refeta:refpt>>h2pteta","xSecWeight*(TMath::Abs(jteta)<3 && rawpt>20 && refpt > 20)");
TFile *f1 = new TFile("1.root","recreate");
 TH1F *h1eta = (TH1F*)h1pteta->ProjectionY("h1eta",h1pteta->GetXaxis()->FindBin(ptmin),h1pteta->GetXaxis()->FindBin(ptmax));
 TH1F *h2eta = (TH1F*)h2pteta->ProjectionY("h2eta",h2pteta->GetXaxis()->FindBin(ptmin),h2pteta->GetXaxis()->FindBin(ptmax));
// h1eta=(TH1F*)h1eta->Rebin(4);
// h2eta=(TH1F*)h2eta->Rebin(4);
// h1eta->Scale(1./h1eta->Integral());
// h2eta->Scale(1./h2eta->Integral());
 h1eta->Scale(1./(ptmax-ptmin));
 h2eta->Scale(1./(ptmax-ptmin));
TGraphErrors *g1,*g2;
double a1eta[1000],a1y[1000],a2eta[1000],a2y[1000],a1yerr[1000],a2yerr[1000];
 NormalizeByBinWidth(h1eta);
 NormalizeByBinWidth(h2eta);
int Nbin = h1eta->GetNbinsX();

 for(int ibin=0;ibin<Nbin;ibin++){
  a1y[ibin] = h1eta->GetBinContent(ibin);
  a1yerr[ibin] = h1eta->GetBinError(ibin);
  a2y[ibin] = h2eta->GetBinContent(ibin);
  a2yerr[ibin] = h2eta->GetBinError(ibin);
  a1eta[ibin]=h1eta->GetBinCenter(ibin);
  a2eta[ibin]=h2eta->GetBinCenter(ibin);
  a1eta[ibin]+=0.465;
  a1eta[ibin]=-a1eta[ibin];
  a2eta[ibin]=h2eta->GetBinCenter(ibin);
  a2eta[ibin]-=0.465;
}
  
  g1= new TGraphErrors(Nbin,a1eta,a1y,0,a1yerr);
  g2= new TGraphErrors(Nbin,a2eta,a2y,0,a2yerr);
  TCanvas *c1 = new TCanvas("c1","c1",0,0,650,650);
  TCanvas *c2 = new TCanvas("c2","c2",0,0,650,650);
  makeMultiPanelCanvas(c1,1,1,0.0,0.0,0.12,0.05,0.06);
  makeMultiPanelCanvas(c2,1,1,0.0,0.0,0.12,0.05,0.06);
  f1->cd();
  gStyle->SetOptStat(0);
  
  //c1->cd(1)->SetLogy();
  c1->cd(1);
  
  TLegend *t1=new TLegend(0.35,0.75,0.8,0.93);
  t1->SetFillColor(0);
  t1->SetBorderSize(0);
  t1->SetFillStyle(0);
  t1->SetTextFont(63);
  t1->SetTextSize(17);

  TH1F *hDum =new TH1F("hdum","",200,-3,3);
  //hDum->SetMaximum(2.00E-4);
  //hDum->SetMaximum(2.00E6);
  hDum->SetMaximum(5e-6);
  hDum->SetMinimum(0);
  hDum->GetXaxis()->SetRangeUser(-3,3);
  //hDum->GetXaxis()->SetTitle("generator #eta_{lab}");
  hDum->GetXaxis()->SetTitle("jet #eta_{CM}");
  //hDum->GetYaxis()->SetTitle("Event Fraction");
  hDum->GetYaxis()->SetTitle("d#sigma/dp_{T}");
  hDum->GetYaxis()->SetTitleSize(0.045);
  hDum->GetYaxis()->SetTitleOffset(1.2);
  hDum->GetYaxis()->SetLabelSize(0.025);
  hDum->GetYaxis()->CenterTitle();
  hDum->GetXaxis()->CenterTitle();
  
  hDum->DrawCopy();
//  h1eta->Scale(1./Lumi1);
//  h2eta->Scale(1./Lumi2);

  g1->SetMarkerStyle(24);
  g1->SetMarkerColor(kBlue);
  g1->SetLineColor(kBlue);
  g1->SetMarkerSize(1.2);
 
  g2->SetMarkerStyle(20);
  g2->SetMarkerColor(kRed);
  g2->SetLineColor(kRed);
  g2->SetMarkerSize(1.2);

  g1->Draw("Psame");
  g2->Draw("Psame");
  
  t1->AddEntry(g1,Form("PPb (%d < p_{T} < %d (GeV/c))",ptmin,ptmax),"p");
  t1->AddEntry(g2,Form("PbP (%d < p_{T} < %d (GeV/c))",ptmin,ptmax),"p");
 
  t1->Draw("same");
  //c1->Update();


//-----------------------------------ratio-----------------------------------
  TGraphErrors *ratio;
  double ratioeta[1000], ratioy[1000];
  c2->cd(1);
  int Nbin_=0;
  for(int ibin=0;ibin<Nbin;ibin++){
        for(int jbin=0;jbin<Nbin;jbin++){
        if(a1eta[ibin] == a2eta[jbin]){
                ratioeta[Nbin_] = a1eta[ibin];
                ratioy[Nbin_] = a1y[ibin]/a2y[jbin];
                Nbin_++;
            }
        }
 }
  ratio= new TGraphErrors(Nbin_,ratioeta,ratioy,0,0);
  ratio->SetTitle("PPb/PbP");
  f1->cd();
// TH1F *ratio1=h1eta->Clone(h1eta->GetName());
//  ratio1->Divide(h2eta);
  ratio->SetMarkerSize(1.2);
  ratio->SetMarkerStyle(21);
  ratio->SetMarkerColor(kBlue);

  TH1F *hDum2 = new TH1F("hDum2","",netabin,etacut);
  hDum2->SetMaximum(1.20);
  hDum2->SetMinimum(0.80);
  //hDum2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hDum->GetXaxis()->SetTitle("jet #eta_{CM}");
  hDum2->GetXaxis()->SetTitleSize(0.055);
  hDum2->GetXaxis()->SetLabelSize(0.045);
  hDum2->GetXaxis()->CenterTitle();

  hDum2->GetYaxis()->SetTitle("pPb/Pbp");
  hDum2->GetYaxis()->SetTitleSize(0.055);
  hDum2->GetYaxis()->SetLabelSize(0.045);
  hDum2->GetYaxis()->CenterTitle();

  hDum2->DrawCopy();
  ratio->Draw("Psame");

 // c2->Update();  

  TLine* l1 = new TLine(-3,1.,3,1.);
  l1->SetLineWidth(1);
  l1->SetLineStyle(3);
  l1->SetLineColor(kBlack);
  l1->Draw();


  
  if(SavePlots) {
    c1->SaveAs(Form("CompareetaSpectrapPb_PbpLabEta_%d_%d.pdf",ptmin,ptmax));
    c2->SaveAs(Form("RatioCompareetaSpectrapPb_PbpLabEta_%d_%d.pdf",ptmin,ptmax));
   }
f1->cd();
h1pteta->Write();
h2pteta->Write();
 g2->Write("gPbPrefeta"); 
 g1->Write("gPPbrefeta"); 
ratio->Write("ratiogeneta");

f1->Close();

}


void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge)
{
  if(canv == 0){
    Error("makeMultiPanelCanvas","GotNullCanvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];
  Float_t Xlow[columns];
   Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
  Float_t PadWidth =
  (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
  (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
  (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
  (1.0/(1.0-edge))+(Float_t)rows-2.0);
  Xlow[0] = leftOffset;
  Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
    }
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
      Xlow[i],Ylow[j],Xup[i],Yup[j]);
      if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
        else pad[i][j]->SetLeftMargin(0);
      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
        else pad[i][j]->SetRightMargin(0);
      if(j==0) pad[i][j]->SetTopMargin(edge);
        else pad[i][j]->SetTopMargin(0);
      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin+0.3);
        else pad[i][j]->SetBottomMargin(0);
      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
  }
 }
}


TH1* NormalizeByBinWidth(TH1 *histo){
  for(int i=1; i<=histo->GetNbinsX();i++){
    double content = histo->GetBinContent(i);
    double error = histo->GetBinError(i);
    histo->SetBinContent(i,content/histo->GetBinWidth(i));
    histo->SetBinError(i,error/histo->GetBinWidth(i));
    }
    return histo;
  }

