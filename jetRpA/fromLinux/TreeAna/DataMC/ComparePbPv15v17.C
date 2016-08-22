#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>

void ComparePbPv15v17(int ieta)
{
  bool SavePlots=kTRUE;
  
  TString InputDir = "/home/xuq/Documents/HI/RpA/output/UnfoldedOutputs";

  double EtaMin[8] = {-2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, -1.0 };
  double EtaMax[8] = {-1.5, -1.0, -0.5,  0.5, 1.0, 1.5, 2.0,  1.0 };

  double RevEtaMin[8];
  double RevEtaMax[8];

  for(int j=0; j<8; j++)
  {
    RevEtaMin[j]=-EtaMax[j];
    RevEtaMax[j]=-EtaMin[j];
  }

  double StartPoint=50.;
  double EndPoint;

  if(TMath::Abs(EtaMin[ieta])==TMath::Abs(EtaMax[ieta]))  EndPoint = 430.;
  else if(ieta>0 && ieta<4 )         EndPoint = 300;
  else EndPoint = 250; 

  TFile *f1 = new TFile(Form("%s/PbP_UnfoPriorGen_akPu3PFKurtForestEtaWeightNoIDCut_jtpt20_EtaBin%.f_%.f_Inc_v13.root",InputDir.Data(),10*EtaMin[ieta],10*EtaMax[ieta]));//akPu3PF PbP v15 JEC

  TFile *f2 = new TFile(Form("%s/PbP_UnfoPriorGen_akPu3PFKurtForestJetPtEtaWeightNoIDCutJECv17_jtpt20_EtaBin%.f_%.f_Inc_v17.root",InputDir.Data(),10*EtaMin[ieta],10*EtaMax[ieta]));//akPu3PF PbP v17 JEC
 
//  TH1F *h1eta = (TH1F*)f1->Get("hReco0");//after unfolding
//  TH1F *h2eta = (TH1F*)f2->Get("hReco0");

  TH1F *h1eta = (TH1F*)f1->Get("hReco0");//After unfolding
  TH1F *h2eta = (TH1F*)f2->Get("hReco0");

  double Lumi1 = 14.0*10E6;
  double Lumi2 = 14.0*10E6;
 
  TCanvas *c1 = new TCanvas("c1","c1",0,0,650,650);
  makeMultiPanelCanvas(c1,1,2,0.0,0.0,0.12,0.01,0.06);
  
  gStyle->SetOptStat(0);
  
  c1->cd(1)->SetLogy();
  
  TLegend *t1=new TLegend(0.65,0.65,0.8,0.93);
  t1->SetFillColor(0);
  t1->SetBorderSize(0);
  t1->SetFillStyle(0);
  t1->SetTextFont(63);
  t1->SetTextSize(17);

  TH1F *hDum =new TH1F("hdum","",50,StartPoint,EndPoint);
  hDum->SetMaximum(8.00E-3);
  hDum->SetMinimum(5.00E-9);
  hDum->GetYaxis()->SetTitle("d#sigma/dp_{T}");
  hDum->GetYaxis()->SetTitleSize(0.055);
  hDum->GetYaxis()->SetLabelSize(0.045);
  hDum->GetYaxis()->CenterTitle();
  
  hDum->DrawCopy();

  h1eta->Scale(1./(Lumi1*0.85));
  h1eta->GetXaxis()->SetRangeUser(StartPoint+5,EndPoint);
  h1eta->SetMarkerStyle(24);
  h1eta->SetMarkerColor(kBlue);
  h1eta->SetLineColor(kBlue);
  h1eta->SetMarkerSize(1.2);
  
  h2eta->Scale(1./(0.85*Lumi2));
  h2eta->GetXaxis()->SetRangeUser(StartPoint+5,EndPoint);
  h2eta->SetMarkerStyle(20);
  h2eta->SetMarkerColor(kRed);
  h2eta->SetLineColor(kRed);
  h2eta->SetMarkerSize(1.2);

  h1eta->DrawCopy("same");
  h2eta->DrawCopy("same");

  t1->SetHeader(Form("%2.1f<#eta_{CM}<%2.1f",EtaMin[ieta],EtaMax[ieta]));  
  t1->AddEntry(h1eta,Form("PbP v15 JEC"),"P");
  t1->AddEntry(h2eta,Form("PbP v17 JEC"),"P");
 
  t1->Draw("same");
  c1->Update();


//-----------------------------------ratio-----------------------------------
  c1->cd(2);
  TH1F *ratio1=h1eta->Clone(h1eta->GetName());
  ratio1->Divide(h2eta);
  ratio1->SetMarkerSize(1.2);
  ratio1->SetMarkerStyle(21);
  ratio1->SetMarkerColor(kBlue);
  TH1F *hDum2 = new TH1F("hDum2","",50,StartPoint,EndPoint);
  hDum2->SetMaximum(1.07);
  hDum2->SetMinimum(0.92);
  hDum2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hDum2->GetXaxis()->SetTitleSize(0.055);
  hDum2->GetXaxis()->SetLabelSize(0.045);
  hDum2->GetXaxis()->CenterTitle();

  hDum2->GetYaxis()->SetTitle("v15 JEC/v17 JEC");
  hDum2->GetYaxis()->SetTitleSize(0.055);
  hDum2->GetYaxis()->SetLabelSize(0.045);
  hDum2->GetYaxis()->CenterTitle();

  hDum2->DrawCopy();
  ratio1->DrawCopy("same");

  c1->Update();  

  TLine* l1 = new TLine(StartPoint,1.,EndPoint,1.);
  l1->SetLineWidth(1);
  l1->SetLineStyle(3);
  l1->SetLineColor(kBlack);
  l1->Draw();


  
  if(SavePlots) {
    c1->SaveAs(Form("Comparev15v17JECEta%.f_%.f.pdf",EtaMin[ieta]*10,EtaMax[ieta]*10));
   }


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

