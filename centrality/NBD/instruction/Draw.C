#include "parameter.h"
#include <iostream>
#include <iomanip>
#include "par.h"
void Draw(){ 
 const int nDil = 8;
 const double centDil[nDil+1] = {1.0,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0};
 double centDilbin[nDil];
for(int i=0;i<nDil;i++)
  centDilbin[i]=0.29+i+N-1-nDil;
 const double NpartDil[nDil] = {8.75,30.51,53.30,86.23,130.06,187.35,261.49,355.45};
 const double NpartDilerr[nDil] = {1.13,3.02,3.95,4.35,4.60,4.44,3.96,2.83};
 const double NcollDil[nDil] = {8.01,38.86,86.85,175.76,326.06,563.21,926.79,1484.49};
 const double NcollDilerr[nDil] = {1.41,6.41,12.48,21.13,34.27,52.66,81.37,120.0};
 TGraphErrors *graphNpartDil = new TGraphErrors(nDil,centDilbin,NpartDil,0,NpartDilerr);
 TGraphErrors *graphNcollDil = new TGraphErrors(nDil,centDilbin,NcollDil,0,NcollDilerr);
 TCanvas *c1 = new TCanvas("c1","c1",1,1,550,460);
 c1->SetLogy();
  c1->SetFillColor(10);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderSize(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  gStyle->SetOptStat(0);
  c1->SetTicks(-1);

	N=N-1;
 TString str="Ncoll";
 TH1D* hist = new TH1D("","",N,0,N);
 hist->GetXaxis()->SetNdivisions(502);
if(method==0)
 hist->SetXTitle("Centrality");
else
 hist->SetXTitle("HF #Sigma E_{T} |#eta|>3");
 hist->SetYTitle(Form("<%s> and systematic errors",str.Data()));
 hist->SetMinimum(1);
 hist->SetMaximum(3999.99);
 hist->GetXaxis()->CenterTitle(0);
 hist->GetYaxis()->CenterTitle(1);
 hist->GetYaxis()->SetTitleOffset(1.1);
 hist->GetXaxis()->SetTitleOffset(1.1);
 hist->GetXaxis()->SetTitleSize(0.056);
 hist->GetYaxis()->SetTitleSize(0.056);
 hist->GetXaxis()->SetLabelSize(0.05);
 hist->GetYaxis()->SetLabelSize(0.05);
hist->GetXaxis()->SetLabelOffset(99);
 hist->Draw();
	
	TFile *f=TFile::Open(outG);
	TGraphErrors* graph = (TGraphErrors*)f->Get(Form("std/%s_graph",str.Data()));
//	TGraphErrors* Gri055_graph = (TGraphErrors*)f->Get(Form("Gri055/%s_graph",str.Data()));
//	TGraphErrors* Gri101_graph = (TGraphErrors*)f->Get(Form("Gri101/%s_graph",str.Data()));
        TVectorD *centbin = (TVectorD*)f->Get(Form("std/G0/centbin"));
        TVectorD *kpoint = (TVectorD*)f->Get(Form("std/G0/kpoint"));

graph->SetTitle("g1");
graph->SetMarkerStyle(20);
graph->SetMarkerColor(1);
graph->SetLineColor(1);
graph->SetLineWidth(2);
graph->SetMarkerSize(1.2);
graph->Draw("Psameez");
/*
Gri055_graph->SetTitle("g2");
Gri055_graph->SetMarkerStyle(33);
Gri055_graph->SetMarkerColor(2);
Gri055_graph->SetLineColor(2);
Gri055_graph->SetLineWidth(2);
Gri055_graph->SetMarkerSize(1.2);
Gri055_graph->Draw("Psameez");

Gri101_graph->SetTitle("g3");
Gri101_graph->SetMarkerStyle(34);
Gri101_graph->SetMarkerColor(4);
Gri101_graph->SetLineColor(4);
Gri101_graph->SetLineWidth(2);
Gri101_graph->SetMarkerSize(1.2);
Gri101_graph->Draw("Psameez");
*/
 graphNpartDil->SetMarkerSize(1.2);
 graphNpartDil->SetLineColor(2);
 graphNpartDil->SetLineWidth(2);
 graphNpartDil->SetMarkerStyle(33);
 graphNpartDil->SetMarkerColor(2);
 graphNcollDil->SetMarkerSize(1.2);
 graphNcollDil->SetLineColor(2);
 graphNcollDil->SetLineWidth(2);
 graphNcollDil->SetMarkerStyle(33);
 graphNcollDil->SetMarkerColor(2);
if(str=="Npart") graphNpartDil->Draw("Psame");
if(str=="Ncoll") graphNcollDil->Draw("Psame");
std::vector<TString> label(N);
for(int i=0;i<N;i++)
        if(method==0)label[i] = Form("%.2f-%.2f%%",(*centbin)[i]*100,(*centbin)[i+1]*100);
        else label[i] = Form("%.2f-%.2f",(*kpoint)[i],(*kpoint)[i+1]);

    TLatex *tex1= new TLatex(0.2,0.9,"CMS Preliminary PbPb #sqrt{s_{NN}} = 2.76 TeV");
    tex1->SetNDC();
    tex1->SetTextColor(1);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.05);
    tex1->Draw();

double y = gPad->GetUymin();
// - 0.2*h->GetYaxis()->GetBinWidth(1);
   TText t;
   t.SetTextAngle(45);
   t.SetTextSize(0.03);
   t.SetTextAlign(33);
   for (int i=0;i<N;i++) {
      double x = hist->GetXaxis()->GetBinCenter(i+1);
      t.DrawText(x,y,label[i]);
   }
TLegend *leg0 = new TLegend(0.18,0.70,0.50,0.85);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0);
    leg0->SetTextFont(42);
    leg0->SetTextSize(0.047);
    leg0->AddEntry(graph,"From fitting","p");
//    leg0->AddEntry(Gri055_graph,"Gribov #Omega=0.55","p");
//    leg0->AddEntry(Gri101_graph,"Gribov #Omega=1.01","p");
 if(str=="Npart")   leg0->AddEntry(graphNpartDil,"Npart from run I","p");
 if(str=="Ncoll")   leg0->AddEntry(graphNcollDil,"Ncoll from run I","p");
	leg0->Draw();	
c1->SaveAs(Form("%sGri.png",str.Data()));
c1->SaveAs(Form("%sGri.pdf",str.Data()));

}
