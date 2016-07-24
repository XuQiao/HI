#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include "par.h"
void Draw(){ 
 TCanvas *c1 = new TCanvas("c1","c1",1,1,550,460);
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
 hist->SetXTitle("Multiplicity");
 hist->SetYTitle(Form("<%s> and systematic errors",str.Data()));
 hist->SetMinimum(-0.001);
 hist->SetMaximum(29.99);
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

TFile *f = TFile::Open("graph.root");
	TGraphErrors* graph = (TGraphErrors*)f->Get(Form("std_%s",str.Data()));
	TGraphErrors* Gri055_graph = (TGraphErrors*)f->Get(Form("Gri055_%s",str.Data()));
	TGraphErrors* Gri101_graph = (TGraphErrors*)f->Get(Form("Gri101_%s",str.Data()));

graph->SetTitle("g1");
graph->SetMarkerStyle(20);
graph->SetMarkerColor(1);
graph->SetLineColor(1);
graph->SetLineWidth(2);
graph->SetMarkerSize(1.2);
graph->Draw("Psameez");

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
//char* label[9]={"60-90%","40-60%","30-40%","20-30%","10-20%","5-10%","1-5%","0-1%","0-90%"};
//char* label[N]={"59.7-100%", "27.4-59.7%","16.6-27.4%","9.41-16.6%", "2.46-9.41%", "0-2.46%"};a
TString label[N];
for(int i=0;i<N;i++)
	if(method==0)label[i] = Form("%.2f-%.2f%%",centbin[N-i-1]*100,centbin[N-i]*100);
	else label[i] = Form("%.2f-%.2f",centbin[N-i-1],centbin[N-i]);

    TLatex *tex1= new TLatex(0.6,26.05,"CMS Preliminary pPb \ #sqrt{s_{NN}} = 5.02 TeV");
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
TLegend *leg0 = new TLegend(0.28,0.63,0.50,0.85);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    leg0->SetTextSize(0.047);
    leg0->AddEntry(graph,"standard","p");
    leg0->AddEntry(Gri055_graph,"Gribov #Omega=0.55","p");
    leg0->AddEntry(Gri101_graph,"Gribov #Omega=1.01","p");
	leg0->Draw();	
c1->SaveAs("Gri_2.png");

}
