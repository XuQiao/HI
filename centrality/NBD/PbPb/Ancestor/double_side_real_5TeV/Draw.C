#include "../parameter.h"
#include <iostream>
#include <iomanip>
#include "par.h"

void Draw(){ 
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
 //TString str="Npart";
 TString str="Ncoll";
 TH1D* hist = new TH1D("","",N,0,N);
 hist->GetXaxis()->SetNdivisions(510);
if(method==0)
 hist->SetXTitle("(100 - Centrality(%))*2");
else
 hist->SetXTitle("HF #Sigma E_{T} |#eta|>3");
 //hist->SetYTitle(Form("<%s> and systematic errors",str.Data()));
 hist->SetYTitle(Form("<%s> and RMS",str.Data()));
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
// hist->GetXaxis()->SetLabelOffset(99);
 hist->Draw();
	
	TFile *f=TFile::Open(outG);
	TGraphErrors* graph = (TGraphErrors*)f->Get(Form("std/%s_graph",str.Data()));
//	TGraphErrors* Gri055_graph = (TGraphErrors*)f->Get(Form("Gri055/%s_graph",str.Data()));
//	TGraphErrors* Gri101_graph = (TGraphErrors*)f->Get(Form("Gri101/%s_graph",str.Data()));
        TVectorD *centbin = (TVectorD*)f->Get(Form("std/G0/centbin"));
        TVectorD *kpoint = (TVectorD*)f->Get(Form("std/G0/kpoint"));
	TVectorD *NcollAver = (TVectorD*)f->Get(Form("std/G0/NcollAver"));
	TVectorD *NpartAver = (TVectorD*)f->Get(Form("std/G0/NpartAver"));
	TVectorD *BAver = (TVectorD*)f->Get(Form("std/G0/BAver"));
	TVectorD *NcollRMS = (TVectorD*)f->Get(Form("std/G0/NcollRMS"));
	TVectorD *NpartRMS = (TVectorD*)f->Get(Form("std/G0/NpartRMS"));
	TVectorD *BRMS = (TVectorD*)f->Get(Form("std/G0/BRMS"));
        TVectorD centbin100 = *centbin;
        TVectorD centbin100reverse = *centbin;
        for(int ibin=0;ibin<N+1;ibin++){
            centbin100reverse[ibin] = centbin100[N-ibin]*100;
        }
        double *acentbin100 = centbin100reverse.GetMatrixArray();
        TH1F* hNcoll = new TH1F("hNcoll","",N,acentbin100);
        TH1F* hNpart = new TH1F("hNpart","",N,acentbin100);
        TH1F* hB = new TH1F("hB","",N,acentbin100);
        for(int ibin=0;ibin<N;ibin++){
        hNcoll->SetBinContent(ibin+1,(*NcollAver)[N-1-ibin]);
        hNcoll->SetBinError(ibin+1,(*NcollRMS)[N-1-ibin]);
        hNpart->SetBinContent(ibin+1,(*NpartAver)[N-1-ibin]);
        hNpart->SetBinError(ibin+1,(*NpartRMS)[N-1-ibin]);
        hB->SetBinContent(ibin+1,(*BAver)[N-1-ibin]);
        hB->SetBinError(ibin+1,(*BRMS)[N-1-ibin]);
        }
        TFile *fout = new TFile("Table_GlauberNBD_v0.root","recreate");
        fout->cd();
        hNcoll->Write();
        hNpart->Write();
        hB->Write();
        fout->Close();

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
std::vector<TString> label(N);
for(int i=0;i<N;i++)
        if(method==0)label[i] = Form("%.2f-%.2f%%",(*centbin)[i]*100,(*centbin)[i+1]*100);
        else label[i] = Form("%.2f-%.2f",(*kpoint)[i],(*kpoint)[i+1]);

    TLatex *tex1= new TLatex(0.2,0.9,"CMS Preliminary PbPb #sqrt{s_{NN}} = 5 TeV");
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
   //   t.DrawText(x,y,label[i]);
   }
TLegend *leg0 = new TLegend(0.18,0.70,0.50,0.85);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0);
    leg0->SetTextFont(42);
    leg0->SetTextSize(0.047);
    leg0->AddEntry(graph,"From Glauber-NBD fitting","p");
//    leg0->AddEntry(Gri055_graph,"Gribov #Omega=0.55","p");
//    leg0->AddEntry(Gri101_graph,"Gribov #Omega=1.01","p");
//	leg0->Draw();	
c1->SaveAs(Form("%sGri.png",str.Data()));
c1->SaveAs(Form("%sGri.pdf",str.Data()));

}
