#include "../parameter.h"
#include <iostream>
#include <iomanip>
#include "par.h"
void Drawmethodcomp(){ 
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

	TFile *f1=TFile::Open(Form("../../OneComp/double_side/%s",outG.Data()));
	TGraphErrors* graph1 = (TGraphErrors*)f1->Get(Form("std/%s_graph",str.Data()));
//	TGraphErrors* Gri055_graph = (TGraphErrors*)f->Get(Form("Gri055/%s_graph",str.Data()));
//	TGraphErrors* Gri101_graph = (TGraphErrors*)f->Get(Form("Gri101/%s_graph",str.Data()));
        TVectorD *centbin1 = (TVectorD*)f1->Get(Form("std/G0/centbin"));
        TVectorD *kpoint1 = (TVectorD*)f1->Get(Form("std/G0/kpoint"));

graph->SetTitle("g1");
graph->SetMarkerStyle(20);
graph->SetMarkerColor(1);
graph->SetLineColor(1);
graph->SetLineWidth(2);
graph->SetMarkerSize(1.2);
graph->Draw("Psameez");

graph1->SetMarkerStyle(24);
graph1->SetMarkerColor(4);
graph1->SetLineColor(4);
graph1->SetLineWidth(2);
graph1->SetMarkerSize(1.2);
graph1->Draw("Psameez");
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
    leg0->AddEntry(graph,"From Ancestor fitting","p");
    leg0->AddEntry(graph1,"From One Comp fitting","p");
//    leg0->AddEntry(Gri055_graph,"Gribov #Omega=0.55","p");
//    leg0->AddEntry(Gri101_graph,"Gribov #Omega=1.01","p");
 if(str=="Npart")   leg0->AddEntry(graphNpartDil,"Npart from run I","p");
 if(str=="Ncoll")   leg0->AddEntry(graphNcollDil,"Ncoll from run I","p");
	leg0->Draw();	
c1->SaveAs(Form("%sGri.png",str.Data()));
c1->SaveAs(Form("%sGri.pdf",str.Data()));


 TCanvas *c2 = new TCanvas("c2","c2",1,1,550,460);
  c2->SetFillColor(10);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderSize(0);
  c2->SetFrameBorderMode(0);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.15);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.02);
  c2->SetTicks(-1);
c2->cd();
TGraphErrors* fdvf1 = (TGraphErrors*)graph->Clone("fdvf1");
TGraphErrors* f1dvDil = (TGraphErrors*)graph->Clone("f1dvDil");
TGraphErrors* fdvDil = (TGraphErrors*)graph->Clone("fdvDil");
TGraphErrors* Dilsys = (TGraphErrors*)graph->Clone("Dilsys");
TGraphErrors* Dilsys2 = (TGraphErrors*)graph->Clone("Dilsys2");
for(int ip = 0;ip<fdvf1->GetN();ip++){
    double x = graph->GetX()[ip];
    double ey = graph->GetEY()[ip];
    double y = graph->GetY()[ip];
    double ey1 = graph1->GetEY()[ip];
    double y1 = graph1->GetY()[ip];
    double eyNcoll = graphNcollDil->GetEY()[ip];
    double yNcoll = graphNcollDil->GetY()[ip];
    double eyNpart = graphNpartDil->GetEY()[ip];
    double yNpart = graphNpartDil->GetY()[ip];
    if(y1!=0 && y!=0){
        fdvf1->SetPoint(ip,x,y/y1);
        fdvf1->SetPointError(ip,0,y/y1*sqrt((ey/y)**2+(ey1/y1)**2));
        if(str=="Ncoll"){
        fdvDil->SetPoint(ip,x,y/yNcoll);
        double yerr = y/yNcoll*sqrt((ey/y)**2+(eyNcoll/yNcoll)**2);
        fdvDil->SetPointError(ip,0,0);
        f1dvDil->SetPoint(ip,x,y1/yNcoll);
        double y1err = y1/yNcoll*sqrt((ey1/y1)**2+(eyNcoll/yNcoll)**2);
        f1dvDil->SetPointError(ip,0,0);
        Dilsys->SetPoint(ip,x,1.);
        Dilsys->SetPointError(ip,0,TMath::Max(fabs(y1/yNcoll-1),fabs(y/yNcoll-1)));
        Dilsys2->SetPoint(ip,x,1.);
        Dilsys2->SetPointError(ip,0,TMath::Max(yerr,y1err));
        }
        else if(str=="Npart"){
        fdvDil->SetPoint(ip,x,y/yNpart);
        double yerr = y/yNpart*sqrt((ey/y)**2+(eyNcoll/yNpart)**2);
        fdvDil->SetPointError(ip,0,0);
        f1dvDil->SetPoint(ip,x,y1/yNpart);
        double y1err = y1/yNpart*sqrt((ey1/y1)**2+(eyNcoll/yNcoll)**2);
        f1dvDil->SetPointError(ip,0,0);
        Dilsys->SetPoint(ip,x,1.);
        Dilsys->SetPointError(ip,0,TMath::Max(fabs(y1/yNpart-1),fabs(y/yNpart-1)));
        Dilsys2->SetPoint(ip,x,1.);
        Dilsys2->SetPointError(ip,0,TMath::Max(yerr,y1err));
        }
}
}
hist->SetMinimum(0.5);
hist->SetMaximum(1.5);
hist->SetYTitle(Form("<%s> ratio",str.Data()));
hist->Draw();
fdvf1->SetMarkerStyle(20);
fdvf1->SetMarkerColor(1);
fdvf1->SetLineColor(1);
fdvf1->SetLineWidth(2);
fdvf1->SetMarkerSize(1.4);
//fdvf1->Draw("Psameez");
f1dvDil->SetMarkerStyle(24);
f1dvDil->SetMarkerColor(2);
f1dvDil->SetLineColor(2);
f1dvDil->SetLineWidth(2);
f1dvDil->SetMarkerSize(1.4);
fdvDil->SetMarkerStyle(20);
fdvDil->SetMarkerColor(4);
fdvDil->SetLineColor(4);
fdvDil->SetLineWidth(2);
fdvDil->SetMarkerSize(1.4);
Dilsys->SetFillColor(kYellow+2);
Dilsys2->SetFillColor(kYellow+1);
//Dilsys2->Draw("e3same");
Dilsys->Draw("e3same");
fdvDil->Draw("Psameez");
f1dvDil->Draw("Psameez");
tex1->Draw();
TLine *l = new TLine(0,1,N,1);
l->SetLineStyle(2);
l->SetLineWidth(3);
l->Draw("same");
TLegend *leg1 = new TLegend(0.38,0.72,0.75,0.86);
    leg1->SetFillColor(10);
    leg1->SetBorderSize(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.04);
    //leg1->AddEntry(fdvf1,"Ancestor fitting/One Comp fitting","p");
    leg1->AddEntry(f1dvDil,"One Component fitting/Run I fitting","p");
    leg1->AddEntry(fdvDil,"Two Component fitting/Run I fitting","p");
    leg1->AddEntry(Dilsys,"RunI fitting systematics","f");
  //  leg1->AddEntry(Dilsys2,"RunI fitting ratio systematics","f");
    double y = gPad->GetUymin()+0.5;
   for (int i=0;i<N;i++) {
      double x = hist->GetXaxis()->GetBinCenter(i+1);
      t.DrawText(x,y,label[i]);
   }
   leg1->Draw();
c2->SaveAs(Form("%sGriratio.png",str.Data()));
c2->SaveAs(Form("%sGriratio.pdf",str.Data()));
}
