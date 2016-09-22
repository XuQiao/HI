#include "fstream.h"
#include "/phenix/u/xuq/util/SimplifyLife.C"

void plotcnComsvsn(){
gStyle->SetLegendFillColor(0);
gStyle->SetLegendBorderSize(0);
TString type = "centIn";
double psize = 1.6;
const int color[10] = {1,2,4,6,4,8,6,2,1,7};
const int style[10] = {20,24,27,28,29,33,34,31,22,21};
const int ncent = 1;
const double bbcmean[ncent] = {5.02};
const double centbin[ncent+1] = {0,100};
const double centmidp[ncent] = {50};
const int npt = 10;
const double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
const double ptmean[npt] = {0.360943, 0.691833, 1.1911, 1.69654, 2.20117, 2.70571, 3.2097, 3.71372, 4.21814, 4.72014};
double c1north[npt];
double c1south[npt];
double c1northerr[npt];
double c1southerr[npt];
double c2north[npt];
double c2south[npt];
double c2northerr[npt];
double c2southerr[npt];
double c3north[npt];
double c3south[npt];
double c3northerr[npt];
double c3southerr[npt];
ifstream fnorth;
ifstream fsouth;
fnorth.open(Form("c1_c2_%s_north.dat",type.Data()));
fsouth.open(Form("c1_c2_%s_south.dat",type.Data()));
for(int ipt=0;ipt<npt;ipt++){
fnorth>>c1north[ipt];
fnorth>>c1northerr[ipt];
fnorth>>c2north[ipt];
fnorth>>c2northerr[ipt];
fnorth>>c3north[ipt];
fnorth>>c3northerr[ipt];
fsouth>>c1south[ipt];
fsouth>>c1southerr[ipt];
fsouth>>c2south[ipt];
fsouth>>c2southerr[ipt];
fsouth>>c3south[ipt];
fsouth>>c3southerr[ipt];
}

TLatex t;
t.SetNDC();
t.SetTextSize(0.04);
TCanvas *c1 = new TCanvas();
TCanvas *c2 = new TCanvas();
TCanvas *c3 = new TCanvas();
TLegend *leg1 = new TLegend(0.22,0.2,0.42,0.38);
TLegend *leg2 = new TLegend(0.22,0.2,0.42,0.38);
TLegend *leg3 = new TLegend(0.22,0.2,0.42,0.38);
TGraphErrors *gr1north = new TGraphErrors(npt,ptmean,c1north,0,c1northerr);
TGraphErrors *gr1south = new TGraphErrors(npt,ptmean,c1south,0,c1southerr);
TGraphErrors *gr2north = new TGraphErrors(npt,ptmean,c2north,0,c2northerr);
TGraphErrors *gr2south = new TGraphErrors(npt,ptmean,c2south,0,c2southerr);
TGraphErrors *gr3north = new TGraphErrors(npt,ptmean,c3north,0,c3northerr);
TGraphErrors *gr3south = new TGraphErrors(npt,ptmean,c3south,0,c3southerr);
SetStyle(*gr1north,psize,1,20,0,0);
SetStyle(*gr1south,psize,2,24,0,0);
SetStyle(*gr2north,psize,1,20,0,0);
SetStyle(*gr2south,psize,2,24,0,0);
SetStyle(*gr3north,psize,1,20,0,0);
SetStyle(*gr3south,psize,2,24,0,0);
SetRange(*gr1north,ptbin[0],-0.10,ptbin[npt],0.001);
SetRange(*gr1south,ptbin[0],-0.10,ptbin[npt],0.001);
SetRange(*gr2north,ptbin[0],-0.001,ptbin[npt],0.015);
SetRange(*gr2south,ptbin[0],-0.001,ptbin[npt],0.015);
SetRange(*gr3north,ptbin[0],-0.015,ptbin[npt],0.015);
SetRange(*gr3south,ptbin[0],-0.015,ptbin[npt],0.015);
gr1north->SetTitle("");
gr1north->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr1north->GetYaxis()->SetTitle("c_{1}");
gr1south->SetTitle("");
gr1south->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr1south->GetYaxis()->SetTitle("c_{1}");
gr2north->SetTitle("");
gr2north->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr2north->GetYaxis()->SetTitle("c_{2}");
gr2south->SetTitle("");
gr2south->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr2south->GetYaxis()->SetTitle("c_{2}");
gr3north->SetTitle("");
gr3north->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr3north->GetYaxis()->SetTitle("c_{3}");
gr3south->SetTitle("");
gr3south->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr3south->GetYaxis()->SetTitle("c_{3}");
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetFillColor(0);
leg2->SetTextSize(0.04);
leg3->SetFillColor(0);
leg3->SetTextSize(0.04);
leg1->AddEntry(gr1north,Form("1.5<|#eta_{asso}|<3.0"),"lp");
leg1->AddEntry(gr1south,Form("1.0<|#eta_{asso}|<3.0"),"lp");
leg2->AddEntry(gr2north,Form("1.5<|#eta_{asso}|<3.0"),"lp");
leg2->AddEntry(gr2south,Form("1.0<|#eta_{asso}|<3.0"),"lp");
leg3->AddEntry(gr3north,Form("1.5<|#eta_{asso}|<3.0"),"lp");
leg3->AddEntry(gr3south,Form("1.0<|#eta_{asso}|<3.0"),"lp");
c1->cd();
gr1north->Draw("AP");leg1->Draw("same");
gr1south->Draw("Psame");

c2->cd();
gr2north->Draw("AP");leg3->Draw("same");
gr2south->Draw("Psame");

c3->cd();
gr3north->Draw("AP");leg3->Draw("same");
gr3south->Draw("Psame");
	
c1->Print(Form("fig/c1svsn_pt_%s.png",type.Data()));
c2->Print(Form("fig/c2svsn_pt_%s.png",type.Data()));
c3->Print(Form("fig/c3svsn_pt_%s.png",type.Data()));
	
}
