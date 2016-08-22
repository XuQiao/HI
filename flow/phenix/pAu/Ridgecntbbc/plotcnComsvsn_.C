#include "fstream.h"
#include "/phenix/u/xuq/util/SimplifyLife.C"

void plotcnComsvsn_(){
gStyle->SetLegendFillColor(0);
gStyle->SetLegendBorderSize(0);
TString type = "ptIn";
double psize = 1.3;
const int color[10] = {1,2,4,6,4,8,6,2,1,7};
const int style[10] = {20,24,27,28,29,33,34,31,22,21};
const int ncent = 9;
const double bbcmean[ncent] = {52.7651, 37.6886, 32.4697, 24.9047, 19.5634, 15.5957, 11.023, 7.34727, 3.46339};
const double centbin[ncent+1] = {0,0.1,0.5,1.0,5.0,10,20,40,60,100};
const double centmidp[ncent] = {0.05,0.3,0.75,3.0,7.5,15,30,50,80};
const int npt = 1;
double ptbin[npt+1] = {1.0,3.0};
double ptmean[npt] = {1.36878};
double c1north[ncent];
double c1south[ncent];
double c1northerr[ncent];
double c1southerr[ncent];
double c2north[ncent];
double c2south[ncent];
double c2northerr[ncent];
double c2southerr[ncent];
double c3north[ncent];
double c3south[ncent];
double c3northerr[ncent];
double c3southerr[ncent];
ifstream fnorth;
ifstream fsouth;
fnorth.open(Form("c1_c2_%s_north.dat",type.Data()));
fsouth.open(Form("c1_c2_%s_south.dat",type.Data()));
for(int icent=0;icent<ncent;icent++){
fnorth>>c1north[icent];
fnorth>>c1northerr[icent];
fnorth>>c2north[icent];
fnorth>>c2northerr[icent];
fnorth>>c3north[icent];
fnorth>>c3northerr[icent];
fsouth>>c1south[icent];
fsouth>>c1southerr[icent];
fsouth>>c2south[icent];
fsouth>>c2southerr[icent];
fsouth>>c3south[icent];
fsouth>>c3southerr[icent];
}

TLatex t;
t.SetNDC();
t.SetTextSize(0.04);
TCanvas *c1 = new TCanvas();
TCanvas *c2 = new TCanvas();
TCanvas *c3 = new TCanvas();
TLegend *leg1 = new TLegend(0.22,0.6,0.42,0.88);
TLegend *leg2 = new TLegend(0.22,0.2,0.42,0.38);
TLegend *leg3 = new TLegend(0.22,0.6,0.42,0.88);
TGraphErrors *gr1north = new TGraphErrors(ncent,bbcmean,c1north,0,c1northerr);
TGraphErrors *gr1south = new TGraphErrors(ncent,bbcmean,c1south,0,c1southerr);
TGraphErrors *gr2north = new TGraphErrors(ncent,bbcmean,c2north,0,c2northerr);
TGraphErrors *gr2south = new TGraphErrors(ncent,bbcmean,c2south,0,c2southerr);
TGraphErrors *gr3north = new TGraphErrors(ncent,bbcmean,c3north,0,c3northerr);
TGraphErrors *gr3south = new TGraphErrors(ncent,bbcmean,c3south,0,c3southerr);
SetStyle(*gr1north,psize,1,20,0,0);
SetStyle(*gr1south,psize,2,24,0,0);
SetStyle(*gr2north,psize,1,20,0,0);
SetStyle(*gr2south,psize,2,24,0,0);
SetStyle(*gr3north,psize,1,20,0,0);
SetStyle(*gr3south,psize,2,24,0,0);
SetRange(*gr1north,bbcmean[0],-0.07,bbcmean[ncent],0.001);
SetRange(*gr1south,bbcmean[0],-0.07,bbcmean[ncent],0.001);
SetRange(*gr2north,bbcmean[0],-0.001,bbcmean[ncent],0.04);
SetRange(*gr2south,bbcmean[0],-0.001,bbcmean[ncent],0.04);
SetRange(*gr3north,bbcmean[0],-0.015,bbcmean[ncent],0.015);
SetRange(*gr3south,bbcmean[0],-0.015,bbcmean[ncent],0.015);
gr1north->SetTitle("");
gr1north->GetXaxis()->SetTitle("Mean hits in vtx");
gr1north->GetYaxis()->SetTitle("c_{1}");
gr1south->SetTitle("");
gr1south->GetXaxis()->SetTitle("Mean hits in vtx");
gr1south->GetYaxis()->SetTitle("c_{1}");
gr2north->SetTitle("");
gr2north->GetXaxis()->SetTitle("Mean hits in vtx");
gr2north->GetYaxis()->SetTitle("c_{2}");
gr2south->SetTitle("");
gr2south->GetXaxis()->SetTitle("Mean hits in vtx");
gr2south->GetYaxis()->SetTitle("c_{2}");
gr3north->SetTitle("");
gr3north->GetXaxis()->SetTitle("Mean hits in vtx");
gr3north->GetYaxis()->SetTitle("c_{3}");
gr3south->SetTitle("");
gr3south->GetXaxis()->SetTitle("Mean hits in vtx");
gr3south->GetYaxis()->SetTitle("c_{3}");
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetFillColor(0);
leg2->SetTextSize(0.04);
leg3->SetFillColor(0);
leg3->SetTextSize(0.04);
leg1->AddEntry(gr1north,Form("north"),"lp");
leg1->AddEntry(gr1south,Form("south"),"lp");
leg2->AddEntry(gr2north,Form("north"),"lp");
leg2->AddEntry(gr2south,Form("south"),"lp");
leg3->AddEntry(gr3north,Form("north"),"lp");
leg3->AddEntry(gr3south,Form("south"),"lp");
c1->cd();
gr1north->Draw("AP");leg1->Draw("same");
gr1south->Draw("Psame");
t.DrawLatex(0.5,0.7,Form("%.1f < p_{T}^{trig} < %.1f (GeV/c)",ptbin[0],ptbin[1]));

c2->cd();
gr2north->Draw("AP");leg3->Draw("same");
gr2south->Draw("Psame");
t.DrawLatex(0.5,0.7,Form("%.1f < p_{T}^{trig} < %.1f (GeV/c)",ptbin[0],ptbin[1]));

c3->cd();
gr3north->Draw("AP");leg3->Draw("same");
gr3south->Draw("Psame");
t.DrawLatex(0.5,0.7,Form("%.1f < p_{T}^{trig} < %.1f (GeV/c)",ptbin[0],ptbin[1]));
	
c1->Print(Form("fig/c1svsn_cent_%s.png",type.Data()));
c2->Print(Form("fig/c2svsn_cent_%s.png",type.Data()));
c3->Print(Form("fig/c3svsn_cent_%s.png",type.Data()));
	
}
