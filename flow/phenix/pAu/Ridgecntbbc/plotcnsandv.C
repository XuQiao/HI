#include "fstream.h"
#include "/phenix/u/xuq/util/SimplifyLife.C"

void plotcnsandv(string t_){
gStyle->SetLegendFillColor(0);
gStyle->SetLegendBorderSize(0);
TString type = t_.c_str();
double psize = 1.4;
const int color[10] = {1,2,4,6,4,8,6,2,1,7};
const int style[10] = {20,24,27,28,29,33,34,31,22,21};
if(type.Contains("pt") && type != "ptccentc"){
const int ncent = 9;
const double bbcmean[ncent] = {30.8494, 23.7815, 20.3412, 15.4924, 11.6095, 8.95753,5.9812,3.9593, 1.85275};
const double centbin[ncent+1] = {0,0.1,0.5,1.0,5.0,10,20,40,60,100};
const double centmidp[ncent] = {0.05,0.3,0.75,3.0,7.5,15,30,50,80};
if(type=="ptIn"){//1.0-3.0
const int npt = 1;
double ptbin[npt+1] = {1.0,3.0};
double ptmean[npt] = {1.36878};
}
else if(type=="ptIn25_4"){
const int npt = 1;
double ptbin[npt+1] = {2.5,4.0};
double ptmean[npt] = {2.92489};
}
else if(type=="ptcoarser"){
const int npt = 4;
double ptbin[npt+1] = {0.2,1.0,2.0,3.0,5.0};
const double ptmean[npt] = {0.519639, 1.29345, 2.32523, 3.51803};
}
else if(type=="ptfiner"){
const int npt = 10;
const double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
const double ptmean[npt] = {0.360943, 0.691833, 1.1911, 1.69654, 2.20117, 2.70571, 3.2097, 3.71372, 4.21814, 4.72014};
}
}
else if(type=="centIn"){
const int ncent = 1;
const double bbcmean[ncent] = {12.0546};
const double centbin[ncent+1] = {0,100};
const double centmidp[ncent] = {50};
const int npt = 10;
const double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
const double ptmean[npt] = {0.360943, 0.691833, 1.1911, 1.69654, 2.20117, 2.70571, 3.2097, 3.71372, 4.21814, 4.72014};
}
else if(type=="ptccentc"){
const int ncent = 4;
const double bbcmean[ncent] = {30.8494, 23.7815, 20.3412, 15.4924};
const double centbin[ncent+1] = {0,0.5,1,10,100};
const double centmidp[ncent] = {0.25,0.75,5.5,55};
const int npt = 5;
const double ptbin[npt+1] = {0.2,0.5,1.0,2.0,3.0,4.0};
const double ptmean[npt] = {0.360943, 0.691833, 1.29345, 2.32523, 3.3542};
}

double c1[ncent][npt];
double c1err[ncent][npt];
double c2[ncent][npt];
double c2err[ncent][npt];
double c3[ncent][npt];
double c3err[ncent][npt];
double c1_[npt][ncent];
double c1err_[npt][ncent];
double c2_[npt][ncent];
double c2err_[npt][ncent];
double c3_[npt][ncent];
double c3err_[npt][ncent];
ifstream f;
f.open(Form("c1_c2_%s.dat",type.Data()));
for(int icent=0;icent<ncent;icent++){
for(int ipt=0;ipt<npt;ipt++){
f>>c1[icent][ipt];
f>>c1err[icent][ipt];
f>>c2[icent][ipt];
f>>c2err[icent][ipt];
f>>c3[icent][ipt];
f>>c3err[icent][ipt];
c1_[ipt][icent] = c1[icent][ipt];
c1err_[ipt][icent] = c1err[icent][ipt];
c2_[ipt][icent] = c2[icent][ipt];
c2err_[ipt][icent] = c2err[icent][ipt];
c3_[ipt][icent] = c3[icent][ipt];
c3err_[ipt][icent] = c3err[icent][ipt];
}
}

TGraphErrors *gr1[npt>ncent?npt:ncent];
TGraphErrors *gr2[npt>ncent?npt:ncent];
TGraphErrors *gr3[npt>ncent?npt:ncent];
TCanvas *c4 = new TCanvas();
TCanvas *c5 = new TCanvas();
TCanvas *c6 = new TCanvas();
TLatex t;
t.SetNDC();
t.SetTextSize(0.04);
TLegend *leg1 = new TLegend(0.62,0.82-npt*0.03,0.82,0.88);
TLegend *leg2 = new TLegend(0.62,0.82-npt*0.03,0.82,0.88);
TLegend *leg3 = new TLegend(0.22,0.82-npt*0.03,0.42,0.88);
for(int ipt = 0;ipt<npt; ipt++){
gr1[ipt] = new TGraphErrors(ncent,bbcmean,c1_[ipt],0,c1err_[ipt]);
gr2[ipt] = new TGraphErrors(ncent,bbcmean,c2_[ipt],0,c2err_[ipt]);
gr3[ipt] = new TGraphErrors(ncent,bbcmean,c3_[ipt],0,c3err_[ipt]);
SetStyle(*gr1[ipt],psize,color[ipt],style[ipt],0,0);
SetStyle(*gr2[ipt],psize,color[ipt],style[ipt],0,0);
SetStyle(*gr3[ipt],psize,color[ipt],style[ipt],0,0);
SetRange(*gr1[ipt],centbin[0],-0.07,centbin[ncent],0.001);
SetRange(*gr2[ipt],centbin[0],-0.001,centbin[ncent],0.02);
SetRange(*gr3[ipt],centbin[0],-0.015,centbin[ncent],0.015);
gr1[ipt]->SetTitle("");
gr1[ipt]->GetXaxis()->SetTitle("Mean BBc charge south");
gr1[ipt]->GetYaxis()->SetTitle("c_{1}");
gr2[ipt]->SetTitle("");
gr2[ipt]->GetXaxis()->SetTitle("Mean BBc charge south");
gr2[ipt]->GetYaxis()->SetTitle("c_{2}");
gr3[ipt]->SetTitle("");
gr3[ipt]->GetXaxis()->SetTitle("Mean BBc charge south");
gr3[ipt]->GetYaxis()->SetTitle("c_{3}");
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetFillColor(0);
leg2->SetTextSize(0.04);
leg3->SetFillColor(0);
leg3->SetTextSize(0.04);
leg1->AddEntry(gr1[ipt],Form("%.1f<p_{T, trig}<%.1f",ptbin[ipt],ptbin[ipt+1]),"lp");
leg2->AddEntry(gr2[ipt],Form("%.1f<p_{T, trig}<%.1f",ptbin[ipt],ptbin[ipt+1]),"lp");
leg3->AddEntry(gr3[ipt],Form("%.1f<p_{T, trig}<%.1f",ptbin[ipt],ptbin[ipt+1]),"lp");
c4->cd();
if(ipt==0){gr1[ipt]->Draw("AP");leg1->Draw("same");t.DrawLatex(0.2,0.2,"c_{1}  arm");}
else gr1[ipt]->Draw("Psame");

c5->cd();
if(ipt==0){gr2[ipt]->Draw("AP");leg2->Draw("same");t.DrawLatex(0.2,0.2,"c_{2}  arm");}
else gr2[ipt]->Draw("Psame");
	
c6->cd();
if(ipt==0){gr3[ipt]->Draw("AP");leg3->Draw("same");t.DrawLatex(0.2,0.2,"c_{3}  arm");}
else gr3[ipt]->Draw("Psame");

}

c4->Print(Form("fig/c1_cent_%s.png",type.Data()));
c5->Print(Form("fig/c2_cent_%s.png",type.Data()));
c6->Print(Form("fig/c3_cent_%s.png",type.Data()));

TCanvas *c4 = new TCanvas();
TCanvas *c5 = new TCanvas();
TCanvas *c6 = new TCanvas();
TLegend *leg1 = new TLegend(0.22,0.82-ncent*0.03,0.42,0.88);
TLegend *leg2 = new TLegend(0.22,0.82-ncent*0.03,0.42,0.88);
TLegend *leg3 = new TLegend(0.22,0.82-ncent*0.03,0.42,0.88);
for(int icent = 0;icent<ncent; icent++){
gr1[icent] = new TGraphErrors(npt,ptmean,c1[icent],0,c1err[icent]);
gr2[icent] = new TGraphErrors(npt,ptmean,c2[icent],0,c2err[icent]);
gr3[icent] = new TGraphErrors(npt,ptmean,c3[icent],0,c3err[icent]);
SetStyle(*gr1[icent],psize,color[icent],style[icent],0,0);
SetStyle(*gr2[icent],psize,color[icent],style[icent],0,0);
SetStyle(*gr3[icent],psize,color[icent],style[icent],0,0);
SetRange(*gr1[icent],ptbin[0],-0.07,ptbin[npt],0.001);
SetRange(*gr2[icent],ptbin[0],-0.001,ptbin[npt],0.01);
SetRange(*gr3[icent],ptbin[0],-0.015,ptbin[npt],0.015);
gr1[icent]->SetTitle("");
gr1[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr1[icent]->GetYaxis()->SetTitle("c_{1}");
gr2[icent]->SetTitle("");
gr2[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr2[icent]->GetYaxis()->SetTitle("c_{2}");
gr3[icent]->SetTitle("");
gr3[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr3[icent]->GetYaxis()->SetTitle("c_{3}");
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetFillColor(0);
leg2->SetTextSize(0.04);
leg3->SetFillColor(0);
leg3->SetTextSize(0.04);
leg1->AddEntry(gr1[icent],Form("%.2f\%<centrality<%.1f\%",centbin[icent],centbin[icent+1]),"lp");
leg2->AddEntry(gr2[icent],Form("%.2f\%<centrality<%.1f\%",centbin[icent],centbin[icent+1]),"lp");
leg3->AddEntry(gr3[icent],Form("%.2f\%<centrality<%.1f\%",centbin[icent],centbin[icent+1]),"lp");
c4->cd();
if(icent==0){gr1[icent]->Draw("AP");leg1->Draw("same");t.DrawLatex(0.2,0.2,"c_{1}  arm");}
else gr1[icent]->Draw("Psame");

c5->cd();
if(icent==0){gr2[icent]->Draw("AP");leg2->Draw("same");t.DrawLatex(0.2,0.2,"c_{2}  arm");}
else gr2[icent]->Draw("Psame");
	
c6->cd();
if(icent==0){gr3[icent]->Draw("AP");leg3->Draw("same");t.DrawLatex(0.2,0.2,"c_{3}  arm");}
else gr3[icent]->Draw("Psame");
}

c4->Print(Form("fig/c1_pt_%s.png",type.Data()));
c5->Print(Form("fig/c2_pt_%s.png",type.Data()));
c6->Print(Form("fig/c3_pt_%s.png",type.Data()));
	
}
