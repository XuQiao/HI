#include "fstream.h"
#include "/home/xuq7/HI/utilities/SimplifyLife.C"

void plotcnbbccsouth(string t_){
gStyle->SetLegendFillColor(0);
gStyle->SetLegendBorderSize(0);
TString type = t_.c_str();
double psize = 1.4;
const int color[10] = {1,2,4,6,4,8,6,2,1,7};
const int style[10] = {20,24,27,28,29,33,34,31,22,21};
if(type.Contains("pt") && type != "ptccentc"){
const int ncent = 9;
double centbin[ncent+1] = {0,0.1,0.5,1,5,10,20,40,60,100}; 
double centmidp[ncent] = {110,50.3333,41.8333,34.3,24,19.6667,12.66667,7.66667,2.8}; //mb
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
const double centbin[ncent+1] = {0,100};
const double centmidp[ncent] = {50};
const int npt = 10;
const double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
const double ptmean[npt] = {0.360943, 0.691833, 1.1911, 1.69654, 2.20117, 2.70571, 3.2097, 3.71372, 4.21814, 4.72014};
}
else if(type=="ptccentc"){
const int ncent = 5;
const double centbin[ncent+1] = {0,1,5,10,40,100};
const double centmidp[ncent] = {0,.5,3,7.5,25,70};
const int npt = 5;
const double ptbin[npt+1] = {0.2,0.5,1.0,2.0,3.0,4.0};
const double ptmean[npt] = {0.360943, 0.691833, 1.29345, 2.32523, 3.3542};
}
double c1north[ncent][npt];
double c1south[ncent][npt];
double c1northerr[ncent][npt];
double c1southerr[ncent][npt];
double c2north[ncent][npt];
double c2south[ncent][npt];
double c2northerr[ncent][npt];
double c2southerr[ncent][npt];
double c3north[ncent][npt];
double c3south[ncent][npt];
double c3northerr[ncent][npt];
double c3southerr[ncent][npt];
double c1north_[npt][ncent];
double c1south_[npt][ncent];
double c1northerr_[npt][ncent];
double c1southerr_[npt][ncent];
double c2north_[npt][ncent];
double c2south_[npt][ncent];
double c2northerr_[npt][ncent];
double c2southerr_[npt][ncent];
double c3north_[npt][ncent];
double c3south_[npt][ncent];
double c3northerr_[npt][ncent];
double c3southerr_[npt][ncent];
ifstream fnorth;
ifstream fsouth;
fnorth.open(Form("c1_c2_%s_north.dat",type.Data()));
fsouth.open(Form("c1_c2_%s_south.dat",type.Data()));
for(int icent=0;icent<ncent;icent++){
for(int ipt=0;ipt<npt;ipt++){
fnorth>>c1north[icent][ipt];
fnorth>>c1northerr[icent][ipt];
fnorth>>c2north[icent][ipt];
fnorth>>c2northerr[icent][ipt];
fnorth>>c3north[icent][ipt];
fnorth>>c3northerr[icent][ipt];
c1north_[ipt][icent] = c1north[icent][ipt];
c1northerr_[ipt][icent] = c1northerr[icent][ipt];
c2north_[ipt][icent] = c2north[icent][ipt];
c2northerr_[ipt][icent] = c2northerr[icent][ipt];
c3north_[ipt][icent] = c3north[icent][ipt];
c3northerr_[ipt][icent] = c3northerr[icent][ipt];
fsouth>>c1south[icent][ipt];
fsouth>>c1southerr[icent][ipt];
fsouth>>c2south[icent][ipt];
fsouth>>c2southerr[icent][ipt];
fsouth>>c3south[icent][ipt];
fsouth>>c3southerr[icent][ipt];
c1south_[ipt][icent] = c1south[icent][ipt];
c1southerr_[ipt][icent] = c1southerr[icent][ipt];
c2south_[ipt][icent] = c2south[icent][ipt];
c2southerr_[ipt][icent] = c2southerr[icent][ipt];
c3south_[ipt][icent] = c3south[icent][ipt];
c3southerr_[ipt][icent] = c3southerr[icent][ipt];
}
}

TGraphErrors *gr1north[npt>ncent?npt:ncent];
TGraphErrors *gr1south[npt>ncent?npt:ncent];
TGraphErrors *gr2north[npt>ncent?npt:ncent];
TGraphErrors *gr2south[npt>ncent?npt:ncent];
TGraphErrors *gr3north[npt>ncent?npt:ncent];
TGraphErrors *gr3south[npt>ncent?npt:ncent];
TCanvas *c1 = new TCanvas();
TCanvas *c2 = new TCanvas();
TCanvas *c3 = new TCanvas();
TCanvas *c4 = new TCanvas();
TCanvas *c5 = new TCanvas();
TCanvas *c6 = new TCanvas();
TLatex t;
t.SetNDC();
t.SetTextSize(0.04);
TLegend *leg1 = new TLegend(0.62,0.82-npt*0.03,0.82,0.88);
TLegend *leg2 = new TLegend(0.62,0.82-npt*0.03,0.82,0.88);
TLegend *leg3 = new TLegend(0.22,0.82-npt*0.03,0.42,0.88);
TLegend *leg4 = new TLegend(0.22,0.82-npt*0.03,0.42,0.88);
TLegend *leg5 = new TLegend(0.62,0.82-npt*0.03,0.82,0.88);
TLegend *leg6 = new TLegend(0.62,0.82-npt*0.03,0.82,0.88);
for(int ipt = 0;ipt<npt; ipt++){
gr1north[ipt] = new TGraphErrors(ncent,centmidp,c1north_[ipt],0,c1northerr_[ipt]);
gr1south[ipt] = new TGraphErrors(ncent,centmidp,c1south_[ipt],0,c1southerr_[ipt]);
gr2north[ipt] = new TGraphErrors(ncent,centmidp,c2north_[ipt],0,c2northerr_[ipt]);
gr2south[ipt] = new TGraphErrors(ncent,centmidp,c2south_[ipt],0,c2southerr_[ipt]);
gr3north[ipt] = new TGraphErrors(ncent,centmidp,c3north_[ipt],0,c3northerr_[ipt]);
gr3south[ipt] = new TGraphErrors(ncent,centmidp,c3south_[ipt],0,c3southerr_[ipt]);
SetStyle(*gr1north[ipt],psize,color[ipt],style[ipt],0,0);
SetStyle(*gr1south[ipt],psize,color[ipt],style[ipt],0,0);
SetStyle(*gr2north[ipt],psize,color[ipt],style[ipt],0,0);
SetStyle(*gr2south[ipt],psize,color[ipt],style[ipt],0,0);
SetStyle(*gr3north[ipt],psize,color[ipt],style[ipt],0,0);
SetStyle(*gr3south[ipt],psize,color[ipt],style[ipt],0,0);
SetRange(*gr1north[ipt],centbin[0],-0.10,centbin[ncent],0.001);
SetRange(*gr1south[ipt],centbin[0],-0.10,centbin[ncent],0.001);
SetRange(*gr2north[ipt],centbin[0],-0.001,centbin[ncent],0.02);
SetRange(*gr2south[ipt],centbin[0],-0.001,centbin[ncent],0.02);
SetRange(*gr3north[ipt],centbin[0],-0.015,centbin[ncent],0.015);
SetRange(*gr3south[ipt],centbin[0],-0.015,centbin[ncent],0.015);
gr1north[ipt]->SetTitle("");
gr1north[ipt]->GetXaxis()->SetTitle("vtxz position(cm)");
gr1north[ipt]->GetYaxis()->SetTitle("c_{1}");
gr1south[ipt]->SetTitle("");
gr1south[ipt]->GetXaxis()->SetTitle("vtxz position(cm)");
gr1south[ipt]->GetYaxis()->SetTitle("c_{1}");
gr2north[ipt]->SetTitle("");
gr2north[ipt]->GetXaxis()->SetTitle("vtxz position(cm)");
gr2north[ipt]->GetYaxis()->SetTitle("c_{2}");
gr2south[ipt]->SetTitle("");
gr2south[ipt]->GetXaxis()->SetTitle("vtxz position(cm)");
gr2south[ipt]->GetYaxis()->SetTitle("c_{2}");
gr3north[ipt]->SetTitle("");
gr3north[ipt]->GetXaxis()->SetTitle("vtxz position(cm)");
gr3north[ipt]->GetYaxis()->SetTitle("c_{3}");
gr3south[ipt]->SetTitle("");
gr3south[ipt]->GetXaxis()->SetTitle("vtxz position(cm)");
gr3south[ipt]->GetYaxis()->SetTitle("c_{3}");
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetFillColor(0);
leg2->SetTextSize(0.04);
leg3->SetFillColor(0);
leg3->SetTextSize(0.04);
leg4->SetFillColor(0);
leg4->SetTextSize(0.04);
leg5->SetFillColor(0);
leg5->SetTextSize(0.04);
leg6->SetFillColor(0);
leg6->SetTextSize(0.04);
leg1->AddEntry(gr1north[ipt],Form("%.1f<p_{T, trig}<%.1f",ptbin[ipt],ptbin[ipt+1]),"lp");
leg2->AddEntry(gr1south[ipt],Form("%.1f<p_{T, trig}<%.1f",ptbin[ipt],ptbin[ipt+1]),"lp");
leg3->AddEntry(gr2north[ipt],Form("%.1f<p_{T, trig}<%.1f",ptbin[ipt],ptbin[ipt+1]),"lp");
leg4->AddEntry(gr2south[ipt],Form("%.1f<p_{T, trig}<%.1f",ptbin[ipt],ptbin[ipt+1]),"lp");
leg5->AddEntry(gr3north[ipt],Form("%.1f<p_{T, trig}<%.1f",ptbin[ipt],ptbin[ipt+1]),"lp");
leg6->AddEntry(gr3south[ipt],Form("%.1f<p_{T, trig}<%.1f",ptbin[ipt],ptbin[ipt+1]),"lp");
c1->cd();
if(ipt==0){gr1north[ipt]->Draw("AP");leg1->Draw("same");t.DrawLatex(0.2,0.2,"c_{1} 1.5<|eta_{asso}|<3.0");}
else gr1north[ipt]->Draw("Psame");

c2->cd();
if(ipt==0) {gr1south[ipt]->Draw("AP");leg2->Draw("same");t.DrawLatex(0.2,0.2,"c_{1} 1.0<|eta_{asso}|<3.0");}
else gr1south[ipt]->Draw("Psame");

c3->cd();
if(ipt==0){gr2north[ipt]->Draw("AP");leg3->Draw("same");t.DrawLatex(0.2,0.2,"c_{2} 1.5<|eta_{asso}|<3.0");}
else gr2north[ipt]->Draw("Psame");
	
c4->cd();
if(ipt==0) {gr2south[ipt]->Draw("AP");leg4->Draw("same");t.DrawLatex(0.2,0.2,"c_{2} 1.0<|eta_{asso}|<3.0");}
else gr2south[ipt]->Draw("Psame");

c5->cd();
if(ipt==0){gr3north[ipt]->Draw("AP");leg5->Draw("same");t.DrawLatex(0.2,0.2,"c_{3} 1.5<|eta_{asso}|<3.0");}
else gr3north[ipt]->Draw("Psame");

c6->cd();
if(ipt==0) {gr3south[ipt]->Draw("AP");leg6->Draw("same");t.DrawLatex(0.2,0.2,"c_{3} 1.0<|eta_{asso}|<3.0");}
else gr3south[ipt]->Draw("Psame");
}

c1->Print(Form("fig/c1north_cent_%s.png",type.Data()));
c2->Print(Form("fig/c1south_cent_%s.png",type.Data()));
c3->Print(Form("fig/c2north_cent_%s.png",type.Data()));
c4->Print(Form("fig/c2south_cent_%s.png",type.Data()));
c5->Print(Form("fig/c3north_cent_%s.png",type.Data()));
c6->Print(Form("fig/c3south_cent_%s.png",type.Data()));

TCanvas *c1 = new TCanvas();
TCanvas *c2 = new TCanvas();
TCanvas *c3 = new TCanvas();
TCanvas *c4 = new TCanvas();
TCanvas *c5 = new TCanvas();
TCanvas *c6 = new TCanvas();
TLegend *leg1 = new TLegend(0.22,0.82-ncent*0.03,0.42,0.88);
TLegend *leg2 = new TLegend(0.22,0.82-ncent*0.03,0.42,0.88);
TLegend *leg3 = new TLegend(0.22,0.82-ncent*0.03,0.42,0.88);
TLegend *leg4 = new TLegend(0.22,0.82-ncent*0.03,0.42,0.88);
TLegend *leg5 = new TLegend(0.22,0.82-ncent*0.03,0.42,0.88);
TLegend *leg6 = new TLegend(0.22,0.82-ncent*0.03,0.42,0.88);
for(int icent = 0;icent<ncent; icent++){
gr1north[icent] = new TGraphErrors(npt,ptmean,c1north[icent],0,c1northerr[icent]);
gr1south[icent] = new TGraphErrors(npt,ptmean,c1south[icent],0,c1southerr[icent]);
gr2north[icent] = new TGraphErrors(npt,ptmean,c2north[icent],0,c2northerr[icent]);
gr2south[icent] = new TGraphErrors(npt,ptmean,c2south[icent],0,c2southerr[icent]);
gr3north[icent] = new TGraphErrors(npt,ptmean,c3north[icent],0,c3northerr[icent]);
gr3south[icent] = new TGraphErrors(npt,ptmean,c3south[icent],0,c3southerr[icent]);
SetStyle(*gr1north[icent],psize,color[icent],style[icent],0,0);
SetStyle(*gr1south[icent],psize,color[icent],style[icent],0,0);
SetStyle(*gr2north[icent],psize,color[icent],style[icent],0,0);
SetStyle(*gr2south[icent],psize,color[icent],style[icent],0,0);
SetStyle(*gr3north[icent],psize,color[icent],style[icent],0,0);
SetStyle(*gr3south[icent],psize,color[icent],style[icent],0,0);
SetRange(*gr1north[icent],ptbin[0],-0.07,ptbin[npt],0.001);
SetRange(*gr1south[icent],ptbin[0],-0.07,ptbin[npt],0.001);
SetRange(*gr2north[icent],ptbin[0],-0.001,ptbin[npt],0.02);
SetRange(*gr2south[icent],ptbin[0],-0.001,ptbin[npt],0.02);
SetRange(*gr3north[icent],ptbin[0],-0.015,ptbin[npt],0.015);
SetRange(*gr3south[icent],ptbin[0],-0.015,ptbin[npt],0.015);
gr1north[icent]->SetTitle("");
gr1north[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr1north[icent]->GetYaxis()->SetTitle("c_{1}");
gr1south[icent]->SetTitle("");
gr1south[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr1south[icent]->GetYaxis()->SetTitle("c_{1}");
gr2north[icent]->SetTitle("");
gr2north[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr2north[icent]->GetYaxis()->SetTitle("c_{2}");
gr2south[icent]->SetTitle("");
gr2south[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr2south[icent]->GetYaxis()->SetTitle("c_{2}");
gr3north[icent]->SetTitle("");
gr3north[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr3north[icent]->GetYaxis()->SetTitle("c_{3}");
gr3south[icent]->SetTitle("");
gr3south[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr3south[icent]->GetYaxis()->SetTitle("c_{3}");
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetFillColor(0);
leg2->SetTextSize(0.04);
leg3->SetFillColor(0);
leg3->SetTextSize(0.04);
leg4->SetFillColor(0);
leg4->SetTextSize(0.04);
leg5->SetFillColor(0);
leg5->SetTextSize(0.04);
leg6->SetFillColor(0);
leg6->SetTextSize(0.04);
leg1->AddEntry(gr1north[icent],Form("%.f<=vtx z<%.f",centbin[icent],centbin[icent+1]),"lp");
leg2->AddEntry(gr1south[icent],Form("%.f<=vtx z<%.f",centbin[icent],centbin[icent+1]),"lp");
leg3->AddEntry(gr2north[icent],Form("%.f<=vtx z<%.f",centbin[icent],centbin[icent+1]),"lp");
leg4->AddEntry(gr2south[icent],Form("%.f<=vtx z<%.f",centbin[icent],centbin[icent+1]),"lp");
leg5->AddEntry(gr3north[icent],Form("%.f<=vtx z<%.f",centbin[icent],centbin[icent+1]),"lp");
leg6->AddEntry(gr3south[icent],Form("%.f<=vtx z<%.f",centbin[icent],centbin[icent+1]),"lp");
c1->cd();
if(icent==0){gr1north[icent]->Draw("AP");leg1->Draw("same");t.DrawLatex(0.2,0.2,"c_{1} 1.5<|eta_{asso}|<3.0");}
else gr1north[icent]->Draw("Psame");

c2->cd();
if(icent==0) {gr1south[icent]->Draw("AP");leg2->Draw("same");t.DrawLatex(0.2,0.2,"c_{1} 1.0<|eta_{asso}|<3.0");}
else gr1south[icent]->Draw("Psame");

c3->cd();
if(icent==0){gr2north[icent]->Draw("AP");leg3->Draw("same");t.DrawLatex(0.2,0.2,"c_{2} 1.5<|eta_{asso}|<3.0");}
else gr2north[icent]->Draw("Psame");
	
c4->cd();
if(icent==0) {gr2south[icent]->Draw("AP");leg4->Draw("same");t.DrawLatex(0.2,0.2,"c_{2} 1.0<|eta_{asso}|<3.0");}
else gr2south[icent]->Draw("Psame");

c5->cd();
if(icent==0){gr3north[icent]->Draw("AP");leg5->Draw("same");t.DrawLatex(0.2,0.2,"c_{3} 1.5<|eta_{asso}|<3.0");}
else gr3north[icent]->Draw("Psame");

c6->cd();
if(icent==0) {gr3south[icent]->Draw("AP");leg6->Draw("same");t.DrawLatex(0.2,0.2,"c_{3} 1.0<|eta_{asso}|<3.0");}
else gr3south[icent]->Draw("Psame");
}

c1->Print(Form("fig/c1north_pt_%s.png",type.Data()));
c2->Print(Form("fig/c1south_pt_%s.png",type.Data()));
c3->Print(Form("fig/c2north_pt_%s.png",type.Data()));
c4->Print(Form("fig/c2south_pt_%s.png",type.Data()));
c5->Print(Form("fig/c3north_pt_%s.png",type.Data()));
c6->Print(Form("fig/c3south_pt_%s.png",type.Data()));
	
}
