#include "root_setting.h"
void plotcorr_4p_pt_compMC(){
//TString dire = "north";
TString dire = "south";
TString coll = "pAu";
//TString corr = "cntbbc150M";
TString corr = "cntbbcmbtrig2B";
//TString corr = "cntbbc";
c1 = new TCanvas("c1"," ",800,460);
makeMultiPanelCanvas(c1,3,2,0,0,0.25,0.2,0.03);
c2 = new TCanvas("c2"," ",800,460);
makeMultiPanelCanvas(c2,3,2,0,0,0.25,0.2,0.03);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
  float const PI = acos(-1.0);
  const int npt = 10;
  const int ncent = 10;
  const int ncentd = 8;
 // TFile *f=TFile::Open(Form("../../MCgen/PHHijing/Nonflow/mergedAnav2%s%s.root",coll.Data(),corr.Data()));
  TFile *f=TFile::Open(Form("../../MCgen/PHHijing/Nonflow/mergedAna%s%s.root",coll.Data(),corr.Data()));
  TFile *fdata=TFile::Open(Form("merged_AnapAumbcentral_PC3Sigma2.root"));
  //TFile *f=TFile::Open(Form("mergedAna%s%s.root",coll.Data(),corr.Data()));
  ofstream fout;
  if(!corr.Contains("bbc"))
  fout.open(Form("Yjet%s_pt.dat",coll.Data()));
TH1F* kforebbcw[npt][ncent];
TH1F* hforebbcw[npt][ncent];
TH1F* kbackbbcw2[npt][ncent];
TH1F* kforebbcw_centIn[npt];
TH1F* hforebbcw_centIn[npt];
TH1F* kbackbbcw2_centIn[npt];
TH1F* kforebbcwd[npt][ncent];
TH1F* hforebbcwd[npt][ncent];
TH1F* kbackbbcw2d[npt][ncent];
TH1F* kforebbcw_centInd[npt];
TH1F* hforebbcw_centInd[npt];
TH1F* kbackbbcw2_centInd[npt];
double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
//double centbin[ncent+1] = {0,1,5,10,20,30,40,60,84}; //fvtx or cent
double centbin[ncent+1] = {100.00,97.51,92.65,85.92,78.00,69.62,53.55,33.84,15.74,5.89,0.00};
//double centbin[ncent+1] = {100.00,0.00};
//double bbccentbin[ncent+1] = {200.0,39.3,27.0,21.4,15.6,9.5,5.8,2.9,0.0}; //fvtx mb cent
//double centbin[ncent] = {108,83.3333,75,56.6667,48.6667,40,29.6667,21.6667,0}; //fvtx or
// double centbin[ncent] = {110,85.3333,77.3333,59.3333,51.6667,43,32.6667,25.3333,0};   //fvtx south
double bbccentbin[ncent+1] = {0,1,2,3,4,5,7,10,14,18,40};
//double bbccentbin[ncent+1] = {0,40};
//double bbccentbin[ncent+1] = {0,1,2,3,5,7,9,12,16,21,100}; //dAu
//double bbccentbin[ncent+1] = {0,2,4,6,8,11,16,20,26,30,100}; //HeAu
 //        
//  double centbin[ncent] = {152.333,118.333,102.333,69.6667,57.3333,45,31,21,0};   //fvtx and
/*double centmin = 1.0; double centmax = 3.0;//has to be boundary
for(int ipt=0; ipt<ncent; ipt++){
if(centmin >= centbin[ipt] && centmin < centbin[ipt+1]){int xcentmin = ipt; continue;}
if(centmax >= centbin[ipt] && centmax < centbin[ipt+1]){int xcentmax = ipt; continue;}
}*/
for(int icent_a=0;icent_a<1;icent_a++){
    if(coll.Contains("Au")){
int xcentmin = icent_a*1+9;
int xcentmind = icent_a*1;
int xcentmax = icent_a*1+10;
int xcentmaxd = icent_a*1+2;
    }
    else{
int xcentmin = icent_a*1;
int xcentmax = (icent_a+1)*10;
int xcentmaxd = (icent_a+1)*6;
    }
double centmin = centbin[xcentmin];
double centmax = centbin[xcentmax];

for(int ipt=0; ipt<npt; ipt++){
for(int icent=0; icent<ncent; icent++){
kforebbcw[ipt][icent] = (TH1F*)f->Get(Form("kfore%sbbc_%d_%d",dire.Data(),icent,ipt));
hforebbcw[ipt][icent] = (TH1F*)f->Get(Form("hfore%sbbc_%d_%d",dire.Data(),icent,ipt));
kbackbbcw2[ipt][icent] = (TH1F*)f->Get(Form("kback%sbbc2_%d_%d",dire.Data(),icent,ipt));
    if(!coll.Contains("Au")){
htemp = (TH1F*)f->Get(Form("kforenorthbbc_%d_%d",icent,ipt));
kforebbcw[ipt][icent]->Add(htemp);
htemp = (TH1F*)f->Get(Form("hforenorthbbc_%d_%d",icent,ipt));
hforebbcw[ipt][icent]->Add(htemp);
htemp = (TH1F*)f->Get(Form("kbacknorthbbc2_%d_%d",icent,ipt));
kbackbbcw2[ipt][icent]->Add(htemp);
    }
}
for(int icent=0; icent<ncentd; icent++){
kforebbcwd[ipt][icent] = (TH1F*)fdata->Get(Form("kfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
hforebbcwd[ipt][icent] = (TH1F*)fdata->Get(Form("hfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
kbackbbcw2d[ipt][icent] = (TH1F*)fdata->Get(Form("kback%sbbcw2_%d_%d",dire.Data(),icent,ipt));
    if(!coll.Contains("Au")){
htemp = (TH1F*)fdata->Get(Form("kforenorthbbcw_%d_%d",icent,ipt));
kforebbcwd[ipt][icent]->Add(htemp);
htemp = (TH1F*)fdata->Get(Form("hforenorthbbcw_%d_%d",icent,ipt));
hforebbcwd[ipt][icent]->Add(htemp);
htemp = (TH1F*)fdata->Get(Form("kbacknorthbbcw2_%d_%d",icent,ipt));
kbackbbcw2d[ipt][icent]->Add(htemp);
    }
}
kforebbcw_centIn[ipt] = (TH1F*)kforebbcw[ipt][xcentmin]->Clone();
hforebbcw_centIn[ipt] = (TH1F*)hforebbcw[ipt][xcentmin]->Clone();
kbackbbcw2_centIn[ipt] = (TH1F*)kbackbbcw2[ipt][xcentmin]->Clone();
kforebbcw_centInd[ipt] = (TH1F*)kforebbcw[ipt][xcentmin]->Clone();
hforebbcw_centInd[ipt] = (TH1F*)hforebbcw[ipt][xcentmin]->Clone();
kbackbbcw2_centInd[ipt] = (TH1F*)kbackbbcw2[ipt][xcentmin]->Clone();
for(int icent=xcentmin+1; icent<xcentmax; icent++){
kforebbcw_centIn[ipt]->Add(kforebbcw[ipt][icent]);
hforebbcw_centIn[ipt]->Add(hforebbcw[ipt][icent]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2[ipt][icent]);
}
for(int icent=xcentmind+1; icent<xcentmaxd; icent++){
kforebbcw_centInd[ipt]->Add(kforebbcwd[ipt][icent]);
hforebbcw_centInd[ipt]->Add(hforebbcwd[ipt][icent]);
kbackbbcw2_centInd[ipt]->Add(kbackbbcw2d[ipt][icent]);
}
kforebbcw_centIn[ipt]->Rebin(2);
hforebbcw_centIn[ipt]->Rebin(2);
kbackbbcw2_centIn[ipt]->Rebin(2);
kforebbcw_centInd[ipt]->Rebin(2);
hforebbcw_centInd[ipt]->Rebin(2);
kbackbbcw2_centInd[ipt]->Rebin(2);
}

TH1F* hpp[npt];
TH1F* hbackpp[npt];
TH1F* hppd[npt];
TH1F* hppratio[npt];
TH1F* hbackppd[npt];
  
int iptstart = 0;
for(int ipt=iptstart; ipt<6; ipt++){
	
if(ipt==iptstart){
for(int ipt_add=0;ipt_add<iptstart;ipt_add++){
kforebbcw_centIn[ipt]->Add(kforebbcw_centIn[ipt_add]);
hforebbcw_centIn[ipt]->Add(hforebbcw_centIn[ipt_add]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2_centIn[ipt_add]);
kforebbcw_centInd[ipt]->Add(kforebbcw_centInd[ipt_add]);
hforebbcw_centInd[ipt]->Add(hforebbcw_centInd[ipt_add]);
kbackbbcw2_centInd[ipt]->Add(kbackbbcw2_centInd[ipt_add]);
}
ptbin[ipt] = ptbin[0];
}

/*
	if(ipt==6){
for(int ipt_add=6;ipt_add<ncent;ipt_add++){
kforebbcw_centIn[ipt]->Add(kforebbcw_centIn[ipt_add]);
hforebbcw_centIn[ipt]->Add(hforebbcw_centIn[ipt_add]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2_centIn[ipt_add]);
}
centbin[ipt] = 800;
}
*/
  hpp[ipt] = (TH1F*)kforebbcw_centIn[ipt]->Clone();
  hbackpp[ipt] = (TH1F*)kbackbbcw2_centIn[ipt]->Clone();
  hppd[ipt] = (TH1F*)kforebbcw_centInd[ipt]->Clone();
  hppratio[ipt] = (TH1F*)kforebbcw_centInd[ipt]->Clone();
  hbackppd[ipt] = (TH1F*)kbackbbcw2_centInd[ipt]->Clone();
  float nbackpp = hbackpp[ipt]->Integral()/2.0/PI;
  float nforepp = hpp[ipt]->Integral()/2.0/PI;
  float nbackppd = hbackppd[ipt]->Integral()/2.0/PI;
  float nforeppd = hppd[ipt]->Integral()/2.0/PI;
  //float ntrig0 = centforedis_0->Integral(11,30);
   for(int i=0; i<20; i++){
     float pp_cont = 1.0*hpp[ipt]->GetBinContent(i+1);
     float pp_contd = 1.0*hppd[ipt]->GetBinContent(i+1);
     //float pp0_err = 1.0*hpp[ipt]->GetBinError(i+1);
     float weight2 = sqrt(1.0*hforebbcw_centIn[ipt]->GetBinContent(i+1));
     float weight2d = sqrt(1.0*hforebbcw_centInd[ipt]->GetBinContent(i+1));

     float backpp_cont = 1.0*hbackpp[ipt]->GetBinContent(i+1);
     float backpp_contd = 1.0*hbackppd[ipt]->GetBinContent(i+1);
    
     float con = pp_cont/backpp_cont*nbackpp/nforepp;
     float err = weight2/backpp_cont*nbackpp/nforepp;
     float cond = pp_contd/backpp_contd*nbackppd/nforeppd;
     float errd = weight2d/backpp_contd*nbackppd/nforeppd;

     hpp[ipt]->SetBinContent(i+1, con);
     hpp[ipt]->SetBinError(i+1, err);
     hppd[ipt]->SetBinContent(i+1, cond);
     hppd[ipt]->SetBinError(i+1, errd);
     hppratio[ipt]->SetBinContent(i+1,con/cond);
     hppratio[ipt]->SetBinError(i+1,con/cond*sqrt((err/con)**2+(errd/cond)**2));
   }

  c1->cd(ipt-iptstart+1);
  gPad->SetTicks(1,1);
if(corr.Contains("bbc")){
    if(coll.Contains("Au")){
  float ymax = 1.1-0.0;//hpp[ipt]->GetMaximum()*1.1;
  float ymin = 0.9+0.0;//hpp[ipt]->GetMinimum()*0.9;
    }
    else{
  float ymax = 1.39-0.0;//hpp[ipt]->GetMaximum()*1.1;
  float ymin = 0.61+0.0;//hpp[ipt]->GetMinimum()*0.9;
    }
}
else{
  float ymax = 7.0-0.0;//hpp[ipt]->GetMaximum()*1.1;
  float ymin = 0.0+0.0;//hpp[ipt]->GetMinimum()*0.9;
}
  hpp[ipt]->SetMinimum(ymin);
  hpp[ipt]->SetMaximum(ymax);

  hpp[ipt]->SetMarkerStyle(20);
  hpp[ipt]->SetMarkerSize(1.3);
  hppd[ipt]->SetMarkerStyle(24);
  hppd[ipt]->SetMarkerSize(1.3);
  
  hpp[ipt]->SetMarkerColor(4);
  hppd[ipt]->SetMarkerColor(2);
  hpp[ipt]->GetYaxis()->SetLabelSize(0.08);
  hpp[ipt]->GetXaxis()->SetLabelSize(0.08);
  hpp[ipt]->GetYaxis()->SetNdivisions(505);
  //hpp[ipt]->GetYaxis()->SetTitleSize(0.6);
  //hpp[ipt]->GetXaxis()->SetTitleSize(0.0);
  hpp[ipt]->Draw("PE");
  hppd[ipt]->Draw("PEsame");
  c2->cd(ipt-iptstart+1);
  hppratio[ipt]->SetMinimum(0.9);
  hppratio[ipt]->SetMaximum(1.1);
  hppratio[ipt]->SetMarkerStyle(20);
  hppratio[ipt]->SetMarkerSize(1.3);
  hppratio[ipt]->Draw("PE");
TLegend *leg1 = new TLegend(0.44442,0.22,0.82,0.32);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.08);
  leg1->SetBorderSize(0);
  leg1->AddEntry(hppratio[ipt],"MC/data","P");
  if(ipt-iptstart==0) leg1->Draw();

  c1->cd(ipt-iptstart+1);
  if(corr.Contains("bbc")){
  TF1 *fun0 = new TF1("fun0","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)+2*[4]*cos(4*x))", -0.5*PI, 1.5*PI);

  fun0->SetLineColor(1);
 // hpp[ipt]->Fit("fun0","R");


  TF1 *fun1 = new TF1("fun1","[0]*(1+2*[1]*cos(x))",   -0.5*PI, 1.5*PI);
  TF1 *fun2 = new TF1("fun2","[0]*(1+2*[1]*cos(2*x))", -0.5*PI, 1.5*PI);
  TF1 *fun3 = new TF1("fun3","[0]*(1+2*[1]*cos(3*x))", -0.5*PI, 1.5*PI);
  TF1 *fun4 = new TF1("fun4","[0]*(1+2*[1]*cos(4*x))", -0.5*PI, 1.5*PI);

//  TF1 *fun0 = new TF1("fun0","pol2",0.1,2);
//  hpp[ipt]->Fit("fun0","R");
TLegend *leg1 = new TLegend(0.22,0.22,0.82,0.32);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.08);
  leg1->SetBorderSize(0);
  leg1->AddEntry(fun0,"1+#Sigma2c_{n}cos(n#Delta#phi)","L");
//  leg1->AddEntry(fun0,"pol2","L");
  //if(ipt-iptstart==1)leg1->Draw();

  TLegend *leg2 = new TLegend(0.12,0.60,0.32,0.82);
  leg2->SetFillColor(10);
  leg2->SetLineStyle(4000);
  leg2->SetLineColor(10);
  leg2->SetLineWidth(0.);
  leg2->SetTextSize(0.080);
  leg2->SetBorderSize(0);
  leg2->AddEntry(fun1,"1+2c_{1}cos(#Delta#phi)","L");
  leg2->AddEntry(fun2,"1+2c_{2}cos(2#Delta#phi)","L");
  leg2->AddEntry(fun3,"1+2c_{3}cos(3#Delta#phi)","L");
  leg2->AddEntry(fun4,"1+2c_{4}cos(4#Delta#phi)","L");
  //if(ipt-iptstart==1)leg2->Draw();
  

  fun1->SetParameters(fun0->GetParameter(0), fun0->GetParameter(1));
  fun2->SetParameters(fun0->GetParameter(0), fun0->GetParameter(2));
  fun3->SetParameters(fun0->GetParameter(0), fun0->GetParameter(3));
  fun4->SetParameters(fun0->GetParameter(0), fun0->GetParameter(4));

  fun1->SetLineColor(3);
  fun2->SetLineColor(2);
  fun3->SetLineColor(6);
  fun4->SetLineColor(9);

  fun1->SetLineStyle(2);
  fun2->SetLineStyle(3);
  fun3->SetLineStyle(4);
  fun4->SetLineStyle(5);
/*
  fun1->Draw("same");
  fun2->Draw("same");
  fun3->Draw("same");
  fun4->Draw("same");
  */
  }
else{
  TLegend *leg3 = new TLegend(0.12,0.60,0.32,0.82);
  leg3->SetFillColor(10);
  leg3->SetLineStyle(4000);
  leg3->SetLineColor(10);
  leg3->SetLineWidth(0.);
  leg3->SetTextSize(0.080);
  leg3->SetBorderSize(0);
  leg3->AddEntry(hpp[ipt],"before zyam", "p");
  //TF1 *fun0 = new TF1("fun0","[0]*x*x+[1]*x+[2]",0.1,2);
  TF1 *fun0 = new TF1("fun0","gaus(0)+gaus(3)",-1.2,PI/2);
  fun0->SetParameters(3,0,1,0.2,0,4);
  hpp[ipt]->Fit("fun0","R");
  hpp[ipt]->DrawCopy();
  TF1 *fun1 = new TF1("fun1","gaus(0)",-1.2,PI/2);
  fun1->SetParameters(fun0->GetParameter(0),fun0->GetParameter(1),fun0->GetParameter(2));
/*  double a=fun0->GetParameter(0);
  double b=fun0->GetParameter(1);
  double c=fun0->GetParameter(2);
  if(-b/2./a<2.)
    double zyam = fun0->Eval(-b/2./a);
  else
    double zyam = fun0->Eval(2.);
    */
  double zyam = fun0->Eval(PI/2);
  hpp_cp = (TH1D*)hpp[ipt]->Clone(Form("hpp_cp_%d",ipt));
   for(int i=0; i<21; i++){
  hpp_cp->SetBinContent(i,hpp[ipt]->GetBinContent(i)-zyam);
   }
   hpp_cp->SetMarkerColor(2);
   hpp_cp->Draw("same");
   //fout<<hpp[ipt]->Integral(hpp[ipt]->FindBin(-1.2),hpp[ipt]->FindBin(PI/2))<<"\t";
   fout<<0.95*fabs(fun1->GetParameter(0)*fun1->GetParameter(2)*sqrt(2*PI))<<endl;
  leg3->AddEntry(hpp_cp,"after zyam", "p");
  if(ipt-iptstart==1)leg3->Draw();
}
  hpp[ipt]->GetXaxis()->SetTitle("#Delta#phi");
  hpp[ipt]->GetXaxis()->SetTitleSize(0.08);
  hpp[ipt]->GetXaxis()->SetLabelSize(0.08);
  hpp[ipt]->GetXaxis()->CenterTitle();
  hpp[ipt]->GetXaxis()->SetTitleOffset(1.1);
  hpp[ipt]->GetYaxis()->SetTitle("C(d#phi)");
  hpp[ipt]->GetYaxis()->SetTitleSize(0.08);
  hpp[ipt]->GetYaxis()->CenterTitle();
  hpp[ipt]->GetYaxis()->SetTitleOffset(1.3);
  hpp[ipt]->GetYaxis()->SetLabelSize(0.08);
  hpp[ipt]->GetYaxis()->SetTitle("C(#Delta#phi)");
  hppratio[ipt]->GetXaxis()->SetTitle("#Delta#phi");
  hppratio[ipt]->GetXaxis()->SetTitleSize(0.08);
  hppratio[ipt]->GetXaxis()->SetLabelSize(0.08);
  hppratio[ipt]->GetXaxis()->CenterTitle();
  hppratio[ipt]->GetXaxis()->SetTitleOffset(1.1);
  hppratio[ipt]->GetYaxis()->SetTitle("C(d#phi)");
  hppratio[ipt]->GetYaxis()->SetTitleSize(0.08);
  hppratio[ipt]->GetYaxis()->CenterTitle();
  hppratio[ipt]->GetYaxis()->SetTitleOffset(1.3);
  hppratio[ipt]->GetYaxis()->SetLabelSize(0.08);
  hppratio[ipt]->GetYaxis()->SetTitle("C(#Delta#phi)");
  TLegend *leg3 = new TLegend(0.12,0.60,0.32,0.82);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->SetTextSize(0.06);
  leg3->AddEntry(hpp[ipt],Form("%s HIJING 200GeV ",coll.Data()));
  leg3->AddEntry(hppd[ipt],Form("%s data 200GeV ",coll.Data()));
  if(ipt-iptstart==1)
  leg3->Draw();
 // TLatex *t=new TLatex(-1.0,0.88*(ymax-ymin)+ymin, Form("%d)",ipt+1));
 // t->SetTextSize(0.08);
 // t->Draw();

  TMarker *m0 = new TMarker(-0.3, 0.91*(ymax-ymin)+ymin, 20);

  if(ipt-iptstart==0){
  TMarker *m0 = new TMarker(-1.1, 0.91*(ymax-ymin)+ymin, 20);
  TLatex *t=new TLatex(-1.1,0.88*(ymax-ymin)+ymin, Form("%s HIJING 200GeV ",coll.Data()));
  t->SetTextSize(0.08);
 // t->Draw();
  }
 // else{
 // TLatex *t=new TLatex(-0.9,0.88*(ymax-ymin)+ymin, Form("%.f <Nch<= %.f",bbccentbin[ipt],bbccentbin[ipt+1]));
 // }
  if(ipt-iptstart==5){
    TLatex *t=new TLatex(-1.1,5.24*(ymax-ymin)+ymin, Form("with weight = particle Energy"));
    t->SetTextSize(0.07);
  //  t->Draw();
  }
  m0->SetMarkerSize(1.3);
  m0->SetMarkerColor(4);
  //m0->Draw();
  
  TLatex *t=new TLatex(-0.4,0.08*(ymax-ymin)+ymin, Form("%.1f < pt < %.1f GeV/c",ptbin[ipt],ptbin[ipt+1]));
  t->SetTextSize(0.08);
  t->Draw();

  //TLatex *t=new TLatex(-1.2,0.78*(ymax-ymin)+ymin, Form("%.1f< cent <%.1f%%",centmax,centmin));
  TLatex *t=new TLatex(-1.2,0.78*(ymax-ymin)+ymin, Form("%.1f< cent <%.f%%",centmax,5.));
  t->SetTextSize(0.08);
  if(ipt-iptstart==2)t->Draw();

  TLatex *t=new TLatex(-1.0,0.66*(ymax-ymin)+ymin, "|#eta_{trig}|<0.35");
  t->SetTextSize(0.08);
  if(ipt-iptstart==0)t->Draw();

if(dire=="north")
  TLatex *t=new TLatex(-1.2,0.24*(ymax-ymin)+ymin, "3.0<#eta_{asso}<4.0");
else
if(corr.Contains("bbc"))
  TLatex *t=new TLatex(-1.2,0.24*(ymax-ymin)+ymin, "-4.0<#eta_{asso}<-3.0");
  else
  TLatex *t=new TLatex(-1.2,0.24*(ymax-ymin)+ymin, "|#eta_{asso}|<0.35");
  t->SetTextSize(0.08);
  t->SetTextColor(2);
  if(ipt-iptstart==0)t->Draw();

}
  c1->Print(Form("fig/%s%sridge_%s_%.1f_%.1f_dataMC.gif",coll.Data(),corr.Data(),dire.Data(),centmin,centmax));
  c1->Print(Form("fig/%s%sridge_%s_%.1f_%.1f_dataMC.pdf",coll.Data(),corr.Data(),dire.Data(),centmin,centmax));
  c2->Print(Form("fig/%s%sridge_%s_%.1f_%.1f_dataMCratio.png",coll.Data(),corr.Data(),dire.Data(),centmin,centmax));
  c2->Print(Form("fig/%s%sridge_%s_%.1f_%.1f_dataMCratio.pdf",coll.Data(),corr.Data(),dire.Data(),centmin,centmax));
}

/*  
  cout<<"*************** v2 ***********"<<endl;
  float v2_0 = funvn0->GetParameter(1)/(zym_pp0 + funvn0->GetParameter(0));
  float v2_1 = funvn1->GetParameter(1)/(zym_pp1 + funvn1->GetParameter(0));
  float v2_2 = funvn2->GetParameter(1)/(zym_pp2 + funvn2->GetParameter(0));
  float v2_3 = funvn3->GetParameter(1)/(zym_pp3 + funvn3->GetParameter(0));
  

  cout<<v2_0<<" "<<v2_1<<" "<<v2_2<<" "<<v2_3<<endl;
 
  
  cout<<funvn0->GetParameter(0)*funvn0->GetParameter(1)<<" "
      <<funvn1->GetParameter(0)*funvn1->GetParameter(1)<<" "
      <<funvn2->GetParameter(0)*funvn2->GetParameter(1)<<" "
      <<funvn3->GetParameter(0)*funvn3->GetParameter(1)<<endl;
  

  cout<<"*************** v2 ***********"<<endl;
  cout<<funvn0->GetParameter(2)<<" "<<funvn1->GetParameter(2)<<" "<<funvn2->GetParameter(2)<<" "<<funvn3->GetParameter(2)<<endl;

  cout<<"*************** v3 ***********"<<endl;
  cout<<funvn0->GetParameter(3)<<" "<<funvn1->GetParameter(3)<<" "<<funvn2->GetParameter(3)<<" "<<funvn3->GetParameter(3)<<endl;
*/
}
