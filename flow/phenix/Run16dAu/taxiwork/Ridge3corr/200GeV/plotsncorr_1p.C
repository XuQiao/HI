#include "root_setting.h"
void plotsncorr_1p(){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
c1 = new TCanvas("c1"," ",500,500);
//makeMultiPanelCanvas(c1,4,1,0,0,0.25,0.2,0.03);
//TString dire = "north";
TString dire = "sn";
  float const PI = acos(-1.0);
const int ncent = 6;
const int npt = 1;
TFile *f=TFile::Open("../../../work/200GeV/output_Ridge.root");
TH1F* kforebbcw[npt][ncent];
TH1F* hforebbcw[npt][ncent];
TH1F* kbackbbcw2[npt][ncent];
TH1F* kforebbcw_centIn[npt];
TH1F* hforebbcw_centIn[npt];
TH1F* kbackbbcw2_centIn[npt];
double ptbin[npt+1] = {0,0};//{0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
//double centbin[ncent+1] = {0,1,5,10,20,30,40,60,100}; //fvtx or cent
//double centbin[ncent+1] = {0,1,5,10,20,30,40,60,100}; //fvtx or cent
double centbin[ncent+1] = {0,5,10,20,40,60,100}; //fvtx or cent
//double centbin[ncent+1] = {0,100}; //fvtx or cent
//double centbin[ncent] = {224,170,146,120,89,65,0};
/*double centmin = 1.0; double centmax = 3.0;//has to be boundary
for(int icent=0; icent<ncent; icent++){
if(centmin >= centbin[icent] && centmin < centbin[icent+1]){int xcentmin = icent; continue;}
if(centmax >= centbin[icent] && centmax < centbin[icent+1]){int xcentmax = icent; continue;}
}*/
for(int icent_a=0;icent_a<ncent;icent_a++){
int xcentmin = icent_a*1;
int xcentmax = (icent_a+1)*1;
double centmin = centbin[xcentmin];
double centmax = centbin[xcentmax];

for(int ipt=0; ipt<npt; ipt++){
for(int icent=0; icent<ncent; icent++){
kforebbcw[ipt][icent] = (TH1F*)f->Get(Form("kfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
//htemp = (TH1F*)f->Get(Form("kforeW%sbbcw_%d_%d",dire.Data(),icent,ipt));
//kforebbcw[ipt][icent]->Add(htemp);
hforebbcw[ipt][icent] = (TH1F*)f->Get(Form("hfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
//htemp = (TH1F*)f->Get(Form("hforeW%sbbcw_%d_%d",dire.Data(),icent,ipt));
//hforebbcw[ipt][icent]->Add(htemp);
kbackbbcw2[ipt][icent] = (TH1F*)f->Get(Form("kback%sbbcw2_%d_%d",dire.Data(),icent,ipt));
//htemp = (TH1F*)f->Get(Form("kbackW%sbbcw_%d_%d",dire.Data(),icent,ipt));
//kbackbbcw2[ipt][icent]->Add(htemp);
}
kforebbcw_centIn[ipt] = (TH1F*)kforebbcw[ipt][xcentmin]->Clone();
hforebbcw_centIn[ipt] = (TH1F*)hforebbcw[ipt][xcentmin]->Clone();
kbackbbcw2_centIn[ipt] = (TH1F*)kbackbbcw2[ipt][xcentmin]->Clone();
for(int icent=xcentmin+1; icent<xcentmax; icent++){
kforebbcw_centIn[ipt]->Add(kforebbcw[ipt][icent]);
hforebbcw_centIn[ipt]->Add(hforebbcw[ipt][icent]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2[ipt][icent]);
}
kforebbcw_centIn[ipt]->Rebin(2);
hforebbcw_centIn[ipt]->Rebin(2);
kbackbbcw2_centIn[ipt]->Rebin(2);
}

TH1F* hpp[npt];
TH1F* hbackpp[npt];
for(int ipt=0; ipt<1; ipt++){
  int ipt_add;
	if(ipt==0){
for(int ipt_add=0;ipt_add<1;ipt_add++){
kforebbcw_centIn[ipt]->Add(kforebbcw_centIn[ipt_add]);
hforebbcw_centIn[ipt]->Add(hforebbcw_centIn[ipt_add]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2_centIn[ipt_add]);
}
ptbin[ipt] = 1.0;
ptbin[ipt+1] = 3.0 ;
}
	if(ipt==1){
kforebbcw_centIn[ipt]->Reset("M");
hforebbcw_centIn[ipt]->Reset("M");
kbackbbcw2_centIn[ipt]->Reset("M");
for(int ipt_add=2;ipt_add<4;ipt_add++){
kforebbcw_centIn[ipt]->Add(kforebbcw_centIn[ipt_add]);
hforebbcw_centIn[ipt]->Add(hforebbcw_centIn[ipt_add]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2_centIn[ipt_add]);
}
ptbin[ipt] = 1.0 ;
ptbin[ipt+1] = 2.0 ;
}
	if(ipt==2){
kforebbcw_centIn[ipt]->Reset("M");
hforebbcw_centIn[ipt]->Reset("M");
kbackbbcw2_centIn[ipt]->Reset("M");
for(int ipt_add=4;ipt_add<6;ipt_add++){
kforebbcw_centIn[ipt]->Add(kforebbcw_centIn[ipt_add]);
hforebbcw_centIn[ipt]->Add(hforebbcw_centIn[ipt_add]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2_centIn[ipt_add]);
}
ptbin[ipt] = 2.0 ;
ptbin[ipt+1] = 3.0 ;
}
	if(ipt==3){
kforebbcw_centIn[ipt]->Reset("M");
hforebbcw_centIn[ipt]->Reset("M");
kbackbbcw2_centIn[ipt]->Reset("M");
for(int ipt_add=6;ipt_add<8;ipt_add++){
kforebbcw_centIn[ipt]->Add(kforebbcw_centIn[ipt_add]);
hforebbcw_centIn[ipt]->Add(hforebbcw_centIn[ipt_add]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2_centIn[ipt_add]);
}
ptbin[ipt] = 3.0 ;
ptbin[ipt+1] = 4.0 ;
}
  hpp[ipt] = (TH1F*)kforebbcw_centIn[ipt]->Clone();
  hbackpp[ipt] = (TH1F*)kbackbbcw2_centIn[ipt]->Clone();
  float nbackpp = hbackpp[ipt]->Integral()/2.0/PI;
  float nforepp = hpp[ipt]->Integral()/2.0/PI;
//float ntrig0 = centforedis_0->Integral(11,30);
   for(int i=0; i<20; i++){
     float pp_cont = 1.0*hpp[ipt]->GetBinContent(i+1);
     //float pp0_err = 1.0*hpp[ipt]->GetBinError(i+1);
     float weight2 = sqrt(1.0*hforebbcw_centIn[ipt]->GetBinContent(i+1));

     float backpp_cont = 1.0*hbackpp[ipt]->GetBinContent(i+1);
    
     float con = pp_cont/backpp_cont*nbackpp/nforepp;
     float err = weight2/backpp_cont*nbackpp/nforepp;

     hpp[ipt]->SetBinContent(i+1, con);
     hpp[ipt]->SetBinError(i+1, err);
   }

  c1->cd(ipt+1);
  gPad->SetTicks(1,1);
  

//  float ymax = 1.032;//hpp[ipt]->GetMaximum()*1.1;
//  float ymin = 0.968;//hpp[ipt]->GetMinimum()*0.9;
  float ymax = 1.04;
  //float ymax = 1.02;//hpp[ipt]->GetMaximum()*1.1;
  float ymin = 0.96;
  //float ymin = 0.98;//hpp[ipt]->GetMinimum()*0.9;

  hpp[ipt]->SetMinimum(ymin);
  hpp[ipt]->SetMaximum(ymax);

  hpp[ipt]->SetMarkerStyle(20);
  hpp[ipt]->SetMarkerSize(1.3);
  
  hpp[ipt]->SetMarkerColor(4);
  hpp[ipt]->GetYaxis()->SetLabelSize(0.04);
  hpp[ipt]->GetXaxis()->SetLabelSize(0.02);
  hpp[ipt]->GetYaxis()->SetNdivisions(505);
  //hpp[ipt]->GetYaxis()->SetTitleSize(0.6);
  //hpp[ipt]->GetXaxis()->SetTitleSize(0.0);
  hpp[ipt]->Draw("PE");
  TF1 *fun0 = new TF1("fun0","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)+2*[4]*cos(4*x))", -0.5*PI, 1.5*PI);

  fun0->SetLineColor(1);
  hpp[ipt]->Fit("fun0","R");


  TF1 *fun1 = new TF1("fun1","[0]*(1+2*[1]*cos(x))",   -0.5*PI, 1.5*PI);
  TF1 *fun2 = new TF1("fun2","[0]*(1+2*[1]*cos(2*x))", -0.5*PI, 1.5*PI);
  TF1 *fun3 = new TF1("fun3","[0]*(1+2*[1]*cos(3*x))", -0.5*PI, 1.5*PI);
  TF1 *fun4 = new TF1("fun4","[0]*(1+2*[1]*cos(4*x))", -0.5*PI, 1.5*PI);

TLegend *leg1 = new TLegend(0.62,0.40,0.82,0.45);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.03);
  leg1->SetBorderSize(0);
  leg1->AddEntry(fun0,"1+#Sigma2c_{n}cos(n#Delta#phi)","L");
  if(ipt==0)leg1->Draw();

  TLegend *leg2 = new TLegend(0.62,0.12,0.82,0.35);
  leg2->SetFillColor(10);
  leg2->SetLineStyle(4000);
  leg2->SetLineColor(10);
  leg2->SetLineWidth(0.);
  leg2->SetTextSize(0.03);
  leg2->SetBorderSize(0);
  leg2->AddEntry(fun1,"1+2c_{1}cos(#Delta#phi)","L");
  leg2->AddEntry(fun2,"1+2c_{2}cos(2#Delta#phi)","L");
  leg2->AddEntry(fun3,"1+2c_{3}cos(3#Delta#phi)","L");
  leg2->AddEntry(fun4,"1+2c_{4}cos(4#Delta#phi)","L");
  if(ipt==0)leg2->Draw();

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

  fun1->Draw("same");
  fun2->Draw("same");
  fun3->Draw("same");
  fun4->Draw("same");

  hpp[ipt]->GetXaxis()->SetTitle("#Delta#phi");
  hpp[ipt]->GetXaxis()->SetTitleSize(0.04);
  hpp[ipt]->GetXaxis()->SetLabelSize(0.04);
  hpp[ipt]->GetXaxis()->CenterTitle();
  hpp[ipt]->GetXaxis()->SetTitleOffset(1.0);
  if(ipt==0)
  hpp[ipt]->GetYaxis()->SetTitle("C(d#phi)");
  else
  hpp[ipt]->GetYaxis()->SetLabelSize(0);
  hpp[ipt]->GetYaxis()->SetTitleSize(0.04);
  hpp[ipt]->GetYaxis()->CenterTitle();
  hpp[ipt]->GetYaxis()->SetTitleOffset(1.);
  hpp[ipt]->GetYaxis()->SetLabelSize(0.04);
  if(ipt==0)
  hpp[ipt]->GetYaxis()->SetTitle("C(#Delta#phi)");
  TLatex *t=new TLatex(-1.0,0.88*(ymax-ymin)+ymin, Form("%d)",icent+1));
  TLatex *t=new TLatex(-1.0,0.88*(ymax-ymin)+ymin, Form("%d)",ipt+1));
  t->SetTextSize(0.06);
 // t->Draw();

  TMarker *m0 = new TMarker(-0.9, 0.91*(ymax-ymin)+ymin, 20);
  m0->SetMarkerSize(1.3);
  m0->SetMarkerColor(4);
  m0->Draw();

  if(ipt==0){
  TLatex *t=new TLatex(-0.8,0.88*(ymax-ymin)+ymin, Form("Run15 p+Al 200GeV"));
  t->SetTextSize(0.04);
  t->Draw();
  TLatex *t=new TLatex(2.4,0.88*(ymax-ymin)+ymin, Form("%.1f < p_{T}^{trig} < %.1f (GeV/c)",ptbin[ipt],ptbin[ipt+1]));
  }
  else
  TLatex *t=new TLatex(-0.4,0.88*(ymax-ymin)+ymin, Form("%.1f < p_{T}^{trig} < %.1f (GeV/c)",ptbin[ipt],ptbin[ipt+1]));
 // TLatex *t=new TLatex(-0.1,0.88*(ymax-ymin)+ymin, Form("p+p %.1f < p_{T} < %.1f (GeV/c)",ptbin[ipt],ptbin[ipt+1]));
  t->SetTextSize(0.03);
  t->Draw();

  TLatex *t=new TLatex(-1.2,0.78*(ymax-ymin)+ymin, Form("%.f - %.f%%",centmin,centmax));
  //TLatex *t=new TLatex(-1.2,0.78*(ymax-ymin)+ymin, Form("Central event with BBC_S charge >= 50"));
  //TLatex *t=new TLatex(-1.2,0.78*(ymax-ymin)+ymin, Form("Minimum Bias"));
  t->SetTextSize(0.04);
  if(ipt==0)t->Draw();
  TLatex *t=new TLatex(2,0.78*(ymax-ymin)+ymin, Form("-3.0 <#eta_{trig} < -1.0"));
  t->SetTextSize(0.04);
  if(ipt==0)t->Draw();

if(dire=="north")
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "3.0<#eta_{asso}<3.9");
  //TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "1.5<#eta_{asso}<2.5");
else
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "-3.9<#eta_{asso}<-3.0");
  //TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "-2.5<#eta_{asso}<-1.5");
  t->SetTextSize(0.04);
  t->SetTextColor(2);
  if(ipt==0)t->Draw();

}
  c1->cd();
/*
  TLatex *tx1=new TLatex(0.30,0.028, "#Delta#phi");
  tx1->SetTextSize(0.04);
  tx1->Draw();

  TLatex *tx2=new TLatex(0.75,0.028, "#Delta#phi");
  tx2->SetTextSize(0.04);
  tx2->Draw();

  TLatex *ty1 = new TLatex(0.052,0.78, "C(#Delta#phi)");
  ty1->SetTextSize(0.04);
  ty1->SetTextAngle(90);
  ty1->Draw();

  TLatex *ty2 = new TLatex(0.052,0.50, "C(#Delta#phi)");
  ty2->SetTextSize(0.04);
  ty2->SetTextAngle(90);
  ty2->Draw();

  TLatex *ty3 = new TLatex(0.052,0.20, "C(#Delta#phi)");
  ty3->SetTextSize(0.04);
  ty3->SetTextAngle(90);
  ty3->Draw();
*/
  c1->Print(Form("fig/ridge_bbcfvtx_%s_%.1f_%.1f.gif",dire.Data(),centmin,centmax));
  c1->Print(Form("fig/ridge_bbcfvtx_%s_%.1f_%.1f.png",dire.Data(),centmin,centmax));
  c1->Print(Form("fig/ridge_bbcfvtx_%s_%.1f_%.1f.pdf",dire.Data(),centmin,centmax));
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
