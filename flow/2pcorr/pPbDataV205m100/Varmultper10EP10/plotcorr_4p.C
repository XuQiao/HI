#include "root_setting.h"
void plotcorr_4p(){
//TString dire = "north";
TString dire = "south";
TString corr = "cntcnt";
c1 = new TCanvas("c1"," ",1400,460);
makeMultiPanelCanvas(c1,5,2,0,0,0.25,0.2,0.03);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
  float const PI = acos(-1.0);
  const int ncent = 10;
  const int npt = 10;
  TFile *f=TFile::Open(Form("mergedAnapAu%s.root",corr.Data()));
  ofstream fout("YjetpAu.dat");
TH1F* kforebbcw[ncent][npt];
TH1F* hforebbcw[ncent][npt];
TH1F* kbackbbcw2[ncent][npt];
TH1F* kforebbcw_ptIn[ncent];
TH1F* hforebbcw_ptIn[ncent];
TH1F* kbackbbcw2_ptIn[ncent];
double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
//double centbin[ncent+1] = {0,1,5,10,20,30,40,60,84}; //fvtx or cent
double centbin[ncent+1] = {100.00,97.51,92.65,85.92,78.00,69.62,53.55,33.84,15.74,5.89,0.00};
//double bbccentbin[ncent+1] = {200.0,39.3,27.0,21.4,15.6,9.5,5.8,2.9,0.0}; //fvtx mb cent
//double centbin[ncent] = {108,83.3333,75,56.6667,48.6667,40,29.6667,21.6667,0}; //fvtx or
// double centbin[ncent] = {110,85.3333,77.3333,59.3333,51.6667,43,32.6667,25.3333,0};   //fvtx south
double bbccentbin[ncent+1] = {0,1,2,3,4,5,7,10,14,18,40};
//  double centbin[ncent] = {152.333,118.333,102.333,69.6667,57.3333,45,31,21,0};   //fvtx and
/*double ptmin = 1.0; double ptmax = 3.0;//has to be boundary
for(int ipt=0; ipt<npt; ipt++){
if(ptmin >= ptbin[ipt] && ptmin < ptbin[ipt+1]){int xptmin = ipt; continue;}
if(ptmax >= ptbin[ipt] && ptmax < ptbin[ipt+1]){int xptmax = ipt; continue;}
}*/
for(int ipt_a=0;ipt_a<1;ipt_a++){
int xptmin = ipt_a*4+2;
int xptmax = (ipt_a+1)*4+2;
double ptmin = ptbin[xptmin];
double ptmax = ptbin[xptmax];

for(int icent=0; icent<ncent; icent++){
for(int ipt=0; ipt<npt; ipt++){
  //  if(icent==0){
kforebbcw[icent][ipt] = (TH1F*)f->Get(Form("kfore%sbbc_%d_%d",dire.Data(),icent,ipt));
hforebbcw[icent][ipt] = (TH1F*)f->Get(Form("hfore%sbbc_%d_%d",dire.Data(),icent,ipt));
kbackbbcw2[icent][ipt] = (TH1F*)f->Get(Form("kback%sbbc2_%d_%d",dire.Data(),icent,ipt));
//    }
//    else{
//kforebbcw[icent][ipt] = (TH1F*)fmb->Get(Form("kfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
//hforebbcw[icent][ipt] = (TH1F*)fmb->Get(Form("hfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
//kbackbbcw2[icent][ipt] = (TH1F*)fmb->Get(Form("kback%sbbcw2_%d_%d",dire.Data(),icent,ipt));
   // }
}
kforebbcw_ptIn[icent] = (TH1F*)kforebbcw[icent][xptmin]->Clone();
hforebbcw_ptIn[icent] = (TH1F*)hforebbcw[icent][xptmin]->Clone();
kbackbbcw2_ptIn[icent] = (TH1F*)kbackbbcw2[icent][xptmin]->Clone();
for(int ipt=xptmin+1; ipt<xptmax; ipt++){
kforebbcw_ptIn[icent]->Add(kforebbcw[icent][ipt]);
hforebbcw_ptIn[icent]->Add(hforebbcw[icent][ipt]);
kbackbbcw2_ptIn[icent]->Add(kbackbbcw2[icent][ipt]);
}
kforebbcw_ptIn[icent]->Rebin(2);
hforebbcw_ptIn[icent]->Rebin(2);
kbackbbcw2_ptIn[icent]->Rebin(2);
}

TH1F* hpp[ncent];
TH1F* hbackpp[ncent];
 
  int icentstart = 0;
for(int icent=icentstart; icent<ncent; icent++){
	
if(icent==icentstart){
for(int icent_add=0;icent_add<icentstart;icent_add++){
kforebbcw_ptIn[icent]->Add(kforebbcw_ptIn[icent_add]);
hforebbcw_ptIn[icent]->Add(hforebbcw_ptIn[icent_add]);
kbackbbcw2_ptIn[icent]->Add(kbackbbcw2_ptIn[icent_add]);
}
centbin[icent] = centbin[0];
bbccentbin[icent] = bbccentbin[0];
}

/*
	if(icent==6){
for(int icent_add=6;icent_add<ncent;icent_add++){
kforebbcw_ptIn[icent]->Add(kforebbcw_ptIn[icent_add]);
hforebbcw_ptIn[icent]->Add(hforebbcw_ptIn[icent_add]);
kbackbbcw2_ptIn[icent]->Add(kbackbbcw2_ptIn[icent_add]);
}
centbin[icent] = 800;
}
*/
  hpp[icent] = (TH1F*)kforebbcw_ptIn[icent]->Clone();
  hbackpp[icent] = (TH1F*)kbackbbcw2_ptIn[icent]->Clone();
  float nbackpp = hbackpp[icent]->Integral()/2.0/PI;
  float nforepp = hpp[icent]->Integral()/2.0/PI;
  //float ntrig0 = ptforedis_0->Integral(11,30);
   for(int i=0; i<20; i++){
     float pp_cont = 1.0*hpp[icent]->GetBinContent(i+1);
     //float pp0_err = 1.0*hpp[icent]->GetBinError(i+1);
     float weight2 = sqrt(1.0*hforebbcw_ptIn[icent]->GetBinContent(i+1));

     float backpp_cont = 1.0*hbackpp[icent]->GetBinContent(i+1);
    
     float con = pp_cont/backpp_cont*nbackpp/nforepp;
     float err = weight2/backpp_cont*nbackpp/nforepp;

     hpp[icent]->SetBinContent(i+1, con);
     hpp[icent]->SetBinError(i+1, err);
   }

  c1->cd(icent-icentstart+1);
  gPad->SetTicks(1,1);

  if(corr.Contains("bbc")){
  float ymax = 1.08-0.0;//hpp[icent]->GetMaximum()*1.1;
  float ymin = 0.92+0.0;//hpp[icent]->GetMinimum()*0.9;
  }
    else{
  float ymax = 4.0-0.0;//hpp[icent]->GetMaximum()*1.1;
  float ymin = 0.0+0.0;//hpp[icent]->GetMinimum()*0.9;
  }

  hpp[icent]->SetMinimum(ymin);
  hpp[icent]->SetMaximum(ymax);

  hpp[icent]->SetMarkerStyle(20);
  hpp[icent]->SetMarkerSize(1.3);
  
  hpp[icent]->SetMarkerColor(4);
  hpp[icent]->GetYaxis()->SetLabelSize(0.08);
  hpp[icent]->GetXaxis()->SetLabelSize(0.08);
  hpp[icent]->GetYaxis()->SetNdivisions(505);
  //hpp[icent]->GetYaxis()->SetTitleSize(0.6);
  //hpp[icent]->GetXaxis()->SetTitleSize(0.0);
  hpp[icent]->Draw("PE");
  if(corr.Contains("bbc")){
  TF1 *fun0 = new TF1("fun0","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)+2*[4]*cos(4*x))", -0.5*PI, 1.5*PI);

  fun0->SetLineColor(1);
  hpp[icent]->Fit("fun0","R");


  TF1 *fun1 = new TF1("fun1","[0]*(1+2*[1]*cos(x))",   -0.5*PI, 1.5*PI);
  TF1 *fun2 = new TF1("fun2","[0]*(1+2*[1]*cos(2*x))", -0.5*PI, 1.5*PI);
  TF1 *fun3 = new TF1("fun3","[0]*(1+2*[1]*cos(3*x))", -0.5*PI, 1.5*PI);
  TF1 *fun4 = new TF1("fun4","[0]*(1+2*[1]*cos(4*x))", -0.5*PI, 1.5*PI);

TLegend *leg1 = new TLegend(0.42,0.18,0.82,0.26);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.08);
  leg1->SetBorderSize(0);
  leg1->AddEntry(fun0,"1+#Sigma2c_{n}cos(n#Delta#phi)","L");
  if(icent-icentstart==1)leg1->Draw();

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
  if(icent-icentstart==1)leg2->Draw();

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
  }
else{
  //TF1 *fun0 = new TF1("fun0","[0]*x*x+[1]*x+[2]",0.1,2);
  TF1 *fun0 = new TF1("fun0","gaus(0)+gaus(3)",-1.2,PI/2);
  fun0->SetParameters(3,0,1,0.2,0,4);
  hpp[icent]->Fit("fun0","R");
  hpp[icent]->DrawCopy();
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
   for(int i=0; i<21; i++){
  hpp[icent]->SetBinContent(i,hpp[icent]->GetBinContent(i)-zyam);
   }
   hpp[icent]->SetMarkerColor(2);
   hpp[icent]->Draw("same");
//   fout<<hpp[icent]->Integral(hpp[icent]->FindBin(-1.2),hpp[icent]->FindBin(PI/2))<<"\t";
   fout<<0.95*fabs(fun1->GetParameter(0)*fun1->GetParameter(2)*sqrt(2*PI))<<endl;
}
  hpp[icent]->GetXaxis()->SetTitle("#Delta#phi");
  hpp[icent]->GetXaxis()->SetTitleSize(0.08);
  hpp[icent]->GetXaxis()->SetLabelSize(0.08);
  hpp[icent]->GetXaxis()->CenterTitle();
  hpp[icent]->GetXaxis()->SetTitleOffset(1.1);
  hpp[icent]->GetYaxis()->SetTitle("C(d#phi)");
  hpp[icent]->GetYaxis()->SetTitleSize(0.08);
  hpp[icent]->GetYaxis()->CenterTitle();
  hpp[icent]->GetYaxis()->SetTitleOffset(1.3);
  hpp[icent]->GetYaxis()->SetLabelSize(0.08);
  hpp[icent]->GetYaxis()->SetTitle("C(#Delta#phi)");
  TLatex *t=new TLatex(-1.0,0.88*(ymax-ymin)+ymin, Form("%d)",icent+1));
  t->SetTextSize(0.08);
 // t->Draw();

  TMarker *m0 = new TMarker(-0.3, 0.91*(ymax-ymin)+ymin, 20);

  if(icent-icentstart==0){
  TMarker *m0 = new TMarker(-1.1, 0.91*(ymax-ymin)+ymin, 20);
  TLatex *t=new TLatex(-1.1,0.88*(ymax-ymin)+ymin, Form("p+Au HIJING 200GeV Nch= %.f",bbccentbin[icent+1]));
  }
  else{
  TLatex *t=new TLatex(-0.9,0.88*(ymax-ymin)+ymin, Form("%.f <Nch<= %.f",bbccentbin[icent],bbccentbin[icent+1]));
  }
  t->SetTextSize(0.08);
  t->Draw();
  if(icent-icentstart==5){
    TLatex *t=new TLatex(-1.1,0.24*(ymax-ymin)+ymin, Form("with weight = particle Energy"));
    t->SetTextSize(0.07);
  //  t->Draw();
  }
  m0->SetMarkerSize(1.3);
  m0->SetMarkerColor(4);
  //m0->Draw();
  
  TLatex *t=new TLatex(-0.4,0.08*(ymax-ymin)+ymin, Form("%.1f < cent < %.1f%%",centbin[icent],centbin[icent+1]));
  t->SetTextSize(0.08);
  t->Draw();

  TLatex *t=new TLatex(-1.2,0.78*(ymax-ymin)+ymin, Form("%.1f<p_{T,trig}<%.1f GeV/c",ptmin,ptmax));
  t->SetTextSize(0.08);
  if(icent-icentstart==0)t->Draw();

  TLatex *t=new TLatex(-1.0,0.66*(ymax-ymin)+ymin, "|#eta_{trig}|<0.35");
  t->SetTextSize(0.08);
  if(icent-icentstart==0)t->Draw();

if(dire=="north")
  TLatex *t=new TLatex(-1.2,0.24*(ymax-ymin)+ymin, "3.0<#eta_{asso}<4.0");
else
  TLatex *t=new TLatex(-1.2,0.24*(ymax-ymin)+ymin, "-4.0<#eta_{asso}<-3.0");
  t->SetTextSize(0.08);
  t->SetTextColor(2);
  if(icent-icentstart==0)t->Draw();

}
  c1->Print(Form("fig/pAuridge%s_%s_%.1f_%.1f.gif",corr.Data(),dire.Data(),ptmin,ptmax));
  c1->Print(Form("fig/pAuridge%s_%s_%.1f_%.1f.pdf",corr.Data(),dire.Data(),ptmin,ptmax));
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
