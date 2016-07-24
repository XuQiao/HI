#include "root_setting.h"
void plotcorr2D_4p_pt(){
//TString dire = "north";
TString dire = "south";
TString coll = "pAu";
TString corr = "full";
c1 = new TCanvas("c1"," ",800,460);
makeMultiPanelCanvas(c1,3,2,0,0,0.25,0.2,0.03);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
  float const PI = acos(-1.0);
  const int npt = 10;
  const int ncent = 10;
  TFile *f=TFile::Open(Form("mergedAna%s%s.root",coll.Data(),corr.Data()));
  ofstream fout;
TH2F* kforebbcw[npt][ncent];
TH2F* hforebbcw[npt][ncent];
TH2F* kbackbbcw2[npt][ncent];
TH2F* kforebbcw_centIn[npt];
TH2F* hforebbcw_centIn[npt];
TH2F* kbackbbcw2_centIn[npt];
double ptbin[ncent+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
//double centbin[ncent+1] = {0,1,5,10,20,30,40,60,84}; //fvtx or cent
double centbin[ncent+1] = {100.00,97.51,92.65,85.92,78.00,69.62,53.55,33.84,15.74,5.89,0.00};
//double bbccentbin[ncent+1] = {200.0,39.3,27.0,21.4,15.6,9.5,5.8,2.9,0.0}; //fvtx mb cent
//double centbin[ncent] = {108,83.3333,75,56.6667,48.6667,40,29.6667,21.6667,0}; //fvtx or
// double centbin[ncent] = {110,85.3333,77.3333,59.3333,51.6667,43,32.6667,25.3333,0};   //fvtx south
double bbccentbin[ncent+1] = {0,1,2,3,4,5,7,10,14,18,40};
//  double centbin[ncent] = {152.333,118.333,102.333,69.6667,57.3333,45,31,21,0};   //fvtx and
/*double centmin = 1.0; double centmax = 3.0;//has to be boundary
for(int ipt=0; ipt<ncent; ipt++){
if(centmin >= centbin[ipt] && centmin < centbin[ipt+1]){int xcentmin = ipt; continue;}
if(centmax >= centbin[ipt] && centmax < centbin[ipt+1]){int xcentmax = ipt; continue;}
}*/
for(int icent_a=0;icent_a<1;icent_a++){
    if(coll=="pAu")
int xcentmin = icent_a*1+9;
    else
int xcentmin = icent_a*1;
int xcentmax = (icent_a+1)*10;
double centmin = centbin[xcentmin];
double centmax = centbin[xcentmax];

for(int ipt=0; ipt<npt; ipt++){
for(int icent=0; icent<ncent; icent++){
kforebbcw[ipt][icent] = (TH2F*)f->Get(Form("kfore%setabbc_%d_%d",dire.Data(),icent,ipt));
hforebbcw[ipt][icent] = (TH2F*)f->Get(Form("kfore%setabbc_%d_%d",dire.Data(),icent,ipt));
kbackbbcw2[ipt][icent] = (TH2F*)f->Get(Form("kback%setabbc2_%d_%d",dire.Data(),icent,ipt));
}
kforebbcw_centIn[ipt] = (TH2F*)kforebbcw[ipt][xcentmin]->Clone();
hforebbcw_centIn[ipt] = (TH2F*)hforebbcw[ipt][xcentmin]->Clone();
kbackbbcw2_centIn[ipt] = (TH2F*)kbackbbcw2[ipt][xcentmin]->Clone();
for(int icent=xcentmin+1; icent<xcentmax; icent++){
kforebbcw_centIn[ipt]->Add(kforebbcw[ipt][icent]);
hforebbcw_centIn[ipt]->Add(hforebbcw[ipt][icent]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2[ipt][icent]);
}
//kforebbcw_centIn[ipt]->RebinX(2);
//hforebbcw_centIn[ipt]->RebinX(2);
//kbackbbcw2_centIn[ipt]->RebinX(2);
}

TH2F* hpp[npt];
TH2F* hbackpp[npt];
  
int iptstart = 0;
for(int ipt=iptstart; ipt<6; ipt++){
	
if(ipt==iptstart){
for(int ipt_add=0;ipt_add<iptstart;ipt_add++){
kforebbcw_centIn[ipt]->Add(kforebbcw_centIn[ipt_add]);
hforebbcw_centIn[ipt]->Add(hforebbcw_centIn[ipt_add]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2_centIn[ipt_add]);
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
  hpp[ipt] = (TH2F*)kforebbcw_centIn[ipt]->Clone();
  hbackpp[ipt] = (TH2F*)kbackbbcw2_centIn[ipt]->Clone();
  float nbackpp = hbackpp[ipt]->Integral()/2.0/PI;
  float nforepp = hpp[ipt]->Integral()/2.0/PI;
  //float ntrig0 = centforedis_0->Integral(11,30);
   for(int i=0; i<40*100; i++){
     float pp_cont = 1.0*hpp[ipt]->GetBinContent(i+1);
     //float pp0_err = 1.0*hpp[ipt]->GetBinError(i+1);
     float weight2 = sqrt(1.0*hforebbcw_centIn[ipt]->GetBinContent(i+1));

     float backpp_cont = 1.0*hbackpp[ipt]->GetBinContent(i+1);
    
     float con = pp_cont/backpp_cont*nbackpp/nforepp;
     //float err = weight2/backpp_cont*nbackpp/nforepp;

     hpp[ipt]->SetBinContent(i+1, con);
     hpp[ipt]->SetBinError(i+1, 0);
   }

  c1->cd(ipt-iptstart+1);
  c1->cd(ipt-iptstart+1)->SetTheta(60.839);
  c1->cd(ipt-iptstart+1)->SetPhi(38.0172);

  gPad->SetTicks(1,1);
if(corr.Contains("bbc")){
    if(coll=="pAu"){
  float ymax = 1.1-0.0;//hpp[ipt]->GetMaximum()*1.1;
  float ymin = 0.9+0.0;//hpp[ipt]->GetMinimum()*0.9;
    }
    else{
  float ymax = 1.3-0.0;//hpp[ipt]->GetMaximum()*1.1;
  float ymin = 0.7+0.0;//hpp[ipt]->GetMinimum()*0.9;
    }
}
else if(corr.Contains("cnt")){
  float ymax = 7.0-0.0;//hpp[ipt]->GetMaximum()*1.1;
  float ymin = 0.0+0.0;//hpp[ipt]->GetMinimum()*0.9;
}
else{
  float ymax = 2.0-0.0;//hpp[ipt]->GetMaximum()*1.1;
  float ymin = 0.8+0.0;//hpp[ipt]->GetMinimum()*0.9;
}
  hpp[ipt]->SetMinimum(ymin);
  hpp[ipt]->SetMaximum(ymax);
  if(corr.Contains("bbc"))
  hpp[ipt]->GetXaxis()->SetRangeUser(2.7,3.7);
  else if(corr.Contains("cnt"))
  hpp[ipt]->GetXaxis()->SetRangeUser(-0.7,0.7);
  else
  hpp[ipt]->GetXaxis()->SetRangeUser(-7,7);
  hpp[ipt]->GetYaxis()->SetRangeUser(-1.5,4.2);

  hpp[ipt]->SetMarkerStyle(20);
  hpp[ipt]->SetMarkerSize(1.3);
  
  hpp[ipt]->SetMarkerColor(4);
  hpp[ipt]->GetYaxis()->SetLabelSize(0.08);
  hpp[ipt]->GetXaxis()->SetLabelSize(0.08);
  hpp[ipt]->GetYaxis()->SetNdivisions(505);
  //hpp[ipt]->GetYaxis()->SetTitleSize(0.6);
  //hpp[ipt]->GetXaxis()->SetTitleSize(0.0);
  hpp[ipt]->Draw("surf1");
  hpp[ipt]->GetXaxis()->SetTitle("#Delta#eta");
  hpp[ipt]->GetYaxis()->SetTitle("#Delta#phi");
  hpp[ipt]->GetXaxis()->SetTitleSize(0.08);
  hpp[ipt]->GetXaxis()->SetLabelSize(0.08);
  hpp[ipt]->GetXaxis()->CenterTitle();
  hpp[ipt]->GetXaxis()->SetTitleOffset(1.1);
  hpp[ipt]->GetZaxis()->SetTitle("C(d#phi)");
  hpp[ipt]->GetYaxis()->SetTitle("#Delta#phi");
  hpp[ipt]->GetYaxis()->SetTitleSize(0.08);
  hpp[ipt]->GetYaxis()->CenterTitle();
  hpp[ipt]->GetYaxis()->SetTitleOffset(1.3);
  hpp[ipt]->GetYaxis()->SetLabelSize(0.08);
  hpp[ipt]->GetZaxis()->SetTitle("C(#Delta#phi)");
 // TLatex *t=new TLatex(-1.0,0.88*(ymax-ymin)+ymin, Form("%d)",ipt+1));
 // t->SetTextSize(0.08);
 // t->Draw();

  TMarker *m0 = new TMarker(-0.3, 0.91*(ymax-ymin)+ymin, 20);

  if(ipt-iptstart==0){
  TMarker *m0 = new TMarker(-1.1, 0.91*(ymax-ymin)+ymin, 20);
  TLatex *t=new TLatex(-1.1,0.88*(ymax-ymin)+ymin, Form("%s HIJING 200GeV ",coll.Data()));
  t->SetTextSize(0.08);
  t->Draw();
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

  TLatex *t=new TLatex(-1.2,0.08*(ymax-ymin)+ymin, Form("%.1f< cent <%.1f%%",centmax,centmin));
  t->SetTextSize(0.08);
  if(ipt-iptstart==2)t->Draw();

  TLatex *t=new TLatex(-1.0,0.08*(ymax-ymin)+ymin, "|#eta_{trig}|<0.35");
  t->SetTextSize(0.08);
  //if(ipt-iptstart==0)t->Draw();

if(dire=="north")
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "3.0<#eta_{asso}<4.0");
else
if(corr.Contains("bbc"))
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "-4.0<#eta_{asso}<-3.0");
  else
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "|#eta_{asso}|<0.35");
  t->SetTextSize(0.08);
  t->SetTextColor(2);
//  if(ipt-iptstart==0)t->Draw();

}
  c1->Print(Form("fig/%s%sridge2D_%s_%.1f_%.1f.gif",coll.Data(),corr.Data(),dire.Data(),centmin,centmax));
  c1->Print(Form("fig/%s%sridge2D_%s_%.1f_%.1f.pdf",coll.Data(),corr.Data(),dire.Data(),centmin,centmax));
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
