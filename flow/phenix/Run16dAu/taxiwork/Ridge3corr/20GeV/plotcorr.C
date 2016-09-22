void plotcorr(){
//TString dire = "north";
TString dire = "south";
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
  float const PI = acos(-1.0);
  const int ncent = 8;
  const int npt = 10;
//  TFile *f=TFile::Open("merged_AnapAufvtxor.root");
//  TFile *f=TFile::Open("merged_AnapAufvtxand.root");
//  TFile *f=TFile::Open("merged_AnapAufvtxsouth.root");
//  TFile *fmb=TFile::Open("merged_PUlow_AnapAumbst.root");
//  TFile *fmb=TFile::Open("merged_PUlow_AnapAumbcentral.root");
  TFile *f=TFile::Open("merged_AnapAumbcentral.root");
//  TFile *f=TFile::Open("merged_AnapAumbst.root");
//  TFile *f=TFile::Open("merged_AnapAufvtxorcent.root");
TH1F* kforebbcw[ncent][npt];
TH1F* hforebbcw[ncent][npt];
TH1F* kbackbbcw2[ncent][npt];
TH1F* kforebbcw_ptIn[ncent];
TH1F* hforebbcw_ptIn[ncent];
TH1F* kbackbbcw2_ptIn[ncent];
double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
double centbin[ncent+1] = {0,1,5,10,20,30,40,60,100}; //fvtx or cent
//double bbccentbin[ncent+1] = {200.0,39.3,27.0,21.4,15.6,9.5,5.8,2.9,0.0}; //fvtx mb cent
//double centbin[ncent] = {108,83.3333,75,56.6667,48.6667,40,29.6667,21.6667,0}; //fvtx or
// double centbin[ncent] = {110,85.3333,77.3333,59.3333,51.6667,43,32.6667,25.3333,0};   //fvtx south
double bbccentbin[ncent+1] = {800,56.3333,44.3333,39.3333,27,21.3333,15.6667,9.66667,5.66667}; //mb
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
kforebbcw[icent][ipt] = (TH1F*)f->Get(Form("kfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
hforebbcw[icent][ipt] = (TH1F*)f->Get(Form("hfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
kbackbbcw2[icent][ipt] = (TH1F*)f->Get(Form("kback%sbbcw2_%d_%d",dire.Data(),icent,ipt));
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
  c1=new TCanvas("c1","c1", 600,750);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
 
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

  float ytop=0.97, ybot=0.07;
  float bot_r = 0.082, top_r = 0.042;
  float ny = 1.0/(1-bot_r-0.02)+1.0/(1-top_r-0.02)+1.0/(1.0-0.04);
  float ygap = (ytop-ybot)/ny;
  float ygap1 = ytop - 1.0/(1-top_r-0.02)*ygap;
  float ygap2 = ytop - 1.0/(1-top_r-0.02)*ygap - 1.0/(1.0-0.04)*ygap;
  float xleft=0.07, xright=0.98, xmiddle=1.0/(0.84/0.85+1.0)*(xright-xleft)+xleft;
  c1->cd();
  icent-=icentstart;
  if(icent==0)TPad *c1_ = new TPad(Form("c1_%d",icent), Form("c1_%d",icent),xleft,ygap1,xmiddle,ytop);
  if(icent==1)TPad *c1_ = new TPad(Form("c1_%d",icent), Form("c1_%d",icent),xmiddle,ygap1,xright,ytop);
  if(icent==2)TPad *c1_ = new TPad(Form("c1_%d",icent), Form("c1_%d",icent),xleft,ygap2,xmiddle,ygap1);
  if(icent==3)TPad *c1_ = new TPad(Form("c1_%d",icent), Form("c1_%d",icent),xmiddle,ygap2,xright,ygap1);
  if(icent==4)TPad *c1_ = new TPad(Form("c1_%d",icent), Form("c1_%d",icent),xleft,ybot,xmiddle,ygap2);
  if(icent==5)TPad *c1_ = new TPad(Form("c1_%d",icent), Form("c1_%d",icent),xmiddle,ybot,xright,ygap2);
  icent+=icentstart;
  c1_->SetFillColor(10);
  c1_->SetFillColor(10);
  c1_->SetBorderMode(0);
  c1_->SetBorderSize(2);
  c1_->SetFrameFillColor(0);
  c1_->SetFrameBorderMode(0);
  c1_->SetFrameBorderMode(0);
  c1_->Draw();
  c1_->cd();
  //c1_1->Range(-0.658827,-4.65723,3.04092,1.94982);
  //c1_1->SetFillColor(10);
  //c1_1->SetBorderMode(0);
  //c1_1->SetBorderSize(2);
  //c1_1->SetLogy();
  gPad->SetLeftMargin(0.130);
  gPad->SetRightMargin(0.02);
 // gPad->SetBottomMargin(0.02);
  gPad->SetTopMargin(top_r);
  gPad->SetTicks(1,1);
  

  float ymax = 1.048;//hpp[icent]->GetMaximum()*1.1;
  float ymin = 0.952;//hpp[icent]->GetMinimum()*0.9;

  hpp[icent]->SetMinimum(ymin);
  hpp[icent]->SetMaximum(ymax);

  hpp[icent]->SetMarkerStyle(20);
  hpp[icent]->SetMarkerSize(1.3);
  
  hpp[icent]->SetMarkerColor(4);
  hpp[icent]->GetYaxis()->SetLabelSize(0.08);
  hpp[icent]->GetXaxis()->SetLabelSize(0.08);
  //hpp[icent]->GetYaxis()->SetTitleSize(0.6);
  //hpp[icent]->GetXaxis()->SetTitleSize(0.0);
  hpp[icent]->Draw("PE");
  TF1 *fun0 = new TF1("fun0","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)+2*[4]*cos(4*x))", -0.5*PI, 1.5*PI);

  fun0->SetLineColor(1);
  hpp[icent]->Fit("fun0","R");


  TF1 *fun1 = new TF1("fun1","[0]*(1+2*[1]*cos(x))",   -0.5*PI, 1.5*PI);
  TF1 *fun2 = new TF1("fun2","[0]*(1+2*[1]*cos(2*x))", -0.5*PI, 1.5*PI);
  TF1 *fun3 = new TF1("fun3","[0]*(1+2*[1]*cos(3*x))", -0.5*PI, 1.5*PI);
  TF1 *fun4 = new TF1("fun4","[0]*(1+2*[1]*cos(4*x))", -0.5*PI, 1.5*PI);

TLegend *leg1 = new TLegend(0.22,0.12,0.82,0.22);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.08);
  leg1->SetBorderSize(0);
  leg1->AddEntry(fun0,"1+#Sigma2c_{n}cos(n#Delta#phi)","L");
  if(icent-icentstart==1)leg1->Draw();

  TLegend *leg2 = new TLegend(0.18,0.52,0.42,0.82);
  leg2->SetFillColor(10);
  leg2->SetLineStyle(4000);
  leg2->SetLineColor(10);
  leg2->SetLineWidth(0.);
  leg2->SetTextSize(0.07);
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

  TLatex *t=new TLatex(-1.0,0.88*(ymax-ymin)+ymin, Form("%d)",icent+1));
  t->SetTextSize(0.08);
  t->Draw();

  TMarker *m0 = new TMarker(-0.3, 0.91*(ymax-ymin)+ymin, 20);
  m0->SetMarkerSize(1.3);
  m0->SetMarkerColor(4);
  m0->Draw();

  TLatex *t=new TLatex(-0.0,0.88*(ymax-ymin)+ymin, Form("p+Au %.1f - %.1f%%",centbin[icent],centbin[icent+1]));
  t->SetTextSize(0.08);
  t->Draw();
  
  TLatex *t=new TLatex(-0.0,0.78*(ymax-ymin)+ymin, Form("%.1f < BBcQs < %.1f",bbccentbin[icent],bbccentbin[icent+1]));
  t->SetTextSize(0.08);
//  t->Draw();

  TLatex *t=new TLatex(-1.2,0.78*(ymax-ymin)+ymin, Form("%.1f<p_{T,trig}<%.1f GeV/c",ptmin,ptmax));
  t->SetTextSize(0.08);
  if(icent-icentstart==0)t->Draw();

  TLatex *t=new TLatex(-1.0,0.66*(ymax-ymin)+ymin, "|#eta_{trig}|<0.35");
  t->SetTextSize(0.08);
  if(icent-icentstart==0)t->Draw();

if(dire=="north")
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "3.1<#eta_{asso}<3.7");
else
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "-3.7<#eta_{asso}<-3.1");
  t->SetTextSize(0.08);
  t->SetTextColor(2);
  if(icent-icentstart==0)t->Draw();

}
  c1->cd();

  TLatex *tx1=new TLatex(0.30,0.048, "#Delta#phi");
  tx1->SetTextSize(0.04);
  tx1->Draw();

  TLatex *tx2=new TLatex(0.75,0.048, "#Delta#phi");
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

  c1->Print(Form("fig/ridge_%s_%.1f_%.1f.png",dire.Data(),ptmin,ptmax));
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
