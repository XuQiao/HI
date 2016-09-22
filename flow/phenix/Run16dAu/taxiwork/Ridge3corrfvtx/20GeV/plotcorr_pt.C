void plotcorr_pt(){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
TString dire = "north";
//TString dire = "south";
  float const PI = acos(-1.0);
const int ncent = 4;
const int npt = 10;
  TFile *f=TFile::Open("../../../work/200GeV/output_cntbbc.root");
TH1F* kforebbcw[npt][ncent];
TH1F* hforebbcw[npt][ncent];
TH1F* kbackbbcw2[npt][ncent];
TH1F* kforebbcw_centIn[npt];
TH1F* hforebbcw_centIn[npt];
TH1F* kbackbbcw2_centIn[npt];
double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
double centbin[ncent+1] = {0,5,10,60,88};
/*double centmin = 1.0; double centmax = 3.0;//has to be boundary
for(int icent=0; icent<ncent; icent++){
if(centmin >= centbin[icent] && centmin < centbin[icent+1]){int xcentmin = icent; continue;}
if(centmax >= centbin[icent] && centmax < centbin[icent+1]){int xcentmax = icent; continue;}
}*/
for(int icent_a=0;icent_a<1;icent_a++){
int xcentmin = icent_a*1;
int xcentmax = (icent_a+1)*1;
double centmin = centbin[xcentmin];
double centmax = centbin[xcentmax];

for(int ipt=0; ipt<npt; ipt++){
for(int icent=0; icent<ncent; icent++){
kforebbcw[ipt][icent] = (TH1F*)f->Get(Form("kfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
hforebbcw[ipt][icent] = (TH1F*)f->Get(Form("hfore%sbbcw_%d_%d",dire.Data(),icent,ipt));
kbackbbcw2[ipt][icent] = (TH1F*)f->Get(Form("kback%sbbcw2_%d_%d",dire.Data(),icent,ipt));
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
  c1=new TCanvas("c1","c1", 600,750);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
for(int ipt=0; ipt<6; ipt++){
	if(ipt==4){
for(int ipt_add=5;ipt_add<6;ipt_add++){
kforebbcw_centIn[ipt]->Add(kforebbcw_centIn[ipt_add]);
hforebbcw_centIn[ipt]->Add(hforebbcw_centIn[ipt_add]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2_centIn[ipt_add]);
}
ptbin[ipt+1] = 3.0 ;
}
	if(ipt==5){
for(int ipt_add=6;ipt_add<10;ipt_add++){
kforebbcw_centIn[ipt]->Add(kforebbcw_centIn[ipt_add]);
hforebbcw_centIn[ipt]->Add(hforebbcw_centIn[ipt_add]);
kbackbbcw2_centIn[ipt]->Add(kbackbbcw2_centIn[ipt_add]);
}
ptbin[ipt+1] = 5.0;
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

  float ytop=0.97, ybot=0.07;
  float bot_r = 0.082, top_r = 0.042;
  float ny = 1.0/(1-bot_r-0.02)+1.0/(1-top_r-0.02)+1.0/(1.0-0.04);
  float ygap = (ytop-ybot)/ny;
  float ygap1 = ytop - 1.0/(1-top_r-0.02)*ygap;
  float ygap2 = ytop - 1.0/(1-top_r-0.02)*ygap - 1.0/(1.0-0.04)*ygap;
  float xleft=0.07, xright=0.98, xmiddle=1.0/(0.84/0.85+1.0)*(xright-xleft)+xleft;
  c1->cd();
  if(ipt==0)TPad *c1_ = new TPad(Form("c1_%d",ipt), Form("c1_%d",ipt),xleft,ygap1,xmiddle,ytop);
  if(ipt==1)TPad *c1_ = new TPad(Form("c1_%d",ipt), Form("c1_%d",ipt),xmiddle,ygap1,xright,ytop);
  if(ipt==2)TPad *c1_ = new TPad(Form("c1_%d",ipt), Form("c1_%d",ipt),xleft,ygap2,xmiddle,ygap1);
  if(ipt==3)TPad *c1_ = new TPad(Form("c1_%d",ipt), Form("c1_%d",ipt),xmiddle,ygap2,xright,ygap1);
  if(ipt==4)TPad *c1_ = new TPad(Form("c1_%d",ipt), Form("c1_%d",ipt),xleft,ybot,xmiddle,ygap2);
  if(ipt==5)TPad *c1_ = new TPad(Form("c1_%d",ipt), Form("c1_%d",ipt),xmiddle,ybot,xright,ygap2);
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
  gPad->SetBottomMargin(0.02);
  //gPad->SetTopMargin(top_r);
  gPad->SetTicks(1,1);
  

  float ymax = 1.02;//hpp[ipt]->GetMaximum()*1.1;
  float ymin = 0.98;//hpp[ipt]->GetMinimum()*0.9;

  hpp[ipt]->SetMinimum(ymin);
  hpp[ipt]->SetMaximum(ymax);

  hpp[ipt]->SetMarkerStyle(20);
  hpp[ipt]->SetMarkerSize(1.3);
  
  hpp[ipt]->SetMarkerColor(4);
  hpp[ipt]->GetYaxis()->SetLabelSize(0.08);
  hpp[ipt]->GetXaxis()->SetLabelSize(0.04);
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

TLegend *leg1 = new TLegend(0.22,0.02,0.82,0.13);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.08);
  leg1->SetBorderSize(0);
  leg1->AddEntry(fun0,"1+#Sigma2c_{n}cos(n#Delta#phi)","L");
  if(ipt==1)leg1->Draw();

  TLegend *leg2 = new TLegend(0.48,0.12,0.82,0.32);
  leg2->SetFillColor(10);
  leg2->SetLineStyle(4000);
  leg2->SetLineColor(10);
  leg2->SetLineWidth(0.);
  leg2->SetTextSize(0.08);
  leg2->SetBorderSize(0);
  leg2->AddEntry(fun1,"1+2c_{1}cos(#Delta#phi)","L");
  leg2->AddEntry(fun2,"1+2c_{2}cos(2#Delta#phi)","L");
  leg2->AddEntry(fun3,"1+2c_{3}cos(3#Delta#phi)","L");
  leg2->AddEntry(fun4,"1+2c_{4}cos(4#Delta#phi)","L");
  if(ipt==1)leg2->Draw();

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

  TLatex *t=new TLatex(-1.0,0.88*(ymax-ymin)+ymin, Form("%d)",ipt+1));
  t->SetTextSize(0.08);
  t->Draw();

  TMarker *m0 = new TMarker(-0.3, 0.91*(ymax-ymin)+ymin, 20);
  m0->SetMarkerSize(1.3);
  m0->SetMarkerColor(4);
  m0->Draw();

  TLatex *t=new TLatex(-0.0,0.88*(ymax-ymin)+ymin, Form("Run16 d+Au 200GeV"));
  t->SetTextSize(0.08);
  //if(ipt==0)
    t->Draw();
  TLatex *t = new TLatex(0,0.78*(ymax-ymin)+ymin,Form("%.1f < p_{T}< %.1f (GeV/c)",ptbin[ipt],ptbin[ipt+1]));
  t->SetTextSize(0.07);
  t->Draw();

  TLatex *t=new TLatex(-1.2,0.78*(ymax-ymin)+ymin, Form("0-5%%",centmin,centmax));
  t->SetTextSize(0.08);
  if(ipt==0)t->Draw();

  TLatex *t=new TLatex(-1.0,0.66*(ymax-ymin)+ymin, "|#eta_{trig}|<0.35");
  t->SetTextSize(0.08);
  if(ipt==0)t->Draw();

if(dire=="north")
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "3.0<#eta_{asso}<3.9");
else
  TLatex *t=new TLatex(-1.2,0.18*(ymax-ymin)+ymin, "-3.9<#eta_{asso}<-3.1");
  t->SetTextSize(0.08);
  t->SetTextColor(2);
  if(ipt==0)t->Draw();

}
  c1->cd();

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

  c1->Print(Form("fig/ridge_6p_pt_%s_%.1f_%.1f.png",dire.Data(),centmin,centmax));
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
