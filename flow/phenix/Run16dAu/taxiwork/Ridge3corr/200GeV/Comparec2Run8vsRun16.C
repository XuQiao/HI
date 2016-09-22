void Comparec2Run8vsRun16(int n = 2){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);

  ifstream finc1("c2_pt_south_dAu_00_05_sys.dat");//dAu 0-5%
  ifstream finc("c1_c2_pt_south_bbc_00_05_dAu_ppg161.dat");
  ifstream finp("c1_c2_central_ptccentc_south.dat");//pp x Et
  float pt[25], ept[25], spt[25];
  float ptp[25];

  float c2_c[25], ec2_c[25], sc2_c[25];
  float c2_p[25], ec2_p[25], sc2_p[25];
  float c1_c[25], ec1_c[25], sc1_c[25];
  float c1_p[25], ec1_p[25], sc1_p[25];
  float tmp;
  
  float ratio[25], eratio[25], sratio[25];

  for(int ipt=0; ipt<25; ipt++){

    ept[ipt]=0;
    spt[ipt]=0.05;
    ptp[ipt]=0.1+0.2*ipt;
    finc1>>pt[ipt]>>c2_c[ipt]>>ec2_c[ipt]>>sc2_c[ipt];
    if(n==0){
   finc>>c1_c[ipt]>>ec1_c[ipt]>>c2_c[ipt]>>ec2_c[ipt]>>tmp>>tmp;
   finp>>c1_p[ipt]>>ec1_p[ipt]>>c2_p[ipt]>>ec2_p[ipt]>>tmp>>tmp;
   if(c1_c[ipt]!=0){
   c2_c[ipt] = -c2_c[ipt]/c1_c[ipt];
   ec2_c[ipt] = c2_c[ipt]*sqrt(TMath::Power(ec2_c[ipt]/c2_c[ipt],2)+TMath::Power(ec1_c[ipt]/c1_c[ipt],2));
   }
   c2_p[ipt] = -c2_p[ipt]/c1_p[ipt];
   ec2_p[ipt] = c2_p[ipt]*sqrt(TMath::Power(ec2_p[ipt]/c2_p[ipt],2)+TMath::Power(ec1_p[ipt]/c1_p[ipt],2));
   }
    if(n==1){
   finc>>c2_c[ipt]>>ec2_c[ipt]>>tmp>>tmp>>tmp>>tmp;
   finp>>c2_p[ipt]>>ec2_p[ipt]>>tmp>>tmp>>tmp>>tmp;
   }
    if(n==2){
    finc>>tmp>>tmp>>c2_c[ipt]>>ec2_c[ipt]>>tmp>>tmp;
    finp>>tmp>>tmp>>c2_p[ipt]>>ec2_p[ipt]>>tmp>>tmp;
    }
    if(n==3){
    finc>>tmp>>tmp>>tmp>>tmp>>c2_c[ipt]>>ec2_c[ipt];
    finp>>tmp>>tmp>>tmp>>tmp>>c2_p[ipt]>>ec2_p[ipt];
    }
    ratio[ipt] = c2_p[ipt]/c2_c[ipt];
    eratio[ipt] = c2_p[ipt]/c2_c[ipt] * sqrt(TMath::Power(ec2_c[ipt]/c2_c[ipt],2)+TMath::Power(ec2_p[ipt]/c2_p[ipt],2));
    //finr>>pt[ipt]>>ratio[ipt]>>eratio[ipt]>>sratio[ipt];

    cout<<c2_c[ipt]<<endl;
  }
  finc.close();
  finp.close();
  //finr.close();

  TGraphErrors * grc = new TGraphErrors(9, pt, c2_c, ept, ec2_c);
  TGraphErrors * grp = new TGraphErrors(25, ptp, c2_p, ept, ec2_p);
  TGraphErrors * grr = new TGraphErrors(9, pt, ratio, ept, eratio);

  TGraphErrors * sgrc = new TGraphErrors(9, pt, c2_c, spt, sc2_c);
  TGraphErrors * sgrp = new TGraphErrors(9, pt, c2_p, spt, sc2_p);
  TGraphErrors * sgrr = new TGraphErrors(9, pt, ratio, spt, sratio);

  c1=new TCanvas("c1","c1", 600,500);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);

  c1->cd();
  TPad *c1_1 = new TPad("c1_1", "c1_1",0.03,0.43,0.97,0.97);
  //TPad *c1_1 = new TPad("c1_1", "c1_1",0.03,ygap1,0.97,ytop);
  c1_1->SetFillColor(10);
  c1_1->SetFillColor(10);
  c1_1->SetBorderMode(0);
  c1_1->SetBorderSize(2);
  c1_1->SetFrameFillColor(0);
  c1_1->SetFrameBorderMode(0);
  c1_1->SetFrameBorderMode(0);
  c1_1->Draw();
  c1_1->cd();
 // c1_1->SetLogy(1);

  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.10);
  gPad->SetBottomMargin(0);
  gPad->SetTopMargin(0.03);
  gPad->SetTicks(1,1);

  TH1F *h = new TH1F("h", "h", 50, 0, 5.0);
  h->GetXaxis()->SetRangeUser(0.5,3.0);
  if(n==0){
  h->SetMaximum(0.8);
  h->SetMinimum(0.2);
  }
  if(n==2){
  h->SetMaximum(0.005);
  h->SetMinimum(0.000);
  }
  if(n==1){
  h->SetMaximum(0.0001);
  h->SetMinimum(-0.01);
  }
  if(n==3){
  h->SetMaximum(0.001);
  h->SetMinimum(-0.001);
  }

  h->GetYaxis()->SetLabelSize(0.07);
  h->GetXaxis()->SetLabelSize(0.00);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleSize(0.09);
  h->GetYaxis()->SetTitleOffset(0.7);
  if(n==0)
  h->GetYaxis()->SetTitle(Form("-c_{2}/c_{1}"));
  else
  h->GetYaxis()->SetTitle(Form("c_{%d}",n));
  h->Draw();

  sgrc->SetMarkerStyle(20);
  sgrc->SetMarkerColor(2);
  sgrc->SetMarkerSize(1.4);
  sgrc->SetFillColor(17);
  //sgrc->Draw("PE2");

  grc->SetMarkerStyle(20);
  grc->SetMarkerSize(1.4);
  grc->SetMarkerColor(2);
  grc->Draw("PE");
  if(n==0)
  TF1 *fpol1= new TF1("fpol1","[0]+[1]*x+[2]*x*x+[3]*x*x",0.5,3.0);
  else
  TF1 *fpol1= new TF1("fpol1","[0]+[1]*x+[2]*x*x+[3]*x*x",0.5,3.0);
  fpol1->SetLineColor(2);
  grc->Fit("fpol1","R");

  sgrp->SetMarkerStyle(22);
  sgrp->SetMarkerColor(4);
  sgrp->SetMarkerSize(1.4);
  sgrp->SetFillColor(8);
  //sgrp->Draw("PE2");

  grp->SetMarkerStyle(22);
  grp->SetMarkerSize(1.4);
  grp->SetMarkerColor(4);
  grp->Draw("PE");
  if(n==0)
  TF1 *fpol2= new TF1("fpol2","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.5,3.0);
  else
  TF1 *fpol2= new TF1("fpol2","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.5,3.0);
  fpol2->SetLineColor(4);
  grp->Fit("fpol2","R");

  TLegend *leg1 = new TLegend(0.28,0.04,0.43,0.28);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.08);
  leg1->SetBorderSize(0);
  //leg1->AddEntry(grc,"PHENIX","");
  leg1->AddEntry(grc,Form("A = c_{%d}^{dAu Run8}; cent:0-5%%",n),"P");
  leg1->AddEntry(grp,Form("B = c_{%d}^{dAu Run16}; cent: 0-5%%",n),"P");
  leg1->Draw();

  TLatex *t1=new TLatex(4.5,0.015, "(a)");
  t1->SetTextSize(0.08);
  t1->Draw();

  c1->cd();
  TPad *c1_2 = new TPad("c1_2", "c1_2",0.03,0.01,0.97,0.43);
  c1_2->SetFillColor(10);
  c1_2->SetFillColor(10);
  c1_2->SetBorderMode(0);
  c1_2->SetBorderSize(2);
  c1_2->SetFrameFillColor(0);
  c1_2->SetFrameBorderMode(0);
  c1_2->SetFrameBorderMode(0);
  c1_2->Draw();
  c1_2->cd();

  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.10);
  gPad->SetBottomMargin(0.30);
  gPad->SetTopMargin(0.00);
  gPad->SetTicks(1,1);

  TH1F *hr = new TH1F("hr", "hr", 50, 0, 5.0);
  hr->GetXaxis()->SetRangeUser(0.5,3.0);
  hr->SetMaximum(2.);
  hr->SetMinimum(0.);
  if(n==3){
  hr->SetMaximum(2);
  hr->SetMinimum(0);
  }
  hr->Draw();

  hr->GetYaxis()->SetLabelSize(0.09);
  hr->GetXaxis()->SetLabelSize(0.10);

  hr->GetYaxis()->SetNdivisions(405);
  hr->GetYaxis()->CenterTitle();
  hr->GetYaxis()->SetTitleSize(0.10);
  hr->GetYaxis()->SetTitleOffset(0.6);
  hr->GetYaxis()->SetTitle("Ratio = B/A");
  //hr->GetYaxis()->SetTitle("(c_{2}^{pp}#times(#SigmaE_{T}^{pp})/(c_{2}^{dAu}#times#SigmaE_{T}^{dAu}))");
  hr->GetXaxis()->CenterTitle();
  hr->GetXaxis()->SetTitleSize(0.11);
  hr->GetXaxis()->SetTitleOffset(1.2);
  hr->GetXaxis()->SetTitle("p_{T}(GeV/c)");

  sgrr->SetMarkerStyle(20);
  sgrr->SetMarkerColor(4);
  sgrr->SetMarkerSize(1.4);
  sgrr->SetFillColor(15);
//  sgrr->Draw("PE2");

  grr->SetMarkerStyle(20);
  grr->SetMarkerSize(1.4);
  grr->SetMarkerColor(4);
//  grr->Draw("PE");

  TF1* fr = new TF1("fr","fpol2/fpol1",0,5);
  fr->Draw("same");
  TLine *tl = new TLine(0.5,1.0,3.0,1.0);
  tl->SetLineStyle(2);
  tl->SetLineWidth(4);
  tl->Draw("same");

  TLatex *t2=new TLatex(4.5,0.2, "(b)");
  t2->SetTextSize(0.10);
//  t2->Draw();
  c1->Print(Form("Compc%dRun8vsRun16_moreptbin.png",n));
}
