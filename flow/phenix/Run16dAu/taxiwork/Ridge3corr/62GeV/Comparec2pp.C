void Comparec2pp(int n=2){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
  TString dire="south";
//  TString dire="north";
  ifstream finp,finc,finper;
  finp.open(Form("c1_c2_central_ptccentc_%s.dat",dire.Data()));//pp x Et
  finper.open(Form("c1_c2_central_per_%s.dat",dire.Data()));//pp x Et
  finc.open(Form("c1_c2_PP_centIn_%s.dat",dire.Data()));
 const int npt = 8;
  float c2_p[npt], ec2_p[npt], ept[npt], spt[npt];
  float c2_per[npt], ec2_per[npt];
  float ptpp[npt], c2_c[npt], ec2_c[npt];
  float c2_pt0[npt], ec2_pt0[npt], spt0[npt];
if(dire=="south"){
float MPPbbc[6] = {3.8,3.8,3.8,3.8,3.8,3.8};
float Mperbbc[6] = {4.21,4.21,4.21,4.21,4.21,4.21};
float Mbbc[6] = {63.93,47.41,38.03,25.73,13.63,4.21};
}
if(dire=="north"){
float MPPbbc[6] = {0.69,0.69,0.69,0.69,0.69,0.69};
float Mbbc[6]  = {11.35,9.66,8.52,6.71,4.55,2.34};
float Mperbbc[6] = {2.34,2.34,2.34,2.34,2.34,2.34};
}
  float tmp;
  
  float ratio[npt], eratio[npt], sratio[npt];
  float ratioper[npt], eratioper[npt], sratioper[npt];
  TGraphErrors * grp;
  TGraphErrors * grc;
  TGraphErrors * grper;
  TGraphErrors * grpr;
  TGraphErrors * grperr;
  float pt[npt]={0.4,0.9,1.5,2.1,2.7};

  for(int ipt=0; ipt<npt; ipt++){

    ept[ipt]=0;
    spt[ipt]=0.05;
    ptpp[ipt]=0.25+0.5*ipt;
    if(n==1){
   finp>>c2_p[ipt]>>ec2_p[ipt]>>tmp>>tmp>>tmp>>tmp;
   }
    else if(n==2){
    finp>>tmp>>tmp>>c2_p[ipt]>>ec2_p[ipt]>>tmp>>tmp;
    finper>>tmp>>tmp>>c2_per[ipt]>>ec2_per[ipt]>>tmp>>tmp;
    finc>>tmp>>tmp>>c2_c[ipt]>>ec2_c[ipt]>>tmp>>tmp;
    }
    else if(n==3){
    finp>>tmp>>tmp>>tmp>>tmp>>c2_p[ipt]>>ec2_p[ipt];
    }
  }
    const int nptpp = 10;
    float ptppmean[nptpp] = {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75};
    TGraphErrors *grPPc2 = new TGraphErrors(nptpp,ptppmean,c2_c,0,ec2_c);
    TF1 *f1 = new TF1("f1","pol3",0,3.5);
    TFitResultPtr r1 = grPPc2->Fit("f1", "S");
    double x[npt]={0.4,0.9,1.5,2.1,2.7};
    double ci1[npt],ci2[npt];
    double cl = 0.683;  // for 1 sigma error
    r1->GetConfidenceIntervals(npt,1,1,x,ci1,cl);
    for(int ipt=0;ipt<npt;ipt++){
        c2_c[ipt] = f1->Eval(x[ipt]);
        ec2_c[ipt] = ci1[ipt];
    c2_c[ipt] *= MPPbbc[0]/Mbbc[0];
    ec2_c[ipt] *= MPPbbc[0]/Mbbc[0];
    c2_per[ipt] *= Mperbbc[0]/Mbbc[0];
    ec2_per[ipt] *= Mperbbc[0]/Mbbc[0];
    ratio[ipt] = c2_c[ipt]/c2_p[ipt];
    eratio[ipt] = c2_c[ipt]/c2_p[ipt] * sqrt(TMath::Power(ec2_c[ipt]/c2_c[ipt],2)+TMath::Power(ec2_p[ipt]/c2_p[ipt],2));
    ratioper[ipt] = c2_per[ipt]/c2_p[ipt];
    eratioper[ipt] = c2_per[ipt]/c2_p[ipt] * sqrt(TMath::Power(ec2_per[ipt]/c2_per[ipt],2)+TMath::Power(ec2_p[ipt]/c2_p[ipt],2));
    }
  finp.close();
  finper.close();
  finc.close();

  grp = new TGraphErrors(npt, pt, c2_p, ept, ec2_p);
  grper = new TGraphErrors(npt, pt, c2_per, ept, ec2_per);
  grc = new TGraphErrors(npt, pt, c2_c, ept, ec2_c);
  grpr = new TGraphErrors(npt, pt, ratio, ept, eratio);
  grperr = new TGraphErrors(npt, pt, ratioper, ept, eratioper);

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
  c1_1->SetLogy(1);

  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.10);
  gPad->SetBottomMargin(0);
  gPad->SetTopMargin(0.03);
  gPad->SetTicks(1,1);

  TH1F *h = new TH1F("h", "h", 50, 0, 5.0);
  h->GetXaxis()->SetRangeUser(0.2,3.0);
  if(dire=="south"){
  if(n==2){
  h->SetMaximum(0.006);
  h->SetMinimum(1e-6);
  }
  else if(n==1){
  h->SetMaximum(0.0001);
  h->SetMinimum(-0.02);
  }
  else if(n==3){
  h->SetMaximum(0.001);
  h->SetMinimum(-0.001);
  }
  }
  if(dire=="north"){
  if(n==2){
  h->SetMaximum(0.015);
  h->SetMinimum(1e-5);
  }
  else if(n==1){
  h->SetMaximum(0.0001);
  h->SetMinimum(-0.04);
  }
  else if(n==3){
  h->SetMaximum(0.001);
  h->SetMinimum(-0.001);
  }
  }
  h->GetYaxis()->SetLabelSize(0.07);
  h->GetXaxis()->SetLabelSize(0.00);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleSize(0.09);
  h->GetYaxis()->SetTitleOffset(0.7);
  h->GetYaxis()->SetTitle(Form("c_{%d}",n));
  h->Draw();
  grp->SetMarkerStyle(20);
  grp->SetMarkerSize(1.4);
  grp->SetMarkerColor(1);
  grp->SetLineColor(1);
  grc->SetMarkerStyle(20);
  grc->SetMarkerSize(1.4);
  grc->SetMarkerColor(4);
  grc->SetLineColor(1);
  grper->SetMarkerStyle(20);
  grper->SetMarkerSize(1.4);
  grper->SetMarkerColor(2);
  grper->SetLineColor(2);
  grp->Draw("PE");
  grper->Draw("PE");
//  grc->Draw("PE");
//  TF1 *fpol2= new TF1("fpol2","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0,5);
//  grp->Fit("fpol2");

if(n==1)
  TLegend *leg1 = new TLegend(0.28,0.04,0.43,0.28);
if(n==2)
  TLegend *leg1 = new TLegend(0.28,0.12,0.43,0.35);
if(n==3)
  TLegend *leg1 = new TLegend(0.28,0.12,0.43,0.35);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.08);
  leg1->SetBorderSize(0);
  leg1->AddEntry(grp,Form("A=c_{%d} d+Au 62 GeV, 0-5%%",n),"P");
  leg1->AddEntry(grper,Form("B=c_{%d} d+Au 62 GeV, 60-100%% * %.2f/%.2f",n, Mbbc[5], Mbbc[0]),"P");
//  leg1->AddEntry(grc,Form("C=c_{%d} minbias p+p 200GeV * %.2f/%.2f", n,MPPbbc[0],Mbbc[0]),"P");
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
  hr->GetXaxis()->SetRangeUser(0.2,3.0);
  hr->SetMaximum(0.5);
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
  hr->GetYaxis()->SetTitle("B/A");
  //hr->GetYaxis()->SetTitle("(c_{2}^{pp}#times(#SigmaE_{T}^{pp})/(c_{2}^{dAu}#times#SigmaE_{T}^{dAu}))");
  hr->GetXaxis()->CenterTitle();
  hr->GetXaxis()->SetTitleSize(0.11);
  hr->GetXaxis()->SetTitleOffset(1.2);
  hr->GetXaxis()->SetTitle("p_{T}(GeV/c)");

  grpr->SetMarkerStyle(20);
  grpr->SetMarkerSize(1.4);
  grpr->SetMarkerColor(4);
  grpr->SetLineColor(4);
  grperr->SetMarkerStyle(20);
  grperr->SetMarkerSize(1.4);
  grperr->SetMarkerColor(2);
  grperr->SetLineColor(2);
//  grpr->Draw("PE");
  grperr->Draw("PE");

  TLatex *t2=new TLatex(4.5,0.2, "(b)");
  t2->SetTextSize(0.10);
//  t2->Draw();
  c1->Print(Form("Compc%dpp_%s.png",n,dire.Data()));
}
