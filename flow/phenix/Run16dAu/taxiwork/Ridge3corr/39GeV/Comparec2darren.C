void Comparec2darren(int n=2){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
//  TString dire="south";
  TString dire="north";
//  TString dire="sn";

  ifstream finp,finc;
  finp.open(Form("c1_c2_central_ptccentc_%s.dat",dire.Data()));//pp x Et
  finc.open(Form("c1_c2_darren_%s.dat",dire.Data()));//pp x Et
  const int npt = 5;
  float c2_p[npt], ec2_p[npt], ept[npt], spt[npt];
  float ptpp[npt], c2_c[7], ec2_c[7];
  float c2_pt0[npt], ec2_pt0[npt], spt0[npt];

  float tmp;
  
  float ratio[npt], eratio[npt], sratio[npt];
  TGraphErrors * grc;
  TGraphErrors * grpr;
  float pt[npt]={0.4,0.9,1.5,2.1,2.7};

  for(int ipt=0; ipt<npt; ipt++){
    ept[ipt]=0;
    spt[ipt]=0.05;
    if(n==1){
   finp>>c2_p[ipt]>>ec2_p[ipt]>>tmp>>tmp>>tmp>>tmp;
   }
    else if(n==2){
    finp>>tmp>>tmp>>c2_p[ipt]>>ec2_p[ipt]>>tmp>>tmp;
    }
    else if(n==3){
    finp>>tmp>>tmp>>tmp>>tmp>>c2_p[ipt]>>ec2_p[ipt];
    }
  }
   c2_p[0]=-9999;
   ec2_p[0]=-9999;
  string tmps;
  for(int ipt=0; ipt<7; ipt++){
    if(n==1){
      finc>>tmps>>tmps>>tmps>>tmps>>tmps>>c2_c[ipt]>>ec2_c[ipt]>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps;
    }
    if(n==2){
      finc>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>c2_c[ipt]>>ec2_c[ipt]>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps;
      cout<<c2_c[ipt]<<endl;
    }
    if(n==3){
      finc>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>tmps>>c2_c[ipt]>>ec2_c[ipt]>>tmps>>tmps;
    }
  }
    float ptd[7] = {0.375,0.625,0.875,1.25,1.75,2.5,4.0};
    grp = new TGraphErrors(7, ptd, c2_c, ept, ec2_c);
    TF1 *f1 = new TF1("f1","pol3",0,3.5);
    TFitResultPtr r1 = grp->Fit("f1", "S");
    double x[7];
    for(int ipt=0;ipt<npt;ipt++){
        x[ipt] = pt[ipt];
    }
    double ci1[npt],ci2[npt];
    double cl = 0.683;  // for 1 sigma error
    r1->GetConfidenceIntervals(npt,1,1,x,ci1,cl);
    for(int ipt=0;ipt<npt;ipt++){
        c2_c[ipt] = f1->Eval(x[ipt]);
        ec2_c[ipt] = ci1[ipt];
    ratio[ipt] = c2_p[ipt]/c2_c[ipt];
    eratio[ipt] = c2_p[ipt]/c2_c[ipt] * sqrt(TMath::Power(ec2_c[ipt]/c2_c[ipt],2)+TMath::Power(ec2_p[ipt]/c2_p[ipt],2));
    }
  finp.close();
  finc.close();

  grc = new TGraphErrors(npt, pt, c2_p, ept, ec2_p);
  grpr = new TGraphErrors(npt, pt, ratio, ept, eratio);

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
  h->GetXaxis()->SetRangeUser(0.,3.5);
  if(dire=="south"){
  if(n==2){
  h->SetMaximum(0.0005);
  h->SetMinimum(0.000);
  }
  else if(n==1){
  h->SetMaximum(0.0001);
  h->SetMinimum(-0.01);
  }
  else if(n==3){
  h->SetMaximum(0.001);
  h->SetMinimum(-0.001);
  }
  }
  if(dire=="north"){
  if(n==2){
  h->SetMaximum(0.015);
  h->SetMinimum(0.000);
  }
  else if(n==1){
  h->SetMaximum(0.0);
  h->SetMinimum(-0.03);
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
  h->DrawCopy();
  grp->SetMarkerStyle(20);
  grp->SetMarkerSize(1.4);
  grp->SetMarkerColor(1);
  grp->SetLineColor(1);
  grc->SetMarkerStyle(20);
  grc->SetMarkerSize(1.4);
  grc->SetMarkerColor(4);
  grc->SetLineColor(1);
  grp->Draw("PEsame");
  grc->Draw("PEsame");
  f1->Draw("same");

if(n==1)
  TLegend *leg1 = new TLegend(0.28,0.04,0.43,0.28);
if(n==2)
  TLegend *leg1 = new TLegend(0.28,0.54,0.43,0.88);
if(n==3)
  TLegend *leg1 = new TLegend(0.28,0.54,0.43,0.88);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.08);
  leg1->SetBorderSize(0);
  if(dire=="south")
  leg1->AddEntry(grc,Form("Run 16 39GeV CNT-BBCs 0-5%%"),"");
  else if(dire=="north")
  leg1->AddEntry(grc,Form("CNT-FVTXs 0-5%%"),"");
  else
  leg1->AddEntry(grc,Form("BBCs-FVTXs 0-5%%"),"");
  leg1->AddEntry(grc,Form("c_{%d} Qiao",n),"P");
  leg1->AddEntry(grp,Form("c_{%d} Darren", n),"P");
  leg1->Draw("same");

  TLatex *t1=new TLatex(4.5,0.015, "(a)");
  t1->SetTextSize(0.08);
  t1->Draw("same");

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
  hr->GetXaxis()->SetRangeUser(0,3.5);
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
  hr->GetYaxis()->SetTitle(Form("c_{%d} Qiao/ c_{%d} Darren",n,n));
  //hr->GetYaxis()->SetTitle("(c_{2}^{pp}#times(#SigmaE_{T}^{pp})/(c_{2}^{dAu}#times#SigmaE_{T}^{dAu}))");
  hr->GetXaxis()->CenterTitle();
  hr->GetXaxis()->SetTitleSize(0.11);
  hr->GetXaxis()->SetTitleOffset(1.2);
  hr->GetXaxis()->SetTitle("p_{T}(GeV/c)");

  grpr->SetMarkerStyle(20);
  grpr->SetMarkerSize(1.4);
  grpr->SetMarkerColor(1);
  grpr->SetLineColor(1);
  grpr->Draw("PE");

  TLatex *t2=new TLatex(4.5,0.2, "(b)");
  t2->SetTextSize(0.10);
//  t2->Draw();
  c1->Print(Form("Compc%dc%d_%s.png",n,n,dire.Data()));
}
