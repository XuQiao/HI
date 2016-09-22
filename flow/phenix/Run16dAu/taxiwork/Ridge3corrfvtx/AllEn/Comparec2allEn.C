void Comparec2allEn(int n=2){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
//  TString dire="south";
  TString dire="north";
  ifstream finp[4];
  finp[0].open(Form("../200GeV/c1_c2_central_ptccentc_%s.dat",dire.Data()));//pp x Et
  finp[1].open(Form("../62GeV/c1_c2_central_ptccentc_%s.dat",dire.Data()));//pp x Et
  finp[2].open(Form("../39GeV/c1_c2_central_ptccentc_%s.dat",dire.Data()));//pp x Et
  finp[3].open(Form("../20GeV/c1_c2_central_ptccentc_%s.dat",dire.Data()));//pp x Et
 //int n = 3;
  float pt[24], c2_p[24], ec2_p[24], ept[24], spt[24];
  float c2_pt0[24], ec2_pt0[24], spt0[24];

  float tmp;
  
  float ratio[24], eratio[24], sratio[24];
  TGraphErrors * grp[4];
  TGraphErrors * grpr[4];

for(int ien=0;ien<4;ien++){
  for(int ipt=0; ipt<24; ipt++){

    ept[ipt]=0;
    spt[ipt]=0.05;
    pt[ipt]=0.1+0.2*ipt;
    if(n==1){
   finp[ien]>>c2_p[ipt]>>ec2_p[ipt]>>tmp>>tmp>>tmp>>tmp;
   }
    else if(n==2){
    finp[ien]>>tmp>>tmp>>c2_p[ipt]>>ec2_p[ipt]>>tmp>>tmp;
    }
    else if(n==3){
    finp[ien]>>tmp>>tmp>>tmp>>tmp>>c2_p[ipt]>>ec2_p[ipt];
    }
    if(ien==0){
    c2_pt0[ipt]=c2_p[ipt];
    ec2_pt0[ipt]=ec2_p[ipt];
    }
    ratio[ipt] = c2_p[ipt]/c2_pt0[ipt];
    eratio[ipt] = c2_p[ipt]/c2_pt0[ipt] * sqrt(TMath::Power(ec2_pt0[ipt]/c2_pt0[ipt],2)+TMath::Power(ec2_p[ipt]/c2_p[ipt],2));
  }
  finp[ien].close();

  grp[ien] = new TGraphErrors(24, pt, c2_p, ept, ec2_p);
  grpr[ien] = new TGraphErrors(24, pt, ratio, ept, eratio);
}

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
  if(dire=="south"){
  if(n==2){
  h->SetMaximum(0.005);
  h->SetMinimum(0.000);
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
  h->SetMinimum(0.000);
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
for(int ien=0;ien<4;ien++){
  grp[ien]->SetMarkerStyle(20+ien);
  grp[ien]->SetMarkerSize(1.4);
  grp[ien]->SetMarkerColor(ien+1);
  grp[ien]->SetLineColor(ien+1);
  grp[ien]->Draw("PE");
}
//  TF1 *fpol2= new TF1("fpol2","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0,5);
//  grp->Fit("fpol2");

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
  //leg1->AddEntry(grc,"PHENIX","");
for(int ien=0;ien<4;ien++){
switch(ien){
    case 0: int en=200; break;
    case 1: int en=62; break;
    case 2: int en=39; break;
    case 3: int en=20; break;
}
  leg1->AddEntry(grp[ien],Form("c_{%d} dAu %d GeV , cent: 0-5%%",n,en),"P");
}
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
  hr->GetYaxis()->SetTitle("Ratio = lower/200GeV");
  //hr->GetYaxis()->SetTitle("(c_{2}^{pp}#times(#SigmaE_{T}^{pp})/(c_{2}^{dAu}#times#SigmaE_{T}^{dAu}))");
  hr->GetXaxis()->CenterTitle();
  hr->GetXaxis()->SetTitleSize(0.11);
  hr->GetXaxis()->SetTitleOffset(1.2);
  hr->GetXaxis()->SetTitle("p_{T}(GeV/c)");

for(int ien=0;ien<4;ien++){
  grpr[ien]->SetMarkerStyle(20+ien);
  grpr[ien]->SetMarkerSize(1.4);
  grpr[ien]->SetMarkerColor(ien+1);
  grpr[ien]->SetLineColor(ien+1);
  grpr[ien]->Draw("PE");
}

  TLatex *t2=new TLatex(4.5,0.2, "(b)");
  t2->SetTextSize(0.10);
//  t2->Draw();
  c1->Print(Form("Compc%dallEn_%s.png",n,dire.Data()));
}
