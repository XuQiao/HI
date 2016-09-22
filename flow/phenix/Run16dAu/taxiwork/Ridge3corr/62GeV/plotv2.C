void plotv2(){
    const int ncent = 6;
    const int npt = 25;
    int scale = 1;
    int centbin[ncent+1] = {0,5,10,20,40,60,100};
    gStyle->SetOptFit(kFALSE);
    gStyle->SetOptStat(kFALSE);

    //-----------------TGraphs--------------------------
//TGraphErrors* gr1[ncent];
//TGraphErrors* gr2[ncent];
TGraphErrors* gr[ncent];
//TGraphErrors* grsub1[ncent];
//TGraphErrors* grsub2[ncent];
TGraphErrors* grsub[ncent];
TGraphErrors* grF[ncent];
    
TGraphErrors *grP[ncent];

    for(int icent=0;icent<ncent;icent++){
    
    gr[icent] = new TGraphErrors(Form("v2_cent%d.dat",icent),"%lg %lg %lg");
    grsub[icent] = new TGraphErrors(Form("v2_cent%d_scale%d.dat",icent,scale));
    grF[icent] = new TGraphErrors(Form("../../Ridge3corrfvtx/62GeV/v2_cent%d.dat",icent),"%lg %lg %lg");

//-----------plotting results-----------------------------------------------------
    if(icent==0 || icent==1){
//    grP[icent] = new TGraphErrors(Form("v2_00_%d_BBCS.dat",icent),"%lg %lg %lg");
    grP[icent] = new TGraphErrors(Form("v2_00_%d_FVTX1S.dat",icent),"%lg %lg %lg");
    grP[icent]->SetName(Form("grP_%d",icent));
    grP[icent]->SetMarkerSize(1.2);
    grP[icent]->SetMarkerStyle(20);
    grP[icent]->SetMarkerColor(1);
    grP[icent]->SetLineColor(1);
 }
    TH1D* h = new TH1D(Form("h_%d",icent),"",100,0,5);
    h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h->GetYaxis()->SetTitle("v_{2}");
    h->GetXaxis()->SetRangeUser(0,3.5);
    h->GetYaxis()->SetRangeUser(0,0.25);

 float ptv2[npt];
  float eptv2[npt], sptv2[npt];
  float v2data[npt], ev2data[npt], sv2data[npt], asv2data1[npt];
  int npoints = 11;
  ifstream finv2("v2_pt_dAu_00_05_sys.dat");

  //Elliptic flow
  for(int i=0; i<npoints; i++)
    {
      eptv2[i]=0;
      sptv2[i]=0.05;
      finv2>>ptv2[i]>>v2data[i]>>ev2data[i]>>sv2data[i];//>>asv2data1[i];
    }
  finv2.close();

  TGraphErrors *grP0 = new TGraphErrors(npoints, ptv2, v2data, eptv2, ev2data);
  grP0->SetMarkerStyle(23);
  grP0->SetMarkerColor(2);
  TF1 *fun = new TF1("fun","pol4",0,5);
 
  /*
  grP[icent]->Fit("fun","Q0");
    TCanvas *c1 = new TCanvas("c1","",500,500);
    h->Draw();
    TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    if(icent==0){
    leg->AddEntry(grP,"Event plane","P");
    leg->AddEntry(fun,"pol 4 fit","L");
    }
    leg->AddEntry(gr1[icent],"use cnt-bbc","P");
    leg->AddEntry(gr2[icent],"use cnt-fvtx","P");
    leg->AddEntry(gr[icent],"use 3-sub","P");
    leg->Draw("same");
    grP[icent]->SetMarkerStyle(20);
    grP[icent]->SetMarkerColor(1);
    gr1[icent]->SetMarkerSize(1.2);
    gr1[icent]->SetMarkerStyle(21);
    gr1[icent]->SetMarkerColor(4);
    gr1[icent]->SetLineColor(4);
    gr2[icent]->SetMarkerSize(1.2);
    gr2[icent]->SetMarkerStyle(20);
    gr2[icent]->SetMarkerColor(2);
    gr2[icent]->SetLineColor(2);
    gr[icent]->SetMarkerSize(1.2);
    gr[icent]->SetMarkerStyle(24);
    gr[icent]->SetMarkerColor(1);
    gr[icent]->SetLineColor(1);
    gr1[icent]->Draw("Psame");
    gr2[icent]->Draw("Psame");
    gr[icent]->Draw("Psame");

//  0-5%
  if(icent==0){
 
  fun->Draw("same");
      for(int i=0; i<npoints; i++)
	{
	  double px1 = ptv2[i]-0.05;
	  double py1 = v2data[i]+sv2data[i];
	  double px2 = ptv2[i]+0.05;
	  double py2 = v2data[i]-sv2data[i];
	  TBox *boxv2 = new TBox(px1,py1,px2,py2);
	  boxv2->SetFillColor(kBlack);
	  boxv2->SetFillStyle(0);
	  boxv2->SetLineColor(kBlack);
	  boxv2->SetLineWidth(1);
	  boxv2->Draw("LSAME");
    }
 
    grP[icent]->Draw("Psame");
  }
    c1->Print(Form("v22pc_cent%d_scale%d.png",icent,scale));

    TCanvas *c2 = new TCanvas("c2","",500,500);
    h->Draw();
    TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    //leg->AddEntry(grP,"Preliminary","P");
    //leg->AddEntry(fun,"pol 4 fit","L");
    leg->AddEntry(gr1[icent],"use cnt-bbc","P");
    leg->AddEntry(grsub1[icent],"use cnt-bbc subtracted","P");
    //leg->AddEntry(gr,"use 3-sub","P");
    leg->Draw();
    gr1[icent]->SetMarkerSize(1.2);
    gr1[icent]->SetMarkerStyle(21);
    gr1[icent]->SetMarkerColor(4);
    gr1[icent]->SetLineColor(4);
    grsub1[icent]->SetMarkerSize(1.2);
    grsub1[icent]->SetMarkerStyle(21);
    grsub1[icent]->SetMarkerColor(2);
    grsub1[icent]->SetLineColor(2);
    gr1[icent]->Draw("Psame");
    grsub1[icent]->Draw("Psame");
  if(icent==0){
  fun->Draw("same");
      for(int i=0; i<npoints; i++)
	{
	  double px1 = ptv2[i]-0.05;
	  double py1 = v2data[i]+sv2data[i];
	  double px2 = ptv2[i]+0.05;
	  double py2 = v2data[i]-sv2data[i];
	  TBox *boxv2 = new TBox(px1,py1,px2,py2);
	  boxv2->SetFillColor(kBlack);
	  boxv2->SetFillStyle(0);
	  boxv2->SetLineColor(kBlack);
	  boxv2->SetLineWidth(1);
	  boxv2->Draw("LSAME");
    }
    grP->Draw("Psame");
  }
    c2->Print(Form("v2cntbbc_cent%d_scale%d.png",icent,scale));

    TCanvas *c3 = new TCanvas("c3","",500,500);
    h->Draw();
    TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(gr1[icent],"use cnt-fvtx","P");
    leg->AddEntry(grsub1[icent],"use cnt-fvtx subtracted","P");
    leg->Draw();
    gr2[icent]->SetMarkerSize(1.2);
    gr2[icent]->SetMarkerStyle(21);
    gr2[icent]->SetMarkerColor(4);
    gr2[icent]->SetLineColor(4);
    grsub2[icent]->SetMarkerSize(1.2);
    grsub2[icent]->SetMarkerStyle(21);
    grsub2[icent]->SetMarkerColor(2);
    grsub2[icent]->SetLineColor(2);
    gr2[icent]->Draw("Psame");
    grsub2[icent]->Draw("Psame");
  if(icent==0){
  fun->Draw("same");
      for(int i=0; i<npoints; i++)
	{
	  double px1 = ptv2[i]-0.05;
	  double py1 = v2data[i]+sv2data[i];
	  double px2 = ptv2[i]+0.05;
	  double py2 = v2data[i]-sv2data[i];
	  TBox *boxv2 = new TBox(px1,py1,px2,py2);
	  boxv2->SetFillColor(kBlack);
	  boxv2->SetFillStyle(0);
	  boxv2->SetLineColor(kBlack);
	  boxv2->SetLineWidth(1);
	  boxv2->Draw("LSAME");
    }
    grP->Draw("Psame");
  }
    c3->Print(Form("v2cntfvtx_cent%d_scale%d.png",icent,scale));
 */   
    TCanvas *c4 = new TCanvas(Form("c4_%d",icent),"",500,500);
    h->Draw();
    TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(gr[icent],"use 3-sub cnt-fvtxs-bbcs","P");
    leg->AddEntry(grF[icent],"use 3-sub cnt-fvtxs-fvtxn","P");
    if(icent==0 || icent==1)
    leg->AddEntry(grP[icent],"use FVTXS EP","P");
    //leg->AddEntry(grsub1[icent],"use 3-sub subtracted","P");
    leg->Draw("same");
    gr[icent]->SetMarkerSize(1.2);
    gr[icent]->SetMarkerStyle(21);
    gr[icent]->SetMarkerColor(4);
    gr[icent]->SetLineColor(4);
   // grsub[icent]->SetMarkerSize(1.2);
   // grsub[icent]->SetMarkerStyle(21);
   // grsub[icent]->SetMarkerColor(2);
   // grsub[icent]->SetLineColor(2);
    grF[icent]->SetMarkerSize(1.2);
    grF[icent]->SetMarkerStyle(21);
    grF[icent]->SetMarkerColor(2);
    grF[icent]->SetLineColor(2);
    gr[icent]->Draw("Psame");
    //grsub[icent]->Draw("Psame");
    grF[icent]->Draw("Psame");
    if(icent==0 || icent==1){
      /*
  fun->Draw("same");
      for(int i=0; i<npoints; i++)
	{
	  double px1 = ptv2[i]-0.05;
	  double py1 = v2data[i]+sv2data[i];
	  double px2 = ptv2[i]+0.05;
	  double py2 = v2data[i]-sv2data[i];
	  TBox *boxv2 = new TBox(px1,py1,px2,py2);
	  boxv2->SetFillColor(kBlack);
	  boxv2->SetFillStyle(0);
	  boxv2->SetLineColor(kBlack);
	  boxv2->SetLineWidth(1);
	  boxv2->Draw("LSAME");
    }
    */
   // grP0->Draw("Psame");
    grP[icent]->Draw("Psame");
  }
    c4->Print(Form("v23sub_cent%d_scale%d.png",icent,scale));
/*
    TCanvas *c5 = new TCanvas("c5","",500,500);
    h->Draw();
    TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(grsub1[icent],"use cnt-bbc subtracted","P");
    leg->AddEntry(grsub2[icent],"use cnt-fvtx subtracted","P");
    leg->AddEntry(grsub[icent],"use 3-sub subtracted","P");
    leg->Draw();
    grsub1[icent]->SetMarkerSize(1.2);
    grsub1[icent]->SetMarkerStyle(21);
    grsub1[icent]->SetMarkerColor(4);
    grsub1[icent]->SetLineColor(4);
    grsub2[icent]->SetMarkerSize(1.2);
    grsub2[icent]->SetMarkerStyle(20);
    grsub2[icent]->SetMarkerColor(2);
    grsub2[icent]->SetLineColor(2);
    grsub[icent]->SetMarkerSize(1.2);
    grsub[icent]->SetMarkerStyle(24);
    grsub[icent]->SetMarkerColor(1);
    grsub[icent]->SetLineColor(1);
    grsub1[icent]->Draw("Psame");
    grsub2[icent]->Draw("Psame");
    grsub[icent]->Draw("Psame");

  if(icent==0){
  fun->Draw("same");
      for(int i=0; i<npoints; i++)
	{
	  double px1 = ptv2[i]-0.05;
	  double py1 = v2data[i]+sv2data[i];
	  double px2 = ptv2[i]+0.05;
	  double py2 = v2data[i]-sv2data[i];
	  TBox *boxv2 = new TBox(px1,py1,px2,py2);
	  boxv2->SetFillColor(kBlack);
	  boxv2->SetFillStyle(0);
	  boxv2->SetLineColor(kBlack);
	  boxv2->SetLineWidth(1);
	  boxv2->Draw("LSAME");
    }
    grP->Draw("Psame");
  }
    c5->Print(Form("v22pc_subtract_cent%d_scale%d.png",icent,scale));
    */
}
}


