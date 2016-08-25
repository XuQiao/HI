void plotc2(){
  TGraphErrors *grv2 = plot("../../EPAnalyzer/pPbDataV205m100/Original",1,1,20);  
  TGraphErrors *grv21 = plot("../../EPAnalyzer/pPbDataV205m100/Varmultsper10EP10",1,6,34);  
  TGraphErrors *grv2obs = plotobs("../../EPAnalyzer/pPbDataV205m100/Original",1,1,20);  
  TGraphErrors *grv2obs1 = plotobs("../../EPAnalyzer/pPbDataV205m100/Varmultsper10EP10",1,6,34);  
  double *EPR = getEPR("../../EPAnalyzer/pPbDataV205m100/Original");
  double *EPR1 = getEPR("../../EPAnalyzer/pPbDataV205m100/Varmultsper10EP10");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
 
  TF1 *V2vsPt = new TF1("V2vsPt","((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);
  TF1 *c2vsPt = new TF1("c2vsPt","TMath::Power(((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600),2)",0.3,6.0);
  ifstream finc("Original/c1c2_centIn_south_cntbbc.dat");
  ifstream finc1("Varmultper10EP10/c1c2_centIn_south_cntbbc.dat");
  
  float pt[10], ept[10];
  float c2_c[10], ec2_c[10];
  float c2_c1[10], ec2_c1[10];
  float c2v2obs[10], ec2v2obs[10];
  float c2v2obs1[10], ec2v2obs1[10];
  float dtc2[10], edtc2[10];
  float dtv2[10], edtv2[10];
  float dtc2mdtv2[10],edtc2mdtv2[10];

  float tmp;
  for(int ipt=0; ipt<9; ipt++){

    pt[ipt] = grv2->GetX()[ipt];
    ept[ipt]=0;

    finc>>tmp>>tmp>>c2_c[ipt]>>ec2_c[ipt]>>tmp>>tmp;//>>sc2_c[ipt];
    finc1>>tmp>>tmp>>c2_c1[ipt]>>ec2_c1[ipt]>>tmp>>tmp;

    float v2_c = grv2->GetY()[ipt];
    float ev2_c = grv2->GetEY()[ipt];
    float v2_c1 = grv21->GetY()[ipt];
    float ev2_c1 = grv21->GetEY()[ipt];
    float v2obs_c = grv2obs->GetY()[ipt];
    float ev2obs_c = grv2obs->GetEY()[ipt];
    float v2obs_c1 = grv2obs1->GetY()[ipt];
    float ev2obs_c1 = grv2obs1->GetEY()[ipt];

    c2v2obs[ipt] = c2_c[ipt]/v2obs_c;
    ec2v2obs[ipt] = sqrt((ec2_c[ipt]*v2obs_c)**2+(ev2obs_c*c2_c[ipt])**2)/v2obs_c;
    c2v2obs1[ipt] = c2_c1[ipt]/v2obs_c1;
    ec2v2obs1[ipt] = sqrt((ec2_c1[ipt]*v2obs_c1)**2+(ev2obs_c1[ipt]*c2_c1)**2)/v2obs_c1;
    
    dtc2[ipt] = (c2_c1[ipt] - c2_c[ipt])/c2_c1[ipt];
    edtc2[ipt] = sqrt(((ec2_c1[ipt]*c2_c[ipt])**2 + (c2_c1[ipt]*ec2_c[ipt])**2)/c2_c[ipt]/c2_c[ipt]);
    dtv2[ipt] = (v2_c1 - v2_c)/v2_c1;
    edtv2[ipt] = sqrt(((ev2_c1*v2_c)**2 + (v2_c1*ev2_c)**2)/v2_c/v2_c);
    dtc2mdtv2[ipt] = (c2_c1[ipt] - c2_c[ipt])/c2_c1[ipt] - (v2_c1 - v2_c)/v2_c1;
//    edtc2mdtv2[ipt] = 0;
    edtc2mdtv2[ipt] = sqrt(((ec2_c1[ipt]*c2_c[ipt])**2 + (c2_c1[ipt]*ec2_c[ipt])**2)/c2_c[ipt]/c2_c[ipt] + ((ev2_c1*v2_c)**2 + (v2_c1*ev2_c)**2)/v2_c/v2_c);

  }   
   
  TLine *l[5];
  finc.close();
  finc1.close();

  TGraphErrors * grc = new TGraphErrors(9, pt, c2_c, ept, ec2_c);
  TGraphErrors * grc1 = new TGraphErrors(9, pt, c2_c1, ept, ec2_c1);
  TGraphErrors * grc2ratio = new TGraphErrors(9, pt, dtc2, ept, edtc2);
  TGraphErrors * grv2ratio = new TGraphErrors(9, pt, dtv2, ept, edtv2);
  TGraphErrors * grratio = new TGraphErrors(9, pt, dtc2mdtv2, ept, edtc2mdtv2);
  TGraphErrors * grc2v2obs = new TGraphErrors(9, pt, c2v2obs, ept, ec2v2obs);
  TGraphErrors * grc2v2obs1 = new TGraphErrors(9, pt, c2v2obs1, ept, ec2v2obs1);

  c1=new TCanvas("c1","c1", 600,600);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);

  c1->cd();
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.10);
  gPad->SetTopMargin(0.01);
  gPad->SetTicks(1,1);

  TH1F *h = new TH1F("h", "h", 80, 0, 8.0);
//  h->SetMaximum(0.032);
//  h->SetMinimum(1.2e-5);
  h->SetMaximum(0.015);
  h->SetMinimum(0);
  h->GetXaxis()->SetRangeUser(0.,6.1);
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetXaxis()->SetLabelSize(0.03);
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitle("c_{2}");
  h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h->Draw();

  grc->SetMarkerStyle(20);
  grc->SetMarkerSize(1.4);
  grc->SetMarkerColor(1);
  grc->SetLineColor(1);
  grc->Draw("PE, same");

  grc1->SetMarkerStyle(34);
  grc1->SetMarkerSize(1.6);
  grc1->SetMarkerColor(6);
  grc1->SetLineColor(6);
  grc1->Draw("PE, same");

  c2vsPt->Draw("same");

  TLegend *leg1 = new TLegend(0.60,0.72,0.87,0.90);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.04);
  leg1->SetBorderSize(0);
  leg1->AddEntry(grc,"c_{2} w/o non-flow","P");
  leg1->AddEntry(grc1,"c_{2} with non-flow","P");
  leg1->AddEntry(c2vsPt,"v_{2,input}^{2}","l");
//  leg1->AddEntry(grcmc,"HIJING c_{2}^{pAu}; cent:0-5%, <N_{bbc}>=18.1","P");
//  leg1->AddEntry(grpmc,"HIJING c_{2}^{pp} #times S, S = (<N_{bbc}^{pp}>/<N_{bbc}^{pAu}>)","P");
  leg1->Draw();

 c1->Print("c2Comp.gif"); 

  c2=new TCanvas("c2","c2", 600,600);
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);

  c2->cd();
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.10);
  gPad->SetTopMargin(0.01);
  gPad->SetTicks(1,1);

//  h->SetMaximum(0.032);
//  h->SetMinimum(1.2e-5);
  h->SetMaximum(0.15);
  h->SetMinimum(0);
  h->GetXaxis()->SetRangeUser(0.,6.1);
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetXaxis()->SetLabelSize(0.03);
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitle("v_{2}");
  h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h->Draw();

  grv2->Draw("PE, same");
  grv21->Draw("PE, same");

  V2vsPt->Draw("same");

  TLegend *leg1 = new TLegend(0.60,0.72,0.87,0.90);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.04);
  leg1->SetBorderSize(0);
  leg1->AddEntry(grv2,"v_{2} w/o non-flow","P");
  leg1->AddEntry(grv21,"v_{2} with non-flow","P");
  leg1->AddEntry(V2vsPt,"v_{2,input}","l");
//  leg1->AddEntry(grcmc,"HIJING c_{2}^{pAu}; cent:0-5%, <N_{bbc}>=18.1","P");
//  leg1->AddEntry(grpmc,"HIJING c_{2}^{pp} #times S, S = (<N_{bbc}^{pp}>/<N_{bbc}^{pAu}>)","P");
  leg1->Draw();

 c2->Print("v2Comp.gif"); 

  c3=new TCanvas("c3","c3", 600,600);
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  c3->SetBorderMode(0);
  c3->SetBorderSize(2);

  c3->cd();
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.10);
  gPad->SetTopMargin(0.01);
  gPad->SetTicks(1,1);

//  h->SetMaximum(0.032);
//  h->SetMinimum(1.2e-5);
  h->SetMaximum(0.5);
  h->SetMinimum(-0.2);
  h->GetXaxis()->SetRangeUser(0.,6.1);
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetXaxis()->SetLabelSize(0.03);
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitle("#frac{#delta c_{2}}{c_{2}} - #frac{#delta v_{2}}{v_{2}}");
  h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h->Draw();

  grratio->SetMarkerStyle(20);
  grratio->SetMarkerSize(1.4);
  grratio->SetMarkerColor(1);
  grratio->SetLineColor(1);
  grratio->Draw("PE, same");

  TLegend *leg1 = new TLegend(0.60,0.72,0.87,0.90);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.04);
  leg1->SetBorderSize(0);
//  leg1->AddEntry(grpmc,"HIJING c_{2}^{pp} #times S, S = (<N_{bbc}^{pp}>/<N_{bbc}^{pAu}>)","P");
const double etap[5]={0,1.0,1.5,2,3};
for(int ieta = 0; ieta<5;ieta++){
//    EPR[ieta] = (*vecDEPR)[ieta];
//    EPR1[ieta] = (*vecDEPR1)[ieta];
    l[ieta] = new TLine(0,(EPR1[ieta]-EPR[ieta])/EPR1[ieta],6.1,(EPR1[ieta]-EPR[ieta])/EPR1[ieta]);
    l[ieta]->SetLineColor(ieta+1);
    l[ieta]->SetLineStyle(ieta);
    l[ieta]->Draw("same");
    l[ieta]->SetLineWidth(2);
    leg1->AddEntry(l[ieta],Form("#eta gap |#eta|>%.1f",etap[ieta]));
}
  leg1->Draw();
  TLatex *t = new TLatex(1,0.3,"Lines are #frac{#delta R}{R}");
  t->Draw();

 c3->Print("deltaratio.gif"); 
  
  c4=new TCanvas("c4","c4", 600,600);
  c4->SetFillColor(10);
  c4->SetFillColor(10);
  c4->SetFillColor(10);
  c4->SetBorderMode(0);
  c4->SetBorderSize(2);

 c4->cd();
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.10);
  gPad->SetTopMargin(0.01);
  gPad->SetTicks(1,1);

//  h->SetMaximum(0.032);
//  h->SetMinimum(1.2e-5);
  h->SetMaximum(1);
  h->SetMinimum(0);
  h->GetXaxis()->SetRangeUser(0.,6.1);
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetXaxis()->SetLabelSize(0.03);
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitle("\"non-flow\" fractions");
  h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h->Draw();

  grc2ratio->SetMarkerStyle(20);
  grc2ratio->SetMarkerSize(1.4);
  grc2ratio->SetMarkerColor(1);
  grc2ratio->SetLineColor(1);
  grc2ratio->Draw("PE, same");
  grv2ratio->SetMarkerStyle(34);
  grv2ratio->SetMarkerSize(1.4);
  grv2ratio->SetMarkerColor(6);
  grv2ratio->SetLineColor(6);
  grv2ratio->Draw("PE, same");
  TLegend *leg1 = new TLegend(0.60,0.72,0.87,0.90);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.04);
  leg1->SetBorderSize(0);
  leg1->AddEntry(grc2ratio,"#frac{#delta c_{2}}{c_{2}}","p");
  leg1->AddEntry(grv2ratio,"#frac{#delta v_{2}}{v_{2}}","p");
  leg1->Draw();
  c4->Print("deltac2v2.gif"); 


  c5=new TCanvas("c5","c5", 600,600);
  c5->SetFillColor(10);
  c5->SetFillColor(10);
  c5->SetFillColor(10);
  c5->SetBorderMode(0);
  c5->SetBorderSize(2);

  c5->cd();
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.10);
  gPad->SetTopMargin(0.01);
  gPad->SetTicks(1,1);

//  h->SetMaximum(0.032);
//  h->SetMinimum(1.2e-5);
  h->SetMaximum(1);
  h->SetMinimum(0);
  h->GetXaxis()->SetRangeUser(0.,6.1);
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetXaxis()->SetLabelSize(0.03);
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitle("#frac{c_{2}}{v_{2}^{obs}}");
  h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h->Draw();

  grc2v2obs->SetMarkerStyle(20);
  grc2v2obs->SetMarkerSize(1.4);
  grc2v2obs->SetMarkerColor(1);
  grc2v2obs->SetLineColor(1);
  grc2v2obs->Draw("PE, same");
  grc2v2obs1->SetMarkerStyle(34);
  grc2v2obs1->SetMarkerSize(1.4);
  grc2v2obs1->SetMarkerColor(6);
  grc2v2obs1->SetLineColor(6);
  grc2v2obs1->Draw("PE, same");
  TLegend *leg1 = new TLegend(0.60,0.72,0.87,0.90);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.04);
  leg1->SetBorderSize(0);
  leg1->AddEntry(grc2v2obs,"w/o non-flow","P");
  leg1->AddEntry(grc2v2obs1,"with non-flow","P");
  leg1->Draw();
 c4->Print("c2v2obs.gif"); 

}

TGraphErrors* plot(TString indir, int ieta, int color,int marker){
	int ibin=0;
	TFile *f = TFile::Open(Form("%s/mergedVobs.root",indir.Data()));
	TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/E_%d/v2",ibin,ieta));
	TVectorD *vecDv2obs = (TVectorD*)f->Get(Form("D_%d/E_%d/v2obs",ibin,ieta));
	TVectorD *vecDv2sub3 = (TVectorD*)f->Get(Form("D_%d/E_%d/v23sub",ibin,ieta));
	//TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/deltavmean",ibin));
	TVectorD *vecDavgpt = (TVectorD*)f->Get(Form("D_%d/avgpt",ibin));
	double *avgpt = vecDavgpt->GetMatrixArray();
	double *v2obs = vecDv2obs->GetMatrixArray();
	double *v2 = vecDv2->GetMatrixArray();
	double *v2sub3 = vecDv2sub3->GetMatrixArray();
	//float *v2err = vecDv2err->GetMatrixArray();
	//TGraphErrors *gr=new TGraphErrors(npt,avgpt,v2,0,0);
	TGraphErrors *gr=new TGraphErrors(9,avgpt,v2,0,0);
	gr->SetTitle("v_{2} vs momentum");
	gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gr->SetMarkerSize(1.6);
	gr->SetMarkerColor(color);
	gr->SetMarkerStyle(marker);
	f->Close();
	return gr;
}
TGraphErrors* plotobs(TString indir, int ieta, int color,int marker){
	int ibin=0;
	TFile *f = TFile::Open(Form("%s/mergedVobs.root",indir.Data()));
	TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/E_%d/v2",ibin,ieta));
	TVectorD *vecDv2obs = (TVectorD*)f->Get(Form("D_%d/E_%d/v2obs",ibin,ieta));
	TVectorD *vecDv2sub3 = (TVectorD*)f->Get(Form("D_%d/E_%d/v23sub",ibin,ieta));
	//TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/deltavmean",ibin));
	TVectorD *vecDavgpt = (TVectorD*)f->Get(Form("D_%d/avgpt",ibin));
	double *avgpt = vecDavgpt->GetMatrixArray();
	double *v2obs = vecDv2obs->GetMatrixArray();
	double *v2 = vecDv2->GetMatrixArray();
	double *v2sub3 = vecDv2sub3->GetMatrixArray();
	//float *v2err = vecDv2err->GetMatrixArray();
	//TGraphErrors *gr=new TGraphErrors(npt,avgpt,v2,0,0);
	TGraphErrors *gr=new TGraphErrors(9,avgpt,v2obs,0,0);
	gr->SetTitle("v_{2} vs momentum");
	gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gr->SetMarkerSize(1.6);
	gr->SetMarkerColor(color);
	gr->SetMarkerStyle(marker);
	f->Close();
	return gr;
}

double* getEPR(TString indir){
	int ibin=0;
	TFile *f = TFile::Open(Form("%s/mergedVobs.root",indir.Data()));
	TVectorD *vecDEPR = (TVectorD*)f->Get(Form("D_%d/EPR",ibin));
	double *EPR = vecDEPR->GetMatrixArray();
	return EPR;
}


