
void PlotV2()
{
  int i,j;
  const int Nbin=10;
  double pt[Nbin];
  double v2[Nbin];
  double v24[Nbin];
  double err[Nbin];
  double err4[Nbin];
  double tmp;
  TF1 *V2vsPt = new TF1("V2vsPt","((x/3.4)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);

  ifstream inR;
  inR.open("V2.txt");
  if(!inR.good())    cout<<" input fail"<<endl;   else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR>>pt[i]; 
    inR>>v2[i];   inR>>err[i];
    inR>>v24[i];  inR>>err4[i];
    cout<<pt[i]<<"  "<<v2[i]<<"  "<<err[i]<<"  "<<v24[i]<<"  "<<err4[i]<<endl;
  } inR.close();

  TCanvas *c1 = new TCanvas("c1","c1",1,1,700,450);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderSize(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.12);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->Divide(2,1,0,0);
  gStyle->SetOptStat(0);
  c1->SetTicks(-1);

 c1->cd(1);
 gPad->SetTicks(-1);
 TH1D* hist = new TH1D("hist","",100,0.,9.09);
 hist->SetXTitle("p_{T} (GeV/c)");
 hist->SetYTitle("v_{2}{Q-Cumulant}");
 hist->SetMinimum(-0.01);
 hist->SetMaximum(0.299);
 hist->GetXaxis()->CenterTitle(1);
 hist->GetYaxis()->CenterTitle(1);
 hist->GetYaxis()->SetTitleOffset(1.4);
 hist->GetXaxis()->SetTitleOffset(0.95);
 hist->GetXaxis()->SetTitleSize(0.056);
 hist->GetYaxis()->SetTitleSize(0.056);
 hist->GetXaxis()->SetLabelSize(0.05);
 hist->GetYaxis()->SetLabelSize(0.05);
 hist->Draw();

  TGraphErrors *graph = new TGraphErrors(Nbin,pt,v2,0,err);
  graph->SetTitle("");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(1);
  graph->SetLineColor(1);
  graph->SetLineStyle(2);
  graph->SetMarkerSize(1.);
  graph->Draw("psameez");
  V2vsPt->SetLineColor(2);
  V2vsPt->SetLineStyle(2);
  V2vsPt->SetLineWidth(1.5);
  V2vsPt->Draw("same");

    TLegend *leg0 = new TLegend(0.22,0.82,0.5,0.95);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    //leg0->SetTextColor(2);
    leg0->SetTextSize(0.051);
    //leg0->SetLineStyle(0.06);
    leg0->AddEntry(V2vsPt,"Input v_{2},  STEG,  10k","l");
    leg0->AddEntry(graph,"v_{2}{2}","p");
    leg0->Draw();


 c1->cd(2);
 gPad->SetRightMargin(0.02);
 gPad->SetTicks(-1);
 TH1D* hist1 = new TH1D("hist1","",100,0.0001,9.09);
 hist1->SetXTitle("p_{T} (GeV/c)");
 hist1->SetYTitle("v_{2}{Q-Cumulant}");
 hist1->SetMinimum(-0.01);
 hist1->SetMaximum(0.299);
 hist1->GetXaxis()->CenterTitle(1);
 hist1->GetYaxis()->CenterTitle(1);
 hist1->GetYaxis()->SetTitleOffset(1.1);
 hist1->GetXaxis()->SetTitleOffset(0.85);
 hist1->GetXaxis()->SetTitleSize(0.066);
 hist1->GetYaxis()->SetTitleSize(0.066);
 hist1->GetXaxis()->SetLabelSize(0.05);
 hist1->GetYaxis()->SetLabelSize(0.05);
 hist1->Draw();

  TGraphErrors *graph2 = new TGraphErrors(Nbin,pt,v24,0,err4);
  graph2->SetTitle("");
  graph2->SetMarkerStyle(21);
  graph2->SetMarkerColor(2);
  graph2->SetLineColor(2);
  graph2->SetLineStyle(2);
  graph2->SetMarkerSize(0.8);
  graph2->Draw("Pezsame");
  V2vsPt->SetLineColor(2);
  V2vsPt->SetLineStyle(2);
  V2vsPt->SetLineWidth(1.5);
  V2vsPt->Draw("same");

    TLegend *leg0 = new TLegend(0.1,0.82,0.42,0.95);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    //leg0->SetTextColor(2);
    leg0->SetTextSize(0.063);
    //leg0->SetLineStyle(0.06);
    leg0->AddEntry(V2vsPt,"Input v_{2},  STEG,  10k","l");
    leg0->AddEntry(graph2,"v_{2}{4}","p");
    leg0->Draw();



  c1->Print("qcValidationSTEG.gif");
  //c1->Print("qcValidationSTEG.pdf");
  //c1->Print("qcValidationSTEG.C");


}

