void V2_gap0_gap1()
{
const int NCentBin=15;
const int NCentBin24=24;
  double pt[NCentBin];
  double pta[NCentBin];
  double pta2[NCentBin];
  double pta3[NCentBin];
  double v2a[NCentBin];
  double v2ae[NCentBin];
  double v21[NCentBin];
  double v21e[NCentBin];
  double v22[NCentBin24];
  double v22e[NCentBin24];
  double v23[NCentBin];
  double v23e[NCentBin];

double temp;

  ofstream sysbn;

  ifstream inv2a;
  ifstream inv21;
  ifstream inv22;
  ifstream inv23;
  inv2a.open("v2.txt");
  if(!inv2a.good())    cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(int i=0;i<NCentBin;i++){
    inv2a>>pt[i];  
    inv2a>>v2a[i];  
    inv2a>>v2ae[i];  
  } 

  for(int i=0;i<NCentBin;i++){
    inv2a>>pta[i];  //inv21>>temp;
    inv2a>>v21[i];
    inv2a>>v21e[i];
  }
  for(int i=0;i<NCentBin;i++){
    inv2a>>pta2[i];  //inv21>>temp;
    inv2a>>v22[i];
    inv2a>>v22e[i];
  }
  for(int i=0;i<NCentBin;i++){
    inv2a>>pta3[i];  //inv21>>temp;
    inv2a>>v23[i];
    inv2a>>v23e[i];
  }

////////////////

////////
  TCanvas *c1 = new TCanvas("c1","c1",1,1,650,550);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderSize(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  //c1->Divide(2,1,0,0);
 gStyle->SetOptStat(0);
  c1->SetTicks(-1);
 TH1D* hist = new TH1D("hist","",100,0.,8.9);
 hist->SetXTitle("p_{T} (GeV/c)");
 hist->SetYTitle("v_{2}{EP}");
 hist->SetMinimum(0.001);
 hist->SetMaximum(0.349);
 hist->GetXaxis()->CenterTitle(1);
 hist->GetYaxis()->CenterTitle(1);
 hist->GetYaxis()->SetTitleOffset(1.1);
 hist->GetXaxis()->SetTitleOffset(0.95);
 hist->GetXaxis()->SetTitleSize(0.066);
 hist->GetYaxis()->SetTitleSize(0.066);
 hist->GetXaxis()->SetLabelSize(0.05);
 hist->GetYaxis()->SetLabelSize(0.05);
 hist->Draw();

  TF1 *V2vsPt = new TF1("V2vsPt","((x/3.4)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);
  //TF1 *V2vsPt = new TF1("V2vsPt","0.000000001*x+0.04",0.2,10);
  V2vsPt->SetLineColor(1);
  V2vsPt->SetLineWidth(0.8);
  V2vsPt->SetLineStyle(2);
  V2vsPt->Draw("same");

  TGraphErrors *gr2 = new TGraphErrors(NCentBin,pt,v2a,0,v2ae);
      gr2->SetTitle("");
      gr2->SetMarkerStyle(20);
      gr2->SetMarkerSize(1.0);
      gr2->SetMarkerColor(1);
      gr2->SetLineWidth(2);
      gr2->SetLineColor(1);
  gr2->Draw("psameez");
 
  TGraphErrors *gr3 = new TGraphErrors(NCentBin,pta,v21,0,v21e);
      gr3->SetMarkerStyle(24);
      gr3->SetMarkerSize(1.);
      gr3->SetMarkerColor(2);
      gr3->SetLineWidth(2);
      gr3->SetLineColor(1);
  gr3->Draw("psameez");

  TGraphErrors *gr4 = new TGraphErrors(NCentBin,pta2,v22,0,v22e);
      gr4->SetMarkerStyle(25);
      gr4->SetMarkerSize(1.3);
      gr4->SetMarkerColor(4);
      gr4->SetLineWidth(2);
      gr4->SetLineColor(1);
  gr4->Draw("psameez");

  TGraphErrors *gr5 = new TGraphErrors(NCentBin,pta3,v23,0,v23e);
      gr5->SetMarkerStyle(28);
      gr5->SetMarkerSize(1);
      gr5->SetMarkerColor(3);
      gr5->SetLineWidth(2);
      gr5->SetLineColor(1);
  gr5->Draw("psameez");

    TLegend *leg = new TLegend(0.19,0.65,0.38,0.93);
    leg->SetFillColor(10);
    leg->SetBorderSize(0.035);
    leg->SetTextFont(42);
    leg->SetTextColor(1);
    leg->SetTextSize(0.04);
    //leg->SetLineStyle(0.06);
    leg->AddEntry(gr2,"|#Delta#eta|>0.0","p");
    leg->AddEntry(gr3,"|#Delta#eta|>1.0","p");
    leg->AddEntry(gr4,"|#Delta#eta|>1.5","p");
    leg->AddEntry(gr5,"|#Delta#eta|>2.0","p");
    leg->AddEntry(V2vsPt,"Input v_{2}","l");
    leg->Draw();
  
    TLatex *tex1= new TLatex(5,0.305,"STEG, 100 events");
    tex1->SetTextColor(1);
    tex1->SetTextSize(0.05);
    tex1->SetTextFont(42);
    tex1->Draw();

    TLatex *tex1= new TLatex(3,0.265,"90>N>=35,  |#eta| < 2.4");
    tex1->SetTextColor(1);
    tex1->SetTextSize(0.05);
    tex1->SetTextFont(42);
    //tex1->Draw();


  c1->Print("eventplanev2.gif");
  //c1->Print("eventplanev2.C");

}

