void Plot2DCorrelation(){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  Int_t nFlagSaveFIGs = 0;
  TFile *f = new TFile("correlationSTEG.root");

  TH2D *hCorrelation = (TH2D*) f->Get("correlation");

   TCanvas *c = new TCanvas("c", "",10,73,530,500);
   gStyle->SetOptStat(0);
   gStyle->SetPalette(1);
   //c->Range(-1.463631,-1.457537,0.9555939,1.145208);
   //TView *view = TView::CreateView(1);
   //view->SetRange(-4.05,-1.472621,1.235,4.05,4.614213,1.55525);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(10);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetTheta(60.839);
   c->SetPhi(38.0172);
   c->SetLeftMargin(0.27);
   c->SetRightMargin(0.06);
   c->SetTopMargin(0.06);
   c->SetBottomMargin(0.18);
   c->SetFrameFillStyle(0);
   c->SetFrameLineStyle(0);
   c->SetFrameBorderMode(0);

   //hCorrelation->SetMaximum(2.6);
   //hCorrelation->SetMinimum(-1.2);

   hCorrelation->SetTitle("");
   hCorrelation->SetLineWidth(3);
   hCorrelation->SetMarkerStyle(20);
   hCorrelation->SetMarkerSize(1.2);
   hCorrelation->GetXaxis()->SetTitle("#Delta#eta");
   hCorrelation->GetXaxis()->SetRange(4,30);
   hCorrelation->GetXaxis()->CenterTitle(true);
   hCorrelation->GetXaxis()->SetNdivisions(505);
   hCorrelation->GetXaxis()->SetLabelFont(42);
   hCorrelation->GetXaxis()->SetLabelSize(0.07);
   hCorrelation->GetXaxis()->SetTitleSize(0.07);
   hCorrelation->GetXaxis()->SetTitleFont(42);
   hCorrelation->GetYaxis()->SetTitle("#Delta#phi");
   hCorrelation->GetYaxis()->CenterTitle(true);
   hCorrelation->GetYaxis()->SetNdivisions(505);
   hCorrelation->GetYaxis()->SetLabelFont(42);
   hCorrelation->GetYaxis()->SetLabelSize(0.07);
   hCorrelation->GetYaxis()->SetTitleSize(0.07);
   hCorrelation->GetYaxis()->SetTitleFont(42);
   hCorrelation->GetZaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}");
   hCorrelation->GetZaxis()->CenterTitle(true);
   hCorrelation->GetZaxis()->SetNdivisions(505);
   hCorrelation->GetZaxis()->SetLabelFont(42);
   hCorrelation->GetZaxis()->SetLabelSize(0.0595);
   hCorrelation->GetZaxis()->SetTitleSize(0.054);
   hCorrelation->GetZaxis()->SetTitleOffset(1.804);
   hCorrelation->GetZaxis()->SetTitleFont(42);
   hCorrelation->Draw("surf1 fb");

   TLatex *   tex = new TLatex(0.7,0.92,"CMS Preliminary");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(20);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.055,0.925,"pp  #sqrt{s} = 7 TeV, N #geq 110");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(20);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.055,0.84,"5 < p_{T}^{trig} < 6 GeV/c");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(20);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.055,0.78,"1 < p_{T}^{assoc} < 2 GeV/c");
tex->SetNDC();
   tex->SetTextFont(43);
   tex->SetTextSize(20);
   tex->SetLineWidth(2);
   tex->Draw();
   c->Modified();
   c->cd();
   c->SetSelected(c);

   if(nFlagSaveFIGs==1) c->SaveAs("hCorrelation.png");







}
