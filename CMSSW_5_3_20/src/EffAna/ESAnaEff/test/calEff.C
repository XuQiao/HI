calEff(){

TFile *fileDS = new TFile("../../AnaEff/test/dsES.root");
TFile *fileES = new TFile("fesDS.root");
TH1D *hMultDS = (TH1D*) fileDS->Get("demo/hMultDS");;
TH1D *hMultES = (TH1D*) fileES->Get("demo/hMultDS");;

//hMultDS->Draw();
//hMultES->SetLineColor(2);
//hMultES->Draw("same");

  TCanvas *c1 = new TCanvas("c1","Efficiency Fraction",500,500);
  TGraphAsymmErrors *gEffMult = new TGraphAsymmErrors(); 
  gEffMult->SetName("gEffMult");
  gEffMult->BayesDivide(hMultES, hMultDS);
  gEffMult->SetMarkerStyle(25);
  gEffMult->SetLineStyle(2);
  gEffMult->SetLineColor(2);
  gEffMult->SetMarkerColor(2);
  gEffMult->Draw();
  //gEffMult->Print("all");

TFile *f = new TFile("ntrkEff.root","RECREATE");
f->cd();
hMultDS->Write();
hMultES->Write();
gEffMult->Write();

}
