TProfile* Getplot(TString sys1, TString sys2, int size, int color, int style, int rebin){

  gROOT->SetStyle("Plain");
  TFile *infile = TFile::Open(Form("glau_%s%s_ntuple_1M.root",sys1.Data(),sys2.Data()));
  infile->cd();
  TProfile *pEcc2Npart = new TProfile("pEcc2Npart","#epsilon_{n} vs Npart in 200 GeV Collisions from Glauber MC;Npart;#epsilon_{n}",252,-0.5,251.5);
  TNtuple *nt = (TNtuple*)infile->Get(Form("nt_%s_%s",sys1.Data(),sys2.Data()));
  nt->Project("pEcc2Npart","Ecc2G:Npart");
  pEcc2Npart->SetMarkerSize(size);
  pEcc2Npart->SetMarkerColor(color);
  pEcc2Npart->SetMarkerStyle(style);
  pEcc2Npart->SetLineColor(color);
  pEcc2Npart->SetStats(0);
  int rebins=rebin;
  pEcc2Npart->Rebin(rebins);
  return pEcc2Npart;
}
void plot(){
  TCanvas *c1 = new TCanvas("c1","c1",100,100,610,500);
  c1->cd();
  TProfile* pdauEcc2Npart = Getplot("d","Au",1.4,1,20,1);
  TProfile* phe3auEcc2Npart = Getplot("He3","Au",1.4,2,20,1);
  TProfile* ppalEcc2Npart = Getplot("p","Al",1.4,4,33,1);
  TProfile* ppal2Ecc2Npart = Getplot("p","Al2",1.4,5,34,1);
  TProfile* ppauEcc2Npart = Getplot("p","Au",1.4,6,29,1);
  pdauEcc2Npart->Draw("p");
  pdauEcc2Npart->GetXaxis()->SetRangeUser(0,55);
  phe3auEcc2Npart->Draw("psame");
  ppalEcc2Npart->Draw("psame");
  ppal2Ecc2Npart->Draw("psame");
  ppauEcc2Npart->Draw("psame");
  TLegend *leg = new TLegend(0.65,0.6,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(pdauEcc2Npart,"d+Au","lp");
  leg->AddEntry(phe3auEcc2Npart,"He3+Au","lp");
  leg->AddEntry(ppalEcc2Npart,"p+Al","lp");
  leg->AddEntry(ppal2Ecc2Npart,"p+Al2","lp");
  leg->AddEntry(ppauEcc2Npart,"p+Au","lp");
  leg->Draw("same");

  /*
  TProfile *pEcc3Npart = new TProfile("pEcc3Npart","#epsilon_{n} vs Npart in 200 GeV CuAu;Npart;#epsilon_{n}",252,-0.5,251.5);
  nt_Cu_Au->Project("pEcc3Npart","Ecc3:Npart");
  pEcc3Npart->SetMarkerSize(1);
  pEcc3Npart->SetMarkerColor(2);
  pEcc3Npart->SetMarkerStyle(25);
  pEcc3Npart->SetLineColor(2);
  pEcc3Npart->SetStats(0);
  pEcc3Npart->Rebin(rebins);
  pEcc3Npart->Draw("psame");
  TLegend *leg = new TLegend(0.7,0.7,0.88,0.86);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.055);
  leg->AddEntry(pEcc2Npart,"#epsilon_{2}","p");
  leg->AddEntry(pEcc3Npart,"#epsilon_{3}","p");
  leg->Draw();
*/
  c1->SaveAs("epsilonNvsNpart_CuAu.png");
}
