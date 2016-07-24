double  Getmean(TString sys1, TString sys2, float cent1, float cent2){

  gROOT->SetStyle("Plain");
  TFile *infile = TFile::Open(Form("glau_%s%s_ntuple_1M.root",sys1.Data(),sys2.Data()));
  infile->cd();
  TH2F *pEcc2B = new TH2F("pEcc2B","#epsilon_{2} vs impact paramter in 200 GeV Collisions from Glauber MC;B;#epsilon_{n}",25200,-0.5,251.5,100,0,1);
  if(sys2.Contains("smeared")){
  TNtuple *nt = (TNtuple*)infile->Get("nt");
  nt->Project("pEcc2B","Ecc3G:B");
  }
  else{
  TNtuple *nt = (TNtuple*)infile->Get(Form("nt_%s_%s",sys1.Data(),sys2.Data()));
  nt->Project("pEcc2B","Ecc2:B");
  }
  TH1D* hB = (TH1D*)pEcc2B->ProjectionX("hB",0,-1);
  for(int ibin=0;ibin<hB->GetNbinsX();ibin++){
      if(hB->Integral(0,ibin) >= cent1 /100. * hB->Integral()){ int bin1 = ibin; break;}
  }
  for(int ibin=0;ibin<hB->GetNbinsX();ibin++){
      if(hB->Integral(0,ibin) >= cent2 /100. * hB->Integral()){ int bin2 = ibin; break;}
  }
  cout<< bin1 << "\t" << bin2 << endl;
  TH1D* hEcc2 = (TH1D*)pEcc2B->ProjectionY("hEcc2",bin1,bin2);
//  hEcc2->Draw();
  cout << hEcc2 -> GetEntries() << endl;
  double mean = hEcc2->GetMean();
  return mean;
}
void calcEcc(){
  double meanEccdAu = Getmean("d","Au_smeared",0,5..);
  double meanEccpAu = Getmean("p","Au_smeared",0,5.);
  double meanEccHe3Au = Getmean("He3","Au_smeared",0,5.);
  double meanEccpAl = Getmean("p","Al",0,5);
  double meanEccpAl2 = Getmean("p","Al2",0,5);
  double meanEccAuAu = Getmean("Au","Au_smeared",0,5.);

  cout << "pAu smeared mean Ecc2 = " << meanEccpAu << endl;
  cout << "dAu smeared mean Ecc2 = " << meanEccdAu << endl;
  cout << "He3Au smeared mean Ecc2 = " << meanEccHe3Au << endl;
  cout << "pAl mean Ecc2 = " << meanEccpAl << endl;
  cout << "pAl2 mean Ecc2 = " << meanEccpAl2 << endl;
  cout << "AuAu smeared mean Ecc2 = " << meanEccAuAu << endl;
}
