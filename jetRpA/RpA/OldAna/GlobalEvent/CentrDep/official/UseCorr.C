void UseCorr(){

TFile *f=TFile::Open("output/weight/pPbEpos_DS_Eff.root","ReadOnly");
TH1F* hEff=(TH1F*)f->Get("hHF_eff_HFEnergyPlus4");
cout<<hEff->GetBinContent(hEff->FindBin(4.50))<<endl;
}

