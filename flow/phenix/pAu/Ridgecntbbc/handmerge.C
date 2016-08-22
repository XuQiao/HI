#include <iostream>
#include <fstream>

void handmerge() {
 
  TH2D *hnvtxntrk_0;

  string PATH = "/phenix/plhf/peng/taxi/Run15pp200MBPro104/6354/data/";
  string file;
  TString filename;
  ifstream fin("run_by_run_output");
  
    filename=Form("%s%s",PATH,"421716.root");
    TFile *read_file = TFile::Open(filename);
  hnvtxntrk_0 = (TH2D*)read_file->Get("hnvtxntrk_0");
  hnvtxntrk_0->Reset();
  hnvtxntrk_0->SetName("hnvtxntrk_0");
  hnvtxntrk_0->SetTitle("hnvtxntrk_0");
  while (getline(fin,file)) {
    filename=Form("%s%s",PATH,file);
    cout << filename << endl;
    TFile *read_file = TFile::Open(filename);
    if (read_file) {
      TH2D *temp = (TH2D*)read_file->Get("hnvtxntrk_0");
      hnvtxntrk_0->Add(temp);
    }
    delete read_file;
  }

  TFile f("output.root","new");
  hnvtxntrk_0->Write();
  f.Close();
}
