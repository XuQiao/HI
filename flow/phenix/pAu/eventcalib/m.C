#include "m.h"
#include <stdlib.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include "TString.h"
#include <fstream>
#include <iostream>

using namespace std;

m::m(TString input, TString input1, TString output):
 OutputFileName(output), 
 InputFileName(input), 
 InputFileName1(input1)
{
}

m::~m()
{
    delete tree;
    d_infile->Close();
    d_infile1->Close();
    d_outfile->Close();
}

int m::Init()
{
  d_infile = TFile::Open(InputFileName,"ReadOnly");
  d_infile1 = TFile::Open(InputFileName1,"ReadOnly");
  if(!d_infile) return 0;
  if(!d_infile1) return 0;
  tree = (TTree*)d_infile->Get("tree");
  tree1 = (TTree*)d_infile1->Get("tree");
  if(!tree) return 0;
  if(!tree1) return 0; 
  d_outfile = new TFile(OutputFileName,"recreate");
  newtree = new TTree("newtree",tree->GetName());
  Long64_t autos = 500000000;  //2G
  newtree -> SetAutoSave(autos);
  tree1 -> SetBranchAddress("cent", &cent);
  tree -> SetBranchAddress("nbbc", &nbbc);
  tree1 -> SetBranchAddress("nbbc", &nbbc1);
  tree1 -> SetBranchAddress("bbct0", &bbct0);
  newtree -> Branch("cent", &cent, "cent/I");
  newtree -> Branch("nbbc", &nbbc1, "nbbc/I");
  newtree -> Branch("bbct0", &bbct0, "bbct0[nbbc]/F");
  return 0;
}

int m::process_event()
{
  int nEvent = tree1->GetEntries();
  cout<<nEvent<<endl;
  int ievent = 0;
  for(int ievent1=0;ievent1 < nEvent; ievent1++){
      tree1->GetEntry(ievent1);
      tree->GetEntry(ievent);
  if(ievent%100000==0){
      std::cout<<"************* ievent= "<<ievent<<"    *************"<<std::endl;
  }

  if(nbbc==nbbc1){
    newtree->Fill();
    ievent++;
  }
}
return 0;
}


int m::End()
{

    std::cout << "InputFileName = " << InputFileName << std::endl;
    std::cout << "InputFileName1 = " << InputFileName1 << std::endl;
    std::cout << "OutputFileName = " << OutputFileName << std::endl;

  if(d_outfile) {
    d_outfile->cd();
    newtree->Write("tree");
  }
  else{
      std::cout<<"ERROR: No output file set!"<<std::endl;
  }

  return 0;
}
