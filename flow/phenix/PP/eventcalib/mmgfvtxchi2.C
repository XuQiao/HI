#include "mmgfvtxchi2.h"
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

mm::mm(TString input, TString input1, TString output):
 OutputFileName(output), 
 InputFileName(input), 
 InputFileName1(input1)
{
}

mm::~mm()
{
    delete tree;
    delete tree1;
    delete newtree;
    d_infile->Close();
    d_infile1->Close();
    d_outfile->Close();
}

int mm::Init()
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
  tree->SetBranchAddress("bbcv", &bbcv);

  tree1 -> SetBranchAddress("bbcv", &bbcv1);
  tree1 -> SetBranchAddress("nfvtxtrack", &nfvtxtrack);
  tree1 -> SetBranchAddress("fvtxchi2", &fvtxchi2);

  newtree -> Branch("nfvtxtrack", &nfvtxtrack, "nfvtxtrack/I");
  newtree -> Branch("fvtxchi2", &fvtxchi2, "fvtxchi2[nfvtxtrack]/F");

  return 0;
}

int mm::process_event()
{
  int nEvent = tree->GetEntries();
  int nEvent1 = tree1->GetEntries();
  int shift=10000;
  cout<<nEvent<<"\t"<<nEvent1<<endl;
  int ievent = 0;
  int ievent1 = 0;
  while(ievent<nEvent && ievent<nEvent1){
      tree1->GetEntry(ievent1);
      tree->GetEntry(ievent);
  if(ievent%100000==0){
      std::cout<<"************* ievent= "<<ievent<<"    *************"<<std::endl;
  }

  if(bbcv==bbcv1){
    newtree->Fill();
  }
  else{
    int i;
    for(i=0;i<shift;i++){
    ievent++;
    tree->GetEntry(ievent);
    if(bbcv==bbcv1) {
        newtree->Fill();
        break;
    }
    }
    if(i==shift){
        ievent-=shift;
        tree->GetEntry(ievent);
        int j;
        for(j=0;j<shift;j++){
            ievent1++;
            tree1->GetEntry(ievent1);
            if(bbcv==bbcv1) {
                newtree->Fill();
                break;
            }
        }
            if(j==shift){
                ievent1-=shift;
            }
          }
    }
        ievent++;
        ievent1++;
}

return 0;
}


int mm::End()
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
