#include "mm.h"
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

mm::mm(TString input, TString input1, TString output, TString output1):
 OutputFileName(output), 
 OutputFileName1(output1), 
 InputFileName(input), 
 InputFileName1(input1)
{
}

mm::~mm()
{
    delete tree;
    delete tree1;
    delete newtree;
    delete newtree1;
    d_infile->Close();
    d_infile1->Close();
    d_outfile->Close();
    d_outfile1->Close();
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
  d_outfile1 = new TFile(OutputFileName1,"recreate");
  newtree1 = new TTree("newtree1",tree->GetName());
  Long64_t autos = 500000000;  //2G
  newtree -> SetAutoSave(autos);
  newtree1 -> SetAutoSave(autos);
  tree->SetBranchAddress("pass_tick", &pass_tick);
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("trig", &trig);
  tree->SetBranchAddress("npc1hits", &npc1hits);
  tree->SetBranchAddress("ZdcEs", &ZdcEs);
  tree->SetBranchAddress("ZdcEn", &ZdcEn);
  tree->SetBranchAddress("bbc_s", &bbc_s);
  tree->SetBranchAddress("bbc_n", &bbc_n);
  tree->SetBranchAddress("bbcv", &bbcv);
  tree->SetBranchAddress("vtxz", &vtxz); 

  tree->SetBranchAddress("ntrack", &ntrack);
  tree->SetBranchAddress("nbbc", &nbbc);
  tree->SetBranchAddress("nvtx", &nvtx);
  tree->SetBranchAddress("nvtxtrack", &nvtxtrack);

  tree->SetBranchAddress("mom", &mom);
  tree->SetBranchAddress("phi0", &phi0);
  tree->SetBranchAddress("the0", &the0);
  tree->SetBranchAddress("alpha", &alpha);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("zed", &zed);
  tree->SetBranchAddress("charge", &charge);
  tree->SetBranchAddress("arm", &arm);
  tree->SetBranchAddress("pc3dphi", &pc3dphi);
  tree->SetBranchAddress("pc3dz", &pc3dz);
  tree->SetBranchAddress("slat", &slat);
  tree->SetBranchAddress("tofdphi", &tofdphi);
  tree->SetBranchAddress("tofdz", &tofdz);
  tree->SetBranchAddress("tofpl", &tofpl);
  tree->SetBranchAddress("qtof", &qtof);
  tree->SetBranchAddress("ttof", &ttof);

  tree->SetBranchAddress("bbccharge", &bbccharge);
  tree->SetBranchAddress("bbct0", &bbct0);
  tree->SetBranchAddress("bbcx", &bbcx);
  tree->SetBranchAddress("bbcy", &bbcy);
  tree->SetBranchAddress("bbcz", &bbcz);

  tree->SetBranchAddress("layer", &layer);
  tree->SetBranchAddress("vtxX", &vtxX);
  tree->SetBranchAddress("vtxY", &vtxY);
  tree->SetBranchAddress("vtxZ", &vtxZ);

  tree->SetBranchAddress("vtxnhits", &vtxnhits);
  tree->SetBranchAddress("vtxpx", &vtxpx);
  tree->SetBranchAddress("vtxpy", &vtxpy);
  tree->SetBranchAddress("vtxpz", &vtxpz);

  tree1 -> SetBranchAddress("bbcv", &bbcv1);
  tree1 -> SetBranchAddress("vtxz", &vtxz);
  tree1 -> SetBranchAddress("fvtxz", &fvtxz);
  tree1 -> SetBranchAddress("nfvtxtrack", &nfvtxtrack);
  tree1 -> SetBranchAddress("farm", &farm);
  tree1 -> SetBranchAddress("fnhits", &fnhits);
  tree1 -> SetBranchAddress("fthe", &fthe);
  tree1 -> SetBranchAddress("feta", &feta);
  tree1 -> SetBranchAddress("fphi", &fphi);
  tree1 -> SetBranchAddress("fvtxX", &fvtxX);
  tree1 -> SetBranchAddress("fvtxY", &fvtxY);
  tree1 -> SetBranchAddress("fvtxZ", &fvtxZ);

  newtree->Branch("pass_tick", &pass_tick, "pass_tick/O");
  newtree->Branch("run", &run, "run/I");
  newtree->Branch("trig", &trig, "trig/I");
  newtree->Branch("npc1hits", &npc1hits, "npc1hits/I");
  newtree->Branch("ZdcEs", &ZdcEs, "ZdcEs/F");
  newtree->Branch("ZdcEn", &ZdcEn, "ZdcEn/F");
  newtree->Branch("bbc_s", &bbc_s, "bbc_s/F");
  newtree->Branch("bbc_n", &bbc_n, "bbc_n/F");
  newtree->Branch("bbcv", &bbcv, "bbcv/F");
  newtree->Branch("vtxz", &vtxz, "vtxz/F"); 
  newtree->Branch("ntrack", &ntrack, "ntrack/I");
  newtree->Branch("nbbc", &nbbc, "nbbc/I");
  newtree->Branch("nvtx", &nvtx, "nvtx/I");
  newtree->Branch("nvtxtrack", &nvtxtrack, "nvtxtrack/I");

  newtree->Branch("mom", &mom, "mom[ntrack]/F");
  newtree->Branch("phi0", &phi0, "phi0[ntrack]/F");
  newtree->Branch("the0", &the0, "the0[ntrack]/F");
  newtree->Branch("alpha", &alpha, "alpha[ntrack]/F");
  newtree->Branch("phi", &phi, "phi[ntrack]/F");
  newtree->Branch("zed", &zed, "zed[ntrack]/F");
  newtree->Branch("charge", &charge, "charge[ntrack]/I");
  newtree->Branch("arm", &arm, "arm[ntrack]/I");
  newtree->Branch("pc3dphi", &pc3dphi, "pc3dphi[ntrack]/F");
  newtree->Branch("pc3dz", &pc3dz, "pc3dz[ntrack]/F");
  newtree->Branch("slat", &slat, "slat[ntrack]/I");
  newtree->Branch("tofdphi", &tofdphi, "tofdphi[ntrack]/F");
  newtree->Branch("tofdz", &tofdz, "tofdz[ntrack]/F");
  newtree->Branch("tofpl", &tofpl, "tofpl[ntrack]/F");
  newtree->Branch("qtof", &qtof, "qtof[ntrack]/F");
  newtree->Branch("ttof", &ttof, "ttof[ntrack]/F");

  newtree->Branch("bbccharge", &bbccharge, "bbccharge[nbbc]/F");
  newtree->Branch("bbct0", &bbct0, "bbct0[nbbc]/F");
  newtree->Branch("bbcx", &bbcx, "bbcx[nbbc]/F");
  newtree->Branch("bbcy", &bbcy, "bbcy[nbbc]/F");
  newtree->Branch("bbcz", &bbcz, "bbcz[nbbc]/F");

  newtree->Branch("layer", &layer, "layer[nvtx]/I");
  newtree->Branch("vtxX", &vtxX, "vtxX[nvtx]/F");
  newtree->Branch("vtxY", &vtxY, "vtxY[nvtx]/F");
  newtree->Branch("vtxZ", &vtxZ, "vtxZ[nvtx]/F");

  newtree->Branch("vtxnhits", &vtxnhits, "vtxnhits[nvtxtrack]/I");
  newtree->Branch("vtxpx", &vtxpx, "vtxpx[nvtxtrack]/F");
  newtree->Branch("vtxpy", &vtxpy, "vtxpy[nvtxtrack]/F");
  newtree->Branch("vtxpz", &vtxpz, "vtxpz[nvtxtrack]/F");

  newtree1 -> Branch("bbcv", &bbcv1,"bbcv/F");
  newtree1 -> Branch("vtxz", &vtxz,"vtxz/F");
  newtree1 -> Branch("fvtxz", &fvtxz,"fvtxz/F");
  newtree1 -> Branch("nfvtxtrack", &nfvtxtrack, "nfvtxtrack/I");
  newtree1 -> Branch("farm", &farm, "farm[nfvtxtrack]/I");
  newtree1 -> Branch("fnhits", &fnhits, "fnhits[nfvtxtrack]/I");
  newtree1 -> Branch("fthe", &fthe, "fthe[nfvtxtrack]/F");
  newtree1 -> Branch("feta", &feta, "feta[nfvtxtrack]/F");
  newtree1 -> Branch("fphi", &fphi, "fphi[nfvtxtrack]/F");
  newtree1 -> Branch("fvtxX", &fvtxX, "fvtxX[nfvtxtrack]/F");
  newtree1 -> Branch("fvtxY", &fvtxY, "fvtxY[nfvtxtrack]/F");
  newtree1 -> Branch("fvtxZ", &fvtxZ, "fvtxZ[nfvtxtrack]/F");

  return 0;
}

int mm::process_event()
{
  int nEvent = tree->GetEntries();
  int nEvent1 = tree1->GetEntries();
  int shift=7000;
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
    newtree1->Fill();
  }
  else{
    ievent++;
    tree->GetEntry(ievent);
      if(bbcv==bbcv1){
        newtree->Fill();
        newtree1->Fill();
      }
      else{
        ievent--;
        ievent1++;
        tree->GetEntry(ievent);
        tree1->GetEntry(ievent1);
      if(bbcv==bbcv1){
        newtree->Fill();
        newtree1->Fill();
      }
      else{
          ievent1--;
          ievent+=2;
          tree->GetEntry(ievent);
          if(bbcv==bbcv1){
             newtree->Fill();
             newtree1->Fill();
          }
          else{
              ievent-=2;
              ievent1+=2;
              tree->GetEntry(ievent);
              tree1->GetEntry(ievent1);
          if(bbcv==bbcv1){
             newtree->Fill();
             newtree1->Fill();
          }
          else{
              ievent1-=2;
              int i;
              for(i=0;i<shift;i++){
                  ievent1+=1;
                  tree1->GetEntry(ievent1);
                  if(bbcv==bbcv1) {
                    newtree->Fill();
                    newtree1->Fill();
                    break;
                  }
              }
              if(i==shift){
                  ievent1-=shift;
                  tree1->GetEntry(ievent1);
                  int j;
                  for(j=0;j<shift;j++){
                      ievent+=1;
                      tree->GetEntry(ievent);
                  if(bbcv==bbcv1) {
                    newtree->Fill();
                    newtree1->Fill();
                    break;
                  }
              }
                  if(j==shift){
                      ievent-=shift;
                  }
          }
          }
      }
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
    std::cout << "OutputFileName1 = " << OutputFileName1 << std::endl;

  if(d_outfile) {
    d_outfile->cd();
    newtree->Write("tree");
    d_outfile1->cd();
    newtree1->Write("tree");
  }
  else{
      std::cout<<"ERROR: No output file set!"<<std::endl;
  }

  return 0;
}
