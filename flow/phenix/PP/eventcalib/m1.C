#include "m1.h"
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

m1::m1(TString input, TString output):
 OutputFileName(output), 
 InputFileName(input)
{
}

m1::~m1()
{
    delete tree;
    d_infile->Close();
    d_outfile->Close();
}

int m1::Init()
{
  d_infile = TFile::Open(InputFileName,"ReadOnly");
  if(!d_infile) return 0;
  tree = (TTree*)d_infile->Get("tree");
  if(!tree) return 0;
  d_outfile = new TFile(OutputFileName,"recreate");
  newtree = new TTree("newtree",tree->GetName());
  Long64_t autos = 500000000;  //2G
  newtree -> SetAutoSave(autos);
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

  return 0;
}

int m1::process_event()
{
  int nEvent = tree->GetEntries();
  cout<<nEvent<<endl;
  for(int ievent=0;ievent < nEvent; ievent++){
      tree->GetEntry(ievent);
  if(ievent%100000==0){
      std::cout<<"************* ievent= "<<ievent<<"    *************"<<std::endl;
  }

  if(ievent!=4807488){
    newtree->Fill();
  }
}
return 0;
}


int m1::End()
{

    std::cout << "InputFileName = " << InputFileName << std::endl;
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
