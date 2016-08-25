#ifndef __mm_H__
#define __mm_H__

#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TProfile.h" 

class mm
{
 public:

  mm(TString input = "", TString input1 = "", TString output="rpanase.root", TString output1="");
  ~mm();

  int Init();
  int process_event();
  int End();
  
 private:
  TString OutputFileName;
  TString OutputFileName1;
  TString InputFileName;
  TString InputFileName1;
  TString coll;
  TFile *d_infile;
  TFile *d_infile1;
  TFile *d_outfile;
  TFile *d_outfile1;

  TTree *tree;
  TTree *tree1;
  TTree *newtree;
  TTree *newtree1;
  bool pass_tick;
  int run;
  int trig;
  int cent;
  int npc1hits;
  float ZdcEs;
  float ZdcEn;
  float bbc_s;
  float bbc_n;
  float bbcv;
  float bbcv1;
  float vtxz;
  float fvtxz;
  int ntrack;
  int nbbc;
  int nvtx;
  int nmpc;
  int nvtxtrack;
  int nfvtxtrack;

  //event plane object
  float Qx[72];
  float Qy[72];
  float Qw[72];

  //cnt track variables
  float mom[1000];
  float phi0[1000];
  float the0[1000];
  float alpha[1000];
  float phi[1000];
  float zed[1000];
  int charge[1000];
  int arm[1000];
  float pc3dphi[1000];
  float pc3dz[1000];
  
  //tof
  int slat[1000];
  float tofdphi[1000];
  float tofdz[1000];
  float tofpl[1000];
  float qtof[1000];
  float ttof[1000];

  //bbc
  float bbccharge[128];
  float bbct0[128];
  float bbcx[128];
  float bbcy[128];
  float bbcz[128];

  //vtx
  int layer[800];
  float vtxX[800];
  float vtxY[800];
  float vtxZ[800];

  //vtx track
  int vtxnhits[800];
  float vtxpx[800];
  float vtxpy[800];
  float vtxpz[800];

  //mpc
  float mpc_e[800];
  float mpc_x[800];
  float mpc_y[800];
  float mpc_z[800];

  //fvtx track
  int farm[1000];
  int fnhits[1000];
  float fthe[1000];
  float feta[1000];
  float fphi[1000];
  float fvtxX[1000];
  float fvtxY[1000];
  float fvtxZ[1000];
};

#endif /* __mm_H__ */
