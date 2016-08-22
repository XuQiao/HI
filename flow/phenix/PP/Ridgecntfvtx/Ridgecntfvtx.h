#ifndef __RIDGECNTFVTX_H__
#define __RIDGECNTFVTX_H__

#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 

class Ridgecntfvtx
{
 public:

  Ridgecntfvtx(TString input = "",TString output="rpanase.root");
  ~Ridgecntfvtx();

  int Init();
  int process_event();
  int End();
  
 private:
  TString OutputFileName;
  TString InputFileName;
  TString coll;
  TFile *d_infile;
  TFile *d_outfile;

  TTree *tree;

//event variables
  int run;
  int trig;
  int cent;
  int npc1hits;
  float ZdcEs;
  float ZdcEn;
  float bbc_s;
  float bbc_n;
  float bbcv;
  float vtxz;
  float fvtxz;
  int ntrack;
  int nbbc;
  int nvtx;
  int nmpc;
  int nvtxtrack;
  int nfvtxtrack;
  
  //QVector
  float Qx[72];
  float Qy[72];
  float Qw[72];

  //cnt track variables
  float mom[1000];
  float phi0[1000];
  float the0[1000];
  float phi[1000];
  float alpha[1000];
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

//----------------------------------------------------------------------------
    
  //BBC
  //normal
  TH1F *hforesouthbbc[10][40];
  TH1F *hbacksouthbbc[10][40];
  
  TH1F *hforenorthbbc[10][40];
  TH1F *hbacknorthbbc[10][40];
//flip
  TH1F *kforesouthbbc[10][40];
  TH1F *kbacksouthbbc[10][40];
  TH2F *kforesouthetabbc[10][40];
  TH2F *kbacksouthetabbc[10][40];

  TH1F *kforenorthbbc[10][40];
  TH1F *kbacknorthbbc[10][40];
  TH2F *kforenorthetabbc[10][40];
  TH2F *kbacknorthetabbc[10][40];
  
 //with weight 
  TH1F *hforesouthbbcw[10][40];
  TH1F *hbacksouthbbcw[10][40];
  
  TH1F *hforenorthbbcw[10][40];
  TH1F *hbacknorthbbcw[10][40];

  TH1F *kforesouthbbcw[10][40];
  TH1F *kbacksouthbbcw[10][40];
  TH2F *kforesouthetabbcw[10][40];
  TH2F *kbacksouthetabbcw[10][40];

  TH1F *kforenorthbbcw[10][40];
  TH1F *kbacknorthbbcw[10][40];
  TH2F *kforenorthetabbcw[10][40];
  TH2F *kbacknorthetabbcw[10][40];

  //tower mixing
  TH1F *hbacksouthbbc2[10][40];
  TH1F *hbacknorthbbc2[10][40];
  TH1F *kbacksouthbbc2[10][40];
  TH1F *kbacknorthbbc2[10][40];
  TH2F *kbacksouthetabbc2[10][40];
  TH2F *kbacknorthetabbc2[10][40];
//with weight
  TH1F *hbacksouthbbcw2[10][40];
  TH1F *hbacknorthbbcw2[10][40];
  TH1F *kbacksouthbbcw2[10][40];
  TH1F *kbacknorthbbcw2[10][40];
  TH2F *kbacksouthetabbcw2[10][40];
  TH2F *kbacknorthetabbcw2[10][40];
  
  float pi;
};

#endif /* __RIDGECNTFVTX_H__ */
