#ifndef __RIDGEDAURUN16_H__
#define __RIDGEDAURUN16_H__

#include <string>
#include "TChain.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 

class RidgedAuRun16
{
 public:

  RidgedAuRun16(std::vector<TString> input, std::vector<TString> input1, const char* output="rpanase.root");
  virtual ~RidgedAuRun16();

  int Init();
  int Inittree();
  int process_event();
  int End();
  
 private:
  const std::string OutputFileName;
  std::vector<TString> InputFileName;
  std::vector<TString> InputFileName1;
  int ievent;
  int jevent;
  int nevent;
  int mbuff;
  int RunNumber;

  TFile *d_outfile;
  TChain *tree;
  TChain *tree1;
  //tree variables
  float        d_bbcz;    // bbcz
  int        centrality; // integer but stored as float in PHGlobal etc
  float        bbc_qn;
  float        bbc_qs;
  unsigned int trigger_scaled;
  unsigned int trigger_live;
  float        bc_x;
  float        bc_y;
//  float        vtx_z;
  float        eventfvtx_x;
  float        eventfvtx_y;
  float        eventfvtx_z;
  float        d_Qx[9];
  float        d_Qy[9];
  float        d_Qw[9];
  float        d_BBC_charge[128];
  float        d_BBC_time0[128];
  bool        d_BBC_valid[128];

  int          npc1;
  int          d_nFVTX_clus;
  float        d_FVTX_x[10000];
  float        d_FVTX_y[10000];
  float        d_FVTX_z[10000];

  int          d_ntrk;
  float        d_mom[1000];
  float        d_phi0[1000];
  float        d_the0[1000];
  float        d_charge[1000];
  float        d_pc3dphi[1000];
  float        d_pc3dz[1000];

  int          d_nfvtxtrk;
  int          d_nfvtxtrk1; //for test
  float        d_fvtxchi[1000];
  int          d_farm[1000];
  int          d_fnhits[1000];
  float        d_feta[1000];
  float        d_fphi[1000];
  float        d_fvtxX[1000];
  float        d_fvtxY[1000];
  float        d_fvtxZ[1000];

  //BBC
  //normal
  TH1D *hforesouthbbc[10][40];
  TH1D *hbacksouthbbc[10][40];

  //flip
  TH1D *kforesouthbbc[10][40];
  TH1D *kbacksouthbbc[10][40];
  
 //with weight 
  TH1D *hforesouthbbcw[10][40];
  TH1D *hbacksouthbbcw[10][40];
  
  TH1D *kforesouthbbcw[10][40];
  TH1D *kbacksouthbbcw[10][40];
  
  TH1D *hbacksouthbbc2[10][40];
  TH1D *kbacksouthbbc2[10][40];
  
  TH1D *hbacksouthbbcw2[10][40];
  TH1D *kbacksouthbbcw2[10][40];

  //north
  TH1D *hforenorthbbc[10][40];
  TH1D *hbacknorthbbc[10][40];
  
  //flip
  TH1D *kforenorthbbc[10][40];
  TH1D *kbacknorthbbc[10][40];

  //with weight
  TH1D *hforenorthbbcw[10][40];
  TH1D *hbacknorthbbcw[10][40];

  TH1D *kforenorthbbcw[10][40];
  TH1D *kbacknorthbbcw[10][40];
  //tower mixing
  
  TH1D *hbacknorthbbc2[10][40];
  TH1D *kbacknorthbbc2[10][40];

  //with weight
  TH1D *hbacknorthbbcw2[10][40];
  TH1D *kbacknorthbbcw2[10][40];

//  TH2F *kforesouthetabbc[10][40];
//  TH2F *kbacksouthetabbc[10][40];

//  TH2F *kforenorthetabbc[10][40];
//  TH2F *kbacknorthetabbc[10][40];

//  TH2F *kforesouthetabbcw[10][40];
//  TH2F *kbacksouthetabbcw[10][40];

//  TH2F *kforenorthetabbcw[10][40];
//  TH2F *kbacknorthetabbcw[10][40];
  
//  TH2F *kbacksouthetabbc2[10][40];
//  TH2F *kbacknorthetabbc2[10][40];

//  TH2F *kbacksouthetabbcw2[10][40];
//  TH2F *kbacknorthetabbcw2[10][40];

  //north
  TH1D *hforesnbbc[10][40];
  TH1D *hbacksnbbc[10][40];
  
  //flip
  TH1D *kforesnbbc[10][40];
  TH1D *kbacksnbbc[10][40];

  //with weight
  TH1D *hforesnbbcw[10][40];
  TH1D *hbacksnbbcw[10][40];

  TH1D *kforesnbbcw[10][40];
  TH1D *kbacksnbbcw[10][40];
  //tower mixing
  
  TH1D *hbacksnbbc2[10][40];
  TH1D *kbacksnbbc2[10][40];

  //with weight
  TH1D *hbacksnbbcw2[10][40];
  TH1D *kbacksnbbcw2[10][40];

  float    pi;
  
  TH2F* hcentbbcs;
  TH2F* hcentnfvtxs;
  TProfile* hrunbbcs[10];
  TProfile* hrunnfvtxs[10];
  TProfile* hrunntrack[10];
  TH2F* hbbcsnfvtxs[10];
//  TH2F* bbcphi;
//  TH2F* cntphi;
  
//fvtx
  TH2D *DCAxydis[2];
  TH2D *DCAxy2dis[2];
  TH2D *DCAcentdis[2];

  TH1F *fvtxdphidis;
  TH1F *fvtxdphidis2;

  TH2D *hvtx0etaz;
  TH2D *hvtx1etaz;

  TH2D *hvtx0etaphi;
  TH2D *hvtx1etaphi;

};

#endif /* __RIDGEDAURUN16_H__ */
