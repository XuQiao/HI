#ifndef __PERFORMTEST_H__
#define __PERFORMTEST_H__

#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include "TChain.h"

#include "RpPar.h"

class PerformTestMB
{
 public:

  PerformTestMB(std::vector<TString> input,std::vector<TString> input1, const char* output="rpanase.root");
  virtual ~PerformTestMB();

  int Init();
  int Inittree();
  int process_event();
  int End();
  
 private:
  std::vector<TString> InputFileName;
  std::vector<TString> InputFileName1;
  const std::string OutputFileName;
  TFile *d_infile;
  int ievent;
  int jevent;
  int RunNumber;

  TFile *d_outfile;
  TChain *tree;
  TChain *tree1;
  //tree variables
  float        d_bbcz;    // bbcz
  int        centrality; // integer but stored as float in PHGlobal etc
  float        bbc_qn;
  float        bbc_qs;
  float        ZdcEn;
  float        ZdcEs;
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
  float        d_BBC_valid[128];

  int          npc1;
  int          d_nFVTX_clus;
  float        d_FVTX_x[10000];
  float        d_FVTX_y[10000];
  float        d_FVTX_z[10000];

  int          d_ntrk;
  float        d_mom[1000];
  float        d_phi0[1000];
  float        d_the0[1000];
  int          d_charge[1000];
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
//---histograms---------------

  TH1F* htrig;
//Tracks
  TH2F *hcntetaphi;
  TH1F *hcntpt;
  TH2F* pc3dphidz_arm0_pos[50];
  TH2F* pc3dphidz_arm1_pos[50];
  TH2F* pc3dphidz_arm0_neg[50];
  TH2F* pc3dphidz_arm1_neg[50];
  TH2F* pc3dphidz_arm0_pos_z[10][50];
  TH2F* pc3dphidz_arm1_pos_z[10][50];
  TH2F* pc3dphidz_arm0_neg_z[10][50];
  TH2F* pc3dphidz_arm1_neg_z[10][50];


//tof
  TH2F* tofdphidz[2];
  TH2F* tofwdphidz[2];
  TH2F* tofsdphisdz;
  TH2F* tofwsdphisdz;
  TH2F* ttofqpratio;
  TH2F* m2qpratio;
  TH2F* ttofp;
  TH2F* m2p;
  TH2F* pinv2chbeta;
  TH2F* deltattofeis;

//vtx
  TH2F* hcluetaphi[4];

//bbc
  TH1F* bbcet;
  TH1F* bbctdevsouth;
  TH1F* bbctdevnorth;
  TH1F* bbctsouth;
  TH1F* bbctnorth;

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

//mpc
 
TH2F* mpcetdis;
TH2F* mpcetetasouth[1];
TH2F* mpcetetanorth[1];
TH2F*  mpc_south_cent;
TH2F*  mpc_north_cent;
TH2F*  mpc_south_north;
TH2F *south_mpc_north_mpc;


//correlate
  TH2F* hcent1cent2[10];
  TH2F* hbbcsZdcEs;
  TH2F* hmixbbcsZdcEs;
  TH2F* hvtxzfvtxz;
  TH2F* hvtxxvtxz;
  TH2F* hvtxyvtxz;
  TH2F* hpc1hitsbbc;
  TH2F* hnpc3hitsntof;
  TH2F* hbbcnbbc;
  TH2F* hbbcsbbcn;
  TH2F* hnvtxnfvtxtrk[4];
  TH2F* hnvtxnmpc[4];
  TH2F* hbbcsnvtx[4];
  TH2F* hbbcnnvtx[4];
  TH2F* hbbcnvtx[4];
  TH2F* hnbbcnclu;
  TH2F* hntracknmpc;
  TH2F *hnfvtxtrkbbc;
  TH2F *hnfvtxtrksnmpcs;
  TH2F *hnfvtxtrknnmpcn;
  TH2F *south_mpc_south_bbc;
  TH2F *north_mpc_north_bbc;
  TH2F* hcentbbct0sigmasouth[10];
  TH2F* hcentbbct0sigmanorth[10];
  TH2F* hcentbbcmixt0sigmasouth[10];
  TH2F* hcentbbcmixt0sigmanorth[10];
  TH2F* hcentbbct0fracsouth[10];
  TH2F* hcentbbct0fracnorth[10];
  TH2F* hcentbbcmixt0fracsouth[10];
  TH2F* hcentbbcmixt0fracnorth[10];
  
  float    pi;
  
};

#endif /* __PERFORMTEST_H__ */
