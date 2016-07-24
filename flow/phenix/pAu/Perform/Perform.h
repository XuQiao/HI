#ifndef __PERFORM_H__
#define __PERFORM_H__

#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TProfile.h" 

class Perform
{
 public:

  Perform(TString input = "", TString input1 = "", TString input2="", TString output="rpanase.root");
  ~Perform();

  int Init();
  int process_event();
  int End();
  
 private:
  TString OutputFileName;
  TString InputFileName;
  TString InputFileName1;
  TString InputFileName2;
  TString coll;
  TFile *d_infile;
  TFile *d_infile1;
  TFile *d_infile2;
  TFile *d_outfile;

  TTree *tree;
  TTree *tree1;
  TTree *tree2;

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
  int nbbc1;
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

//-----------------------------------------------------------------------------------------

//Tracks
  TH2F *hcntetaphi;
  TH2F *hcentcnteta;
  TH2F *hcentcntpt;
  TH1F *hcntpt;
  TH2F* pc3dphidz_arm0_pos[50];
  TH2F* pc3dphidz_arm1_pos[50];
  TH2F* pc3dphidz_arm0_neg[50];
  TH2F* pc3dphidz_arm1_neg[50];
  TH2F* pc3dphidz_arm0_pos_z[4][50];
  TH2F* pc3dphidz_arm1_pos_z[4][50];
  TH2F* pc3dphidz_arm0_neg_z[4][50];
  TH2F* pc3dphidz_arm1_neg_z[4][50];


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
  TH2F* deltattofwis;

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
 
TH2F* mpcetetasouth[10];
TH2F* mpcetetanorth[10];
TH2F*  mpc_south_cent;
TH2F*  mpc_north_cent;
TH2F*  mpc_south_north;
TH2F *south_mpc_north_mpc;


//correlate
  TH2F* hcent1cent2[10];
  TH2F* hbbcsZdcEs;
  TH2F* hmixbbcsZdcEs;
  TH2F* hcentbbcs;
  TH2F* hcentntrk;
  TH2F* hcentbbct0sigmasouth[10];
  TH2F* hcentbbct0sigmanorth[10];
  TH2F* hcentbbcmixt0sigmasouth[10];
  TH2F* hcentbbcmixt0sigmanorth[10];
  TH2F* hcentbbct0fracsouth[10];
  TH2F* hcentrabbct0fracsouth[10];
  TH2F* hcentbbct0fracnorth[10];
  TH2F* hcentbbcmixt0fracsouth[10];
  TH2F* hcentbbcmixt0fracnorth[10];

  TH2F* hvtxzfvtxz;
  TH2F* hpc1hitsbbc;
  TH2F* hnpc3hitsntof;
  TH2F* hbbcnbbc;
  TH2F* hbbcsbbcn;
  TH2F* hnvtxnfvtxtrk[4];
  TH2F* hnvtxnmpc[4];
  TH2F* hnvtxntrk[4];
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
  
  //run QA
  TH2F *hrunbbcs;
  TH2F *hrunbbcn;
  TH2F *hrunntrack[5];
  
};

#endif /* __PERFORM_H__ */
