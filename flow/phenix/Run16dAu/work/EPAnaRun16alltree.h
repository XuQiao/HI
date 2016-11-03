#ifndef __EPANARUN16TREE_H__
#define __EPANARUN16TREE_H__

#include <vector>
#include <list>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"

#include "RpPar.h"
#include "PPBBCpmtweight.h"

using namespace std;

class EPAnaRun16alltree
{
 public:

  EPAnaRun16alltree(std::vector<TString> input, std::vector<TString> input1, const char* output="rpanase.root");
  virtual ~EPAnaRun16alltree();

  int Init();
  int Inittree();
  int process_event();
  int End();
  int SetcalFlag(int _flag);
  void Getrec();
  void Getflt();
    

 protected:

 private:
  std::string OutputFileName;
  std::vector<TString> InputFileName;
  std::vector<TString> InputFileName1;
  TFile *d_infile;
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

  int          d_nbbc;
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


  int ievent;
  int jevent;
 // int nevnt;
  int calFlag;
 //std::vector<int> nevent;
  //int iseg;
  //int ndst;
  int RunNumber;
  //bool seqcal;
  std::vector<std::list<std::string> > filelist;
  //std::vector<std::string> fnames;

  TH1F* sinflt[ncent][nbbcz][nhar][nsub][nord];
  TH1F* cosflt[ncent][nbbcz][nhar][nsub][nord];
  
//  std::vector<std::vector<TH1F*> > sinflttmp[ncent][nbbcz][nhar][nsub][nord];
//  std::vector<std::vector<TH1F*> > cosflttmp[ncent][nbbcz][nhar][nsub][nord];

  TH1F* q[ncent][nbbcz][nhar][nsub][nxy];
//  std::vector<std::vector<TH1F*> > qtmp[ncent][nbbcz][nhar][nsub][nxy];  //invidual calibration
  TH1F* qRec[ncent][nbbcz][nhar][nsub][nxy];
  TH1F* psi[ncent][nbbcz][nhar][nsub];
//  std::vector<std::vector<TH1F*> > psitmp[ncent][nbbcz][nhar][nsub];
  TH1F* psiFla[ncent][nbbcz][nhar][nsub];

   TH1F* phiweight[ncent][nbbcz][nr][nhar][nsub];
   TProfile* phiweightbbc[ncent][nbbcz][nhar][nsub];

  float meanx[ncent][nbbcz][nhar][nsub];
  float meany[ncent][nbbcz][nhar][nsub];
  float rmsx[ncent][nbbcz][nhar][nsub];
  float rmsy[ncent][nbbcz][nhar][nsub];

  float cosfltarr[ncent][nbbcz][nhar][nsub][nord];
  float sinfltarr[ncent][nbbcz][nhar][nsub][nord];

  /*
  TH1F* EPRCNTBBCS[ncent][nhar];
  TH1F* EPRCNTFVTX1S[ncent][nhar];
  TH1F* EPRBBCSFVTX1S[ncent][nhar];
  TH1F* EPRCNTFVTX1S[ncent][nhar];
  TH1F* EPRCNTFVTX2S[ncent][nhar];
  TH1F* EPRCNTFVTX1LS[ncent][nhar];
  TH1F* EPRCNTFVTX2LS[ncent][nhar];
  TH1F* EPRCNTFVTX3LS[ncent][nhar];
  TH1F* EPRCNTFVTX4LS[ncent][nhar];
  TH1F* EPRBBCSFVTX2S[ncent][nhar];
  */
  TH1F* EPRFVTX1NFVTX1S[ncent][nhar];
  TH1F* EPRBBCSFVTX1S[ncent][nhar];
  TH1F* EPRBBCSFVTX1LS[ncent][nhar];
  TH1F* EPRBBCSFVTX2LS[ncent][nhar];
  TH1F* EPRBBCSFVTX3LS[ncent][nhar];
  TH1F* EPRBBCSFVTX4LS[ncent][nhar];
  TH1F* EPRBBCSFVTX1p2p3LS[ncent][nhar];
  TH1F* EPRBBCSFVTX1p2p4LS[ncent][nhar];
  TH1F* EPRBBCSFVTXtrkS[ncent][nhar];
  /*
  TH1F* EPRCNTFVTX1trkS[ncent][nhar];
  TH1F* EPRCNTFVTX2trkS[ncent][nhar];

  TH1F* EPRCNTcBBCSc[ncent][nhar];
  TH1F* EPRCNTcFVTX1S[ncent][nhar];

  TH1F* EPRBBCSFVTX1trkS[ncent][nhar];
  TH1F* EPRBBCSFVTX2trkS[ncent][nhar];
          
  TH1F* EPRBBCScFVTX1S[ncent][nhar];
  TH1F* EPRCNTcBBCSc[ncent][nhar];
  TH1F* EPRCNTcFVTX1Sc[ncent][nhar];
  TH1F* EPRBBCScFVTX1Sc[ncent][nhar];
  */
/*
  TH2F* vobsFVTX2S[ncent][nhar][nphi];
  TH2F* vobsFVTX2Ssq[ncent][nhar][nphi];
 */ 

  TH2F* vobsBBCS[ncent][nhar][nphi];
  TH2F* vobsFVTX1N[ncent][nhar][nphi];
  TH2F* vobsFVTX1S[ncent][nhar][nphi];
  TH2F* vobsFVTX1LS[ncent][nhar][nphi];
  TH2F* vobsFVTX2LS[ncent][nhar][nphi];
  TH2F* vobsFVTX3LS[ncent][nhar][nphi];
  TH2F* vobsFVTX4LS[ncent][nhar][nphi];
  TH2F* vobsFVTX1p2p3LS[ncent][nhar][nphi];
  TH2F* vobsFVTX1p2p4LS[ncent][nhar][nphi];
  TH2F* vobsFVTXtrkS[ncent][nhar][nphi];
  TH2F* vobsFVTX1LN[ncent][nhar][nphi];
  TH2F* vobsFVTX2LN[ncent][nhar][nphi];
  TH2F* vobsFVTX3LN[ncent][nhar][nphi];
  TH2F* vobsFVTX4LN[ncent][nhar][nphi];
  
  TH2F* vBBCS[ncent][nhar][nphi];
  TH2F* vFVTX1N[ncent][nhar][nphi];
  TH2F* vFVTX1S[ncent][nhar][nphi];
  TH2F* vFVTX1LS[ncent][nhar][nphi];
  TH2F* vFVTX2LS[ncent][nhar][nphi];
  TH2F* vFVTX3LS[ncent][nhar][nphi];
  TH2F* vFVTX4LS[ncent][nhar][nphi];
  TH2F* vFVTX1p2p3LS[ncent][nhar][nphi];
  TH2F* vFVTX1p2p4LS[ncent][nhar][nphi];
  TH2F* vFVTXtrkS[ncent][nhar][nphi];
  TH2F* vFVTX1LN[ncent][nhar][nphi];
  TH2F* vFVTX2LN[ncent][nhar][nphi];
  TH2F* vFVTX3LN[ncent][nhar][nphi];
  TH2F* vFVTX4LN[ncent][nhar][nphi];

  TH2F* vnBBCS[ncent][nhar][nphi];
  TH2F* vnFVTX1N[ncent][nhar][nphi];
  TH2F* vnFVTX1S[ncent][nhar][nphi];
  TH2F* vnFVTX1LS[ncent][nhar][nphi];
  TH2F* vnFVTX2LS[ncent][nhar][nphi];
  TH2F* vnFVTX3LS[ncent][nhar][nphi];
  TH2F* vnFVTX4LS[ncent][nhar][nphi];
  TH2F* vnFVTX1p2p3LS[ncent][nhar][nphi];
  TH2F* vnFVTX1p2p4LS[ncent][nhar][nphi];
  TH2F* vnFVTXtrkS[ncent][nhar][nphi];
  TH2F* vnFVTX1LN[ncent][nhar][nphi];
  TH2F* vnFVTX2LN[ncent][nhar][nphi];
  TH2F* vnFVTX3LN[ncent][nhar][nphi];
  TH2F* vnFVTX4LN[ncent][nhar][nphi];
 /* 
  TH2F* vobsBBCS[ncent][nhar][nphi];
  TH2F* vobsBBCSsq[ncent][nhar][nphi];
  TH2F* vobsFVTX1S[ncent][nhar][nphi];
  TH2F* vobsFVTX1Ssq[ncent][nhar][nphi];
  TH2F* vobscBBCSc[ncent][nhar][nphi];
  TH2F* vobscBBCScsq[ncent][nhar][nphi];
  TH2F* vobscFVTX1Sc[ncent][nhar][nphi];
  TH2F* vobscFVTX1Scsq[ncent][nhar][nphi];
  TH2F* vobsFVTX1trkS[ncent][nhar][nphi];
  TH2F* vobsFVTX1trkSsq[ncent][nhar][nphi];
  TH2F* vobsFVTX2trkS[ncent][nhar][nphi];
  TH2F* vobsFVTX2trkSsq[ncent][nhar][nphi];
  */
  TH2F* hfvtx[4];
  TH2F* hbbc;
  TH2F* hfvtxtrk;
  TH2F* hbbcsnfvtxtrk;
  TH1F* hnfvtxtrk;
};

#endif /* __EPANARUN16TREE_H__ */
