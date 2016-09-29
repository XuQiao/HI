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

  std::vector<std::vector<TH1F*> > sinflt[ncent][nbbcz][nhar][nsub][nord];
  std::vector<std::vector<TH1F*> > cosflt[ncent][nbbcz][nhar][nsub][nord];
  
//  std::vector<std::vector<TH1F*> > sinflttmp[ncent][nbbcz][nhar][nsub][nord];
//  std::vector<std::vector<TH1F*> > cosflttmp[ncent][nbbcz][nhar][nsub][nord];

  std::vector<std::vector<TH1F*> > q[ncent][nbbcz][nhar][nsub][nxy];
//  std::vector<std::vector<TH1F*> > qtmp[ncent][nbbcz][nhar][nsub][nxy];  //invidual calibration
  std::vector<std::vector<TH1F*> > qRec[ncent][nbbcz][nhar][nsub][nxy];
  std::vector<std::vector<TH1F*> > psi[ncent][nbbcz][nhar][nsub];
//  std::vector<std::vector<TH1F*> > psitmp[ncent][nbbcz][nhar][nsub];
  std::vector<std::vector<TH1F*> > psiFla[ncent][nbbcz][nhar][nsub];

   std::vector<std::vector<TH1F*> > phiweight[ncent][nbbcz][nhar][nsub];
   std::vector<std::vector<TProfile*> > phiweightbbc[ncent][nbbcz][nhar][nsub];

  float meanx[ncent][nbbcz][nhar][nsub][nangle1][nagle2];
  float meany[ncent][nbbcz][nhar][nsub][nangle1][nagle2];
  float rmsx[ncent][nbbcz][nhar][nsub][nangle1][nagle2];
  float rmsy[ncent][nbbcz][nhar][nsub][nangle1][nagle2];

  float cosfltarr[ncent][nbbcz][nhar][nsub][nord][nangle1][nagle2];
  float sinfltarr[ncent][nbbcz][nhar][nsub][nord][nangle1][nagle2];

  /*
  TH1F* EPRCNTBBCS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTFVTX1S[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCSFVTX1S[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTFVTX1S[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTFVTX2S[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTFVTX1LS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTFVTX2LS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTFVTX3LS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTFVTX4LS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCSFVTX2S[nangle1][nangle2][ncent][nhar];
  */
  TH1F* EPRBBCSFVTX1S[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCSFVTX1LS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCSFVTX2LS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCSFVTX3LS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCSFVTX4LS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCSFVTX1p2p3LS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCSFVTXtrkS[nangle1][nangle2][ncent][nhar];
  /*
  TH1F* EPRCNTFVTX1trkS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTFVTX2trkS[nangle1][nangle2][ncent][nhar];

  TH1F* EPRCNTcBBCSc[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTcFVTX1S[nangle1][nangle2][ncent][nhar];

  TH1F* EPRBBCSFVTX1trkS[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCSFVTX2trkS[nangle1][nangle2][ncent][nhar];
          
  TH1F* EPRBBCScFVTX1S[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTcBBCSc[nangle1][nangle2][ncent][nhar];
  TH1F* EPRCNTcFVTX1Sc[nangle1][nangle2][ncent][nhar];
  TH1F* EPRBBCScFVTX1Sc[nangle1][nangle2][ncent][nhar];
  */
/*
  TH2F* vobsFVTX2S[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX2Ssq[nangle1][nangle2][ncent][nhar][nphi];
 */ 

  TH2F* vobsBBCS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX1S[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX1LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX2LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX3LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX4LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX1p2p3LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTXtrkS[nangle1][nangle2][ncent][nhar][nphi];
  
  TH2F* vBBCS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vFVTX1S[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vFVTX1LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vFVTX2LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vFVTX3LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vFVTX4LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vFVTX1p2p3LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vFVTXtrkS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vnBBCS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vnFVTX1S[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vnFVTX1LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vnFVTX2LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vnFVTX3LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vnFVTX4LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vnFVTX1p2p3LS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vnFVTXtrkS[nangle1][nangle2][ncent][nhar][nphi];
 /* 
  TH2F* vobsBBCS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsBBCSsq[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX1S[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX1Ssq[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobscBBCSc[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobscBBCScsq[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobscFVTX1Sc[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobscFVTX1Scsq[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX1trkS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX1trkSsq[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX2trkS[nangle1][nangle2][ncent][nhar][nphi];
  TH2F* vobsFVTX2trkSsq[nangle1][nangle2][ncent][nhar][nphi];
  */
  TH2F* hfvtx[4];
  TH2F* hbbc;
  TH2F* hfvtxtrk;
  TH2F* hbbcsnfvtxtrk;
  TH1F* hnfvtxtrk;
};

#endif /* __EPANARUN16TREE_H__ */
