
#ifndef __QCUMULANTONLINE_H__
#define __QCUMULANTONLINE_H__

#include <vector>
#include <list>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include <TVectorD.h>
#include <TStopwatch.h>
#include "TLorentzVector.h"
#include "RpPar.h"

// To force loading of libraries
class TF1;
class TH1D;
class TProfile;
class TFile;
class TCanvas;

class QCumulantOnline
{
 public:

  QCumulantOnline(std::vector<TString> input,const char* output="rpanase.root");
  virtual ~QCumulantOnline();

  int Init();
  int Inittree();
  int process_event();
  int End();
  int SetRun(int _run);
  int SetcalFlag(int _flag);
  int Setweight(bool _isweight);

 private:
  std::string OutputFileName;
  std::vector<TString> InputFileName;
  int ievent;
  int jevent;
  int calFlag;
  bool isweight;
  int RunNumber;
  TFile *d_infile;
  TFile *d_outfile;
  TChain *tree;
  //tree variables
  float        event;
  float        d_bbcz;    // bbcz
  float        centrality; // integer but stored as float in PHGlobal etc
  float        bbc_qn;
  float        bbc_qs;
  unsigned int trigger_scaled;
  unsigned int trigger_live;
  float        bc_x;
  float        bc_y;
  float        vtx_z;
  float        eventfvtx_x;
  float        eventfvtx_y;
  float        eventfvtx_z;
  float        d_Qx[9];
  float        d_Qy[9];
  float        d_Qw[9];
  float        d_BBC_charge[64];

  int          npc1;
  int          d_nFVTX_clus;
  int          d_nFVTXN_clus;
  int          d_nFVTXS_clus;
  float        d_FVTX_x[10000];
  float        d_FVTX_y[10000];
  float        d_FVTX_z[10000];

  int          d_ntrk;
  float        d_px[1000];
  float        d_py[1000];
  float        d_pz[1000];
  float        d_pc3sdphi[1000];
  float        d_pc3sdz[1000];

  // --- Some histograms ---------------------------------------------
  TH1D* hd[ncent][nsub][nhar][ncorr]; 
  TH1D* hc[ncent][nsub][nhar][ncorr];
  TH1D* hv[ncent][nsub][nhar][ncorr];
  
  TH1D* hwd[ncent][nsub][nhar][ncorr];
  TH1D* hwc[ncent][nsub][nhar][ncorr];
  TH1D* hwv[ncent][nsub][nhar][ncorr];

  TH1D* hcos1[ncent][nsub][nhar]; 
  TH1D* hsin1[ncent][nsub][nhar];
  TH1D* hcos1p2[ncent][nsub][nhar]; 
  TH1D* hsin1p2[ncent][nsub][nhar];
  TH1D* hcos1m2m3[ncent][nsub][nhar]; 
  TH1D* hsin1m2m3[ncent][nsub][nhar];
  
  TH1D* hwcos1[ncent][nsub][nhar];
  TH1D* hwsin1[ncent][nsub][nhar];
  TH1D* hwcos1p2[ncent][nsub][nhar]; 
  TH1D* hwsin1p2[ncent][nsub][nhar];
  TH1D* hwcos1m2m3[ncent][nsub][nhar]; 
  TH1D* hwsin1m2m3[ncent][nsub][nhar];

  TH1D* hdpr[ncent][nsub][nhar][npt][ncorr];
  TH1D* hcpr[ncent][nsub][nhar][npt][ncorr];
  TH1D* hvpr[ncent][nsub][nhar][npt][ncorr];
  
  TH1D* hwdpr[ncent][nsub][nhar][npt][ncorr];
  TH1D* hwcpr[ncent][nsub][nhar][npt][ncorr];
  TH1D* hwvpr[ncent][nsub][nhar][npt][ncorr];
  
  TH1D* hcos1pr[ncent][nsub][nhar][npt]; 
  TH1D* hsin1pr[ncent][nsub][nhar][npt];
  TH1D* hcos1p2pr[ncent][nsub][nhar][npt]; 
  TH1D* hsin1p2pr[ncent][nsub][nhar][npt];
  TH1D* hcos1p2m3pr[ncent][nsub][nhar][npt]; 
  TH1D* hsin1p2m3pr[ncent][nsub][nhar][npt];
  TH1D* hcos1m2m3pr[ncent][nsub][nhar][npt]; 
  TH1D* hsin1m2m3pr[ncent][nsub][nhar][npt];
  
  TH1D* hwcos1pr[ncent][nsub][nhar][npt];
  TH1D* hwsin1pr[ncent][nsub][nhar][npt];
  TH1D* hwcos1p2pr[ncent][nsub][nhar][npt]; 
  TH1D* hwsin1p2pr[ncent][nsub][nhar][npt];
  TH1D* hwcos1p2m3pr[ncent][nsub][nhar][npt]; 
  TH1D* hwsin1p2m3pr[ncent][nsub][nhar][npt];
  TH1D* hwcos1m2m3pr[ncent][nsub][nhar][npt]; 
  TH1D* hwsin1m2m3pr[ncent][nsub][nhar][npt];
  
  TH1D* phiweight[ncent][nbbcz][nhar][nsub];
};
#endif /* __QCUMULANTONLINE_H__ */
