#ifndef __M_H__
#define __M_H__

#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TProfile.h" 

class m
{
 public:

  m(TString input = "", TString input1 = "",  TString output="rpanase.root");
  ~m();

  int Init();
  int process_event();
  int End();
  
 private:
  TString OutputFileName;
  TString InputFileName;
  TString InputFileName1;
  TString coll;
  TFile *d_infile;
  TFile *d_infile1;
  TFile *d_outfile;

  TTree *tree;
  TTree *tree1;
  TTree *newtree;
  float bbcv;
  float bbcv1;
  float vtxz;
  float fvtxz;
  int nfvtxtrack;
  int farm[1000];
  int fnhits[1000];
  float fthe[1000];
  float feta[1000];
  float fphi[1000];
  float fvtxX[1000];
  float fvtxY[1000];
  float fvtxZ[1000];
};

#endif /* __M_H__ */
