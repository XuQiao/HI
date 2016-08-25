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

  mm(TString input = "", TString input1 = "", TString output="rpanase.root");
  ~mm();

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
  int nfvtxtrack;
  float fvtxchi2[1000];
};

#endif /* __mm_H__ */
