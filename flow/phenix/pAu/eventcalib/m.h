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
  int cent;
  int nbbc;
  int nbbc1;
  float bbct0[128];
};

#endif /* __M_H__ */
