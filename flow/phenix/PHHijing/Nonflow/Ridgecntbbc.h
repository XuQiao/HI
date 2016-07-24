#ifndef __RIDGECNTBBC_H__
#define __RIDGECNTBBC_H__

#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 

class Ridgecntbbc
{
 public:

  Ridgecntbbc(TString input = "", TString output="rpanase.root");
  ~Ridgecntbbc();

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
    int n;

//particle variables
    float px[10000];
    float py[10000];
    float pz[10000];
    float E[10000];
    int KS[10000];
    int KF[10000];
    float M[10000];

//----------------------------------------------------------------------------
  //BBC
  TH1F* hnch;
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

#endif /* __RIDGECNTBBC_H__ */
