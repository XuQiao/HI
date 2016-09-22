#include "PerformTestMB.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include "TString.h"
#include "TVectorD.h"
#include <fstream>
#include <string>
#include <TROOT.h>
#include "TSystem.h"
#include <iostream>

//#include "Run16dAupc3dphidzcalibsmoothpass1.h"
#include "func.h"

using namespace std;
int get_fvtx_layer(float);
void initialize_pmt_position();
float d_pmt_z = -1443.5; // same for all tubes

const double mpion = 0.139570;
const double phbeta = 29.9792458;

const int centbin = 1;
const int ntrks = 40;
const int nbbcs = 800;
const int nvtxs = 800;
const int nfvtxs = 800;
const int nmpcs = 400;
const int ncluss = 64;
const float mpc_e_cut = 10.0;

const int ntry = 1;
const int nwidth = 10;
//const float bbct_sigma[ntry] = {0.5,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,-1};
const float bbct_sigma[ntry] = {-1};
const float bbct_width[nwidth] = {0.1,0.2,0.3,0.4,0.5,0.7,0.8,0.9,1.0,1.2};

const int nbuff = 10;
float buff_cent[nbuff] = {};
float buff_bbcs[nbuff] = {};
float buff_ZdcEs[nbuff] = {};
float buff_bbct0meansouth[nbuff] = {};
float buff_bbct0meannorth[nbuff] = {};
float buff_bbct0sqsouth[nbuff] = {};
float buff_bbct0sqnorth[nbuff] = {};
float buff_nbbcau[nbuff] = {};
float buff_nbbcde[nbuff] = {};
TH1F buff_bbctsouthtmp[nbuff];
TH1F buff_bbctnorthtmp[nbuff];
vector<float> buff_bbctau[nbuff];
vector<float> buff_bbctde[nbuff];
int ibuff=0;

bool isgood(TH1F h, float sigma);
float isgood1(vector<float> bbct0, float bbct0mean, int nbbc_, float sigma);

//_____________________________________________________________________________________________________________________________
PerformTestMB::PerformTestMB(std::vector<TString> input, std::vector<TString> input1, const char* output) :
  InputFileName(input), InputFileName(input1), OutputFileName(output), ievent(0), jevent(0), RunNumber(0)
{
  d_outfile=NULL;
  htrig=NULL;
//Tracks
  hcntetaphi = NULL;
  hcntpt = NULL;
for(int i=0;i<50;i++){
pc3dphidz_arm0_pos[i] = NULL;
pc3dphidz_arm1_pos[i] = NULL;
pc3dphidz_arm0_neg[i] = NULL;
pc3dphidz_arm1_neg[i] = NULL;
}
for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
    for(int i=0;i<50;i++){
        pc3dphidz_arm0_pos_z[ibbcz][i]=NULL;
        pc3dphidz_arm1_pos_z[ibbcz][i]=NULL;
        pc3dphidz_arm0_neg_z[ibbcz][i]=NULL;
        pc3dphidz_arm1_neg_z[ibbcz][i]=NULL;
    }
}


//tof
for(int ich=0;ich<2;ich++){
  tofdphidz[ich] = NULL;
  tofwdphidz[ich] = NULL;
}
  tofsdphisdz = NULL;
  tofwsdphisdz = NULL;
  ttofqpratio = NULL;
  m2qpratio = NULL;
  ttofp = NULL;
  m2p = NULL;
  pinv2chbeta = NULL;
  deltattofeis = NULL;

//vtx
  for (int ilayer=0; ilayer<4; ilayer++) {
  hcluetaphi[ilayer]  = NULL;
}
//bbc
  bbcet = NULL;

//fvtx
  for (int iarm=0; iarm<2; iarm++) {
  DCAxydis[iarm] = NULL;
  DCAxy2dis[iarm] = NULL;
  DCAcentdis[iarm] = NULL;
}
  fvtxdphidis = NULL;
  fvtxdphidis2 = NULL;

  hvtx0etaz = NULL;
  hvtx1etaz = NULL;

  hvtx0etaphi = NULL;
  hvtx1etaphi = NULL;

//mpc
 
mpcetdis  = NULL;
for (int i=0; i<centbin; i++) {
mpcetetasouth[i] = NULL;
mpcetetanorth[i] = NULL;
}
 mpc_south_cent = NULL;
 mpc_north_cent = NULL;
 mpc_south_north = NULL;
south_mpc_north_mpc = NULL;


//correlate
  hvtxzfvtxz = NULL;
  hvtxxvtxz = NULL;
  hvtxyvtxz = NULL;
  hpc1hitsbbc = NULL;
  hnpc3hitsntof = NULL;
  hbbcnbbc = NULL;
  hbbcsbbcn = NULL;
  for (int ilayer=0; ilayer<4; ilayer++) {
  hnvtxnfvtxtrk[ilayer] = NULL;
  hnvtxnmpc[ilayer] = NULL;
  hbbcsnvtx[ilayer] = NULL;
  hbbcnnvtx[ilayer] = NULL;
  hbbcnvtx[ilayer] = NULL;
}
  hnbbcnclu = NULL;
  hntracknmpc = NULL;
  hnfvtxtrkbbc = NULL;
  hnfvtxtrksnmpcs = NULL;
  hnfvtxtrknnmpcn = NULL;
  south_mpc_south_bbc = NULL;
  north_mpc_north_bbc = NULL;
for(int itry=0;itry<ntry;itry++){
  hcentbbct0sigmasouth[itry]=NULL;
  hcentbbct0sigmanorth[itry]=NULL;
  hcentbbcmixt0sigmasouth[itry]=NULL;
  hcentbbcmixt0sigmanorth[itry]=NULL;
}
for(int itry=0;itry<nwidth;itry++){
  hcentbbct0fracsouth[itry]=NULL;
  hcentbbct0fracnorth[itry]=NULL;
  hcentbbcmixt0fracsouth[itry]=NULL;
  hcentbbcmixt0fracnorth[itry]=NULL;
}
  
  pi=0.0;

}

//_____________________________________________________________________________________________________________________________
PerformTestMB::~PerformTestMB()
{
  cout << " PerformTestMB::~PerformTestMB " << endl;
	}

//_____________________________________________________________________________________________________________________________
int PerformTestMB::Init()
{
TH1::SetDefaultSumw2(); 
  cout << " PerformTestMB::Init " << endl;
  char name[200];
  pi=acos(-1.0);
  
  d_outfile = new TFile(OutputFileName.c_str(),"recreate");
  sprintf(name,"htrig");    htrig = new TH1F(name,name,10000,0,10000);
//Tracks
sprintf(name,"hcntetaphi"); hcntetaphi = new TH2F(name,name,200,-1,1,200,-2*pi,2*pi);
sprintf(name,"hcntpt"); hcntpt = new TH1F(name,name,200,0,10);
for(int i=0;i<50;i++){
sprintf(name,"pc3dphidz_arm0_pos_%d",i),pc3dphidz_arm0_pos[i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
sprintf(name,"pc3dphidz_arm1_pos_%d",i),pc3dphidz_arm1_pos[i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
sprintf(name,"pc3dphidz_arm0_neg_%d",i),pc3dphidz_arm0_neg[i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
sprintf(name,"pc3dphidz_arm1_neg_%d",i),pc3dphidz_arm1_neg[i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
}
for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
for(int i=0;i<50;i++){
sprintf(name,"pc3dphidz_arm0_pos_z%d_%d",ibbcz,i),pc3dphidz_arm0_pos_z[ibbcz][i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
sprintf(name,"pc3dphidz_arm1_pos_z%d_%d",ibbcz,i),pc3dphidz_arm1_pos_z[ibbcz][i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
sprintf(name,"pc3dphidz_arm0_neg_z%d_%d",ibbcz,i),pc3dphidz_arm0_neg_z[ibbcz][i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
sprintf(name,"pc3dphidz_arm1_neg_z%d_%d",ibbcz,i),pc3dphidz_arm1_neg_z[ibbcz][i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
}
}


//tof
for(int ich=0;ich<2;ich++){
sprintf(name,"tofdphidz_%d",ich);  tofdphidz[ich] = new TH2F(name,name,200,-0.2,0.2,200,-10,10);
sprintf(name,"tofwdphidz_%d",ich);  tofwdphidz[ich] = new TH2F(name,name,200,-0.2,0.2,200,-10,10);
}
sprintf(name,"tofsdphisdz");  tofsdphisdz = new TH2F(name,name,200,-10,10,200,-10,10);
sprintf(name,"tofwsdphisdz");  tofwsdphisdz = new TH2F(name,name,200,-10,10,200,-10,10);
sprintf(name,"ttofqpratio");  ttofqpratio = new TH2F(name,name,500,0,100,500,-8,8);
sprintf(name,"m2qpratio");  m2qpratio = new TH2F(name,name,500,0,2,500,-8,8);
sprintf(name,"ttofp");  ttofp = new TH2F(name,name,500,0,100,500,0,10);
sprintf(name,"pinv2chbeta"); pinv2chbeta = new TH2F(name,name,500,0,10,500,-0.8,0.8);
sprintf(name,"m2p");  m2p = new TH2F(name,name,500,0,2,500,0,10);
sprintf(name,"deltattofeis"); deltattofeis = new TH2F(name,name,960,0,960,1600,-20,60);

//vtx
for(int i=0;i<4;i++){
sprintf(name,"hcluetaphi_%d",i); hcluetaphi[i] = new TH2F(name,name,300,-3.0,3.0,100,-4,4);
}

//bbc
sprintf(name,"bbcet"); bbcet = new TH1F(name,name,500,0,20);
sprintf(name,"bbctdevsouth"); bbctdevsouth = new TH1F(name,name,2000,-10,10);
sprintf(name,"bbctdevnorth"); bbctdevnorth = new TH1F(name,name,2000,-10,10);
sprintf(name,"bbctsouth"); bbctsouth = new TH1F(name,name,2000,0,20);
sprintf(name,"bbctnorth"); bbctnorth = new TH1F(name,name,2000,0,20);

//fvtx
  for(int iarm=0; iarm<2; iarm++){
    sprintf(name,"DCAxydis_%d",iarm); DCAxydis[iarm] = new TH2D(name,name, 100, -2, 2, 100, -2, 2);
    sprintf(name,"DCAxy2dis_%d",iarm); DCAxy2dis[iarm] = new TH2D(name,name, 100, -2, 2, 100, -2, 2);
    sprintf(name,"DCAcentdis_%d",iarm); DCAcentdis[iarm] = new TH2D(name,name, 10, 0, 10, 100, 0, 10);
  }
 // sprintf(name,"fvtxdphidis"); fvtxdphidis = new TH1F(name, name, 100, -10.0, 10.0);
 // sprintf(name,"fvtxdphidis2"); fvtxdphidis2 = new TH1F(name, name, 100, -10.0, 10.0);
  sprintf(name,"hvtx0etaz"); hvtx0etaz = new TH2D(name,name,240, -12, 12, 600, -3.0, 3.0);
  sprintf(name,"hvtx1etaz"); hvtx1etaz = new TH2D(name,name,240, -12, 12, 600, -3.0, 3.0);
  sprintf(name,"hvtx0etaphi"); hvtx0etaphi = new TH2D(name,name, 240, -3.14*1.1, 3.14*1.1, 600, -3.0, 3.0);
  sprintf(name,"hvtx1etaphi"); hvtx1etaphi = new TH2D(name,name, 240, -3.14*1.1, 3.14*1.1, 600, -3.0, 3.0);

//mpc
  sprintf(name,"mpcetdis"); mpcetdis = new TH2F(name,name,600,-0.5,599.5,60,0.0,6.0);
  sprintf(name,"mpc_south_cent"); mpc_south_cent = new TH2F(name, name, 100, -0.5, 99.5, 250, -0.5, 249.5);
  sprintf(name,"mpc_north_cent"); mpc_north_cent = new TH2F(name, name, 100, -0.5, 99.5, 250, -0.5, 249.5);
  sprintf(name,"mpc_south_north"); mpc_south_north = new TH2F(name, name, 600, -0.5, 599.5, 600, -0.5, 599.5);
  sprintf(name,"south_mpc_north_mpc");  south_mpc_north_mpc = new TH2F(name, name, 600, -0.5, 19.5, 600, -0.5, 19.5);
for (int i=0; i<centbin; i++) {
 sprintf(name,"mpcetetasouth_%d", i);    mpcetetasouth[i] = new TH2F(name,name,30, 3.0, 4.0, 60,0.0,6.0);
 sprintf(name,"mpcetetanorth_%d", i);    mpcetetanorth[i] = new TH2F(name,name,30, 3.0, 4.0, 60,0.0,6.0);
}

//correlate
  sprintf(name,"hbbcsZdcEs"); hbbcsZdcEs = new TH2F(name,name,600,0,300,500,0,5000);
  sprintf(name,"hmixbbcsZdcEs"); hmixbbcsZdcEs = new TH2F(name,name,600,0,600,500,0,10000);
  sprintf(name,"hnpc3hitsntof"); hnpc3hitsntof = new TH2F(name,name,50,0,50,50,0,50);
  sprintf(name,"hbbcnbbc"); hbbcnbbc = new TH2F(name,name,600,0,200,100,0,100);
  sprintf(name,"hpc1hitsbbc"); hpc1hitsbbc = new TH2F(name, name, 20, 0, 20, 600, 0, 200);
  sprintf(name,"hvtxzfvtxz"); hvtxzfvtxz = new TH2F(name,name,100,-25.,25.,100,-25.,25.);
  sprintf(name,"hvtxxvtxz"); hvtxxvtxz = new TH2F(name,name,20,-2.,2.,100,-10.,10.);
  sprintf(name,"hvtxyvtxz"); hvtxyvtxz = new TH2F(name,name,20,-2.,2.,100,-10.,10.);
  for(int i=0;i<4;i++){
	sprintf(name,"hbbcsnvtx_%d",i);  hbbcsnvtx[i] = new TH2F(name,name,600,0,200,500,0,500);
	sprintf(name,"hbbcnnvtx_%d",i);  hbbcnnvtx[i] = new TH2F(name,name,600,0,200,500,0,500);
	sprintf(name,"hbbcnvtx_%d",i);  hbbcnvtx[i] = new TH2F(name,name,600,0,200,500,0,500);
        sprintf(name,"hnvtxnmpc_%d",i); hnvtxnmpc[i] = new TH2F(name,name,500,0,500,100,0,100);
  	sprintf(name,"hnvtxnfvtxtrk_%d",i); hnvtxnfvtxtrk[i] = new TH2F(name,name,500,0,500,100,0,100);
  }
  sprintf(name,"hnbbcnclu"); hnbbcnclu = new TH2F(name,name,100,0,100,50,0,50);
  sprintf(name,"hbbcsbbcn"); hbbcsbbcn = new TH2F(name,name,600,0,200,600,0,200);
  sprintf(name,"hntracknmpc"); hntracknmpc = new TH2F(name,name,50,0,50,100,0,100);
  sprintf(name,"hnfvtxtrkbbc"); hnfvtxtrkbbc = new TH2F(name,name,50,0,50,600,0,200);
  sprintf(name,"hnfvtxtrksnmpcs"); hnfvtxtrksnmpcs = new TH2F(name,name,50,0,50,50,0,50);
  sprintf(name,"hnfvtxtrksnmpcn"); hnfvtxtrknnmpcn = new TH2F(name,name,50,0,50,50,0,50);
  sprintf(name,"south_mpc_south_bbc");  south_mpc_south_bbc = new TH2F(name, name, 600, -0.5, 19.5, 600, -0.5, 199.5);
  sprintf(name,"north_mpc_north_bbc");  north_mpc_north_bbc = new TH2F(name, name, 600, -0.5, 19.5, 600, -0.5, 199.5);

for(int itry=0;itry<ntry;itry++){
  sprintf(name,"hcentbbct0sigmasouth_%d",itry); hcentbbct0sigmasouth[itry] = new TH2F(name,name,600,0,300,500,0,10);
  sprintf(name,"hcentbbct0sigmanorth_%d",itry); hcentbbct0sigmanorth[itry] = new TH2F(name,name,600,0,300,500,0,10);
  sprintf(name,"hcentbbcmixt0sigmasouth_%d",itry); hcentbbcmixt0sigmasouth[itry] = new TH2F(name,name,600,0,300,500,0,10);
  sprintf(name,"hcentbbcmixt0sigmanorth_%d",itry); hcentbbcmixt0sigmanorth[itry] = new TH2F(name,name,600,0,300,500,0,10);
}
for(int itry=0;itry<nwidth;itry++){
  sprintf(name,"hcent1cent2_%d",itry); hcent1cent2[itry] = new TH2F(name,name,100,0,100,100,0,100);
  sprintf(name,"hcentbbct0fracsouth_%d",itry); hcentbbct0fracsouth[itry] = new TH2F(name,name,600,0,300,110,0,1.1);
  sprintf(name,"hcentbbct0fracnorth_%d",itry); hcentbbct0fracnorth[itry] = new TH2F(name,name,600,0,300,110,0,1.1);
  sprintf(name,"hcentbbcmixt0fracsouth_%d",itry); hcentbbcmixt0fracsouth[itry] = new TH2F(name,name,600,0,300,110,0,1.1);
  sprintf(name,"hcentbbcmixt0fracnorth_%d",itry); hcentbbcmixt0fracnorth[itry] = new TH2F(name,name,600,0,300,110,0,1.1);
}

  cout<<"finish of initialize"<<endl;
  return 0;
}

int PerformTestMB::Inittree(){
  //------------------------------------------------------------//
  //               Initializing Tree Variables                  //
  //------------------------------------------------------------//
  tree = NULL;
  tree1 = NULL;
  tree = new TChain("tree");
  tree1 = new TChain("tree");
  for(unsigned int itree=0;itree<InputFileName.size();itree++){
     cout<<InputFileName[itree]<<endl;
     try{
         tree->Add(InputFileName[itree]);}
     catch(int e){ cout << "An exception occurred. Exception Nr. " << e << '\n';}
  }
  for(unsigned int itree=0;itree<InputFileName1.size();itree++){
     cout<<InputFileName1[itree]<<endl;
     try{
         tree1->Add(InputFileName1[itree]);}
     catch(int e){ cout << "An exception occurred. Exception Nr. " << e << '\n';}
  }
  cout << "Now getting ready to read in the tree branch addresses and stuff...." << endl;

  // List of branches
  TBranch* b_run;
//  TBranch* b_event;   //!
  TBranch* b_bbc_z;   //!
  TBranch* b_centrality;   //!
  TBranch* b_bbc_qn;   //!
  TBranch* b_bbc_qs;   //!
  TBranch* b_zdc_qn;   //!
  TBranch* b_zdc_qs;   //!
  TBranch* b_npc1;   //!
  TBranch* b_trigger_scaled;   //!
//  TBranch* b_trigger_live;   //!
//  TBranch* b_d_Qx;   //!
//  TBranch* b_d_Qy;   //!
//  TBranch* b_d_Qw;   //!
//  TBranch* b_bc_x;   //!
//  TBranch* b_bc_y;   //!
//  TBranch* b_vtx_z;   //!
  TBranch* b_fvtx_x;   //!
  TBranch* b_fvtx_y;   //!
  TBranch* b_fvtx_z;   //!
  TBranch* b_d_BBC_charge;   //!
  TBranch* b_d_BBC_time0;   //!
  TBranch* b_d_BBC_valid;   //!
  TBranch* b_d_nFVTX_clus;   //!
//  TBranch* b_d_nFVTXN_clus;   //!
//  TBranch* b_d_nFVTXS_clus;   //!
  TBranch* b_d_FVTX_x;   //!
  TBranch* b_d_FVTX_y;   //!
  TBranch* b_d_FVTX_z;   //!
  TBranch* b_ntrk;   //!
  TBranch* b_mom;   //!
  TBranch* b_phi0;   //!
  TBranch* b_the0;   //!
  TBranch* b_charge;   //!
  TBranch* b_pc3dphi;   //!
  TBranch* b_pc3dz;   //!
  // TBranch* b_d_ntrk;   //!
  // TBranch* b_d_cntpx;   //!
  // TBranch* b_d_cntpy;   //!
  // TBranch* b_d_cntpz;   //!
  TBranch* b_d_nfvtxtrk;   //!
  TBranch* b_d_nfvtxtrk1;  //!
  TBranch* b_fvtxchi;
  TBranch* b_farm;
  TBranch* b_fnhits;
  TBranch* b_feta;
  TBranch* b_fphi;
  TBranch* b_fvtxX;
  TBranch* b_fvtxY;
  TBranch* b_fvtxZ;


  tree->SetBranchAddress("run",&RunNumber,&b_run);
  tree->SetBranchAddress("bbcv",&d_bbcz,&b_bbc_z);
  tree->SetBranchAddress("cent",&centrality,&b_centrality);
  tree->SetBranchAddress("bbc_s",&bbc_qs,&b_bbc_qs);
  tree->SetBranchAddress("bbc_n",&bbc_qn,&b_bbc_qn);
  tree->SetBranchAddress("ZdcEs",&ZdcEs,&b_zdc_qs);
  tree->SetBranchAddress("ZdcEn",&ZdcEn,&b_zdc_qn);
  tree->SetBranchAddress("npc1hits",&npc1,&b_npc1);
  tree->SetBranchAddress("trig",&trigger_scaled,&b_trigger_scaled);

  tree->SetBranchAddress("fvtxx",&eventfvtx_x,&b_fvtx_x);
  tree->SetBranchAddress("fvtxy",&eventfvtx_y,&b_fvtx_y);
  tree->SetBranchAddress("fvtxz",&eventfvtx_z,&b_fvtx_z);

  tree->SetBranchAddress("bbccharge",d_BBC_charge,&b_d_BBC_charge);
  tree->SetBranchAddress("bbct0",d_BBC_time0,&b_d_BBC_time0);
  tree1->SetBranchAddress("bbcvalid",d_BBC_valid,&b_d_BBC_valid);

  tree->SetBranchAddress("nclus",&d_nFVTX_clus,&b_d_nFVTX_clus);
  tree->SetBranchAddress("fclusX",d_FVTX_x,&b_d_FVTX_x);
  tree->SetBranchAddress("fclusY",d_FVTX_y,&b_d_FVTX_y);
  tree->SetBranchAddress("fclusZ",d_FVTX_z,&b_d_FVTX_z);

  tree->SetBranchAddress("ntrack",&d_ntrk,&b_ntrk);
  tree->SetBranchAddress("mom",d_mom,&b_mom);
  tree->SetBranchAddress("phi0",d_phi0,&b_phi0);
  tree->SetBranchAddress("the0",d_the0,&b_the0);
  tree->SetBranchAddress("charge",d_charge,&b_charge);
  tree->SetBranchAddress("pc3dphi",d_pc3dphi,&b_pc3dphi);
  tree->SetBranchAddress("pc3dz",d_pc3dz,&b_pc3dz);

  tree->SetBranchAddress("nfvtxtrack",&d_nfvtxtrk,&b_d_nfvtxtrk);
  tree1->SetBranchAddress("nfvtxtrack",&d_nfvtxtrk1,&b_d_nfvtxtrk1);
  tree1->SetBranchAddress("fvtxchi2",d_fvtxchi,&b_fvtxchi);
  tree->SetBranchAddress("farm",d_farm,&b_farm);
  tree->SetBranchAddress("fnhits",d_fnhits,&b_fnhits);
  tree->SetBranchAddress("feta",d_feta,&b_feta);
  tree->SetBranchAddress("fphi",d_fphi,&b_fphi);
  tree->SetBranchAddress("fvtxX",d_fvtxX,&b_fvtxX);
  tree->SetBranchAddress("fvtxY",d_fvtxY,&b_fvtxY);
  tree->SetBranchAddress("fvtxZ",d_fvtxZ,&b_fvtxZ);
  return 0;
}

//_____________________________________________________________________________________________________________________________
int PerformTestMB::process_event()
{
  TH1F bbctsouthtmp(Form("bbctsouthtmp_%d",ievent),"bbctsouthtmp",401,0.05,20.05);
  TH1F bbctnorthtmp(Form("bbctnorthtmp_%d",ievent),"bbctnorthtmp",401,0.05,20.05);
  
  int nEvent = tree->GetEntries();
  cout<<nEvent<<endl;
  for(ievent=0;ievent < nEvent; ievent++){
      tree->GetEntry(ievent);
      tree1->GetEntry(ievent);

  if(ievent%10000==0) {
    cout<<"************* ievent= "<<ievent<<"    *************"<<endl;
  }

  //global
    unsigned int trigger_FVTXNSBBCScentral = 0x00100000;
    unsigned int trigger_FVTXNSBBCS        = 0x00400000;
    unsigned int trigger_BBCLL1narrowcent  = 0x00000008;
    unsigned int trigger_BBCLL1narrow      = 0x00000010;
    unsigned int accepted_triggers = 0;
      // --- Run16dAu200                                                                                                    
   if ( RunNumber >= 454774 && RunNumber <= 455639 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
      // --- Run16dAu62                                                                                                                     
   if ( RunNumber >= 455792 && RunNumber <= 456283 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
      // --- Run16dAu20                                                                   
   if ( RunNumber >= 456652 && RunNumber <= 457298 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;
      // --- Run16dAu39                                                                        
   if ( RunNumber >= 457634 && RunNumber <= 458167 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;

  unsigned int passes_trigger = trigger_scaled & accepted_triggers;
  if ( passes_trigger == 0 )      continue;
  float bbcv = d_bbcz;
  int cent = centrality;
//  float bbc_s = bbc_qs;
  float vvertex = bbcv;

  
//  float fvtxx = eventfvtx_x;
//  float fvtxy = eventfvtx_y;
  float fvtxz = eventfvtx_z;
  if(fvtxz != fvtxz){
      fvtxz = -9999;
  }
  
  if(RunNumber>=455792 && RunNumber<=458167){ // low energies use fvtx as vertex if there is one
  if(fvtxz != -9999)
    vvertex = fvtxz;
  }

 int ibbcz = -1;
 ibbcz = nbbcz*(vvertex+10)/20;
 if(ibbcz<0||ibbcz>=nbbcz) continue;

 int icent = -9999;
 if(cent<0) continue;
 icent = 0;
 if(icent<0) continue;

 // --- all numbers from Darren 2016-06-23
      const float x_off = 0.3;
      const float beam_angle = 0.001;
      float vtx_z = vvertex;
      float vtx_x = x_off + atan(beam_angle)*vtx_z;
      float vtx_y = 0.02;

    hvtxxvtxz->Fill(eventfvtx_x,eventfvtx_z);
    hvtxyvtxz->Fill(eventfvtx_y,eventfvtx_z);
    jevent++;

 // Tracks 
  int ntrack=0;
  for(int itrk=0; itrk< d_ntrk; itrk++){
      float mom    = d_mom[itrk];
      float phi0    = d_phi0[itrk];
      float the0    = d_the0[itrk];
      int charge = d_charge[itrk];
      float pc3dphi   = d_pc3dphi[itrk];
      float pc3dz   = d_pc3dz[itrk];
      float pt        = mom*sin(the0);
      float pz        = mom*cos(the0);
      float px = pt * cos(phi0);
      float py = pt * sin(phi0);
      float eta       = -log(tan(0.5*the0));
      px = pz*sin(-beam_angle) + px*cos(-beam_angle);
      float phi = atan2(py,px);
      pt = sqrt(px*px+py*py);
      int dcarm=0;
      if(px>0) dcarm=1;
      
      hcntetaphi->Fill(eta,phi);
      hcntpt->Fill(pt);      
	for(int ipt = 0;ipt<50;ipt++){
      if(pt>=0.1*ipt && pt<0.1*(ipt+1)){
	if(charge>0 && dcarm == 0)	pc3dphidz_arm0_pos[ipt]->Fill(pc3dphi,pc3dz);
	if(charge>0 && dcarm == 1)	pc3dphidz_arm1_pos[ipt]->Fill(pc3dphi,pc3dz);
	if(charge<0 && dcarm == 0)	pc3dphidz_arm0_neg[ipt]->Fill(pc3dphi,pc3dz);
	if(charge<0 && dcarm == 1)	pc3dphidz_arm1_neg[ipt]->Fill(pc3dphi,pc3dz);
	if(charge>0 && dcarm == 0)	pc3dphidz_arm0_pos_z[ibbcz][ipt]->Fill(pc3dphi,pc3dz);
	if(charge>0 && dcarm == 1)	pc3dphidz_arm1_pos_z[ibbcz][ipt]->Fill(pc3dphi,pc3dz);
	if(charge<0 && dcarm == 0)	pc3dphidz_arm0_neg_z[ibbcz][ipt]->Fill(pc3dphi,pc3dz);
	if(charge<0 && dcarm == 1)	pc3dphidz_arm1_neg_z[ibbcz][ipt]->Fill(pc3dphi,pc3dz);
	}
      }
      if(pt>0.2&&pt<5.0){
	ntrack++;
      }

    }
  
  // beam beam counter (bbc) r.p. // -----------------------------------
  int nbbcau=0;//bbc south
  int nbbcde=0;//bbc north
  
  float bbcde_t0[nbbcs];
  float bbcau_t0[nbbcs];

  float bbct0meansouth = 0, bbct0sqsouth = 0;
  float bbct0meannorth = 0, bbct0sqnorth = 0;

    for(int ipmt = 0; ipmt < 128; ipmt++){
    float bbcx, bbcy, bbcz;
    if(ipmt < 64){
    bbcx      = d_pmt_x[ipmt];
    bbcy      = d_pmt_y[ipmt];
    bbcz      = d_pmt_z;
    }
    else{
    bbcx      = d_pmt_x[ipmt-64];
    bbcy      = d_pmt_y[ipmt-64];
    bbcz      = -d_pmt_z;
    }
    float charge = d_BBC_charge[ipmt];
    float time0 = d_BBC_time0[ipmt];
    bbcx = bbcx - vtx_x*10.0;
    bbcy = bbcy - vtx_y*10.0;
    bbcz = bbcz - vtx_z*10.0;
    // --- rotation
    bbcx = bbcz*sin(-beam_angle) + bbcx*cos(-beam_angle);

    if (time0>0 && charge>0) {
      int iarm = 0;
      if (bbcz > 0) iarm = 1;
      float val=charge;
      
	if(iarm==0){
	  bbcau_t0[nbbcau] = time0;
          bbct0meansouth += time0;
          bbct0sqsouth += time0*time0;
	  nbbcau++;
	}
	else{
	  bbcde_t0[nbbcde] = time0;
          bbct0meannorth += time0;
          bbct0sqnorth += time0*time0;
          nbbcde++;
	}
      bbcet->Fill(val);
      }
  }
bbctsouthtmp.Reset("M");
bbctnorthtmp.Reset("M");

for(int ipmt=0;ipmt<nbbcau;ipmt++){
   //bbctsouthtmp.Fill(bbcau_t0[ipmt]-bbct0meansouth/nbbcau);
   bbctsouthtmp.Fill(bbcau_t0[ipmt]);
}
for(int ipmt=0;ipmt<nbbcde;ipmt++){
   //bbctnorthtmp.Fill(bbcde_t0[ipmt]-bbct0meannorth/nbbcde);
   bbctnorthtmp.Fill(bbcde_t0[ipmt]);
}

float bbct0sigmasouth = TMath::Sqrt(bbct0sqsouth/nbbcau - (bbct0meansouth/nbbcau) * (bbct0meansouth/nbbcau));
float bbct0sigmanorth = TMath::Sqrt(bbct0sqnorth/nbbcde - (bbct0meannorth/nbbcde) * (bbct0meannorth/nbbcde));
for(int ipmt = 0; ipmt < 128; ipmt++){
    float bbcx, bbcy, bbcz;
    if(ipmt < 64){
    bbcx      = d_pmt_x[ipmt];
    bbcy      = d_pmt_y[ipmt];
    bbcz      = d_pmt_z;
    }
    else{
    bbcx      = d_pmt_x[ipmt-64];
    bbcy      = d_pmt_y[ipmt-64];
    bbcz      = -d_pmt_z;
    }
    float charge = d_BBC_charge[ipmt];
    float time0 = d_BBC_time0[ipmt];
    bbcx = bbcx - vtx_x*10.0;
    bbcy = bbcy - vtx_y*10.0;
    bbcz = bbcz - vtx_z*10.0;
    // --- rotation
    bbcx = bbcz*sin(-beam_angle) + bbcx*cos(-beam_angle);

  if (charge>0 && time0>0) {
    int iarm = 0;
    if (bbcz > 0) iarm = 1;
      if(iarm==0){
          bbctdevsouth->Fill(time0-bbct0meansouth/nbbcau);
          bbctsouth->Fill(time0);
      }
      else{
          bbctdevnorth->Fill(time0-bbct0meannorth/nbbcde);
          bbctnorth->Fill(time0);
      }
  }
}

//FVTX tracks
  int nfvtxtrk=0;
//  int nfvtxtrks=0;
//  int nfvtxtrkn=0;

//Correlate histograms
  hbbcsZdcEs->Fill(bbc_qs,ZdcEs);
  hvtxzfvtxz->Fill(vtx_z,fvtxz);
  hpc1hitsbbc->Fill(npc1,bbc_qs+bbc_qn);
  hbbcnbbc->Fill(bbc_qs+bbc_qn,nbbcau+nbbcde);
  hbbcsbbcn->Fill(bbc_qs,bbc_qn);
  hnfvtxtrkbbc->Fill(nfvtxtrk,bbc_qs+bbc_qn);
/*
  if(jevent % 2*nbuff >= nbuff){
      hmixbbcsZdcEs->Fill(bbc_qs+buff_bbcs[ibuff],ZdcEs+buff_ZdcEs[ibuff]);
      buff_bbctau[ibuff].insert(buff_bbctau[ibuff].begin()+buff_nbbcau[ibuff],bbcau_t0,bbcau_t0+nbbcau);
      buff_bbctde[ibuff].insert(buff_bbctde[ibuff].begin()+buff_nbbcde[ibuff],bbcde_t0,bbcde_t0+nbbcde);
      float bbcmixt0sigmasouth = TMath::Sqrt((bbct0sqsouth+buff_bbct0sqsouth[ibuff])/(nbbcau+buff_nbbcau[ibuff])-TMath::Power((bbct0meansouth+buff_bbct0meansouth[ibuff])/(nbbcau+buff_nbbcau[ibuff]),2));
      float bbcmixt0sigmanorth = TMath::Sqrt((bbct0sqnorth+buff_bbct0sqnorth[ibuff])/(nbbcde+buff_nbbcde[ibuff])-TMath::Power((bbct0meannorth+buff_bbct0meannorth[ibuff])/(nbbcde+buff_nbbcde[ibuff]),2));
  //    bool mixPileUpFilterSouthFlag = false;
  //    bool mixPileUpFilterNorthFlag = false;
 //     TH1F hmixbbctsouthtmp = *(TH1F*)bbctsouthtmp.Clone("hmixbbctsouthtmp");
      TH1F hmixbbctsouthtmp(bbctsouthtmp);
      hmixbbctsouthtmp.SetName("hmixbbctsouthtmp");
      if(hmixbbctsouthtmp.GetNbinsX()==buff_bbctsouthtmp[ibuff].GetNbinsX()){
      hmixbbctsouthtmp.Add(&buff_bbctsouthtmp[ibuff]);
      }
 //     TH1F hmixbbctnorthtmp = *(TH1F*)bbctnorthtmp.Clone("hmixbbctnorthtmp");
      TH1F hmixbbctnorthtmp(bbctnorthtmp);
      hmixbbctnorthtmp.SetName("hmixbbctnorthtmp");
      if(hmixbbctnorthtmp.GetNbinsX()==buff_bbctnorthtmp[ibuff].GetNbinsX()){
      hmixbbctnorthtmp.Add(&buff_bbctnorthtmp[ibuff]);
      }

      for(int itry=0;itry<ntry;itry++){
      if(isgood(hmixbbctsouthtmp, bbct_sigma[itry]) && isgood(hmixbbctnorthtmp, bbct_sigma[itry])){
      hcentbbcmixt0sigmasouth[itry]->Fill(bbc_qs+buff_bbcs[ibuff],bbcmixt0sigmasouth);
      hcentbbcmixt0sigmanorth[itry]->Fill(bbc_qs+buff_bbcs[ibuff],bbcmixt0sigmanorth);
      }
      }
      for(int itry=0;itry<nwidth;itry++){
      float mixfracsouth= isgood1(buff_bbctau[ibuff], (bbct0meansouth+buff_bbct0meansouth[ibuff])/(nbbcau+buff_nbbcau[ibuff]),nbbcau+buff_nbbcau[ibuff],bbct_width[itry]);
      float mixfracnorth= isgood1(buff_bbctde[ibuff], (bbct0meannorth+buff_bbct0meannorth[ibuff])/(nbbcde+buff_nbbcde[ibuff]),nbbcde+buff_nbbcde[ibuff],bbct_width[itry]);
      hcentbbcmixt0fracsouth[itry]->Fill(bbc_qs+buff_bbcs[ibuff],mixfracsouth);
      hcentbbcmixt0fracnorth[itry]->Fill(bbc_qs+buff_bbcs[ibuff],mixfracnorth);
      if(bbct_width[itry]==0.5){
          int fracbin = mixfracsouth/0.1;
          if(fracbin == 10) fracbin-=1;
          hcent1cent2[fracbin]->Fill(cent,buff_cent[ibuff]);
      }
      }
  }
  buff_cent[ibuff] = cent;
  buff_bbcs[ibuff] = bbc_qs;
  buff_ZdcEs[ibuff] = ZdcEs;
  buff_bbct0meansouth[ibuff] = bbct0meansouth;
  buff_bbct0meannorth[ibuff] = bbct0meannorth;
  buff_bbct0sqsouth[ibuff] = bbct0sqsouth;
  buff_bbct0sqnorth[ibuff] = bbct0sqnorth;
  buff_nbbcau[ibuff] = nbbcau;
  buff_nbbcde[ibuff] = nbbcde;
  buff_bbctsouthtmp[ibuff] = bbctsouthtmp;
  buff_bbctnorthtmp[ibuff] = bbctnorthtmp;
  vector<float>().swap(buff_bbctde[ibuff]);
  vector<float>().swap(buff_bbctau[ibuff]);
  buff_bbctau[ibuff].insert(buff_bbctau[ibuff].begin() , bbcau_t0 , bbcau_t0 + nbbcau ) ; 
  buff_bbctde[ibuff].insert(buff_bbctde[ibuff].begin() , bbcde_t0 , bbcde_t0 + nbbcde ) ; 
  ibuff ++;
  ibuff=ibuff % nbuff;
  
  vector<float> bbct0au;
  vector<float> bbct0de;
  vector<float>().swap(bbct0au);
  vector<float>().swap(bbct0de);
  bbct0au.insert(bbct0au.begin() , bbcau_t0 , bbcau_t0 + nbbcau ) ; 
  bbct0de.insert(bbct0de.begin() , bbcde_t0 , bbcde_t0 + nbbcde ) ; 
  for(int itry=0;itry<ntry;itry++){
  if (isgood(bbctsouthtmp, bbct_sigma[itry]) && isgood(bbctnorthtmp, bbct_sigma[itry])){
  hcentbbct0sigmasouth[itry]->Fill(bbc_qs, bbct0sigmasouth);
  hcentbbct0sigmanorth[itry]->Fill(bbc_qs, bbct0sigmanorth);
  }
  }
  for(int itry=0;itry<nwidth;itry++){
  float fracsouth= isgood1(bbct0au, bbct0meansouth/nbbcau,nbbcau,bbct_width[itry]);
  float fracnorth= isgood1(bbct0de, bbct0meannorth/nbbcde,nbbcde,bbct_width[itry]);
  hcentbbct0fracsouth[itry]->Fill(bbc_qs,fracsouth);
  hcentbbct0fracnorth[itry]->Fill(bbc_qs,fracnorth);
  }
*/
 // delete bbctsouthtmp;
 // delete bbctnorthtmp;
 // delete hmixbbctsouthtmp;
 // delete hmixbbctnorthtmp;
}
return 0;
}   

//_____________________________________________________________________________________________________________________________
int PerformTestMB::End()
{
  cout << "End of PerformTestMB for Run " << RunNumber << endl;
  cout << "Total # of events = " << ievent << " " << jevent << " " << endl;
  cout << "OutputFileName = " << OutputFileName << endl;

  //HistoManager->dumpHistos(OutputFileName);
  if(d_outfile) {
    d_outfile->cd();
    //  Here I will write the output file...
    
    //hRun->Write();
    //htrig->Write();
//Tracks
  //hcntetaphi->Write();
  //hcntpt->Write();
for(int i=0;i<50;i++){
pc3dphidz_arm0_pos[i]->Write();
pc3dphidz_arm1_pos[i]->Write();
pc3dphidz_arm0_neg[i]->Write();
pc3dphidz_arm1_neg[i]->Write();
}
for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
    for(int i=0;i<50;i++){
        pc3dphidz_arm0_pos_z[ibbcz][i]->Write();
        pc3dphidz_arm1_pos_z[ibbcz][i]->Write();
        pc3dphidz_arm0_neg_z[ibbcz][i]->Write();
        pc3dphidz_arm1_neg_z[ibbcz][i]->Write();
    }
}

/*
//tof
for(int i=0;i<2;i++){
  tofdphidz[i]->Write();
  tofwdphidz[i]->Write();
}
  tofsdphisdz->Write();
  tofwsdphisdz->Write();
  ttofqpratio->Write();
  m2qpratio->Write();
  m2p->Write();
  ttofp->Write();
  pinv2chbeta->Write();
  deltattofeis->Write();

//vtx
for(int i=0;i<4;i++){
  hcluetaphi[i]->Write();
}

//bbc
  bbcet->Write();
  bbctdevsouth->Write();
  bbctdevnorth->Write();
  bbctsouth->Write();
    bbctnorth->Write();

//fvtx
for(int arm=0;arm<2;arm++){
  DCAxydis[arm]->Write();
  DCAxy2dis[arm]->Write();
  DCAcentdis[arm]->Write();
}
//  TH1F *fvtxdphidis->Write();
//  TH1F *fvtxdphidis2->Write();

  hvtx0etaz->Write();
  hvtx1etaz->Write();

  hvtx0etaphi->Write();
  hvtx1etaphi->Write();

//mpc
 
mpcetdis->Write();
for (int i=0; i<centbin; i++) {
mpcetetasouth[i]->Write();
mpcetetanorth[i]->Write();
}
 mpc_south_cent->Write();
 mpc_north_cent->Write();
 south_mpc_north_mpc->Write();
 mpc_south_north->Write();

*/
//correlate
  hvtxzfvtxz->Write();
  hvtxxvtxz->Write();
  hvtxyvtxz->Write();
  hbbcsZdcEs->Write();
  hmixbbcsZdcEs->Write();
  hpc1hitsbbc->Write();
  hnpc3hitsntof->Write();
  hbbcnbbc->Write();
  hbbcsbbcn->Write();
  /*
for(int i=0;i<4;i++){
  hnvtxnfvtxtrk[i]->Write();
  hnvtxnmpc[i]->Write();
  hbbcsnvtx[i]->Write();
  hbbcnnvtx[i]->Write();
  hbbcnvtx[i]->Write();
}
  hnbbcnclu->Write();
  hntracknmpc->Write();
  hnfvtxtrkbbc->Write();
  hnfvtxtrksnmpcs->Write();
  hnfvtxtrknnmpcn->Write();
  south_mpc_south_bbc->Write();
  north_mpc_north_bbc->Write();
  */
  /*
for(int itry=0;itry<ntry;itry++){
  hcentbbct0sigmasouth[itry]->Write();
  hcentbbct0sigmanorth[itry]->Write();
  hcentbbcmixt0sigmasouth[itry]->Write();
  hcentbbcmixt0sigmanorth[itry]->Write();
}
for(int itry=0;itry<nwidth;itry++){
  hcent1cent2[itry]->Write();
  hcentbbct0fracsouth[itry]->Write();
  hcentbbct0fracnorth[itry]->Write();
  hcentbbcmixt0fracsouth[itry]->Write();
  hcentbbcmixt0fracnorth[itry]->Write();
}
 */ 
    //rxngain ->Write();

    //htree->Write();
    
    d_outfile->Close();
    //delete d_outfile;
    
  }else{

    cout<<"ERROR: No output file set!"<<endl;

  }
  return 0;
}

bool isgood(TH1F h, float sigma){
  if(sigma == -1) return true;
  bool Flag = true;
  int imaxbin = h.GetMaximumBin();
  float bincenter = h.GetBinCenter(imaxbin);
  for(int ibin=0;ibin<h.GetNbinsX();ibin++){
      if((h.GetBinCenter(ibin)<(bincenter-sigma)) || (h.GetBinCenter(ibin)>(bincenter+sigma))){
          if(h.GetBinContent(ibin)!=0){
              Flag = false;
              break;
          }
      }
  }

  if(Flag){
  h.SetBinContent(imaxbin,0);
  imaxbin = h.GetMaximumBin();
  bincenter = h.GetBinCenter(imaxbin);
  for(int ibin=0;ibin<h.GetNbinsX();ibin++){
      if((h.GetBinCenter(ibin)<(bincenter-sigma)) || (h.GetBinCenter(ibin)>(bincenter+sigma))){
          if(h.GetBinContent(ibin)!=0){
              Flag = false;
              break;
          }
      }
  }
  }
  return Flag;
}

float isgood1(vector<float> bbct0, float bbct0mean, int nbbc_, float sigma){
  if(sigma == -1) return true;
  int countin=0;
  for(int i=0;i<nbbc_;i++){
      if(fabs(bbct0[i]-bbct0mean)<=sigma) countin++;
  }
  return 1.*countin/nbbc_;
}
