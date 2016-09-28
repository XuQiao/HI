#include "RidgedAuRun16.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include "TString.h"
#include "TRandom.h"
#include <fstream>
#include <iostream>

#include "Run16dAupc3dphidzcalibsmoothpass1.h"
#include "func.h"

using namespace std;

int get_fvtx_layer(float);
void initialize_pmt_position();
float d_pmt_z = -1443.5; // same for all tubes

const int bbcz_cut = 10;
const int vtxz_cut = 10;

const int ntrks = 200;
const int nbbc = 800;

const int centbin = 6;
const int ptbin = 25;
const int vzbin = 10;
const int nbuff = 40;

const int nbufftrks = centbin*vzbin*nbuff*ntrks;

const int nbuffbbcs = centbin*vzbin*nbuff*nbbc;

//typedef vector<float> track_buff;
float trackbuff_pt[nbufftrks];
float trackbuff_phi[nbufftrks];
float trackbuff_eta[nbufftrks];

//typedef vector<float> bbcs_buff;
float bbcsbuff_et[nbuffbbcs];
float bbcsbuff_phi[nbuffbbcs];
float bbcsbuff_eta[nbuffbbcs];

//typedef vector<float> bbcn_buff;
float bbcnbuff_et[nbuffbbcs];
float bbcnbuff_phi[nbuffbbcs];
float bbcnbuff_eta[nbuffbbcs];

int ntrack_buff[centbin][vzbin];
int dtrack_buff[centbin][vzbin];
int buff_ntrack[centbin][vzbin][nbuff];

int nbbcs_buff[centbin][vzbin];
int dbbcs_buff[centbin][vzbin];
int buff_nbbcs[centbin][vzbin][nbuff];

int nbbcn_buff[centbin][vzbin];
int dbbcn_buff[centbin][vzbin];
int buff_nbbcn[centbin][vzbin][nbuff];


//_____________________________________________________________________________________________________________________________
RidgedAuRun16::RidgedAuRun16(std::vector<TString> input, std::vector<TString> input1, const char* output) :
  OutputFileName(output), InputFileName(input), InputFileName1(input1), ievent(0), jevent(0), nevent(0), RunNumber(0)
{

  d_outfile=NULL;
 // htree=NULL;
  
  hcentbbcs = NULL;
  hcentnfvtxs = NULL;
for(int i=0;i<centbin;i++){
  hrunbbcs[i] = NULL;
  hrunnfvtxs[i] = NULL;
  hrunntrack[i] = NULL;
  hbbcsnfvtxs[i] = NULL;

     for(int j=0;j<40;j++){

       //BBC
       //normal
        hforesouthbbc[i][j]=NULL;
        hbacksouthbbc[i][j]=NULL;
        //flip
        kforesouthbbc[i][j]=NULL;
        kbacksouthbbc[i][j]=NULL;
        //with weight
        hforesouthbbcw[i][j]=NULL;
        hbacksouthbbcw[i][j]=NULL;
        //flip
        kforesouthbbcw[i][j]=NULL;
        kbacksouthbbcw[i][j]=NULL;
        //tower mixing
        hbacksouthbbc2[i][j]=NULL;
        kbacksouthbbc2[i][j]=NULL;
        
        hbacksouthbbcw2[i][j]=NULL;
        kbacksouthbbcw2[i][j]=NULL;
  //north
        hforenorthbbc[i][j]=NULL;
        hbacknorthbbc[i][j]=NULL;

        kforenorthbbc[i][j]=NULL;
        kbacknorthbbc[i][j]=NULL;
  
        hforenorthbbcw[i][j]=NULL;
        hbacknorthbbcw[i][j]=NULL;

        kforenorthbbcw[i][j]=NULL;
        kbacknorthbbcw[i][j]=NULL;
        
        hbacknorthbbc2[i][j]=NULL;
        kbacknorthbbc2[i][j]=NULL;
        
        hbacknorthbbcw2[i][j]=NULL;
        kbacknorthbbcw2[i][j]=NULL;
  //south-north
        hforesnbbc[i][j]=NULL;
        hbacksnbbc[i][j]=NULL;

        kforesnbbc[i][j]=NULL;
        kbacksnbbc[i][j]=NULL;
  
        hforesnbbcw[i][j]=NULL;
        hbacksnbbcw[i][j]=NULL;

        kforesnbbcw[i][j]=NULL;
        kbacksnbbcw[i][j]=NULL;
        
        hbacksnbbc2[i][j]=NULL;
        kbacksnbbc2[i][j]=NULL;
        
        hbacksnbbcw2[i][j]=NULL;
        kbacksnbbcw2[i][j]=NULL;
//2D
        //kforesouthetabbc[i][j]=NULL;
        //kbacksouthetabbc[i][j]=NULL;
        
//        kforenorthetabbc[i][j]=NULL;
//        kbacknorthetabbc[i][j]=NULL;

//        kforesouthetabbcw[i][j]=NULL;
//        kbacksouthetabbcw[i][j]=NULL;

//        kforenorthetabbcw[i][j]=NULL;
//        kbacknorthetabbcw[i][j]=NULL;
  
//        kbacksouthetabbc2[i][j]=NULL;
//        kbacknorthetabbc2[i][j]=NULL;

//        kbacksouthetabbcw2[i][j]=NULL;
//        kbacknorthetabbcw2[i][j]=NULL;
     }
  }
    
    for(int iarm=0;iarm<2;iarm++){
        DCAxydis[iarm]=NULL;
        DCAxy2dis[iarm]=NULL;
        DCAcentdis[iarm]=NULL;
    }
      fvtxdphidis=NULL;
      fvtxdphidis2=NULL;
    hvtx0etaz=NULL;
    hvtx1etaz=NULL;

    hvtx0etaphi=NULL;
    hvtx1etaphi=NULL;


  pi=0.0;

}

//_____________________________________________________________________________________________________________________________
RidgedAuRun16::~RidgedAuRun16()
{
  cout << " RidgedAuRun16::~RidgedAuRun16 " << endl;
	}

//_____________________________________________________________________________________________________________________________
int RidgedAuRun16::Init()
{
  cout << " RidgedAuRun16::Init " << endl;
  char name[80];
  pi=acos(-1.0);

  d_outfile = new TFile(OutputFileName.c_str(),"recreate");

  sprintf(name,"hcentbbcs"); hcentbbcs = new TH2F(name,name,100,0,100,300,0,300);
  sprintf(name,"hcentnfvtxs"); hcentnfvtxs = new TH2F(name,name,100,0,100,100,0,100);
  for (int icent=0; icent<centbin; icent++) {
  sprintf(name,"hrunbbcs_%d",icent); hrunbbcs[icent] = new TProfile(name,name,460000-454000,454000,460000,0,300);
  sprintf(name,"hrunnfvtxs_%d",icent); hrunnfvtxs[icent] = new TProfile(name,name,460000-454000,454000,460000,0,100);
  sprintf(name,"hrunntrack_%d",icent); hrunntrack[icent] = new TProfile(name,name,460000-454000,454000,460000,0,50);
  //sprintf(name,"hbbcsnfvtxs_%d",icent); hbbcsnfvtxs[icent] = new TH2F(name,name,600,0,300,50,0,50);
  sprintf(name,"hbbcsnfvtxs_%d",icent); hbbcsnfvtxs[icent] = new TH2F(name,name,50,0,50,50,0,50);

    //bbc correlation
    for(int ipt=0; ipt<ptbin; ipt++){
      //cnt with bbc
      sprintf(name,"hforesouthbbc_%d_%d",icent,ipt);
      hforesouthbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacksouthbbc_%d_%d",icent,ipt);
      hbacksouthbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      //flip
      sprintf(name,"kforesouthbbc_%d_%d",icent,ipt);
      kforesouthbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksouthbbc_%d_%d",icent,ipt);
      kbacksouthbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      //with weight
      sprintf(name,"hforesouthbbcw_%d_%d",icent,ipt);
      hforesouthbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacksouthbbcw_%d_%d",icent,ipt);
      hbacksouthbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      //flip
      sprintf(name,"kforesouthbbcw_%d_%d",icent,ipt);
      kforesouthbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksouthbbcw_%d_%d",icent,ipt);
      kbacksouthbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      //mixing with bbc
      sprintf(name,"hbacksouthbbc2_%d_%d",icent,ipt);
      hbacksouthbbc2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);
      
      //flip
      sprintf(name,"kbacksouthbbc2_%d_%d",icent,ipt);
      kbacksouthbbc2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      //with weight
      sprintf(name,"hbacksouthbbcw2_%d_%d",icent,ipt);
      hbacksouthbbcw2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      //flip
      sprintf(name,"kbacksouthbbcw2_%d_%d",icent,ipt);
      kbacksouthbbcw2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      //north
      sprintf(name,"hforenorthbbc_%d_%d",icent,ipt);
      hforenorthbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacknorthbbc_%d_%d",icent,ipt);
      hbacknorthbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kforenorthbbc_%d_%d",icent,ipt);
      kforenorthbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthbbc_%d_%d",icent,ipt);
      kbacknorthbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hforenorthbbcw_%d_%d",icent,ipt);
      hforenorthbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacknorthbbcw_%d_%d",icent,ipt);
      hbacknorthbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kforenorthbbcw_%d_%d",icent,ipt);
      kforenorthbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthbbcw_%d_%d",icent,ipt);
      kbacknorthbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);
      
      sprintf(name,"hbacknorthbbc2_%d_%d",icent,ipt);
      hbacknorthbbc2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthbbc2_%d_%d",icent,ipt);
      kbacknorthbbc2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacknorthbbcw2_%d_%d",icent,ipt);
      hbacknorthbbcw2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthbbcw2_%d_%d",icent,ipt);
      kbacknorthbbcw2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      //south - north
      sprintf(name,"hforesnbbc_%d_%d",icent,ipt);
      hforesnbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacksnbbc_%d_%d",icent,ipt);
      hbacksnbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kforesnbbc_%d_%d",icent,ipt);
      kforesnbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksnbbc_%d_%d",icent,ipt);
      kbacksnbbc[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hforesnbbcw_%d_%d",icent,ipt);
      hforesnbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacksnbbcw_%d_%d",icent,ipt);
      hbacksnbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kforesnbbcw_%d_%d",icent,ipt);
      kforesnbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksnbbcw_%d_%d",icent,ipt);
      kbacksnbbcw[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);
      
      sprintf(name,"hbacksnbbc2_%d_%d",icent,ipt);
      hbacksnbbc2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksnbbc2_%d_%d",icent,ipt);
      kbacksnbbc2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacksnbbcw2_%d_%d",icent,ipt);
      hbacksnbbcw2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksnbbcw2_%d_%d",icent,ipt);
      kbacksnbbcw2[icent][ipt] = new TH1D(name,name,40, -0.5*pi, 1.5*pi);

/*
      //2D
      sprintf(name,"kforesouthetabbc_%d_%d",icent,ipt);
      kforesouthetabbc[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksouthetabbc_%d_%d",icent,ipt);
      kbacksouthetabbc[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kforenorthetabbc_%d_%d",icent,ipt);
      kforenorthetabbc[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthetabbc_%d_%d",icent,ipt);
      kbacknorthetabbc[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksouthetabbc2_%d_%d",icent,ipt);
      kbacksouthetabbc2[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthetabbc2_%d_%d",icent,ipt);
      kbacknorthetabbc2[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      //with weight
      sprintf(name,"kforesouthetabbcw_%d_%d",icent,ipt);
      kforesouthetabbcw[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksouthetabbcw_%d_%d",icent,ipt);
      kbacksouthetabbcw[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kforenorthetabbcw_%d_%d",icent,ipt);
      kforenorthetabbcw[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthetabbcw_%d_%d",icent,ipt);
      kbacknorthetabbcw[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksouthetabbcw2_%d_%d",icent,ipt);
      kbacksouthetabbcw2[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthetabbcw2_%d_%d",icent,ipt);
      kbacknorthetabbcw2[icent][ipt] = new TH2F(name, name, 100, -5.0, 5.0, 40, -0.5*pi, 1.5*pi);
      */
    }
  }

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

  
  memset(ntrack_buff, 0, sizeof(ntrack_buff));
  memset(dtrack_buff, 0, sizeof(dtrack_buff));
  memset(buff_ntrack, 0, sizeof(buff_ntrack));

  memset(nbbcs_buff, 0, sizeof(nbbcs_buff));
  memset(dbbcs_buff, 0, sizeof(dbbcs_buff));
  memset(buff_nbbcs, 0, sizeof(buff_nbbcs));

  memset(nbbcn_buff, 0, sizeof(nbbcn_buff));
  memset(dbbcn_buff, 0, sizeof(dbbcn_buff));
  memset(buff_nbbcn, 0, sizeof(buff_nbbcn));

  memset(trackbuff_pt, 0, sizeof(trackbuff_pt));
  memset(trackbuff_phi, 0, sizeof(trackbuff_phi));
  memset(trackbuff_eta, 0, sizeof(trackbuff_eta));
  
  memset(bbcnbuff_et, 0, sizeof(bbcnbuff_et));
  memset(bbcnbuff_phi, 0, sizeof(bbcnbuff_phi));
  memset(bbcnbuff_eta, 0, sizeof(bbcnbuff_eta));

  memset(bbcsbuff_et, 0, sizeof(bbcsbuff_et));
  memset(bbcsbuff_phi, 0, sizeof(bbcsbuff_phi));
  memset(bbcsbuff_eta, 0, sizeof(bbcsbuff_eta));

  cout<<"finish of initialize"<<endl;
  return 0;
}

int RidgedAuRun16::Inittree(){
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
//  TBranch* b_event;   //!
  TBranch* b_bbc_z;   //!
  TBranch* b_centrality;   //!
  TBranch* b_bbc_qn;   //!
  TBranch* b_bbc_qs;   //!
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

  
  tree->SetBranchAddress("run",&RunNumber,&b_bbc_z);
  tree->SetBranchAddress("bbcv",&d_bbcz,&b_bbc_z);
  tree->SetBranchAddress("cent",&centrality,&b_centrality);
  tree->SetBranchAddress("bbc_n",&bbc_qn,&b_bbc_qn);
  tree->SetBranchAddress("bbc_s",&bbc_qs,&b_bbc_qs);
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
int RidgedAuRun16::process_event(){
  int nEvent = tree->GetEntries();
  cout<<nEvent<<endl;
  for(ievent=0;ievent < nEvent; ievent++){
      tree->GetEntry(ievent);
      tree1->GetEntry(ievent);
  if(ievent%10000==0) {
    cout<<"************* ievent= "<<ievent<<"    *************"<<endl;
  }

//event check
   if(d_nfvtxtrk != d_nfvtxtrk1) {cout << d_nfvtxtrk << " EVENT NOT SYNC ! " << d_nfvtxtrk1 << endl; continue;}

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
  float bbc_s = bbc_qs;
  float vvertex = bbcv;
  
  float fvtxz = eventfvtx_z;
  if(fvtxz != fvtxz){
      fvtxz = -9999;
  }
  
  if(RunNumber>=455792 && RunNumber<=458167){ // low energies use fvtx as vertex if there is one
  if(fvtxz != -9999)
    vvertex = fvtxz;
  }

 int ivz = -1;
 ivz = vzbin*(vvertex+10)/20;
 if(ivz<0||ivz>=vzbin) continue;

 int icent = -9999;
 if(cent<0) continue;
 if(cent<=5) icent = 0;
 else if(cent<=10) icent = 1;
 else if(cent<=20) icent = 2;
 else if(cent<=40) icent = 3;
 else if(cent<=60) icent = 4;
 else icent = 5;
 if(icent>=centbin||icent<0) continue;

 // --- all numbers from Darren 2016-06-23
      const float x_off = 0.3;
      const float beam_angle = 0.001;
      float vtx_z = vvertex;
      float vtx_x = x_off + atan(beam_angle)*vtx_z;
      float vtx_y = 0.02;
 if ( RunNumber >= 456652 && RunNumber <= 457298 && d_nFVTX_clus > 300) continue;
 if ( RunNumber >= 457634 && RunNumber <= 458167 && d_nFVTX_clus > 500) continue;

  jevent++;

  int ntrack = 0;
  float track_pt[ntrks];
  float track_phi[ntrks];
  float track_eta[ntrks];

   for(int itrk=0; itrk< d_ntrk; itrk++){
      float mom    = d_mom[itrk];
      float phi0    = d_phi0[itrk];
      float the0    = d_the0[itrk];
      float eta       = -log(tan(0.5*the0));
      int charge = d_charge[itrk];
      float pc3dphi   = d_pc3dphi[itrk];
      float pc3dz   = d_pc3dz[itrk];
      float pt        = mom*sin(the0);
      float pz        = mom*cos(the0);
      float px = pt * cos(phi0);
      float py = pt * sin(phi0);
      px = pz*sin(-beam_angle) + px*cos(-beam_angle);
      float phi = atan2(py,px);
      int dcarm=0;
      if(px>0) dcarm=1;

      double sdphi = calcsdphi(pc3dphi,dcarm,charge,mom,RunNumber);
      double sdz =  calcsdz(pc3dz,dcarm,charge,mom,RunNumber);
      if(fabs(sdphi)<2.0 && fabs(sdz)<2.0){

     // if(m2tof>0.6 && m2tof<1.2)   //proton
      //if(m2tof>-0.2 && m2tof<0.1) //pion
      //if(m2tof>0.15 && m2tof<0.45) //kaon

      if(pt>0.2&&pt<5.0){
	track_pt[ntrack]=pt;
	track_phi[ntrack]=phi;
	track_eta[ntrack]=eta;
	ntrack++;
	}
      }
    } //track loop
 

  // beam beam counter (bbc) r.p. // -----------------------------------
  int nbbcs=0;//bbc south
  float bbcs_et[nbbc];
  float bbcs_phi[nbbc];
  float bbcs_eta[nbbc];

 // int nbbcn=0;//bbc north
 // float bbcn_et[nbbc];
 // float bbcn_phi[nbbc];
 // float bbcn_eta[nbbc];
    int jbbc = 0;
 for(int ipmt = 0; ipmt < 128; ipmt++){
        float bbcx, bbcy, bbcz;
        if(!d_BBC_valid[ipmt]) continue;
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
    float charge = d_BBC_charge[jbbc];
    float time0 = d_BBC_time0[jbbc];
    jbbc++;
    bbcx = bbcx - vtx_x*10.0;
    bbcy = bbcy - vtx_y*10.0;
    bbcz = bbcz - vtx_z*10.0;
    // --- rotation
    bbcx = bbcz*sin(-beam_angle) + bbcx*cos(-beam_angle);
    int iarm = 0;
    if(bbcz>0) iarm = 1;
    if (time0>0 && charge>0) {
      float phi=atan2(bbcy,bbcx);
      float val = charge;
      float rad = sqrt(bbcx*bbcx+bbcy*bbcy);
      float the = atan2(rad,bbcz - vtx_z*10.0); // bbcz-bbcv*10.0);
      float eta = -log(tan(0.5*the));

	if(iarm==0){
	  bbcs_et[nbbcs] = val;
	  bbcs_phi[nbbcs] = phi; 
	  bbcs_eta[nbbcs] = eta; 
	  nbbcs++;
	}
	else{
	 // bbcn_et[nbbcn] = val;
         // bbcn_phi[nbbcn] = phi;
         // bbcn_eta[nbbcn] = eta;
         // nbbcn++;
	}
      }
    }

  
//FVTX

  int nfvtxs=0;//fvtxau in the et cut
  float fvtxs_phi[nbbc];
  float fvtxs_eta[nbbc];
  float fvtxs_et[nbbc];

//  int nfvtxn=0;//fvtxau in the et cut
//  float fvtxn_phi[nbbc];
//  float fvtxn_eta[nbbc];
//  float fvtxn_et[nbbc];
  
    for(int iftrk = 0; iftrk < d_nfvtxtrk; iftrk++){
      int fvtx_arm = d_farm[iftrk];
      int nhits = d_fnhits[iftrk];
      float fvtx_eta = d_feta[iftrk];
      float fvtx_the = 2.*atan(exp(-fvtx_eta));
      float fvtx_phi = d_fphi[iftrk];

      float     fvtx_x = d_fvtxX[iftrk];
      float     fvtx_y = d_fvtxY[iftrk];
      float     fvtx_z = d_fvtxZ[iftrk];
      float     chi2   = d_fvtxchi[iftrk];
      if(fvtx_arm) continue;
      if(fvtx_the==0) continue;
      if(chi2>4) continue;
 
      float DCA_x = fvtx_x + tan(fvtx_the)*cos(fvtx_phi)*(fvtxz - fvtx_z);
      float DCA_y = fvtx_y + tan(fvtx_the)*sin(fvtx_phi)*(fvtxz - fvtx_z);

     bool dcacut = fabs(DCA_x)<2.0 && fabs(DCA_y)<2.0; //fabs(sigma_dcax)<2.0 && fabs(sigma_dcay)<2.0;
      if(fvtx_phi<10 && fvtx_phi > -10 && fabs(fvtx_eta)<3.5 && dcacut && nhits>=3){
          if(fvtx_eta<-1.0 && fvtx_eta>-3.0){
            fvtxs_phi[nfvtxs] = fvtx_phi;
            fvtxs_eta[nfvtxs] = fvtx_eta;
            fvtxs_et[nfvtxs] = 1.;
            nfvtxs++;
          }
 //         if(fvtx_eta>1.0 && fvtx_eta<3.0){
 //           bbcs_phi[nbbcs] = fvtx_phi;
 //           bbcs_eta[nbbcs] = fvtx_eta;
 //           bbcs_et[nbbcs] = 1.;
 //           nbbcs++;
 //         }
        }
  }

  hrunntrack[icent]->Fill(RunNumber,ntrack);
  hrunnfvtxs[icent]->Fill(RunNumber,nfvtxs);
  hrunbbcs[icent]->Fill(RunNumber,bbc_s);
  hbbcsnfvtxs[icent]->Fill(bbc_s,nfvtxs);
  hcentnfvtxs->Fill(cent,nfvtxs);
  hcentbbcs->Fill(cent,bbc_s);

  //correlation
  int sign =(int)(gRandom->Rndm()*2);//put the flip sign here since per track/evt
  //1) //foreground
  
//  float sum_mpcsouth[ptbin]={0,0,0,0,0,0,0,0,0,0};
//  float sum_mpcnorth[ptbin]={0,0,0,0,0,0,0,0,0,0};

// float sum_bbcsouth[ptbin]={0,0,0,0,0,0,0,0,0,0};
//  float sum_bbcnorth[ptbin]={0,0,0,0,0,0,0,0,0,0};

//  float ptmean[ptbin]={0.35, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75};
  

  for(int itrk=0; itrk<ntrack; itrk++){
    int ipt_cor = track_pt[itrk]/0.2;
    //BBC
    for(int ipmt=0; ipmt<nbbcs; ipmt++){
      float val = bbcs_et[ipmt];
      float bbc_phi = bbcs_phi[ipmt];
      //float bbc_eta = bbcs_eta[ipmt];

      float dphi = track_phi[itrk] - bbc_phi;
      //float deta = track_eta[itrk] - bbc_eta;
      //int sign =(int)(gRandom->Rndm()*2);

      if(dphi<-0.5*pi) dphi+=2*pi;
      if(dphi>1.5*pi) dphi-=2*pi;

      float dphia = dphi;
      //if(sign) {
      //if(dphi<0.5*pi) dphia  = -dphi;
      //else dphia = 2*pi - dphi;
      //}

      if(sign) dphia = -dphia;
      if(dphia<-0.5*pi) dphia+=2*pi;
      if(dphia>1.5*pi) dphia-=2*pi;

      //sum_bbcsouth[ipt_cor]+=val*cos(2*dphi);

      hforesouthbbc[icent][ipt_cor]->Fill(dphi);
      kforesouthbbc[icent][ipt_cor]->Fill(dphia);
//      kforesouthetabbc[icent][ipt_cor]->Fill(deta, dphia);

      hforesouthbbcw[icent][ipt_cor]->Fill(dphi,val*val);
      kforesouthbbcw[icent][ipt_cor]->Fill(dphia,val);
//    kforesouthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);
    }//bbc south

    for(int ipmt=0; ipmt<nfvtxs; ipmt++){
      float val = fvtxs_et[ipmt];
      float bbc_phi = fvtxs_phi[ipmt];
      //float bbc_eta = fvtxs_eta[ipmt];

      float dphi = track_phi[itrk] - bbc_phi;
      //float deta = track_eta[itrk] - bbc_eta;
      //int sign =(int)(gRandom->Rndm()*2);

      if(dphi<-0.5*pi) dphi+=2*pi;
      if(dphi>1.5*pi) dphi-=2*pi;

      float dphia = dphi;
      //if(sign) {
      //if(dphi<0.5*pi) dphia  = -dphi;
      //else dphia = 2*pi - dphi;
      //}

      if(sign) dphia = -dphia;
      if(dphia<-0.5*pi) dphia+=2*pi;
      if(dphia>1.5*pi) dphia-=2*pi;

      //sum_bbcnorth[ipt_cor]+=val*cos(2*dphi);

      hforenorthbbc[icent][ipt_cor]->Fill(dphi);
      kforenorthbbc[icent][ipt_cor]->Fill(dphia);
//      kforenorthetabbc[icent][ipt_cor]->Fill(deta, dphia);

      hforenorthbbcw[icent][ipt_cor]->Fill(dphi,val*val);
      kforenorthbbcw[icent][ipt_cor]->Fill(dphia,val);
//      kforenorthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);
    }//bbc north
  }
  
  for(int iftrk=0; iftrk<nfvtxs; iftrk++){
    int ifpt_cor = 0;
    //BBC
    for(int ipmt=0; ipmt<nbbcs; ipmt++){
      float val = bbcs_et[ipmt]*fvtxs_et[iftrk];
      float bbc_phi = bbcs_phi[ipmt];
      //float bbc_eta = bbcs_eta[ipmt];

      float dphi = fvtxs_phi[iftrk] - bbc_phi;
      //float deta = fvtxs_eta[iftrk] - bbc_eta;
      //int sign =(int)(gRandom->Rndm()*2);

      if(dphi<-0.5*pi) dphi+=2*pi;
      if(dphi>1.5*pi) dphi-=2*pi;

      float dphia = dphi;
      //if(sign) {
      //if(dphi<0.5*pi) dphia  = -dphi;
      //else dphia = 2*pi - dphi;
      //}

      if(sign) dphia = -dphia;
      if(dphia<-0.5*pi) dphia+=2*pi;
      if(dphia>1.5*pi) dphia-=2*pi;

      //sum_bbcsouth[ipt_cor]+=val*cos(2*dphi);

      hforesnbbc[icent][ifpt_cor]->Fill(dphi);
      kforesnbbc[icent][ifpt_cor]->Fill(dphia);
//      kforesouthetabbc[icent][ipt_cor]->Fill(deta, dphia);

      hforesnbbcw[icent][ifpt_cor]->Fill(dphi,val*val);
      kforesnbbcw[icent][ifpt_cor]->Fill(dphia,val);
//    kforesouthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);
    }//south-north
  }
  
  //2) associate tower + mixing track
  if(dtrack_buff[icent][ivz]>0){
    
    for(int ibuff=0;ibuff<dtrack_buff[icent][ivz];ibuff++){
	
      for(int itrk_buff=0; itrk_buff<buff_ntrack[icent][ivz][ibuff]; itrk_buff++){
	
	int ipt_cor = trackbuff_pt[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff]/0.2;
	
	for(int ipmt=0; ipmt<nbbcs; ipmt++){
          float val = bbcs_et[ipmt];
          float bbc_phi = bbcs_phi[ipmt];
          //float bbc_eta = bbcs_eta[ipmt];

          float dphi = trackbuff_phi[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff] - bbc_phi;
	  //float deta = trackbuff_eta[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff] - bbc_eta;
          //int sign =(int)(gRandom->Rndm()*2);

          if(dphi<-0.5*pi) dphi+=2*pi;
          if(dphi>1.5*pi) dphi-=2*pi;

          float dphia = dphi;
          //if(sign) {
	  //if(dphi<0.5*pi) dphia  = -dphi;
	  //else dphia = 2*pi - dphi;
          //}

	  if(sign) dphia = -dphia;
	  if(dphia<-0.5*pi) dphia+=2*pi;
	  if(dphia>1.5*pi) dphia-=2*pi;

          hbacksouthbbc[icent][ipt_cor]->Fill(dphi);
          kbacksouthbbc[icent][ipt_cor]->Fill(dphia);
//          kbacksouthetabbc[icent][ipt_cor]->Fill(deta, dphia);

          hbacksouthbbcw[icent][ipt_cor]->Fill(dphi,val*val);
          kbacksouthbbcw[icent][ipt_cor]->Fill(dphia,val);
//          kbacksouthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);

        }//bbcsouth

	for(int ipmt=0; ipmt<nfvtxs; ipmt++){
          float val = fvtxs_et[ipmt];
          float bbc_phi = fvtxs_phi[ipmt];
          //float bbc_eta = fvtxs_eta[ipmt];

          float dphi = trackbuff_phi[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff] - bbc_phi;
          //float deta = trackbuff_eta[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff] - bbc_eta;
          //int sign =(int)(gRandom->Rndm()*2);

          if(dphi<-0.5*pi) dphi+=2*pi;
          if(dphi>1.5*pi) dphi-=2*pi;

          float dphia = dphi;
          //if(sign) {
	  //if(dphi<0.5*pi) dphia  = -dphi;
	  //else dphia = 2*pi - dphi;
          //}

	  if(sign) dphia = -dphia;
	  if(dphia<-0.5*pi) dphia+=2*pi;
	  if(dphia>1.5*pi) dphia-=2*pi;

          hbacknorthbbc[icent][ipt_cor]->Fill(dphi);
          kbacknorthbbc[icent][ipt_cor]->Fill(dphia);
 //         kbacknorthetabbc[icent][ipt_cor]->Fill(deta, dphia);
          
          hbacknorthbbcw[icent][ipt_cor]->Fill(dphi,val*val);
          kbacknorthbbcw[icent][ipt_cor]->Fill(dphia,val);
 //         kbacknorthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);
        }//bbcnorth

      }
    }
  }
  
  if(dbbcn_buff[icent][ivz]>0){
    
    for(int ibuff=0;ibuff<dbbcn_buff[icent][ivz];ibuff++){
	
      for(int iftrk_buff=0; iftrk_buff<buff_nbbcn[icent][ivz][ibuff]; iftrk_buff++){
	
	int ifpt_cor = 0;
	
	for(int ipmt=0; ipmt<nbbcs; ipmt++){
          float val = bbcs_et[ipmt]*bbcnbuff_et[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+iftrk_buff];
          float bbc_phi = bbcs_phi[ipmt];
          //float bbc_eta = bbcs_eta[ipmt];

          float dphi = bbcnbuff_phi[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+iftrk_buff]- bbc_phi;
	  //float deta = bbcnbuff_eta[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+iftrk_buff] - bbc_eta;
          //int sign =(int)(gRandom->Rndm()*2);

          if(dphi<-0.5*pi) dphi+=2*pi;
          if(dphi>1.5*pi) dphi-=2*pi;

          float dphia = dphi;
          //if(sign) {
	  //if(dphi<0.5*pi) dphia  = -dphi;
	  //else dphia = 2*pi - dphi;
          //}

	  if(sign) dphia = -dphia;
	  if(dphia<-0.5*pi) dphia+=2*pi;
	  if(dphia>1.5*pi) dphia-=2*pi;

          hbacksnbbc[icent][ifpt_cor]->Fill(dphi);
          kbacksnbbc[icent][ifpt_cor]->Fill(dphia);
 //         kbacksnetabbc[icent][ipt_cor]->Fill(deta, dphia);

          hbacksnbbcw[icent][ifpt_cor]->Fill(dphi,val*val);
          kbacksnbbcw[icent][ifpt_cor]->Fill(dphia,val);
 //         kbacksouthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);

        }//bbcsouth-north
      }
    }
  }

  //3) trigger track + mixing tower
  for(int itrk=0; itrk<ntrack; itrk++){
    int ipt_cor = track_pt[itrk]/0.2;

    if(dbbcs_buff[icent][ivz]>0){
      for(int ibuff=0;ibuff<dbbcs_buff[icent][ivz];ibuff++){
        for(int ipmt_buff=0; ipmt_buff<buff_nbbcs[icent][ivz][ibuff]; ipmt_buff++){
          float val = bbcsbuff_et[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ipmt_buff];
          float dphi = track_phi[itrk] - bbcsbuff_phi[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ipmt_buff];
          //float deta = track_eta[itrk] - bbcsbuff_eta[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ipmt_buff];
          //int sign =(int)(gRandom->Rndm()*2);
          if(dphi<-0.5*pi) dphi+=2*pi;
          if(dphi>1.5*pi) dphi-=2*pi;

          float dphia = dphi;
          //if(sign) {
	  //if(dphi<0.5*pi) dphia  = -dphi;
	  //else dphia = 2*pi - dphi;
          //}

	  if(sign) dphia = -dphia;
	  if(dphia<-0.5*pi) dphia+=2*pi;
	  if(dphia>1.5*pi) dphia-=2*pi;

          hbacksouthbbc2[icent][ipt_cor]->Fill(dphi);
          kbacksouthbbc2[icent][ipt_cor]->Fill(dphia);
//          kbacksouthetabbc2[icent][ipt_cor]->Fill(deta, dphia);

          hbacksouthbbcw2[icent][ipt_cor]->Fill(dphi,val*val);
          kbacksouthbbcw2[icent][ipt_cor]->Fill(dphia,val);
//          kbacksouthetabbcw2[icent][ipt_cor]->Fill(deta, dphia, val);
        }
      }
    }//end bbc south

    if(dbbcn_buff[icent][ivz]>0){
      for(int ibuff=0;ibuff<dbbcn_buff[icent][ivz];ibuff++){
        for(int ipmt_buff=0; ipmt_buff<buff_nbbcn[icent][ivz][ibuff]; ipmt_buff++){
          float val = bbcnbuff_et[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ipmt_buff];
          float dphi = track_phi[itrk] - bbcnbuff_phi[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ipmt_buff];
          //float deta = track_eta[itrk] - bbcnbuff_eta[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ipmt_buff];
          //int sign =(int)(gRandom->Rndm()*2);

          if(dphi<-0.5*pi) dphi+=2*pi;
          if(dphi>1.5*pi) dphi-=2*pi;

          float dphia = dphi;
          //if(sign) {
	  //if(dphi<0.5*pi) dphia  = -dphi;
	  //else dphia = 2*pi - dphi;
          //}

	  if(sign) dphia = -dphia;
	  if(dphia<-0.5*pi) dphia+=2*pi;
	  if(dphia>1.5*pi) dphia-=2*pi;

          hbacknorthbbc2[icent][ipt_cor]->Fill(dphi);
          kbacknorthbbc2[icent][ipt_cor]->Fill(dphia);
//          kbacknorthetabbc2[icent][ipt_cor]->Fill(deta, dphia);
          
	  hbacknorthbbcw2[icent][ipt_cor]->Fill(dphi,val*val);
          kbacknorthbbcw2[icent][ipt_cor]->Fill(dphia,val);
//          kbacknorthetabbcw2[icent][ipt_cor]->Fill(deta, dphia, val);
        }
      }
    }//end bbc north
    
  }//end of 3 

  for(int iftrk=0; iftrk<nfvtxs; iftrk++){
    int ifpt_cor = 0;

    if(dbbcs_buff[icent][ivz]>0){
      for(int ibuff=0;ibuff<dbbcs_buff[icent][ivz];ibuff++){
        for(int ipmt_buff=0; ipmt_buff<buff_nbbcs[icent][ivz][ibuff]; ipmt_buff++){
          float val = bbcsbuff_et[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ipmt_buff]*fvtxs_et[iftrk];
          float dphi = fvtxs_phi[iftrk] - bbcsbuff_phi[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ipmt_buff];
          //float deta = fvtxs_eta[iftrk] - bbcsbuff_eta[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ipmt_buff];
          //int sign =(int)(gRandom->Rndm()*2);

          if(dphi<-0.5*pi) dphi+=2*pi;
          if(dphi>1.5*pi) dphi-=2*pi;

          float dphia = dphi;
          //if(sign) {
	  //if(dphi<0.5*pi) dphia  = -dphi;
	  //else dphia = 2*pi - dphi;
          //}

	  if(sign) dphia = -dphia;
	  if(dphia<-0.5*pi) dphia+=2*pi;
	  if(dphia>1.5*pi) dphia-=2*pi;

          hbacksnbbc2[icent][ifpt_cor]->Fill(dphi);
          kbacksnbbc2[icent][ifpt_cor]->Fill(dphia);
//          kbacksouthetabbc2[icent][ipt_cor]->Fill(deta, dphia);

          hbacksnbbcw2[icent][ifpt_cor]->Fill(dphi,val*val);
          kbacksnbbcw2[icent][ifpt_cor]->Fill(dphia,val);
//          kbacksouthetabbcw2[icent][ipt_cor]->Fill(deta, dphia, val);
        }
      }
    }
  }

  nevent++;

  //4) buff track
  if(ntrack>0){//rebuff
    
    int ibuff=ntrack_buff[icent][ivz];
    
    int itrk_buff=0;//track number in each buff

    if(ntrack>0){//buffing lambda
      for(int itrk=0;itrk<ntrack;itrk++){
	
        trackbuff_pt[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff]=track_pt[itrk];
	trackbuff_phi[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff]=track_phi[itrk];
	trackbuff_eta[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff]=track_eta[itrk];
        itrk_buff++;
      }
      if(itrk_buff>0){
        buff_ntrack[icent][ivz][ibuff]=itrk_buff;
        ntrack_buff[icent][ivz]++;
        if(dtrack_buff[icent][ivz]<nbuff) dtrack_buff[icent][ivz]++;
        if(ntrack_buff[icent][ivz]>=nbuff) ntrack_buff[icent][ivz]=0;
      }
    }
  }//end of 4

  //6) buff bbcs
  if(nbbcs>0){//rebuff
    int ibuff=nbbcs_buff[icent][ivz];
    
    int ibbcs_buff=0;//bbcs number in each buff

    if(nbbcs>0){//buffing lambda
      for(int ibbc=0;ibbc<nbbcs;ibbc++){
	bbcsbuff_et[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ibbcs_buff]=bbcs_et[ibbc];
	bbcsbuff_phi[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ibbcs_buff]=bbcs_phi[ibbc];
        bbcsbuff_eta[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ibbcs_buff]=bbcs_eta[ibbc];

        ibbcs_buff++;
      }
      if(ibbcs_buff>0){
        buff_nbbcs[icent][ivz][ibuff]=ibbcs_buff;
        nbbcs_buff[icent][ivz]++;
        if(dbbcs_buff[icent][ivz]<nbuff) dbbcs_buff[icent][ivz]++;
        if(nbbcs_buff[icent][ivz]>=nbuff) nbbcs_buff[icent][ivz]=0;
      }
    }
  }

  if(nfvtxs>0){//rebuff
    int ibuff=nbbcn_buff[icent][ivz];

    int ibbcn_buff=0;//bbcn number in each buff

    if(nfvtxs>0){//buffing lambda
      for(int ibbc=0;ibbc<nfvtxs;ibbc++){
        bbcnbuff_et[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ibbcn_buff]=fvtxs_et[ibbc];
        bbcnbuff_phi[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ibbcn_buff]=fvtxs_phi[ibbc];
        bbcnbuff_eta[icent*vzbin*nbuff*nbbc+ivz*nbuff*nbbc+ibuff*nbbc+ibbcn_buff]=fvtxs_eta[ibbc];

        ibbcn_buff++;
      }
      if(ibbcn_buff>0){
        buff_nbbcn[icent][ivz][ibuff]=ibbcn_buff;
        nbbcn_buff[icent][ivz]++;
        if(dbbcn_buff[icent][ivz]<nbuff) dbbcn_buff[icent][ivz]++;
        if(nbbcn_buff[icent][ivz]>=nbuff) nbbcn_buff[icent][ivz]=0;
      }
    }
  }

  }//event loop

return 0;
}   
//_____________________________________________________________________________________________________________________________
int RidgedAuRun16::End()
{
  cout << "End of RidgedAuRun16 for Run " << RunNumber << endl;
  cout << "Total # of events = " << ievent << " " << jevent << " " << nevent << endl;
  cout << "OutputFileName = " << OutputFileName << endl;

  if(d_outfile) {
   
    d_outfile->cd();
    //  Here I will write the output file...

    hcentbbcs->Write();
    hcentnfvtxs->Write();
    for(int icent=0; icent<centbin; icent++){
    hrunbbcs[icent]->Write();
    hrunnfvtxs[icent]->Write();
    hrunntrack[icent]->Write();
    hbbcsnfvtxs[icent]->Write();
      for(int ipt=0; ipt<ptbin; ipt++){
       //BBC
       //normal
        hforesouthbbc[icent][ipt]->Write();
        hbacksouthbbc[icent][ipt]->Write();
        //flip
        kforesouthbbc[icent][ipt]->Write();
        kbacksouthbbc[icent][ipt]->Write();
        //with weight
        hforesouthbbcw[icent][ipt]->Write();
        hbacksouthbbcw[icent][ipt]->Write();
        //flip
        kforesouthbbcw[icent][ipt]->Write();
        kbacksouthbbcw[icent][ipt]->Write();
        //tower mixing
        hbacksouthbbc2[icent][ipt]->Write();
        kbacksouthbbc2[icent][ipt]->Write();
        
        hbacksouthbbcw2[icent][ipt]->Write();
        kbacksouthbbcw2[icent][ipt]->Write();
  //north
        hforenorthbbc[icent][ipt]->Write();
        hbacknorthbbc[icent][ipt]->Write();

        kforenorthbbc[icent][ipt]->Write();
        kbacknorthbbc[icent][ipt]->Write();
  
        hforenorthbbcw[icent][ipt]->Write();
        hbacknorthbbcw[icent][ipt]->Write();

        kforenorthbbcw[icent][ipt]->Write();
        kbacknorthbbcw[icent][ipt]->Write();
        
        hbacknorthbbc2[icent][ipt]->Write();
        kbacknorthbbc2[icent][ipt]->Write();
        
        hbacknorthbbcw2[icent][ipt]->Write();
        kbacknorthbbcw2[icent][ipt]->Write();
  //south-north
        if(ipt==0){
        hforesnbbc[icent][ipt]->Write();
        hbacksnbbc[icent][ipt]->Write();

        kforesnbbc[icent][ipt]->Write();
        kbacksnbbc[icent][ipt]->Write();
  
        hforesnbbcw[icent][ipt]->Write();
        hbacksnbbcw[icent][ipt]->Write();

        kforesnbbcw[icent][ipt]->Write();
        kbacksnbbcw[icent][ipt]->Write();
        
        hbacksnbbc2[icent][ipt]->Write();
        kbacksnbbc2[icent][ipt]->Write();
        
        hbacksnbbcw2[icent][ipt]->Write();
        kbacksnbbcw2[icent][ipt]->Write();
      }
    }
    }
/*
    for(int iarm=0;iarm<2;iarm++){
      DCAxydis[iarm]->Write();
      DCAxy2dis[iarm]->Write();
      DCAcentdis[iarm]->Write();
    }
      //  TH1F *fvtxdphidis->Write();
      //  TH1F *fvtxdphidis2->Write();
    hvtx0etaz->Write();
    hvtx1etaz->Write();

    hvtx0etaphi->Write();
    hvtx1etaphi->Write();
*/

    //htree->Write();
    
    d_outfile->Close();
    //delete d_outfile;
    
  }else{

    cout<<"ERROR: No output file set!"<<endl;

  }

  return 0;
}
