#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TMath.h>
#include <TString.h>
#include <TLatex.h>
#include "TComplex.h"
#include "QCumulantOnline.h"
#include <iostream>

using namespace std;
int get_fvtx_layer(float);
void initialize_pmt_position();
float d_pmt_x[64];
float d_pmt_y[64];
float d_pmt_z = -1443.5; // same for all tubes

QCumulantOnline::QCumulantOnline(std::vector<TString> input, const char* output) :
  OutputFileName(output), InputFileName(input), ievent(0), jevent(0), calFlag(0)
{
  d_outfile = NULL;
  isweight = true;
  for(int icent=0;icent<ncent;icent++){
    for(int isub=0;isub<nsub;isub++){
      for(int ihar=0;ihar<nhar;ihar++){
        hcos1[icent][isub][ihar]=NULL;
        hsin1[icent][isub][ihar]=NULL;
        hcos1p2[icent][isub][ihar]=NULL; 
        hsin1p2[icent][isub][ihar]=NULL;
        hcos1m2m3[icent][isub][ihar]=NULL; 
        hsin1m2m3[icent][isub][ihar]=NULL;

        hwcos1[icent][isub][ihar]=NULL;
        hwsin1[icent][isub][ihar]=NULL;
        hwcos1p2[icent][isub][ihar]=NULL; 
        hwsin1p2[icent][isub][ihar]=NULL;
        hwcos1m2m3[icent][isub][ihar]=NULL; 
        hwsin1m2m3[icent][isub][ihar]=NULL;

        for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
            phiweight[icent][ibbcz][ihar][isub]=NULL;
        }
        
        for(int ipt=0;ipt<npt;ipt++){
            hcos1pr[icent][isub][ihar][ipt]=NULL;
            hsin1pr[icent][isub][ihar][ipt]=NULL;
            hcos1p2pr[icent][isub][ihar][ipt]=NULL; 
            hsin1p2pr[icent][isub][ihar][ipt]=NULL;
            hcos1p2m3pr[icent][isub][ihar][ipt]=NULL; 
            hsin1p2m3pr[icent][isub][ihar][ipt]=NULL;
            hcos1m2m3pr[icent][isub][ihar][ipt]=NULL; 
            hsin1m2m3pr[icent][isub][ihar][ipt]=NULL;

            hwcos1pr[icent][isub][ihar][ipt]=NULL;
            hwsin1pr[icent][isub][ihar][ipt]=NULL;
            hwcos1p2pr[icent][isub][ihar][ipt]=NULL; 
            hwsin1p2pr[icent][isub][ihar][ipt]=NULL;
            hwcos1p2m3pr[icent][isub][ihar][ipt]=NULL; 
            hwsin1p2m3pr[icent][isub][ihar][ipt]=NULL;
            hwcos1m2m3pr[icent][isub][ihar][ipt]=NULL; 
            hwsin1m2m3pr[icent][isub][ihar][ipt]=NULL;
        }
        for(int icorr=0;icorr<ncorr;icorr++){
        hd[icent][isub][ihar][icorr]=NULL;
        hc[icent][isub][ihar][icorr]=NULL;
        hv[icent][isub][ihar][icorr]=NULL;
        hwd[icent][isub][ihar][icorr]=NULL;
        hwc[icent][isub][ihar][icorr]=NULL;
        hwv[icent][isub][ihar][icorr]=NULL;
        
        for(int ipt=0;ipt<npt;ipt++){
            hdpr[icent][isub][ihar][ipt][icorr]=NULL;
            hcpr[icent][isub][ihar][ipt][icorr]=NULL;
            hvpr[icent][isub][ihar][ipt][icorr]=NULL;
            hwdpr[icent][isub][ihar][ipt][icorr]=NULL;
            hwcpr[icent][isub][ihar][ipt][icorr]=NULL;
            hwvpr[icent][isub][ihar][ipt][icorr]=NULL;
        }
    }
    }
    }
  }
}

QCumulantOnline::~QCumulantOnline()
{
//  d_outfile->Close();
//  delete bbccalib;
//  delete bbcgeo;
  cout << " QCumulantOnline::~QCumulantOnline " << endl;
}

int QCumulantOnline::Inittree(){
  //------------------------------------------------------------//
  //               Initializing Tree Variables                  //
  //------------------------------------------------------------//
  tree = NULL;
  tree = new TChain("ntp_event");
  for(unsigned int itree=0;itree<InputFileName.size();itree++){
     cout<<InputFileName[itree]<<endl;
     try{
         tree->Add(InputFileName[itree]);}
     catch(int e){ cout << "An exception occurred. Exception Nr. " << e << '\n';}
  }

  cout << "Now getting ready to read in the tree branch addresses and stuff...." << endl;

  // List of branches
  TBranch* b_event;   //!
  TBranch* b_bbc_z;   //!
  TBranch* b_centrality;   //!
  TBranch* b_bbc_qn;   //!
  TBranch* b_bbc_qs;   //!
  TBranch* b_npc1;   //!
  TBranch* b_trigger_scaled;   //!
  TBranch* b_trigger_live;   //!
  TBranch* b_d_Qx;   //!
  TBranch* b_d_Qy;   //!
  TBranch* b_d_Qw;   //!
  TBranch* b_bc_x;   //!
  TBranch* b_bc_y;   //!
  TBranch* b_vtx_z;   //!
  TBranch* b_fvtx_x;   //!
  TBranch* b_fvtx_y;   //!
  TBranch* b_fvtx_z;   //!
  TBranch* b_d_BBC_charge;   //!
  TBranch* b_d_nFVTX_clus;   //!
  TBranch* b_d_nFVTXN_clus;   //!
  TBranch* b_d_nFVTXS_clus;   //!
  TBranch* b_d_FVTX_x;   //!
  TBranch* b_d_FVTX_y;   //!
  TBranch* b_d_FVTX_z;   //!
  TBranch* b_ntrk;   //!
  TBranch* b_px;   //!
  TBranch* b_py;   //!
  TBranch* b_pz;   //!
  TBranch* b_pc3sdphi;   //!
  TBranch* b_pc3sdz;   //!
  // TBranch* b_d_ntrk;   //!
  // TBranch* b_d_cntpx;   //!
  // TBranch* b_d_cntpy;   //!
  // TBranch* b_d_cntpz;   //!

  tree->SetBranchAddress("bbc_z",&d_bbcz,&b_bbc_z);
  tree->SetBranchAddress("centrality",&centrality,&b_centrality);
  tree->SetBranchAddress("bbc_qn",&bbc_qn,&b_bbc_qn);
  tree->SetBranchAddress("bbc_qs",&bbc_qs,&b_bbc_qs);
  tree->SetBranchAddress("npc1",&npc1,&b_npc1);
  tree->SetBranchAddress("event",&event,&b_event);
  tree->SetBranchAddress("trigger_scaled",&trigger_scaled,&b_trigger_scaled);
  tree->SetBranchAddress("trigger_live",&trigger_live,&b_trigger_live);
  tree->SetBranchAddress("bc_x",&bc_x,&b_bc_x);
  tree->SetBranchAddress("bc_y",&bc_y,&b_bc_y);
  tree->SetBranchAddress("vtx_z",&vtx_z,&b_vtx_z);

  tree->SetBranchAddress("fvtx_x",&eventfvtx_x,&b_fvtx_x);
  tree->SetBranchAddress("fvtx_y",&eventfvtx_y,&b_fvtx_y);
  tree->SetBranchAddress("fvtx_z",&eventfvtx_z,&b_fvtx_z);

  tree->SetBranchAddress("d_BBC_charge",d_BBC_charge,&b_d_BBC_charge);
  tree->SetBranchAddress("d_Qx",d_Qx,&b_d_Qx);
  tree->SetBranchAddress("d_Qy",d_Qy,&b_d_Qy);
  tree->SetBranchAddress("d_Qw",d_Qw,&b_d_Qw);

  tree->SetBranchAddress("d_nFVTX_clus",&d_nFVTX_clus,&b_d_nFVTX_clus);
  tree->SetBranchAddress("d_nFVTXN_clus",&d_nFVTXN_clus,&b_d_nFVTXN_clus);
  tree->SetBranchAddress("d_nFVTXS_clus",&d_nFVTXS_clus,&b_d_nFVTXS_clus);
  tree->SetBranchAddress("d_FVTX_x",d_FVTX_x,&b_d_FVTX_x);
  tree->SetBranchAddress("d_FVTX_y",d_FVTX_y,&b_d_FVTX_y);
  tree->SetBranchAddress("d_FVTX_z",d_FVTX_z,&b_d_FVTX_z);

  tree->SetBranchAddress("d_ntrk",&d_ntrk,&b_ntrk);
  tree->SetBranchAddress("d_cntpx",d_px,&b_px);
  tree->SetBranchAddress("d_cntpy",d_py,&b_py);
  tree->SetBranchAddress("d_cntpz",d_pz,&b_pz);
  tree->SetBranchAddress("d_cntpc3sdphi",d_pc3sdphi,&b_pc3sdphi);
  tree->SetBranchAddress("d_cntpc3sdz",d_pc3sdz,&b_pc3sdz);

  return 0;
}

int QCumulantOnline::Init(){
  initialize_pmt_position();
  float pi = acos(-1.0);
    TH1::SetDefaultSumw2();
    for(int icent=0;icent<ncent;icent++){
      for(int isub=0;isub<nsub;isub++){
        for(int ihar=0;ihar<nhar;ihar++){
        hcos1[icent][isub][ihar]=new TH1D(Form("hcos1_%d_%d_%d",icent,isub,ihar),Form("hcos1_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hsin1[icent][isub][ihar]=new TH1D(Form("hsin1_%d_%d_%d",icent,isub,ihar),Form("hsin1_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hcos1p2[icent][isub][ihar]=new TH1D(Form("hcos1p2_%d_%d_%d",icent,isub,ihar),Form("hcos1p2_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hsin1p2[icent][isub][ihar]=new TH1D(Form("hsin1p2_%d_%d_%d",icent,isub,ihar),Form("hsin1p2_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hcos1m2m3[icent][isub][ihar]=new TH1D(Form("hcos1m2m3_%d_%d_%d",icent,isub,ihar),Form("hcos1m2m3_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hsin1m2m3[icent][isub][ihar]=new TH1D(Form("hsin1m2m3_%d_%d_%d",icent,isub,ihar),Form("hsin1m2m3_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        
        hwcos1[icent][isub][ihar]=new TH1D(Form("hwcos1_%d_%d_%d",icent,isub,ihar),Form("hwcos1_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hwsin1[icent][isub][ihar]=new TH1D(Form("hwsin1_%d_%d_%d",icent,isub,ihar),Form("hwsin1_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hwcos1p2[icent][isub][ihar]=new TH1D(Form("hwcos1p2_%d_%d_%d",icent,isub,ihar),Form("hwcos1p2_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hwsin1p2[icent][isub][ihar]=new TH1D(Form("hwsin1p2_%d_%d_%d",icent,isub,ihar),Form("hwsin1p2_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hwcos1m2m3[icent][isub][ihar]=new TH1D(Form("hwcos1m2m3_%d_%d_%d",icent,isub,ihar),Form("hwcos1m2m3_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        hwsin1m2m3[icent][isub][ihar]=new TH1D(Form("hwsin1m2m3_%d_%d_%d",icent,isub,ihar),Form("hwsin1m2m3_%d_%d_%d",icent,isub,ihar),220,-1.1,1.1);
        for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
          phiweight[icent][ibbcz][ihar][isub]=new TH1D(Form("phiweight_%d_%d_%d_%d",icent,ibbcz,ihar,isub),Form("phiweight_%d_%d_%d_%d",icent,ibbcz,ihar,isub),50,-pi,pi);
        }
        for(int ipt=0;ipt<npt;ipt++){
        hcos1pr[icent][isub][ihar][ipt]=new TH1D(Form("hcos1pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hcos1pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hsin1pr[icent][isub][ihar][ipt]=new TH1D(Form("hsin1pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hsin1pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hcos1p2pr[icent][isub][ihar][ipt]=new TH1D(Form("hcos1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hcos1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hsin1p2pr[icent][isub][ihar][ipt]=new TH1D(Form("hsin1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hsin1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hcos1p2m3pr[icent][isub][ihar][ipt]=new TH1D(Form("hcos1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hcos1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hsin1p2m3pr[icent][isub][ihar][ipt]=new TH1D(Form("hsin1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hsin1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hcos1m2m3pr[icent][isub][ihar][ipt]=new TH1D(Form("hcos1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hcos1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hsin1m2m3pr[icent][isub][ihar][ipt]=new TH1D(Form("hsin1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hsin1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hwcos1pr[icent][isub][ihar][ipt]=new TH1D(Form("hwcos1pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hwcos1pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hwsin1pr[icent][isub][ihar][ipt]=new TH1D(Form("hwsin1pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hwsin1pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hwcos1p2pr[icent][isub][ihar][ipt]=new TH1D(Form("hwcos1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hwcos1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hwsin1p2pr[icent][isub][ihar][ipt]=new TH1D(Form("hwsin1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hwsin1p2pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hwcos1p2m3pr[icent][isub][ihar][ipt]=new TH1D(Form("hwcos1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hwcos1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hwsin1p2m3pr[icent][isub][ihar][ipt]=new TH1D(Form("hwsin1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hwsin1p2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hwcos1m2m3pr[icent][isub][ihar][ipt]=new TH1D(Form("hwcos1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hwcos1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        hwsin1m2m3pr[icent][isub][ihar][ipt]=new TH1D(Form("hwsin1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),Form("hwsin1m2m3pr_%d_%d_%d_%d",icent,isub,ihar,ipt),220,-1.1,1.1);
        }
          for(int icorr=0;icorr<ncorr;icorr++){
            hd[icent][isub][ihar][icorr]=new TH1D(Form("hd_%d_%d_%d_%d",icent,isub,ihar,icorr),Form("hd_%d_%d_%d_%d",icent,isub,ihar,icorr),220,-1.1,1.1);
            hc[icent][isub][ihar][icorr]=new TH1D(Form("hc_%d_%d_%d_%d",icent,isub,ihar,icorr),Form("hc_%d_%d_%d_%d",icent,isub,ihar,icorr),220,-1.1,1.1);
            hv[icent][isub][ihar][icorr]=new TH1D(Form("hv_%d_%d_%d_%d",icent,isub,ihar,icorr),Form("hv_%d_%d_%d_%d",icent,isub,ihar,icorr),220,-1.1,1.1);
            hwd[icent][isub][ihar][icorr]=new TH1D(Form("hwd_%d_%d_%d_%d",icent,isub,ihar,icorr),Form("hwd_%d_%d_%d_%d",icent,isub,ihar,icorr),220,-1.1,1.1);
            hwc[icent][isub][ihar][icorr]=new TH1D(Form("hwc_%d_%d_%d_%d",icent,isub,ihar,icorr),Form("hwc_%d_%d_%d_%d",icent,isub,ihar,icorr),220,-1.1,1.1);
            hwv[icent][isub][ihar][icorr]=new TH1D(Form("hwv_%d_%d_%d_%d",icent,isub,ihar,icorr),Form("hwv_%d_%d_%d_%d",icent,isub,ihar,icorr),220,-1.1,1.1);
        for(int ipt=0;ipt<npt;ipt++){
            hdpr[icent][isub][ihar][ipt][icorr]=new TH1D(Form("hdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),Form("hdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),220,-1.1,1.1);
            hcpr[icent][isub][ihar][ipt][icorr]=new TH1D(Form("hcpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),Form("hcpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),220,-1.1,1.1);
            hvpr[icent][isub][ihar][ipt][icorr]=new TH1D(Form("hvpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),Form("hvpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),220,-1.1,1.1);
            hwdpr[icent][isub][ihar][ipt][icorr]=new TH1D(Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),Form("hwdpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),220,-1.1,1.1);
            hwcpr[icent][isub][ihar][ipt][icorr]=new TH1D(Form("hwcpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),Form("hwcpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),220,-1.1,1.1);
            hwvpr[icent][isub][ihar][ipt][icorr]=new TH1D(Form("hwvpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),Form("hwvpr_%d_%d_%d_%d_%d",icent,isub,ihar,ipt,icorr),220,-1.1,1.1);
        }
          }
        }
      }
    }

  d_outfile = new TFile(OutputFileName.c_str(),"recreate");
  cout<<"finish of initialize"<<endl;
  return 0;
}

int QCumulantOnline::process_event(){
  int nEvent = tree->GetEntries();
  cout<<nEvent<<endl;
  for(ievent=0;ievent < nEvent; ievent++){
      tree->GetEntry(ievent);
  if(ievent%10000==0) {
    cout<<"QCumulant calFlag = "<< calFlag << "************* ievent= "<<ievent<<"    *************"<<endl;
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
  int cent = (int)centrality;
//  float bbc_s = d_global->getBbcChargeS();
  float vvertex = bbcv;
  /*
  PHPoint vertex = vtxout->get_Vertex();
  double vtxx = vertex.getX();
  double vtxy = vertex.getY();
  double vtxz = vertex.getZ();
  if(TMath::Abs(vtxx) > 10 || TMath::Abs(vtxy) > 10 || TMath::Abs(vtxz) > 30) continue;
  */
  float fvtxx = eventfvtx_x;
  float fvtxy = eventfvtx_y;
  float fvtxz = eventfvtx_z;
  if(fvtxz != fvtxz){
      fvtxz = -9999;
  }
  if(RunNumber>=455792 && RunNumber<=458167){ // low energies use fvtx as vertex if there is one
  if(fvtxz != -9999)
    vvertex = fvtxz;
  }

  if(fabs(vvertex)>10) continue;
 
  const int vzbin = 10;
  int ibbcz = -1;
  ibbcz = vzbin*(vvertex+10)/20;
  if(ibbcz<0||ibbcz>=vzbin) continue;

 int icent = -9999;
 if(cent<0) continue;
 if(cent<=5) icent = 0;
 else if(cent<=10) icent = 1;
 if(icent<0) continue;

 // --- all numbers from Darren 2016-06-23
      const float x_off = 0.3;
      const float beam_angle = 0.001;
      float vtx_z = vvertex;
      float vtx_x = x_off + atan(beam_angle)*vtx_z;
      float vtx_y = 0.02;

      // --- radius cut using FVTX coordinates
      if ( sqrt(pow(fvtxx-vtx_x,2.0) +  pow(fvtxy-vtx_y,2.0)) >= 0.15 ) continue;

      int nclus = 0;
      int ncluss = 0;
      int nclusn = 0;
    for(int iclus = 0; iclus < d_nFVTX_clus; iclus++){
        float fvtx_x      = d_FVTX_x[iclus];
        float fvtx_y      = d_FVTX_y[iclus];
        float fvtx_z      = d_FVTX_z[iclus];
    int iarm = fvtx_z<0?0:1;//get_fvtx_layer(fvtx_z);
    fvtx_x = fvtx_x-vtx_x;
    fvtx_y = fvtx_y-vtx_y;
    fvtx_z = fvtx_z-vtx_z;
    fvtx_x = fvtx_z*sin(-beam_angle) + fvtx_x*cos(-beam_angle);
    if( (fabs(fvtx_x)>999) ||(fabs(fvtx_y)>999) || (fabs(fvtx_z)>999)) continue;
      nclus++;
      if(iarm==0)
      ncluss++;
      else
      nclusn++;
    }
   if ( RunNumber >= 454774 && RunNumber <= 455639 ) if(nclus>800||ncluss>400||nclusn>400) continue;
   if ( RunNumber >= 455792 && RunNumber <= 456283 ) if(nclus>800||ncluss>400||nclusn>400) continue; 
   if ( RunNumber >= 456652 && RunNumber <= 457298 ) if(nclus>600||ncluss>300||nclusn>300) continue;
   if ( RunNumber >= 457634 && RunNumber <= 458167 ) if(nclus>800||ncluss>400||nclusn>400) continue;

  jevent++;

  int iharE=0;
  if(nhar==1 || nhar==2) iharE = 1.0;
  

 double Qx[nsub];
 double Qy[nsub];
 double Qxk3[nsub];
 double Qyk3[nsub];
 double Qx2[nsub];
 double Qy2[nsub];
 double Qx2k2[nsub];
 double Qy2k2[nsub];
 double Qx3[nsub];
 double Qy3[nsub];
 double Qw[nsub];
 double qx[npt];
 double qy[npt];
 double qxk2[npt];
 double qyk2[npt];
 double qx2[npt];
 double qy2[npt];
 double qx2k1[npt];
 double qy2k1[npt];
 double qw[npt];
 double QS1[nsub];
 double QS2[nsub];
 double QS3[nsub];
 double QS4[nsub];
 double qs1[npt];
 double qs2[npt];
 double qs3[npt];

for(int ihar=0;ihar<nhar;ihar++){
for(int isub=0;isub<nsub;isub++){
      Qx[isub] = 0;
      Qy[isub] = 0;
      Qxk3[isub] = 0;
      Qyk3[isub] = 0;
      Qx2[isub] = 0;
      Qx2[isub] = 0;
      Qx2k2[isub] = 0;
      Qx2k2[isub] = 0;
      Qy2[isub] = 0;
      Qx3[isub] = 0;
      Qy3[isub] = 0;
      Qw[isub] = 0;
      QS1[isub] = 0;
      QS2[isub] = 0;
      QS3[isub] = 0;
      QS4[isub] = 0;
    }
    for(int ipt=0;ipt<npt;ipt++){
        qx[ipt] = 0;
        qy[ipt] = 0;
        qxk2[ipt] = 0;
        qyk2[ipt] = 0;
        qx2[ipt] = 0;
        qy2[ipt] = 0;
        qx2k1[ipt] = 0;
        qy2k1[ipt] = 0;
        qw[ipt] = 0;
        qs1[ipt] = 0;
        qs2[ipt] = 0;
        qs3[ipt] = 0;
    }

  int n = ihar + 1 + iharE;

    for(int iclus = 0; iclus < d_nFVTX_clus; iclus++){
        float fvtx_x      = d_FVTX_x[iclus];
        float fvtx_y      = d_FVTX_y[iclus];
        float fvtx_z      = d_FVTX_z[iclus];
    int istation = get_fvtx_layer(fvtx_z);
    fvtx_x = fvtx_x-vtx_x;
    fvtx_y = fvtx_y-vtx_y;
    fvtx_z = fvtx_z-vtx_z;
    fvtx_x = fvtx_z*sin(-beam_angle) + fvtx_x*cos(-beam_angle);
    double weight = 1.0;
//    double fvtx_r = sqrt(pow(fvtx_x,2.0)+pow(fvtx_y,2.0));
    if( (fabs(fvtx_x)>999) ||(fabs(fvtx_y)>999) || (fabs(fvtx_z)>999)) continue;
//    double fvtx_the = atan2(fvtx_r,fvtx_z-vvertex);
    double fvtx_phi = atan2(fvtx_y,fvtx_x);
    if(fvtx_z>0) continue;
//    double fvtx_eta = -log(tan(0.5*fvtx_the));
    if(isweight && calFlag == 0){
        phiweight[icent][ibbcz][ihar][istation]->Fill(fvtx_phi);
        phiweight[icent][ibbcz][ihar][5]->Fill(fvtx_phi);
    }
    else{
        if(isweight){
        if(phiweight[icent][ibbcz][ihar][istation]->GetBinContent(phiweight[icent][ibbcz][ihar][istation]->FindBin(fvtx_phi))!=0){
        if(phiweight[icent][ibbcz][ihar][istation]->FindBin(fvtx_phi)>0 && phiweight[icent][ibbcz][ihar][istation]->FindBin(fvtx_phi)<=50)
        weight = phiweight[icent][ibbcz][ihar][istation]->Integral()/phiweight[icent][ibbcz][ihar][istation]->GetNbinsX()/phiweight[icent][ibbcz][ihar][istation]->GetBinContent(phiweight[icent][ibbcz][ihar][istation]->FindBin(fvtx_phi));
     //if(fabs(weight-1.)>0.2) weight = 0.;
        else weight = 0.;
    }
        else weight = 0.;
        }
    Qx[istation]+=weight*cos(n*fvtx_phi);
    Qy[istation]+=weight*sin(n*fvtx_phi);
    Qxk3[istation]+=weight*weight*weight*cos(n*fvtx_phi);
    Qyk3[istation]+=weight*weight*weight*sin(n*fvtx_phi);
    Qx2[istation]+=weight*cos(2*n*fvtx_phi);
    Qy2[istation]+=weight*sin(2*n*fvtx_phi);
    Qx2k2[istation]+=weight*weight*cos(2*n*fvtx_phi);
    Qy2k2[istation]+=weight*weight*sin(2*n*fvtx_phi);
    Qx3[istation]+=weight*cos(3*n*fvtx_phi);
    Qy3[istation]+=weight*sin(3*n*fvtx_phi);
    Qw[istation]+=weight;
    QS1[istation]+=weight;
    QS2[istation]+=weight*weight;
    QS3[istation]+=weight*weight*weight;
    QS4[istation]+=weight*weight*weight*weight;

        if(isweight){
        if(phiweight[icent][ibbcz][ihar][5]->GetBinContent(phiweight[icent][ibbcz][ihar][5]->FindBin(fvtx_phi))!=0){
        if(phiweight[icent][ibbcz][ihar][5]->FindBin(fvtx_phi)>0 && phiweight[icent][ibbcz][ihar][5]->FindBin(fvtx_phi)<=50)
        weight = phiweight[icent][ibbcz][ihar][5]->Integral()/phiweight[icent][ibbcz][ihar][5]->GetNbinsX()/phiweight[icent][ibbcz][ihar][5]->GetBinContent(phiweight[icent][ibbcz][ihar][5]->FindBin(fvtx_phi));
        else weight = 0.;
        //if(fabs(weight-1.)>0.2) weight = 0.;
        }
        else weight = 0.;
        }
    Qx[5]+=weight*cos(n*fvtx_phi);
    Qy[5]+=weight*sin(n*fvtx_phi);
    Qxk3[5]+=weight*weight*weight*cos(n*fvtx_phi);
    Qyk3[5]+=weight*weight*weight*sin(n*fvtx_phi);
    Qx2[5]+=weight*cos(2*n*fvtx_phi);
    Qy2[5]+=weight*sin(2*n*fvtx_phi);
    Qx2k2[5]+=weight*weight*cos(2*n*fvtx_phi);
    Qy2k2[5]+=weight*weight*sin(2*n*fvtx_phi);
    Qx3[5]+=weight*cos(3*n*fvtx_phi);
    Qy3[5]+=weight*sin(3*n*fvtx_phi);
    Qw[5]+=weight;
    QS1[5]+=weight;
    QS2[5]+=weight*weight;
    QS3[5]+=weight*weight*weight;
    QS4[5]+=weight*weight*weight*weight;
    }//calFlag
    }
  

  // beam beam counter (bbc) r.p. // -----------------------------------
    for(int ipmt = 0; ipmt < 64; ipmt++){
    float bbcx      = d_pmt_x[ipmt];
    float bbcy      = d_pmt_y[ipmt];
    float charge = d_BBC_charge[ipmt];
    if(charge <= 0) continue;
    float bbcz      = d_pmt_z;
    bbcx = bbcx - vtx_x*10.0;
    bbcy = bbcy - vtx_y*10.0;
    bbcz = bbcz - vtx_z*10.0;
    // --- rotation
    bbcx = bbcz*sin(-beam_angle) + bbcx*cos(-beam_angle);

    if (charge>0 && bbcz<0) {
      double phi=atan2(bbcy,bbcx);
      double weight = charge;
    if(isweight && calFlag == 0){
        phiweight[icent][ibbcz][ihar][4]->Fill(phi);
    }
    else{
        if(isweight){
        if(phiweight[icent][ibbcz][ihar][4]->GetBinContent(phiweight[icent][ibbcz][ihar][4]->FindBin(phi))!=0){
        if(phiweight[icent][ibbcz][ihar][4]->FindBin(phi)>0 && phiweight[icent][ibbcz][ihar][4]->FindBin(phi)<=50)
        weight = phiweight[icent][ibbcz][ihar][4]->Integral()/phiweight[icent][ibbcz][ihar][4]->GetNbinsX()/phiweight[icent][ibbcz][ihar][4]->GetBinContent(phiweight[icent][ibbcz][ihar][4]->FindBin(phi));
        else weight = 0.;
        //if(fabs(weight-1.)>0.2) weight = 0.;
        }
        else weight = 0.;
        weight = charge * weight;
        }
      Qx[4] += weight * cos(n*phi);
      Qy[4] += weight * sin(n*phi);
      Qxk3[4]+=weight*weight*weight*cos(n*phi);
      Qyk3[4]+=weight*weight*weight*sin(n*phi);
      Qx2[4] += weight*cos(2*n*phi);
      Qy2[4] += weight*sin(2*n*phi);
      Qx2k2[4]+=weight*weight*cos(2*n*phi);
      Qy2k2[4]+=weight*weight*sin(2*n*phi);
      Qx3[4] += weight*cos(3*n*phi);
      Qy3[4] += weight*sin(3*n*phi);
      Qw[4] += weight;
      QS1[4]+=weight;
      QS2[4]+=weight*weight;
      QS3[4]+=weight*weight*weight;
      QS4[4]+=weight*weight*weight*weight;
    }//calFlag
    }
  }
    
if(calFlag>0){
// Tracks 
    for(int itrk=0; itrk< d_ntrk; itrk++){
      float px    = d_px[itrk];
      float py    = d_py[itrk];
      float pz    = d_pz[itrk];
      int dcarm=0;
      if(px>0) dcarm=1;
      px = pz*sin(-beam_angle) + px*cos(-beam_angle);
      float pt = sqrt(px*px+py*py);
      float phi = atan2(py,px);
      double weight = 1.0;

      double sdphi = d_pc3sdphi[itrk];//calcsdphi(d_scnt->get_pc3dphi(),d_scnt->get_dcarm(),d_scnt->get_charge(),d_scnt->get_mom());
      double sdz = d_pc3sdz[itrk];//calcsdz(d_scnt->get_pc3dz(),d_scnt->get_dcarm(),d_scnt->get_charge(),d_scnt->get_mom());
      if(fabs(sdphi)<2.0 && fabs(sdz)<2.0){
        if(isweight){
            weight = pt;
        }
            Qx[6] += weight*cos(n*phi);  
            Qy[6] += weight*sin(n*phi);       
            Qxk3[6] += weight*weight*weight*cos(n*phi);  
            Qyk3[6] += weight*weight*weight*sin(n*phi);       
            Qx2[6] += weight*cos(2*n*phi);
            Qy2[6] += weight*sin(2*n*phi);
            Qx2k2[6] += weight*weight*cos(2*n*phi);
            Qy2k2[6] += weight*weight*sin(2*n*phi);
            Qx3[6] += weight*cos(3*n*phi);
            Qy3[6] += weight*sin(3*n*phi);
            Qw[6] += weight;
            QS1[6]+=weight;
            QS2[6]+=weight*weight;
            QS3[6]+=weight*weight*weight;
            QS4[6]+=weight*weight*weight*weight;
            for(int ipt=0;ipt<npt;ipt++){
            if(pt>=ptbin[ipt] && pt<ptbin[ipt+1]){
            qx[ipt] += 1.0 * cos(n*phi); //particle of interest
            qy[ipt] += 1.0 * sin(n*phi); //particle of interest
            qxk2[ipt] += weight * weight * cos(n*phi); //POI & RP
            qyk2[ipt] += weight * weight * sin(n*phi); //POI & RP
            qx2[ipt] += 1.0 * cos(2*n*phi); //particle of interest
            qy2[ipt] += 1.0 * sin(2*n*phi); //particle of interest
            qx2k1[ipt] += weight * cos(2*n*phi); //POI & RP
            qy2k1[ipt] += weight * sin(2*n*phi); //POI & RP
            qw[ipt] += 1.0; //particle of interest
            qs1[ipt] += weight;
            qs2[ipt] += weight*weight;
            qs3[ipt] += weight*weight*weight;
            }
            }
      }
  }

for(int isub=0;isub<nsub;isub++){
//-------------------Integrated Qcumulant--------------------

  //Q vectors are ready
  TComplex Q = TComplex(Qx[isub],Qy[isub]); //Qn
  TComplex Qstar = TComplex(Qx[isub],-Qy[isub]); //Qn*
  TComplex Qk3 = TComplex(Qxk3[isub],Qyk3[isub]); //Qn
  TComplex Qk3star = TComplex(Qxk3[isub],-Qyk3[isub]); //Qn
  TComplex Q2 = TComplex(Qx2[isub],Qy2[isub]); //Q2n
  TComplex Q2star = TComplex(Qx2[isub],-Qy2[isub]); //Q2n*
  TComplex Q2k2 = TComplex(Qx2k2[isub],Qy2k2[isub]); //Q2n
  TComplex Q2k2star = TComplex(Qx2k2[isub],-Qy2k2[isub]); //Q2n
  TComplex Q3 = TComplex(Qx3[isub],Qy3[isub]); //Q3n
  double QW = Qw[isub];
  double S1 = QS1[isub]; 
  double S2 = QS2[isub]; 
  double S3 = QS3[isub]; 
  double S4 = QS4[isub]; 
//  if(QW<2) continue;
  
//  double W1 = QW;
//  double W2 = QW*(QW-1);
//  double W3 = QW*(QW-1)*(QW-2);
//  double W4 = QW*(QW-1)*(QW-2)*(QW-3);
  double W6 = QW*(QW-1)*(QW-2)*(QW-3)*(QW-4)*(QW-5);
  
  double M1 = S1;
  double M11 = S1*S1-S2;
  double M111 = S1*S1*S1-3.*S2*S1+2.*S3;
  double M1111 = S1*S1*S1*S1-6.*S2*S1*S1+8.*S3*S1+3.*S2*S2-6.*S4;

  //---2nd order Q-cumulant----
  double sqrt2Q = Q.Rho2(); //|Qn|^2
  
//  double d2 = W2>0?(sqrt2Q-QW)/(QW*(QW-1)):-9999; //<2>(n|n)
  double d2 = M11!=0?(sqrt2Q-S2)/M11:-9999; //<2>(n|n)
  
  //Detector non-uniform acceptance
  //double cos1 = Q.Re()/QW;
  //double sin1 = Q.Im()/QW;
  double cos1 = M1!=0?Q.Re()/M1:-9999;
  double sin1 = M1!=0?Q.Im()/M1:-9999;
  
  //---4st order Q-cumulant----
  double sqrt4Q = Q.Rho2()*Q.Rho2(); //|Qn|^4
  double sqrt2Q2 = Q2.Rho2(); //|Q2n|^2
  //double d4 = W4>0?(sqrt4Q+sqrt2Q2-2.*(Q2*Qstar*Qstar).Re()-4.*(QW-2)*sqrt2Q)/(QW*(QW-1)*(QW-2)*(QW-3))+2./((QW-1)*(QW-2)):-9999; //<4>(n,n|n,n)
  double d4 = M1111!=0?(sqrt4Q+Q2k2.Rho2()-2.*(Q2k2*Qstar*Qstar).Re()+8.*(Qk3*Qstar).Re()-4.*S2*sqrt2Q-6.*S4+2.*S2*S2)/M1111:-9999; //<4>(n,n|n,n)
  
  //Detector non-uniform acceptance
  //double cos1p2 = W2>0?(Q*Q-Q2).Re()/(QW*(QW-1)):-9999;
  //double sin1p2 = W2>0?(Q*Q-Q2).Im()/(QW*(QW-1)):-9999;
  double cos1p2 = M11!=0?(Q*Q-Q2k2).Re()/M11:-9999;
  double sin1p2 = M11!=0?(Q*Q-Q2k2).Im()/M11:-9999;
  //double cos1m2m3 = W3>0?(Q*Qstar*Qstar-Q*Q2star-2.*(QW-1)*Qstar).Re()/(QW*(QW-1)*(QW-2)):-9999;
  //double sin1m2m3 = W3>0?(Q*Qstar*Qstar-Q*Q2star-2.*(QW-1)*Qstar).Im()/(QW*(QW-1)*(QW-2)):-9999;
  double cos1m2m3 = M111!=0?(Q*Qstar*Qstar-Q*Q2k2star-2.*S2*Qstar+2.*Qk3star).Re()/M111:-9999;
  double sin1m2m3 = M111!=0?(Q*Qstar*Qstar-Q*Q2k2star-2.*S2*Qstar+2.*Qk3star).Im()/M111:-9999;

  //---6st order Q-cumulant----
  double sqrt6Q = Q.Rho2()*Q.Rho2()*Q.Rho2(); //|Qn|^6
  double sqrt2Q3 = Q3.Rho2(); //|Qn|^6
  double d6 = W6>0?(sqrt6Q+9.*sqrt2Q2*sqrt2Q-6.*(Q2*Q*Qstar*Qstar*Qstar).Re()+4.*(Q3*Qstar*Qstar*Qstar).Re()-12.*(Q3*Q2star*Qstar).Re()+18.*(QW-4)*(Q2*Qstar*Qstar).Re()+4.*sqrt2Q3)/(QW*(QW-1)*(QW-2)*(QW-3)*(QW-4)*(QW-5))-9.*(sqrt4Q+sqrt2Q2)/(QW*(QW-1)*(QW-2)*(QW-3)*(QW-5))+18.*sqrt2Q/(QW*(QW-1)*(QW-3)*(QW-4))-6./((QW-1)*(QW-2)*(QW-3)):-9999; //<6>(n,n,n|n,n,n)
  

  //----correlators------
  //double c2 = d2; //cn{2}
  //double c4 = d4-2*d2*d2; //cn{4}
  //double c6 = d6-9*d2*d4+12*c2*c2*c2; //cn{6}
  double c2 = d2-cos1*cos1-sin1*sin1;
  double c4 = d4-2*d2*d2-4*cos1*cos1m2m3+4*sin1*sin1m2m3-cos1p2*cos1p2-sin1p2*sin1p2+4*cos1p2*(cos1*cos1-sin1*sin1)+8*sin1p2*sin1*cos1+8*d2*(cos1*cos1+sin1*sin1)-6*(cos1*cos1+sin1*sin1)*(cos1*cos1+sin1*sin1);
  double c6 = d6-9*d2*d4+12*c2*c2*c2; //cn{6}
  
  //---flow harmonics------
  double v2 = c2>=0?sqrt(c2):-9999; //vn{2}
  double v4 = c4<=0?TMath::Power(-c4,1./4):-9999; //vn{4}
  double v6 = c6>=0?TMath::Power(1./4*c6,1./6):-9999; //vn{6}

  //----event histograms-----
  hd[icent][isub][ihar][0]->Fill(d2);
  hd[icent][isub][ihar][1]->Fill(d4);
  hd[icent][isub][ihar][2]->Fill(d6);
  hc[icent][isub][ihar][0]->Fill(c2);
  hc[icent][isub][ihar][1]->Fill(c4);
  hc[icent][isub][ihar][2]->Fill(c6);
  hv[icent][isub][ihar][0]->Fill(v2);
  hv[icent][isub][ihar][1]->Fill(v4);
  hv[icent][isub][ihar][2]->Fill(v6);

  hcos1[icent][isub][ihar]->Fill(cos1);
  hsin1[icent][isub][ihar]->Fill(sin1);
  hcos1p2[icent][isub][ihar]->Fill(cos1p2);
  hsin1p2[icent][isub][ihar]->Fill(sin1p2);
  hcos1m2m3[icent][isub][ihar]->Fill(cos1m2m3);
  hsin1m2m3[icent][isub][ihar]->Fill(sin1m2m3);
/*
  hwd[icent][isub][ihar][0]->Fill(d2,W2);
  hwd[icent][isub][ihar][1]->Fill(d4,W4);
  hwd[icent][isub][ihar][2]->Fill(d6,W6);
  hwc[icent][isub][ihar][0]->Fill(c2,W2);
  hwc[icent][isub][ihar][1]->Fill(c4,W4);
  hwc[icent][isub][ihar][2]->Fill(c6,W6);
  hwv[icent][isub][ihar][0]->Fill(v2,W2);
  hwv[icent][isub][ihar][1]->Fill(v4,W4);
  hwv[icent][isub][ihar][2]->Fill(v6,W6);
  */
  hwd[icent][isub][ihar][0]->Fill(d2,M11);
  hwd[icent][isub][ihar][1]->Fill(d4,M1111);
  hwd[icent][isub][ihar][2]->Fill(d6,W6);
  hwc[icent][isub][ihar][0]->Fill(c2,M11);
  hwc[icent][isub][ihar][1]->Fill(c4,M1111);
  hwc[icent][isub][ihar][2]->Fill(c6,W6);
  hwv[icent][isub][ihar][0]->Fill(v2,M11);
  hwv[icent][isub][ihar][1]->Fill(v4,M1111);
  hwv[icent][isub][ihar][2]->Fill(v6,W6);
  /*
  hwcos1[icent][isub][ihar]->Fill(cos1,W1);
  hwsin1[icent][isub][ihar]->Fill(sin1,W1);
  hwcos1p2[icent][isub][ihar]->Fill(cos1p2,W2);
  hwsin1p2[icent][isub][ihar]->Fill(sin1p2,W2);
  hwcos1m2m3[icent][isub][ihar]->Fill(cos1m2m3,W3);
  hwsin1m2m3[icent][isub][ihar]->Fill(sin1m2m3,W3);
  */
  hwcos1[icent][isub][ihar]->Fill(cos1,M1);
  hwsin1[icent][isub][ihar]->Fill(sin1,M1);
  hwcos1p2[icent][isub][ihar]->Fill(cos1p2,M11);
  hwsin1p2[icent][isub][ihar]->Fill(sin1p2,M11);
  hwcos1m2m3[icent][isub][ihar]->Fill(cos1m2m3,M111);
  hwsin1m2m3[icent][isub][ihar]->Fill(sin1m2m3,M111);

  for(int ipt=0;ipt<npt;ipt++){
//-------------------Differential Qcumulant-------------------
//q vectors are ready
  TComplex q = TComplex(qx[ipt],qy[ipt]); //qn
  TComplex qstar = TComplex(qx[ipt],-qy[ipt]); //qn*
  TComplex qk2 = TComplex(qxk2[ipt],qyk2[ipt]); //qn
  TComplex qk2star = TComplex(qxk2[ipt],-qyk2[ipt]); //qn
  TComplex q2 = TComplex(qx2[ipt],qy2[ipt]); //q2n
  TComplex q2k1 = TComplex(qx2k1[ipt],qy2k1[ipt]); //q2n
  TComplex q2star = TComplex(qx2[ipt],-qy2[ipt]); //q2n
 // TComplex q3 = TComplex(qx3,qy3); //q3n
  double qW = qw[ipt];
//  if(qW<1) continue;

  TComplex p=0; //overlap
  TComplex pstar=0; //overlap
  TComplex pk2;
  TComplex pk2star;
  TComplex p2k1;
  TComplex p2=0; //overlap
  double s1=0;
  double s2=0;
  double s3=0;
  double pW=0;

  if(isub == 6){
    p = q;
    pstar = qstar;
    p2 = q2;
    pk2 = qk2;
    pk2star = qk2star;
    p2k1 = q2k1;
    pW = qW;
    s1 = qs1[ipt];
    s2 = qs2[ipt];
    s3 = qs3[ipt];
  }
  else{
    p = 0;
    pstar = 0;
    p2 = 0;
    pk2 = 0;
    pk2star = 0;
    p2k1 = 0;
    pW = 0;
    s1 = 0;
    s2 = 0;
    s3 = 0;
  }
  
  double w1 = qW;
//  double w2 = qW*QW-pW;
//  double w3 = (qW*QW-2.*pW)*(QW-1);
//  double w4 = (qW*QW-3.*pW)*(QW-1)*(QW-2);
 
  double m01 = qW*S1-s1;
  double m011 = qW*(S1*S1-S2)-2.*(s1*S1-s2);
  double m0111 = qW*(S1*S1*S1-3.*S1*S2+2.*S3)-3.*(s1*(S1*S1-S2)+2.*(s3-s2*S1));


  //---2nd order Q-cumulant----
  //double d2pr = w2>0?((q*Qstar-pW)/(qW*QW-pW)).Re():-9999; //<2'>(n_|n__)
  double d2pr = m01!=0?((q*Qstar-s1)/m01).Re():-9999; //<2'>(n_|n__)
  
  //Detector non-uniform acceptance
  double cos1pr = q.Re()/qW;
  double sin1pr = q.Im()/qW;
  
  //---4st order Q-cumulant----
//  double d4pr = ((q*Q*Qstar*Qstar-q2*Qstar*Qstar-q*Q*Q2star+q2*Q2star-2.*qW*sqrt2Q-2.*(QW-3)*q*Qstar+2.*Q*qstar+2.*qW*(QW-3))/(qW*(QW-1)*(QW-2)*(QW-3))).Re(); //<4'>(n_n__|n__n__)
  //double d4pr = w4>0?((q*Q*Qstar*Qstar-p2*Qstar*Qstar-q*Q*Q2star+p2*Q2star-2.*pW*sqrt2Q-2.*QW*q*Qstar+7.*p*Qstar-Q*pstar+2.*q*Qstar+2.*pW*(QW-3))/((QW-1)*(QW-2)*(qW*QW-3.*pW))).Re():-9999; //<4'>(n_n__|n__n__)
  double d4pr = m0111!=0?((q*Q*Qstar*Qstar-p2k1*Qstar*Qstar-q*Q*Q2k2star+p2k1*Q2k2star-2.*s1*sqrt2Q-2.*S2*q*Qstar+7.*pk2*Qstar-Q*pk2star+2.*q*Qk3star+2.*s1*S2-6.*s3)/m0111).Re():-9999; //<4'>(n_n__|n__n__)
  
  //Detector non-uniform acceptance
  //double cos1p2pr = w2>0?(q*Q-p2).Re()/(qW*QW-pW):-9999;
  //double sin1p2pr = w2>0?(q*Q-p2).Im()/(qW*QW-pW):-9999;
  double cos1p2pr = m01!=0?(q*Q-p2k1).Re()/m01:-9999;
  double sin1p2pr = m01!=0?(q*Q-p2k1).Im()/m01:-9999;
  //double cos1p2m3pr = w3>0?(q*sqrt2Q-q*QW-p2*Qstar-pW*Q+2.*p).Re()/((qW*QW-2.*pW)*(QW-1)):-9999;
  //double sin1p2m3pr = w3>0?(q*sqrt2Q-q*QW-p2*Qstar-pW*Q+2.*p).Im()/((qW*QW-2.*pW)*(QW-1)):-9999;
  double cos1p2m3pr = m011!=0?(q*sqrt2Q-q*S2-p2k1*Qstar-s1*Q+2.*pk2).Re()/m011:-9999;
  double sin1p2m3pr = m011!=0?(q*sqrt2Q-q*S2-p2k1*Qstar-s1*Q+2.*pk2).Im()/m011:-9999;
  //double cos1m2m3pr = w3>0?(q*Qstar*Qstar-q*Q2star-2.*pW*Qstar+2.*pstar).Re()/((qW*QW-2.*pW)*(QW-1)):-9999;
  //double sin1m2m3pr = w3>0?(q*Qstar*Qstar-q*Q2star-2.*pW*Qstar+2.*pstar).Im()/((qW*QW-2.*pW)*(QW-1)):-9999;
  double cos1m2m3pr = m011>0?(q*Qstar*Qstar-q*Q2k2star-2.*s1*Qstar+2.*pk2star).Re()/m011:-9999;
  double sin1m2m3pr = m011>0?(q*Qstar*Qstar-q*Q2k2star-2.*s1*Qstar+2.*pk2star).Im()/m011:-9999;

  //----correlators------
  //double c2pr = d2pr;
  //double c4pr = d4pr-2*d2pr*d2;
    double c2pr = d2pr-cos1pr*cos1-sin1pr*sin1;
    double c4pr = d2pr-2*d2pr*d2
-cos1pr*cos1m2m3
+sin1pr*sin1m2m3
-cos1*cos1m2m3pr
+sin1*sin1m2m3pr
-2*cos1*cos1p2m3pr
-2*sin1*sin1p2m3pr
-cos1p2pr*cos1p2
-sin1p2pr*sin1p2
+2*cos1p2*(cos1pr*cos1-sin1pr*sin1)+2*sin1p2*(cos1pr*sin1+sin1pr*cos1)
+4*d2*(cos1pr*cos1+sin1pr*sin1)
+2*cos1p2pr*(cos1*cos1-sin1*sin1)
+4*sin1p2pr*cos1*sin1
+4*d2pr*(cos1*cos1+sin1*sin1)
-6*(cos1*cos1-sin1*sin1)*
 (cos1pr*cos1-sin1pr*sin1)
-12*cos1*sin1*
 (sin1pr*cos1-cos1pr*sin1);

  //----flow harmonics------
  double v2pr = c2pr/v2;
  double v4pr = -c4pr/v4/v4/v4;

  //----event histograms----------
  hdpr[icent][isub][ihar][ipt][0]->Fill(d2pr);
  hdpr[icent][isub][ihar][ipt][1]->Fill(d4pr);
  hcpr[icent][isub][ihar][ipt][0]->Fill(c2pr);
  hcpr[icent][isub][ihar][ipt][1]->Fill(c4pr);
  hvpr[icent][isub][ihar][ipt][0]->Fill(v2pr);
  hvpr[icent][isub][ihar][ipt][1]->Fill(v4pr);
  
  hcos1pr[icent][isub][ihar][ipt]->Fill(cos1pr);
  hsin1pr[icent][isub][ihar][ipt]->Fill(sin1pr);
  hcos1p2pr[icent][isub][ihar][ipt]->Fill(cos1p2pr);
  hsin1p2pr[icent][isub][ihar][ipt]->Fill(sin1p2pr);
  hcos1p2m3pr[icent][isub][ihar][ipt]->Fill(cos1p2m3pr);
  hsin1p2m3pr[icent][isub][ihar][ipt]->Fill(sin1p2m3pr);
  hcos1m2m3pr[icent][isub][ihar][ipt]->Fill(cos1m2m3pr);
  hsin1m2m3pr[icent][isub][ihar][ipt]->Fill(sin1m2m3pr);
/*
  hwdpr[icent][isub][ihar][ipt][0]->Fill(d2pr,w2);
  hwdpr[icent][isub][ihar][ipt][1]->Fill(d4pr,w4);
  hwcpr[icent][isub][ihar][ipt][0]->Fill(c2pr,w2);
  hwcpr[icent][isub][ihar][ipt][1]->Fill(c4pr,w4);
  hwvpr[icent][isub][ihar][ipt][0]->Fill(v2pr,w2);
  hwvpr[icent][isub][ihar][ipt][1]->Fill(v4pr,w4);
  */
  hwdpr[icent][isub][ihar][ipt][0]->Fill(d2pr,m01);
  hwdpr[icent][isub][ihar][ipt][1]->Fill(d4pr,m0111);
  hwcpr[icent][isub][ihar][ipt][0]->Fill(c2pr,m01);
  hwcpr[icent][isub][ihar][ipt][1]->Fill(c4pr,m0111);
  hwvpr[icent][isub][ihar][ipt][0]->Fill(v2pr,m01);
  hwvpr[icent][isub][ihar][ipt][1]->Fill(v4pr,m0111);
  /*
  hwcos1pr[icent][isub][ihar][ipt]->Fill(cos1pr,w1);
  hwsin1pr[icent][isub][ihar][ipt]->Fill(sin1pr,w1);
  hwcos1p2pr[icent][isub][ihar][ipt]->Fill(cos1p2pr,w2);
  hwsin1p2pr[icent][isub][ihar][ipt]->Fill(sin1p2pr,w2);
  hwcos1p2m3pr[icent][isub][ihar][ipt]->Fill(cos1p2m3pr,w3);
  hwsin1p2m3pr[icent][isub][ihar][ipt]->Fill(sin1p2m3pr,w3);
  hwcos1m2m3pr[icent][isub][ihar][ipt]->Fill(cos1m2m3pr,w3);
  hwsin1m2m3pr[icent][isub][ihar][ipt]->Fill(sin1m2m3pr,w3);
  */
  hwcos1pr[icent][isub][ihar][ipt]->Fill(cos1pr,w1);
  hwsin1pr[icent][isub][ihar][ipt]->Fill(sin1pr,w1);
  hwcos1p2pr[icent][isub][ihar][ipt]->Fill(cos1p2pr,m01);
  hwsin1p2pr[icent][isub][ihar][ipt]->Fill(sin1p2pr,m01);
  hwcos1p2m3pr[icent][isub][ihar][ipt]->Fill(cos1p2m3pr,m011);
  hwsin1p2m3pr[icent][isub][ihar][ipt]->Fill(sin1p2m3pr,m011);
  hwcos1m2m3pr[icent][isub][ihar][ipt]->Fill(cos1m2m3pr,m011);
  hwsin1m2m3pr[icent][isub][ihar][ipt]->Fill(sin1m2m3pr,m011);

  }//ipt
}//isub
}//calFlag
}//ihar

}//event loop
  return 0;
}

int QCumulantOnline::End()
{

  cout << "End of QCumulantOnline for Run " << RunNumber << endl;
  cout << "Total # of events = " << ievent << ", Total # of passed events = "<< jevent << endl;
  cout << "OutputFileName = " << OutputFileName << endl;
  
  if(d_outfile && d_outfile->IsOpen()) {
    d_outfile->cd();

  for(int icent=0;icent<ncent;icent++){
    for(int isub=0;isub<nsub;isub++){
      for(int ihar=0;ihar<nhar;ihar++){
        hcos1[icent][isub][ihar]->Write();
        hsin1[icent][isub][ihar]->Write();
//        hcos1p2[icent][isub][ihar]->Write(); 
//        hsin1p2[icent][isub][ihar]->Write();
//        hcos1m2m3[icent][isub][ihar]->Write(); 
//        hsin1m2m3[icent][isub][ihar]->Write();

        hwcos1[icent][isub][ihar]->Write();
        hwsin1[icent][isub][ihar]->Write();
        hwcos1p2[icent][isub][ihar]->Write(); 
        hwsin1p2[icent][isub][ihar]->Write();
        hwcos1m2m3[icent][isub][ihar]->Write(); 
        hwsin1m2m3[icent][isub][ihar]->Write();
        
        for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
            phiweight[icent][ibbcz][ihar][isub]->Write();
        }
        
        for(int ipt=0;ipt<npt;ipt++){
            /*
            hcos1pr[icent][isub][ihar][ipt]->Write();
            hsin1pr[icent][isub][ihar][ipt]->Write();
            hcos1p2pr[icent][isub][ihar][ipt]->Write(); 
            hsin1p2pr[icent][isub][ihar][ipt]->Write();
            hcos1p2m3pr[icent][isub][ihar][ipt]->Write(); 
            hsin1p2m3pr[icent][isub][ihar][ipt]->Write();
            hcos1m2m3pr[icent][isub][ihar][ipt]->Write(); 
            hsin1m2m3pr[icent][isub][ihar][ipt]->Write();
            */
            hwcos1pr[icent][isub][ihar][ipt]->Write();
            hwsin1pr[icent][isub][ihar][ipt]->Write();
            hwcos1p2pr[icent][isub][ihar][ipt]->Write(); 
            hwsin1p2pr[icent][isub][ihar][ipt]->Write();
            hwcos1p2m3pr[icent][isub][ihar][ipt]->Write(); 
            hwsin1p2m3pr[icent][isub][ihar][ipt]->Write();
            hwcos1m2m3pr[icent][isub][ihar][ipt]->Write(); 
            hwsin1m2m3pr[icent][isub][ihar][ipt]->Write();
        }
        for(int icorr=0;icorr<ncorr;icorr++){
//        hd[icent][isub][ihar][icorr]->Write();
//        hc[icent][isub][ihar][icorr]->Write();
//        hv[icent][isub][ihar][icorr]->Write();
        hwd[icent][isub][ihar][icorr]->Write();
        hwc[icent][isub][ihar][icorr]->Write();
        hwv[icent][isub][ihar][icorr]->Write();
        for(int ipt=0;ipt<npt;ipt++){
//            hdpr[icent][isub][ihar][ipt][icorr]->Write();
//            hcpr[icent][isub][ihar][ipt][icorr]->Write();
//            hvpr[icent][isub][ihar][ipt][icorr]->Write();
            hwdpr[icent][isub][ihar][ipt][icorr]->Write();
            hwcpr[icent][isub][ihar][ipt][icorr]->Write();
            hwvpr[icent][isub][ihar][ipt][icorr]->Write();
        }
    }
    }
    }
  }
  }
    return 0;
}

int QCumulantOnline::SetRun(int _run){
    RunNumber = _run;
    return 0;
  }

int QCumulantOnline::SetcalFlag(int _flag){
    calFlag = _flag;
    return 0;
  }

int QCumulantOnline::Setweight(bool _isweight){
    isweight = _isweight;
    return 0;
  }

int get_fvtx_layer(float z)
{
  // --- south side
  if ( z < -18 && z > -24 ) return 0;
  if ( z < -24 && z > -30 ) return 1;
  if ( z < -30 && z > -35 ) return 2;
  if ( z < -35 )            return 3;
  // --- north side
  if ( z > 18 && z < 24 ) return 0;
  if ( z > 24 && z < 30 ) return 1;
  if ( z > 30 && z < 35 ) return 2;
  if ( z > 35 )           return 3;
  // --- invalid numbers...
  cout<<"get_fvtx_layer::invalid z =  "<<z<<endl;
  return -1;
}

void initialize_pmt_position()
{

  d_pmt_x[0] = -123;
  d_pmt_y[0] = 42.6;
  d_pmt_x[1] = -123;
  d_pmt_y[1] = 14.2;
  d_pmt_x[2] = -98.4;
  d_pmt_y[2] = 85.2;
  d_pmt_x[3] = -98.4;
  d_pmt_y[3] = 56.8;
  d_pmt_x[4] = -98.4;
  d_pmt_y[4] = 28.4;
  d_pmt_x[5] = -73.8;
  d_pmt_y[5] = 99.4;
  d_pmt_x[6] = -73.8;
  d_pmt_y[6] = 71;
  d_pmt_x[7] = -73.8;
  d_pmt_y[7] = 42.6;
  d_pmt_x[8] = -73.8;
  d_pmt_y[8] = 14.2;
  d_pmt_x[9] = -49.2;
  d_pmt_y[9] = 113.6;
  d_pmt_x[10] = -49.2;
  d_pmt_y[10] = 85.2;
  d_pmt_x[11] = -49.2;
  d_pmt_y[11] = 56.8;
  d_pmt_x[12] = -24.6;
  d_pmt_y[12] = 127.8;
  d_pmt_x[13] = -24.6;
  d_pmt_y[13] = 99.4;
  d_pmt_x[14] = -24.6;
  d_pmt_y[14] = 71;
  d_pmt_x[15] = 0;
  d_pmt_y[15] = 113.6;
  d_pmt_x[16] = 0;
  d_pmt_y[16] = 85.2;
  d_pmt_x[17] = 24.6;
  d_pmt_y[17] = 127.8;
  d_pmt_x[18] = 24.6;
  d_pmt_y[18] = 99.4;
  d_pmt_x[19] = 24.6;
  d_pmt_y[19] = 71;
  d_pmt_x[20] = 49.2;
  d_pmt_y[20] = 113.6;
  d_pmt_x[21] = 49.2;
  d_pmt_y[21] = 85.2;
  d_pmt_x[22] = 49.2;
  d_pmt_y[22] = 56.8;
  d_pmt_x[23] = 73.8;
  d_pmt_y[23] = 99.4;
  d_pmt_x[24] = 73.8;
  d_pmt_y[24] = 71;
  d_pmt_x[25] = 73.8;
  d_pmt_y[25] = 42.6;
  d_pmt_x[26] = 73.8;
  d_pmt_y[26] = 14.2;
  d_pmt_x[27] = 98.4;
  d_pmt_y[27] = 85.2;
  d_pmt_x[28] = 98.4;
  d_pmt_y[28] = 56.8;
  d_pmt_x[29] = 98.4;
  d_pmt_y[29] = 28.4;
  d_pmt_x[30] = 123;
  d_pmt_y[30] = 42.6;
  d_pmt_x[31] = 123;
  d_pmt_y[31] = 14.2;
  d_pmt_x[32] = 123;
  d_pmt_y[32] = -42.6;
  d_pmt_x[33] = 123;
  d_pmt_y[33] = -14.2;
  d_pmt_x[34] = 98.4;
  d_pmt_y[34] = -85.2;
  d_pmt_x[35] = 98.4;
  d_pmt_y[35] = -56.8;
  d_pmt_x[36] = 98.4;
  d_pmt_y[36] = -28.4;
  d_pmt_x[37] = 73.8;
  d_pmt_y[37] = -99.4;
  d_pmt_x[38] = 73.8;
  d_pmt_y[38] = -71;
  d_pmt_x[39] = 73.8;
  d_pmt_y[39] = -42.6;
  d_pmt_x[40] = 73.8;
  d_pmt_y[40] = -14.2;
  d_pmt_x[41] = 49.2;
  d_pmt_y[41] = -113.6;
  d_pmt_x[42] = 49.2;
  d_pmt_y[42] = -85.2;
  d_pmt_x[43] = 49.2;
  d_pmt_y[43] = -56.8;
  d_pmt_x[44] = 24.6;
  d_pmt_y[44] = -127.8;
  d_pmt_x[45] = 24.6;
  d_pmt_y[45] = -99.4;
  d_pmt_x[46] = 24.6;
  d_pmt_y[46] = -71;
  d_pmt_x[47] = -0;
  d_pmt_y[47] = -113.6;
  d_pmt_x[48] = -0;
  d_pmt_y[48] = -85.2;
  d_pmt_x[49] = -24.6;
  d_pmt_y[49] = -127.8;
  d_pmt_x[50] = -24.6;
  d_pmt_y[50] = -99.4;
  d_pmt_x[51] = -24.6;
  d_pmt_y[51] = -71;
  d_pmt_x[52] = -49.2;
  d_pmt_y[52] = -113.6;
  d_pmt_x[53] = -49.2;
  d_pmt_y[53] = -85.2;
  d_pmt_x[54] = -49.2;
  d_pmt_y[54] = -56.8;
  d_pmt_x[55] = -73.8;
  d_pmt_y[55] = -99.4;
  d_pmt_x[56] = -73.8;
  d_pmt_y[56] = -71;
  d_pmt_x[57] = -73.8;
  d_pmt_y[57] = -42.6;
  d_pmt_x[58] = -73.8;
  d_pmt_y[58] = -14.2;
  d_pmt_x[59] = -98.4;
  d_pmt_y[59] = -85.2;
  d_pmt_x[60] = -98.4;
  d_pmt_y[60] = -56.8;
  d_pmt_x[61] = -98.4;
  d_pmt_y[61] = -28.4;
  d_pmt_x[62] = -123;
  d_pmt_y[62] = -42.6;
  d_pmt_x[63] = -123;
  d_pmt_y[63] = -14.2;

}



