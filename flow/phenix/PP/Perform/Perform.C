#include "Perform.h"
#include <stdlib.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include "TString.h"
#include <fstream>
#include <iostream>

using namespace std;

const float bbcz_cut = 10.0;
const int centbin = 1;
const int nbbcs = 800;
const int ncluss = 800;
const int nmpcs = 800;
const int nfvtxs = 800;
const float mpc_e_cut = 10.0;

const double mpion = 0.139570;
const double phbeta = 29.9792458;
float buff_bbcs = 0;
float buff_ZdcEs = 0;

Perform::Perform(TString input, TString output):
 OutputFileName(output), 
 InputFileName(input), 
  run(-9999),
  trig(-9999),
  cent(-9999),
  npc1hits(-9999),
  ZdcEs(-9999),
  ZdcEn(-9999),
  bbc_s(-9999),
  bbc_n(-9999),
  bbcv(-9999),
  vtxz(-9999),
  fvtxz(-9999),
  ntrack(-9999),
  nbbc(-9999),
  nvtx(-9999),
  nmpc(-9999),
  nvtxtrack(-9999),
  nfvtxtrack(-9999),

  //cnt track variables
  mom(),
  phi0(),
  the0(),
  phi(),
  alpha(),
  zed(),
  charge(),
  arm(),
  pc3dphi(),
  pc3dz(),
  
  //tof
  slat(),
  tofdphi(),
  tofdz(),
  tofpl(),
  qtof(),
  ttof(),

  //bbc
  bbccharge(),
  bbct0(),
  bbcx(),
  bbcy(),
  bbcz(),

  //vtx
  layer(),
  vtxX(),
  vtxY(),
  vtxZ(),

  //vtx track
  vtxnhits(),
  vtxpx(),
  vtxpy(),
  vtxpz(),

  //mpc
  mpc_e(),
  mpc_x(),
  mpc_y(),
  mpc_z(),

  //fvtx track
  farm(),
  fnhits(),
  fthe(),
  feta(),
  fphi(),
  fvtxX(),
  fvtxY(),
  fvtxZ()
{
}

Perform::~Perform()
{
    delete tree;
    d_infile->Close();
    d_outfile->Close();
}

int Perform::Init()
{
  TH1::SetDefaultSumw2();   
  coll = "";
  TString temp = OutputFileName;
  temp.ToLower();
  if(temp.Contains("pp") && temp.Contains("mb")) coll = "ppminbias";
  else if(temp.Contains("pal") && temp.Contains("mbst")) coll = "pAlminbias";
  else if(temp.Contains("pal") && temp.Contains("mbcentral")) coll = "pAlcentral";
  else if(temp.Contains("pp") && temp.Contains("fvtxand")) coll = "ppfvtxand";
  else if(temp.Contains("pal") && temp.Contains("fvtxand")) coll = "pAlfvtxand";
  else if(temp.Contains("pp") && temp.Contains("fvtxor")) coll = "ppfvtxor";
  else if(temp.Contains("pal") && temp.Contains("fvtxor")) coll = "pAlfvtxor";
  else if(temp.Contains("pal") && temp.Contains("fvtxsouth")) coll = "pAlfvtxsouth";
  else if(temp.Contains("pal") && temp.Contains("fvtxnorth")) coll = "pAlfvtxnorth";
  else return 0;

//  OutputFileName.Insert(OutputFileName.Length()-5,coll);

  d_outfile = new TFile(OutputFileName,"recreate");
  d_infile = TFile::Open(InputFileName,"ReadOnly");
  if(!d_infile) return 0;
  tree = (TTree*)d_infile->Get("tree");
  if(!tree) return 0;
  
  tree -> SetBranchAddress("run", &run);
  tree -> SetBranchAddress("trig", &trig);
  tree -> SetBranchAddress("npc1hits", &npc1hits);
  tree->SetBranchAddress("ZdcEs",&ZdcEs);
  tree->SetBranchAddress("ZdcEn",&ZdcEn);
  if(coll.Contains("pA"))
  tree->SetBranchAddress("cent",&cent);
  tree -> SetBranchAddress("bbc_s", &bbc_s);
  tree -> SetBranchAddress("bbc_n", &bbc_n);
  tree -> SetBranchAddress("bbcv", &bbcv);
  if(!coll.Contains("pAl")){
  tree -> SetBranchAddress("vtxz", &vtxz);
  }
//  tree -> SetBranchAddress("fvtxz", &fvtxz);
  tree -> SetBranchAddress("ntrack", &ntrack);
  tree -> SetBranchAddress("nbbc", &nbbc);
  if(!coll.Contains("pAl")){
  tree -> SetBranchAddress("nvtx", &nvtx);
  tree -> SetBranchAddress("nvtxtrack", &nvtxtrack);
  }
//  tree -> SetBranchAddress("nmpc", &nmpc);
//  tree -> SetBranchAddress("nfvtxtrack", &nfvtxtrack);

  tree -> SetBranchAddress("mom", &mom);
  tree -> SetBranchAddress("phi0", &phi0);
  tree -> SetBranchAddress("the0", &the0);
  tree -> SetBranchAddress("alpha", &alpha);
  tree -> SetBranchAddress("phi", &phi);
  tree -> SetBranchAddress("zed", &zed);
  tree -> SetBranchAddress("charge", &charge);
  tree -> SetBranchAddress("arm", &arm);
  tree -> SetBranchAddress("pc3dphi", &pc3dphi);
  tree -> SetBranchAddress("pc3dz", &pc3dz);
  tree -> SetBranchAddress("slat", &slat);
  tree -> SetBranchAddress("tofdphi", &tofdphi);
  tree -> SetBranchAddress("tofdz", &tofdz);
  tree -> SetBranchAddress("tofpl", &tofpl);
  tree -> SetBranchAddress("qtof", &qtof);
  tree -> SetBranchAddress("ttof", &ttof);

  tree -> SetBranchAddress("bbccharge", &bbccharge);
  tree -> SetBranchAddress("bbct0", &bbct0);
  tree -> SetBranchAddress("bbcx", &bbcx);
  tree -> SetBranchAddress("bbcy", &bbcy);
  tree -> SetBranchAddress("bbcz", &bbcz);

  if(!coll.Contains("pAl")){
  tree -> SetBranchAddress("layer", &layer);
  tree -> SetBranchAddress("vtxX", &vtxX);
  tree -> SetBranchAddress("vtxY", &vtxY);
  tree -> SetBranchAddress("vtxZ", &vtxZ);

  tree -> SetBranchAddress("vtxnhits", &vtxnhits);
  tree -> SetBranchAddress("vtxpx", &vtxpx);
  tree -> SetBranchAddress("vtxpy", &vtxpy);
  tree -> SetBranchAddress("vtxpz", &vtxpz);
  }
 /*
  tree -> SetBranchAddress("mpc_e", &mpc_e);
  tree -> SetBranchAddress("mpc_x", &mpc_x);
  tree -> SetBranchAddress("mpc_y", &mpc_y);
  tree -> SetBranchAddress("mpc_z", &mpc_z);
*/
  /*
  tree -> SetBranchAddress("farm", &farm);
  tree -> SetBranchAddress("fnhits", &fnhits);
  tree -> SetBranchAddress("fthe", &fthe);
  tree -> SetBranchAddress("feta", &feta);
  tree -> SetBranchAddress("fphi", &fphi);
  tree -> SetBranchAddress("fvtxX", &fvtxX);
  tree -> SetBranchAddress("fvtxY", &fvtxY);
  tree -> SetBranchAddress("fvtxZ", &fvtxZ);
*/
char name[512];
float pi = acos(-1);

//Tracks
sprintf(name,"hcntetaphi"); hcntetaphi = new TH2F(name,name,200,-1,1,200,-2*pi,2*pi);
sprintf(name,"hcntpt"); hcntpt = new TH1F(name,name,2000,0,10);
for(int i=0;i<50;i++){
sprintf(name,"pc3dphidz_arm0_pos_%d",i),pc3dphidz_arm0_pos[i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
sprintf(name,"pc3dphidz_arm1_pos_%d",i),pc3dphidz_arm1_pos[i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
sprintf(name,"pc3dphidz_arm0_neg_%d",i),pc3dphidz_arm0_neg[i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
sprintf(name,"pc3dphidz_arm1_neg_%d",i),pc3dphidz_arm1_neg[i] = new TH2F(name,name,200,-0.1,0.1,100,-10,10);
}

//tof
for(int ich=0;ich<2;ich++){
sprintf(name,"tofdphidz_%d",ich);  tofdphidz[ich] = new TH2F(name,name,200,-0.2,0.2,200,-10,10);
sprintf(name,"tofwdphidz_%d",ich);  tofwdphidz[ich] = new TH2F(name,name,200,-0.2,0.2,200,-10,10);
}
sprintf(name,"ttofqpratio");  ttofqpratio = new TH2F(name,name,500,0,100,500,-8,8);
sprintf(name,"m2qpratio");  m2qpratio = new TH2F(name,name,500,0,2,500,-8,8);
sprintf(name,"ttofp");  ttofp = new TH2F(name,name,500,0,100,500,0,10);
sprintf(name,"pinv2chbeta"); pinv2chbeta = new TH2F(name,name,500,0,10,500,-0.8,0.8);
sprintf(name,"m2p");  m2p = new TH2F(name,name,500,0,2,500,0,10);
sprintf(name,"deltattofeis"); deltattofeis = new TH2F(name,name,960,0,960,1600,-20,60);
sprintf(name,"deltattofwis"); deltattofwis = new TH2F(name,name,512,0,512,1600,-20,60);

//vtx
for(int i=0;i<4;i++){
sprintf(name,"hcluetaphi_%d",i); hcluetaphi[i] = new TH2F(name,name,300,-3.0,3.0,100,-4,4);
}

//bbc
sprintf(name,"bbcet"); bbcet = new TH1F(name,name,500,0,20);

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
  sprintf(name,"mpc_south_cent"); mpc_south_cent = new TH2F(name, name, 100, -0.5, 99.5, 250, -0.5, 249.5);
  sprintf(name,"mpc_north_cent"); mpc_north_cent = new TH2F(name, name, 100, -0.5, 99.5, 250, -0.5, 249.5);
  sprintf(name,"mpc_south_north"); mpc_south_north = new TH2F(name, name, 600, -0.5, 599.5, 600, -0.5, 599.5);
  sprintf(name,"south_mpc_north_mpc");  south_mpc_north_mpc = new TH2F(name, name, 600, -0.5, 19.5, 600, -0.5, 19.5);
for (int i=0; i<centbin; i++) {
 sprintf(name,"mpcetetasouth_%d", i);    mpcetetasouth[i] = new TH2F(name,name,30, 3.0, 4.0, 60,0.0,6.0);
 sprintf(name,"mpcetetanorth_%d", i);    mpcetetanorth[i] = new TH2F(name,name,30, 3.0, 4.0, 60,0.0,6.0);
}

//correlate
  sprintf(name,"hbbcsZdcEs"); hbbcsZdcEs = new TH2F(name,name,600,0,300,500,0,2500);
  sprintf(name,"hmixbbcsZdcEs"); hmixbbcsZdcEs = new TH2F(name,name,600,0,300,500,0,5000);
  sprintf(name,"hnpc3hitsntof"); hnpc3hitsntof = new TH2F(name,name,200,0,200,100,0,100);
  sprintf(name,"hbbcnbbc"); hbbcnbbc = new TH2F(name,name,600,0,300,600,0,600);
  sprintf(name,"hpc1hitsbbc"); hpc1hitsbbc = new TH2F(name, name, 100, 0, 100, 600, 0, 300);
  sprintf(name,"hvtxzfvtxz"); hvtxzfvtxz = new TH2F(name,name,100,-25.,25.,100,-25.,25.);
  for(int i=0;i<4;i++){
	sprintf(name,"hbbcsnvtx_%d",i);  hbbcsnvtx[i] = new TH2F(name,name,600,0,300,500,0,500);
	sprintf(name,"hbbcnnvtx_%d",i);  hbbcnnvtx[i] = new TH2F(name,name,600,0,300,500,0,500);
	sprintf(name,"hbbcnvtx_%d",i);  hbbcnvtx[i] = new TH2F(name,name,600,0,300,500,0,500);
        sprintf(name,"hnvtxnmpc_%d",i); hnvtxnmpc[i] = new TH2F(name,name,500,0,500,500,0,500);
        sprintf(name,"hnvtxntrk_%d",i); hnvtxntrk[i] = new TH2F(name,name,500,0,500,500,0,500);
  	sprintf(name,"hnvtxnfvtxtrk_%d",i); hnvtxnfvtxtrk[i] = new TH2F(name,name,500,0,500,500,0,500);
  }
  sprintf(name,"hnbbcnclu"); hnbbcnclu = new TH2F(name,name,600,0,300,500,0,500);
  sprintf(name,"hbbcsbbcn"); hbbcsbbcn = new TH2F(name,name,600,0,300,600,0,300);
  sprintf(name,"hntracknmpc"); hntracknmpc = new TH2F(name,name,100,0,100,500,0,500);
  sprintf(name,"hnfvtxtrkbbc"); hnfvtxtrkbbc = new TH2F(name,name,500,0,500,600,0,300);
  sprintf(name,"hnfvtxtrksnmpcs"); hnfvtxtrksnmpcs = new TH2F(name,name,200,0,200,500,0,500);
  sprintf(name,"hnfvtxtrksnmpcn"); hnfvtxtrknnmpcn = new TH2F(name,name,200,0,200,500,0,500);
  sprintf(name,"south_mpc_south_bbc");  south_mpc_south_bbc = new TH2F(name, name, 600, -0.5, 19.5, 600, -0.5, 199.5);
  sprintf(name,"north_mpc_north_bbc");  north_mpc_north_bbc = new TH2F(name, name, 600, -0.5, 19.5, 600, -0.5, 199.5);

  //run QA
  sprintf(name,"hrunbbcs"); hrunbbcs = new TH2F(name,name,432008-422070+21,422070-10.5,432008+10.5,600,0,300);
  sprintf(name,"hrunbbcn"); hrunbbcn = new TH2F(name,name,432008-422070+21,422070-10.5,432008+10.5,600,0,300);
  sprintf(name,"hrunntrack_arm0_pos"); hrunntrack[1] = new TH2F(name,name,432008-422070+21,422070-10.5,432008+10.5,50,0,50);
  sprintf(name,"hrunntrack_arm0_neg"); hrunntrack[2] = new TH2F(name,name,432008-422070+21,422070-10.5,432008+10.5,50,0,50);
  sprintf(name,"hrunntrack_arm1_pos"); hrunntrack[3] = new TH2F(name,name,432008-422070+21,422070-10.5,432008+10.5,50,0,50);
  sprintf(name,"hrunntrack_arm1_neg"); hrunntrack[4] = new TH2F(name,name,432008-422070+21,422070-10.5,432008+10.5,50,0,50);
  sprintf(name,"hrunntrack"); hrunntrack[0] = new TH2F(name,name,432008-422070+21,422070-10.5,432008+10.5,50,0,50);

  return 0;
}

int Perform::process_event()
{
  int nEvent = tree->GetEntries();
  cout << nEvent << endl;
  for(int ievent=0;ievent < nEvent; ievent++){
      tree->GetEntry(ievent);

  if(ievent%100000==0) {
      std::cout<<"************* ievent= "<<ievent<<"    *************"<<std::endl;
  }

//global
  //std::cout << "RunNumber = "<<run<<std::endl;
  bool isMB = (trig & 0x00000010);
  bool isCentral = (trig & 0x00000008);
  bool isFVtxsouth = 0;
  bool isFVtxnorth = 0;
  bool isFVtxAnd = 0;
  bool isFVtxOr = 0;
  if(coll.Contains("pA")){
      isFVtxsouth = trig & 0x00000800;
      isFVtxnorth = trig & 0x00000400;
      isFVtxAnd = (trig & 0x00000400) && (trig & 0x00000800);
      isFVtxOr = (trig & 0x00000400) || (trig & 0x00000800);
  }
  else if(coll.Contains("pp")){
      if(run>425926){
          isFVtxAnd = (trig & 0x00000400);
          isFVtxOr = (trig & 0x00000800);
      }
      else{
          isFVtxAnd = (trig & 0x00000400) && (trig & 0x00000800);
          isFVtxOr = (trig & 0x00000400) || (trig & 0x00000800);
      }
  }
  else return 0;
  if(coll.Contains("minbias")){
      if(!isMB) continue;   //MinBias trig
  }
  else if(coll.Contains("central")){
      if(!isCentral) continue;   //Central trig
  }
  else if(coll.Contains("fvtxor")){
      if(! isFVtxOr) continue;   //FVtx trig or 
  }
  else if(coll.Contains("fvtxand")){
      if(! isFVtxAnd) continue;   //FVtx trig or 
  }
  else if(coll.Contains("fvtxsouth")){
      if(! isFVtxsouth) continue;   //FVtx trig south
  }
  else if(coll.Contains("fvtxnorth")){
      if(! isFVtxnorth) continue;   //FVtx trig south
  }
  else return 0;

  int icent = -9999;
  if(bbc_s>0) icent=0;
  if(icent>=centbin || icent <0) continue;
//cnt track
  int ntrk = 0;
  int ntof = 0;
  int ntrkdiff[5]={};
  for(int itrk=0;itrk<ntrack;itrk++){
      float cntphi       = phi0[itrk];
      float pt        = mom[itrk] * sin(the0[itrk]);
      float cnteta       = -log(tan(0.5*the0[itrk]));

      hcntetaphi->Fill(cnteta,cntphi);      
      hcntpt->Fill(pt);
      for(int ipt = 0;ipt<50;ipt++){
          if(mom[itrk]>=0.1*ipt && mom[itrk]<0.1*(ipt+1)){
              if(charge[itrk] > 0 && arm[itrk] == 0)	pc3dphidz_arm0_pos[ipt]->Fill(pc3dphi[itrk],pc3dz[itrk]);
              if(charge[itrk] > 0 && arm[itrk] == 1)	pc3dphidz_arm1_pos[ipt]->Fill(pc3dphi[itrk],pc3dz[itrk]);
              if(charge[itrk] < 0 && arm[itrk] == 0)	pc3dphidz_arm0_neg[ipt]->Fill(pc3dphi[itrk],pc3dz[itrk]);
              if(charge[itrk] < 0 && arm[itrk] == 1)	pc3dphidz_arm1_neg[ipt]->Fill(pc3dphi[itrk],pc3dz[itrk]);
	}
      }
 //     float emcsdphi = d_scnt->get_emcsdphi();
 //     float emcsdz = d_scnt->get_emcsdz();
 //     bool emc_matching = fabs(emcsdphi)<2.0 && fabs(emcsdz)<2.0;
      
      /*
       if(pt>0.2&&pt<5.0)&&emc_matching)
      if(arm[itrk]==1){
	if(striptofw>=0 && tofw_matching && qtofw>60 && qtofw<600){
	  good_tofw = true;
	}
      }
      bool good_tofe=false;
      if(arm[itrk]==0){
	if(slat>=0 && tofe_matching && qtofe>0.002){
	  good_tofe=true;
	}
      }
      */

      if(slat[itrk]>0&&pt>0.02&&pt<5.0&&tofpl[itrk]>400){
	if(arm[itrk]==0 && charge[itrk] > 0) 	tofdphidz[0]->Fill(tofdphi[itrk],tofdz[itrk]);
	if(arm[itrk]==0 && charge[itrk] < 0)	tofdphidz[1]->Fill(tofdphi[itrk],tofdz[itrk]);
      	if(arm[itrk]==1 && charge[itrk] > 0) 	tofwdphidz[0]->Fill(tofdphi[itrk],tofdz[itrk]);
	if(arm[itrk]==1 && charge[itrk] < 0)	tofwdphidz[1]->Fill(tofdphi[itrk],tofdz[itrk]);
      float deltat = ttof[itrk] - tofpl[itrk]/phbeta*sqrt(mpion*mpion/mom[itrk]/mom[itrk]+1);
      float m2tof = mom[itrk]*mom[itrk]*(ttof[itrk]*ttof[itrk]*phbeta*phbeta/(tofpl[itrk]*tofpl[itrk])-1);
      //if(m2tof>0.6 && m2tof<1.2)//proton
      //if(m2tof>-0.2 && m2tof<0.1)//pion
      ttofqpratio->Fill(ttof[itrk],charge[itrk]/mom[itrk]);
      m2qpratio->Fill(m2tof,charge[itrk]/mom[itrk]);
      ttofp->Fill(ttof[itrk],mom[itrk]);
      m2p->Fill(m2tof,mom[itrk]);
      pinv2chbeta->Fill(1./mom[itrk],charge[itrk]*(1.-mom[itrk]/TMath::Sqrt(mom[itrk]*mom[itrk]+m2tof)));
       if(mom[itrk]>0.5&&mom[itrk]<2.0){
       if(arm[itrk] == 0) deltattofeis->Fill(slat[itrk],deltat);
       if(arm[itrk] == 1) deltattofwis->Fill(slat[itrk],deltat);
       }
      ntof++;
      }
      if(pt>0.2&&pt<5.0){
    if(charge[itrk] > 0 && arm[itrk] == 0)     ntrkdiff[1]++; 
    if(charge[itrk] < 0 && arm[itrk] == 0)     ntrkdiff[2]++;
    if(charge[itrk] > 0 && arm[itrk] == 1)     ntrkdiff[3]++;
    if(charge[itrk] < 0 && arm[itrk] == 1)     ntrkdiff[4]++;
	ntrk++;
      }
  }
  /*
  // beam beam counter (bbc) r.p. // -----------------------------------
  int nbbcau=0;//bbc south
  float bbcau_et[nbbcs];
  float bbcau_phi[nbbcs];

  int nbbcde=0;//bbc north
  float bbcde_et[nbbcs];
  float bbcde_phi[nbbcs];

  for (int ipmt=0; ipmt<nbbc; ipmt++) {
    if (bbccharge[ipmt]>0) {
      int iarm = 0;
      if (bbcz[ipmt] > 0) iarm = 1;
      float bbcphi=atan2(bbcy[ipmt],bbcx[ipmt]);
      float val=bbccharge[ipmt];
      
      //float rad = sqrt(bbcx*bbcx+bbcy*bbcy);
      //float the = atan2(rad,bbcz - bbcv); // bbcz-bbcv*10.0);
      //float eta = -log(tan(0.5*the));

      if (val>0) {
	if(iarm==0){
	  bbcau_et[nbbcau] = val;
	  bbcau_phi[nbbcau] = bbcphi;
	  nbbcau++;
	}
	else{
	  bbcde_et[nbbcde] = val;
          bbcde_phi[nbbcde] = bbcphi;
          nbbcde++;
	}
      bbcet->Fill(val);
      }
    }
  }
  */
/*
//Clusters
  int nclu[4]={};

  int ncluau=0;//clu south
  int nclude=0;//clu north

  if(coll.Contains("pp")){
  float cluau_eta[ncluss];
  float cluau_phi[ncluss];

  float clude_eta[ncluss];
  float clude_phi[ncluss];

  for (int ij=0; ij<nvtx; ij++) {
    float X = vtxX[ij];
    float Y = vtxY[ij];
    float z = vtxZ[ij];
    int   ilayer = layer[ij];

    if (ilayer>=0 && ilayer<4) {
      float cluphi=atan2(Y,X);
      float rad=sqrt(X*X+Y*Y);
      float the=atan2(rad,z-vtxz);
      float clueta=-log(tan(0.5*the));

      if(X!=X||Y!=Y||z!=z||fabs(clueta)>=3.5) continue;

	nclu[ilayer]++;
	if(ilayer==0 && nclu[0]<=80)
        hcluetaphi[ilayer]->Fill(clueta,cluphi);
	if(ilayer!=0)
        hcluetaphi[ilayer]->Fill(clueta,cluphi);

      if (ilayer==0) {
	if(fabs(clueta)>1.0 && fabs(clueta)<3.0){
	  cluau_phi[ncluau] = cluphi;
	  cluau_eta[ncluau] = clueta;
	  ncluau++;
	}
	if(fabs(clueta)<1.0){
	  clude_phi[nclude] = cluphi;
	  clude_eta[nclude] = clueta;
	  nclude++;
	}
      }
    }  
  }
  }

//  if(nvtx>200) return 1;
//  if(ncluau<=4) return 1;
*/
  /*
//FVTX

  int nfvtxau=0;//fvtxau in the et cut
  float fvtxau_phi[nfvtxs];
  float fvtxau_eta[nfvtxs];

  int nfvtxde=0;//fvtxau in the et cut
  float fvtxde_phi[nfvtxs];
  float fvtxde_eta[nfvtxs];
  
  int nfvtxtrk=0;
  int nfvtxtrks=0;
  int nfvtxtrkn=0;

  int nfvtxtrkall=0;
  int nfvtxtrkalls=0;
  int nfvtxtrkalln=0;

  for (int ij=0; ij<nfvtxtrack; ij++) {
      short fvtx_arm = farm[ij];
      short nhits = fnhits[ij];
      float fvtx_the = fthe[ij];
      float fvtx_eta = feta[ij];
      float fvtx_phi = fphi[ij];

      float 	fvtx_x = fvtxX[ij];
      float 	fvtx_y = fvtxY[ij];
      float 	fvtx_z = fvtxZ[ij];
 	
      if(fvtx_the==0) continue;
 
      nfvtxtrkall++;
      if(fvtx_arm==0) nfvtxtrkalls++;
      if(fvtx_arm==1) nfvtxtrkalln++;
      //short_chi2 = fvtx_trk->get_short_chi2_ndf();
      //chi2 = fvtx_trk->get_chi2_ndf();
      float DCA_x = fvtx_x + tan(fvtx_the)*cos(fvtx_phi)*(fvtxz - fvtx_z);
      float DCA_y = fvtx_y + tan(fvtx_the)*sin(fvtx_phi)*(fvtxz - fvtx_z);
      float DCA_R = sqrt((DCA_x*DCA_x) + (DCA_y*DCA_y));

      if(nhits>=3){
	DCAxydis[fvtx_arm]->Fill(DCA_x, DCA_y);
	DCAcentdis[fvtx_arm]->Fill(icent, DCA_R);
      }
    
      float sigma_dcay = 999.9;
      float sigma_dcax = 999.9;

      if(fvtx_arm==0){
	sigma_dcax= fabs(DCA_x+0.1425)/0.1585;
	sigma_dcay= fabs(DCA_y-0.0527)/0.1668;
      }
      else{
	sigma_dcax= fabs(DCA_x+0.1591)/0.1643;
	sigma_dcay= fabs(DCA_y-0.0533)/0.1990;
      }

     bool dcacut = true;//sigma_dcax<2.0 && sigma_dcay<2.0;
      if(fvtx_phi<10 && fvtx_phi > -10 && fabs(fvtx_eta)<3.5 && dcacut && nhits>=3){
      //if(fvtx_phi<10 && fvtx_phi > -10 && fabs(fvtx_eta)<3.5 && DCA_R<2.0)
	if(icent==0){
	  hvtx0etaz->Fill(fvtxz, fvtx_eta);
	  hvtx0etaphi->Fill(fvtx_phi, fvtx_eta);
	} 
	if(icent==0 && fabs(fvtx_eta)<2.5 && fabs(fvtx_eta)>1.5){
          hvtx1etaz->Fill(fvtxz, fvtx_eta);
          hvtx1etaphi->Fill(fvtx_phi, fvtx_eta);
        }

	DCAxy2dis[fvtx_arm]->Fill(DCA_x, DCA_y);

	//nfvtxtrk++;
	//if(arm==0) nfvtxtrks++;
	//else nfvtxtrkn++;
	//old 1.2 2.7

	//bool eta_south = (arm==0)&&(fvtx_eta<-0.034*fvtxz-1.5) && (fvtx_eta>-0.034*fvtxz-2.5);
	//bool eta_north = (arm==1)&&(fvtx_eta<-0.034*fvtxz+2.5) && (fvtx_eta>-0.034*fvtxz+1.5);
	
	//if(!eta_south && !eta_north) continue;
	
	nfvtxtrk++;
        if(fvtx_arm==0) nfvtxtrks++;
        else nfvtxtrkn++;
	//cout<<"nfvtxtrk "<<nfvtxtrk<<" nfvtxtrks "<<nfvtxtrks<<" nfvtxtrkn "<<nfvtxtrkn<<endl;

	if(nfvtxtrk<150&&nfvtxtrks<100&&nfvtxtrkn<60){
	  if(fvtx_eta<-1.0){
	    fvtxau_phi[nfvtxau] = fvtx_phi;
	    fvtxau_eta[nfvtxau] = fvtx_eta;
	    nfvtxau++;
	  }
	  else if(fvtx_eta>1.0){
	    fvtxde_phi[nfvtxde] = fvtx_phi;
	    fvtxde_eta[nfvtxde] = fvtx_eta;
	    nfvtxde++;
	  }
	}
      }
    }
*/
/*
 // muon piston calorimeter (mpc) r.p. // -----------------------------------
  float mpc_southe=0, mpc_northe=0;
 
  int nmpc_south = 0;
  int nmpc_north = 0;
 
  int nmpcau=0;//mpcau in the et cut
  float mpcau_et[nmpcs];
  float mpcau_phi[nmpcs];
  float mpcau_eta[nmpcs];

  int nmpcde=0;//mpcau in the et cut
  float mpcde_et[nmpcs];
  float mpcde_phi[nmpcs];
  float mpcde_eta[nmpcs];

  int itower;
  for (itower=0; itower<nmpc; itower++){
    float mpce = mpc_e[itower];
    float mpcx = mpc_x[itower];       // get x position of tower
    float mpcy = mpc_y[itower];       // get y position of tower
    float mpcz = mpc_z[itower];       // get z position of tower
    float mpc_phi = atan2(mpcy,mpcx); // + 0.015*dumflt->GetRandom();
    float mpc_rad = sqrt(mpcx*mpcx+mpcy*mpcy);
    float mpc_the = atan2(mpc_rad,mpcz-bbcv); // mpc_z-bbcv);

    float mpc_eta = -log(tan(0.5*mpc_the));

    if (fabs(mpcz) < 1.) continue;          // skip non-existent channels

    int iarm = 0; //south
    if (mpcz > 0) iarm = 1;//north
    int iring = 0;
    mpc_rad=fabs(mpc_eta);
    iring = (int)((mpc_rad-3.0)/0.25);
    if (iring<0) iring=0;
    if (iring>3) iring=3;
    float val = mpce*sin(mpc_the);

    //mpcgain->Fill(fee_ch+0.0, val);
    //mpcmapp->Fill(fee_ch+0.0, val);

    if(mpce>mpc_e_cut){
      //if(iarm==0) mpcsouthgain->Fill(fee_ch+0.0, val);
      //else if(iarm==1) mpcnorthgain->Fill(fee_ch+0.0, val);

      if(iarm==0) {
	nmpc_south++;
	mpc_southe += val;
      }
      else if(iarm==1) {
	nmpc_north++;
	mpc_northe += val;
      }
    }

  //reaction plane
    for (int ih=0; ih<4; ih++) {
      float vc=val*cos((ih+1.0)*mpc_phi);
      float vs=val*sin((ih+1.0)*mpc_phi);
      int ioff=iring*4;
      if (iarm==0) {sumxy[ih][12][0]+=vc; sumxy[ih][12][1]+=vs; sumxy[ih][12][2]+=val;}
      if (iarm==1) {sumxy[ih][13][0]+=vc; sumxy[ih][13][1]+=vs; sumxy[ih][13][2]+=val;}
      if (iarm<2 ) {sumxy[ih][14][0]+=vc; sumxy[ih][14][1]+=vs; sumxy[ih][14][2]+=val;}
      if (iarm==0) {sumxy[ih][22+ioff][0]+=vc; sumxy[ih][22+ioff][1]+=vs; sumxy[ih][22+ioff][2]+=val;}
      if (iarm==1) {sumxy[ih][23+ioff][0]+=vc; sumxy[ih][23+ioff][1]+=vs; sumxy[ih][23+ioff][2]+=val;}
      if (iarm<2 ) {sumxy[ih][24+ioff][0]+=vc; sumxy[ih][24+ioff][1]+=vs; sumxy[ih][24+ioff][2]+=val;}
      if (ih%2==0 && iarm==1)  {vc*=-1.0; vs*=-1.0;}
      if (iarm<2 ) {sumxy[ih][15][0]+=vc; sumxy[ih][15][1]+=vs; sumxy[ih][15][2]+=val;}
      if (iarm<2 ) {sumxy[ih][25+ioff][0]+=vc; sumxy[ih][25+ioff][1]+=vs; sumxy[ih][25+ioff][2]+=val;}
    }
    
    if(mpce>mpc_e_cut){
      if(iarm==0) {
	mpcetetasouth[icent]->Fill(mpc_eta, val);
      }
      else if(iarm==1) {
	mpcetetanorth[icent]->Fill(mpc_eta, val);
      }
    }

    //if(val>0.5 && val<1.0)
    if(mpce>mpc_e_cut){
      if(iarm==0){
	mpcau_et[nmpcau] = val;
	mpcau_phi[nmpcau] = mpc_phi;
	mpcau_eta[nmpcau] = mpc_eta;
	nmpcau++;
      }
      else {
	mpcde_et[nmpcde] = val;
        mpcde_phi[nmpcde] = mpc_phi;
	mpcde_eta[nmpcde] = mpc_eta;
        nmpcde++;
      }
    }  
  }
*/
 
//  mpc_south_north->Fill(nmpc_south, nmpc_north);
//  south_mpc_north_mpc->Fill(mpc_southe, mpc_northe);
  
//Correlate histograms
//  hbbcsZdcEs->Fill(bbc_s,ZdcEs);
//  if(ievent % 2 == 1) hmixbbcsZdcEs->Fill(bbc_s+buff_bbcs,ZdcEs+buff_ZdcEs);
//  buff_bbcs = bbc_s;
//  buff_ZdcEs = ZdcEs;
/*
  hvtxzfvtxz->Fill(vtxz,fvtxz);
  hpc1hitsbbc->Fill(npc1hits,bbc_s+bbc_n);
  hnpc3hitsntof->Fill(ntrack,ntof);
  hbbcnbbc->Fill(bbc_s+bbc_n,nbbcau+nbbcde);
  */
  hbbcsbbcn->Fill(bbc_s,bbc_n);
/*
  for (int ilayer=0; ilayer<4; ilayer++) {
	hnvtxnfvtxtrk[ilayer]->Fill(nclu[ilayer],nfvtxtrk);
	hnvtxnmpc[ilayer]->Fill(nclu[ilayer],nmpc_south+nmpc_north);
	hnvtxntrk[ilayer]->Fill(nclu[ilayer],ntrk);
  	hbbcsnvtx[ilayer]->Fill(bbc_s,nclu[ilayer]);
  	hbbcnnvtx[ilayer]->Fill(bbc_n,nclu[ilayer]);
  	hbbcnvtx[ilayer]->Fill(bbc_s+bbc_n,nclu[ilayer]);
}
  hnbbcnclu->Fill(nbbcau+nbbcde,ncluau+nclude);
  hntracknmpc->Fill(ntrk,nmpc_south+nmpc_north);
  hnfvtxtrkbbc->Fill(nfvtxtrk,bbc_s+bbc_n);
  hnfvtxtrksnmpcs->Fill(nfvtxtrks,nmpc_south);
  hnfvtxtrknnmpcn->Fill(nfvtxtrkn,nmpc_north);
  mpc_south_cent->Fill(icent, nmpc_south);
  mpc_north_cent->Fill(icent, nmpc_north);
  south_mpc_south_bbc->Fill(mpc_southe, bbc_s);
  north_mpc_north_bbc->Fill(mpc_northe, bbc_n);
*/
  hrunbbcs->Fill(run,bbc_s);
  hrunbbcn->Fill(run,bbc_n);
  ntrkdiff[0] = ntrk;
  for(int i=0;i<4;i++){
      hrunntrack[i]->Fill(run,ntrkdiff[i]);
 }
}
return 0;
}

int Perform::End()
{

    std::cout << "InputFileName = " << InputFileName << std::endl;
    std::cout << "OutputFileName = " << OutputFileName << std::endl;


  if(d_outfile) {
    d_outfile->cd();

//Tracks
//  hcntetaphi->Write();
//  hcntpt->Write();
for(int i=0;i<50;i++){
pc3dphidz_arm0_pos[i]->Write();
pc3dphidz_arm1_pos[i]->Write();
pc3dphidz_arm0_neg[i]->Write();
pc3dphidz_arm1_neg[i]->Write();
}
/*
//tof
for(int i=0;i<2;i++){
  tofdphidz[i]->Write();
  tofwdphidz[i]->Write();
}
  ttofqpratio->Write();
  m2qpratio->Write();
  m2p->Write();
  ttofp->Write();
  pinv2chbeta->Write();
  deltattofeis->Write();
  deltattofwis->Write();
*/
  /*  
//vtx
for(int i=0;i<4;i++){
  hcluetaphi[i]->Write();
}

//bbc
  bbcet->Write();

//fvtx
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
//mpc
 /*
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
//  hvtxzfvtxz->Write();
//  hbbcsZdcEs->Write();
//  hmixbbcsZdcEs->Write();
//  hpc1hitsbbc->Write();
//  hnpc3hitsntof->Write();
//  hbbcnbbc->Write();
  hbbcsbbcn->Write();
/*
  for(int i=0;i<4;i++){
  hnvtxnfvtxtrk[i]->Write();
  hnvtxnmpc[i]->Write();
  hnvtxntrk[i]->Write();
  hbbcsnvtx[i]->Write();
  hbbcnnvtx[i]->Write();
  hbbcnvtx[i]->Write();
}
*/
//  hnbbcnclu->Write();
//  hntracknmpc->Write();
//  hnfvtxtrkbbc->Write();
//  hnfvtxtrksnmpcs->Write();
//  hnfvtxtrknnmpcn->Write();
//  south_mpc_south_bbc->Write();
//  north_mpc_north_bbc->Write();
  hrunbbcs->Write();
  hrunbbcn->Write();
  for(int i=0;i<5;i++){
    hrunntrack[i]->Write();
  }
  }
  else{

      std::cout<<"ERROR: No output file set!"<<std::endl;

  }

  return 0;
}
