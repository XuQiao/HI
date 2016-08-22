#include "Ridgecntbbc.h"
#include <stdlib.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include "TString.h"
#include <TRandom.h>
#include <fstream>
#include <iostream>

#include "Run15pppc3dphidzcalibsmoothpass1.h"

using namespace std;

const int bbcz_cut = 10;
const int vtxz_cut = 10;

const int ntrks = 200;
const int nbbcs = 800;
const int ncluss = 800;

const int centbin = 8;
const int ptbin = 10;
const int vzbin = 10;
const int nbuff = 40;

const int nbufftrks = centbin*vzbin*nbuff*ntrks;

const int nbuffbbcs = centbin*vzbin*nbuff*nbbcs;

typedef vector<float> track_buff;
track_buff trackbuff_pt(nbufftrks);
track_buff trackbuff_phi(nbufftrks);
track_buff trackbuff_eta(nbufftrks);

typedef vector<float> bbcau_buff;
bbcau_buff bbcaubuff_et(nbuffbbcs);
bbcau_buff bbcaubuff_phi(nbuffbbcs);
bbcau_buff bbcaubuff_eta(nbuffbbcs);

typedef vector<float> bbcde_buff;
bbcde_buff bbcdebuff_et(nbuffbbcs);
bbcde_buff bbcdebuff_phi(nbuffbbcs);
bbcde_buff bbcdebuff_eta(nbuffbbcs);

int ntrack_buff[centbin][vzbin];
int dtrack_buff[centbin][vzbin];
int buff_ntrack[centbin][vzbin][nbuff];

int nbbcau_buff[centbin][vzbin];
int dbbcau_buff[centbin][vzbin];
int buff_nbbcau[centbin][vzbin][nbuff];

int nbbcde_buff[centbin][vzbin];
int dbbcde_buff[centbin][vzbin];
int buff_nbbcde[centbin][vzbin][nbuff];

Ridgecntbbc::Ridgecntbbc(TString input, TString output):
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

  //QVector
  Qx(),
  Qy(),
  Qw(),

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

Ridgecntbbc::~Ridgecntbbc()
{
    delete tree;
    d_infile->Close();
    d_outfile->Close();
}

int Ridgecntbbc::Init()
{
  coll = "";
  TString temp = OutputFileName;
  temp.ToLower();
  //fetchdphidz();
  if(temp.Contains("pp") && temp.Contains("mb")) coll = "ppminbias";
  else if(temp.Contains("pal") && temp.Contains("mbst")) coll = "pAlminbias";
  else if(temp.Contains("pal") && temp.Contains("mbcentral")) coll = "pAlcentral";
  else if(temp.Contains("pau") && temp.Contains("mbst")) coll = "pAuminbias";
  else if(temp.Contains("pau") && temp.Contains("mbcentral")) coll = "pAucentral";
  else if(temp.Contains("pp") && temp.Contains("fvtxand")) coll = "ppfvtxand";
  else if(temp.Contains("pal") && temp.Contains("fvtxand")) coll = "pAlfvtxand";
  else if(temp.Contains("pau") && temp.Contains("fvtxand")) coll = "pAufvtxand";
  else if(temp.Contains("pp") && temp.Contains("fvtxor")) coll = "ppfvtxor";
  else if(temp.Contains("pal") && temp.Contains("fvtxor")) coll = "pAlfvtxor";
  else if(temp.Contains("pal") && temp.Contains("fvtxsouth")) coll = "pAlfvtxsouth";
  else if(temp.Contains("pal") && temp.Contains("fvtxnorth")) coll = "pAlfvtxnorth";
  else if(temp.Contains("pau") && temp.Contains("fvtxor")) coll = "pAufvtxor";
  else if(temp.Contains("pau") && temp.Contains("fvtxsouth")) coll = "pAufvtxsouth";
  else if(temp.Contains("pau") && temp.Contains("fvtxnorth")) coll = "pAufvtxnorth";
  else return 0;

//  OutputFileName.Insert(OutputFileName.Length()-5,coll);

  d_outfile = new TFile(OutputFileName,"recreate");
      d_infile = TFile::Open(InputFileName,"ReadOnly");
  if(!d_infile) return 0;
  tree = (TTree*)d_infile->Get("tree");
  if(!tree) return 0;
  tree -> SetBranchAddress("run", &run);
  tree -> SetBranchAddress("trig", &trig);
  if(coll.Contains("pA"))
  tree -> SetBranchAddress("cent", &cent);
  tree -> SetBranchAddress("npc1hits", &npc1hits);
  tree -> SetBranchAddress("bbc_s", &bbc_s);
  tree -> SetBranchAddress("bbc_n", &bbc_n);
  tree -> SetBranchAddress("bbcv", &bbcv);
  if(!coll.Contains("pAl"))
  tree -> SetBranchAddress("vtxz", &vtxz);
//  tree -> SetBranchAddress("fvtxz", &fvtxz);
  tree -> SetBranchAddress("ntrack", &ntrack);
  tree -> SetBranchAddress("nbbc", &nbbc);
  if(!coll.Contains("pAl")){
  tree -> SetBranchAddress("nvtx", &nvtx);
  tree -> SetBranchAddress("nvtxtrack", &nvtxtrack);
  }
//  tree -> SetBranchAddress("nmpc", &nmpc);
// tree -> SetBranchAddress("nfvtxtrack", &nfvtxtrack);

  tree -> SetBranchAddress("mom", &mom);
  tree -> SetBranchAddress("phi0", &phi0);
  tree -> SetBranchAddress("the0", &the0);
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
  tree -> SetBranchAddress("bbcx", &bbcx);
  tree -> SetBranchAddress("bbcy", &bbcy);
  tree -> SetBranchAddress("bbcz", &bbcz);
/*
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
  */
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
 pi = acos(-1);


  for (int icent=0; icent<centbin; icent++) {
    //bbc correlation
    for(int ipt=0; ipt<ptbin; ipt++){
      //cnt with bbc
      sprintf(name,"hforesouthbbc_%d_%d",icent,ipt);
      hforesouthbbc[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacksouthbbc_%d_%d",icent,ipt);
      hbacksouthbbc[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hforenorthbbc_%d_%d",icent,ipt);
      hforenorthbbc[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacknorthbbc_%d_%d",icent,ipt);
      hbacknorthbbc[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      //flip
      sprintf(name,"kforesouthbbc_%d_%d",icent,ipt);
      kforesouthbbc[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksouthbbc_%d_%d",icent,ipt);
      kbacksouthbbc[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kforenorthbbc_%d_%d",icent,ipt);
      kforenorthbbc[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthbbc_%d_%d",icent,ipt);
      kbacknorthbbc[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);
//with weight
      sprintf(name,"hforesouthbbcw_%d_%d",icent,ipt);
      hforesouthbbcw[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacksouthbbcw_%d_%d",icent,ipt);
      hbacksouthbbcw[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hforenorthbbcw_%d_%d",icent,ipt);
      hforenorthbbcw[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacknorthbbcw_%d_%d",icent,ipt);
      hbacknorthbbcw[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      //flip
      sprintf(name,"kforesouthbbcw_%d_%d",icent,ipt);
      kforesouthbbcw[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacksouthbbcw_%d_%d",icent,ipt);
      kbacksouthbbcw[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kforenorthbbcw_%d_%d",icent,ipt);
      kforenorthbbcw[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthbbcw_%d_%d",icent,ipt);
      kbacknorthbbcw[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);



      //mixing with bbc
      sprintf(name,"hbacksouthbbc2_%d_%d",icent,ipt);
      hbacksouthbbc2[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacknorthbbc2_%d_%d",icent,ipt);
      hbacknorthbbc2[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      //flip
      sprintf(name,"kbacksouthbbc2_%d_%d",icent,ipt);
      kbacksouthbbc2[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthbbc2_%d_%d",icent,ipt);
      kbacknorthbbc2[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      //with weight
      sprintf(name,"hbacksouthbbcw2_%d_%d",icent,ipt);
      hbacksouthbbcw2[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"hbacknorthbbcw2_%d_%d",icent,ipt);
      hbacknorthbbcw2[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      //flip
      sprintf(name,"kbacksouthbbcw2_%d_%d",icent,ipt);
      kbacksouthbbcw2[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);

      sprintf(name,"kbacknorthbbcw2_%d_%d",icent,ipt);
      kbacknorthbbcw2[icent][ipt] = new TH1F(name,name,40, -0.5*pi, 1.5*pi);


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
    }
  }

  memset(ntrack_buff, 0, sizeof(ntrack_buff));
  memset(dtrack_buff, 0, sizeof(dtrack_buff));
  memset(buff_ntrack, 0, sizeof(buff_ntrack));

  memset(nbbcau_buff, 0, sizeof(nbbcau_buff));
  memset(dbbcau_buff, 0, sizeof(dbbcau_buff));
  memset(buff_nbbcau, 0, sizeof(buff_nbbcau));

  memset(nbbcde_buff, 0, sizeof(nbbcde_buff));
  memset(dbbcde_buff, 0, sizeof(dbbcde_buff));
  memset(buff_nbbcde, 0, sizeof(buff_nbbcde));

  return 0;
}

int Ridgecntbbc::process_event()
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
  if(fabs(bbcv)>bbcz_cut) continue;
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

  float bbc_add = bbc_s + bbc_n;
  
  int ivz = 0;
  int icent=  -9999;
  if(bbc_add >= 26.3) icent = 0;
  else if(bbc_add >= 19.3) icent = 1;
  else if(bbc_add >= 16.3) icent = 2;
  else if(bbc_add >= 13) icent = 3;
  else if(bbc_add >= 11) icent = 4;
  else if(bbc_add >= 9.6) icent = 5;
  else if(bbc_add >= 7) icent = 6;
  else icent = 7;
  if(icent>=centbin||icent<0) continue;
 
//cnt track
  int ntrk = 0;
  float track_pt[ntrks];
  float track_phi[ntrks];
  float track_eta[ntrks];
  for(int itrk=0;itrk<ntrack;itrk++){
      float cntphi       = phi0[itrk];
      float pt        = mom[itrk] * sin(the0[itrk]);
      float cnteta       = -log(tan(0.5*the0[itrk]));
      double sdphi = calcsdphi(pc3dphi[itrk],arm[itrk],charge[itrk],mom[itrk]);
      double sdz = calcsdz(pc3dz[itrk],arm[itrk],charge[itrk],mom[itrk]);
      if(fabs(sdphi)<3.0 && fabs(sdz)<3.0){
          if(pt>0.2&&pt<5.0){
              track_pt[ntrk]=pt;
              track_phi[ntrk]=cntphi;
              track_eta[ntrk]=cnteta;
              ntrk++;
          }
      }
  }

  // beam beam counter (bbc) r.p. // -----------------------------------
  int nbbcau=0;//bbc south
  float bbcau_et[nbbcs];
  float bbcau_phi[nbbcs];
  float bbcau_eta[nbbcs];

  int nbbcde=0;//bbc north
  float bbcde_et[nbbcs];
  float bbcde_phi[nbbcs];
  float bbcde_eta[nbbcs];

  for (int ipmt=0; ipmt<nbbc; ipmt++) {
    if (bbccharge[ipmt]>0) {
      int iarm = 0;
      if (bbcz[ipmt] > 0) iarm = 1;
      float bbcphi=atan2(bbcy[ipmt],bbcx[ipmt]);
      float val=bbccharge[ipmt];
      
      float rad = sqrt(bbcx[ipmt]*bbcx[ipmt]+bbcy[ipmt]*bbcy[ipmt]);
      float the = atan2(rad,bbcz[ipmt] - bbcv); // bbcz-bbcv*10.0);
      float eta = -log(tan(0.5*the));

      if (val>0) {
	if(iarm==0){
	  bbcau_et[nbbcau] = val;
	  bbcau_phi[nbbcau] = bbcphi; 
	  bbcau_eta[nbbcau] = eta; 
	  nbbcau++;
	}
	else{
	  bbcde_et[nbbcde] = val;
          bbcde_phi[nbbcde] = bbcphi;
          bbcde_eta[nbbcde] = eta;
          nbbcde++;
	}
      }
    }
  }


  //correlation
  int sign =(int)(gRandom->Rndm()*2);//put the flip sign here since per track/evt
  //1) //foreground
  
  for(int itrk=0; itrk<ntrk; itrk++){
    int ipt_cor = track_pt[itrk]/0.5;
    //BBC
    for(int ipmt=0; ipmt<nbbcau; ipmt++){
      float val = bbcau_et[ipmt];
      float bbc_phi = bbcau_phi[ipmt];
      float bbc_eta = bbcau_eta[ipmt];

      float dphi = track_phi[itrk] - bbc_phi;
      float deta = track_eta[itrk] - bbc_eta;
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
      kforesouthetabbc[icent][ipt_cor]->Fill(deta, dphia);

      hforesouthbbcw[icent][ipt_cor]->Fill(dphi,val*val);
      kforesouthbbcw[icent][ipt_cor]->Fill(dphia,val);
      kforesouthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);
    }//bbc south

    for(int ipmt=0; ipmt<nbbcde; ipmt++){
      float val = bbcde_et[ipmt];
      float bbc_phi = bbcde_phi[ipmt];
      float bbc_eta = bbcde_eta[ipmt];

      float dphi = track_phi[itrk] - bbc_phi;
      float deta = track_eta[itrk] - bbc_eta;
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
      kforenorthetabbc[icent][ipt_cor]->Fill(deta, dphia);

      hforenorthbbcw[icent][ipt_cor]->Fill(dphi,val*val);
      kforenorthbbcw[icent][ipt_cor]->Fill(dphia,val);
      kforenorthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);
    }//bbc north
  }
  
  //2) associate tower + mixing track
  if(dtrack_buff[icent][ivz]>0){
    
    for(int ibuff=0;ibuff<dtrack_buff[icent][ivz];ibuff++){
	
      for(int itrk_buff=0; itrk_buff<buff_ntrack[icent][ivz][ibuff]; itrk_buff++){
	
	int ipt_cor = trackbuff_pt[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff]/0.5;
	
	for(int ipmt=0; ipmt<nbbcau; ipmt++){
          float val = bbcau_et[ipmt];
          float bbc_phi = bbcau_phi[ipmt];
          float bbc_eta = bbcau_eta[ipmt];

          float dphi = trackbuff_phi[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff] - bbc_phi;
	  float deta = trackbuff_eta[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff] - bbc_eta;
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
          kbacksouthetabbc[icent][ipt_cor]->Fill(deta, dphia);

          hbacksouthbbcw[icent][ipt_cor]->Fill(dphi,val*val);
          kbacksouthbbcw[icent][ipt_cor]->Fill(dphia,val);
          kbacksouthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);

        }//bbcsouth

	for(int ipmt=0; ipmt<nbbcde; ipmt++){
          float val = bbcde_et[ipmt];
          float bbc_phi = bbcde_phi[ipmt];
          float bbc_eta = bbcde_eta[ipmt];

          float dphi = trackbuff_phi[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff] - bbc_phi;
          float deta = trackbuff_eta[icent*vzbin*nbuff*ntrks+ivz*nbuff*ntrks+ibuff*ntrks+itrk_buff] - bbc_eta;
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
          kbacknorthetabbc[icent][ipt_cor]->Fill(deta, dphia);
          
          hbacknorthbbcw[icent][ipt_cor]->Fill(dphi,val*val);
          kbacknorthbbcw[icent][ipt_cor]->Fill(dphia,val);
          kbacknorthetabbcw[icent][ipt_cor]->Fill(deta, dphia, val);

        }//bbcnorth

      }
    }
  }
  
  //3) trigger track + mixing tower
  for(int itrk=0; itrk<ntrk; itrk++){
    int ipt_cor = track_pt[itrk]/0.5;

    if(dbbcau_buff[icent][ivz]>0){
      for(int ibuff=0;ibuff<dbbcau_buff[icent][ivz];ibuff++){
        for(int ipmt_buff=0; ipmt_buff<buff_nbbcau[icent][ivz][ibuff]; ipmt_buff++){
          float val = bbcaubuff_et[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ipmt_buff];
          float dphi = track_phi[itrk] - bbcaubuff_phi[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ipmt_buff];
          float deta = track_eta[itrk] - bbcaubuff_eta[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ipmt_buff];
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
          kbacksouthetabbc2[icent][ipt_cor]->Fill(deta, dphia);

          hbacksouthbbcw2[icent][ipt_cor]->Fill(dphi,val*val);
          kbacksouthbbcw2[icent][ipt_cor]->Fill(dphia,val);
          kbacksouthetabbcw2[icent][ipt_cor]->Fill(deta, dphia, val);
        }
      }
    }//end bbc south

    if(dbbcde_buff[icent][ivz]>0){
      for(int ibuff=0;ibuff<dbbcde_buff[icent][ivz];ibuff++){
        for(int ipmt_buff=0; ipmt_buff<buff_nbbcde[icent][ivz][ibuff]; ipmt_buff++){
          float val = bbcdebuff_et[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ipmt_buff];
          float dphi = track_phi[itrk] - bbcdebuff_phi[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ipmt_buff];
          float deta = track_eta[itrk] - bbcdebuff_eta[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ipmt_buff];
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
          kbacknorthetabbc2[icent][ipt_cor]->Fill(deta, dphia);
          
	  hbacknorthbbcw2[icent][ipt_cor]->Fill(dphi,val*val);
          kbacknorthbbcw2[icent][ipt_cor]->Fill(dphia,val);
          kbacknorthetabbcw2[icent][ipt_cor]->Fill(deta, dphia, val);
        }
      }
    }//end bbc north
    
  }//end of 3 

  //4) buff track
  if(ntrk>0){//rebuff
    
    int ibuff=ntrack_buff[icent][ivz];
    
    int itrk_buff=0;//track number in each buff

    if(ntrk>0){//buffing lambda
      for(int itrk=0;itrk<ntrk;itrk++){
	
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

  //6) buff bbcau
  if(nbbcau>0){//rebuff
    int ibuff=nbbcau_buff[icent][ivz];
    
    int ibbc_buff=0;//bbcau number in each buff

    if(nbbcau>0){//buffing lambda
      for(int ibbc=0;ibbc<nbbcau;ibbc++){
	bbcaubuff_et[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ibbc_buff]=bbcau_et[ibbc];
	bbcaubuff_phi[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ibbc_buff]=bbcau_phi[ibbc];
        bbcaubuff_eta[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ibbc_buff]=bbcau_eta[ibbc];

        ibbc_buff++;
      }
      if(ibbc_buff>0){
        buff_nbbcau[icent][ivz][ibuff]=ibbc_buff;
        nbbcau_buff[icent][ivz]++;
        if(dbbcau_buff[icent][ivz]<nbuff) dbbcau_buff[icent][ivz]++;
        if(nbbcau_buff[icent][ivz]>=nbuff) nbbcau_buff[icent][ivz]=0;
      }
    }
  }

  if(nbbcde>0){//rebuff
    int ibuff=nbbcde_buff[icent][ivz];

    int ibbc_buff=0;//bbcde number in each buff

    if(nbbcde>0){//buffing lambda
      for(int ibbc=0;ibbc<nbbcde;ibbc++){
        bbcdebuff_et[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ibbc_buff]=bbcde_et[ibbc];
        bbcdebuff_phi[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ibbc_buff]=bbcde_phi[ibbc];
        bbcdebuff_eta[icent*vzbin*nbuff*nbbcs+ivz*nbuff*nbbcs+ibuff*nbbcs+ibbc_buff]=bbcde_eta[ibbc];

        ibbc_buff++;
      }
      if(ibbc_buff>0){
        buff_nbbcde[icent][ivz][ibuff]=ibbc_buff;
        nbbcde_buff[icent][ivz]++;
        if(dbbcde_buff[icent][ivz]<nbuff) dbbcde_buff[icent][ivz]++;
        if(nbbcde_buff[icent][ivz]>=nbuff) nbbcde_buff[icent][ivz]=0;
      }
    }
  }
}
  return 0;
}   
//_____________________________________________________________________________________________________________________________
int Ridgecntbbc::End()
{
    std::cout << "InputFileName = " << InputFileName << std::endl;
    std::cout << "OutputFileName = " << OutputFileName << std::endl;
    

  if(d_outfile) {
    d_outfile->cd();

    for(int icent=0; icent<centbin; icent++){
      //BBC
      for(int ipt=0; ipt<ptbin; ipt++){
	hforesouthbbc[icent][ipt]->Write();
	hbacksouthbbc[icent][ipt]->Write();
	
	hforenorthbbc[icent][ipt]->Write();
	hbacknorthbbc[icent][ipt]->Write();
	
	kforesouthbbc[icent][ipt]->Write();
        kbacksouthbbc[icent][ipt]->Write();
	kforesouthetabbc[icent][ipt]->Write();
        kbacksouthetabbc[icent][ipt]->Write();

        kforenorthbbc[icent][ipt]->Write();
        kbacknorthbbc[icent][ipt]->Write();
        kforenorthetabbc[icent][ipt]->Write();
        kbacknorthetabbc[icent][ipt]->Write();


	hforesouthbbcw[icent][ipt]->Write();
	hbacksouthbbcw[icent][ipt]->Write();
	
	hforenorthbbcw[icent][ipt]->Write();
	hbacknorthbbcw[icent][ipt]->Write();
	
	kforesouthbbcw[icent][ipt]->Write();
        kbacksouthbbcw[icent][ipt]->Write();
	kforesouthetabbcw[icent][ipt]->Write();
        kbacksouthetabbcw[icent][ipt]->Write();

        kforenorthbbcw[icent][ipt]->Write();
        kbacknorthbbcw[icent][ipt]->Write();
        kforenorthetabbcw[icent][ipt]->Write();
        kbacknorthetabbcw[icent][ipt]->Write();

	//tow mixing
	hbacksouthbbc2[icent][ipt]->Write();
	hbacknorthbbc2[icent][ipt]->Write();
        kbacksouthbbc2[icent][ipt]->Write();
        kbacknorthbbc2[icent][ipt]->Write();
        kbacksouthetabbc2[icent][ipt]->Write();
        kbacknorthetabbc2[icent][ipt]->Write();

	hbacksouthbbcw2[icent][ipt]->Write();
	hbacknorthbbcw2[icent][ipt]->Write();
        kbacksouthbbcw2[icent][ipt]->Write();
        kbacknorthbbcw2[icent][ipt]->Write();
        kbacksouthetabbcw2[icent][ipt]->Write();
        kbacknorthetabbcw2[icent][ipt]->Write();
      }
    }

    //htree->Write();
    
    //d_outfile->Close();
    //delete d_outfile;
    
  }
  else{

      cout<<"ERROR: No output file set!"<<endl;

  }

  return 0;
}
