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
#include "TBranchElement.h"
#include <fstream>
#include <iostream>

using namespace std;

const int bbcz_cut = 10;
const int vtxz_cut = 10;

const int ntrks = 200;
const int nbbcs = 800;
const int ncluss = 800;

const int centbin = 10;
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
 n(0),
 px(),
 py(),
 pz(),
 E(),
 KS(),
 KF(),
 M()
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
  d_outfile = new TFile(OutputFileName,"recreate");
      d_infile = TFile::Open(InputFileName,"ReadOnly");
  if(!d_infile) return 0;
  tree = (TTree*)d_infile->Get("T");
  if(!tree) return 0;
  TBranchElement *array = (TBranchElement*)tree->GetBranch("part_array");
  tree->SetMakeClass(1);
    tree->SetBranchAddress("part_array.fKS",KS);
    tree->SetBranchAddress("part_array.fKF",KF);
    tree->SetBranchAddress("part_array.fPx",px);
    tree->SetBranchAddress("part_array.fPy",py);
    tree->SetBranchAddress("part_array.fPz",pz);
    tree->SetBranchAddress("part_array.fEnergy",E);
    tree->SetBranchAddress("part_array.fMass",M);
    array->SetAddress(&n);
  TH1::SetDefaultSumw2(kTRUE);
  
char name[512];
 pi = acos(-1);

      sprintf(name,"hnch");
      hnch = new TH1F(name,name, 500, 0,500);
    
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
  if(n>10000)  {cout<<n<<"!"<<endl;  continue;}
  int icent = -1;
  int ntrk = 0;
  int ivz = 0;
  int nch = 0;

  float track_pt[ntrks];
  float track_phi[ntrks];
  float track_eta[ntrks];
  int nbbcau=0;//bbc south
  float bbcau_et[nbbcs];
  float bbcau_phi[nbbcs];
  float bbcau_eta[nbbcs];

  int nbbcde=0;//bbc north
  float bbcde_et[nbbcs];
  float bbcde_phi[nbbcs];
  float bbcde_eta[nbbcs];

        for(int iparticle = 0;iparticle < n; iparticle++){
            float pt = sqrt(px[iparticle]*px[iparticle]+py[iparticle]*py[iparticle]);
            float theta = atan(pt/pz[iparticle]);
            if(theta<0) theta = TMath::Pi()+theta;
            float eta = -log(tan(theta/2));
            float phi;
            if(px[iparticle]>0){
                phi = atan(py[iparticle]/px[iparticle]);
            if(py[iparticle]<0)
                phi = 2*TMath::Pi()+atan(py[iparticle]/px[iparticle]);
            }
            else
                phi = TMath::Pi()+atan(py[iparticle]/px[iparticle]);
            float Energy = E[iparticle];
        if(!(fabs(KF[iparticle]) == 211 || fabs(KF[iparticle]) == 213 || fabs(KF[iparticle]) == 321 || fabs(KF[iparticle]) == 323 ||fabs(KF[iparticle]) == 2212)) continue;
           // if(KS[iparticle]!=1) continue;
          if(fabs(eta)<0.35 && pt>0.2 && pt< 5.0){
              track_pt[ntrk] = pt;
              track_phi[ntrk] = phi;
              track_eta[ntrk] = eta;
              ntrk++;
          }
            if(fabs(eta)>3 && fabs(eta)<4 && pt>0.05){
      int iarm = 0;
      if (eta > 0) iarm = 1;
	if(iarm==0){
	  bbcau_et[nbbcau] = Energy;
	  bbcau_phi[nbbcau] = phi; 
	  bbcau_eta[nbbcau] = eta; 
	  nbbcau++;
	}
	else{
	  bbcde_et[nbbcde] = Energy;
          bbcde_phi[nbbcde] = phi;
          bbcde_eta[nbbcde] = eta;
          nbbcde++;
	}
         if(iarm==0)   nch++;
    }
        //if((fabs(KF[iparticle]) == 211 || fabs(KF[iparticle]) == 321 || fabs(KF[iparticle]) == 2212)) 
    }

   // int chbin[centbin+1] = {0,1,2,3,4,5,7,10,14,18,40}; //both side
   // int chbin[centbin+1] = {0,1,2,3,4,5,7,10,12,16,40}; //pAu
   // int chbin[centbin+1] = {0,1,2,3,5,7,9,12,16,21,100}; //dAu
    int chbin[centbin+1] = {0,2,4,6,8,11,16,20,26,30,100}; //HeAu
    hnch->Fill(nch);
    for(int ich=0;ich<centbin;ich++){
        if(nch>chbin[ich] && nch<=chbin[ich+1]){
            icent = ich;
            break;
        }
    }
    if(icent<0 || icent>=centbin) continue;

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

    hnch->Write();

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
    
    d_outfile->Close();
    //delete d_outfile;
    
  }
  else{

      cout<<"ERROR: No output file set!"<<endl;

  }

  return 0;
}
