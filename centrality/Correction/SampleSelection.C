#include <stdio.h>
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TTree.h"
#include <TMath.h>
#define Pi 3.14159265359

  const Int_t maxn = 900000;
  Float_t b, sumethf, sume3gev;
  Float_t npart, ncoll;
  Int_t nhreco, nhgen, nhf;
  Float_t ptreco[maxn], etareco[maxn], phireco[maxn];
  Float_t ptgen[maxn], etagen[maxn], phigen[maxn];
  Float_t ethf[maxn], ehf[maxn], etahf[maxn], phihf[maxn];
  Float_t gb, gsumethf, gsume3gev;
  Float_t gnpart, gncoll;
  Int_t gnhreco, gnhgen, gnhf;
  Float_t gptreco[maxn], getareco[maxn], gphireco[maxn];
  Float_t gptgen[maxn], getagen[maxn], gphigen[maxn];
  Float_t gethf[maxn], gehf[maxn], getahf[maxn], gphihf[maxn];

Bool_t highPurity[maxn];
Float_t a1[maxn], a2[maxn],a3[maxn],a4[maxn],a5[maxn],a6[maxn];
Int_t nhNtrk;
Float_t sume3gev;

Int_t HLT_PAZeroBiasPixel_SingleTrack_v1=0;
Int_t phltPixelClusterShapeFilter=0;
Int_t phfPosFilter1=0;
Int_t phfNegFilter1=0;
Int_t pprimaryvertexFilter=0;
Int_t Run=0;
Int_t CentBin=-1;
Int_t lumi=0;
Float_t vz;

  long Nevt;

void centTree1(){

  //TFile * file = TFile::Open("PAPhysics-pt04-test.root");
  TFile * file = TFile::Open("/tmp/tuos/pPb_MinBiasTree_v4_210498-210906_json.root");
  TTree *tree;
  tree=(TTree*)hiEvtAnalyzer->Get("HiTree");
  //tree->AddFriend("evttree=hiEvtAnalyzer/HiTree");
  //tree->AddFriend("pftree=pfTowers/hf");
  tree->AddFriend("pptrack=pptracks/trackTree");
  tree->AddFriend("skimTree=skimanalysis/HltTree");
  tree->AddFriend("hltTree=hltanalysis/HltTree");

//     tree->SetBranchAddress("mult",&nhgen);
     //tree->SetBranchAddress("b",&b);
     tree->SetBranchAddress("vz",&vz);
//     tree->SetBranchAddress("Ncoll",&ncoll);
     //tree->SetBranchAddress("Npart",&npart);
//     tree->SetBranchAddress("pt",&ptgen);
//     tree->SetBranchAddress("eta",&etagen);
//     tree->SetBranchAddress("phi",&phigen);
Int_t xxxx;
     tree->SetBranchAddress("nTrk",&nhreco);
     tree->SetBranchAddress("trkPt",&ptreco);
     tree->SetBranchAddress("trkEta",&etareco);
     tree->SetBranchAddress("trkPhi",&phireco);
     tree->SetBranchAddress("trkDz1",&a1);
     tree->SetBranchAddress("trkDzError1",&a2);
     tree->SetBranchAddress("trkDxy1",&a3);
     tree->SetBranchAddress("trkDxyError1",&a4);
     tree->SetBranchAddress("trkPtError",&a5);
     tree->SetBranchAddress("highPurity",&highPurity);
     //tree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
     tree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&xxxx);
     tree->SetBranchAddress("phltPixelClusterShapeFilter",&phltPixelClusterShapeFilter);
     tree->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
     tree->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
     tree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
     //tree->SetBranchAddress("Run",&Run);
     tree->SetBranchAddress("hiBin",&CentBin);
     tree->SetBranchAddress("lumi",&lumi);

     Nevt=tree->GetEntries();  cout<<"Entries= "<<Nevt<<endl;

const int Ncentb=100;
const int Ncutb=10;
double bb[Ncutb+1]={0,5,10,15,20,30,40,50,60,70,100};
double Ntrks[Ncutb]={0};
int tmpb=-1;
int allntrks=0;
long int evtused=0;
     for(long ne=0; ne<Nevt; ne++){
     //for(long ne=0; ne<5000; ne++){
       if(ne%3000==0)  cout<<"Have run "<<ne<<" events. "<<endl;
       tree->GetEntry(ne);

//if(!(phltPixelClusterShapeFilter&&phfPosFilter1&&phfNegFilter1&&pprimaryvertexFilter&&HLT_PAZeroBiasPixel_SingleTrack_v1)) continue;
if(!(phltPixelClusterShapeFilter&&phfPosFilter1&&phfNegFilter1&&pprimaryvertexFilter&&xxxx&&fabs(vz)<15)) continue;
evtused++;
tmpb=-1;
      for(int i=0;i<Ncutb;i++){
        if(bb[i]<=CentBin&&bb[i+1]>CentBin) tmpb=i;
      }
     
       int nhtmp=0;
       for(int ni=0;ni<nhreco;ni++){
       if(highPurity[ni]&&fabs(a1[ni]/a2[ni])<3&&fabs(a3[ni]/a4[ni])<3&&a5[ni]/ptreco[ni]<0.1&&ptreco[ni]>10.0&&fabs(etareco[ni])<2.4){
         //gptreco[nhtmp]=ptreco[ni]; getareco[nhtmp]=etareco[ni]; gphireco[nhtmp]=phireco[ni];
          nhtmp++;
          Ntrks[tmpb]=Ntrks[tmpb]+1;
        }
       }
       //nhNtrk=nhtmp;
       allntrks+=nhtmp;
     }

for(int i=0;i<Ncutb;i++)
cout<<""<<Ntrks[i]<<endl;

cout<<endl<<allntrks<<endl;
cout<<endl<<evtused<<endl;

     //fctree->Close();
     delete file;
}






void RecentTree(){

centTree1();
//centTree2();

}


