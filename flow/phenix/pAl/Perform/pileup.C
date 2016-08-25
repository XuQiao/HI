#include <stdlib.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include "TRandom.h" 
#include "TRandom3.h" 
#include "TString.h"
#include <fstream>
#include <iostream>

std::string readline(char* name, int iline);
int getnumber(char* name);

void pileup()
{
  TH1::SetDefaultSumw2();   
  int trig;
  float bbcv;
  int cent;
  int nbbc;
  float bbct0[128];
  float bbccharge[128];
  float bbcz[128];
  
  TChain *ch = new TChain("tree","A Chain");
  TChain *ch1 = new TChain("tree","A Chain");
  char infile[] = "../filelistmb.dat";
  char infile1[] = "../filelistmb1.dat";
  cout << "Adding Chain..." << endl;
 // for(int i=0;i < getnumber(infile);i++){
  for(int i=0;i < 10;i++){
      if(i % 10 == 0) cout << i << " file;" << endl;
      TString filename(readline(infile,i));
      TString filename1(readline(infile1,i));
      ch->Add(Form("%s",filename.Data()));
      ch1->Add(Form("%s",filename1.Data()));
  }
  cout<<ch->GetEntries()<<endl;
  int nEvent = ch->GetEntries();
  ch -> SetBranchAddress("trig", &trig);
  ch -> SetBranchAddress("bbcv", &bbcv);
  ch1 -> SetBranchAddress("cent", &cent);
  ch -> SetBranchAddress("nbbc", &nbbc);
  ch1 -> SetBranchAddress("bbct0", &bbct0);
  ch -> SetBranchAddress("bbccharge",&bbccharge);
  ch -> SetBranchAddress("bbcz",&bbcz);

  const int ncent = 3;
  const int nhisto = 20;
  TH1F* bbcdeltat0s[ncent][nhisto];
  TH1F* bbcdeltat0n[ncent][nhisto];
  TH1F* bbcmixdeltat0s[nhisto];
  TH1F* bbcmixdeltat0n[nhisto];
  char name[128];
  for(int ihisto=0;ihisto<nhisto;ihisto++){
      sprintf(name,"bbcmixdeltat0s_%d",ihisto);
      bbcmixdeltat0s[ihisto] = new TH1F(name,name,200,-5,5);
      sprintf(name,"bbcmixdeltat0n_%d",ihisto);
      bbcmixdeltat0n[ihisto] = new TH1F(name,name,200,-5,5);
  for(int icent=0;icent<ncent;icent++){
      sprintf(name,"bbcdeltat0s_%d_%d",icent,ihisto);
      bbcdeltat0s[icent][ihisto] = new TH1F(name,name,200,-5,5);
      sprintf(name,"bbcdeltat0n_%d_%d",icent,ihisto);
      bbcdeltat0n[icent][ihisto] = new TH1F(name,name,200,-5,5);
  }
  }
  int nr[ncent]={};
  int nmix=0;
  TRandom3 *r = new TRandom3(0);
  int ievent = r->Integer(nEvent);
  int istep = r->Integer(100);
  int ibuff = 0;
  int nbuff = 10;
  vector<float> buff_bbct0s;
  vector<float> buff_bbct0n;
  vector<float> buff_bbccharge;
  vector<float> buff_bbcz;
  while(nr[0]<nhisto || nr[1]<nhisto || nr[2]<nhisto || nmix < nhisto){
//    int ievent = r->Integer(nEvent);
    ch -> GetEntry(ievent);
    ch1 -> GetEntry(ievent);
    ievent += istep;
    bool isMB = (trig & 0x00000010);
    bool iscentral = (trig & 0x00000008);
    //if (!isMB) continue;
    if (!iscentral) continue;
    if (fabs(bbcv)>10) continue;
    cout<<ievent<<"\t"<<nr[0]<<"\t"<<nr[1]<<"\t"<<nr[2]<<"\t"<<nmix<<endl;
    float bbctmeans=0, bbctmeann=0;
    if(nmix < nhisto){
        if(ibuff % nbuff == 0){
        bbctmeans=0;
        int nbbcau=0;
        bbctmeann=0;
        int nbbcde=0;
        for(int ipmt=0;ipmt<nbbc;ipmt++){
            if (bbccharge[ipmt]>0) {
            int iarm = 0;
            if (bbcz[ipmt] > 0) iarm = 1;
            if(iarm==0){
            bbctmeans+=bbct0[ipmt];
            nbbcau++;
            }
            else{
            bbctmeann+=bbct0[ipmt];
            nbbcde++;
            }
            }
        }
        for(int ipmt=0;ipmt<buff_bbct0s.size();ipmt++){
            bbctmeans+=buff_bbct0s[ipmt];
        }
        for(int ipmt=0;ipmt<buff_bbct0n.size();ipmt++){
            bbctmeann+=buff_bbct0s[ipmt];
        }
        bbctmeans/=(buff_bbct0s.size()+nbbcau);
        bbctmeann/=(buff_bbct0n.size()+nbbcde);
        for(int ipmt=0;ipmt<buff_bbct0s.size();ipmt++){
           // bbcmixdeltat0s[nmix]->Fill(buff_bbct0s[ipmt]-bbct0[0]);
            bbcmixdeltat0s[nmix]->Fill(buff_bbct0s[ipmt]-bbctmeans);
        }
        for(int ipmt=0;ipmt<buff_bbct0n.size();ipmt++){
           // bbcmixdeltat0n[nmix]->Fill(buff_bbct0n[ipmt]-bbct0[0]);
            bbcmixdeltat0n[nmix]->Fill(buff_bbct0n[ipmt]-bbctmeann);
        }
        buff_bbct0s.clear();
        buff_bbct0n.clear();
        for(int ipmt=0;ipmt<nbbc;ipmt++){
            if (bbccharge[ipmt]>0) {
            int iarm = 0;
            if (bbcz[ipmt] > 0) iarm = 1;
            if(iarm==0){
                //bbcmixdeltat0s[nmix]->Fill(bbct0[ipmt]-bbct0[0]);
                bbcmixdeltat0s[nmix]->Fill(bbct0[ipmt]-bbctmeans);
                buff_bbct0s.push_back(bbct0[ipmt]);
            }
            else{
                //bbcmixdeltat0n[nmix]->Fill(bbct0[ipmt]-bbct0[0]);
                bbcmixdeltat0n[nmix]->Fill(bbct0[ipmt]-bbctmeann);
                buff_bbct0n.push_back(bbct0[ipmt]);
            }
            }
        }
        nmix++;
        }
        ibuff ++;
    }
    if(nr[0] < nhisto && cent<5 && nbbc>55){
        bbctmeans=0;
        nbbcau=0;
        bbctmeann=0;
        nbbcde=0;
        for(int ipmt=0;ipmt<nbbc;ipmt++){
            if (bbccharge[ipmt]>0) {
            int iarm = 0;
            if (bbcz[ipmt] > 0) iarm = 1;
            if(iarm==0){
            bbctmeans+=bbct0[ipmt];
            nbbcau++;
            }
            else{
            bbctmeann+=bbct0[ipmt];
            nbbcde++;
            }
            }
        }
        for(int ipmt=0;ipmt<nbbc;ipmt++){
            if (bbccharge[ipmt]>0) {
            int iarm = 0;
            if (bbcz[ipmt] > 0) iarm = 1;
            if(iarm==0)
                //bbcdeltat0s[0][nr[0]]->Fill(bbct0[ipmt]-bbct0[0]);
                bbcdeltat0s[0][nr[0]]->Fill(bbct0[ipmt]-bbctmeans/nbbcau);
            else 
                //bbcdeltat0n[0][nr[0]]->Fill(bbct0[ipmt]-bbct0[0]);
                bbcdeltat0n[0][nr[0]]->Fill(bbct0[ipmt]-bbctmeann/nbbcde);
            }
        }
        nr[0]++;
    }

    if(nr[1] < nhisto && cent>5 && cent<10 && nbbc > 40){
        bbctmeans=0;
        nbbcau=0;
        bbctmeann=0;
        nbbcde=0;
        for(int ipmt=0;ipmt<nbbc;ipmt++){
            if (bbccharge[ipmt]>0) {
            int iarm = 0;
            if (bbcz[ipmt] > 0) iarm = 1;
            if(iarm==0){
            bbctmeans+=bbct0[ipmt];
            nbbcau++;
            }
            else{
            bbctmeann+=bbct0[ipmt];
            nbbcde++;
            }
            }
        }
        for(int ipmt=0;ipmt<nbbc;ipmt++){
            if (bbccharge[ipmt]>0) {
            int iarm = 0;
            if (bbcz[ipmt] > 0) iarm = 1;
            if(iarm==0)
                //bbcdeltat0s[1][nr[1]]->Fill(bbct0[ipmt]-bbct0[0]);
                bbcdeltat0s[1][nr[1]]->Fill(bbct0[ipmt]-bbctmeans/nbbcau);
            else 
                //bbcdeltat0n[1][nr[1]]->Fill(bbct0[ipmt]-bbct0[0]);
                bbcdeltat0n[1][nr[1]]->Fill(bbct0[ipmt]-bbctmeann/nbbcde);
            }
        }
        nr[1]++;
    }

    if(nr[2] < nhisto && cent>10 && cent<20 && nbbc > 30){
        bbctmeans=0;
        nbbcau=0;
        bbctmeann=0;
        nbbcde=0;
        for(int ipmt=0;ipmt<nbbc;ipmt++){
            if (bbccharge[ipmt]>0) {
            int iarm = 0;
            if (bbcz[ipmt] > 0) iarm = 1;
            if(iarm==0){
            bbctmeans+=bbct0[ipmt];
            nbbcau++;
            }
            else{
            bbctmeann+=bbct0[ipmt];
            nbbcde++;
            }
            }
        }
        for(int ipmt=0;ipmt<nbbc;ipmt++){
            if (bbccharge[ipmt]>0) {
            int iarm = 0;
            if (bbcz[ipmt] > 0) iarm = 1;
            if(iarm==0)
            //    bbcdeltat0s[2][nr[2]]->Fill(bbct0[ipmt]-bbct0[0]);
                bbcdeltat0s[2][nr[2]]->Fill(bbct0[ipmt]-bbctmeans/nbbcau);
            else 
           //     bbcdeltat0n[2][nr[2]]->Fill(bbct0[ipmt]-bbct0[0]);
                bbcdeltat0n[2][nr[2]]->Fill(bbct0[ipmt]-bbctmeann/nbbcde);
            }
        }
        nr[2]++;
    }
  }

 TFile *fout = new TFile("RndmEventdeltaT0_Mean.root","Recreate");
 fout->cd();
  for(int ihisto=0;ihisto<nhisto;ihisto++){
    bbcmixdeltat0s[ihisto]->Write();
    bbcmixdeltat0n[ihisto]->Write();
  for(int icent=0;icent<ncent;icent++){
    bbcdeltat0s[icent][ihisto]->Write();
    bbcdeltat0n[icent][ihisto]->Write();
  }
  }
}

std::string readline(char* name, int iline){
        std::ifstream backstory(name);
        std::string line;
        if (backstory.is_open())
                if(backstory.good()){
                        for(int i = 0; i < iline+1; ++i)
                           getline(backstory, line);
                }
        return line;
}

int getnumber(char* name){
        std::ifstream backstory(name);
        std::string line;
        int i=0;
        while (!backstory.eof())
            if(backstory.good()){
                getline(backstory, line);
                i++;
            }
        return i-1;
}
