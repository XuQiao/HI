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

void pileupm2()
{
  TH1::SetDefaultSumw2();   
  TChain *ch = new TChain("tree","A Chain");
  char infile[] = "../filelistmb.dat";
  cout << "Adding Chain..." << endl;
  for(int i=0;i < getnumber(infile);i++){
      if(i % 10 == 0) cout << i << " file;" << endl;
      TString filename(readline(infile,i));
      ch->Add(Form("%s",filename.Data()));
  }

  int trig;
  float ZdcEs;
  float ZdcEn;
  vector<float> buff_Zdcs;
  vector<float> buff_bbcs;
  double trigbin[30];
  for(int i=0;i<30;i++){
      trigbin[i] = 2**i-1+0.5;
  }
//  TH2F *hmixbbcsZdcEs = new TH2F("hmixbbcsZdcEs","hmixbbcsZdcEs",600,0,600,2500,0,5000);
  TH1F* htrig = new TH1F("htrig","htrig",29,trigbin);

  ch -> SetBranchAddress("trig", &trig);
  ch -> SetBranchAddress("ZdcEs", &ZdcEs);
  ch -> SetBranchAddress("ZdcEn", &ZdcEn);
  long nEvent = ch->GetEntries();
  cout << nEvent << endl;
  /*
  for(int ievent=0;ievent < nEvent; ievent++){
      ArrayZdcEs[ievent] = ZdcEs;
      ArrayZdcEn[ievent] = ZdcEn;
      Arraytrig[ievent] = trig;
  }
*/
  for(long ievent=0;ievent < 10000000; ievent++){
      ch->GetEntry(ievent);
  htrig -> Fill(trig);
 // if(ievent > 10000000) continue;
 
      if(ievent % 1000000 == 0 ) cout<< "Processing " << ievent << " events" << endl;
      /*
      float ZdcEs1 = ArrayZdcEs[ievent];
      float ZdcEn1 = ArrayZdcEn[ievent];
      //cout << "Zdc south energy = "<< ZdcEs1 << endl;
      bool isMB = (Arraytrig[ievent] & 0x00000010);
      bool isCentral = (Arraytrig[ievent] & 0x00000008);
      
      if(isCentral){
          hZdcEs -> Fill(ZdcEs1);
          hZdcEn -> Fill(ZdcEn1);
      }

      int snEvent = r->Integer(nEvent);
          //nEvent - ievent; 
      float ZdcEs2 = ArrayZdcEs[snEvent];
      //cout << "Zdc south energy other = "<< ZdcEs2 << endl;

      if(isMB) 
          hmixZdcEs -> Fill(ZdcEs1 + ZdcEs2);
        */
    }


    TFile *fout = new TFile("test.root","Recreate");
    fout->cd();
    htrig->Write();
 //   hZdcEs->Write();
 //   hZdcEn->Write();
 //   hmixZdcEs->Write();
    fout->Close();
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
