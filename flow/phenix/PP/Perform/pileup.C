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
  TChain *ch = new TChain("tree","A Chain");
  char infile[] = "../filelistfvtx.dat";
  cout << "Adding Chain..." << endl;
  for(int i=0;i < getnumber(infile);i++){
      if(i % 10 == 0) cout << i << " file;" << endl;
      TString filename(readline(infile,i));
      ch->Add(Form("%s",filename.Data()));
  }
  cout<<ch->GetEntries()<<endl;
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
