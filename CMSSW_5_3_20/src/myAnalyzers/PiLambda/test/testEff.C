#include <cstdlib> 
#include <iostream>
#include "seleXi.h"
#include <TApplication.h>

// #### Weights 7TeV: 1 -> D6T, 2 -> pythia8, 3 -> p0 
// #### Weights 900GeV: 11 -> D6T, 12 -> pythia8, 13 -> p0  

int main(int argc, char **argv){
   if(argc != 10){
     std::cout << "Missing parameters. Exiting" << std::endl;
     std::cout << "Usage: doEff outputDir weight srqtS(900GeV/7TeV) xi3dCut batpi3dCut v0pi3dCut v0proton3dCut xisepCut vtxClCut \n";
     exit(0);
   }
   int theWeight = atoi(argv[2]);
   
   float xi3dCut = atof(argv[4]);
   float batpi3dCut = atof(argv[5]);
   float v0pi3dCut = atof(argv[6]);
   float v0proton3dCut = atof(argv[7]);
   float xisepCut = atof(argv[8]);
   float vtxClCut = atof(argv[9]);

  if (argv[3]==TString("7TeV")) {
    TFile* the7TeVMCFile = new TFile("/nfs/data37/cms/v0ntuple/mc_7TeV_tuneD6T_START36_V10_PARTIAL_latest.root");  
    //TFile* the7TeVMCFile = new TFile("/nfs/data37/cms/v0ntuple/mc_7TeV_pythia8_START36_V10.root");
    //TFile* the7TeVMCFile = new TFile("/nfs/data37/cms/v0ntuple/mc_7TeV_tuneP0_START36_V10.root");
    TTree* the7TeVMCTree = (TTree*)the7TeVMCFile->Get("ntupleXiOm/ntuple");
    seleXi eff7TeV(the7TeVMCTree);
    eff7TeV.SetNEvents(-1);
    eff7TeV.SetDataType("mc", argv[3]);
    eff7TeV.SetWeightChoice(theWeight);  // 1 = 7TeV, D6T, 2 = 7TeV p8, 3 = 7TeV P0, 11 = 900 GeV D6T, 
    eff7TeV.SetStyle();
    eff7TeV.SetOutputDir(argv[1]);
    eff7TeV.SetCutValues(xi3dCut, batpi3dCut, v0pi3dCut, v0proton3dCut, xisepCut, vtxClCut);
    eff7TeV.Loop();
    eff7TeV.WriteHistos();
  }  
  
  if (argv[3]==TString("900GeV")) {
    TFile* the900GeVMCFile = new TFile("/nfs/data37/cms/v0ntuple/mc_900GeV_pythia8_START36_V10A_PARTIAL.root");
    TTree* the900GeVMCTree = (TTree*)the900GeVMCFile->Get("ntupleXiOm/ntuple");
    seleXi eff900GeV(the900GeVMCTree);
    eff900GeV.SetNEvents(-1);
    eff900GeV.SetDataType("mc", argv[3]);
    eff900GeV.SetWeightChoice(theWeight);    //14 = mixed MC
    eff900GeV.SetStyle();
    eff900GeV.SetOutputDir(argv[1]);
    eff900GeV.SetCutValues(xi3dCut, batpi3dCut, v0pi3dCut, v0proton3dCut, xisepCut, vtxClCut);
    eff900GeV.Loop();
    eff900GeV.WriteHistos();
  }
  return 0;
}
