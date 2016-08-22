#include <cstdlib> 
#include <iostream>
#include "seleXiMC.h"
#include <TApplication.h>
#include "/usr/users/eluiggi/bin/include/config.h"

// #### Weights 7TeV: 1 -> D6T, 2 -> pythia8, 3 -> p0 
// #### Weights 900GeV: 11 -> D6T, 12 -> pythia8, 13 -> p0  

int main(int argc, char **argv){
   Config theConfig(argv[1], argv);

   std::string theInFile = theConfig.pString("InputRootFile");
   std::string theOutFile = theConfig.pString("OutputRootDir");
   std::string theEnergy = theConfig.pString("SqrtS");
   int nEvents = theConfig.pInt("NEvents");
   int theWeight = theConfig.pInt("Weight");
   std::string theWeightType = theConfig.pString("WeightType");
   double xi3dCut = theConfig.pDouble("Xi3dCut");
   double batpi3dCut = theConfig.pDouble("BatPi3dCut");
   double v0pi3dCut = theConfig.pDouble("V0Pi3dCut");
   double v0proton3dCut = theConfig.pDouble("V0Proton3dCut");
   double xisepCut = theConfig.pDouble("XiSepCut");
   double vtxClCut = theConfig.pDouble("VtxClCut");

   TFile* theMCFile = new TFile(theInFile.c_str());  					       
   TTree* theTree = (TTree*)theMCFile->Get("ntupleXiOm/ntuple");			       
   seleXiMC eff7TeV(theTree);								       
   eff7TeV.SetNEvents(nEvents); 								       
   eff7TeV.SetDataType("mc", theEnergy);							       
   eff7TeV.SetWeightChoice(theWeight, theWeightType);  // 1 = 7TeV, D6T, 2 = 7TeV p8, 3 = 7TeV P0, 11 = 900 GeV D6T,  
   eff7TeV.SetStyle();										       
   eff7TeV.SetOutputDir(theOutFile);								       
   eff7TeV.SetCutValues(xi3dCut, batpi3dCut, v0pi3dCut, v0proton3dCut, xisepCut, vtxClCut);	       
   eff7TeV.Loop();										       
   eff7TeV.WriteHistos();									       

  return 0;
}
