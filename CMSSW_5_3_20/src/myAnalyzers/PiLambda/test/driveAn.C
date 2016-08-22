#include <cstdlib> 
#include <iostream>
#include "seleXiData.h"
#include "/usr/users/eluiggi/bin/include/config.h"

int main(int argc, char **argv){
  Config theConfig(argv[1], argv);
  
  std::string theInFile = theConfig.pString("InputRootFile");															     
  std::string theEffFile = theConfig.pString("EffFile");
  std::string theOutRootDir = theConfig.pString("OutputRootDir");
  int nEvents = theConfig.pInt("NEvents");
  double xi3dCut = theConfig.pDouble("Xi3dCut");												     
  double batpi3dCut = theConfig.pDouble("BatPi3dCut");
  double v0pi3dCut = theConfig.pDouble("V0Pi3dCut");
  double v0proton3dCut = theConfig.pDouble("V0Proton3dCut");
  double xisepCut = theConfig.pDouble("XiSepCut");
  double vtxClCut = theConfig.pDouble("VtxClCut"); 


  TFile* the7TeVFile = new TFile(theInFile.c_str());						     
  TTree* the7TeVTree = (TTree*)the7TeVFile->Get("ntupleXiOm/ntuple");				     
  
  seleXiData the7TeVXiPlotter(the7TeVTree);							     
  the7TeVXiPlotter.SetNEvents(nEvents); 							     
  the7TeVXiPlotter.SetEffFile(theEffFile);							     
  the7TeVXiPlotter.SetOutputDir(theOutRootDir); 						     
  the7TeVXiPlotter.SetStyle();									     
  the7TeVXiPlotter.SetCutValues(xi3dCut, batpi3dCut, v0pi3dCut, v0proton3dCut, xisepCut, vtxClCut);  
  the7TeVXiPlotter.Loop();									     
  the7TeVXiPlotter.WriteHistos();								     
  return 0;
}
