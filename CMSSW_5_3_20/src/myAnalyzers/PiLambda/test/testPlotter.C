#include <cstdlib> 
#include <iostream>
#include "xiPlotter.h"
#include <TApplication.h>

int main(int argc, char **argv){
  //TApplication* theApp = new TApplication("Rint", &argc, argv);
  
  if(argc != 9){														     
   std::cout << "Missing parameters. Exiting" << std::endl;									     
   std::cout << "Usage: doEff  outputDir xi3dCut batpi3dCut v0pi3dCut v0proton3dCut xisepCut vtxClCut \n";  
   exit(0);															     
  }																     
  																     
  float xi3dCut = atof(argv[3]);												     
  float batpi3dCut = atof(argv[4]);												     
  float v0pi3dCut = atof(argv[5]);												     
  float v0proton3dCut = atof(argv[6]);												     
  float xisepCut = atof(argv[7]);												     
  float vtxClCut = atof(argv[8]);												     

 
   TFile* the900GeVFile = new TFile("/nfs/data37/cms/v0ntuple/data_900GeV_june19.root");					     
   TTree* the900GeVTree = (TTree*)the900GeVFile->Get("ntupleXiOm/ntuple");							     
  
   xiPlotter theXi900GeVPlotter(the900GeVTree); 										     
   theXi900GeVPlotter.SetNEvents(-1);												     
   theXi900GeVPlotter.SetOutputDir(argv[1]); 								     
   theXi900GeVPlotter.SetEffFile(argv[2]+std::string("/efficienciesXi_900GeV.root"));					     
   theXi900GeVPlotter.SetStyle();												     
   theXi900GeVPlotter.SetCutValues(xi3dCut, batpi3dCut, v0pi3dCut, v0proton3dCut, xisepCut, vtxClCut);
   theXi900GeVPlotter.Loop();													     
   theXi900GeVPlotter.WriteHistos();												     

/*

    TFile* the7TeVFile = new TFile("/nfs/data37/cms/v0ntuple/data_7TeV_june22.root");
    TTree* the7TeVTree = (TTree*)the7TeVFile->Get("ntupleXiOm/ntuple");
  
    xiPlotter the7TeVXiPlotter(the7TeVTree);
    the7TeVXiPlotter.SetNEvents(-1);
    the7TeVXiPlotter.SetEffFile(argv[2]+std::string("/efficienciesXi_7TeV.root"));
    the7TeVXiPlotter.SetOutputDir(argv[1]);
    the7TeVXiPlotter.SetStyle();
    the7TeVXiPlotter.SetCutValues(xi3dCut, batpi3dCut, v0pi3dCut, v0proton3dCut, xisepCut, vtxClCut);
    the7TeVXiPlotter.Loop();
    the7TeVXiPlotter.WriteHistos();
*/
 

 
  return 0;
}
