#include <vector>
using namespace std;
void Writehisto(){
    IntVar *pl = new IntVar();
    vector<TString> sin;
    for(int ifile = 0; ifile < 191; ifile ++){
        TString sfile = Form("root://eoscms//store/group/phys_heavyions/kjung/pPb_EposMinBias_5TeV_8022_Forest_corrL1/EPOS5TeV_GEN_SIM/crab_pPb_EposMinBias_5TeV_8022_Forest_corrL1/161106_174912/0000/HiForestAOD_%d.root",ifile+1);
        sin.push_back(sfile);
    }
//    sin.push_back("/afs/cern.ch/user/k/kjung/public/forestTests/HiForestAOD_run284755_lumi414.root");
    pl -> Init(sin);
    pl -> Fill();
    pl -> Write("output_AODcut.root");
}
