/*
Global and forward
Centrality
HF energy vs η, ET
Number of vertices
Correlation of vertices of the foreground and background
Vertex distribution
Track multiplicity distributions
Event content
*/
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include <vector>

using namespace std;
class IntVar{
 public:
    IntVar();
    ~IntVar();
    void Init(vector<TString> sin);
    void Fill();
    void Write(TString sout);

 private:
    TChain *tree;
    float HFp3;
    float HFp4;
    float HFm3;
    float HFm4;
    float ET;
    float eta;
    float vtz;
    int Bin;
    int Ntrack;
    float vz;
    float HFhit;
    int HLT_PAZeroBiasPixel_SingleTrack_v1;
    int phltPixelClusterShapeFilter;
    int phfPosFilter1;
    int phfNegFilter1;
    int pprimaryvertexFilter;
    int pBeamScrapingFilter;
    int pVertexFilterCutGplus;

    vector<TH1F*> h1;
    vector<TH2F*> h2;
};
