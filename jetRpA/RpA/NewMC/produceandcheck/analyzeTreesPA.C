#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <set>
#include "TH1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TChain.h"
#include "TF1.h"
#include "TMath.h"
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
//#include "/afs/cern.ch/work/y/ymao/public/RpA/Corrections/TrackCorrector2D.C"
//#include "/afs/cern.ch/work/y/ymao/public/RpA/Corrections/Correction.C"

//These includes cause some complications in CMSSW_5_3_8_HI_patchX. Commented out for pp.
//If you want to recalculate the JECs on the fly again, just uncomment everything in the updateJEC blocks
//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;
typedef std::vector<trigger::TriggerObject> trigO;

TStopwatch timer;

// ******* GLOBAL DECLARATIONS **********
const int QCDpthatBins = 8;
//const int QCDpthatBins = 9;
//const int QCDpthatBins = 1;
const int dataFiles = 10292;
int startfile ;
int endfile ;
int combinationMethod ;

const TString algo = "akPu3PF" ;
const double deta[]={-2.2, -1.2, -0.7, -0.3, 0.3, 0.7,1.2,2.2} ;
const int netabin = sizeof(deta)/sizeof(Double_t)-1 ;
const Double_t jetPtBin[]={3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;
const int pthatbin[10] = {15,30,50,80,120,170,220,280,370, 9999};
const double wght[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 1.088E-07, 0};

//***************************************


//**********************************************************
// Count the MC events to appropriately weight the pthat bins
//**********************************************************

int *countMCevents(std::string infile, int isMC){

  TChain *ch = NULL;
//  if(usePUsub) 
    ch = new TChain(Form("%sJetAnalyzer/t", algo.Data()));
//  else ch = new TChain("ak3PFJetAnalyzer/t");
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  int *MCentries = new int[QCDpthatBins];
  int entrySum[QCDpthatBins]; 
  for(int ifile=0; ifile<QCDpthatBins; ifile++){
    instr >> filename;
    ch->Add(filename.c_str());
  //  if(ifile==0) MCentries[ifile] = ch->GetEntries();
    if(ifile<QCDpthatBins-1) 
         entrySum[ifile] = ch->GetEntries(Form("pthat>=%d && pthat<%d",pthatbin[ifile], pthatbin[ifile+1]));
     else entrySum[ifile] = ch->GetEntries(Form("pthat>=%d",pthatbin[ifile]));  
    cout << "Summed MCentries[" << ifile << "]: " << entrySum[ifile] << endl;
    }
      for(int ifile=0; ifile<QCDpthatBins; ifile++){
     if(ifile==0) MCentries[ifile] = entrySum[ifile] ;
     else MCentries[ifile] = entrySum[ifile]; 
    //     MCentries[ifile] = ch->GetEntries(Form("pthat<%d",pthatbin[ifile+1]));
 //   MCentries[ifile] = ch->GetEntries();
  }
//  MCentries[0] = ch->GetEntries("pthat<15");
/*  MCentries[0] = ch->GetEntries("pthat>=15 && pthat<30");
  MCentries[1] = ch->GetEntries("pthat>=30 && pthat<50");
  MCentries[2] = ch->GetEntries("pthat>=50 && pthat<80");
  MCentries[3] = ch->GetEntries("pthat>=80 && pthat<120");
  MCentries[4] = ch->GetEntries("pthat>=120 && pthat<170");
  MCentries[5] = ch->GetEntries("pthat>=170 && pthat<220");
  MCentries[6] = ch->GetEntries("pthat>=220 && pthat<280");

  MCentries[0] = ch->GetEntries("pthat<30");
  MCentries[1] = ch->GetEntries("pthat<50");
  MCentries[2] = ch->GetEntries("pthat<80");
  MCentries[3] = ch->GetEntries("pthat<120");
  MCentries[4] = ch->GetEntries("pthat<170");
  MCentries[5] = ch->GetEntries("pthat<220");
  MCentries[6] = ch->GetEntries("pthat<280");
  if(QCDpthatBins==8)MCentries[7] = ch->GetEntries("pthat>=280");
  else {
   MCentries[7] = ch->GetEntries("pthat>=280 && pthat<370");
//   MCentries[7] = ch->GetEntries("pthat<370");
   MCentries[8] = ch->GetEntries("pthat>=370");
  }
*/  for(int i=0; i<QCDpthatBins; i++){
    cout << "QCD MCentries[" << i << "]: " << MCentries[i] << endl;
  }


  
  return MCentries;
}

//**********************************************************
// Trigger-Combine the data in order to unfold properly later
// need to adjust here to use the leading jet binning combination
//**********************************************************

//[0] = MB, [1] = Jet20, [2] = Jet40, [3] = Jet60, [4] = Jet80, [5] = Jet100
double trigComb(bool *trg, double *pscl, double pt, int combinationMethod){
  double weight=0;
//old method used by HIN-12-003
/*  //new method can be used for 5 triggers: [0] = Jet20, [4]= Jet100
    if(trg[0] && !trg[1] && !trg[2] && !trg[3] && !trg[4]) weight = 1./(1./pscl[0]);
    else if(trg[1] && !trg[2] && !trg[3] && !trg[4]) weight = 1.;
    else if(trg[2] && !trg[3] && !trg[4]) weight = 1.;
    else if(trg[3] && !trg[4]) weight = 1.;
    else weight =1. ;
  //  if(trg[4]) weight = 1.;

    if(trg[0] && !trg[1] && !trg[2] ) weight = 1./(1./pscl[0]);
    if(trg[1] && !trg[2]) weight = 1./(1./pscl[1]);
    else weight = 1. ; ///(1./pscl[2]) ;

  //HIN-12-017 (charged part RpA) combination method - solid but loses events that slip through the pt bins
    if(trg[0] && pt>40 && pt<60) weight = pscl[0];
    if(trg[1] && pt>60 && pt<75) weight = pscl[1];
    if(trg[2] && pt>75 && pt<95) weight = pscl[2];
    if(trg[3] && pt>95 && pt<120) weight = pscl[3];
    if(trg[4] && pt>120) weight = 1.;
*/
 if(combinationMethod==1){ 
   //Final CMS Methodology - jet-by-jet weighting method for CMS trigger-type setup.
        if(trg[4] && pt>=100) weight = 1; //pscl[4];
        if(trg[3] && pt>=80 && pt<100) weight = pscl[3];
        if(trg[2] && pt>=60 && pt<80) weight = pscl[2];
        if(trg[1] && pt>=40 && pt<60) weight = pscl[1];
        if(trg[0] && pt>=20 && pt<40) weight = pscl[0];
 }
 else {
 //HIN-12-017 (charged part RpA) combination method - solid but loses events that slip through the pt bins
    if(trg[0] && pt>40 && pt<=60) weight = pscl[0];
    if(trg[1] && pt>60 && pt<=75) weight = pscl[1];
    if(trg[2] && pt>75 && pt<=95) weight = pscl[2];
    if(trg[3] && pt>95 && pt<=120) weight = pscl[3];
    if(trg[4] && pt>120) weight = 1.;
 }   
  return weight;
}

//**********************************************************
// "get" the trigger prescales by counting trigger overlap
//**********************************************************

double* getPscls(std::string infile, int nFiles){
  
  TChain *dataCH = NULL;
//  if(usePUsub){
    dataCH = new TChain(Form("%sJetAnalyzer/t", algo.Data()));
//  }
//  else dataCH = new TChain("ak3PFJetAnalyzer/t");
  TChain *dataCH2 = new TChain("hltanalysis/HltTree");
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  for(int ifile=0; ifile<nFiles; ifile++){
    instr >> filename;
    dataCH->Add(filename.c_str());
    dataCH2->Add(filename.c_str());
  }
  dataCH->AddFriend(dataCH2, "hltanalysis/HltTree");
//    TFile *fin=NULL;
//    for(int ifile=0; ifile<nFiles; ifile++){
//        fin = TFile::Open(filename.c_str(), "readonly");
//        Long64_t nentries = dataCH->GetEntries();
//        for (Long64_t i=0; i<nentries;i++) {
//            dataCH->GetEntry(i);
//
//        }
//    }
//    bool jet40 = false ;
//    bool jet60 = false ;
//    bool jet75 = false ;
//    bool jet95 = false ;
//    bool jet120 = false ;
//    bool combined = false ;
    
  //Set up trigger combination prescales for data
  double ov1, ov2, ov3, ov4;
  ov1 = dataCH->GetEntries("HLT_PAJet20_NoJetID_v1 && HLT_PAJet80_NoJetID_v1");
  ov2 = dataCH->GetEntries("HLT_PAJet40_NoJetID_v1 && HLT_PAJet80_NoJetID_v1");
  ov3 = dataCH->GetEntries("HLT_PAJet60_NoJetID_v1 && HLT_PAJet80_NoJetID_v1");
  ov4 = dataCH->GetEntries("HLT_PAJet80_NoJetID_v1");
  double *pscls = new double[4];
  pscls[0] = ov4/ov1;
  pscls[1] = ov4/ov2;
  pscls[2] = ov4/ov3;
  pscls[3] = 1.;
  return pscls;
}

// *************************************************************
float triggerMatch(int trgObjSize, float* trigPhi, float* trigEta, float *trigPt, float jtphi, float jteta, float jtpt){
  float triggerPt=0;
  float closestMatch=999999.;
  
  for(int iObj=0; iObj<trgObjSize; iObj++){
    if(abs(trigPhi[iObj]-jtphi)<0.2 && abs(trigEta[iObj]-jteta)<0.2 && abs(trigPt[iObj]-jtpt)/jtpt<1. && abs(trigPt[iObj]-TMath::Floor(trigPt[iObj]))>0.0001){
      if(sqrt(pow((trigPhi[iObj]-jtphi),2)+pow((trigEta[iObj]-jteta),2))<closestMatch){
closestMatch = sqrt(pow((trigPhi[iObj]-jtphi),2)+pow((trigEta[iObj]-jteta),2));
triggerPt = trigPt[iObj];
      }
    }
  }
  return triggerPt;
}

// ****************************************************************


//**********************************************************
// ~~~ MAIN PROGRAM ~~~
//**********************************************************

void analyzeTreesPA(int isRecopp=1, int ppPbPb=0, int isMC=1, int doNtuples=1, int doJets=1, int doTracks=1, int updateJEC=0, int cbin=-1, bool ExpandedTree=false, bool useJetTrgAssociation=1)
{
  // isMC=0 --> Real data, ==1 --> QCD, ==2 --> bJet, ==3 --> cJet
  Float_t minJetPt=0.;
  
  Float_t maxJetEta=3.;
  Float_t minMuPt=5;
  TString coll = "PPb";
   if(coll=="PbP") updateJEC = 1 ;
  if(useJetTrgAssociation)  combinationMethod = 1 ;
  else combinationMethod = 0 ;
   if(isMC){
     startfile=0;
     endfile=QCDpthatBins ;
    }
   else {
      startfile=atoi(getenv("FIRST"));
    endfile=atoi(getenv("LAST"));
    }
   cout <<"analysis  "  << "  in file " << startfile << "  and " <<endfile<<endl ;

  timer.Start();  
  // cbin = -1 --> 0-100%
  // cbin = 0 --> 0-20%
  // cbin = 1 --> 20-50%
  // cbin =2 --> 50-100%
  if(!ppPbPb) cbin=-1;
  int useWeight=1;
//  double         triggerPt;
  trigO *HLT_PAJet_NoJetID_v1_trigObject[6];
    int trigObjSize[5];

  float trigObjPt[5][1000];
  float trigObjEta[5][1000];
  float trigObjPhi[5][1000];

  for(int i=0; i<6; i++){
    HLT_PAJet_NoJetID_v1_trigObject[i] = new trigO;
  } 

  //Get vz weight histograms for MC if needed
//  TH1D *hMCvz[2], hDatavz;
//  TFile *fMCvz, *fDatavz, *fBvz;
//  if(isMC){  
      TF1 * fCen = new TF1("fCen","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x+[5]*x*x*x*x+[6]*x*x*x*x*x)", 0., 100.);
      TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);a
if(coll=="PPb"){
      // fCen->SetParameters(8.68073e-03, 5.09356e+00, -1.33053e-02, 1.46904e-03, -6.99681e-05, 1.06721e-06, -5.21398e-09);
         fCen->SetParameters(1.20916e-02, 5.02157e+00, -3.38300e-02, 1.87647e-03, -6.76442e-05, 9.08602e-07, -4.01536e-09);//parameterize on 05.03 after approval
      fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);
	}
       else if(coll=="PbP"){
         fCen->SetParameters(1.05408e-02, 5.27477e+00, -8.03382e-02, 3.51669e-03, -8.85332e-05, 1.08917e-06, -4.90091e-09);
	 fVz->SetParameters(1.54398e+00, -8.56155e-03, -1.40026e-02, 4.01020e-05, 3.47683e-05); //latest parameterization
        //fVz->SetParameters(1.66731e+00,-2.43367e-03,-1.42488e-02,7.40147e-06,3.22477e-05);
}    
                    else{
                        fVz->SetParameters(1.,0,0,0,0);	}
//    fMCvz = new TFile("MCvzDistr.root");
//    hMCvz[0] = (TH1D*)(fMCvz->Get("hvz"))->Clone("hMCvz_0");
//    fDatavz = new TFile("DatavzDistro.root");
//    hDatavz = (TH1D*)(fMCvz->Get("hvzData"))->Clone("hDatavz");
//    fBvz = new TFile("BvzDistr.root");
//    hMCvz[1] = (TH1D*)(fBvz->Get("hvzB"))->Clone("hMCvz_1");
//    
//    hMCvz[0]->Scale(1./hMCvz[0]->Integral());
//    hMCvz[1]->Scale(1./hMCvz[1]->Integral());
//    hDatavz->Scale(1./hDatavz->Integral());
//  }

//      TString CorrfileName;
//    CorrfileName="/afs/cern.ch/work/y/ymao/public/RpA/Corrections/pPbHijing_EffCorr_VRw.root";
 //   Corrector corr_Qiao(CorrfileName);

      TFile *fcrel3 = NULL ;
     TH1D *C_rel= NULL ;
 
     fcrel3 = TFile::Open(Form("/afs/cern.ch/work/y/ymao/public/RpA/Corrections/Casym_pPb_double_hcalbins_algo_%s_pt100_140_jet80_alphahigh_20_phicut250.root", algo.Data()), "readonly");
     if(fcrel3)  C_rel=(TH1D*)fcrel3->Get("C_asym"); 
  cout << "Analyzing Trees! Assuming " << QCDpthatBins  << endl;

 // int pthatbin[10] = {15,30,50,80,120,170,220,280,370, 9999};
//  double wght[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 1.088E-07, 0};
/*  if(QCDpthatBins==8) {
    pthatbin[]= {15,30,50,80,120,170,220,280, 9999};
    wght[]={5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 0};
  }
 else {
   pthatbin[]= {15,30,50,80,120,170,220,280,370, 9999};
    wght[]={5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 1.088E-07, 0};
  } 
*/ //  int pthatbin[QCDpthatBins+1] = {80, 9999};
//  double wght[QCDpthatBins+1]={4.412E-04, 0};

  TFile *fin=NULL;
  std::string infile;
  int *MCentr = NULL;
  double *pscls = NULL;
 
  //File load
    if(!isMC){
  //    infile = "pPbDataForest.txt";
  //   infile = "DataForestTrig60to80to100.txt";
     infile = "pAFileListKurtv4.txt";
    }
    else {
    //  infile = "pPbMCKurtForest.txt";
   //   infile = "PbpForestList.txt";
    //  infile = "pPbForestListProd25.txt";
     if(coll=="PPb") infile = "pPbMCKurtForest.txt";
        else infile = "ppMCKurtForest.txt";
    //  infile = "test.txt";
    }

  if(!isMC){
   // pscls = getPscls(infile,QCDpthatBins,usePUsub);
  }
  
  int dupRuns[6] = {181912,181913,181938,181950,181985,182124};
  
  std::vector<int> usedEvents[6];
  int nDup=0;

  //Declaration of leaves types
  Int_t evt;
  Int_t run;
  Int_t bin;
  Int_t lumi;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
  Float_t hf;
  Float_t vz;
  Int_t nref;
  Float_t rawpt[1000];
  Float_t jtpt[1000];
  Float_t jteta[1000];
  Float_t jty[1000];
  Float_t jtphi[1000];
	Int_t chargedN[1000];
Int_t photonN[1000];
Int_t neutralN[1000];
  Float_t jtpu[1000];
  Float_t discr_ssvHighEff[1000];
  Float_t discr_ssvHighPur[1000];
  //Float_t discr_csvMva[1000];
  Float_t discr_csvSimple[1000];
  // Float_t discr_muByIp3[1000];
  Float_t discr_muByPt[1000];
  Float_t discr_prob[1000];
  Float_t discr_probb[1000];
  Float_t discr_tcHighEff[1000];
  Float_t discr_tcHighPur[1000];
  Int_t nsvtx[1000];
  Int_t svtxntrk[1000];
  Float_t svtxdl[1000];
  Float_t svtxdls[1000];
  Float_t svtxm[1000];
  Float_t svtxpt[1000];
  Int_t nIPtrk[1000];
  Int_t nselIPtrk[1000];
  Int_t nIP;
  Int_t nTrk;
  Float_t trkPt[1000];
  Float_t trkPhi[1000];
  Float_t trkEta[1000];
  Int_t ipJetIndex[10000];
  Float_t ipPt[10000];
  Float_t ipProb0[10000];
  Float_t ipProb1[10000];
  Float_t ip2d[10000];
  Float_t ip2dSig[10000];
  Float_t ip3d[10000];
  Float_t ip3dSig[10000];
  Float_t ipDist2Jet[10000];
  Float_t ipDist2JetSig[10000];
  Float_t ipClosest2Jet[10000];
  //Float_t mue[1000];
  Float_t mupt[1000];
  Float_t muptPF[1000];
  Float_t mueta[1000];
  Float_t muphi[1000];
  //Float_t mudr[1000];
  Float_t muptrel[1000];
  //Int_t muchg[1000];
  Float_t pthat;
  Int_t beamId1;
  Int_t beamId2;
  Int_t subid[1000];
  Float_t refpt[1000];
  Float_t refeta[1000];
  Float_t refy[1000];
  Float_t refphi[1000];
  Float_t refdphijt[1000];
  Float_t refdrjt[1000];
  Float_t refparton_pt[1000];
  Int_t refparton_flavor[1000];
  Int_t refparton_flavorForB[1000];
  Bool_t refparton_isGSP[1000];
  Int_t          HLT_PAZeroBiasPixel_SingleTrack_v1;
  Int_t          HLT_PAJet20_NoJetID_v1;
  Int_t          HLT_PAJet40_NoJetID_v1;
  Int_t          HLT_PAJet60_NoJetID_v1;
  Int_t          HLT_PAJet80_NoJetID_v1;
  Int_t          HLT_PAJet100_NoJetID_v1;
  Int_t          HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl;
  Int_t          HLT_PAJet20_NoJetID_v1_Prescl;
  Int_t          HLT_PAJet40_NoJetID_v1_Prescl;
  Int_t          HLT_PAJet60_NoJetID_v1_Prescl;
  Int_t          HLT_PAJet80_NoJetID_v1_Prescl;
  Int_t          HLT_PAJet100_NoJetID_v1_Prescl;
  Int_t pVertexFilterCutGplusUpsPP;
  Int_t pVertexFilterCutGplus;
  Int_t pPAcollisionEventSelectionPA;
  Int_t pBeamScrapingFilter;
  Int_t pprimaryvertexFilter;
  Int_t phfPosFilter1;
  Int_t phfNegFilter1;
  Int_t pHBHENoiseFilter;

  /*
Int_t ngen;
Int_t genmatchindex[1000];
Float_t genpt[1000];
Float_t geneta[1000];
Float_t geny[1000];
Float_t genphi[1000];
Float_t gendphijt[1000];
Float_t gendrjt[1000];
*/

  float chargedMax[1000];
  float photonMax[1000];
  float neutralMax[1000];
  float chargedSum[1000];
  float photonSum[1000];
  float neutralSum[1000];
  float muSum[1000];
  float eSum[1000];

  Double_t weight, xSecWeight, centWeight, vzWeight;
 
  centWeight = 1. ; 
  int nHLTBit;
  bool hltBit[12];

  Int_t pvSel;
  Int_t hbheNoiseSel;
  Int_t spikeSel;
  Int_t collSell;

  TFile *fout=NULL;

   if ( isRecopp && useJetTrgAssociation ) { // pp reco, jet triggered
      if(isMC) fout=new TFile(Form("/afs/cern.ch/work/q/qixu/private/RpA/skimTree/%sMCKurtForest_ppReco_%s_QCDjetTrig_JetPt%.fnoIPupperCut.root",coll.Data(), algo.Data(), minJetPt),"recreate");
      else fout=new TFile(Form("/afs/cern.ch/work/q/qixu/private/RpA/skimTree/%sdata_ppReco_%s_AlljetTrigKurtTrCombFile%d_%d_JetPt%.fnoIPupperCut.root", coll.Data(),algo.Data(),startfile,endfile, minJetPt),"recreate");
    }
   else if ( isRecopp && !useJetTrgAssociation ) { // pp reco, jet triggered
      if(isMC) fout=new TFile(Form("/afs/cern.ch/work/q/qixu/private/RpA/skimTree/%sMCKurtForest_ppReco_%s_QCDjetTrig_JetPt%.fnoIPupperCut.root",coll.Data(), algo.Data(), minJetPt),"recreate");
      else fout=new TFile(Form("/afs/cern.ch/work/q/qixu/private/RpA/skimTree/%sdata_ppReco_%s_AlljetTrigEricTrCombFile%d_%d_JetPt%.fnoIPupperCut.root", coll.Data(), algo.Data(),startfile,endfile, minJetPt),"recreate");
    }
    else { // hi reco, jet triggered
      if(isMC)fout=new TFile("/afs/cern.ch/work/q/qixu/private/RpA/skimTree/ppMC_hiReco_jetTrig_addGSP_up.root","recreate");
      else fout=new TFile("/afs/cern.ch/work/q/qixu/private/RpA/skimTree/ppdata_hiReco_jetTrig_regPFforJets.root","recreate");
    }


  TH1D *hbin = new TH1D("hbin","hbin",100,0.5,100.5);
  TH1D *hbinw = new TH1D("hbinw","hbinw",100,0.5,100.5);
  hbin->Sumw2(); hbinw->Sumw2();

  TH1D *hvz = new TH1D("hvz","hvz",120,-15.,15.);
  TH1D *hvzw = new TH1D("hvzw","hvzw",120,-15.,15.);
  hvz->Sumw2(); hvzw->Sumw2();
  TH1F * hJetPtEtaBin[netabin];
  for(Int_t ieta = 0 ; ieta < netabin ; ieta++){
     hJetPtEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), nJetPtBin, jetPtBin);
     hJetPtEtaBin[ieta]->Sumw2();
   }

//  TH1D *hjtpt = new TH1D("hjtpt","hjtpt",1000,0,1000);
  TH1D *hjtpt = new TH1D("hjtpt","hjtpt", nJetPtBin, jetPtBin);
  TH1D *hrefpt = new TH1D("hrefpt","hrefpt", nJetPtBin, jetPtBin);
  TH1F *hpthat = new TH1F("hpthat","hpthat", 400, 0., 400.);
  TH1F *hpthatw = new TH1F("hpthatw","hpthatw", 400, 0., 400.);
  TH2F *hpthatcs = new TH2F("pthatweight","pthatweight", 400, 0., 400., 200, -15., -5.);
   hrefpt->Sumw2(); hpthat->Sumw2();hpthatw->Sumw2();hpthatcs->Sumw2();

  TH1D *hjtptB = new TH1D("hjtptB","hjtptB",68,80,330);
  TH1D *hjtptC = new TH1D("hjtptC","hjtptC",68,80,330);
  TH1D *hjtptL = new TH1D("hjtptL","hjtptL",68,80,330);
  TH1D *hjtptU = new TH1D("hjtptU","hjtptU",68,80,330);
  hjtpt->Sumw2(); hjtptB->Sumw2(); hjtptC->Sumw2(); hjtptL->Sumw2(); hjtptU->Sumw2();

  TH1D *hrawpt = new TH1D("hrawpt","hrawpt",500,0,500);
  TH1D *hrawptB = new TH1D("hrawptB","hrawptB",68,80,330);
  TH1D *hrawptC = new TH1D("hrawptC","hrawptC",68,80,330);
  TH1D *hrawptL = new TH1D("hrawptL","hrawptL",68,80,330);
  hrawpt->Sumw2(); hrawptB->Sumw2(); hrawptC->Sumw2(); hrawptL->Sumw2();

  TH1F *hjteta = new TH1F("hjteta","hjteta",100,-5,5);
  TH1D *hjtetaB = new TH1D("hjtetaB","hjtetaB",40,-2,2);
  TH1D *hjtetaC = new TH1D("hjtetaC","hjtetaC",40,-2,2);
  TH1D *hjtetaL = new TH1D("hjtetaL","hjtetaL",40,-2,2);
    TH1F *hrefeta = new TH1F("hrefeta","hrefeta",100,-5,5);
  TProfile *hetares = new TProfile("hetares","hetares",100,-5,5);
  hjteta->Sumw2(); hjtetaB->Sumw2(); hjtetaC->Sumw2(); hjtetaL->Sumw2();
  hrefeta->Sumw2();  hetares->Sumw2();

  TH1D *hjtphi = new TH1D("hjtphi","hjtphi",40,-1.*acos(-1.),acos(-1.));
  TH1D *hjtphiB = new TH1D("hjtphiB","hjtphiB",40,-1.*acos(-1.),acos(-1.));
  TH1D *hjtphiC = new TH1D("hjtphiC","hjtphiC",40,-1.*acos(-1.),acos(-1.));
  TH1D *hjtphiL = new TH1D("hjtphiL","hjtphiL",40,-1.*acos(-1.),acos(-1.));
  hjtphi->Sumw2(); hjtphiB->Sumw2(); hjtphiC->Sumw2(); hjtphiL->Sumw2();

  TH1D *hdiscr_csvSimple = new TH1D("hdiscr_csvSimple","hdiscr_csvSimple",25,0,1);
  TH1D *hdiscr_csvSimpleB = new TH1D("hdiscr_csvSimpleB","hdiscr_csvSimpleB",25,0,1);
  TH1D *hdiscr_csvSimpleC = new TH1D("hdiscr_csvSimpleC","hdiscr_csvSimpleC",25,0,1);
  TH1D *hdiscr_csvSimpleL = new TH1D("hdiscr_csvSimpleL","hdiscr_csvSimpleL",25,0,1);
  hdiscr_csvSimple->Sumw2(); hdiscr_csvSimpleB->Sumw2(); hdiscr_csvSimpleC->Sumw2(); hdiscr_csvSimpleL->Sumw2();
  
  TH1D *hdiscr_prob = new TH1D("hdiscr_prob","hdiscr_prob",25,0,2.5);
  TH1D *hdiscr_probB = new TH1D("hdiscr_probB","hdiscr_probB",25,0,2.5);
  TH1D *hdiscr_probC = new TH1D("hdiscr_probC","hdiscr_probC",25,0,2.5);
  TH1D *hdiscr_probL = new TH1D("hdiscr_probL","hdiscr_probL",25,0,2.5);
  hdiscr_prob->Sumw2(); hdiscr_probB->Sumw2(); hdiscr_probC->Sumw2(); hdiscr_probL->Sumw2();
  
  TH1D *hdiscr_ssvHighEff = new TH1D("hdiscr_ssvHighEff","hdiscr_ssvHighEff",25,1,6);
  TH1D *hdiscr_ssvHighEffB = new TH1D("hdiscr_ssvHighEffB","hdiscr_ssvHighEffB",25,1,6);
  TH1D *hdiscr_ssvHighEffC = new TH1D("hdiscr_ssvHighEffC","hdiscr_ssvHighEffC",25,1,6);
  TH1D *hdiscr_ssvHighEffL = new TH1D("hdiscr_ssvHighEffL","hdiscr_ssvHighEffL",25,1,6);
  hdiscr_ssvHighEff->Sumw2(); hdiscr_ssvHighEffB->Sumw2(); hdiscr_ssvHighEffC->Sumw2(); hdiscr_ssvHighEffL->Sumw2();
  
  TH1D *hdiscr_ssvHighPur = new TH1D("hdiscr_ssvHighPur","hdiscr_ssvHighPur",25,1,6);
  TH1D *hdiscr_ssvHighPurB = new TH1D("hdiscr_ssvHighPurB","hdiscr_ssvHighPurB",25,1,6);
  TH1D *hdiscr_ssvHighPurC = new TH1D("hdiscr_ssvHighPurC","hdiscr_ssvHighPurC",25,1,6);
  TH1D *hdiscr_ssvHighPurL = new TH1D("hdiscr_ssvHighPurL","hdiscr_ssvHighPurL",25,1,6);
  hdiscr_ssvHighPur->Sumw2(); hdiscr_ssvHighPurB->Sumw2(); hdiscr_ssvHighPurC->Sumw2(); hdiscr_ssvHighPurL->Sumw2();

  TH1D *hdiscr_tcHighEff = new TH1D("hdiscr_tcHighEff","hdiscr_tcHighEff",25,1,6);
  TH1D *hdiscr_tcHighEffB = new TH1D("hdiscr_tcHighEffB","hdiscr_tcHighEffB",25,1,6);
  TH1D *hdiscr_tcHighEffC = new TH1D("hdiscr_tcHighEffC","hdiscr_tcHighEffC",25,1,6);
  TH1D *hdiscr_tcHighEffL = new TH1D("hdiscr_tcHighEffL","hdiscr_tcHighEffL",25,1,6);
  hdiscr_tcHighEff->Sumw2(); hdiscr_tcHighEffB->Sumw2(); hdiscr_tcHighEffC->Sumw2(); hdiscr_tcHighEffL->Sumw2();
  
  TH1D *hdiscr_tcHighPur = new TH1D("hdiscr_tcHighPur","hdiscr_tcHighPur",25,1,6);
  TH1D *hdiscr_tcHighPurB = new TH1D("hdiscr_tcHighPurB","hdiscr_tcHighPurB",25,1,6);
  TH1D *hdiscr_tcHighPurC = new TH1D("hdiscr_tcHighPurC","hdiscr_tcHighPurC",25,1,6);
  TH1D *hdiscr_tcHighPurL = new TH1D("hdiscr_tcHighPurL","hdiscr_tcHighPurL",25,1,6);
  hdiscr_tcHighPur->Sumw2(); hdiscr_tcHighPurB->Sumw2(); hdiscr_tcHighPurC->Sumw2(); hdiscr_tcHighPurL->Sumw2();
  
  TH1D *hnsvtx = new TH1D("hnsvtx","hnsvtx",6,-0.5,5.5);
  TH1D *hnsvtxB = new TH1D("hnsvtxB","hnsvtxB",6,-0.5,5.5);
  TH1D *hnsvtxC = new TH1D("hnsvtxC","hnsvtxC",6,-0.5,5.5);
  TH1D *hnsvtxL = new TH1D("hnsvtxL","hnsvtxL",6,-0.5,5.5);
  hnsvtx->Sumw2(); hnsvtxB->Sumw2(); hnsvtxC->Sumw2(); hnsvtxL->Sumw2();
  
  TH1D *hsvtxntrk = new TH1D("hsvtxntrk","hsvtxntrk",12,-0.5,11.5);
  TH1D *hsvtxntrkB = new TH1D("hsvtxntrkB","hsvtxntrkB",12,-0.5,11.5);
  TH1D *hsvtxntrkC = new TH1D("hsvtxntrkC","hsvtxntrkC",12,-0.5,11.5);
  TH1D *hsvtxntrkL = new TH1D("hsvtxntrkL","hsvtxntrkL",12,-0.5,11.5);
  hsvtxntrk->Sumw2(); hsvtxntrkB->Sumw2(); hsvtxntrkC->Sumw2(); hsvtxntrkL->Sumw2();

  TH1D *hsvtxdl = new TH1D("hsvtxdl","hsvtxdl",20,0,10);
  TH1D *hsvtxdlB = new TH1D("hsvtxdlB","hsvtxdlB",20,0,10);
  TH1D *hsvtxdlC = new TH1D("hsvtxdlC","hsvtxdlC",20,0,10);
  TH1D *hsvtxdlL = new TH1D("hsvtxdlL","hsvtxdlL",20,0,10);
  hsvtxdl->Sumw2(); hsvtxdlB->Sumw2(); hsvtxdlC->Sumw2(); hsvtxdlL->Sumw2();

  TH1D *hsvtxdls = new TH1D("hsvtxdls","hsvtxdls",40,0,80);
  TH1D *hsvtxdlsB = new TH1D("hsvtxdlsB","hsvtxdlsB",40,0,80);
  TH1D *hsvtxdlsC = new TH1D("hsvtxdlsC","hsvtxdlsC",40,0,80);
  TH1D *hsvtxdlsL = new TH1D("hsvtxdlsL","hsvtxdlsL",40,0,80);
  hsvtxdls->Sumw2(); hsvtxdlsB->Sumw2(); hsvtxdlsC->Sumw2(); hsvtxdlsL->Sumw2();
  
  TH1D *hsvtxm = new TH1D("hsvtxm","hsvtxm",32,0,8);
  TH1D *hsvtxmB = new TH1D("hsvtxmB","hsvtxmB",32,0,8);
  TH1D *hsvtxmC = new TH1D("hsvtxmC","hsvtxmC",32,0,8);
  TH1D *hsvtxmL = new TH1D("hsvtxmL","hsvtxmL",32,0,8);
  hsvtxm->Sumw2(); hsvtxmB->Sumw2(); hsvtxmC->Sumw2(); hsvtxmL->Sumw2();
  
  TH1D *hsvtxmSV3 = new TH1D("hsvtxmSV3","hsvtxmSV3",32,0,8);
  TH1D *hsvtxmSV3B = new TH1D("hsvtxmSV3B","hsvtxmSV3B",32,0,8);
  TH1D *hsvtxmSV3C = new TH1D("hsvtxmSV3C","hsvtxmSV3C",32,0,8);
  TH1D *hsvtxmSV3L = new TH1D("hsvtxmSV3L","hsvtxmSV3L",32,0,8);
  hsvtxmSV3->Sumw2(); hsvtxmSV3B->Sumw2(); hsvtxmSV3C->Sumw2(); hsvtxmSV3L->Sumw2();
  
  TH1D *hsvtxpt = new TH1D("hsvtxpt","hsvtxpt",20,0,100);
  TH1D *hsvtxptB = new TH1D("hsvtxptB","hsvtxptB",20,0,100);
  TH1D *hsvtxptC = new TH1D("hsvtxptC","hsvtxptC",20,0,100);
  TH1D *hsvtxptL = new TH1D("hsvtxptL","hsvtxptL",20,0,100);
  hsvtxpt->Sumw2(); hsvtxptB->Sumw2(); hsvtxptC->Sumw2(); hsvtxptL->Sumw2();
  
  TH1D *hsvtxptSV3 = new TH1D("hsvtxptSV3","hsvtxptSV3",20,0,100);
  TH1D *hsvtxptSV3B = new TH1D("hsvtxptSV3B","hsvtxptSV3B",20,0,100);
  TH1D *hsvtxptSV3C = new TH1D("hsvtxptSV3C","hsvtxptSV3C",20,0,100);
  TH1D *hsvtxptSV3L = new TH1D("hsvtxptSV3L","hsvtxptSV3L",20,0,100);
  hsvtxptSV3->Sumw2(); hsvtxptSV3B->Sumw2(); hsvtxptSV3C->Sumw2(); hsvtxptSV3L->Sumw2();
  
  TH1D *hnIPtrk = new TH1D("hnIPtrk","hnIPtrk",100,0,100);
  TH1D *hnIPtrkB = new TH1D("hnIPtrkB","hnIPtrkB",100,0,100);
  TH1D *hnIPtrkC = new TH1D("hnIPtrkC","hnIPtrkC",100,0,100);
  TH1D *hnIPtrkL = new TH1D("hnIPtrkL","hnIPtrkL",100,0,100);
  hnIPtrk->Sumw2(); hnIPtrkB->Sumw2(); hnIPtrkC->Sumw2(); hnIPtrkL->Sumw2();
  
  TH1D *hnselIPtrk = new TH1D("hnselIPtrk","hnselIPtrk",100,0,100);
  TH1D *hnselIPtrkB = new TH1D("hnselIPtrkB","hnselIPtrkB",100,0,100);
  TH1D *hnselIPtrkC = new TH1D("hnselIPtrkC","hnselIPtrkC",100,0,100);
  TH1D *hnselIPtrkL = new TH1D("hnselIPtrkL","hnselIPtrkL",100,0,100);
  hnselIPtrk->Sumw2(); hnselIPtrkB->Sumw2(); hnselIPtrkC->Sumw2(); hnselIPtrkL->Sumw2();
  
  TH1D *hmuptrel = new TH1D("hmuptrel","hmuptrel",40,0,4);
  TH1D *hmuptrelB = new TH1D("hmuptrelB","hmuptrelB",40,0,4);
  TH1D *hmuptrelC = new TH1D("hmuptrelC","hmuptrelC",40,0,4);
  TH1D *hmuptrelL = new TH1D("hmuptrelL","hmuptrelL",40,0,4);
  hmuptrel->Sumw2(); hmuptrelB->Sumw2(); hmuptrelC->Sumw2(); hmuptrelL->Sumw2();
  
  TH1D *hmuptrelSV2 = new TH1D("hmuptrelSV2","hmuptrelSV2",40,0,4);
  TH1D *hmuptrelSV2B = new TH1D("hmuptrelSV2B","hmuptrelSV2B",40,0,4);
  TH1D *hmuptrelSV2C = new TH1D("hmuptrelSV2C","hmuptrelSV2C",40,0,4);
  TH1D *hmuptrelSV2L = new TH1D("hmuptrelSV2L","hmuptrelSV2L",40,0,4);
  hmuptrelSV2->Sumw2(); hmuptrelSV2B->Sumw2(); hmuptrelSV2C->Sumw2(); hmuptrelSV2L->Sumw2();
  
  TH1D *hmuptrelSV3 = new TH1D("hmuptrelSV3","hmuptrelSV3",40,0,4);
  TH1D *hmuptrelSV3B = new TH1D("hmuptrelSV3B","hmuptrelSV3B",40,0,4);
  TH1D *hmuptrelSV3C = new TH1D("hmuptrelSV3C","hmuptrelSV3C",40,0,4);
  TH1D *hmuptrelSV3L = new TH1D("hmuptrelSV3L","hmuptrelSV3L",40,0,4);
  hmuptrelSV3->Sumw2(); hmuptrelSV3B->Sumw2(); hmuptrelSV3C->Sumw2(); hmuptrelSV3L->Sumw2();
  
  TH1D *hipPt = new TH1D("hipPt","hipPt",40,0,40);
  TH1D *hipPtB = new TH1D("hipPtB","hipPtB",40,0,40);
  TH1D *hipPtC = new TH1D("hipPtC","hipPtC",40,0,40);
  TH1D *hipPtL = new TH1D("hipPtL","hipPtL",40,0,40);
  hipPt->Sumw2(); hipPtB->Sumw2(); hipPtC->Sumw2(); hipPtL->Sumw2();
  
  TH1D *hipProb0 = new TH1D("hipProb0","hipProb0",40,-1,1);
  TH1D *hipProb0B = new TH1D("hipProb0B","hipProb0B",40,-1,1);
  TH1D *hipProb0C = new TH1D("hipProb0C","hipProb0C",40,-1,1);
  TH1D *hipProb0L = new TH1D("hipProb0L","hipProb0L",40,-1,1);
  hipProb0->Sumw2(); hipProb0B->Sumw2(); hipProb0C->Sumw2(); hipProb0L->Sumw2();
  
  TH1D *hipProb1 = new TH1D("hipProb1","hipProb1",40,-1,1);
  TH1D *hipProb1B = new TH1D("hipProb1B","hipProb1B",40,-1,1);
  TH1D *hipProb1C = new TH1D("hipProb1C","hipProb1C",40,-1,1);
  TH1D *hipProb1L = new TH1D("hipProb1L","hipProb1L",40,-1,1);
  hipProb1->Sumw2(); hipProb1B->Sumw2(); hipProb1C->Sumw2(); hipProb1L->Sumw2();
  
  TH1D *hip2d = new TH1D("hip2d","hip2d",40,-0.1,0.1);
  TH1D *hip2dB = new TH1D("hip2dB","hip2dB",40,-0.1,0.1);
  TH1D *hip2dC = new TH1D("hip2dC","hip2dC",40,-0.1,0.1);
  TH1D *hip2dL = new TH1D("hip2dL","hip2dL",40,-0.1,0.1);
  hip2d->Sumw2(); hip2dB->Sumw2(); hip2dC->Sumw2(); hip2dL->Sumw2();
  
  TH1D *hip2dSig = new TH1D("hip2dSig","hip2dSig",70,-35,35);
  TH1D *hip2dSigB = new TH1D("hip2dSigB","hip2dSigB",70,-35,35);
  TH1D *hip2dSigC = new TH1D("hip2dSigC","hip2dSigC",70,-35,35);
  TH1D *hip2dSigL = new TH1D("hip2dSigL","hip2dSigL",70,-35,35);
  hip2dSig->Sumw2(); hip2dSigB->Sumw2(); hip2dSigC->Sumw2(); hip2dSigL->Sumw2();

  TH1D *hip2d1 = new TH1D("hip2d1","hip2d1",40,-0.1,0.1);
  TH1D *hip2d1B = new TH1D("hip2d1B","hip2d1B",40,-0.1,0.1);
  TH1D *hip2d1C = new TH1D("hip2d1C","hip2d1C",40,-0.1,0.1);
  TH1D *hip2d1L = new TH1D("hip2d1L","hip2d1L",40,-0.1,0.1);
  hip2d1->Sumw2(); hip2d1B->Sumw2(); hip2d1C->Sumw2(); hip2d1L->Sumw2();

  TH1D *hip2dSig1 = new TH1D("hip2dSig1","hip2dSig1",70,-35,35);
  TH1D *hip2dSig1B = new TH1D("hip2dSig1B","hip2dSig1B",70,-35,35);
  TH1D *hip2dSig1C = new TH1D("hip2dSig1C","hip2dSig1C",70,-35,35);
  TH1D *hip2dSig1L = new TH1D("hip2dSig1L","hip2dSig1L",70,-35,35);
  hip2dSig1->Sumw2(); hip2dSig1B->Sumw2(); hip2dSig1C->Sumw2(); hip2dSig1L->Sumw2();

  TH1D *hip2d2 = new TH1D("hip2d2","hip2d2",40,-0.1,0.1);
  TH1D *hip2d2B = new TH1D("hip2d2B","hip2d2B",40,-0.1,0.1);
  TH1D *hip2d2C = new TH1D("hip2d2C","hip2d2C",40,-0.1,0.1);
  TH1D *hip2d2L = new TH1D("hip2d2L","hip2d2L",40,-0.1,0.1);
  hip2d2->Sumw2(); hip2d2B->Sumw2(); hip2d2C->Sumw2(); hip2d2L->Sumw2();

  TH1D *hip2dSig2 = new TH1D("hip2dSig2","hip2dSig2",70,-35,35);
  TH1D *hip2dSig2B = new TH1D("hip2dSig2B","hip2dSig2B",70,-35,35);
  TH1D *hip2dSig2C = new TH1D("hip2dSig2C","hip2dSig2C",70,-35,35);
  TH1D *hip2dSig2L = new TH1D("hip2dSig2L","hip2dSig2L",70,-35,35);
  hip2dSig2->Sumw2(); hip2dSig2B->Sumw2(); hip2dSig2C->Sumw2(); hip2dSig2L->Sumw2();

  TH1D *hip2d3 = new TH1D("hip2d3","hip2d3",40,-0.1,0.1);
  TH1D *hip2d3B = new TH1D("hip2d3B","hip2d3B",40,-0.1,0.1);
  TH1D *hip2d3C = new TH1D("hip2d3C","hip2d3C",40,-0.1,0.1);
  TH1D *hip2d3L = new TH1D("hip2d3L","hip2d3L",40,-0.1,0.1);
  hip2d3->Sumw2(); hip2d3B->Sumw2(); hip2d3C->Sumw2(); hip2d3L->Sumw2();

  TH1D *hip2dSig3 = new TH1D("hip2dSig3","hip2dSig3",70,-35,35);
  TH1D *hip2dSig3B = new TH1D("hip2dSig3B","hip2dSig3B",70,-35,35);
  TH1D *hip2dSig3C = new TH1D("hip2dSig3C","hip2dSig3C",70,-35,35);
  TH1D *hip2dSig3L = new TH1D("hip2dSig3L","hip2dSig3L",70,-35,35);
  hip2dSig3->Sumw2(); hip2dSig3B->Sumw2(); hip2dSig3C->Sumw2(); hip2dSig3L->Sumw2();
  
  TH1D *hip3d = new TH1D("hip3d","hip3d",40,-0.1,0.1);
  TH1D *hip3dB = new TH1D("hip3dB","hip3dB",40,-0.1,0.1);
  TH1D *hip3dC = new TH1D("hip3dC","hip3dC",40,-0.1,0.1);
  TH1D *hip3dL = new TH1D("hip3dL","hip3dL",40,-0.1,0.1);
  hip3d->Sumw2(); hip3dB->Sumw2(); hip3dC->Sumw2(); hip3dL->Sumw2();
  
  TH1D *hip3dSig = new TH1D("hip3dSig","hip3dSig",70,-35,35);
  TH1D *hip3dSigB = new TH1D("hip3dSigB","hip3dSigB",70,-35,35);
  TH1D *hip3dSigC = new TH1D("hip3dSigC","hip3dSigC",70,-35,35);
  TH1D *hip3dSigL = new TH1D("hip3dSigL","hip3dSigL",70,-35,35);
  hip3dSig->Sumw2(); hip3dSigB->Sumw2(); hip3dSigC->Sumw2(); hip3dSigL->Sumw2();
  
  TH1D *hip3d1 = new TH1D("hip3d1","hip3d1",40,-0.1,0.1);
  TH1D *hip3d1B = new TH1D("hip3d1B","hip3d1B",40,-0.1,0.1);
  TH1D *hip3d1C = new TH1D("hip3d1C","hip3d1C",40,-0.1,0.1);
  TH1D *hip3d1L = new TH1D("hip3d1L","hip3d1L",40,-0.1,0.1);
  hip3d1->Sumw2(); hip3d1B->Sumw2(); hip3d1C->Sumw2(); hip3d1L->Sumw2();

  TH1D *hip3dSig1 = new TH1D("hip3dSig1","hip3dSig1",70,-35,35);
  TH1D *hip3dSig1B = new TH1D("hip3dSig1B","hip3dSig1B",70,-35,35);
  TH1D *hip3dSig1C = new TH1D("hip3dSig1C","hip3dSig1C",70,-35,35);
  TH1D *hip3dSig1L = new TH1D("hip3dSig1L","hip3dSig1L",70,-35,35);
  hip3dSig1->Sumw2(); hip3dSig1B->Sumw2(); hip3dSig1C->Sumw2(); hip3dSig1L->Sumw2();

  TH1D *hip3d2 = new TH1D("hip3d2","hip3d2",40,-0.1,0.1);
  TH1D *hip3d2B = new TH1D("hip3d2B","hip3d2B",40,-0.1,0.1);
  TH1D *hip3d2C = new TH1D("hip3d2C","hip3d2C",40,-0.1,0.1);
  TH1D *hip3d2L = new TH1D("hip3d2L","hip3d2L",40,-0.1,0.1);
  hip3d2->Sumw2(); hip3d2B->Sumw2(); hip3d2C->Sumw2(); hip3d2L->Sumw2();

  TH1D *hip3dSig2 = new TH1D("hip3dSig2","hip3dSig2",70,-35,35);
  TH1D *hip3dSig2B = new TH1D("hip3dSig2B","hip3dSig2B",70,-35,35);
  TH1D *hip3dSig2C = new TH1D("hip3dSig2C","hip3dSig2C",70,-35,35);
  TH1D *hip3dSig2L = new TH1D("hip3dSig2L","hip3dSig2L",70,-35,35);
  hip3dSig2->Sumw2(); hip3dSig2B->Sumw2(); hip3dSig2C->Sumw2(); hip3dSig2L->Sumw2();

  TH1D *hip3d3 = new TH1D("hip3d3","hip3d3",40,-0.1,0.1);
  TH1D *hip3d3B = new TH1D("hip3d3B","hip3d3B",40,-0.1,0.1);
  TH1D *hip3d3C = new TH1D("hip3d3C","hip3d3C",40,-0.1,0.1);
  TH1D *hip3d3L = new TH1D("hip3d3L","hip3d3L",40,-0.1,0.1);
  hip3d3->Sumw2(); hip3d3B->Sumw2(); hip3d3C->Sumw2(); hip3d3L->Sumw2();

  TH1D *hip3dSig3 = new TH1D("hip3dSig3","hip3dSig3",70,-35,35);
  TH1D *hip3dSig3B = new TH1D("hip3dSig3B","hip3dSig3B",70,-35,35);
  TH1D *hip3dSig3C = new TH1D("hip3dSig3C","hip3dSig3C",70,-35,35);
  TH1D *hip3dSig3L = new TH1D("hip3dSig3L","hip3dSig3L",70,-35,35);
  hip3dSig3->Sumw2(); hip3dSig3B->Sumw2(); hip3dSig3C->Sumw2(); hip3dSig3L->Sumw2();

  TH1D *hipDist2Jet = new TH1D("hipDist2Jet","hipDist2Jet",40,-0.1,0);
  TH1D *hipDist2JetB = new TH1D("hipDist2JetB","hipDist2JetB",40,-0.1,0);
  TH1D *hipDist2JetC = new TH1D("hipDist2JetC","hipDist2JetC",40,-0.1,0);
  TH1D *hipDist2JetL = new TH1D("hipDist2JetL","hipDist2JetL",40,-0.1,0);
  hipDist2Jet->Sumw2(); hipDist2JetB->Sumw2(); hipDist2JetC->Sumw2(); hipDist2JetL->Sumw2();
  
  TH1D *hipDist2JetSig = new TH1D("hipDist2JetSig","hipDist2JetSig",40,-0.1,0.1);
  TH1D *hipDist2JetSigB = new TH1D("hipDist2JetSigB","hipDist2JetSigB",40,-0.1,0.1);
  TH1D *hipDist2JetSigC = new TH1D("hipDist2JetSigC","hipDist2JetSigC",40,-0.1,0.1);
  TH1D *hipDist2JetSigL = new TH1D("hipDist2JetSigL","hipDist2JetSigL",40,-0.1,0.1);
  hipDist2JetSig->Sumw2(); hipDist2JetSigB->Sumw2(); hipDist2JetSigC->Sumw2(); hipDist2JetSigL->Sumw2();
  
  TH1D *hipClosest2Jet = new TH1D("hipClosest2Jet","hipClosest2Jet",40,0,1);
  TH1D *hipClosest2JetB = new TH1D("hipClosest2JetB","hipClosest2JetB",40,0,1);
  TH1D *hipClosest2JetC = new TH1D("hipClosest2JetC","hipClosest2JetC",40,0,1);
  TH1D *hipClosest2JetL = new TH1D("hipClosest2JetL","hipClosest2JetL",40,0,1);
  hipClosest2Jet->Sumw2(); hipClosest2JetB->Sumw2(); hipClosest2JetC->Sumw2(); hipClosest2JetL->Sumw2();

  Double_t t_jtpt, t_jteta, t_jtphi, t_rawpt, t_refpt,  t_refeta, t_subid, t_discr_prob, t_discr_ssvHighEff, t_discr_ssvHighPur, t_discr_csvSimple, t_svtxm;
Int_t t_chargedN, t_neutralN, t_photonN;
Double_t t_chargedSum, t_photonSum, t_neutralSum, t_chargedMax, t_neutralMax,t_photonMax, t_eSum, t_muSum;
  Double_t t_pthat, t_weight, t_vz;
  Int_t t_refparton_flavorForB;
  Int_t trigIndex, t_bin, t_run,  t_nref;
  Double_t         t_hiHFplusEta4;
  Double_t         t_hiHFminusEta4;
  Int_t t_nIP;
  Double_t t_ipPt[100], t_ipProb0[100];
 Int_t t_ipJetIndex[100];
  Int_t t_nTrk;
  Double_t t_trkPt, t_trkEta, t_trkPhi;

  TTree *nt = new TTree("nt","");

  nt->Branch("nref",&t_nref, "nref/I");
  nt->Branch("jtpt",&t_jtpt,"jtpt/D");
  nt->Branch("jteta",&t_jteta,"jteta/D");
  nt->Branch("jtphi",&t_jtphi,"jtphi/D");
  nt->Branch("chargedN",&t_chargedN,"chargedN/I");
  nt->Branch("photonN",&t_photonN,"photonN/I");
  nt->Branch("neutralN",&t_neutralN,"neutralN/I");
  nt->Branch("chargedSum",&t_chargedSum,"chargedSum/D");
  nt->Branch("neutralSum",&t_neutralSum,"neutralSum/D");
  nt->Branch("photonSum",&t_photonSum,"photonSum/D");
  nt->Branch("eSum",&t_eSum,"eSum/D");
  nt->Branch("muSum",&t_muSum,"muSum/D");
  nt->Branch("chargedMax",&t_chargedMax,"chargedMax/D");
  nt->Branch("photonMax",&t_photonMax,"photonMax/D");
  nt->Branch("neutralMax",&t_neutralMax,"neutralMax/D");
  nt->Branch("subid",&t_subid,"subid/I");
  nt->Branch("rawpt",&t_rawpt,"rawpt/D");
 if(isMC) nt->Branch("refpt",&t_refpt,"refpt/D");
 if(isMC) nt->Branch("refeta",&t_refeta,"refeta/D");
  nt->Branch("refparton_flavorForB",&t_refparton_flavorForB,"refparton_flavorForB/I");
  nt->Branch("discr_prob",&t_discr_prob,"discr_prob/D");
  nt->Branch("discr_ssvHighEff",&t_discr_ssvHighEff,"discr_ssvHighEff/D");
  nt->Branch("discr_ssvHighPur",&t_discr_ssvHighPur,"discr_ssvHighPur/D");
  nt->Branch("discr_csvSimple",&t_discr_csvSimple,"discr_csvSimple/D");
  nt->Branch("svtxm",&t_svtxm,"svtxm/D");
//  nt->Branch("bin",&t_bin,"bin/I");
  nt->Branch("hiBin",&t_bin,"hiBin/I");
  nt->Branch("evt",&evt, "evt/I");
  nt->Branch("run",&t_run,"run/I");
  nt->Branch("vz",&t_vz,"vz/D");
  nt->Branch("hiHFplusEta4",&t_hiHFplusEta4,"hiHFplusEta4/D");
  nt->Branch("hiHFminusEta4",&t_hiHFminusEta4,"hiHFminusEta4/D");
  if(ppPbPb){
    nt->Branch("trigIndex",&trigIndex,"trigIndex/I");
  }
  if(ExpandedTree){
    nt->Branch("nIP",&t_nIP);
    nt->Branch("ipPt",t_ipPt,"ipPt[nIP]/D");
    nt->Branch("ipProb0",t_ipProb0,"ipProb0[nIP]/D");
    nt->Branch("ipJetIndex",t_ipJetIndex,"ipJetIndex[nIP]/I");
  }
    if(doTracks){
    nt->Branch("nTrk",&t_nTrk, "nTrk/I");
    nt->Branch("trkPt",&t_trkPt,"trkPt/D");
    nt->Branch("trkPhi",&t_trkPhi,"trkPhi/D");
    nt->Branch("trkEta",&t_trkEta,"trkEta/D");
  }
  if(!ppPbPb && !isMC){
      nt->Branch("HLT_PAZeroBiasPixel_SingleTrack_v1 ",&HLT_PAZeroBiasPixel_SingleTrack_v1,"HLT_PAZeroBiasPixel_SingleTrack_v1/I");
    nt->Branch("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_NoJetID_v1,"HLT_PAJet20_noJetID_v1/I");
    nt->Branch("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_NoJetID_v1,"HLT_PAJet40_noJetID_v1/I");
    nt->Branch("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_NoJetID_v1,"HLT_PAJet60_noJetID_v1/I");
    nt->Branch("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_NoJetID_v1,"HLT_PAJet80_noJetID_v1/I");
    nt->Branch("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_NoJetID_v1,"HLT_PAJet100_noJetID_v1/I");
    nt->Branch("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl ",&HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl,"HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl/I");
    nt->Branch("HLT_PAJet20_noJetID_v1_Prescl",&HLT_PAJet20_NoJetID_v1_Prescl,"HLT_PAJet20_noJetID_v1_Prescl/I");
    nt->Branch("HLT_PAJet40_noJetID_v1_Prescl",&HLT_PAJet40_NoJetID_v1_Prescl,"HLT_PAJet40_noJetID_v1_Prescl/I");
    nt->Branch("HLT_PAJet60_noJetID_v1_Prescl",&HLT_PAJet60_NoJetID_v1_Prescl,"HLT_PAJet60_noJetID_v1_Prescl/I");
    nt->Branch("HLT_PAJet80_noJetID_v1_Prescl",&HLT_PAJet80_NoJetID_v1_Prescl,"HLT_PAJet80_noJetID_v1_Prescl/I");
    nt->Branch("HLT_PAJet100_noJetID_v1_Prescl",&HLT_PAJet100_NoJetID_v1_Prescl,"HLT_PAJet100_noJetID_v1_Prescl/I");	
    nt->Branch("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA,"pPAcollisionEventSelectionPA/I");
    nt->Branch("pVertexFilterCutGplusUpsPP",&pVertexFilterCutGplusUpsPP,"pVertexFilterCutGplusUpsPP/I");
    nt->Branch("pVertexFilterCutGplus",&pVertexFilterCutGplus,"pVertexFilterCutGplus/I");
    nt->Branch("pBeamScrapingFilter",&pBeamScrapingFilter,"pBeamScrapingFilter/I");
    nt->Branch("pprimaryvertexFilter",&pprimaryvertexFilter,"pprimaryvertexFilter/I");
    nt->Branch("phfPosFilter1",&phfPosFilter1,"phfPosFilter1/I");
    nt->Branch("phfNegFilter1",&phfNegFilter1,"phfNegFilter1/I");
    nt->Branch("pHBHENoiseFilter",&pHBHENoiseFilter,"pHBHENoiseFilter/I");
  }

  if(isMC) nt->Branch("pthat",&t_pthat,"pthat/D");
  nt->Branch("weight",&t_weight,"weight/D");

  // grab the JEC's
   
/*  if(updateJEC){
   string L2Name, L3Name; 
    JetCorrectorParameters* parHI44x_l2, * parHI44x_l3;
    vector<JetCorrectorParameters> vpar_HI44x;
    FactorizedJetCorrector *_JEC_HI44X = new FactorizedJetCorrector(vpar_HI44x);//JR

     
    //cout<<" updating the JECs, USING REGPF "<<endl;

    //string L2Name = "JEC/JEC_regPF_L2Relative_AK3PF.txt";
    //string L3Name = "JEC/JEC_regPF_L3Absolute_AK3PF.txt";
    if(algo=="akPu3PF"){
      L2Name = "/afs/cern.ch/work/y/ymao/pA/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/txt/PbpJEC/JEC_dijet_L2Relative_AKPu3PF.txt";
      L3Name = "/afs/cern.ch/work/y/ymao/pA/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/txt/PbpJEC/JEC_dijet_L3Absolute_AKPu3PF.txt";
   }
  else { 
     L2Name = "/afs/cern.ch/work/y/ymao/pA/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/txt/PbpJEC/JEC_dijet_L2Relative_AK3PF.txt";
     L3Name = "/afs/cern.ch/work/y/ymao/pA/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/txt/PbpJEC/JEC_dijet_L3Absolute_AK3PF.txt";
   }  
 
      parHI44x_l2 = new JetCorrectorParameters(L2Name.c_str());
      parHI44x_l3 = new JetCorrectorParameters(L3Name.c_str());

      vpar_HI44x.push_back(*parHI44x_l2);
      vpar_HI44x.push_back(*parHI44x_l3);
      _JEC_HI44X = new FactorizedJetCorrector(vpar_HI44x);

          
  }
  */
 
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  vector<string> listVector;
 if(!isMC){
    for(int ifile=0; ifile<10292; ifile++){
       instr >> filename;
       listVector.push_back(filename);
    }
 } 

  int nFiles=0;
  if(ppPbPb) nFiles=1;
  else if(isMC){
    endfile=QCDpthatBins ;
   startfile=0 ; 
    nFiles=QCDpthatBins;
  }
  else{
//    nFiles=dataFiles;
   if(endfile>0 && endfile<dataFiles ) nFiles = endfile ; 
  else nFiles = dataFiles ; 
  }
     float w = 1.;

    TTree *t;
    TTree *tEvt = NULL;
    TTree *tHlt = NULL;
    TTree *tTrk = NULL;
      TBranch* tweight;
//   vector<string> listVector;
  for(int ifile = startfile; ifile<nFiles; ifile++){
    //Add b/c statistics to the HF statistics
   // if(!ppPbPb){
   if(isMC){ 
     if(( ifile<QCDpthatBins) || !isMC){
        instr >> filename;
      }
      fin = TFile::Open(filename.c_str(), "readonly");
   std::cout << "File: " << filename.c_str() << std::endl;
   }
  else { 
      fin = TFile::Open(listVector.at(ifile).c_str(), "readonly");
   std::cout << "File: " << listVector.at(ifile).c_str() << std::endl;
    }
  //  if(usePUsub) 
       t = (TTree*) fin->Get(Form("%sJetAnalyzer/t", algo.Data()));
  //  else t = (TTree*) fin->Get("ak3PFJetAnalyzer/t");
    TTree *tSkim = (TTree*) fin->Get("skimanalysis/HltTree");
    if(!ppPbPb) tEvt = (TTree*) fin->Get("hiEvtAnalyzer/HiTree");
    if(!ppPbPb) tHlt = (TTree*) fin->Get("hltanalysis/HltTree");
    if(!ppPbPb) tTrk = (TTree*) fin->Get("ppTrack/trackTree");
    if(!t || !tSkim || (!tEvt&&!ppPbPb) || (!tHlt&&!ppPbPb)){ cout << "Error! Can't find one of the trees!" << endl; exit(0);}
     
    if(tEvt) t->AddFriend("hiEvtAnalyzer/HiTree");
    if(tSkim) t->AddFriend("skimanalysis/HltTree");
    if(tHlt) t->AddFriend("hltanalysis/HltTree");
    if(tTrk) t->AddFriend("ppTrack/trackTree");
     
    Long64_t nentries = t->GetEntries();

    t->SetBranchAddress("evt",&evt);
    t->SetBranchAddress("lumi",&lumi);
   // t->SetBranchAddress("run",&run);
  //  if(cbin != -1 || ppPbPb) 
     t->SetBranchAddress("hiBin",&bin);
    if(!isMC) t->SetBranchAddress("run",&run);
    if(ppPbPb) t->SetBranchAddress("hf",&hf);
    t->SetBranchAddress("vz",&vz);
    if(isMC) t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nref",&nref);
    t->SetBranchAddress("rawpt",&rawpt);
    t->SetBranchAddress("jtpt",&jtpt);
    t->SetBranchAddress("jteta",&jteta);
    t->SetBranchAddress("jtphi",&jtphi);
    if(!ppPbPb){
      t->SetBranchAddress("jtpu",&jtpu);
      t->SetBranchAddress("jty",&jty);
    }
    t->SetBranchAddress("discr_ssvHighEff",discr_ssvHighEff);
    t->SetBranchAddress("discr_ssvHighPur",discr_ssvHighPur);
    //t->SetBranchAddress("discr_csvMva",discr_csvMva);
    t->SetBranchAddress("discr_csvSimple",discr_csvSimple);
    //t->SetBranchAddress("discr_muByIp3",discr_muByIp3);
    t->SetBranchAddress("discr_muByPt",discr_muByPt);
    t->SetBranchAddress("discr_prob",discr_prob);

    t->SetBranchAddress("discr_probb",discr_probb);
    t->SetBranchAddress("discr_tcHighEff",discr_tcHighEff);
    t->SetBranchAddress("discr_tcHighPur",discr_tcHighPur);
    t->SetBranchAddress("nsvtx",nsvtx);
    t->SetBranchAddress("svtxntrk",svtxntrk);
    t->SetBranchAddress("svtxdl",svtxdl);
    t->SetBranchAddress("svtxdls",svtxdls);
    t->SetBranchAddress("svtxpt",svtxpt);
    
    t->SetBranchAddress("svtxm",svtxm);

    t->SetBranchAddress("nIPtrk",nIPtrk);
    t->SetBranchAddress("nselIPtrk",nselIPtrk);
    t->SetBranchAddress("nIP",&nIP);

    if(doTracks){
      t->SetBranchAddress("ipJetIndex",ipJetIndex);
      t->SetBranchAddress("ipPt",ipPt);
      t->SetBranchAddress("ipProb0",ipProb0);
      //t->SetBranchAddress("ipProb1",ipProb1);
      t->SetBranchAddress("ip2d",ip2d);
      t->SetBranchAddress("ip2dSig",ip2dSig);
      t->SetBranchAddress("ip3d",ip3d);
      t->SetBranchAddress("ip3dSig",ip3dSig);
      t->SetBranchAddress("ipDist2Jet",ipDist2Jet);
      //t->SetBranchAddress("ipDist2JetSig",ipDist2JetSig);
      t->SetBranchAddress("ipClosest2Jet",ipClosest2Jet);

      t->SetBranchAddress("nTrk",&nTrk);
      t->SetBranchAddress("trkPt",&trkPt);
      t->SetBranchAddress("trkPhi",&trkPhi);
      t->SetBranchAddress("trkEta",&trkEta);     
    }

    t->SetBranchAddress("mupt",mupt);
    if(ppPbPb) t->SetBranchAddress("muptPF",muptPF);
    if(!ppPbPb && !isMC) {
   //    t->SetBranchAddress("pVertexFilterCutGplusUpsPP",&pVertexFilterCutGplusUpsPP);
       t->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
       t->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
       t->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
       t->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
       t->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
       t->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    }
    /*
t->SetBranchAddress("mue",mue);
t->SetBranchAddress("mueta",mueta);
t->SetBranchAddress("muphi",muphi);
t->SetBranchAddress("mudr",mudr);
t->SetBranchAddress("muptrel",muptrel);
t->SetBranchAddress("muchg",muchg);
*/
        tweight = t->GetBranch("weight");
    if(isMC){
      t->SetBranchAddress("pthat",&pthat);
      t->SetBranchAddress("beamId1",&beamId1);
      t->SetBranchAddress("beamId2",&beamId2);
      t->SetBranchAddress("refpt",refpt);
      t->SetBranchAddress("subid",subid);
      t->SetBranchAddress("refeta",refeta);
      t->SetBranchAddress("refy",refy);
      t->SetBranchAddress("refphi",refphi);
      t->SetBranchAddress("refdphijt",refdphijt);
      t->SetBranchAddress("refdrjt",refdrjt);
      t->SetBranchAddress("refparton_pt",refparton_pt);
      t->SetBranchAddress("refparton_flavor",refparton_flavor);
      t->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);
      t->SetBranchAddress("refparton_isGSP",refparton_isGSP);
      

    //  if(isMC){
        if(!tweight){
         if(ifile==0){
         cout << "Weight not found in Tree. Calculating..." << endl;
         useWeight=0;
         }
        }
        if(!ppPbPb && !useWeight && ifile==0){
       //  MCentr = countMCevents(infile, usePUsub, isMC);
         MCentr = countMCevents(infile, isMC);
        // MCentr = (int) nentries;
         // if(isMC>1){
         // for(int lm=HFpthatBins+2; lm<QCDpthatBins+1; lm++){
         // MCentr[HFpthatBins] += MCentr[lm]; //hack because we go to pthat bin 540 in QCD jet and only pthat bin 170 in b/c jet MC
         // }
         // }
         for(int i=0; i<QCDpthatBins; i++){
         cout << "MCentr["<<i<<"]: " << *(MCentr+i) << endl;
         }
        }
  //    }
    }
    t->SetBranchAddress("chargedN",chargedN);
    t->SetBranchAddress("photonN",photonN);
    t->SetBranchAddress("neutralN",neutralN);
    t->SetBranchAddress("chargedMax",chargedMax);
    t->SetBranchAddress("photonMax",photonMax);
    t->SetBranchAddress("neutralMax",neutralMax);
    t->SetBranchAddress("chargedSum",chargedSum);
    t->SetBranchAddress("photonSum",photonSum);
    t->SetBranchAddress("neutralSum",neutralSum);
    t->SetBranchAddress("muSum",muSum);
    t->SetBranchAddress("eSum",eSum);
    
    if(isMC&&useWeight){
      t->SetBranchAddress("xSecWeight",&xSecWeight);
      t->SetBranchAddress("centWeight",&centWeight);
      t->SetBranchAddress("vzWeight",&vzWeight);
    }

    if(ppPbPb){
      t->SetBranchAddress("nHLTBit",&nHLTBit);
      t->SetBranchAddress("hltBit",hltBit);
      
      tSkim->SetBranchAddress("pvSel",&pvSel);
      tSkim->SetBranchAddress("hbheNoiseSel",&hbheNoiseSel);
      tSkim->SetBranchAddress("spikeSel",&spikeSel);
      tSkim->SetBranchAddress("collSell",&collSell);
    }
    if(!ppPbPb){
      t->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
      t->SetBranchAddress("HLT_PAJet20_NoJetID_v1",&HLT_PAJet20_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&HLT_PAJet40_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&HLT_PAJet60_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&HLT_PAJet80_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet100_NoJetID_v1",&HLT_PAJet100_NoJetID_v1);
      t->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl);
      t->SetBranchAddress("HLT_PAJet20_NoJetID_v1_Prescl",&HLT_PAJet20_NoJetID_v1_Prescl);
      t->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&HLT_PAJet40_NoJetID_v1_Prescl);
      t->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&HLT_PAJet60_NoJetID_v1_Prescl);
      t->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&HLT_PAJet80_NoJetID_v1_Prescl);
      t->SetBranchAddress("HLT_PAJet100_NoJetID_v1_Prescl",&HLT_PAJet100_NoJetID_v1_Prescl);
      t->SetBranchAddress("HLT_PAJet20_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[0]);
            t->SetBranchAddress("HLT_PAJet40_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[1]);
            t->SetBranchAddress("HLT_PAJet60_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[2]);
            t->SetBranchAddress("HLT_PAJet80_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[3]);
            t->SetBranchAddress("HLT_PAJet100_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[4]);
      t->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
      t->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
      t->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
      t->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4);
      t->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4);
    }


    int gspCounter=0;
    //nentries=10;
    for (Long64_t i=0; i<nentries;i++) {
       
      if (i%1000000==0) cout<<" i = "<<i<<" out of "<<nentries<<" ("<<(int)(100*(float)i/(float)nentries)<<"%)"<<endl;
      
      tSkim->GetEntry(i);
      t->GetEntry(i);
      tEvt->GetEntry(i);
      if(ppPbPb && isMC){
        // temporarily remove cuts from MC
        if(!pvSel||!spikeSel) continue; //hbheNoise doesn't work in mixed events
      }
      if(ppPbPb){
        //if(!pvSel||!hbheNoiseSel||!spikeSel) continue;
        // turn off spike and on coll Sel
        if(!pvSel||!hbheNoiseSel||!collSell){
         //cout<<" selection failed, pvSel="<<pvSel<<", hbheNoiseSel="<<hbheNoiseSel<<" , collSell="<<collSell<<endl;
         continue;
        }
      }

      //Cut to remove events that correspond to the twiki "good events" but not the golden lumi filter
      if(!isMC){
      //  if(((int)run>211256 ) || ((int)run<210676 )) continue;
       if(coll=="PPb"){
        if(((int)run>211256 ) ) continue; //! pPb run
       }
      else {
        if(((int)run<=211256 ) ) continue; //! Pbp run
       }  
      }

       double trgPrescl[5] = {(double)HLT_PAJet20_NoJetID_v1_Prescl, (double)HLT_PAJet40_NoJetID_v1_Prescl, (double)HLT_PAJet60_NoJetID_v1_Prescl, (double)HLT_PAJet80_NoJetID_v1_Prescl, (double)HLT_PAJet100_NoJetID_v1_Prescl};
        bool trgDec[5] = {(bool)HLT_PAJet20_NoJetID_v1, (bool)HLT_PAJet40_NoJetID_v1, (bool)HLT_PAJet60_NoJetID_v1, (bool)HLT_PAJet80_NoJetID_v1, (bool)HLT_PAJet100_NoJetID_v1};
      if(!ppPbPb){
        if(!isMC){
        //  if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA || !pVertexFilterCutGplus) continue;
         if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA ) continue;
    //  if(!HLT_PAZeroBiasPixel_SingleTrack_v1 && !HLT_PAJet20_NoJetID_v1 && !HLT_PAJet40_NoJetID_v1 && !HLT_PAJet60_NoJetID_v1 && !HLT_PAJet80_NoJetID_v1 && !HLT_PAJet100_NoJetID_v1 ) continue;
      if(!HLT_PAJet20_NoJetID_v1 && !HLT_PAJet40_NoJetID_v1 && !HLT_PAJet60_NoJetID_v1 && !HLT_PAJet80_NoJetID_v1 && !HLT_PAJet100_NoJetID_v1 ) continue;
        }
     /*   else{
          if(!pHBHENoiseFilter || !pPAcollisionEventSelectionPA) continue;
        //  if(!pHBHENoiseFilter ) continue;
        }
       */
       //! for 0-90% centrality bins selection using HF sum energy
       t_hiHFplusEta4 = hiHFplusEta4 ;
       t_hiHFminusEta4 = hiHFminusEta4 ;
     // if((t_hiHFplusEta4+t_hiHFminusEta4)<2.87) continue ;
      }

      for(int j=0; j<5; j++){ trigObjSize[j] = HLT_PAJet_NoJetID_v1_trigObject[j]->size();}
         //Fill the trigger Pt/Eta/Phi from the TriggerObjects
         for(int ii=0; ii<5; ii++){
            for(unsigned int iObj=0; iObj<trigObjSize[ii]; iObj++){
               trigObjPt[ii][iObj] = HLT_PAJet_NoJetID_v1_trigObject[ii]->at(iObj).pt();
               trigObjEta[ii][iObj] = HLT_PAJet_NoJetID_v1_trigObject[ii]->at(iObj).eta();
               trigObjPhi[ii][iObj] = HLT_PAJet_NoJetID_v1_trigObject[ii]->at(iObj).phi();
            }
         }

/*         //Protection against random trigger benefit
                int maxtrg= -1;
             for(int ii=4; ii>=0; ii--){
                 if(trgDec[ii]==1){
                      maxtrg=ii;
                      break;
                    }
                }
          float maximumTrgPt=0;
                bool rtBenefitFlag=0;
                for(int ii=0; ii<5; ii++){
                    for(unsigned int iObj=0; iObj<trigObjSize[ii]; iObj++){
                        if(trigObjPt[ii][iObj]>=(20+maxtrg*20)) rtBenefitFlag=1;
                        if(trigObjPt[ii][iObj]>maximumTrgPt) maximumTrgPt = trigObjPt[ii][iObj];
                    }
                }
                if(!isMC && (maxtrg==-1 || !rtBenefitFlag)){
                    cout << "event skipped from RT!" << endl;
                    cout << "maxtrg: " << (20+20*maxtrg) << endl;
                    cout << "trgPt: " << maximumTrgPt << endl;
                    continue;
                }
                //end rtbenefit protection
  */        std::set<int> usedObjNos;
          std::set<int> usedtrgNos;

 
 /*     if(ppPbPb){
        if(cbin==-1){
         // do nothing
         t_bin=bin;
        }
        else if(cbin==0){
         if(bin>=8) continue;
         else t_bin=bin;
        }
        else if(cbin==1) {
         if(bin<8||bin>=20) continue;
         else t_bin=bin;
        }
        else if(cbin==2){
         if(bin<20) continue;
         else t_bin=bin;
        }
        else {
         cout<<" bin not defined "<<endl;
         return;
        }
      }
      else 
 */      t_bin=bin;
       t_run = run ;
       t_vz = vz ;
 
//      if(isMC&&!ppPbPb){
//        if(beamId1==2112 || beamId2==2112) continue;
//      }
      

   if(coll=="PPb" || coll=="PbP" ){  if(fabs(vz)>15.) continue;}
      
      
        
 /*     
      if(updateJEC){
        
        for(int ij=0; ij<nref; ij++){        
          _JEC_HI44X->setJetEta(jteta[ij]);
          _JEC_HI44X->setJetPt(rawpt[ij]);
          jtpt[ij] = rawpt[ij]*_JEC_HI44X->getCorrection();
        }
      }
   */   
    /*  if(useWeight){
        if(isMC)w=weight;
      }
   */   //trigger weighting in pp data
 //     if(!ppPbPb && !isMC){
      //  bool trgDec[6] = {(bool)HLT_PAZeroBiasPixel_SingleTrack_v1, (bool)HLT_PAJet20_NoJetID_v1, (bool)HLT_PAJet40_NoJetID_v1, (bool)HLT_PAJet60_NoJetID_v1, (bool)HLT_PAJet80_NoJetID_v1, (bool)HLT_PAJet100_NoJetID_v1};
 //     }

      
      //Do the weighting = x-sec / Nentries, where Nentries is weighted differently for B/C jets and QCD jets
      if(isMC){
        bool isFiltered=0;
        vzWeight=1;
        t_pthat=pthat;
        int j=0;
         
        while(pthat>pthatbin[j] && j<QCDpthatBins) {j++;
     //   while(j<QCDpthatBins) {
     //      j++;
     //  cout <<"j==" << j <<endl ; 
         if(j==QCDpthatBins) 
           xSecWeight = ((wght[j-1])/MCentr[j-1]); //wght[0] = pthat>15, MCentr[0] = pthat<15. I know it's dumb - bear with me.
         else 
          xSecWeight = ((wght[j-1]-wght[j])/MCentr[j-1]);
      //  }
       }
     
        w= xSecWeight ;
/*     // if((ifile<(QCDpthatBins-1) && pthat>pthatbin[ifile+1])) continue ;
      if(pthat<pthatbin[ifile]) continue ;
       if(ifile==QCDpthatBins-1) 
          // w = ((wght[ifile]-wght[ifile+1])/MCentr[ifile]); //wght[0] = pthat>15, MCentr[0] = pthat<15. I know it's dumb - bear with me.
           w = ((wght[ifile])/MCentr[ifile]); //wght[0] = pthat>15, MCentr[0] = pthat<15. I know it's dumb - bear with me.
         else { 
          w = ((wght[ifile]-wght[ifile+1])/MCentr[ifile]); 
          if((pthat>pthatbin[ifile+1])) continue ;
        }
  */      int vzbin = (int) TMath::Ceil(vz+15.+0.4); // 0.4 is the pixel detector shift
      //  if(vzbin>0&&vzbin<=30)
       //     vzWeight = hDatavz->GetBinContent(vzbin)/hMCvz[isFiltered]->GetBinContent(vzbin);
        if(coll=="PPb" || coll=="PbP" )  {
            vzWeight = fVz->Eval(vz);
             centWeight = fCen->Eval(bin);
          }
        t_weight=w*vzWeight*centWeight;
     //   t_weight=1;
          if(i%100000==0)cout <<" pthat =" << pthatbin[ifile] << "  cross section cs= " << wght[ifile] << " Nevt = " << MCentr[ifile] << " weight calculated cs/Nevt = " <<  xSecWeight <<endl ;
        hpthat->Fill(pthat);
        hpthatw->Fill(pthat, xSecWeight);
        hpthatcs->Fill(pthat, TMath::Log10(xSecWeight));
      }
    //  else t_weight=w;
   t_nref = nref;
   t_nTrk =  nTrk;

      int useEvent=0;
      
      int trackPosition =0;
  //    triggerPt=0;
       //lets get the leading jet information here
      int leadJetIndex = -1;
      float leadJetPt = minJetPt ;
     for(int ij=0;ij<nref;ij++){
        if(rawpt[ij]<22.||fabs(jteta[ij])>3.0) continue ;
        	float jetpt = jtpt[ij];
                if(jetpt>leadJetPt){
                 leadJetPt = jetpt;
                leadJetIndex = ij;
              }

   } //! end of leading jet search

  //old implementation with trigger Obj mathcing to offline jets
/*   for(int ij=0;ij<nref;ij++){
        if(rawpt[ij]<20.&&fabs(jteta[ij])>maxJetEta) continue ; 
    //to assigne the trigger pt according to its trigger combination method
    if(useJetTrgAssociation){
      for(int ii=0; ii<5; ii++){
        if(trgDec[ii]){
          triggerPt = triggerMatch(trigObjSize[ii], trigObjPhi[ii], trigObjEta[ii], trigObjPt[ii], jtphi[ij], jteta[ij], jtpt[ij]);
     //if you find a trigger match that has the right pt window for that particular trigger, break out of the for loop!
           if((triggerPt>(ii+1)*20 && triggerPt<(ii+2)*20) || (triggerPt>100 && ii==4)){
                  break;
             }
          }
        }
      }
  else{
    triggerPt = leadJetPt;
    }

  //      w = trigComb(trgDec, trgPrescl);
    w = trigComb(trgDec, trgPrescl, triggerPt, combinationMethod);
 } //end of leading jet search
//  triggerPt = leadJetPt;
*/
 // new implementation with max trigger Obj for weight calculation
    double trgPtWeight=0;
    double maxTrgPt =0;
      double maxpt=0;
      int maxtrg=-1;
  if(!isMC){
   if(useJetTrgAssociation){
    for(int ii=0; ii<5; ii++){
       if(trgDec[ii]){
       	  for(int isize=0; isize<trigObjSize[ii]; isize++){
	    double triggerPt = HLT_PAJet_NoJetID_v1_trigObject[ii]->at(isize).pt();
	    if(triggerPt > maxpt && (triggerPt-TMath::Floor(triggerPt)) > 0.0001){
	 //   if(triggerPt > maxpt ){
	      maxpt = triggerPt;
       //      cout <<" Kurt maxPt found " << maxpt << endl ;
	      maxtrg = ii;
	    }
	  }
	}      
      }
        
     if(maxtrg>=0 && maxtrg<5){
         maxTrgPt = maxpt;
     //   cout <<" maxTrgPt found " << maxTrgPt << endl ;
        trgPtWeight = trgPrescl[maxtrg];
        combinationMethod = 1 ;
    }
/*   else {
       maxTrgPt = leadJetPt;
   }
*/  }
   else {
       maxTrgPt = leadJetPt;
       if(maxTrgPt >40 && maxTrgPt <=60) trgPtWeight = trgPrescl[0];
       if(maxTrgPt >60 && maxTrgPt <=75) trgPtWeight = trgPrescl[1];
       if(maxTrgPt >75 && maxTrgPt <=95) trgPtWeight = trgPrescl[2];
       if(maxTrgPt >95 && maxTrgPt <=120) trgPtWeight = trgPrescl[3];
       if(maxTrgPt >120 ) trgPtWeight = trgPrescl[4];
      combinationMethod = 0 ;
   }
    w = trigComb(trgDec, trgPrescl, maxTrgPt, combinationMethod);
  // w*=corr_Qiao.getEventWeightHFPlus4bak(hiHFplusEta4,kTRUE); 
   t_weight=w;
    //   cout <<" Eric maxTrgPt found " << maxTrgPt <<" weight =" << t_weight << endl ;
 } //only for data weight
 
         //inclusive track 
        if(doTracks && doNtuples){
         for(int itrk=0; itrk<nTrk; itrk++){
                if(TMath::Abs(trkEta[itrk]) <2.4){
                 t_trkPt=(double)trkPt[itrk];
                 t_trkPhi=(double)trkPhi[itrk];
                 t_trkEta=(double)trkEta[itrk];
                }
         }
         }

   for(int ij=0;ij<nref;ij++){
      
        trackPosition+=nselIPtrk[ij];
      
//        if(useGSP==2){
//         if(refparton_isGSP[ij]==1){
//         gspCounter++;
//         if(gspCounter%2==0) continue;
//         }        
//        }
//        if(useGSP==3){
//         if(refparton_isGSP[ij]==0){
//         gspCounter++;
//         if(gspCounter%2==0) continue;
//         }        
//        }
          if(rawpt[ij]<22.||fabs(jteta[ij])>maxJetEta) continue ;
          double jetweight = 1; 
        //  if(isMC && subid[ij]!=0) continue ;
        //  if(isMC && jtpt[ij]>4*pthat) continue ;
      //  if(jtpt[ij]>minJetPt && fabs(jteta[ij])<maxJetEta){
        if(jtpt[ij]>minJetPt ){
            
         if(doNtuples){
        
         t_jtpt=jtpt[ij];
         t_jteta=jteta[ij];
         t_jtphi=jtphi[ij];
         t_rawpt=rawpt[ij];
         t_chargedN=chargedN[ij];
         t_photonN=photonN[ij];
         t_neutralN=neutralN[ij];
         t_neutralMax=neutralMax[ij];
         t_photonMax=photonMax[ij];
         t_chargedMax=chargedMax[ij];
         t_neutralSum=neutralSum[ij];
         t_photonSum=photonSum[ij];
         t_chargedSum=chargedSum[ij];
         t_eSum=eSum[ij];
         t_muSum=muSum[ij];
       //  if(coll=="PPb") jetweight*=C_rel->GetBinContent(C_rel->FindBin(jteta[ij]));

         if(isMC)t_refpt=refpt[ij];
         if(isMC)t_refeta=refeta[ij];
         t_subid=subid[ij];
         t_refparton_flavorForB=refparton_flavorForB[ij];
         t_discr_prob=discr_prob[ij];
         t_discr_ssvHighEff=discr_ssvHighEff[ij];
         t_discr_ssvHighPur=discr_ssvHighPur[ij];
         t_discr_csvSimple=discr_csvSimple[ij];
         t_svtxm=svtxm[ij];
        
         //Find jet tracks that correspond to the jet & apply proximity cuts
         if(ExpandedTree){
         t_nIP=nselIPtrk[ij];
         int counter=0;
         for(int itrk=0; itrk<nIP; itrk++){
                if(ipJetIndex[itrk] == ij){
                 t_ipProb0[counter] = ipProb0[itrk];
                 t_ipPt[counter] = ipPt[itrk];
                 t_ipJetIndex[counter] = ij;
                 counter++;
                }
         }
         }
         nt->Fill();


         }

        
         if(!doJets) continue;
        
        
         useEvent=1;
        
       //  hjtpt->Fill(jtpt[ij],w);
       if(isMC){
        if(coll=="PPb"){
        if(TMath::Abs((jteta[ij]+0.465))<=1.)hjtpt->Fill(jtpt[ij],w);
        if(TMath::Abs((refeta[ij]+0.465))<=1.)hrefpt->Fill(refpt[ij], w);
        }
       else {
          if(TMath::Abs((jteta[ij]))<=1.)hjtpt->Fill(jtpt[ij],w);
        if(TMath::Abs((refeta[ij]))<=1.)hrefpt->Fill(refpt[ij], w);
       }
      }
     else {
       if(coll=="PPb"){
        if(TMath::Abs((jteta[ij]+0.465))<=1.)hjtpt->Fill(jtpt[ij],w*jetweight);
       }
      else {
        if(TMath::Abs((jteta[ij]-0.465))<=1.)hjtpt->Fill(jtpt[ij],w);
       } 
        }
          Int_t dEtaBin = -1 ;
          for(Int_t ieta = 0 ; ieta <netabin; ieta++){
           if(isMC){
             if(coll=="PPb"){
                if((refeta[ij]+0.465)>deta[ieta]&&(refeta[ij]+0.465)<=deta[ieta+1]) dEtaBin = ieta ;
             }
            else {
            if((refeta[ij])>deta[ieta]&&(refeta[ij])<=deta[ieta+1]) dEtaBin = ieta ;
               } 
           }            
           else{
           if(coll=="PPb"){ 
             if((jteta[ij]+0.465)>deta[ieta]&&(jteta[ij]+0.465)<=deta[ieta+1]) dEtaBin = ieta ;
            }
            else {
             if((jteta[ij]-0.465)>deta[ieta]&&(jteta[ij]-0.465)<=deta[ieta+1]) dEtaBin = ieta ;
            }
           }
          }
         if(dEtaBin!=-1) hJetPtEtaBin[dEtaBin]->Fill(jtpt[ij],w*jetweight);

         hrawpt->Fill(rawpt[ij],w);

         hjteta->Fill(jteta[ij],w);
        if(isMC){
         hrefeta->Fill(refeta[ij],w);
         hetares->Fill(refeta[ij], (refeta[ij]-jteta[ij])/refeta[ij], w);
        } 
         hjtphi->Fill(jtphi[ij],w);

         //*
         hdiscr_csvSimple->Fill(discr_csvSimple[ij],w);

        
         hdiscr_prob->Fill(discr_prob[ij],w);

         hdiscr_ssvHighEff->Fill(discr_ssvHighEff[ij],w);
        
         hdiscr_ssvHighPur->Fill(discr_ssvHighPur[ij],w);

         hdiscr_tcHighEff->Fill(discr_tcHighEff[ij],w);
        
         hdiscr_tcHighPur->Fill(discr_tcHighPur[ij],w);

         //*
         hnsvtx->Fill(nsvtx[ij],w);
        
            if(nsvtx[ij]>0){
                
                hsvtxntrk->Fill(svtxntrk[ij],w);
                
                
                // require at least 1 tracks as in btagging @ 7 TeV note
                if(svtxntrk[ij]>1){
                    
                    hsvtxdl->Fill(svtxdl[ij],w);
                    
                    hsvtxdls->Fill(svtxdls[ij],w);
                    
                    
                    hsvtxm->Fill(svtxm[ij],w);
                    
                    
                    hsvtxpt->Fill(svtxpt[ij],w);
                    
                    
                    if(svtxntrk[ij]>=3) {
                        
                        hsvtxmSV3->Fill(svtxm[ij],w);
                        
                        
                        hsvtxptSV3->Fill(svtxpt[ij],w);
                        
                    }
                }
            }
        
         hnIPtrk->Fill(nIPtrk[ij],w);

        
        
         hnselIPtrk->Fill(nselIPtrk[ij],w);
        

        }  //! jet pt cuts
      }  //! end of jet loop
      if(useEvent){
        if(isMC){
         hbinw->Fill(bin,w);
         hbin->Fill(bin,xSecWeight*vzWeight);
        }
        else hbin->Fill(bin);
        
        if(isMC){
         hvzw->Fill(vz,w);
         if(ppPbPb)hvz->Fill(vz,xSecWeight*centWeight);
         else hvz->Fill(vz,xSecWeight);
        }
        else hvz->Fill(vz);
      }
      //reset event-level containers
      usedObjNos.clear();
      usedtrgNos.clear();
     for(int ii=0; ii<5; ii++){ HLT_PAJet_NoJetID_v1_trigObject[ii]->clear(); }       
    }  //! events loop
    //close file before moving to the next one!!
    fin->Close();
    t = NULL ;
    tEvt = NULL;
    tHlt = NULL;
    tweight = NULL; 
  } //! file loop
  fout->cd();

  hbin->Write(); hbinw->Write(); hvz->Write(); hvzw->Write();
  
  hjtpt->Write();
  if(isMC)hrefpt->Write();
  hpthat->Write();
  hpthatw->Write();
  hpthatcs->Write();
  for(Int_t ieta = 0 ; ieta <netabin; ieta++){
      hJetPtEtaBin[ieta]->Write();
    } 
  hrawpt->Write();
  
  hjteta->Write();
   hrefeta->Write();
  hetares->Write();
 
  hjtphi->Write();
  
  hdiscr_csvSimple->Write();
  
  hdiscr_prob->Write();
  
  hdiscr_ssvHighEff->Write();

  hdiscr_ssvHighPur->Write();

  hdiscr_tcHighEff->Write();

  hdiscr_tcHighPur->Write();

  hnsvtx->Write();

  hsvtxntrk->Write();
  
  hsvtxdl->Write();

  hsvtxdls->Write();

  hsvtxm->Write();

  hsvtxmSV3->Write();

  hsvtxpt->Write();

  hsvtxptSV3->Write();

  hnIPtrk->Write();

  hnselIPtrk->Write();

  hipPt->Write();

  hipProb0->Write();

  hipProb1->Write();

  hip2d->Write();

  hip2dSig->Write();

  hip2d1->Write();

  hip2dSig1->Write();

  hip2d2->Write();

  hip2dSig2->Write();

  hip2d3->Write();

  hip2dSig3->Write();

  hip3d->Write();

  hip3dSig->Write();

  hip3d1->Write();

  hip3dSig1->Write();

  hip3d2->Write();

  hip3dSig2->Write();

  hip3d3->Write();

  hip3dSig3->Write();


  hipDist2Jet->Write();

  hipDist2JetSig->Write();

  hipClosest2Jet->Write();

  nt->Write();
  
  fout->Close();

}



