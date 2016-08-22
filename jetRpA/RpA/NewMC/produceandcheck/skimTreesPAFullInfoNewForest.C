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

using namespace std;

//TStopwatch timer;

// ******* GLOBAL DECLARATIONS **********
const int QCDpthatBins = 10;
//const int QCDpthatBins = 9;
//const int QCDpthatBins = 1;
const int dataFiles = 7210;
//const int dataFiles = 4828;

int startfile ;
int endfile ;
int combinationMethod ;

TString coll = "PPb";
const TString algo = "akPu3Calo" ;
//const double deta[]={-2.5, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5} ;
const double deta[]={-1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5} ;
const int netabin = sizeof(deta)/sizeof(Double_t)-1 ;
const double jetPtBin[]={ 3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;
const int pthatbin[] = {15,30,50,80,120,170,220,280,370,460, 9999};
const double wght[] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 1.088E-07, 2.527E-08,0};

//**********************************************************
// Count the MC events to appropriately weight the pthat bins
//**********************************************************

int *countMCevents(std::string infile){

  TChain *ch = NULL;
    ch = new TChain(Form("%sJetAnalyzer/t", algo.Data()));
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  int *MCentries = new int[QCDpthatBins];
  int entrySum[QCDpthatBins]; 
  for(int ifile=0; ifile<QCDpthatBins; ifile++){
    instr >> filename;
    ch->Add(filename.c_str());
    if(ifile<QCDpthatBins-1) 
         entrySum[ifile] = ch->GetEntries(Form("pthat>=%d && pthat<%d",pthatbin[ifile], pthatbin[ifile+1]));
     else entrySum[ifile] = ch->GetEntries(Form("pthat>=%d",pthatbin[ifile]));  
    cout << "Summed MCentries[" << ifile << "]: " << entrySum[ifile] << endl;
    }
      for(int ifile=0; ifile<QCDpthatBins; ifile++){
     if(ifile==0) MCentries[ifile] = entrySum[ifile] ;
     else MCentries[ifile] = entrySum[ifile]; 
  }
 for(int i=0; i<QCDpthatBins; i++){
    cout << "QCD MCentries[" << i << "]: " << MCentries[i] << endl;
  }
  return MCentries;
}


//**********************************************************
// Trigger-Combine the data in order to unfold properly later
// need to adjust here to use the leading jet binning combination
//**********************************************************

//[0] = MB, [1] = Jet20, [2] = Jet40, [3] = Jet60, [4] = Jet80, [5] = Jet100
double trigComb(bool *trg, double *pscl, double pt, int combMethod){
  double weight=0;
 if(combMethod==1){
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

// ****************************************************************


void skimTreesPAFullInfoNewForest(int isMC=0)
 // isMC=0 --> Real data, ==1 --> QCD,
{
    if(isMC){
     startfile=0;
     endfile=QCDpthatBins ;
    }
	else{
    startfile=atoi(getenv("FIRST"));
    endfile=atoi(getenv("LAST"));}
   cout <<"analysis  "  << "  in file " << startfile << "  and " <<endfile<<endl ;
  int useWeight=1;

  //timer.Start();

    TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
    if(coll=="PPb")
   //     fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);
          fVz->SetParameters(1.47442e+00, -2.83100e-03, -1.19295e-02, 1.10312e-05, 2.64814e-05); //! new official MC >QM14
    else   if(coll=="PbP")
 //       fVz->SetParameters(1.54398e+00, -8.56155e-03, -1.40026e-02, 4.01020e-05, 3.47683e-05); //latest parameterization
  //      fVz->SetParameters(1.66731e+00,-2.43367e-03,-1.42488e-02,7.40147e-06,3.22477e-05);
           fVz->SetParameters(1.49736e+00, -6.93060e-03, -1.26864e-02, 2.98693e-05, 2.89538e-05); //! new official MC after QM14 
    else
        fVz->SetParameters(1.,0,0,0,0);
    
    TF1 * fCen = new TF1("fCen","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x+[5]*x*x*x*x+[6]*x*x*x*x*x)", 0., 100.);
    if(coll=="PPb")
     //   fCen->SetParameters(8.68073e-03, 5.09356e+00, -1.33053e-02, 1.46904e-03, -6.99681e-05, 1.06721e-06, -5.21398e-09);
//  fCen->SetParameters(1.20916e-02, 5.02157e+00, -3.38300e-02, 1.87647e-03, -6.76442e-05, 9.08602e-07, -4.01536e-09);//parameterize on 05.03 after approval
    //     fCen->SetParameters(7.92204e-03, 4.52005e+00, 9.77340e-02, -5.00362e-03, 9.74735e-05, -8.93897e-07, 3.39375e-09);//parameterize on new official MC after QM on 12/03/14
           fCen->SetParameters(6.06918e-03, 4.84487e+00, 4.26255e-02, -1.30682e-03, 1.94753e-05, -2.53606e-07, 1.61323e-09); //! parameterize on new official MC using double side HF on 02/02/15
    else if(coll=="PbP")
    //    fCen->SetParameters(1.05408e-02, 5.27477e+00, -8.03382e-02, 3.51669e-03, -8.85332e-05, 1.08917e-06, -4.90091e-09);
    //   fCen->SetParameters(2.89263e-02, 3.43643e+00, 5.62562e-02, -1.86555e-03, 1.97924e-06, 3.15416e-07, -1.97946e-09);//parameterize on new official MC after QM on 12/03/14
     fCen->SetParameters(5.10893e-03,4.88698e+00,8.37930e-02,-3.77127e-03, 7.90191e-05,-9.04877e-07, 4.26221e-09); //! parameterize on new official MC using double side HF on 02/02/15
    else
        fCen->SetParameters(1, 0, 0, 0, 0, 0, 0);
    

  TFile *fin=NULL;
  std::string infile;
   int *MCentr = NULL;

    if(!isMC){
             if(coll=="PPb")   infile = "ForwardPbDataFile.txt";
      // if(coll=="PPb")   infile = "test.txt";
              else  infile = "BackwardPbDataFile.txt";
    }
    else{
          if(coll=="PPb") infile = "OfficialMCForwardPb.txt";
      else if(coll=="PbP") infile = "OfficialMCBackwardPb.txt";
    //  else if(coll=="PbP") infile = "test.txt";
              else infile = "ppMCKurtForest.txt"; 
    }
 
  //Declaration of leaves types
  Int_t evt;
  Int_t run;
  Int_t bin;
  Int_t lumi;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
  Float_t vz;
  Int_t nref;
  Int_t ngen;
  Float_t rawpt[1000];
  Float_t jtpt[1000];
  Float_t jteta[1000];
  Float_t jty[1000];
  Float_t jtphi[1000];
  Float_t jtpu[1000];
  Int_t subid[1000];
  Float_t pthat;
  Float_t refpt[1000];
  Float_t refeta[1000];
  Float_t genpt[1000];
  Float_t geneta[1000];
  Float_t refy[1000];
  Float_t refphi[1000];
    Float_t pt;
    Float_t eta;
    Float_t phi;

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
  Int_t pprimaryVertexFilter;
  Int_t phfPosFilter3;
  Int_t phfNegFilter3;
  Int_t pHBHENoiseFilter;

  Int_t chargedN[1000];
  Int_t photonN[1000];
  Int_t neutralN[1000];
  Float_t chargedMax[1000];
  Float_t photonMax[1000];
  Float_t neutralMax[1000];
  Float_t chargedSum[1000];
  Float_t photonSum[1000];
  Float_t neutralSum[1000];
  Float_t muSum[1000];
  Float_t eSum[1000];
 
  // track variables
  Int_t t_nTrk;
  Int_t nTrk;
  Float_t t_trkPt[1000];
  Float_t trkPt[1000];
  Float_t t_trkPhi[1000];
  Float_t trkPhi[1000];
  Float_t t_trkEta[1000];
  Float_t trkEta[1000];
  Float_t         trkPtError[1000];   //[nTrk]
  Bool_t          highPurity[1000];   //[nTrk]
  Float_t         trkDxy1[1000];   //[nTrk]
  Float_t         trkDxyError1[1000];   //[nTrk]
  Float_t         trkDz1[1000];   //[nTrk]
  Float_t         trkDzError1[1000];   //[nTrk]
  Float_t weight, xSecWeight, centWeight, vzWeight;
 
  centWeight = 1. ; 

  TFile *fout=NULL;
  TString datatype ;
  if(isMC) datatype = "MC";
  else datatype = "DATA";
  if(isMC)
   fout=new TFile(Form("/cms/store/user/ymao/pA5TEV/Mixing/STARTHI53V27/merged/%s%s%sskimFullInfoLowerpthatCutfile%d_%d.root",datatype.Data(), coll.Data(), algo.Data(),startfile,endfile),"recreate");
   else 
   fout=new TFile(Form("/cms/store/user/ymao/pA5TEV/PAHighPt/HIRun2013-PromptReco-v1/skimfile/ForwardPb/%s%s%sskimFullInfofile%d_%d.root",datatype.Data(), coll.Data(), algo.Data(),startfile,endfile),"recreate");

  TTree *nt = new TTree("nt","");

//Jet Variables
  nt->Branch("nref",&nref, "nref/I");
  nt->Branch("ngen",&ngen, "ngen/I");
  nt->Branch("jtpt",jtpt,"jtpt[nref]/F");
  nt->Branch("jty",jteta,"jty[nref]/F");
  nt->Branch("jteta",jteta,"jteta[nref]/F");
  nt->Branch("jtphi",jtphi,"jtphi[nref]/F");
  nt->Branch("rawpt",rawpt,"rawpt[nref]/F");
if(isMC){  
     nt->Branch("subid",subid,"subid[nref]/I");
     nt->Branch("refpt",refpt,"refpt[nref]/F");
     nt->Branch("refeta",refeta,"refeta[nref]/F");
     nt->Branch("genpt",genpt,"genpt[ngen]/F");
     nt->Branch("geneta",geneta,"geneta[ngen]/F");
}

//JetID Variables
  nt->Branch("chargedN",chargedN,"chargedN[nref]/I");
  nt->Branch("neutralN",neutralN,"neutralN[nref]/I");
  nt->Branch("photonN",photonN,"photonN[nref]/I");
  nt->Branch("chargedMax",chargedMax,"chargedMax[nref]/F");
  nt->Branch("neutralMax",neutralMax,"neutralMax[nref]/F");
  nt->Branch("chargedMax",chargedMax,"chargedMax[nref]/F");
  nt->Branch("photonMax",photonMax,"photonMax[nref]/F");
  nt->Branch("chargedSum",chargedSum,"chargedSum[nref]/F");
  nt->Branch("neutralSum",neutralSum,"neutralSum[nref]/F");
  nt->Branch("photonSum",photonSum,"photonSum[nref]/F");
  nt->Branch("muSum",muSum,"muSum[nref]/F");
  nt->Branch("eSum",eSum,"eSum[nref]/F");

  //track variables

    nt->Branch("nTrk",&nTrk, "nTrk/I");
    nt->Branch("trkPt",trkPt,"trkPt[nTrk]/F");
    nt->Branch("trkPhi",trkPhi,"trkPhi[nTrk]/F");
    nt->Branch("trkEta",trkEta,"trkEta[nTrk]/F");

//Event Variables
  nt->Branch("hiBin",&bin,"hiBin/I");
  nt->Branch("evt",&evt, "evt/I");
  nt->Branch("run",&run,"run/I");
  nt->Branch("vz",&vz,"vz/F");
  nt->Branch("hiHFplusEta4",&hiHFplusEta4,"hiHFplusEta4/F");
  nt->Branch("hiHFminusEta4",&hiHFminusEta4,"hiHFminusEta4/F");
 if(!isMC){
    nt->Branch("pt",&pt,"pt/F");
    nt->Branch("eta",&eta,"eta/F");
    nt->Branch("phi",&phi,"phi/F");
    }
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
    nt->Branch("pVertexFilterCutGplus",&pVertexFilterCutGplus,"pVertexFilterCutGplus/I");
    nt->Branch("pBeamScrapingFilter",&pBeamScrapingFilter,"pBeamScrapingFilter/I");
    nt->Branch("pprimaryvertexFilter",&pprimaryVertexFilter,"pprimaryVertexFilter/I");
    nt->Branch("phfPosFilter1",&phfPosFilter3,"phfPosFilter3/I");
    nt->Branch("phfNegFilter1",&phfNegFilter3,"phfNegFilter3/I");
    nt->Branch("pHBHENoiseFilter",&pHBHENoiseFilter,"pHBHENoiseFilter/I");
    nt->Branch("weight",&weight,"weight/F");
   if(isMC) {
	nt->Branch("pthat",&pthat,"pthat/F");
	nt->Branch("xSecWeight",&xSecWeight,"xSecWeight/F");
	nt->Branch("vzWeight",&vzWeight,"vzWeight/F");
	nt->Branch("centWeight",&centWeight,"centWeight/F");
   }


  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  vector<string> listVector;
 if(!isMC){
    for(int ifile=0; ifile<dataFiles; ifile++){
       instr >> filename;
       listVector.push_back(filename);
    }
 } 
else {
   for(int ifile=0; ifile<QCDpthatBins; ifile++){
       instr >> filename;
       listVector.push_back(filename);
    }
 } 
  int nFiles=0;
  if(isMC){
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
    TTree *tObj = NULL;
    TTree *tTrk = NULL;
    TTree *tSkim = NULL;
      TBranch* tweight;

  for(int ifile = startfile; ifile<nFiles; ifile++){
/*      if(isMC){ 
     if(( ifile<QCDpthatBins)){
        instr >> filename;
      }
      fin = TFile::Open(filename.c_str(), "readonly");
   std::cout << "File: " << filename.c_str() << std::endl;
   }
  else { 
*/      fin = TFile::Open(listVector.at(ifile).c_str(), "readonly");
   std::cout << "File: " << listVector.at(ifile).c_str() << std::endl;
//    }

       t = (TTree*) fin->Get(Form("%sJetAnalyzer/t", algo.Data()));
    tSkim = (TTree*) fin->Get("skimanalysis/HltTree");
    tEvt = (TTree*) fin->Get("hiEvtAnalyzer/HiTree");
    tHlt = (TTree*) fin->Get("hltanalysis/HltTree");
      if(!isMC)tObj = (TTree*) fin->Get("hltobject/jetObjTree");
    tTrk = (TTree*) fin->Get("ppTrack/trackTree");
    if(!t || !tSkim || (!tEvt) || (!tHlt)){ cout << "Error! Can't find one of the trees!" << endl; exit(0);}
     
    if(tEvt) t->AddFriend("hiEvtAnalyzer/HiTree");
    if(tSkim) t->AddFriend("skimanalysis/HltTree");
    if(tHlt) t->AddFriend("hltanalysis/HltTree");
      if(tObj) t->AddFriend("hltobject/jetObjTree");
   if(tTrk) t->AddFriend("ppTrack/trackTree");
     
    Long64_t nentries = t->GetEntries();
    cout<<nentries<<endl;
    t->SetBranchAddress("evt",&evt);
    t->SetBranchAddress("lumi",&lumi);
   t->SetBranchAddress("hiBin",&bin);
 if(!isMC)  t->SetBranchAddress("run",&run);
 if(isMC){ 
      t->SetBranchAddress("ngen",&ngen);
      t->SetBranchAddress("weight",&weight);
      t->SetBranchAddress("pthat",&pthat);
      t->SetBranchAddress("refpt",refpt);
      t->SetBranchAddress("genpt",genpt);
      t->SetBranchAddress("subid",subid);
      t->SetBranchAddress("refeta",refeta);
      t->SetBranchAddress("geneta",geneta);
      t->SetBranchAddress("refy",refy);
      t->SetBranchAddress("refphi",refphi);
 }
    t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("nref",&nref);
    t->SetBranchAddress("rawpt",rawpt);
    t->SetBranchAddress("jtpt",jtpt);
    t->SetBranchAddress("jteta",jteta);
    t->SetBranchAddress("jtphi",jtphi);
      t->SetBranchAddress("jtpu",jtpu);
      t->SetBranchAddress("jty",jty);

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
    
      t->SetBranchAddress("nTrk",&t_nTrk);
      t->SetBranchAddress("trkPt",t_trkPt);
      t->SetBranchAddress("trkPhi",t_trkPhi);
      t->SetBranchAddress("trkEta",t_trkEta);
      t->SetBranchAddress("trkPtError",trkPtError);     	   
      t->SetBranchAddress("highPurity",highPurity);     	   
      t->SetBranchAddress("trkDz1",trkDz1);     	   
      t->SetBranchAddress("trkDzError1",trkDzError1);     	   
      t->SetBranchAddress("trkDxy1",trkDxy1);     	   
      t->SetBranchAddress("trkDxyError1",trkDxyError1);     	   
      if(!isMC){

      t->SetBranchAddress("pt",&pt);
      t->SetBranchAddress("eta",&eta);
      t->SetBranchAddress("phi",&phi);
      }
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
      t->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
       t->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
       t->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
       t->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
       t->SetBranchAddress("phfPosFilter3",&phfPosFilter3);
       t->SetBranchAddress("phfNegFilter3",&phfNegFilter3);
       t->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
      t->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4);
      t->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4);
	 tweight = t->GetBranch("weight");
      if(isMC){
        if(!tweight){
         if(ifile==0){
         cout << "Weight not found in Tree. Calculating..." << endl;
         useWeight=0;
         }
        }
        if(!useWeight && ifile==0){
         MCentr = countMCevents(infile);
         for(int i=0; i<QCDpthatBins; i++){
         cout << "MCentr["<<i<<"]: " << *(MCentr+i) << endl;
         }
   	}
	}

    for (Long64_t i=0; i<nentries;i++) {
       
      if (i%100000==0) cout <<" ifile ="<< ifile <<" i = "<<i<<" out of "<<nentries<<" ("<<(Float_t)(100*(Float_t)i/(Float_t)nentries)<<"%)"<<endl;
      
      tSkim->GetEntry(i);
      t->GetEntry(i);
      tEvt->GetEntry(i);
      //  if(!isMC)tObj->GetEntry(i);
 double trgPrescl[5] = {(double)HLT_PAJet20_NoJetID_v1_Prescl, (double)HLT_PAJet40_NoJetID_v1_Prescl, (double)HLT_PAJet60_NoJetID_v1_Prescl, (double)HLT_PAJet80_NoJetID_v1_Prescl, (double)HLT_PAJet100_NoJetID_v1_Prescl};
 bool trgDec[5] = {(bool)HLT_PAJet20_NoJetID_v1, (bool)HLT_PAJet40_NoJetID_v1, (bool)HLT_PAJet60_NoJetID_v1, (bool)HLT_PAJet80_NoJetID_v1, (bool)HLT_PAJet100_NoJetID_v1};
  //   if(!isMC){
       if(!HLT_PAJet20_NoJetID_v1 && !HLT_PAJet40_NoJetID_v1 && !HLT_PAJet60_NoJetID_v1 && !HLT_PAJet80_NoJetID_v1 && !HLT_PAJet100_NoJetID_v1 ) continue;
 //     }

    if(isMC){
        vzWeight=1;
        int j=0;
        
        while(pthat>pthatbin[j] && j<QCDpthatBins) {
          j++;
         if(j==QCDpthatBins) 
           xSecWeight = ((wght[j-1])/MCentr[j-1]);
         else 
          xSecWeight = ((wght[j-1]-wght[j])/MCentr[j-1]);
       }
          if(pthat<pthatbin[ifile]) continue ;  
	 // if(ifile<QCDpthatBins-1 && pthat>pthatbin[ifile+1]) continue ;
          w= xSecWeight ;
            vzWeight = fVz->Eval(vz);
             centWeight = fCen->Eval(bin);
      weight=w*vzWeight*centWeight;
   //   weight=w;
    }

		  // new implementation with max trigger Obj for weight calculation
    double maxTrgPt =0;
  if(!isMC){
      maxTrgPt = pt ;
      
    w = trigComb(trgDec, trgPrescl, maxTrgPt, 1);
  // w*=corr_Qiao.getEventWeightHFPlus4bak(hiHFplusEta4,kTRUE); 
   weight=w;
 } //only for data weight
	  nTrk=0;
         for(int itrk=0; itrk<t_nTrk; itrk++){
             if(((highPurity[itrk])
             && (fabs(trkDz1[itrk]/trkDzError1[itrk])<3)
             && (fabs(trkDxy1[itrk]/trkDxyError1[itrk])<3)
             && (trkPtError[itrk]/t_trkPt[itrk]<0.1)
            )){
                if(t_trkPt[itrk]>1. && fabs(t_trkEta[itrk])<2.4){
                 trkPt[nTrk]=t_trkPt[itrk];
                 trkPhi[nTrk]=t_trkPhi[itrk];
                 trkEta[nTrk]=t_trkEta[itrk];
		 nTrk++;
                }
            }
         }

     nt->Fill();
    }  //! events loop
    fin->Close();
    t = NULL ;
    tEvt = NULL;
    tSkim = NULL;
    tTrk = NULL;
    tObj= NULL;
    tHlt = NULL;
    tweight = NULL; 
  } //! file loop
  fout->cd();
  
//  nt->Write();
   fout->Write();
 
  fout->Close();
  cout << "work done!!!" <<endl ;
  nt = NULL ;
}
