#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TCut.h"

using namespace std;

TString coll = "PPb";
TString TrigName; //"JetMB","Jet20","Jet40","Jet60","Jet80","Jet100"
double weight;

class hist_class{
public:
    hist_class(TString the_desc);
    void Delete();
    void Write();
   
    TString desc;

    bool IsMC;

        TH1F* hEv;
        TH1F* hEvRun;
        TH1I* hRun;
        TH1F* hEvRunHLT;
        TH1F* hEvRunHLTBeam;
        TH1F* hEvRunHLTBeamppr;
        TH1F* hEvRunHLTBeampprHFp;
        TH1F* hEvRunHLTBeampprHFpHFn;
        TH1F* hEvRunHLTBeampprHFpHFnNoi;
        TH1F* hEvRunHLTBeampprHFpHFnNoivz;
        TH1F* hVz_cut1;
        TH1F* hVz_cut2;
        TH1F* hVz_cut3;
        TH1F* hVz_cut4;
        TH1F* hVz_cut7;

};

hist_class::hist_class(TString the_desc)
{
  desc = the_desc;
  IsMC = kFALSE; 
   
 hEv=new TH1F("hEv","hEv",4,-2,2);
 hEvRun=new TH1F("hEvRun","hEvRun",4,-2,2);
 hRun=new TH1I("hRun","hRun",2000,210000,212000);
 hEvRunHLT=new TH1F("hEvRunHLT","hEvRunHLT",4,-2,2);
 hEvRunHLTBeam=new TH1F("hEvRunHLTBeam","hEvRunHLTBeam",4,-2,2);
 hEvRunHLTBeamppr=new TH1F("hEvRunHLTBeamppr","hEvRunHLTBeamppr",4,-2,2);
 hEvRunHLTBeampprHFp=new TH1F("hEvRunHLTBeampprHFp","hEvRunHLTBeampprHFp",4,-2,2);
 hEvRunHLTBeampprHFpHFn=new TH1F("hEvRunHLTBeampprHFpHFn","hEvRunHLTBeampprHFpHFn",4,-2,2);
 hEvRunHLTBeampprHFpHFnNoi=new TH1F("hEvRunHLTBeampprHFpHFnNoi","hEvRunHLTBeampprHFpHFnNoi",4,-2,2);
 hEvRunHLTBeampprHFpHFnNoivz=new TH1F("hEvRunHLTBeampprHFpHFnNoivz","hEvRunHLTBeampprHFpHFnNoivz",4,-2,2);
         hVz_cut1=new TH1F("hVz_cut1","hVz_cut1",1000,-50.,50.);        hVz_cut1->Sumw2();
         hVz_cut2=new TH1F("hVz_cut2","hVz_cut2",1000,-50.,50.);        hVz_cut2->Sumw2();
         hVz_cut3=new TH1F("hVz_cut3","hVz_cut3",1000,-50.,50.);        hVz_cut3->Sumw2();
         hVz_cut4=new TH1F("hVz_cut4","hVz_cut4",1000,-50.,50.);        hVz_cut4->Sumw2();
         hVz_cut7=new TH1F("hVz_cut7","hVz_cut7",1000,-50.,50.);        hVz_cut7->Sumw2();
	
}


void hist_class::Delete()
{
         delete hEv;
         delete hEvRun;
         delete hRun;
         delete hEvRunHLT;
         delete hEvRunHLTBeam;
         delete hEvRunHLTBeamppr;
         delete hEvRunHLTBeampprHFp;
         delete hEvRunHLTBeampprHFpHFn;
         delete hEvRunHLTBeampprHFpHFnNoi;
         delete hEvRunHLTBeampprHFpHFnNoivz;
         delete hVz_cu11;
         delete hVz_cut2;
         delete hVz_cut3;
         delete hVz_cut4;
         delete hVz_cut7;
}

void hist_class::Write()
{
      TString out_name;
         out_name=Form("%sEvSel_runcut.root",TrigName.Data());

       TFile *out_file = new TFile(Form("/scratch/xuq7/RpA/GlobalEvent/EventSel/%s",out_name.Data()), "RECREATE");

         hEv->Write();
         hEvRun->Write();
         hRun->Write();
         hEvRunHLT->Write();
         hEvRunHLTBeam->Write();
         hEvRunHLTBeamppr->Write();
         hEvRunHLTBeampprHFp->Write();
         hEvRunHLTBeampprHFpHFn->Write();
         hEvRunHLTBeampprHFpHFnNoi->Write();
         hEvRunHLTBeampprHFpHFnNoivz->Write();
         hVz_cut1->Write();
         hVz_cut2->Write();
         hVz_cut3->Write();
         hVz_cut4->Write();
         hVz_cut7->Write();

  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void anaEv()
{
  TrigName = getenv("TRIG");
//  TrigName = "Jet20";   
  std::cout<<"start working\n";
  cout<<"analyze trig = "<<TrigName.Data()<<endl; 
  hist_class *my_hists = new hist_class("pfjet");

  if(my_hists->IsMC==kTRUE) cout<<"Analyzing MC!"<<endl;  
  if(my_hists->IsMC==kFALSE) cout<<"Analyzing DATA!"<<endl;

  if(my_hists->IsMC==kTRUE)  TFile *f = new TFile("ppMCKurtForestNoBoost_ppReco_akPu3PF_QCDjetTrig_JetPt0noIPupperCut.root");
  if(my_hists->IsMC==kFALSE){
/*	if(TrigName=="MB")	
	TFile *f = new TFile("pPbdata_ppreco_akPu3PF_MinBiasTrig_Merged.root");
else*/
	TFile *f = new TFile("/scratch/xuq7/RpA/TreeAna/PPbdata_ppReco_akPu3PF_AlljetTrigKurtTrCombFile_JetPt0noIPupperCut.root");  
}

  TTree *nt = (TTree*)f->Get("nt");

  Double_t jtpt[1000],jteta[1000],jtphi[1000];
  Double_t hiHFplusEta4;
  Int_t HLT_PAJet20_noJetID_v1;
  Int_t HLT_PAJet40_noJetID_v1;
  Int_t HLT_PAJet60_noJetID_v1;
  Int_t HLT_PAJet80_noJetID_v1;
  Int_t HLT_PAJet100_noJetID_v1;
  Int_t HLT_PAZeroBiasPixel_SingleTrack_v1;
  Int_t pPAcollisionEventSelectionPA;
  Int_t pprimaryvertexFilter;
  Int_t phfPosFilter1;
  Int_t phfNegFilter1;
  Int_t pHBHENoiseFilter;
  Int_t pBeamScrapingFilter;
  Int_t nref;
  Int_t run;
  Int_t bin;
  Int_t nTrk;
  Double_t MCweight,pthat;
  Double_t trkPt,trkPhi,trkEta;
  Double_t vz;
  Double_t rawpt,refpt,refeta;
  
  nt->SetBranchAddress("run",&run);
  nt->SetBranchAddress("bin",&bin);
  nt->SetBranchAddress("vz",&vz); 
  nt->SetBranchAddress("pthat",&pthat);  

  nt->SetBranchAddress("rawpt",&rawpt);
  nt->SetBranchAddress("refpt",&refpt);
  nt->SetBranchAddress("refeta",&refeta);  

  nt->SetBranchAddress("nTrk",&nTrk);
  nt->SetBranchAddress("trkPt",&trkPt);
  nt->SetBranchAddress("trkPhi",&trkPhi);
  nt->SetBranchAddress("trkEta",&trkEta);
  
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("jtpt",jtpt);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("jtphi",jtphi);
  nt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);	
  nt->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
  nt->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
  nt->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
  nt->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
  nt->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  nt->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
  nt->SetBranchAddress("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_noJetID_v1);  
  nt->SetBranchAddress("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_noJetID_v1);
//  nt->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);  
  nt->SetBranchAddress("weight",&MCweight);

  Int_t Nevents = 0;
  Int_t Nevt_40_60 = 0;
  Int_t Nevt_60_75 = 0;
  Int_t Nevt_75_95 = 0;
  Int_t Nevt_95_120 = 0;
  Int_t Nevt_120 = 0;

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
//    for(int i=0; i<10000; i++){
    nt->GetEntry(i);
my_hists->hEv->Fill(1);
    int jetTr2;
    if(TrigName=="Jet20") jetTr2 = HLT_PAJet20_noJetID_v1;
    else if (TrigName=="Jet40") jetTr2 = HLT_PAJet40_noJetID_v1;
    else if (TrigName=="Jet60") jetTr2 = HLT_PAJet60_noJetID_v1;
    else if (TrigName=="Jet80") jetTr2 = HLT_PAJet80_noJetID_v1;
    else if (TrigName=="Jet100")jetTr2 = HLT_PAJet100_noJetID_v1;
    else jetTr2 = HLT_PAZeroBiasPixel_SingleTrack_v1;
      
my_hists->hRun->Fill(run);
if(run>211256) continue;//only for pPb runs at moment
if(run<210676) continue ;
my_hists->hEvRun->Fill(1);
        my_hists->hVz_cut1->Fill(vz);

    if(!jetTr2) continue;
	my_hists->hEvRunHLT->Fill(1);
        my_hists->hVz_cut2->Fill(vz);

    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
 
           if(!(pBeamScrapingFilter)) continue ;
        my_hists->hEvRunHLTBeam->Fill(1);
        my_hists->hVz_cut3->Fill(vz);
           if(!(pprimaryvertexFilter)) continue ;
        my_hists->hVz_cut4->Fill(vz);
        my_hists->hEvRunHLTBeamppr->Fill(1);
           if(!(phfPosFilter1)) continue ;
        my_hists->hEvRunHLTBeampprHFp->Fill(1);
           if(!(phfNegFilter1)) continue ;
  //	if(!(pPAcollisionEventSelectionPA)) continue;
        my_hists->hEvRunHLTBeampprHFpHFn->Fill(1);
          if(!pHBHENoiseFilter) continue ;
        my_hists->hEvRunHLTBeampprHFpHFnNoi->Fill(1);
        my_hists->hVz_cut7->Fill(vz);
        if(TMath::Abs(vz)>15.) continue ;
        my_hists->hEvRunHLTBeampprHFpHFnNoivz->Fill(1);
}
  my_hists->Write();
  cout<<"working done\n";
 
}  
