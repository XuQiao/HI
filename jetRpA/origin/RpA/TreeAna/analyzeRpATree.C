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
#include "/home/nih/CMSSW_5_3_8_HI_patch2/src/NewTree/Corrections/Correction.C"
using namespace std;
const double deta[] = {-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,2.2};
const int netabin = sizeof(deta)/sizeof(double)-1;

TString algo="ak3PF";//"akPu3PF"
TString coll = "PPb";
TString TrigName = "Jet20" ; //"JetMB","Jet20","Jet40","Jet60","Jet80","Jet100"
double weight;

class hist_class{
public:
    hist_class(TString the_desc);
    void Delete();
    void Write();
   
    TString desc;

    bool IsMC;

    TH1D * NEvents;
    TH1D * NevtCounter;
    TH1D * JetAbove40;
    TH1D * JetAbove60;
    TH1D * JetAbove75;
    TH1D * JetAbove95;
    TH1D * JetAbove120;
  
    TH1F * jetpt;

    TH1F * jetptEtaBin[7];
};

hist_class::hist_class(TString the_desc)
{
  desc = the_desc;
  IsMC = kTRUE; 
   
  NEvents = new TH1D("Nevents_0-100%","Nevents_0-100%",100,0,2.); 
  NEvents -> Sumw2();
   
  NevtCounter = new TH1D("NeventsCounter_0-100%","NeventsCounter_0-100%",1,0,1);
  NevtCounter -> Sumw2();

  JetAbove40 = new TH1D("NeventsWithJetAbove40_0-100%","NeventsWithJetAbove40_0-100%",1,0,1); 
  JetAbove60 = new TH1D("NeventsWithJetAbove60_0-100%","NeventsWithJetAbove60_0-100%",1,0,1);
  JetAbove75 = new TH1D("NevntsWithJetAbove75_0-100%","NeventsWithJetAbove75_0-100%",1,0,1);
  JetAbove95 = new TH1D("NevntsWithJetAbove95_0-100%","NeventsWithJetAbove95_0-100%",1,0,1);
  JetAbove120 = new TH1D("NevntsWithJetAbove120_0-100%","NeventsWithJetAbove120_0-100%",1,0,1);
 
  jetpt = new TH1F("jetpt_0-100%","jetpt_0-100%",1000,0.,1000.);
  jetpt -> Sumw2();

  for(int ieta=0; ieta < netabin; ieta++){
    jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_Cen0-100%%",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f_Cen0-100%%",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
    jetptEtaBin[ieta]->Sumw2();
    }
}


void hist_class::Delete()
{
  delete NEvents;
  delete NevtCounter;
  delete JetAbove40;
  delete JetAbove60;
  delete JetAbove75;
  delete JetAbove95;
  delete JetAbove120;
  delete jetpt;
  for(int ieta=0; ieta < netabin; ieta++) delete jetptEtaBin[ieta];
}

void hist_class::Write()
{
  TString dataType;
  TString out_name;
  if(IsMC) dataType = "MC";
  else dataType = "DATA";

  out_name=Form("%s%s%s%sDiJetMassJet.root",dataType.Data(),coll.Data(),TrigName.Data(),algo.Data());
  
  TFile *out_file = new TFile(Form("/home/nih/CMSSW_5_3_8_HI_patch2/src/NewTree/%s",out_name.Data()),"RECREATE");  

  NEvents->Write();
  NevtCounter->Write();
  JetAbove40->Write();
  JetAbove60->Write();
  JetAbove75->Write();
  JetAbove95->Write();
  JetAbove120->Write(); 
  jetpt->Write();
  for(int ieta=0; ieta<netabin; ieta++) jetptEtaBin[ieta]->Write(); 
  
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void analyzeRpATree()
{
  std::cout<<"start working\n";
  cout<<"analyze trig = "<<TrigName<<endl; 
  hist_class *my_hists = new hist_class("pfjet");

  if(my_hists->IsMC==kTRUE) cout<<"DATA Type is MC!"<<endl;  
  if(my_hists->IsMC==kFALSE) cout<<"DATA Type is DATA!"<<endl;

  if(my_hists->IsMC==kTRUE)  TFile *f = new TFile("ppMCKurtForestNoBoost_ppReco_akPu3PF_QCDjetTrig_JetPt0noIPupperCut.root");
  if(my_hists->IsMC==kFALSE) TFile *f = new TFile("pPbdata_ppReco_akPu3PF_AlljetTrigKurtTrCombFile0_400_JetPt0noIPupperCut.root");  

  TTree *nt = (TTree*)f->Get("nt");

  Double_t jtpt[1000],jteta[1000],jtphi[1000];
  Double_t hiHFplusEta4;
  Int_t HLT_PAJet20_noJetID_v1;
  Int_t HLT_PAJet40_noJetID_v1;
  Int_t HLT_PAJet60_noJetID_v1;
  Int_t HLT_PAJet80_noJetID_v1;
  Int_t HLT_PAJet100_noJetID_v1;
  Int_t HLT_PAZeroBiasPixel_SingleTrack_v1;
  Int_t nref;
  Double_t MCweight;

  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("jtpt",jtpt);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("jtphi",jtphi);
  nt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
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

  TString CorrfileName;
  CorrfileName="/home/nih/CMSSW_5_3_8_HI_patch2/src/NewTree/Corrections/pPbHijing_EffCorr_VRw.root";
  Corrector corr_Qiao(CorrfileName);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
    nt->GetEntry(i);

    weight =1.;

    double HFbin = hiHFplusEta4;

  if(my_hists->IsMC==kFALSE){
    int jetTr2;
    if(TrigName=="Jet20") jetTr2 = HLT_PAJet20_noJetID_v1;
    else if (TrigName=="Jet40") jetTr2 = HLT_PAJet40_noJetID_v1;
    else if (TrigName=="Jet60") jetTr2 = HLT_PAJet60_noJetID_v1;
    else if (TrigName=="Jet80") jetTr2 = HLT_PAJet80_noJetID_v1;
    else if (TrigName=="Jet100")jetTr2 = HLT_PAJet100_noJetID_v1;
    else jetTr2 = HLT_PAZeroBiasPixel_SingleTrack_v1;
    if(!jetTr2) continue;
    } 
  
    double leadingJet = -999.;
    for(int j4i = 0; j4i < nref; j4i++) {
      if(jtpt[j4i] > leadingJet ) leadingJet = jtpt[j4i]; //find leading jet  
    }   

    bool jetAbove = false;
    bool jetAbove40 = false;
    bool jetAbove60 = false;
    bool jetAbove75 = false;
    bool jetAbove95 = false;
    bool jetAbove120 = false;
  
    if(leadingJet>40.) jetAbove40 = true;
    if(leadingJet>60.) jetAbove60 = true;
    if(leadingJet>75.) jetAbove75 = true;
    if(leadingJet>95.) jetAbove95 = true;
    if(leadingJet>120.) jetAbove120 = true; 

    if(!jetAbove40) Nevents++;
    else if(jetAbove40 && !jetAbove60) Nevt_40_60++;
    else if(jetAbove60 && !jetAbove75) Nevt_60_75++;
    else if(jetAbove75 && !jetAbove95) Nevt_75_95++;
    else if(jetAbove95 && !jetAbove120) Nevt_95_120++;
    else Nevt_120++;
  
    if(TrigName=="Jet20") jetAbove = jetAbove40 && !jetAbove60;
    else if(TrigName=="Jet40") jetAbove = jetAbove60 && !jetAbove75;
    else if(TrigName=="Jet60") jetAbove = jetAbove75 && !jetAbove95;
    else if(TrigName=="Jet80") jetAbove = jetAbove95 && !jetAbove120;
    else if(TrigName=="Jet100") jetAbove = jetAbove120;
    else jetAbove = !jetAbove40; 
    
    if(my_hists->IsMC==kFALSE)   if(!jetAbove) continue;//jet trigger cut

    if(i%1000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
 
    if(my_hists->IsMC==kFALSE)  weight*=corr_Qiao.getEventWeightHFPlus4bak(HFbin,kTRUE);//efficiency weight
    
    if(my_hists->IsMC==kTRUE) weight*=MCweight; // weight for MC


    my_hists->NEvents->Fill(1,weight);

    for(int j4i = 0; j4i < nref; j4i++){
      double jet_pt = jtpt[j4i];
      double jet_eta = jteta[j4i];
      int dEtaBin = -1.;
      if(coll=="PPb"){
        if(TMath::Abs(jet_eta+0.465)<=1.) my_hists->jetpt->Fill(jet_pt,weight);
      }

      for(int ieta =0; ieta<netabin;ieta++){
        if(coll=="PPb"){
          if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
        }
      }//assign the eta bin for jets
      if(dEtaBin!=-1) my_hists->jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
   }//jet loop

  my_hists->NevtCounter->SetBinContent(1,Nevents*weight);
  my_hists->JetAbove40->SetBinContent(1,Nevt_40_60*weight);
  my_hists->JetAbove60->SetBinContent(1,Nevt_60_75*weight);
  my_hists->JetAbove75->SetBinContent(1,Nevt_75_95*weight);
  my_hists->JetAbove95->SetBinContent(1,Nevt_95_120*weight);
  my_hists->JetAbove120->SetBinContent(1,Nevt_120*weight);

  }//loop over tree

  my_hists->Write();
  cout<<"working done\n";
 
}  









  











