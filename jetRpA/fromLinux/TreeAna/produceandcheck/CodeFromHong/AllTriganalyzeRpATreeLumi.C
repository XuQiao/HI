// this code is modified to use data only
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
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/Correction.C"
using namespace std;
const double deta[] = {-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,2.2};
const int netabin = sizeof(deta)/sizeof(double)-1;
TString TrigName[] = {"Jet20","Jet40","Jet60","Jet80","Jet100"};
const int trigRes[]={20,55,75,95,120,999};

TString algo="akPu3PF";//"akPu3PF"
TString coll = "PPb";
double weight;

class hist_class{
public:
    hist_class();
    void Delete();
    void Write();

    TH1F * jetpt[5];
    TH2F * jetptEta[5];
    TH2F * jetptphi[5];
    TH2F * rawptJES[5];
    TH2F * rawptJESEtaBin[5][7];
	
    TH1F * jetptEtaBin[5][7];
};

hist_class::hist_class()
{
  for(int itrig=0; itrig< 5; itrig++){
    jetpt[itrig] = new TH1F(Form("%sjetpt_0-100%%",TrigName[itrig].Data()),Form("%sjetpt_0-100%%",TrigName[itrig].Data()),1000,0.,1000.);
    jetpt[itrig] -> Sumw2(); 
    jetptEta[itrig] = new TH2F(Form("%sjetptEta_0-100%%",TrigName[itrig].Data()),Form("%sjetptEta_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,100,-5.,5.);
    jetptEta[itrig] -> Sumw2(); 
    jetptphi[itrig] = new TH2F(Form("%sjetptphi_0-100%%",TrigName[itrig].Data()),Form("%sjetptphi_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi[itrig] -> Sumw2(); 
    rawptJES[itrig] = new TH2F(Form("%srawptJES_0-100%%",TrigName[itrig].Data()),Form("%srawptJES_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,200,0.,2.);
    rawptJES[itrig] -> Sumw2(); 
    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[itrig][ieta] = new TH1F(Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[itrig][ieta]->Sumw2();
      rawptJESEtaBin[itrig][ieta] = new TH2F(Form("%srawptJESEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%srawptJESEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[itrig][ieta]->Sumw2();
      }
    }
}


void hist_class::Delete()
{
  for(int itrig=0; itrig < 5; itrig++){
    delete jetpt[itrig];
    delete jetptEta[itrig];
    delete jetptphi[itrig];
    delete rawptJES[itrig];
    for(int ieta=0; ieta < netabin; ieta++){
	 delete jetptEtaBin[itrig][ieta];
	 delete rawptJESEtaBin[itrig][ieta];
	}
    }
}

void hist_class::Write()
{
  TString dataType;
  TString out_name;
  dataType = "DATA";

  out_name=Form("AllTrigLumi%s%s%sDiJetMass.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/scratch/xuq7/RpA/TreeAna/%s",out_name.Data()),"RECREATE");  

  for(int itrig=0; itrig < 5; itrig++){ 
    jetpt[itrig]->Write();
    jetptEta[itrig]->Write();
    jetptphi[itrig]->Write();
    rawptJES[itrig]->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[itrig][ieta]->Write(); 
	 rawptJESEtaBin[itrig][ieta]->Write(); 
 	 }
	}
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void AllTriganalyzeRpATreeLumi()
{
  std::cout<<"start working\n";
  hist_class *my_hists = new hist_class();

  cout<<"Analyzing DATA!"<<endl;

  TFile *f = new TFile("/scratch/xuq7/RpA/TreeAna/data_ppReco_akPu3PF_AlljetTrigKurtCombinationWeight_JetPt0noIPupperCut.root");  
  
  TTree *nt = (TTree*)f->Get("nt");

  Double_t jtpt,jteta,jtphi,rawpt;
  Int_t HLT_PAJet20_noJetID_v1;
  Int_t HLT_PAJet40_noJetID_v1;
  Int_t HLT_PAJet60_noJetID_v1;
  Int_t HLT_PAJet80_noJetID_v1;
  Int_t HLT_PAJet100_noJetID_v1;

  Int_t HLT_PAJet20_noJetID_v1_Prescl;
  Int_t HLT_PAJet40_noJetID_v1_Prescl;
  Int_t HLT_PAJet60_noJetID_v1_Prescl;
  Int_t HLT_PAJet80_noJetID_v1_Prescl;
  Int_t HLT_PAJet100_noJetID_v1_Prescl;
  Int_t run;
  Double_t hiHFplusEta4;

  nt->SetBranchAddress("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_noJetID_v1);
  
  nt->SetBranchAddress("HLT_PAJet20_noJetID_v1_Prescl",&HLT_PAJet20_noJetID_v1_Prescl);
  nt->SetBranchAddress("HLT_PAJet40_noJetID_v1_Prescl",&HLT_PAJet40_noJetID_v1_Prescl);
  nt->SetBranchAddress("HLT_PAJet60_noJetID_v1_Prescl",&HLT_PAJet60_noJetID_v1_Prescl);
  nt->SetBranchAddress("HLT_PAJet80_noJetID_v1_Prescl",&HLT_PAJet80_noJetID_v1_Prescl);
  nt->SetBranchAddress("HLT_PAJet100_noJetID_v1_Prescl",&HLT_PAJet100_noJetID_v1_Prescl); 
  nt->SetBranchAddress("run",&run);
  nt->SetBranchAddress("rawpt",&rawpt);
  nt->SetBranchAddress("jtpt",&jtpt);
  nt->SetBranchAddress("jteta",&jteta);
  nt->SetBranchAddress("jtphi",&jtphi);
  nt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);

  TString CorrfileName;
  CorrfileName="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/pPbHijing_EffCorr_VRw.root";
  Corrector corr_Qiao(CorrfileName);
	fcrel3 = TFile::Open("Casym_pPb_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root", "readonly");
   if(fcrel3)  C_rel=(TH1D*)fcrel3->Get("C_asym");

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
//    for(int i=0; i<50000; i++){
    nt->GetEntry(i);

    double HFbin = hiHFplusEta4;

    weight =1.;

    if(i%1000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
 
     weight*=corr_Qiao.getEventWeightHFPlus4bak(HFbin,kTRUE);//efficiency weight
    
    double jet_pt = jtpt;
    double jet_eta = jteta;  
    double jet_phi = jtphi; 
    double raw_pt = rawpt;
	double jetweight=1.;
  //  cout << jet_pt << jet_eta <<endl;
  
    int jetTr2[5];
    double Prescl[5];
    for(int itrig=0; itrig<5; itrig++){
    if(itrig==0) {
      jetTr2[itrig] = HLT_PAJet20_noJetID_v1;
      Prescl[itrig]=HLT_PAJet20_noJetID_v1_Prescl;}
    else if (itrig==1) {
      jetTr2[itrig] = HLT_PAJet40_noJetID_v1; 
      Prescl[itrig]=HLT_PAJet40_noJetID_v1_Prescl;}
    else if (itrig==2) {
      jetTr2[itrig] = HLT_PAJet60_noJetID_v1;
      Prescl[itrig]=HLT_PAJet60_noJetID_v1_Prescl;}
    else if (itrig==3) {
      jetTr2[itrig] = HLT_PAJet80_noJetID_v1;
      Prescl[itrig]=HLT_PAJet80_noJetID_v1_Prescl;}
    else if (itrig==4) {
      jetTr2[itrig] = HLT_PAJet100_noJetID_v1;
      Prescl[itrig]=HLT_PAJet100_noJetID_v1_Prescl;}
    }

    for(int itrig=0; itrig < 5; itrig++){
      int dEtaBin = -1.;
	// jetweight*=C_rel->GetBinContent(C_rel->FindBin(jet_eta));
	
	if(TMath::Abs(jet_eta)<=3){
         if(jetTr2[itrig]==1 && jet_pt>trigRes[itrig] && jet_pt<=trigRes[itrig+1]){
	my_hists->jetptEta[itrig]->Fill(jet_pt*jetweight,jet_eta,weight*Prescl[itrig]);
         my_hists->jetptphi[itrig]->Fill(jet_pt*jetweight,jet_phi,weight*Prescl[itrig]);
	}
	}
	if(rawpt>=22){
      if(TMath::Abs(jet_eta+0.465)<=1. && jetTr2[itrig]==1 ) {
		if( jet_pt>trigRes[itrig] && jet_pt<=trigRes[itrig+1]){
		my_hists->jetpt[itrig]->Fill(jet_pt*jetweight,weight*Prescl[itrig]);
        	my_hists->rawptJES[itrig]->Fill(raw_pt,jet_pt/raw_pt,weight*Prescl[itrig]);
         }
	}

      for(int ieta =0; ieta<netabin;ieta++){
          if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1 && jetTr2[itrig]==1){
	if( jet_pt>trigRes[itrig] && jet_pt<=trigRes[itrig+1]){
	 my_hists->jetptEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,weight*Prescl[itrig]);
	 my_hists->rawptJESEtaBin[itrig][dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,weight*Prescl[itrig]);
		}
	}
	}
    } //loop over itriag
  }//loop over tree

  my_hists->Write();
  cout<<"working done\n";
 
}  
