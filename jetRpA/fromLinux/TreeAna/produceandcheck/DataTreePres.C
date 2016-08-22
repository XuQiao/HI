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
///#include "/scratch/xuq7/RpA/OldAna/Corrections/Correction.C"
using namespace std;
//const double deta[] = {-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,2.2};
const double deta[] = {-2.0, -1.5, -1.0, -0.5, 0.5, 1.0,1.5,2.0};
const int netabin = sizeof(deta)/sizeof(double)-1;
TString TrigName[] = {"Jet20","Jet40","Jet60","Jet80","Jet100"};
const int trigRes[]={20,55,75,95,120,999};

TString algo="akPu3PF";//"akPu3PF"
TString coll = "PPb";
double weight;

void DataTreePres(){
    TH1F * jetpt[5];
 /*   TH2F * jetptEta[5];
    TH2F * jetptphi[5];
    TH2F * jetEtaphi[5];
    TH2F * rawptJES[5];
    TH2F * rawptJESEtaBin[5][7];
    TH1F * jetptEtaBin[5][7];
  */
for(int itrig=0; itrig< 5; itrig++){
    jetpt[itrig] = new TH1F(Form("%sjetpt_0-100%%",TrigName[itrig].Data()),Form("%sjetpt_0-100%%",TrigName[itrig].Data()),1000,0.,1000.);
    jetpt[itrig] -> Sumw2(); 
/*    jetptEta[itrig] = new TH2F(Form("%sjetptEta_0-100%%",TrigName[itrig].Data()),Form("%sjetptEta_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,100,-5.,5.);
    jetptEta[itrig] -> Sumw2(); 
    jetptphi[itrig] = new TH2F(Form("%sjetptphi_0-100%%",TrigName[itrig].Data()),Form("%sjetptphi_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi[itrig] -> Sumw2(); 
    jetEtaphi[itrig] = new TH2F(Form("%sjetEtaphi_0-100%%",TrigName[itrig].Data()),Form("%sjetEtaphi_0-100%%",TrigName[itrig].Data()),100,-5.,5.,200, -TMath::Pi(), TMath::Pi());
    jetEtaphi[itrig] -> Sumw2(); 

    rawptJES[itrig] = new TH2F(Form("%srawptJES_0-100%%",TrigName[itrig].Data()),Form("%srawptJES_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,200,0.,2.);
    rawptJES[itrig] -> Sumw2(); 
    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[itrig][ieta] = new TH1F(Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[itrig][ieta]->Sumw2();
      rawptJESEtaBin[itrig][ieta] = new TH2F(Form("%srawptJESEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%srawptJESEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[itrig][ieta]->Sumw2();
      }*/
    }

  std::cout<<"start working\n";

  cout<<"Analyzing DATA!"<<endl;
if(coll=="PPb")
  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/PPbdata_ppReco_akPu3PF_AlljetTrigKurtTrCombFile_merged.root");
else if(coll=="PbP")
  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/PbPdataKurt_ppReco_akPu3PF_AlljetTrigKurtTrCombFile0_6815_JetPt0noIPupperCut.root");
      TFile *fcrel3 = NULL ;
     TH1D *C_rel= NULL ;

     if(coll=="PPb")
        fcrel3 = TFile::Open(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/Corrections/Casym_pPb_double_hcalbins_algo_%s_pt100_140_jet80_alphahigh_20_phicut250.root", algo.Data()), "readonly");
     if(coll=="PbP")
        fcrel3 = TFile::Open(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/Corrections/Casym_Pbp_double_hcalbins_algo_%s_pt100_140_jet80_alphahigh_20_phicut250.root", algo.Data()), "readonly");

     if(fcrel3)  C_rel=(TH1D*)fcrel3->Get("C_asym");

   TF1 * fUE = new TF1("fUE","1-[0]/pow(x,[1])",20,600);
  if(algo=="ak3PF")  fUE->SetParameters(0.8648,0.8167);
  else if (algo=="akPu3PF") fUE->SetParameters(0.3015,0.8913);
  else fUE->SetParameters(0.,0.);

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
  Double_t weight;
  Double_t chargedMax,chargedSum,photonSum,neutralSum,eSum,muSum;

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
  nt->SetBranchAddress("jtpt",&jtpt);
  nt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
  nt->SetBranchAddress("chargedMax",&chargedMax);
  nt->SetBranchAddress("chargedSum",&chargedSum);
  nt->SetBranchAddress("photonSum",&photonSum);
  nt->SetBranchAddress("neutralSum",&neutralSum);
  nt->SetBranchAddress("muSum",&muSum);
  nt->SetBranchAddress("eSum",&eSum);
  nt->SetBranchAddress("weight",&weight);

//  TString CorrfileName;
//  CorrfileName="/scratch/xuq7/RpA/OldAna/Corrections/pPbHijing_EffCorr_VRw.root";
//  Corrector corr_Qiao(CorrfileName);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
    nt->GetEntry(i);

    double HFbin = hiHFplusEta4;

   weight =1.;

    if(i%100000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
 
//    weight*=corr_Qiao.getEventWeightHFPlus4bak(HFbin,kTRUE);//efficiency weight
    
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
//if(chargedMax/jet_pt<=0.05 || chargedMax/jet_pt>=1.00) continue; //Jet ID cut
if((chargedSum+photonSum+neutralSum+eSum+muSum)/jet_pt>1.01) continue; //Jet ID cut

    for(int itrig=0; itrig < 5; itrig++){
      int dEtaBin = -1.;
//	jetweight*=(fUE->Eval(jet_pt))*C_rel->GetBinContent(C_rel->FindBin(jet_eta));        
//   if(raw_pt<22) continue;
//   if(TMath::Abs(jet_eta+0.465)>3.) continue;
//   if((chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt>1.01) continue;

//	 if(jetTr2[itrig]==1 && jet_pt>trigRes[itrig] && jet_pt<=trigRes[itrig+1]){
         if(jetTr2[itrig]==1){
	double	evtweight = weight * Prescl[itrig];
//	jetptEta[itrig]->Fill(jet_pt*jetweight,jet_eta,evtweight);
//         jetptphi[itrig]->Fill(jet_pt*jetweight,jet_phi,evtweight);
//         jetEtaphi[itrig]->Fill(jet_eta,jet_phi,evtweight);

      if(TMath::Abs(jet_eta+0.465)<=1.) {
		jetpt[itrig]->Fill(jet_pt*jetweight,evtweight);
      //  	rawptJES[itrig]->Fill(raw_pt,jet_pt/raw_pt,evtweight);
	}
/*
      for(int ieta =0; ieta<netabin;ieta++){
          if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 jetptEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,evtweight);
	 rawptJESEtaBin[itrig][dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,evtweight);
	}
*/
	} //trig selection and pt range selection
    } //loop over itrig
  }//loop over tree
 
 TString dataType;
  TString out_name;
  dataType = "DATA";
  out_name=Form("%s%s%sJetID2_Univ_woCom.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/%s",out_name.Data()),"RECREATE");

  for(int itrig=0; itrig < 5; itrig++){
    jetpt[itrig]->Write();
/*    jetptEta[itrig]->Write();
    jetptphi[itrig]->Write();
    jetEtaphi[itrig]->Write();
    rawptJES[itrig]->Write();
    for(int ieta=0; ieta<netabin; ieta++){
         jetptEtaBin[itrig][ieta]->Write();
         rawptJESEtaBin[itrig][ieta]->Write();
         }*/
        }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;

  cout<<"working done\n";
 
}  

