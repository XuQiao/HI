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
//#include "/scratch/xuq7/RpA/OldAna/Corrections/Correction.C"
using namespace std;
//const double deta[] = {-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,2.2};
const double deta[] = {-2.0, -1.5, -1.0, -0.5, 0.5, 1.0,1.5,2.0};
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
    TH2F * jetEtaphi[5];
    TH2F * jetptchMaxpt[5];      //Added
    TH2F * jetptchMax[5];    //Added
    TH2F * jetptneuMaxpt[5];     //Added
    TH2F * jetptneuMaxr[5];   //Added
    TH2F * jetptneuMax[5];       //Added
    TH2F * jetptphoMaxpt[5]; //Added 
    TH2F * jetptphoMax[5];   //Added
    TH2F * jetptchSum[5];    //Added
    TH2F * jetptneuSum[5];    //Added
    TH2F * jetptphoSum[5];    //Added
    TH2F * jetptchSumpt[5];    //Added
    TH2F * jetptneuSumpt[5];    //Added
    TH2F * jetptphoSumpt[5];    //Added
    TH2F* jetptchN[5];       //Added
    TH2F* jetptneuN[5];      //Added
    TH2F* jetptphoN[5];      //Added
    TH2F * jetptSumSum[5];   //Added
    TH2F * jetptSumSumrawpt[5];   //Added
    TH2F * jetptSumSumpt[5];   //Added
    TH2F * rawptJES[5];
    TH2F * rawptJESEtaBin[5][7];
    TH2F * jetptchMaxptEtaBin[5][7];   //Added
    TH2F * jetptSumSumptEtaBin[5][7];   //Added
    TH2F * jetptneuMaxrEtaBin[5][7];
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
    jetEtaphi[itrig] = new TH2F(Form("%sjetEtaphi_0-100%%",TrigName[itrig].Data()),Form("%sjetEtaphi_0-100%%",TrigName[itrig].Data()),100,-5.,5.,200, -TMath::Pi(), TMath::Pi());
    jetEtaphi[itrig] -> Sumw2(); 

    jetptchMax[itrig] = new TH2F(Form("%sjetptchMax_0-100%%",TrigName[itrig].Data()), Form("jetptchMax_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);       //Added
        jetptchMax[itrig]->Sumw2();      //Added

        jetptchMaxpt[itrig] = new TH2F(Form("%sjetptchMaxpt_0-100%%",TrigName[itrig].Data()), Form("jetptchMaxpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 2);    //Added
    jetptchMaxpt[itrig]->Sumw2();      //Added

        jetptneuMaxpt[itrig] = new TH2F(Form("%sjetptneuMaxpt_0-100%%",TrigName[itrig].Data()), Form("jetptneuMaxpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 2.);    //Added
    jetptneuMaxpt[itrig]->Sumw2();      //Added

        jetptneuMax[itrig] = new TH2F(Form("%sjetptneuMax_0-100%%",TrigName[itrig].Data()), Form("jetptneuMax_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
    jetptneuMax[itrig]->Sumw2();      //Added

        jetptneuMaxr[itrig] = new TH2F(Form("%sjetptneuMaxr_0-100%%",TrigName[itrig].Data()), Form("jetptneuMaxr_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 1000, 0., 5.);    //Added
        jetptneuMaxr[itrig]->Sumw2();      //Added

       jetptphoMaxpt[itrig] = new TH2F(Form("%sjetptphoMaxpt_0-100%%",TrigName[itrig].Data()), Form("jetptphoMaxpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptphoMaxpt[itrig]->Sumw2();      //Added

        jetptphoMax[itrig] = new TH2F(Form("%sjetptphoMax_0-100%%",TrigName[itrig].Data()), Form("jetptphoMax_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoMax[itrig]->Sumw2();      //Added

        jetptchSum[itrig] = new TH2F(Form("%sjetptchSum_0-100%%",TrigName[itrig].Data()), Form("jetptchSum_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptchSum[itrig]->Sumw2();      //Added

        jetptneuSum[itrig] = new TH2F(Form("%sjetptneuSum_0-100%%",TrigName[itrig].Data()), Form("jetptneuSum_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptneuSum[itrig]->Sumw2();      //Added

        jetptphoSum[itrig] = new TH2F(Form("%sjetptphoSum_0-100%%",TrigName[itrig].Data()), Form("jetptphoSum_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoSum[itrig]->Sumw2();      //Added

        jetptchSumpt[itrig] = new TH2F(Form("%sjetptchSumpt_0-100%%",TrigName[itrig].Data()), Form("jetptchSumpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptchSumpt[itrig]->Sumw2();      //Added

        jetptneuSumpt[itrig] = new TH2F(Form("%sjetptneuSumpt_0-100%%",TrigName[itrig].Data()), Form("jetptneuSumpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptneuSumpt[itrig]->Sumw2();      //Added

        jetptphoSumpt[itrig] = new TH2F(Form("%sjetptphoSumpt_0-100%%",TrigName[itrig].Data()), Form("jetptphoSumpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptphoSumpt[itrig]->Sumw2();      //Added

        jetptchN[itrig] = new TH2F(Form("%sjetptchN_0-100%%",TrigName[itrig].Data()), Form("jetptchN_0-100%%",TrigName[itrig].Data()),1000,0., 1000., 100, 0, 100);
        jetptchN[itrig]->Sumw2();

        jetptneuN[itrig] = new TH2F(Form("%sjetptneuN_0-100%%",TrigName[itrig].Data()), Form("jetptneuN_0-100%%",TrigName[itrig].Data()),1000,0., 1000., 100, 0, 100);
        jetptneuN[itrig]->Sumw2();

        jetptphoN[itrig] = new TH2F(Form("%sjetptphoN_0-100%%",TrigName[itrig].Data()), Form("jetptphoN_0-100%%",TrigName[itrig].Data()),1000,0., 1000., 100, 0, 100);
        jetptphoN[itrig]->Sumw2();

        jetptSumSum[itrig] = new TH2F(Form("%sjetptSumSum_0-100%%",TrigName[itrig].Data()), Form("jetptSumSum_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 1000, 0., 1000);    //Added
        jetptSumSum[itrig]->Sumw2();      //Added

        jetptSumSumrawpt[itrig] = new TH2F(Form("%sjetptSumSumrawpt_0-100%%",TrigName[itrig].Data()), Form("jetptSumSumrawpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumrawpt[itrig]->Sumw2();      //Added
        jetptSumSumpt[itrig] = new TH2F(Form("%sjetptSumSumpt_0-100%%",TrigName[itrig].Data()), Form("jetptSumSumpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumpt[itrig]->Sumw2();      //Added

    rawptJES[itrig] = new TH2F(Form("%srawptJES_0-100%%",TrigName[itrig].Data()),Form("%srawptJES_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,200,0.,2.);
    rawptJES[itrig] -> Sumw2(); 
    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[itrig][ieta] = new TH1F(Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[itrig][ieta]->Sumw2();
      jetptchMaxptEtaBin[itrig][ieta] = new TH2F(Form("%sjetptchMaxptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptchMaxptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200, 0.,2.);
      jetptchMaxptEtaBin[itrig][ieta]->Sumw2();
      jetptSumSumptEtaBin[itrig][ieta] = new TH2F(Form("%sjetptSumSumptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptSumSumptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,2000, 0.,10);
      jetptSumSumptEtaBin[itrig][ieta]->Sumw2();
      jetptneuMaxrEtaBin[itrig][ieta] = new TH2F(Form("%sjetptneuMaxrEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptneuMaxrEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,1000, 0. ,5.);
      jetptneuMaxrEtaBin[itrig][ieta]->Sumw2();
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
    delete jetEtaphi[itrig];
            delete jetptchMax[itrig];        //Added
        delete  jetptchMaxpt[itrig];       //Added
        delete  jetptneuMaxpt[itrig]; //Added
        delete    jetptneuMaxr[itrig];   //Added
        delete    jetptneuMax[itrig];   //Added
        delete    jetptphoMaxpt[itrig]; //Added 
        delete    jetptphoMax[itrig];   //Added
        delete    jetptchSum[itrig];
        delete    jetptneuSum[itrig];
        delete    jetptphoSum[itrig];
        delete    jetptchSumpt[itrig];
        delete    jetptneuSumpt[itrig];
        delete    jetptphoSumpt[itrig];
        delete    jetptchN[itrig];
        delete    jetptneuN[itrig];
        delete    jetptphoN[itrig];
        delete    jetptSumSum[itrig];   //Added
        delete    jetptSumSumrawpt[itrig];   //Added
        delete    jetptSumSumpt[itrig];   //Added
    delete rawptJES[itrig];
    for(int ieta=0; ieta < netabin; ieta++){
	 delete jetptEtaBin[itrig][ieta];
	 delete jetptSumSumptEtaBin[itrig][ieta];
	 delete jetptchMaxptEtaBin[itrig][ieta];
	 delete jetptneuMaxrEtaBin[itrig][ieta];
	 delete rawptJESEtaBin[itrig][ieta];
	}
    }
}

void hist_class::Write()
{
  TString dataType;
  TString out_name;
  dataType = "DATA";
  out_name=Form("%s%s%sUniv_woCom.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/TreeAna/%s",out_name.Data()),"RECREATE");  

  for(int itrig=0; itrig < 5; itrig++){ 
    jetpt[itrig]->Write();
    jetptEta[itrig]->Write();
    jetptphi[itrig]->Write();
    jetEtaphi[itrig]->Write();
        jetptchMax[itrig]->Write();      //Added
jetptchMaxpt[itrig]->Write();       //Added
jetptneuMaxpt[itrig]->Write(); //Added
jetptneuMaxr[itrig]->Write();   //Added
jetptneuMax[itrig]->Write();   //Added
jetptphoMaxpt[itrig]->Write(); //Added 
jetptphoMax[itrig]->Write();   //Added
jetptchSum[itrig]->Write();
jetptneuSum[itrig]->Write();
jetptphoSum[itrig]->Write();
jetptchSumpt[itrig]->Write();
jetptneuSumpt[itrig]->Write();
jetptphoSumpt[itrig]->Write();
jetptchN[itrig]->Write();
jetptneuN[itrig]->Write();
jetptphoN[itrig]->Write();
jetptSumSum[itrig]->Write();   //Added
jetptSumSumrawpt[itrig]->Write();   //Added
jetptSumSumpt[itrig]->Write();   //Added
    rawptJES[itrig]->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[itrig][ieta]->Write(); 
	 jetptSumSumptEtaBin[itrig][ieta]->Write();
	 jetptchMaxptEtaBin[itrig][ieta]->Write();
	 jetptneuMaxrEtaBin[itrig][ieta]->Write();
	 rawptJESEtaBin[itrig][ieta]->Write(); 
 	 }
	}
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void DataanalyzerRpATree()
{
  std::cout<<"start working\n";
  hist_class *my_hists = new hist_class();

  cout<<"Analyzing DATA!"<<endl;

  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/PPbdata_ppReco_akPu3PF_AlljetTrigKurtTrCombFile_merged.root");
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
	Int_t chargedN, neutralN, photonN;
Double_t chargedSum, photonSum, neutralSum, chargedMax, neutralMax,photonMax, eSum, muSum;
  Double_t weight;

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
    nt->SetBranchAddress("chargedN",&chargedN);
    nt->SetBranchAddress("photonN",&photonN);
    nt->SetBranchAddress("neutralN",&neutralN);
    nt->SetBranchAddress("chargedMax",&chargedMax);
    nt->SetBranchAddress("photonMax",&photonMax);
    nt->SetBranchAddress("neutralMax",&neutralMax);
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
//    for(int i=0; i<50000; i++){
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

    for(int itrig=0; itrig < 5; itrig++){
      int dEtaBin = -1.;
	jetweight*=(fUE->Eval(jet_pt))*C_rel->GetBinContent(C_rel->FindBin(jet_eta));        
   if(raw_pt<22) continue;
   if(TMath::Abs(jet_eta+0.465)>3.) continue;
//   if((chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt>1.01) continue;

	 if(jetTr2[itrig]==1 && jet_pt>trigRes[itrig] && jet_pt<=trigRes[itrig+1]){
//        if(jetTr2[itrig]==1){
	double	evtweight = weight * Prescl[itrig];
	my_hists->jetptEta[itrig]->Fill(jet_pt*jetweight,jet_eta,evtweight);
         my_hists->jetptphi[itrig]->Fill(jet_pt*jetweight,jet_phi,evtweight);
         my_hists->jetEtaphi[itrig]->Fill(jet_eta,jet_phi,evtweight);

      if(TMath::Abs(jet_eta+0.465)<=1.) {
	  my_hists->jetptchMax[itrig]->Fill(jet_pt*jetweight, chargedMax, evtweight);  //Added
    	  my_hists->jetptneuMax[itrig]->Fill(jet_pt*jetweight, neutralMax, evtweight);
    	  my_hists->jetptphoMax[itrig]->Fill(jet_pt*jetweight, photonMax, evtweight);
    	  my_hists->jetptchMaxpt[itrig]->Fill(jet_pt*jetweight, chargedMax/jet_pt, evtweight);
    	  my_hists->jetptneuMaxpt[itrig]->Fill(jet_pt*jetweight, neutralMax/jet_pt, evtweight);
    	  my_hists->jetptphoMaxpt[itrig]->Fill(jet_pt*jetweight, photonMax/jet_pt, evtweight);
          my_hists->jetptneuMaxr[itrig]->Fill(jet_pt*jetweight, neutralMax/TMath::Max(chargedSum,neutralSum) ,evtweight);
          my_hists->jetptchSum[itrig]->Fill(jet_pt*jetweight, chargedSum, evtweight);
          my_hists->jetptneuSum[itrig]->Fill(jet_pt*jetweight, neutralSum, evtweight);
          my_hists->jetptphoSum[itrig]->Fill(jet_pt*jetweight, photonSum, evtweight);
          my_hists->jetptchSumpt[itrig]->Fill(jet_pt*jetweight, chargedSum/jet_pt, evtweight);
          my_hists->jetptneuSumpt[itrig]->Fill(jet_pt*jetweight, neutralSum/jet_pt, evtweight);
          my_hists->jetptphoSumpt[itrig]->Fill(jet_pt*jetweight, photonSum/jet_pt, evtweight);

          my_hists->jetptchN[itrig]->Fill(jet_pt*jetweight, chargedN, evtweight);
          my_hists->jetptneuN[itrig]->Fill(jet_pt*jetweight, neutralN, evtweight);
          my_hists->jetptphoN[itrig]->Fill(jet_pt*jetweight, photonN, evtweight);

          my_hists->jetptSumSum[itrig]->Fill(jet_pt*jetweight, chargedSum+photonSum+neutralSum+muSum+eSum,evtweight);
          my_hists->jetptSumSumpt[itrig]->Fill(jet_pt*jetweight, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,evtweight);
          my_hists->jetptSumSumrawpt[itrig]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,evtweight);

		my_hists->jetpt[itrig]->Fill(jet_pt*jetweight,evtweight);
        	my_hists->rawptJES[itrig]->Fill(raw_pt,jet_pt/raw_pt,evtweight);
	}

      for(int ieta =0; ieta<netabin;ieta++){
          if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 my_hists->jetptEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,evtweight);
	 my_hists->jetptSumSumptEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,evtweight);
	 my_hists->jetptchMaxptEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,chargedMax/jet_pt,evtweight);
	 my_hists->jetptneuMaxrEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,neutralMax/TMath::Max(chargedSum,neutralSum),evtweight);
	 my_hists->rawptJESEtaBin[itrig][dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,evtweight);
	}

	} //trig selection and pt range selection
    } //loop over itrig
  }//loop over tree

  my_hists->Write();
  cout<<"working done\n";
 
}  


