// this code is modified to use MC only
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
//const double deta[] = {-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,2.2};
const double deta[] = {-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0};
const int netabin = sizeof(deta)/sizeof(double)-1;

TString algo="akPu3PF";//"akPu3PF"
TString coll = "PPb";
TString trigName[5] = {"Jet20","Jet40","Jet60","Jet80","Jet100"} ;   

void DataskimTree_temp(){
    TH1F * jetpt[5];
    TH2F * jetptEta;
    TH2F * jetptEta_woRes;
    TH2F * jetptphi;
    TH2F * jetEtaphi;
    TH2F * rawptJES;
    TH1F * Vz;
    TH1F * Cent;

    TH2F * rawptJESEtaBin[7];
    TH1F * jetptEtaBin[7];
for(int itrig=0; itrig<5 ;itrig++){
    jetpt[itrig] = new TH1F(Form("%sjetpt",trigName[itrig].Data()),Form("%sjetpt",trigName[itrig].Data()),1000,0.,1000.);
    jetpt[itrig] -> Sumw2(); 
}
    jetptEta = new TH2F(Form("jetptEta"),Form("jetptEta"),1000,0.,1000.,100,-5.,5.);
    jetptEta -> Sumw2(); 
    jetptEta_woRes = new TH2F(Form("jetptEta_woRes"),Form("jetptEta_woRes"),1000,0.,1000.,100,-5.,5.);
    jetptEta_woRes -> Sumw2(); 
    jetptphi = new TH2F(Form("jetptphi"),Form("jetptphi"),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi -> Sumw2(); 
    jetEtaphi = new TH2F(Form("jetEtaphi"),Form("jetEtaphi"),100,-5.,5.,200, -TMath::Pi(), TMath::Pi());
    jetEtaphi -> Sumw2(); 

    rawptJES = new TH2F(Form("rawptJES"),Form("rawptJES"),1000,0.,1000.,200,0.,2.);
    rawptJES -> Sumw2(); 
    Vz = new TH1F(Form("Vz"),Form("Vz"),600,-30.,30.);
    Vz-> Sumw2(); 
    Cent = new TH1F(Form("Cent"),Form("Cent"),100,0.,100.);
    Cent-> Sumw2(); 

    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[ieta]->Sumw2();
      rawptJESEtaBin[ieta] = new TH2F(Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[ieta]->Sumw2();
      }

  std::cout<<"start working\n";

  cout<<"Analyzing Data!"<<endl;

  TFile *f = new TFile(Form("/home/xuq/Documents/HI/RpA/output/skimTree/DATA%s%sskimfile.root",coll.Data(),algo.Data()));
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

  Float_t jtpt[1000],jteta[1000],jtphi[1000],rawpt[1000];
  Float_t t_chargedMax[1000],t_chargedSum[1000],t_photonSum[1000],t_neutralSum[1000],t_eSum[1000],t_muSum[1000];
  Int_t run,nref,hiBin;
  Int_t HLT_PAJet20_noJetID_v1;
  Int_t HLT_PAJet40_noJetID_v1;
  Int_t HLT_PAJet60_noJetID_v1;
  Int_t HLT_PAJet80_noJetID_v1;
  Int_t HLT_PAJet100_noJetID_v1;

  Int_t pPAcollisionEventSelectionPA;
  Int_t pHBHENoiseFilter;

  Float_t weight;
  Float_t vz;

  nt->SetBranchAddress("run",&run);
  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("hiBin",&hiBin);
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("rawpt",rawpt);
  nt->SetBranchAddress("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_noJetID_v1);

  nt->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
  nt->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  nt->SetBranchAddress("jtpt",jtpt);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("jtphi",jtphi);
    nt->SetBranchAddress("chargedMax",&t_chargedMax);
  nt->SetBranchAddress("chargedSum",&t_chargedSum);
  nt->SetBranchAddress("photonSum",&t_photonSum);
  nt->SetBranchAddress("neutralSum",&t_neutralSum);
  nt->SetBranchAddress("muSum",&t_muSum);
  nt->SetBranchAddress("eSum",&t_eSum);

  nt->SetBranchAddress("weight",&weight);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
    nt->GetEntry(i);
if(run>211256) continue;//only for pPb runs at moment
if((TMath::Abs(vz)>15) || (!pPAcollisionEventSelectionPA) || (!pHBHENoiseFilter)) continue;
 if(!HLT_PAJet20_noJetID_v1 && !HLT_PAJet40_noJetID_v1 && !HLT_PAJet60_noJetID_v1 && !HLT_PAJet80_noJetID_v1 && !HLT_PAJet100_noJetID_v1 ) continue;

	Vz->Fill(vz);
	Cent->Fill(hiBin);

      int jetTr2[5];
    for(int itrig=0; itrig<5; itrig++){
    if(itrig==0) 
      jetTr2[itrig] = HLT_PAJet20_noJetID_v1;
    else if (itrig==1) 
      jetTr2[itrig] = HLT_PAJet40_noJetID_v1;
    else if (itrig==2) 
      jetTr2[itrig] = HLT_PAJet60_noJetID_v1;
    else if (itrig==3) 
      jetTr2[itrig] = HLT_PAJet80_noJetID_v1;
    else if (itrig==4) 
      jetTr2[itrig] = HLT_PAJet100_noJetID_v1;
    }
    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
for(int j4i = 0; j4i < nref; j4i++){
    double jet_pt = jtpt[j4i];
    double jet_eta = jteta[j4i];  
    double jet_phi = jtphi[j4i]; 
    double raw_pt = rawpt[j4i];
       double chargedMax = t_chargedMax[j4i];
    double chargedSum = t_chargedSum[j4i];
    double photonSum = t_photonSum[j4i];
    double neutralSum = t_neutralSum[j4i];
    double muSum = t_muSum[j4i];
    double eSum = t_eSum[j4i];

    double jetweight = 1.;
//	jetweight*=(fUE->Eval(jet_pt))*C_rel->GetBinContent(C_rel->FindBin(jet_eta));        
	//if(chargedMax/jet_pt<=0.05 || chargedMax/jet_pt>=1.00) continue; //Jet ID cut
if((chargedSum+photonSum+neutralSum+eSum+muSum)/jet_pt>1.01) continue; //Jet ID cut
 
       if(raw_pt<20) continue;
 /*       if(TMath::Abs(jet_eta)<=3){
                jetptEta->Fill(jet_pt*jetweight,jet_eta,weight);
                jetptEta_woRes->Fill(jet_pt,jet_eta,weight);
                jetptphi->Fill(jet_pt*jetweight,jet_phi,weight);
                jetEtaphi->Fill(jet_eta,jet_phi,weight);
        } 
*/
      int dEtaBin = -1.;
	if(coll=="PPb")	jet_eta = jet_eta+0.465;
	else if(coll=="PbP")	jet_eta = jet_eta-0.465;
     if(TMath::Abs(jet_eta)<=1.) {
    for(int itrig=0; itrig < 5; itrig++){
  if(jetTr2[itrig]==1)
		jetpt[itrig]->Fill(jet_pt*jetweight,weight);
//        	rawptJES->Fill(raw_pt,jet_pt/raw_pt,weight);
	}
}
/*
      for(int ieta =0; ieta<netabin;ieta++){
          if(jet_eta>deta[ieta]&&jet_eta<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 jetptEtaBin[dEtaBin]->Fill(jet_pt*jetweight,weight);
	 rawptJESEtaBin[dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,weight);
	}*/
      } //loop over jet
}	//loop over tree

  TString dataType;
  TString out_name;
  dataType = "DATA";
  out_name=Form("%s%s%sJetID2skimUniv_Trig.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/%s",out_name.Data()),"RECREATE");  

    for(int itrig=0; itrig < 5; itrig++){
    jetpt[itrig]->Write();}
    jetptEta->Write();
    jetptEta_woRes->Write();
    jetptphi->Write();
    jetEtaphi->Write();
    rawptJES->Write();
    Vz->Write();
    Cent->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[ieta]->Write(); 
	 rawptJESEtaBin[ieta]->Write(); 
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;

  cout<<"working done\n";
 
} 
