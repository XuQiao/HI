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
    
void MCTree(){
    TH1F * refjetpt;
    TH1F * jetpt;
    TH2F * jetptEta;
    TH2F * refptEta;
    TH2F * jetptphi;
    TH2F * jetEtaphi;
    TH2F * rawptJES;
    TH2F * refptJES;

    TH2F * rawptJESEtaBin[7];
    TH1F * refjetptEtaBin[7];
    TH2F * refptJESEtaBin[7];
    TH1F * jetptEtaBin[7];

    refjetpt = new TH1F(Form("refjetpt"),Form("refjetpt"),1000,0.,1000.);
    refjetpt -> Sumw2(); 
    jetpt = new TH1F(Form("jetpt"),Form("jetpt"),1000,0.,1000.);
    jetpt -> Sumw2(); 
    jetptEta = new TH2F(Form("jetptEta"),Form("jetptEta"),1000,0.,1000.,100,-5.,5.);
    jetptEta -> Sumw2(); 
    refptEta = new TH2F(Form("refptEta"),Form("refptEta"),1000,0.,1000.,100,-5.,5.);
    refptEta -> Sumw2(); 
    jetptphi = new TH2F(Form("jetptphi"),Form("jetptphi"),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi -> Sumw2(); 
    jetEtaphi = new TH2F(Form("jetEtaphi"),Form("jetEtaphi"),100,-5.,5.,200, -TMath::Pi(), TMath::Pi());
    jetEtaphi -> Sumw2(); 

    rawptJES = new TH2F(Form("rawptJES"),Form("rawptJES"),1000,0.,1000.,200,0.,2.);
    rawptJES -> Sumw2(); 
    refptJES = new TH2F(Form("refptJES"),Form("refptJES"),1000,0.,1000.,200,0.,2.);
    refptJES -> Sumw2(); 

    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[ieta]->Sumw2();
      refjetptEtaBin[ieta] = new TH1F(Form("refjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("refjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      refjetptEtaBin[ieta]->Sumw2();
      refptJESEtaBin[ieta] = new TH2F(Form("refptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("refptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      refptJESEtaBin[ieta]->Sumw2();
      rawptJESEtaBin[ieta] = new TH2F(Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[ieta]->Sumw2();
      }

  std::cout<<"start working\n";

  cout<<"Analyzing MC!"<<endl;

if(coll=="PPb")
  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/PPbMCKurtForest_ppReco_akPu3PF_QCDjetQiaoIDCut_JetPt0noIPupperCut.root");  
else if(coll=="PbP")
  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/PbPMCKurtForestVzWeight_ppReco_akPu3PF_QCDjetQiaoIDCut_JetPt0noIPupperCut.root");  
  
  TTree *nt = (TTree*)f->Get("nt");

Double_t jtpt,jteta,geneta,jtphi,rawpt,refpt,refeta,genpt;
Int_t subid;
Int_t nref,hiBin;
  Double_t weight;
  Double_t vz,pthat;

  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("hiBin",&hiBin);
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("rawpt",&rawpt);
  nt->SetBranchAddress("jtpt",&jtpt);
  nt->SetBranchAddress("refpt",&refpt);
  nt->SetBranchAddress("refeta",&refeta);
  nt->SetBranchAddress("subid",&subid);
  nt->SetBranchAddress("jteta",&jteta);
  nt->SetBranchAddress("jtphi",&jtphi);
    nt->SetBranchAddress("weight",&weight);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
    nt->GetEntry(i);
	if(TMath::Abs(vz)>15) continue;
        
    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
    double jet_pt = jtpt;
    double ref_pt = refpt;
    double ref_eta = refeta;
    double jet_eta = jteta;  
    double jet_phi = jtphi; 
    double raw_pt = rawpt;
    int sub_id = subid;
	refptEta->Fill(ref_pt,ref_eta,weight);

   if(coll=="PPb") ref_eta = ref_eta+0.465;
        else if(coll=="PbP")    ref_eta = ref_eta-0.465;

	if(TMath::Abs(ref_eta) <= 1.) {
        refjetpt->Fill(ref_pt, weight);
       }
    int dEtaBin = -1. ;
        for(Int_t ieta = 0 ; ieta <netabin; ieta++){
               if(ref_eta>deta[ieta]&&ref_eta<=deta[ieta+1]) dEtaBin = ieta ;
         } //assign the eta bin for jets
      if(dEtaBin!=-1){
        refjetptEtaBin[dEtaBin]->Fill(ref_pt, weight);
       }

	        if(raw_pt<22) continue;
        if(TMath::Abs(jet_eta)<=3){
                jetptEta->Fill(jet_pt,jet_eta,weight);
                jetptphi->Fill(jet_pt,jet_phi,weight);
                jetEtaphi->Fill(jet_eta,jet_phi,weight);
        } 

      int dEtaBin = -1.;
	   if(coll=="PPb") jet_eta = jet_eta+0.465;
        else if(coll=="PbP")    jet_eta = jet_eta-0.465;
     if(TMath::Abs(jet_eta)<=1.) {
		jetpt->Fill(jet_pt,weight);
		refptJES->Fill(ref_pt,jet_pt/ref_pt,weight);
        	rawptJES->Fill(raw_pt,jet_pt/raw_pt,weight);
	}
      for(int ieta =0; ieta<netabin;ieta++){
          if(jet_eta>deta[ieta]&&jet_eta<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
	 rawptJESEtaBin[dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,weight);
	 refptJESEtaBin[dEtaBin]->Fill(ref_pt,jet_pt/ref_pt,weight);
	}
      } //loop over jet
	
  TString dataType;
  TString out_name;
  dataType = "MC";
  out_name=Form("%s%s%sUniv.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/%s",out_name.Data()),"RECREATE");  

	
    refjetpt->Write();
    jetpt->Write();
    jetptEta->Write();
    refptEta->Write();
    jetptphi->Write();
    jetEtaphi->Write();
    rawptJES->Write();
    refptJES->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[ieta]->Write(); 
	 refjetptEtaBin[ieta]->Write(); 
	 refptJESEtaBin[ieta]->Write(); 
	 rawptJESEtaBin[ieta]->Write(); 
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;

  cout<<"working done\n";
 
} 
