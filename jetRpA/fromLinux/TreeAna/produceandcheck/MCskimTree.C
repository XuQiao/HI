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
TString wtype = "";    

void MCskimTree(){
    TH1F * refjetpt;
    TH1F * genjetpt;
    TH1F * jetpt;
    TH2F * jetptEta;
    TH2F * refptEta;
    TH2F * jetptphi;
    TH2F * jetEtaphi;
    TH2F * rawptJES;
    TH2F * refptJES;
    TH1F * Vz;
    TH1F * VzW;
    TH1F * Cent;
    TH1F * CentvzW;
    TH1F * CentW;
    TH1F * Pthat;
    TH1F * PthatAfXw;

    TH2F * rawptJESEtaBin[7];
    TH1F * refjetptEtaBin[7];
    TH1F * genjetptEtaBin[7];
    TH2F * refptJESEtaBin[7];
    TH1F * jetptEtaBin[7];

    refjetpt = new TH1F(Form("refjetpt"),Form("refjetpt"),1000,0.,1000.);
    refjetpt -> Sumw2(); 
    genjetpt = new TH1F(Form("genjetpt"),Form("genjetpt"),1000,0.,1000.);
    genjetpt -> Sumw2(); 
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
    Vz = new TH1F(Form("Vz"),Form("Vz"),600,-30.,30.);
    Vz-> Sumw2(); 
    VzW = new TH1F(Form("VzW"),Form("VzW"),600,-30.,30.);
    VzW-> Sumw2(); 
    Cent = new TH1F(Form("Cent"),Form("Cent"),100,0.,100.);
    Cent-> Sumw2(); 
    CentvzW = new TH1F(Form("CentvzW"),Form("CentvzW"),100,0.,100.);
    CentvzW-> Sumw2(); 
    CentW = new TH1F(Form("CentW"),Form("CentW"),100,0.,100.);
    CentW-> Sumw2(); 
    Pthat = new TH1F(Form("Pthat"),Form("Pthat"),1000,0.,1000.);
    Pthat-> Sumw2(); 
    PthatAfXw = new TH1F(Form("PthatAfXw"),Form("PthatAfXw"),1000,0.,1000.);
    PthatAfXw-> Sumw2(); 

    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[ieta]->Sumw2();
      refjetptEtaBin[ieta] = new TH1F(Form("refjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("refjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      refjetptEtaBin[ieta]->Sumw2();
      genjetptEtaBin[ieta] = new TH1F(Form("genjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("genjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      genjetptEtaBin[ieta]->Sumw2();
      refptJESEtaBin[ieta] = new TH2F(Form("refptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("refptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      refptJESEtaBin[ieta]->Sumw2();
      rawptJESEtaBin[ieta] = new TH2F(Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[ieta]->Sumw2();
      }

  std::cout<<"start working\n";

  cout<<"Analyzing MC!"<<endl;

if(coll=="PPb")
  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/MCPPbakPu3PFskimfile0_8.root");  
else if(coll=="PP")
  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/MCPPakPu3PFskimfile0_9.root");  
else if(coll=="PbP")
  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/MCPbPakPu3PFskimfile0_8.root");  
  
  TTree *nt = (TTree*)f->Get("nt");

  Float_t jtpt[100],jteta[100],geneta[100],jtphi[100],rawpt[100],refpt[100],refeta[100],genpt[100];
Int_t subid[100];
Int_t nref,ngen,hiBin;
  Float_t xSecWeight,weight,vzWeight,centWeight;
  Float_t vz,pthat;

  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("hiBin",&hiBin);
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("rawpt",rawpt);
  nt->SetBranchAddress("jtpt",jtpt);
  nt->SetBranchAddress("refpt",refpt);
  nt->SetBranchAddress("refeta",refeta);
  nt->SetBranchAddress("ngen",&ngen);
  nt->SetBranchAddress("genpt",genpt);
  nt->SetBranchAddress("geneta",geneta);
  nt->SetBranchAddress("subid",subid);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("jtphi",jtphi);
    nt->SetBranchAddress("pthat",&pthat);
    nt->SetBranchAddress("xSecWeight",&xSecWeight);
    nt->SetBranchAddress("vzWeight",&vzWeight);
    nt->SetBranchAddress("centWeight",&centWeight);
    nt->SetBranchAddress("weight",&weight);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
    nt->GetEntry(i);
	if(TMath::Abs(vz)>15) continue;
        
	Vz->Fill(vz,xSecWeight);
	VzW->Fill(vz,xSecWeight*vzWeight);
	Cent->Fill(hiBin,xSecWeight);
	CentvzW->Fill(hiBin,xSecWeight*vzWeight);
	CentW->Fill(hiBin,weight);
	Pthat->Fill(pthat);
	PthatAfXw->Fill(pthat,xSecWeight);

    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
for(int j4i = 0; j4i < nref; j4i++){
    double jet_pt = jtpt[j4i];
    double ref_pt = refpt[j4i];
    double ref_eta = refeta[j4i];
    double jet_eta = jteta[j4i];  
    double jet_phi = jtphi[j4i]; 
    double raw_pt = rawpt[j4i];
    int sub_id = subid[j4i];

	if(wtype=="xSec") weight = xSecWeight;

	refptEta->Fill(ref_pt,ref_eta,xSecWeight);	
	if(coll=="PPb") ref_eta = ref_eta+0.465;
        else if(coll=="PbP")    ref_eta = ref_eta-0.465;

	if(TMath::Abs(ref_eta) <= 1.) {
        refjetpt->Fill(ref_pt, xSecWeight);
       }
    int dEtaBin = -1. ;
        for(Int_t ieta = 0 ; ieta <netabin; ieta++){
               if(ref_eta>deta[ieta]&&ref_eta<=deta[ieta+1]) dEtaBin = ieta ;
         } //assign the eta bin for jets
      if(dEtaBin!=-1){
        refjetptEtaBin[dEtaBin]->Fill(ref_pt, xSecWeight);
       }
	        if(raw_pt<22) continue;
        if(TMath::Abs(jet_eta)<=3 && sub_id==0){
                jetptEta->Fill(jet_pt,jet_eta,weight);
                jetptphi->Fill(jet_pt,jet_phi,weight);
                jetEtaphi->Fill(jet_eta,jet_phi,weight);
	}
        if(coll=="PPb") jet_eta = jet_eta+0.465;
        else if(coll=="PbP")    jet_eta = jet_eta-0.465;
	else;

      int dEtaBin = -1.;
     if(TMath::Abs(jet_eta)<=1.) {
		jetpt->Fill(jet_pt,weight);
		refptJES->Fill(ref_pt,jet_pt/ref_pt,xSecWeight);
        	rawptJES->Fill(raw_pt,jet_pt/raw_pt,xSecWeight);
	}

      for(int ieta =0; ieta<netabin;ieta++){
          if(jet_eta>deta[ieta]&&jet_eta<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
	 rawptJESEtaBin[dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,xSecWeight);
	 refptJESEtaBin[dEtaBin]->Fill(ref_pt,jet_pt/ref_pt,xSecWeight);
	}
      } //loop over jet
	
	for(int j5i = 0; j5i < ngen ; j5i++) {
    double gen_pt=genpt[j5i];
    double gen_eta=geneta[j5i];
        if(coll=="PPb") gen_eta = gen_eta+0.465;
        else if(coll=="PbP")    gen_eta = gen_eta-0.465;
        if(TMath::Abs(gen_eta)<=1.) {
          genjetpt->Fill(gen_pt, xSecWeight);
       }
     int    dEtaBin=-1;
        for(Int_t ieta = 0 ; ieta <netabin; ieta++){
               if(gen_eta>deta[ieta]&&gen_eta<=deta[ieta+1]) dEtaBin = ieta ;
         } //assign the eta bin for jets
      if(dEtaBin!=-1){
        genjetptEtaBin[dEtaBin]->Fill(gen_pt, xSecWeight);
        }
	} //end of gen jet loop

}	//loop over tree

  TString dataType;
  TString out_name;
  dataType = "MC";
  out_name=Form("%s%s%s%sskimUniv.root",dataType.Data(),coll.Data(),algo.Data(),wtype.Data());

  TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/%s",out_name.Data()),"RECREATE");  

	
    jetpt->Write();
    jetptEta->Write();
    genjetpt->Write();
    refjetpt->Write();
    refptEta->Write();
    refptJES->Write();
    jetptphi->Write();
    jetEtaphi->Write();
    rawptJES->Write();
    Vz->Write();
    VzW->Write();
    Cent->Write();
    CentvzW->Write();
    CentW->Write();
    Pthat->Write();
    PthatAfXw->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[ieta]->Write(); 
	 genjetptEtaBin[ieta]->Write(); 
	 refjetptEtaBin[ieta]->Write(); 
	 refptJESEtaBin[ieta]->Write(); 
	 rawptJESEtaBin[ieta]->Write(); 
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;

  cout<<"working done\n";
}
