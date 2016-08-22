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
const double deta[] = {-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,2.2};
const int netabin = sizeof(deta)/sizeof(double)-1;
const TString JetIDName[]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN","PPcut"};
const int nJetID = sizeof(JetIDName)/sizeof(TString);

TString algo="akPu3PF";//"akPu3PF"
TString coll = "PPb";

class hist_class{
public:
    hist_class();
    void Write();

    TH1F * jetpt;
    TH2F * jetptEta;
    TH2F * jetptphi;
    TH2F * jetEtaphi;
    TH2F * jetptjetid[nJetID];
    TH2F * jetptjetid_real[nJetID];
    TH2F * jetptjetid_fake[nJetID];
    TH1F * Pthat;
    TH1F * PthatAfXw;
    TH2F * rawptJES;

    TH2F * rawptJESEtaBin[netabin];
    TH2F * jetptjetidEtaBin[netabin][nJetID];
    TH1F * jetptEtaBin[netabin];
};

hist_class::hist_class()
{
    jetpt = new TH1F(Form("jetpt"),Form("jetpt"),1000,0.,1000.);
    jetpt -> Sumw2(); 
    jetptEta = new TH2F(Form("jetptEta"),Form("jetptEta"),1000,0.,1000.,100,-5.,5.);
    jetptEta -> Sumw2(); 
    jetptphi = new TH2F(Form("jetptphi"),Form("jetptphi"),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi -> Sumw2(); 
    jetEtaphi = new TH2F(Form("jetEtaphi"),Form("jetEtaphi"),100,-5.,5.,200, -TMath::Pi(), TMath::Pi());
    jetEtaphi -> Sumw2(); 

       Pthat = new TH1F(Form("Pthat"),Form("Pthat"),1000,0.,1000.);
        Pthat->Sumw2();
        PthatAfXw = new TH1F(Form("PthatAfXw"),Form("PthatAfXw"),1000,0.,1000.);
        PthatAfXw->Sumw2();

//----------------------------------------------------------------------------------------------

        for(int ijetid=0;ijetid<nJetID;ijetid++){
        if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr")){
        jetptjetid[ijetid] = new TH2F(Form("jetpt%s",JetIDName[ijetid].Data()), Form("jetpt%s",JetIDName[ijetid].Data()), 1000, 0., 1000., 200, 0., 2.);      //Added
        jetptjetid_real[ijetid] = new TH2F(Form("jetpt%s_real",JetIDName[ijetid].Data()), Form("jetpt%s_real",JetIDName[ijetid].Data()), 1000, 0., 1000., 200, 0., 2.);       //Added
        jetptjetid_fake[ijetid] = new TH2F(Form("jetpt%s_fake",JetIDName[ijetid].Data()), Form("jetpt%s_fake",JetIDName[ijetid].Data()), 1000, 0., 1000., 200, 0., 2.);       //Added
        }
        else if(JetIDName[ijetid].Contains("N")){
        jetptjetid[ijetid] = new TH2F(Form("jetpt%s",JetIDName[ijetid].Data()), Form("jetpt%s",JetIDName[ijetid].Data()), 1000, 0., 1000., 100, 0., 100.);    //Added
        jetptjetid_real[ijetid] = new TH2F(Form("jetpt%s_real",JetIDName[ijetid].Data()), Form("jetpt%s_real",JetIDName[ijetid].Data()), 1000, 0., 1000., 100, 0., 100.);     //Added
        jetptjetid_fake[ijetid] = new TH2F(Form("jetpt%s_fake",JetIDName[ijetid].Data()), Form("jetpt%s_fake",JetIDName[ijetid].Data()), 1000, 0., 1000., 100, 0., 100.);     //Added
        }
        else{
        jetptjetid[ijetid] = new TH2F(Form("jetpt%s",JetIDName[ijetid].Data()), Form("jetpt%s",JetIDName[ijetid].Data()), 1000, 0., 1000., 3000, 0., 300.);   //Added
        jetptjetid_real[ijetid] = new TH2F(Form("jetpt%s_real",JetIDName[ijetid].Data()), Form("jetpt%s_real",JetIDName[ijetid].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptjetid_fake[ijetid] = new TH2F(Form("jetpt%s_fake",JetIDName[ijetid].Data()), Form("jetpt%s_fake",JetIDName[ijetid].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        }
        jetptjetid[ijetid]->Sumw2();      //Added
        jetptjetid_real[ijetid]->Sumw2(); //Added
        jetptjetid_fake[ijetid]->Sumw2(); //Added
        }

//----------------------------------------------------------------

    rawptJES = new TH2F(Form("rawptJES"),Form("rawptJES"),1000,0.,1000.,200,0.,2.);
    rawptJES -> Sumw2(); 
    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[ieta]->Sumw2();
	for(int ijetid=0;ijetid<nJetID;ijetid++){
        if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr")){
        jetptjetidEtaBin[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptjetidEtaBin[ieta][ijetid]->Sumw2();      //Added
        }
        }
      rawptJESEtaBin[ieta] = new TH2F(Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[ieta]->Sumw2();
      }
}

void hist_class::Write()
{
  TString dataType;
  TString out_name;
  dataType = "MC";
  out_name=Form("%s%s%s.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/scratch/xuq7/RpA/NewMC/%s",out_name.Data()),"RECREATE");  
    jetpt->Write();
    jetptEta->Write();
    jetptphi->Write();
    jetEtaphi->Write();
	 for(int ijetid=0;ijetid<nJetID;ijetid++){
        jetptjetid[ijetid]->Write();
        jetptjetid_real[ijetid]->Write();
        jetptjetid_fake[ijetid]->Write();
        }
    Pthat->Write();
    PthatAfXw->Write();
    rawptJES->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[ieta]->Write(); 
	for(int ijetid=0;ijetid<nJetID;ijetid++){
        if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr")){
                 jetptjetidEtaBin[ieta][ijetid]->Write();
        }
        }
	 rawptJESEtaBin[ieta]->Write(); 
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void PPbanalyzeRpATree()
{
  std::cout<<"start working\n";
  hist_class *my_hists = new hist_class();

  cout<<"Analyzing MC!"<<endl;

  TFile *f = new TFile("/home/maoy/working/pA/RpA/TreeProduce/CMSSW_5_3_20/src/PPbMCOfficialForestVzWeight_ppReco_akPu3PF_QCDjetTrigJECv8_JetPt0noIPupperCut.root");
  
  TTree *nt = (TTree*)f->Get("nt");

  Double_t jtpt,jteta,jtphi,rawpt,refpt;
	Int_t chargedN, neutralN, photonN,subid;
Double_t chargedSum, photonSum, neutralSum, chargedMax, neutralMax,photonMax, eSum, muSum;
  Double_t weight;

  nt->SetBranchAddress("rawpt",&rawpt);
  nt->SetBranchAddress("jtpt",&jtpt);
  nt->SetBranchAddress("refpt",&refpt);
  nt->SetBranchAddress("subid",&subid);
  nt->SetBranchAddress("jteta",&jteta);
  nt->SetBranchAddress("jtphi",&jtphi);
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

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
//    for(int i=0; i<50000; i++){
    nt->GetEntry(i);

    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
 
    
    double jet_pt = jtpt;
    double ref_pt = refpt;
    double jet_eta = jteta;  
    double jet_phi = jtphi; 
    double raw_pt = rawpt;
  
      int dEtaBin = -1.;

	if(subid==0 && fabs(jet_eta)<2 ){
		my_hists->jetptEta->Fill(jet_pt,jet_eta,weight);
		my_hists->jetptphi->Fill(jet_pt,jet_phi,weight);
		my_hists->jetEtaphi->Fill(jet_eta,jet_phi,weight);
	}
	double jetidv[nJetID]={chargedMax,chargedSum,neutralMax,neutralSum,photonMax,photonSum,chargedMax/jet_pt,chargedSum/jet_pt,neutralMax/jet_pt,neutralSum/jet_pt,photonMax/jet_pt,photonSum/jet_pt,(chargedSum+neutralSum+photonSum+muSum+eSum)/jet_pt,(chargedSum+neutralSum+photonSum+muSum+eSum)/raw_pt,neutralMax/TMath::Max(chargedSum,neutralSum),(double)chargedN,(double)neutralN,(double)photonN};
     if(TMath::Abs(jet_eta+0.465)<=1.) {
	for(int ijetid=0;ijetid<nJetID;ijetid++){
	my_hists->jetptjetid[ijetid]->Fill(jet_pt, jetidv[ijetid], weight);	//Added
	if(ref_pt<0){	//Added
	  my_hists->jetptjetid_fake[ijetid]->Fill(jet_pt, jetidv[ijetid], weight);	//Added
	}
	else{
	  my_hists->jetptjetid_real[ijetid]->Fill(jet_pt, jetidv[ijetid], weight);	//Added
	}	
	}	
		my_hists->jetpt->Fill(jet_pt,weight);
        	my_hists->rawptJES->Fill(raw_pt,jet_pt/raw_pt,weight);
	}

      for(int ieta =0; ieta<netabin;ieta++){
          if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 my_hists->jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
	for(int ijetid=0;ijetid<nJetID;ijetid++){
        if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr")){
        my_hists->jetptjetidEtaBin[dEtaBin][ijetid]->Fill(jet_pt,jetidv[ijetid],weight);
        }
        }
	 my_hists->rawptJESEtaBin[dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,weight);
	}
  }//loop over tree

  my_hists->Write();
  cout<<"working done\n";
 
}  
