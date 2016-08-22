// this code is modified to use Data only
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
const TString JetIDName[]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","eSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN","PPcut","PPcutTight","PPcutTighter","SumSumpt1"};
const int nJetID = sizeof(JetIDName)/sizeof(TString);


TString algo="akPu3PF";//"akPu3PF"
TString coll = "PbP";
    
void DataJetIDskimTree(){
    TH2F* jetptjetid[nJetID];
    TH1F* jetpt;
    TH2F* Spikeetaphi;
    TH2F* Spikeetapt;
        jetpt = new TH1F(Form("jetpt"),Form("jetpt"),1000,0.,1000.);
        jetpt -> Sumw2(); 
        Spikeetaphi = new TH2F(Form("Spikeetaphi"),Form("Spikeetaphi"),200,-5.,5.,200, -TMath::Pi(), TMath::Pi());
        Spikeetaphi -> Sumw2();
        Spikeetapt = new TH2F(Form("Spikeetapt"),Form("Spikeetapt"),200,-5.,5.,1000,0.,1000.);
        Spikeetapt -> Sumw2();
   for(int ijetid=0;ijetid<nJetID;ijetid++){
        if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr")){
        jetptjetid[ijetid] = new TH2F(Form("jetpt%s",JetIDName[ijetid].Data()), Form("jetpt%s",JetIDName[ijetid].Data()), 1000, 0., 1000., 200, 0., 2.);      //Added
        }
        else if(JetIDName[ijetid].Contains("N") || JetIDName[ijetid].Contains("PP")){
        jetptjetid[ijetid] = new TH2F(Form("jetpt%s",JetIDName[ijetid].Data()), Form("jetpt%s",JetIDName[ijetid].Data()), 1000, 0., 1000., 100, 0., 100.);    //Added
        }
        else{
        jetptjetid[ijetid] = new TH2F(Form("jetpt%s",JetIDName[ijetid].Data()), Form("jetpt%s",JetIDName[ijetid].Data()), 1000, 0., 1000., 3000, 0., 300.);   //Added
        }
        jetptjetid[ijetid]->Sumw2();      //Added
        }

    TH2F * jetptjetidEtaBin[netabin][nJetID];
    TH1F * jetptEtaBin[netabin];
		
    for(int ieta=0; ieta < netabin; ieta++){
        jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
        jetptEtaBin[ieta]->Sumw2();
       for(int ijetid=0;ijetid<nJetID;ijetid++){
        if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr")){
        jetptjetidEtaBin[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptjetidEtaBin[ieta][ijetid]->Sumw2();      //Added
        }
        if(JetIDName[ijetid].Contains("PP") || JetIDName[ijetid].Contains("N")){
        jetptjetidEtaBin[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 100, 0., 100.);    //Added
        jetptjetidEtaBin[ieta][ijetid]->Sumw2();      //Added
        }
        }
    }

  std::cout<<"start working\n";

  cout<<"Analyzing DATA!"<<endl;

  TFile *f = new TFile(Form("/cms/store/user/qixu/jetRpA/skimTree/DATA%s%sskimfile.root",coll.Data(),algo.Data()));
     
     TFile *fcrel3 = NULL ;
     TH1D *C_rel= NULL ;

        fcrel3 = TFile::Open(Form("/cms/store/user/qixu/jetRpA/RpA/NewMC/Corrections/CasymYaxian_%s_double_hcalbins_algo_%s_pt100_140_jet50_alphahigh_20_phicut250.root",coll.Data(),algo.Data()), "readonly");
     if(fcrel3)  C_rel=(TH1D*)fcrel3->Get("C_asym");

   TF1 * fUE = new TF1("fUE","1-[0]/pow(x,[1])",20,600);
  if(algo=="ak3PF")  fUE->SetParameters(0.8648,0.8167);
  else if (algo=="akPu3PF") fUE->SetParameters(0.3015,0.8913);
  else fUE->SetParameters(0.,0.);

  TTree *nt = (TTree*)f->Get("nt");

  Float_t jtpt[1000],jteta[1000],jtphi[1000],rawpt[1000];
        Int_t t_chargedN[100], t_neutralN[100], t_photonN[100];
Float_t t_chargedSum[100], t_photonSum[100], t_neutralSum[100], t_chargedMax[100], t_neutralMax[100], t_photonMax[100], t_eSum[100], t_muSum[100];
  Int_t run,nref,hiBin;
  Int_t HLT_PAJet20_noJetID_v1;
  Int_t HLT_PAJet40_noJetID_v1;
  Int_t HLT_PAJet60_noJetID_v1;
  Int_t HLT_PAJet80_noJetID_v1;
  Int_t HLT_PAJet100_noJetID_v1;

  Int_t pPAcollisionEventSelectionPA;
  Int_t pprimaryVertexFilter;
  Int_t pHBHENoiseFilter;

  Float_t weight;
  Float_t vz,pt;

  nt->SetBranchAddress("run",&run);
  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("pt",&pt);
  nt->SetBranchAddress("hiBin",&hiBin);
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("rawpt",rawpt);
  nt->SetBranchAddress("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_noJetID_v1);
  nt->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
  if(coll=="PPb")
  nt->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
  if(coll=="PbP")
  nt->SetBranchAddress("pprimaryvertexFilter",&pprimaryVertexFilter);
  nt->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  nt->SetBranchAddress("jtpt",jtpt);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("jtphi",jtphi);    
    nt->SetBranchAddress("chargedN",t_chargedN);
    nt->SetBranchAddress("photonN",t_photonN);
    nt->SetBranchAddress("neutralN",t_neutralN);
    nt->SetBranchAddress("chargedMax",t_chargedMax);
    nt->SetBranchAddress("photonMax",t_photonMax);
    nt->SetBranchAddress("neutralMax",t_neutralMax);
    nt->SetBranchAddress("chargedSum",t_chargedSum);
    nt->SetBranchAddress("photonSum",t_photonSum);
    nt->SetBranchAddress("neutralSum",t_neutralSum);
    nt->SetBranchAddress("muSum",t_muSum);
    nt->SetBranchAddress("eSum",t_eSum);
    nt->SetBranchAddress("weight",&weight);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
    nt->GetEntry(i);
if(coll == "PPb" && run>211256) continue;//only for pPb runs at moment
if(coll == "PbP" && run<=211256) continue;//only for pPb runs at moment
if((TMath::Abs(vz)>15) || (!pPAcollisionEventSelectionPA) || (!pprimaryVertexFilter)|| (!pHBHENoiseFilter)) continue;
 if(!HLT_PAJet20_noJetID_v1 && !HLT_PAJet40_noJetID_v1 && !HLT_PAJet60_noJetID_v1 && !HLT_PAJet80_noJetID_v1 && !HLT_PAJet100_noJetID_v1 ) continue;

    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
for(int j4i = 0; j4i < nref; j4i++){
    double jet_pt = jtpt[j4i];
    double jet_eta = jteta[j4i];  
    double jet_phi = jtphi[j4i];  
    double raw_pt = rawpt[j4i];
    double jetweight = 1.;
 //       jetweight*=(fUE->Eval(jet_pt))*C_rel->GetBinContent(C_rel->FindBin(jet_eta));
        jetweight*=C_rel->GetBinContent(C_rel->FindBin(jet_eta));

		double chargedMax = t_chargedMax[j4i];
		double chargedSum = t_chargedSum[j4i];
		double neutralMax = t_neutralMax[j4i];
		double neutralSum = t_neutralSum[j4i];
		double photonMax = t_photonMax[j4i];
		double photonSum = t_photonSum[j4i];
		int chargedN = t_chargedN[j4i];
		int neutralN = t_neutralN[j4i];
		int photonN = t_photonN[j4i];
		double muSum = t_muSum[j4i];
		double eSum = t_eSum[j4i];
  //      if((chargedN == 0 || chargedSum == 0) && TMath::Abs(jet_eta)< 2.4) continue;  // jet id selection
	double PPTighter0 = (double)(neutralSum/jet_pt < 0.8 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 1.0 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4));
	double PPTighter1 = (double)(neutralSum/jet_pt < 0.8 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 1.0 && chargedSum/jet_pt < 0.95 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
	double PPTighter2 = (double)(neutralSum/jet_pt < 0.9 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 1.0 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
	double PPTighter3 = (double)(neutralSum/jet_pt < 0.9 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 1.0 && chargedSum/jet_pt < 0.95 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
	double PPTighter4 = (double)(neutralSum/jet_pt < 0.9 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 0.95 && chargedSum/jet_pt < 0.95 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
	double PPTighter5 = (double)(neutralSum/jet_pt < 0.8 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 0.95 && chargedSum/jet_pt < 0.95 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
        double PPTighter = PPTighter0*TMath::Power(2,5)+PPTighter1*TMath::Power(2,4)+PPTighter2*TMath::Power(2,3)+PPTighter3*TMath::Power(2,2)+PPTighter4*TMath::Power(2,1)+PPTighter5*TMath::Power(2,0)+0.5;
	double jetidv[nJetID]={chargedMax,chargedSum,neutralMax,neutralSum,photonMax,photonSum,chargedMax/jet_pt,chargedSum/jet_pt,neutralMax/jet_pt,neutralSum/jet_pt,photonMax/jet_pt,photonSum/jet_pt,eSum/jet_pt,(chargedSum+neutralSum+photonSum+muSum+eSum)/jet_pt,(chargedSum+neutralSum+photonSum+muSum+eSum)/raw_pt,neutralMax/TMath::Max(chargedSum,neutralSum),(double)chargedN,(double)neutralN,(double)photonN,(double)(neutralSum/jet_pt<1.0 && eSum/jet_pt<1.0 && photonSum/jet_pt<1.0 && ((chargedSum>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) ), (double)(neutralSum/jet_pt<0.9 && eSum/jet_pt<1.0 && photonSum/jet_pt<0.9 && ((chargedSum>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) ),PPTighter,(chargedSum+neutralSum+muSum+eSum)/jet_pt};
       if(raw_pt<22 || fabs(jet_eta)>5) continue;
       if(jet_pt>4*pt) continue;
    //if(neutralSum/jet_pt > 0.8 && chargedSum/jet_pt < 0.1) {
    if(neutralSum/jet_pt > 0.8) {
        Spikeetaphi->Fill(jet_eta,jet_phi, weight);
        Spikeetapt->Fill(jet_eta,jet_pt, weight);
    }
       
        jet_pt = jet_pt*jetweight;
        
      int dEtaBin = -1;
	if(coll=="PPb")	jet_eta=jet_eta+0.465;
	if(coll=="PbP")	jet_eta=jet_eta-0.465;


     if(TMath::Abs(jet_eta)<=1.) {
         jetpt->Fill(jet_pt,weight);
	for(int ijetid=0;ijetid<nJetID;ijetid++){
            jetptjetid[ijetid]->Fill(jet_pt, jetidv[ijetid],weight);     //Added
	}	
    }
      for(int ieta =0; ieta<netabin;ieta++){
          if(jet_eta>deta[ieta]&&jet_eta<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
         jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
	for(int ijetid=0;ijetid<nJetID;ijetid++){
        if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr") || JetIDName[ijetid].Contains("PP") || JetIDName[ijetid].Contains("N")){
        jetptjetidEtaBin[dEtaBin][ijetid]->Fill(jet_pt,jetidv[ijetid],weight);
	}
        }
      }
      } //loop over jet
	
}	//loop over tree

  TString dataType;
  TString out_name;
  dataType = "DATA";
  out_name=Form("%s%s%sskimJetID.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/tmp/xuq7/%s",out_name.Data()),"RECREATE");  

    jetpt->Write();
    Spikeetaphi->Write();
    Spikeetapt->Write();
    for(int ijetid=0;ijetid<nJetID;ijetid++){
        jetptjetid[ijetid]->Write();
    }
    for(int ieta=0; ieta<netabin; ieta++){
        for(int ijetid=0;ijetid<nJetID;ijetid++){
          if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr") || JetIDName[ijetid].Contains("PP") || JetIDName[ijetid].Contains("N")){
         jetptjetidEtaBin[ieta][ijetid]->Write();
        }
        }
         jetptEtaBin[ieta]->Write();
     }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;

  cout<<"working done\n";
 
} 
