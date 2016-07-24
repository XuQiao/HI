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

void DataskimTreePres(){
    TH1F * jetpt[5];
    TH2F * jetptEta[5];
    TH2F * jetptphi[5];
    TH1F * jetptEtaBin[5][7];
  
for(int itrig=0; itrig< 5; itrig++){
    jetpt[itrig] = new TH1F(Form("%sjetpt_0-100%%",TrigName[itrig].Data()),Form("%sjetpt_0-100%%",TrigName[itrig].Data()),1000,0.,1000.);
    jetpt[itrig] -> Sumw2(); 
    jetptEta[itrig] = new TH2F(Form("%sjetptEta_0-100%%",TrigName[itrig].Data()),Form("%sjetptEta_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,100,-5.,5.);
    jetptEta[itrig] -> Sumw2(); 
    jetptphi[itrig] = new TH2F(Form("%sjetptphi_0-100%%",TrigName[itrig].Data()),Form("%sjetptphi_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi[itrig] -> Sumw2(); 
    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[itrig][ieta] = new TH1F(Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[itrig][ieta]->Sumw2();
      }
    }


  cout<<"start working\n";

  cout<<"Analyzing DATA!"<<endl;
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

  Float_t jtpt[100],jteta[100],jtphi[100],rawpt[100];
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
  Int_t run,nref;
  Int_t pPAcollisionEventSelectionPA;
  Int_t pHBHENoiseFilter;
  Float_t hiHFplusEta4;
  Float_t weight, evtweight, vz;
  Float_t t_chargedMax[100], t_chargedSum[100], t_photonSum[100], t_neutralSum[100], t_muSum[100], t_eSum[100];

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
  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
  nt->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("rawpt",rawpt);
  nt->SetBranchAddress("jtpt",jtpt);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("jtphi",jtphi);
  nt->SetBranchAddress("chargedMax",t_chargedMax);
  nt->SetBranchAddress("chargedSum",t_chargedSum);
  nt->SetBranchAddress("photonSum",t_photonSum);
  nt->SetBranchAddress("neutralSum",t_neutralSum);
  nt->SetBranchAddress("muSum",t_muSum);
  nt->SetBranchAddress("eSum",t_eSum);

  nt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
//nt->SetBranchAddress("evtweight",&evtweight);

//  TString CorrfileName;
//  CorrfileName="/scratch/xuq7/RpA/OldAna/Corrections/pPbHijing_EffCorr_VRw.root";
//  Corrector corr_Qiao(CorrfileName);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
    nt->GetEntry(i);
    if(i%100000==1)    cout<<"analyzing "<< i <<" th event"<<endl;
if(run>211256) continue;//only for pPb runs at moment
//if(run<=211256) continue;//only for Pbp runs at moment
if((TMath::Abs(vz)>15) || (!pPAcollisionEventSelectionPA) || (!pHBHENoiseFilter)) continue;
if(!HLT_PAJet20_noJetID_v1 && !HLT_PAJet40_noJetID_v1 && !HLT_PAJet60_noJetID_v1 && !HLT_PAJet80_noJetID_v1 && !HLT_PAJet100_noJetID_v1 ) continue;

    double HFbin = hiHFplusEta4;

 
//    weight*=corr_Qiao.getEventWeightHFPlus4bak(HFbin,kTRUE);//efficiency evtweight
	   weight =1.;
    
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
//        if(jetTr2[itrig]==1)
  if(jetTr2[itrig]==1)
{
	evtweight = weight * Prescl[itrig];
  for(int j4i=0; j4i<nref;j4i++){
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
    double jetweight=1.;
if(chargedMax/jet_pt<=0.05 || chargedMax/jet_pt>=1.00) continue; //Jet ID cut
//if((chargedSum+photonSum+neutralSum+eSum+muSum)/jet_pt>1.01) continue; //Jet ID cut
//	jetevtweight*=(fUE->Eval(jet_pt))*C_rel->GetBinContent(C_rel->FindBin(jet_eta));        

//   if(jet_pt>trigRes[itrig] && jet_pt<=trigRes[itrig+1]){
	   
//	if(raw_pt<22) continue;
//   	if(TMath::Abs(jet_eta+0.465)>3.) continue;

//	jetptEta[itrig]->Fill(jet_pt*jetweight,jet_eta,evtweight);
//         jetptphi[itrig]->Fill(jet_pt*jetweight,jet_phi,evtweight);
      
	int dEtaBin = -1.;
      if(TMath::Abs(jet_eta+0.465)<=1.) {
	jetpt[itrig]->Fill(jet_pt*jetweight,evtweight);
	}

/*      for(int ieta =0; ieta<netabin;ieta++){
          if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 jetptEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,evtweight);
	}*/
//      	} //pt range selection
      } //loop over jet
    }//trig selection
       } //loop over itrig
  }//loop over tree


  TString dataType;
  TString out_name;
  dataType = "DATA";
  out_name=Form("%s%s%sJetID1_skimUniv_woCom.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/%s",out_name.Data()),"RECREATE");

  for(int itrig=0; itrig < 5; itrig++){
    jetpt[itrig]->Write();
    jetptEta[itrig]->Write();
    jetptphi[itrig]->Write();
    for(int ieta=0; ieta<netabin; ieta++){
         jetptEtaBin[itrig][ieta]->Write();
         }
        }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;

  cout<<"working done\n";
}  


