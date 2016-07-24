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

TString algo="akPu3PF";//"akPu3PF"
TString coll = "PPb";
    
void DataJetIDskimTree(){
    TH2F * jetptchMaxpt;      //Added
    TH2F * jetptchMax;    //Added
    TH2F * jetptneuMaxpt;     //Added
    TH2F * jetptneuMaxr;   //Added
    TH2F * jetptneuMax;       //Added
    TH2F * jetptphoMaxpt; //Added 
    TH2F * jetptphoMax;   //Added
    TH2F * jetptchSum;    //Added
    TH2F * jetptneuSum;    //Added
    TH2F * jetptphoSum;    //Added
    TH2F * jetptchSumpt;    //Added
    TH2F * jetptneuSumpt;    //Added
    TH2F * jetptphoSumpt;    //Added
    TH2F * jetptchN;       //Added
    TH2F * jetptneuN;      //Added
    TH2F * jetptphoN;      //Added
    TH2F * jetptSumSum;   //Added
    TH2F * jetptSumSumrawpt;   //Added
    TH2F * jetptSumSumpt;   //Added

TH2F * jetptchMaxEtaBin[7];	//Added
TH2F * jetptchMaxptEtaBin[7];	//Added
TH2F * jetptneuMaxptEtaBin[7];	//Added
TH2F * jetptneuMaxrEtaBin[7];	//Added
TH2F * jetptphoMaxptEtaBin[7];	//Added
TH2F * jetptchSumptEtaBin[7];	//Added
TH2F * jetptneuSumptEtaBin[7];	//Added
TH2F * jetptphoSumptEtaBin[7];	//Added
TH2F * jetptSumSumptEtaBin[7];	//Added
TH2F * jetptSumSumrawptEtaBin[7];		//Added

    jetptchMax = new TH2F(Form("jetptchMaxEtaBin-10_10"), Form("jetptchMaxEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);       //Added
        jetptchMax->Sumw2();      //Added

        jetptchMaxpt = new TH2F(Form("jetptchMaxptEtaBin-10_10"), Form("jetptchMaxptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
    jetptchMaxpt->Sumw2();      //Added

        jetptneuMaxpt = new TH2F(Form("jetptneuMaxptEtaBin-10_10"), Form("jetptneuMaxptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
    jetptneuMaxpt->Sumw2();      //Added

        jetptneuMax = new TH2F(Form("jetptneuMaxEtaBin-10_10"), Form("jetptneuMaxEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
    jetptneuMax->Sumw2();      //Added

        jetptneuMaxr = new TH2F(Form("jetptneuMaxrEtaBin-10_10"), Form("jetptneuMaxrEtaBin-10_10"), 1000, 0., 1000., 1000, 0., 5.);    //Added
        jetptneuMaxr->Sumw2();      //Added

       jetptphoMaxpt = new TH2F(Form("jetptphoMaxptEtaBin-10_10"), Form("jetptphoMaxptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptphoMaxpt->Sumw2();      //Added

        jetptphoMax = new TH2F(Form("jetptphoMaxEtaBin-10_10"), Form("jetptphoMaxEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoMax->Sumw2();      //Added

        jetptchSum = new TH2F(Form("jetptchSumEtaBin-10_10"), Form("jetptchSumEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptchSum->Sumw2();      //Added

        jetptneuSum = new TH2F(Form("jetptneuSumEtaBin-10_10"), Form("jetptneuSumEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptneuSum->Sumw2();      //Added

        jetptphoSum = new TH2F(Form("jetptphoSumEtaBin-10_10"), Form("jetptphoSumEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoSum->Sumw2();      //Added

        jetptchSumpt = new TH2F(Form("jetptchSumptEtaBin-10_10"), Form("jetptchSumptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptchSumpt->Sumw2();      //Added

        jetptneuSumpt = new TH2F(Form("jetptneuSumptEtaBin-10_10"), Form("jetptneuSumptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptneuSumpt->Sumw2();      //Added

        jetptphoSumpt = new TH2F(Form("jetptphoSumptEtaBin-10_10"), Form("jetptphoSumptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptphoSumpt->Sumw2();      //Added

        jetptchN = new TH2F(Form("jetptchNEtaBin-10_10"), Form("jetptchNEtaBin-10_10"),1000,0., 1000., 100, 0, 100);
        jetptchN->Sumw2();

        jetptneuN = new TH2F(Form("jetptneuNEtaBin-10_10"), Form("jetptneuNEtaBin-10_10"),1000,0., 1000., 100, 0, 100);
        jetptneuN->Sumw2();

        jetptphoN = new TH2F(Form("jetptphoNEtaBin-10_10"), Form("jetptphoNEtaBin-10_10"),1000,0., 1000., 100, 0, 100);
        jetptphoN->Sumw2();

        jetptSumSum = new TH2F(Form("jetptSumSumEtaBin-10_10"), Form("jetptSumSumEtaBin-10_10"), 1000, 0., 1000., 1000, 0., 1000);    //Added
        jetptSumSum->Sumw2();      //Added

        jetptSumSumrawpt = new TH2F(Form("jetptSumSumrawptEtaBin-10_10"), Form("jetptSumSumrawptEtaBin-10_10"), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumrawpt->Sumw2();      //Added
        jetptSumSumpt = new TH2F(Form("jetptSumSumptEtaBin-10_10"), Form("jetptSumSumptEtaBin-10_10"), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumpt->Sumw2();      //Added
		
    for(int ieta=0; ieta < netabin; ieta++){

	jetptchMaxEtaBin[ieta] = new TH2F(Form("jetptchMaxEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptchMaxEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 3000,0.,300.);    //Added
        jetptchMaxEtaBin[ieta]->Sumw2();      //Added
	jetptchMaxptEtaBin[ieta] = new TH2F(Form("jetptchMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptchMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptchMaxptEtaBin[ieta]->Sumw2();      //Added
        jetptphoMaxptEtaBin[ieta] = new TH2F(Form("jetptphoMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptphoMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptphoMaxptEtaBin[ieta]->Sumw2();      //Added
        jetptneuMaxptEtaBin[ieta] = new TH2F(Form("jetptneuMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptneuMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptneuMaxptEtaBin[ieta]->Sumw2();      //Added
        jetptchSumptEtaBin[ieta] = new TH2F(Form("jetptchSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptchSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptchSumptEtaBin[ieta]->Sumw2();      //Added
        jetptneuSumptEtaBin[ieta] = new TH2F(Form("jetptneuSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptneuSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptneuSumptEtaBin[ieta]->Sumw2();      //Added
        jetptphoSumptEtaBin[ieta] = new TH2F(Form("jetptphoSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptphoSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptphoSumptEtaBin[ieta]->Sumw2();      //Added
        jetptSumSumptEtaBin[ieta] = new TH2F(Form("jetptSumSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptSumSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 2000, 0., 10.);    //Added
        jetptSumSumptEtaBin[ieta]->Sumw2();        
	jetptSumSumrawptEtaBin[ieta] = new TH2F(Form("jetptSumSumrawptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptSumSumrawptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 2000, 0., 10.);    //Added
        jetptSumSumrawptEtaBin[ieta]->Sumw2();
	jetptneuMaxrEtaBin[ieta] = new TH2F(Form("jetptneuMaxrEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptneuMaxrEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 1000, 0., 5.);    //Added
        jetptneuMaxrEtaBin[ieta]->Sumw2();

      }

  std::cout<<"start working\n";

  cout<<"Analyzing DATA!"<<endl;

if(coll=="PPb")
  TFile *f = new TFile(Form("/home/xuq/Documents/HI/RpA/output/skimTree/DATA%s%sskimfile.root",coll.Data(),algo.Data()));
if(coll=="PbP")
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
        Int_t t_chargedN[100], t_neutralN[100], t_photonN[100];
Float_t t_chargedSum[100], t_photonSum[100], t_neutralSum[100], t_chargedMax[100], t_neutralMax[100], t_photonMax[100], t_eSum[100], t_muSum[100];
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
if(run>211256) continue;//only for pPb runs at moment
if((TMath::Abs(vz)>15) || (!pPAcollisionEventSelectionPA) || (!pHBHENoiseFilter)) continue;
 if(!HLT_PAJet20_noJetID_v1 && !HLT_PAJet40_noJetID_v1 && !HLT_PAJet60_noJetID_v1 && !HLT_PAJet80_noJetID_v1 && !HLT_PAJet100_noJetID_v1 ) continue;

    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
for(int j4i = 0; j4i < nref; j4i++){
    double jet_pt = jtpt[j4i];
    double jet_eta = jteta[j4i];  
    double raw_pt = rawpt[j4i];
    double jetweight = 1.;
        jetweight*=(fUE->Eval(jet_pt))*C_rel->GetBinContent(C_rel->FindBin(jet_eta));

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
	 if(raw_pt<22) continue;
      int dEtaBin = -1;
	if(coll=="PPb")	jet_eta=jet_eta+0.465;
	if(coll=="PbP")	jet_eta=jet_eta-0.465;
     if(TMath::Abs(jet_eta)<=1.) {
	  jetptchMax->Fill(jet_pt*jetweight, chargedMax, weight);  //Added
    	  jetptneuMax->Fill(jet_pt*jetweight, neutralMax, weight);
    	  jetptphoMax->Fill(jet_pt*jetweight, photonMax, weight);
    	  jetptchMaxpt->Fill(jet_pt*jetweight, chargedMax/jet_pt, weight);
    	  jetptneuMaxpt->Fill(jet_pt*jetweight, neutralMax/jet_pt, weight);
    	  jetptphoMaxpt->Fill(jet_pt*jetweight, photonMax/jet_pt, weight);
          jetptneuMaxr->Fill(jet_pt*jetweight, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          jetptchSum->Fill(jet_pt*jetweight, chargedSum, weight);
          jetptneuSum->Fill(jet_pt*jetweight, neutralSum, weight);
          jetptphoSum->Fill(jet_pt*jetweight, photonSum, weight);
          jetptchSumpt->Fill(jet_pt*jetweight, chargedSum/jet_pt, weight);
          jetptneuSumpt->Fill(jet_pt*jetweight, neutralSum/jet_pt, weight);
          jetptphoSumpt->Fill(jet_pt*jetweight, photonSum/jet_pt, weight);

          jetptchN->Fill(jet_pt*jetweight, chargedN, weight);
          jetptneuN->Fill(jet_pt*jetweight, neutralN, weight);
          jetptphoN->Fill(jet_pt*jetweight, photonN, weight);

          jetptSumSum->Fill(jet_pt*jetweight, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          jetptSumSumpt->Fill(jet_pt*jetweight, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          jetptSumSumrawpt->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);
	
	}
	
      for(int ieta =0; ieta<netabin;ieta++){
          if(jet_eta>deta[ieta]&&jet_eta<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 jetptchMaxEtaBin[dEtaBin]->Fill(jet_pt*jetweight,chargedMax,weight);     //Added
	 jetptchMaxptEtaBin[dEtaBin]->Fill(jet_pt*jetweight,chargedMax/jet_pt,weight);     //Added
//	 jetptneuMaxptEtaBin[dEtaBin]->Fill(jet_pt*jetweight,neutralMax/jet_pt,weight);     //Added
//	 jetptphoMaxptEtaBin[dEtaBin]->Fill(jet_pt*jetweight,photonMax/jet_pt,weight);     //Added
//	 jetptchSumptEtaBin[dEtaBin]->Fill(jet_pt*jetweight,chargedSum/jet_pt,weight);     //Added
//	 jetptneuSumptEtaBin[dEtaBin]->Fill(jet_pt*jetweight,neutralSum/jet_pt,weight);     //Added
//	 jetptphoSumptEtaBin[dEtaBin]->Fill(jet_pt*jetweight,photonSum/jet_pt,weight);     //Added
//	 jetptSumSumrawptEtaBin[dEtaBin]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);
	 jetptSumSumptEtaBin[dEtaBin]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
	 jetptneuMaxrEtaBin[dEtaBin]->Fill(jet_pt*jetweight,neutralMax/TMath::Max(chargedSum,neutralSum),weight);
	}
      } //loop over jet
	
}	//loop over tree

  TString dataType;
  TString out_name;
  dataType = "DATA";
  out_name=Form("%s%s%sskimJetID.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/%s",out_name.Data()),"RECREATE");  

jetptchMax->Write();       //Added
jetptchMaxpt->Write();       //Added
jetptneuMaxpt->Write(); //Added
jetptneuMaxr->Write();   //Added
jetptneuMax->Write();   //Added
jetptphoMaxpt->Write(); //Added 
jetptphoMax->Write();   //Added
jetptchSum->Write();
jetptneuSum->Write();
jetptphoSum->Write();
jetptchSumpt->Write();
jetptneuSumpt->Write();
jetptphoSumpt->Write();
jetptchN->Write();
jetptneuN->Write();
jetptphoN->Write();
jetptSumSum->Write();   //Added
jetptSumSumrawpt->Write();   //Added
jetptSumSumpt->Write();   //Added
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptchMaxEtaBin[ieta]->Write();     //Added
	 jetptchMaxptEtaBin[ieta]->Write();     //Added
	 jetptneuMaxptEtaBin[ieta]->Write();    //Added
	 jetptneuMaxrEtaBin[ieta]->Write();     //Added
	 jetptphoMaxptEtaBin[ieta]->Write();    //Added
	 jetptchSumptEtaBin[ieta]->Write();     //Added
  	 jetptneuSumptEtaBin[ieta]->Write();    //Added
	 jetptphoSumptEtaBin[ieta]->Write();    //Added
	 jetptSumSumptEtaBin[ieta]->Write();    //Added
	 jetptSumSumrawptEtaBin[ieta]->Write();         //Added
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;

  cout<<"working done\n";
 
} 
