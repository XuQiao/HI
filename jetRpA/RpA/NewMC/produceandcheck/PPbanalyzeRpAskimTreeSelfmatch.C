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

class hist_class{
public:
    hist_class();
    void Write();
    TH1F * jetpt;
    TH1F * jetpt_fake;
    TH1F * jetpt_real;
    TH1F * jetpt_selfmatch_fake;
    TH1F * jetpt_selfmatch_real;
    TH1F * jetpt_selfmatch1_fake;
    TH1F * jetpt_selfmatch1_real;
    TH1D * hNmatched;
    TH1D * hNmatched1;
    TH1D * hdeltaR;
    TH2F * hptdeltaptfrc;
    TH2F * hptdeltaptfrcmatch1;
    TH2F * hptdeltaptfrcmatch20;
    TH2F * hptdeltaptfrcmatch21;
    TH2F * hptgenpt;
    TH2F * hptgenptmatch1;
    TH2F * hptgenptmatch20;
    TH2F * hptgenptmatch21;
    TH2F * hdeltaRdeltaptfrcmatch1;
    TH2F * hdeltaRdeltaptfrcmatch20;
    TH2F * hdeltaRdeltaptfrcmatch21;
    TH1F * jetptEtaBin[netabin];
    TH1F * jetptEtaBin_fake[netabin];
    TH1F * jetptEtaBin_real[netabin];
    TH1F * jetptEtaBin_selfmatch_fake[netabin];
    TH1F * jetptEtaBin_selfmatch_real[netabin];
    TH1F * jetptEtaBin_selfmatch1_fake[netabin];
    TH1F * jetptEtaBin_selfmatch1_real[netabin];
};

hist_class::hist_class()
{
    jetpt = new TH1F(Form("jetpt"),Form("jetpt"),1000,0.,1000.);
    jetpt -> Sumw2(); 
    jetpt_fake = new TH1F(Form("jetpt_fake"),Form("jetpt_fake"),1000,0.,1000.);
    jetpt_fake -> Sumw2(); 
    jetpt_real = new TH1F(Form("jetpt_real"),Form("jetpt_real"),1000,0.,1000.);
    jetpt_real -> Sumw2(); 
    jetpt_selfmatch_fake = new TH1F(Form("jetpt_selfmatch_fake"),Form("jetpt_selfmatch_fake"),1000,0.,1000.);
    jetpt_selfmatch_fake -> Sumw2(); 
    jetpt_selfmatch_real = new TH1F(Form("jetpt_selfmatch_real"),Form("jetpt_selfmatch_real"),1000,0.,1000.);
    jetpt_selfmatch_real -> Sumw2(); 
    jetpt_selfmatch1_fake = new TH1F(Form("jetpt_selfmatch1_fake"),Form("jetpt_selfmatch1_fake"),1000,0.,1000.);
    jetpt_selfmatch1_fake -> Sumw2(); 
    jetpt_selfmatch1_real = new TH1F(Form("jetpt_selfmatch1_real"),Form("jetpt_selfmatch1_real"),1000,0.,1000.);
    jetpt_selfmatch1_real -> Sumw2(); 

    hNmatched = new TH1D(Form("hNmatched"),Form("hNmatched"),20,0.,20.);
    hNmatched -> Sumw2();
    hNmatched1 = new TH1D(Form("hNmatched1"),Form("hNmatched1"),20,0.,20.);
    hNmatched1 -> Sumw2();
    hdeltaR = new TH1D(Form("hdeltaR"),Form("hdeltaR"),1000,0.,10.);
    hdeltaR -> Sumw2();
    hptdeltaptfrc = new TH2F(Form("hptdeltaptfrc"),Form("hptdeltaptfrc"),1000,0.,1000.,1000,0.,10.);
    hptdeltaptfrc -> Sumw2();
    hptdeltaptfrcmatch1 = new TH2F(Form("hptdeltaptfrcmatch1"),Form("hdeltaptfrcmatch1"),1000,0.,1000.,1000,0.,10.);
    hptdeltaptfrcmatch1 -> Sumw2();
    hptdeltaptfrcmatch20 = new TH2F(Form("hptdeltaptfrcmatch20"),Form("hptdeltaptfrcmatch20"),1000,0.,1000.,1000,0.,10.);
    hptdeltaptfrcmatch20 -> Sumw2();
    hptdeltaptfrcmatch21 = new TH2F(Form("hptdeltaptfrcmatch21"),Form("hptdeltaptfrcmatch21"),1000,0.,1000.,1000,0.,10.);
    hptdeltaptfrcmatch21 -> Sumw2();
    hptgenpt = new TH2F(Form("hptgenpt"),Form("hptgenpt"),1000,0.,1000.,1000,0.,1000.);
    hptgenpt -> Sumw2();
    hptgenptmatch1 = new TH2F(Form("hptgenptmatch1"),Form("hptgenptmatch1"),1000,0.,1000.,1000,0.,1000.);
    hptgenptmatch1 -> Sumw2();
    hptgenptmatch20 = new TH2F(Form("hptgenptmatch20"),Form("hptgenptmatch20"),1000,0.,1000.,1000,0.,1000.);
    hptgenptmatch20 -> Sumw2();
    hptgenptmatch21 = new TH2F(Form("hptgenptmatch21"),Form("hptgenptmatch21"),1000,0.,1000.,1000,0.,1000.);
    hptgenptmatch21 -> Sumw2();
    hdeltaRdeltaptfrcmatch1 = new TH2F(Form("hdeltaRdeltaptfrcmatch1"),Form("hdeltaRdeltaptfrcmatch1"),1000,0.,10.,1000,0.,10.);
    hdeltaRdeltaptfrcmatch1 -> Sumw2();
    hdeltaRdeltaptfrcmatch20 = new TH2F(Form("hdeltaRdeltaptfrcmatch20"),Form("hdeltaRdeltaptfrcmatch20"),1000,0.,10.,1000,0.,10.);
    hdeltaRdeltaptfrcmatch20 -> Sumw2();
    hdeltaRdeltaptfrcmatch21 = new TH2F(Form("hdeltaRdeltaptfrcmatch21"),Form("hdeltaRdeltaptfrcmatch21"),1000,0.,10.,1000,0.,10.);
    hdeltaRdeltaptfrcmatch21 -> Sumw2();

    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[ieta]->Sumw2();
      jetptEtaBin_fake[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_fake",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin_fake%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_fake[ieta]->Sumw2();
      jetptEtaBin_real[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_real",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin_real%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_real[ieta]->Sumw2();
      jetptEtaBin_selfmatch_fake[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_selfmatch_fake",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f_selfmatch_fake",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_selfmatch_fake[ieta]->Sumw2();
      jetptEtaBin_selfmatch_real[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_selfmatch_real",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f_selfmatch_real",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_selfmatch_real[ieta]->Sumw2();
      jetptEtaBin_selfmatch1_fake[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_selfmatch1_fake",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f_selfmatch_fake",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_selfmatch1_fake[ieta]->Sumw2();
      jetptEtaBin_selfmatch1_real[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_selfmatch1_real",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f_selfmatch_real",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_selfmatch1_real[ieta]->Sumw2();
      }
}

void hist_class::Write()
{
  TString dataType;
  TString out_name;
  dataType = "MC";
  out_name=Form("%s%s%s_useskimSelfmatch.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/cms/store/user/qixu/jetRpA/RpA/NewMC/%s",out_name.Data()),"RECREATE");  
    jetpt->Write();
    jetpt_fake->Write();
    jetpt_real->Write();
    jetpt_selfmatch_fake->Write();
    jetpt_selfmatch_real->Write();
    jetpt_selfmatch1_fake->Write();
    jetpt_selfmatch1_real->Write();
    hNmatched->Write();
    hNmatched1->Write();
    hdeltaR->Write();
    hptdeltaptfrc->Write();
    hptdeltaptfrcmatch1->Write();
    hptdeltaptfrcmatch20->Write();
    hptdeltaptfrcmatch21->Write();
    hptgenpt->Write();
    hptgenptmatch1->Write();
    hptgenptmatch20->Write();
    hptgenptmatch21->Write();
    hdeltaRdeltaptfrcmatch1->Write();
    hdeltaRdeltaptfrcmatch20->Write();
    hdeltaRdeltaptfrcmatch21->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[ieta]->Write(); 
         jetptEtaBin_fake[ieta]->Write();
         jetptEtaBin_real[ieta]->Write();
         jetptEtaBin_selfmatch_fake[ieta]->Write();
         jetptEtaBin_selfmatch1_fake[ieta]->Write();
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void PPbanalyzeRpAskimTree()
{
  std::cout<<"start working\n";
  hist_class *my_hists = new hist_class();

  cout<<"Analyzing MC!"<<endl;

//  TFile *f = new TFile(Form("/store/user/qixu/jetRpA/skimTree/MC%sakPu3PFskimFullInfofile0_10.root",coll.Data()));
  TFile *f = new TFile(Form("/cms/store/user/qixu/jetRpA/skimTree/MC%s%sskimfile0_10final.root",coll.Data(),algo.Data()));
  
  TTree *nt = (TTree*)f->Get("nt");

  Float_t jtpt[100],jteta[100],geneta[100],genphi[100],jtphi[100],rawpt[100],refpt[100],refeta[100],genpt[100];
	Int_t t_chargedN[100], t_neutralN[100], t_photonN[100],subid[100];
Float_t t_chargedSum[100], t_photonSum[100], t_neutralSum[100], t_chargedMax[100], t_neutralMax[100], t_photonMax[100], t_eSum[100], t_muSum[100];
Int_t nref,ngen,hiBin;
  Float_t weight,xSecWeight,vzWeight,centWeight;
  Float_t vz,pthat;
  Int_t HLT_PAJet20_noJetID_v1;
  Int_t HLT_PAJet40_noJetID_v1;
  Int_t HLT_PAJet60_noJetID_v1;
  Int_t HLT_PAJet80_noJetID_v1;
  Int_t HLT_PAJet100_noJetID_v1;

  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("hiBin",&hiBin);
  nt->SetBranchAddress("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_noJetID_v1);
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("ngen",&ngen);
  nt->SetBranchAddress("rawpt",rawpt);
  nt->SetBranchAddress("jtpt",jtpt);
  nt->SetBranchAddress("refpt",refpt);
  nt->SetBranchAddress("refeta",refeta);
  nt->SetBranchAddress("genpt",genpt);
  nt->SetBranchAddress("subid",subid);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("geneta",geneta);
  nt->SetBranchAddress("genphi",genphi);
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
    nt->SetBranchAddress("pthat",&pthat);
    nt->SetBranchAddress("xSecWeight",&xSecWeight);
    nt->SetBranchAddress("vzWeight",&vzWeight);
    nt->SetBranchAddress("centWeight",&centWeight);
    nt->SetBranchAddress("weight",&weight);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
//    for(int i=0; i<50000; i++){
    nt->GetEntry(i);
if(TMath::Abs(vz)>15) continue;
 if(!HLT_PAJet20_noJetID_v1 && !HLT_PAJet40_noJetID_v1 && !HLT_PAJet60_noJetID_v1 && !HLT_PAJet80_noJetID_v1 && !HLT_PAJet100_noJetID_v1 ) continue;
        
    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
for(int j4i = 0; j4i < nref; j4i++){
    double jet_pt = jtpt[j4i];
    double ref_pt = refpt[j4i];
    double ref_eta = refeta[j4i];
    double jet_eta = jteta[j4i];  
    double jet_phi = jtphi[j4i]; 
    double raw_pt = rawpt[j4i];
	if(raw_pt<22.0 || fabs(jet_eta)>5.0) continue;
        if(jet_pt>4*pthat) continue;

        int matchflag=0;
        int matchflag1=0;
        double deltaptfrc[10];
        double genpt_match[10];
        double deltaR_[10];
        if(TMath::Abs(jet_eta)<=2){
	for(int j5i = 0; j5i < ngen ; j5i++) {
        double gen_pt=genpt[j5i];
        double gen_eta=geneta[j5i];
        double gen_phi=genphi[j5i];
        TVector3 jetVec, genVec;
        jetVec.SetPtEtaPhi(jet_pt,jet_eta,jet_phi);
        genVec.SetPtEtaPhi(gen_pt,gen_eta,gen_phi);
        double deltaR = jetVec.DeltaR(genVec);
        my_hists->hdeltaR->Fill(deltaR,weight);
        my_hists->hptdeltaptfrc->Fill(jet_pt,TMath::Abs(jet_pt-gen_pt)/jet_pt,weight);
        my_hists->hptgenpt->Fill(jet_pt,gen_pt,weight);
        if(deltaR<0.3){ 
        deltaR_[matchflag]=deltaR;
        genpt_match[matchflag]=gen_pt;
        deltaptfrc[matchflag]=TMath::Abs(jet_pt-gen_pt)/jet_pt;
        matchflag++;
        }
        if(deltaR<0.3 && TMath::Abs(jet_pt-gen_pt)/jet_pt < 0.5){
           matchflag1++; 
        }
        }
        my_hists->hNmatched->Fill(matchflag,weight);
        my_hists->hNmatched1->Fill(matchflag1,weight);
        if(matchflag==1){
        my_hists->hptdeltaptfrcmatch1->Fill(jet_pt,deltaptfrc[0],weight);
        my_hists->hdeltaRdeltaptfrcmatch1->Fill(deltaR_[0],deltaptfrc[0],weight);
        my_hists->hptgenptmatch1->Fill(jet_pt,genpt_match[0],weight);
        }
        if(matchflag==2){
          int  i0 = deltaR_[0]<deltaR_[1]?0:1;
          int  i1 = deltaR_[0]<deltaR_[1]?1:0;
        my_hists->hptdeltaptfrcmatch20->Fill(jet_pt,deltaptfrc[i0],weight);
        my_hists->hptdeltaptfrcmatch21->Fill(jet_pt,deltaptfrc[i1],weight);
        my_hists->hptgenptmatch20->Fill(jet_pt,genpt_match[i0],weight);
        my_hists->hptgenptmatch21->Fill(jet_pt,genpt_match[i1],weight);
        my_hists->hdeltaRdeltaptfrcmatch20->Fill(deltaR_[i0],deltaptfrc[i0],weight);
        my_hists->hdeltaRdeltaptfrcmatch21->Fill(deltaR_[i1],deltaptfrc[i1],weight);
        }
        }

        if(coll =="PPb") {
            ref_eta = ref_eta+0.465;
            jet_eta = jet_eta+0.465;
        }
        if(coll == "PbP"){
            ref_eta = ref_eta-0.465;
            jet_eta = jet_eta-0.465;
        }

      int dEtaBin = -1.;
     if(TMath::Abs(jet_eta)<=1.) {
	if(ref_pt<0)	//Added
          my_hists->jetpt_fake->Fill(jet_pt,weight);
	else
          my_hists->jetpt_real->Fill(jet_pt,weight);
        if(matchflag==0)    //self-matching fake
          my_hists->jetpt_selfmatch_fake->Fill(jet_pt,weight);
        else
          my_hists->jetpt_selfmatch_real->Fill(jet_pt,weight);
        if(matchflag1==0)    //self-matching fake
          my_hists->jetpt_selfmatch1_fake->Fill(jet_pt,weight);
        else
          my_hists->jetpt_selfmatch1_real->Fill(jet_pt,weight);
	my_hists->jetpt->Fill(jet_pt,weight);
	}	

      for(int ieta =0; ieta<netabin;ieta++){
          if(jet_eta>deta[ieta]&&jet_eta<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 my_hists->jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
         if(ref_pt<0)
            my_hists->jetptEtaBin_fake[dEtaBin]->Fill(jet_pt,weight);
         else
            my_hists->jetptEtaBin_real[dEtaBin]->Fill(jet_pt,weight);
         if(matchflag==0)
            my_hists->jetptEtaBin_selfmatch_fake[dEtaBin]->Fill(jet_pt,weight);
         else
            my_hists->jetptEtaBin_selfmatch_real[dEtaBin]->Fill(jet_pt,weight);
         if(matchflag1==0)
            my_hists->jetptEtaBin_selfmatch1_fake[dEtaBin]->Fill(jet_pt,weight);
         else
            my_hists->jetptEtaBin_selfmatch1_real[dEtaBin]->Fill(jet_pt,weight);
	}
      } //loop over jet
}	//loop over tree

  my_hists->Write();
  cout<<"working done\n";
 
}  
