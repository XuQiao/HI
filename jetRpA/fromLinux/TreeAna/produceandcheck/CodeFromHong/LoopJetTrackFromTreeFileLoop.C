#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TCut.h"
#include "/home/xuq/Documents/HI/RpA/output/TreeAna/Corrections/TrackCorrector2D.C"
#include <stdio.h>

using namespace std;

// ******* GLOBAL DECLARATIONS **********
const int QCDpthatBins = 8;
//const int QCDpthatBins = 9;
//const int QCDpthatBins = 1;
const int dataFiles = 10292;
int startfile ;
int endfile ;

const int sigma = 70 ;  //units in mb
const double Ncoll = 6.9 ; //for pPb Ncoll
const double dphicut = 2*TMath::Pi()/3. ;
const double leadingjetcut = 120. ;
const double subleadjetcut = 30. ;
double conesize = 0.3; 
//***************************************

//const double deta[]={-2.2, -1.2, -0.7, -0.3, 0.3, 0.7,1.2,2.2} ;
const double deta[]={-2., -1.5, -1.0, -0.5, 0.5, 1.0,1.5,2.0} ;
const int netabin = sizeof(deta)/sizeof(Double_t)-1 ;
//const Double_t jetPtBin[]={3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
const Double_t jetPtBin[]={3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};

const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;
const int pthatbin[10] = {15,30,50,80,120,170,220,280,370, 9999};
const double wght[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 1.088E-07, 0};
const Double_t TrkBin[]={0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 138., 155.2, 172.4, 189.6, 206.8};
const int nTrkBin = sizeof(TrkBin)/sizeof(Double_t)-1 ;
const int ntrketabin = 7 ;
const double dtrketa[]={-1.8, -1.3, -0.8, -0.3, 0.3, 0.8,1.3,1.8} ;

const float pPbLumi=20.7; //unit in nb, need to convert to mb for cross sections
const float PbpLumi=14.0; //unit in nb, need to convert to mb for cross sections
const bool SavePlot=kFALSE;
const bool IsMC = kFALSE ;
const TString coll = "PPb";
const TString algo = "akPu3PF" ;
TString TrigName = "Jet20" ; //"JetMB","Jet20","Jet40","Jet60","Jet80","Jet100"
const int nbin = netabin+1 ;
const int ntrkbin = ntrketabin+1  ;

//do tracking corrections
const bool doTrackCorrections = kTRUE; //for tracking efficiency correction

class Histos{
public:
    Histos();
    void Delete();
    void Write();
    
    TString desc;
    
    TH1D * NEvents;
    TH1D * JetAbove20;
    TH1D * JetAbove40;
    TH1D * JetAbove60;
    TH1D * JetAbove75;
    TH1D * JetAbove95;
    TH1D * JetAbove120;
    
    TH1F * jetpt[nbin];
    TH1F * trkpt[ntrkbin];
    TH1F * jettrkpt[ntrkbin];

    TH1F * jet20 ;
    TH1F * jet40 ;
    TH1F * jet60 ;
    TH1F * jet80 ;
    TH1F * jet100 ;

    TH1F * leadjetpt;
    TH2F * leadjetptDiEta;
    TH2F * leadjetptDPhi;
    TH2F * leadjetptMjj;
    TH2F * leadjetptratio;
    TH2F * MjjDiEta;
    TH2F * AvgPtBalanceMid;
    TH2F * AvgPtBalanceFor;
    TProfile * ProbeEtaMeanB;
    TH1F * subjetpt;
    TH1F * dijetMass;
   
    
};
Histos::Histos()
{
    NEvents = new TH1D("Nevents_0-100%","Nevents_0-100%",100,0,2.);
    NEvents -> Sumw2();
    jet20 = new TH1F(Form("jet20EtaBin%.f_%.f",-1.0*10, 1.0*10),Form("jet20EtaBin%.f_%.f",-1.0*10, 1.0*10),1000,0.,1000.);
    jet20->Sumw2();
    jet40 = new TH1F(Form("jet40EtaBin%.f_%.f",-1.0*10, 1.0*10),Form("jet40EtaBin%.f_%.f",-1.0*10, 1.0*10),1000,0.,1000.);
    jet40->Sumw2();
    jet60 = new TH1F(Form("jet60EtaBin%.f_%.f",-1.0*10, 1.0*10),Form("jet60EtaBin%.f_%.f",-1.0*10, 1.0*10),1000,0.,1000.);
    jet60->Sumw2();
    jet80 = new TH1F(Form("jet80EtaBin%.f_%.f",-1.0*10, 1.0*10),Form("jet80EtaBin%.f_%.f",-1.0*10, 1.0*10),1000,0.,1000.);
    jet80->Sumw2();
    jet100 = new TH1F(Form("jet100EtaBin%.f_%.f",-1.0*10, 1.0*10),Form("jet100EtaBin%.f_%.f",-1.0*10, 1.0*10),1000,0.,1000.);
    jet100->Sumw2();
    
    JetAbove20 = new TH1D("NeventsWithJetBelow40_0-100%","NeventsWithJetBelow40_0-100%",1,0,1);
    JetAbove20 -> Sumw2();
    
    JetAbove40 = new TH1D("NeventsWithJetAbove40_0-100%","NeventsWithJetAbove40_0-100%",1,0,1);
    JetAbove60 = new TH1D("NeventsWithJetAbove60_0-100%","NeventsWithJetAbove60_0-100%",1,0,1);
    JetAbove75 = new TH1D("NevntsWithJetAbove75_0-100%","NeventsWithJetAbove75_0-100%",1,0,1);
    JetAbove95 = new TH1D("NevntsWithJetAbove95_0-100%","NeventsWithJetAbove95_0-100%",1,0,1);
    JetAbove120 = new TH1D("NevntsWithJetAbove120_0-100%","NeventsWithJetAbove120_0-100%",1,0,1);
    
    for(int ieta=0; ieta <=netabin; ieta++){
        if(ieta==netabin)
          //  jetpt[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",-1.0*10, 1.0*10),Form("jetptEtaBin%.f_%.f",-1.0*10, 1.0*10),1000, 0., 1000); //nJetPtBin, jetPtBin);
            jetpt[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",-1.0*10, 1.0*10),Form("jetptEtaBin%.f_%.f",-1.0*10, 1.0*10),nJetPtBin, jetPtBin);
            else
                jetpt[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), nJetPtBin, jetPtBin);
              //  jetpt[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000, 0., 1000.) ;  //nJetPtBin, jetPtBin);
        
        jetpt[ieta]->Sumw2();
    }

      for(int ieta=0; ieta <=ntrketabin; ieta++){
        if(ieta==ntrketabin){
            trkpt[ieta] = new TH1F(Form("trackptEtaBin%.f_%.f",-1.0*10, 1.0*10),Form("trackptEtaBin%.f_%.f",-1.0*10, 1.0*10),nTrkBin, TrkBin);
            jettrkpt[ieta] = new TH1F(Form("jettrackptEtaBin%.f_%.f",-1.0*10, 1.0*10),Form("jettrackptEtaBin%.f_%.f",-1.0*10, 1.0*10),nTrkBin, TrkBin);
           }
          else{
                trkpt[ieta] = new TH1F(Form("trackptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("trackptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), nTrkBin, TrkBin);
                jettrkpt[ieta] = new TH1F(Form("jettrackptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jettrackptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), nTrkBin, TrkBin);
        }
       trkpt[ieta]->Sumw2();
       jettrkpt[ieta]->Sumw2();
    }
          leadjetpt = new TH1F(Form("leadingjetpt"), Form("leadingjetpt"), 1000, 0., 1000.);
       leadjetpt->Sumw2();
        subjetpt = new TH1F(Form("subleadingjetpt"), Form("subleadingjetpt"), 1000, 0., 1000.);
       subjetpt->Sumw2();
       dijetMass= new TH1F(Form("dijetMass"), Form("dijetMass"), 1000, 0., 1000.);
       dijetMass->Sumw2();
        leadjetptDiEta = new TH2F(Form("leadjetptDiEta"), Form("leadjetptDiEta"), 1000, 0., 1000., 100, -5., 5.);
       leadjetptDiEta->Sumw2();
       leadjetptDiEta->GetXaxis()->SetTitle("p_{T}^{leading} (GeV/c)");
        leadjetptDiEta->GetYaxis()->SetTitle("(#eta_{1}+#eta_{2})/2");
        leadjetptDPhi = new TH2F(Form("leadjetptDeltaPhi"), Form("leadjetptDeltaPhi"), 1000, 0., 1000., 80, 0., 8.);
       leadjetptDPhi->Sumw2();
       leadjetptDPhi->GetXaxis()->SetTitle("p_{T}^{leading} (GeV/c)");
        leadjetptDPhi->GetYaxis()->SetTitle("|(#phi_{1}-#phi_{2})|");
        leadjetptratio = new TH2F(Form("leadjetptratio"), Form("leadjetptratio"), 1000, 0., 1000., 200, 0., 2.);
       leadjetptratio->Sumw2();
       leadjetptratio->GetXaxis()->SetTitle("p_{T}^{leading} (GeV/c)");
        leadjetptratio->GetYaxis()->SetTitle("p_{T}^{2}/p_{T}^{1}");
        leadjetptMjj = new TH2F(Form("leadjetptMjj"), Form("leadjetptMjj"), 1000, 0., 1000., 1000, 0., 1000.);
       leadjetptMjj->Sumw2();
       leadjetptMjj->GetXaxis()->SetTitle("p_{T}^{leading} (GeV/c)");
        leadjetptMjj->GetYaxis()->SetTitle("M_{jj} (GeV/c)^{2}");
        MjjDiEta = new TH2F(Form("MjjDiEta"), Form("MjjDiEta"), 1000, 0., 1000., 100, -5., 5.);
        MjjDiEta->Sumw2();
        MjjDiEta->GetYaxis()->SetTitle("(#eta_{1}+#eta_{2})/2");
        MjjDiEta->GetXaxis()->SetTitle("M_{jj} (GeV/c)^{2}");
  
       AvgPtBalanceMid = new TH2F(Form("AvgPtBalanceMid"), Form("AvgPtBalanceMid"), 1000, 0., 1000., 100, -2., 2.);
       AvgPtBalanceMid->Sumw2();
       AvgPtBalanceMid->GetYaxis()->SetTitle("2*(p_{T}^{probe}-p_{T}^{ref})/(p_{T}^{ref}+p_{T}^{probe})");
       AvgPtBalanceMid->GetXaxis()->SetTitle("(p_{T}^{ref}+p_{T}^{probe})/2");

       AvgPtBalanceFor = new TH2F(Form("AvgPtBalanceForward"), Form("AvgPtBalanceForward"), 1000, 0., 1000., 100, -2., 2.);
       AvgPtBalanceFor->Sumw2();
       AvgPtBalanceFor->GetYaxis()->SetTitle("2*(p_{T}^{probe}-p_{T}^{ref})/(p_{T}^{ref}+p_{T}^{probe})");
       AvgPtBalanceFor->GetXaxis()->SetTitle("(p_{T}^{ref}+p_{T}^{probe})/2");

       ProbeEtaMeanB = new TProfile(Form("ProbeEtaMeanB"), Form("ProbeEtaMeanB"), 600, -3., 3.);
       ProbeEtaMeanB->Sumw2();
       ProbeEtaMeanB->GetXaxis()->SetTitle("#eta^{probe}");
       ProbeEtaMeanB->GetYaxis()->SetTitle("Relative Response");
}
void Histos::Delete()
{
    for(int ieta=0; ieta<=netabin; ieta++) delete jetpt[ieta];
    for(int ieta=0; ieta<=ntrketabin; ieta++) delete trkpt[ieta];
    for(int ieta=0; ieta<=ntrketabin; ieta++) delete jettrkpt[ieta];

    delete    jet20;
    delete jet40;
    delete jet60;
    delete jet80;
    delete jet100;

    delete leadjetpt;
    delete leadjetptDiEta;
    delete leadjetptDPhi;
    delete leadjetptMjj;
   delete leadjetptratio;
    delete MjjDiEta;
    delete subjetpt;
    delete dijetMass;
    delete ProbeEtaMeanB;
    delete AvgPtBalanceMid;
    delete AvgPtBalanceFor;
}
void Histos::Write()
{
    TString dataType;
    TString out_name;
    if(IsMC) dataType = "MC";
    else dataType = "DATA";
    
 //   out_name=Form("%s%s%s%sJetSpectra.root",dataType.Data(),coll.Data(),TrigName.Data(),algo.Data());
 //   out_name=Form("%s%s%sJetSpectraCombineTriggerQiaoTreeJetPtEtaWeightNoJetIDCut.root",dataType.Data(),coll.Data(),algo.Data());
 //   out_name=Form("%s%s%sJetSpectraCombineTriggerQiaoTreeNoWeightNoJetIDCut.root",dataType.Data(),coll.Data(),algo.Data());
 //   out_name=Form("%s%s%sJetSpectraCombineTriggerYueShiJECNoWeightNoJetIDCut.root",dataType.Data(),coll.Data(),algo.Data());
    out_name=Form("%s%s%sJetSpectraCombineTriggerJetPtNoJetIDCutNoResidualTrkCorrFile%d_%d.root",dataType.Data(),coll.Data(),algo.Data(), startfile, endfile);
 //   out_name=Form("%s%s%sJetSpectraCombineTriggerJetPtEtaWeightNewBin.root",dataType.Data(),coll.Data(),algo.Data());
//    out_name=Form("%s%s%sJetSpectraCombineTriggerNoWeightNoJetIDCutRawPt.root",dataType.Data(),coll.Data(),algo.Data());
    
    TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/%s",out_name.Data()),"RECREATE");
    
//    NEvents->Write();
//    JetAbove20->Write();
//    JetAbove40->Write();
//    JetAbove60->Write();
//    JetAbove75->Write();
//    JetAbove95->Write();
//    JetAbove120->Write();
    for(int ieta=0; ieta<=netabin; ieta++) jetpt[ieta]->Write();
    for(int ieta=0; ieta<=ntrketabin; ieta++) trkpt[ieta]->Write(); 
    for(int ieta=0; ieta<=ntrketabin; ieta++) jettrkpt[ieta]->Write(); 

    jet20->Write();
    jet40->Write();
    jet60->Write();
    jet80->Write();
    jet100->Write();
    
    leadjetpt->Write();
    leadjetptDiEta->Write();
    leadjetptDPhi->Write();
    leadjetptMjj->Write();
    leadjetptratio->Write();
    MjjDiEta->Write();
    subjetpt->Write();
    dijetMass->Write();
    AvgPtBalanceFor->Write();
    AvgPtBalanceMid->Write();
    ProbeEtaMeanB->Write();
    out_file->Close();
    cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void LoopJetTrackFromTreeFileLoop()
{
/*  
  gStyle->SetErrorX(0.);
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetOptLogz(1);
  gStyle->SetPadRightMargin(0.13);	
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
*/

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

     if(IsMC){
     startfile=0;
     endfile=QCDpthatBins ;
    }
   else {
//      startfile=atoi(getenv("FIRST"));
//    endfile=atoi(getenv("LAST"));
    }
  
  // input files
    std::cout<<"start working\n";
//    cout<<"analyze  file = "<<startfile << " to " <<endfile <<endl;
    cout <<"MC = " << IsMC <<endl;
    
    // Initialize Histograms
    Histos *hist = new Histos();
    
    if (IsMC) cout<<"This is a MC analysis"<<endl;
    else cout<< "This is a data analysis"<<endl;
    
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
     
        TrackCorrector2D corr("/home/xuq/Documents/HI/RpA/output/TreeAna/Corrections/trackCorrections_HIN12017v5_XSecWeighted.root");
   if(doTrackCorrections) corr.load("trkCorr_HIN12017");

      TString fileName ;
    
    TFile *tFile;
    TTree *tJet;
    
    if(IsMC==kTRUE){
    /*    if(coll=="PbP")
          fileName = Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/%sMCKurtForestVzWeight_ppReco_%s_QCDjetQiaoIDCut_JetPt0noIPupperCut.root", coll.Data(), algo.Data());
       else
   */     fileName = Form("/home/xuq/Documents/HI/RpA/output/skimTree/MC%s%sskimfile0_8.root", coll.Data(), algo.Data());
       }
    else {
        if(coll=="PbP")
          //  fileName = Form("/afs/cern.ch/work/y/ymao/CMSSW_6_2_0_pre7/src/histos/pPbdata_ppreco_%s_jetTrig_etashift_FinalTC_FullPUForest_AddSvtxInfo_NoEvtSelection_WithChgMax_EtaResCorr.root", algo.Data());
         //   fileName = Form("/afs/cern.ch/work/y/ymao/public/RpA/ForQiao/PbPdataKurt_ppReco_%s_AlljetTrigKurtTrCombFile0_6815_JetPt0noIPupperCut.root", algo.Data());
         //   fileName = Form("/afs/cern.ch/work/y/ymao/CMSSW_6_2_0_pre7/src/histos/PbPdata_ppReco_%s_AlljetTrigKurtTrCombFile0_6815_JetPt0noIPupperCut.root", algo.Data());
            fileName = Form("DATAPbP%sskimFullInfofile0_6815.root", algo.Data());
        else
          //  fileName = Form("/Users/ymao/group/CMS/anaOutputs/pPb/RpPb/data_ppReco_%s_AlljetTrigKurtCombinationWeight_JetPt0noIPupperCut.root", algo.Data());
     //       fileName = Form("/afs/cern.ch/work/y/ymao/CMSSW_6_2_0_pre7/src/histos/%sdata_ppReco_%s_AlljetTrigKurtCombinationWeight_JetPt0noIPupperCut.root", coll.Data(), algo.Data());
     //       fileName = Form("/afs/cern.ch/work/q/qixu/public/RpA//%sdata_ppReco_%s_AlljetTrigKurtTrCombFile_merged.root", coll.Data(), algo.Data());
         //   fileName = Form("/afs/cern.ch/work/y/ymao/CMSSW_6_2_0_pre7/src/histos/DATA%s%sskimFullInfoCombined.root", coll.Data(), algo.Data());
            fileName = Form("/home/xuq/Documents/HI/RpA/output/skimTree/DATA%s%sskimfile.root", coll.Data(), algo.Data());
    }
    
    tFile = TFile::Open(fileName.Data(),"readonly");
    tJet = (TTree*)tFile->Get("nt");

   
/*    double rawpt;
    double jtpt;
    double refpt;
    double refeta;
    double jteta;
*/    
    float weight;
    float pthat;
    int nref;
/*
     float rawpt;
    float jtpt;
    float refpt;
    float refeta;
    float jteta;
*/
  Float_t rawpt[1000];
  Float_t jtpt[1000];
  Float_t jteta[1000];
  Float_t jty[1000];
  Float_t jtphi[1000];
  Float_t refpt[1000];
  Float_t refeta[1000];

  Int_t chargedN[1000];
  Int_t photonN[1000];
  Int_t neutralN[1000];
  Float_t chargedMax[1000];
  Float_t photonMax[1000];
  Float_t neutralMax[1000];
  Float_t chargedSum[1000];
  Float_t photonSum[1000];
  Float_t neutralSum[1000];
  Float_t muSum[1000];
  Float_t eSum[1000];


    int nTrk;
/*    double trkPt;
    double trkEta;
    double trkPhi;
*/  //  float trkPt;
  //  float trkEta;
  //  float trkPhi;
    Float_t trkPt[1000];
    Float_t trkPhi[1000];
    Float_t trkEta[1000];

    int bin;
    int run;
    int evtSelection;
    Int_t evt;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
  Float_t vz;
  Int_t pPAcollisionEventSelectionPA;
  Int_t pprimaryvertexFilter;
  Int_t pHBHENoiseFilter;
    int HLT_PAJet20_noJetID_v1;
    int HLT_PAJet40_noJetID_v1;
    int HLT_PAJet60_noJetID_v1;
    int HLT_PAJet80_noJetID_v1;
    int HLT_PAJet100_noJetID_v1;
 //   int HLT_PAZeroBiasPixel_SingleTrack_v1;
    int HLT_PAJet20_noJetID_v1_Prescl;
    int HLT_PAJet40_noJetID_v1_Prescl;
    int HLT_PAJet60_noJetID_v1_Prescl;
    int HLT_PAJet80_noJetID_v1_Prescl;
    int HLT_PAJet100_noJetID_v1_Prescl;
//    int HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl;

    tJet->SetBranchAddress("nref",&nref);
    tJet->SetBranchAddress("rawpt" , rawpt );
    tJet->SetBranchAddress("jtpt" , jtpt );
   if(IsMC){
    tJet->SetBranchAddress("refpt", refpt);
    tJet->SetBranchAddress("refeta", refeta);
   }
    tJet->SetBranchAddress("jteta", jteta);
    
    tJet->SetBranchAddress("chargedN",chargedN);
    tJet->SetBranchAddress("photonN",photonN);
    tJet->SetBranchAddress("neutralN",neutralN);

    tJet->SetBranchAddress("chargedMax",chargedMax);
    tJet->SetBranchAddress("photonMax",photonMax);
    tJet->SetBranchAddress("neutralMax",neutralMax);
    tJet->SetBranchAddress("chargedSum",chargedSum);
    tJet->SetBranchAddress("photonSum",photonSum);
    tJet->SetBranchAddress("neutralSum",neutralSum);
    tJet->SetBranchAddress("muSum",muSum);
    tJet->SetBranchAddress("eSum",eSum);
 
    tJet->SetBranchAddress("nTrk",&nTrk);
    tJet->SetBranchAddress("trkPt",trkPt);
    tJet->SetBranchAddress("trkPhi",trkPhi);
    tJet->SetBranchAddress("trkEta",trkEta);
    
    tJet->SetBranchAddress("evt",&evt);
    tJet->SetBranchAddress("pthat",&pthat);
    tJet->SetBranchAddress("weight",&weight);
    tJet->SetBranchAddress("vz",&vz);
    tJet->SetBranchAddress("run",&run);
    tJet->SetBranchAddress("evtSelection",&evtSelection);
    tJet->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
    tJet->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    tJet->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
    tJet->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4);
    tJet->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4);

    tJet->SetBranchAddress("hibin",&bin);
    tJet->SetBranchAddress("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_noJetID_v1);
    tJet->SetBranchAddress("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_noJetID_v1);
    tJet->SetBranchAddress("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_noJetID_v1);
    tJet->SetBranchAddress("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_noJetID_v1);
    tJet->SetBranchAddress("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_noJetID_v1);
//    tJet->SetBranchAddress("HLT_ZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
    tJet->SetBranchAddress("HLT_PAJet20_noJetID_v1_Prescl",&HLT_PAJet20_noJetID_v1_Prescl);
    tJet->SetBranchAddress("HLT_PAJet40_noJetID_v1_Prescl",&HLT_PAJet40_noJetID_v1_Prescl);
    tJet->SetBranchAddress("HLT_PAJet60_noJetID_v1_Prescl",&HLT_PAJet60_noJetID_v1_Prescl);
    tJet->SetBranchAddress("HLT_PAJet80_noJetID_v1_Prescl",&HLT_PAJet80_noJetID_v1_Prescl);
    tJet->SetBranchAddress("HLT_PAJet100_noJetID_v1_Prescl",&HLT_PAJet100_noJetID_v1_Prescl);
//    tJet->SetBranchAddress("HLT_ZeroBiasPixel_SingleTrack_v1_Prescl",&HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl);

    
    Int_t Nevents = 0;
    Int_t Nevt_40_60 = 0;
    Int_t Nevt_60_75 = 0;
    Int_t Nevt_75_95 = 0;
    Int_t Nevt_95_120 = 0;
    Int_t Nevt_120 = 0;

       TVector3 jet1_vec;
       TVector3 jet2_vec;
    TVector3 jet_vec;
    TVector3 track_vec;

//    TString CorrfileName;
//    CorrfileName="~/Dropbox/CMS/work/RpA/HIN_14-001/TreeAnalysis/Corrections/pPbHijing_EffCorr_VRw.root";
//    Corrector corr_Qiao(CorrfileName);
    

    Int_t nentries = tJet->GetEntries();
    
      //  weight =1.;
    for(int i=0; i<nentries; i++){
        tJet->GetEntry(i);
             //Cut to remove events that correspond to the twiki "good events" but not the golden lumi filter
      if(!IsMC){
      //  if(((int)run>211256 ) || ((int)run<210676 )) continue;
       if(coll=="PPb"){
        if(((int)run>211256 ) ) continue; //! pPb run
       }
      else {
        if(((int)run<=211256 ) ) continue; //! Pbp run
       }
        if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA ) continue;
      }
      if(coll=="PPb" || coll=="PbP" ){  if(fabs(vz)>15.) continue;}

       double jetweight = 1;     
        
        double HFbin = hiHFplusEta4;
        
       //lets get the leading jet information here
      int leadJetIndex = -1;
      float leadJetPt = 0. ;
     for(int ij=0;ij<nref;ij++){
        if(rawpt[ij]<22.|| fabs(jteta[ij])>3.0) continue ;
                float jetpt = jtpt[ij];
                if(jetpt>leadJetPt){
                 leadJetPt = jetpt;
                leadJetIndex = ij;
              }

   } //! end of leading jet search

  
        bool jetAbove = false;
        bool jetAbove40 = false;
        bool jetAbove60 = false;
        bool jetAbove75 = false;
        bool jetAbove95 = false;
        bool jetAbove120 = false;
        
        if(HLT_PAJet20_noJetID_v1 && leadJetPt>40.) jetAbove40 = true;
        if(HLT_PAJet20_noJetID_v1 && leadJetPt>60.) jetAbove60 = true;
        if(HLT_PAJet20_noJetID_v1 && leadJetPt>75.) jetAbove75 = true;
        if(HLT_PAJet20_noJetID_v1 && leadJetPt>95.) jetAbove95 = true;
        if(HLT_PAJet20_noJetID_v1 && leadJetPt>120.) jetAbove120 = true;
        
        if(!jetAbove40) Nevents++;
        else if(jetAbove40 && !jetAbove60) Nevt_40_60++;
        else if(jetAbove60 && !jetAbove75) Nevt_60_75++;
        else if(jetAbove75 && !jetAbove95) Nevt_75_95++;
        else if(jetAbove95 && !jetAbove120) Nevt_95_120++;
        else Nevt_120++;
        
        if(TrigName=="Jet20") jetAbove = jetAbove40 && !jetAbove60;
        else if(TrigName=="Jet40") jetAbove = jetAbove60 && !jetAbove75;
        else if(TrigName=="Jet60") jetAbove = jetAbove75 && !jetAbove95;
        else if(TrigName=="Jet80") jetAbove = jetAbove95 && !jetAbove120;
        else if(TrigName=="Jet100") jetAbove = jetAbove120;
        else jetAbove = !jetAbove40;
        
  //      if(!jetAbove) continue;//jet trigger cut
        
        if(i%100000==1)    cout<<"analyzing "<< i <<" th event"<<endl;
        
//        if(IsMC==kFALSE)
//            weight*=corr_Qiao.getEventWeightHFPlus4bak(HFbin,kTRUE);//efficiency weight
//        else
//          weight*=weight; // weight for MC
        
        hist->NEvents->Fill(1,weight);

        for(int j4i = 0; j4i < nref; j4i++){
            double jet_pt = jtpt[j4i];
            double jet_eta = jteta[j4i];
            jetweight =1. ;
           if(!IsMC){
             if(coll=="PPb" || coll=="PbP") jetweight*=(fUE->Eval(jet_pt))*C_rel->GetBinContent(C_rel->FindBin(jet_eta));
         //    if(coll=="PPb") jetweight*=C_rel->GetBinContent(C_rel->FindBin(jet_eta));
            }
            int dEtaBin = -1.;
              if (rawpt[j4i]<22) continue;
         //   if((chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt>1.01) continue ;
            if(coll=="PPb"){
                if(TMath::Abs(jet_eta+0.465)<=1.) {
           //     if(TMath::Abs(jet_eta)<=1.) {
            //       hist->jetpt[netabin]->Fill(jet_pt,weight);
            //        hist->jetpt[netabin]->Fill(rawpt,weight);
         //           hist->jetpt[netabin]->Fill(rawpt*jetweight,weight);
                    hist->jetpt[netabin]->Fill(jet_pt*jetweight,weight);
                    if(HLT_PAJet20_noJetID_v1)hist->jet20->Fill(jet_pt,HLT_PAJet20_noJetID_v1_Prescl);
                    if(HLT_PAJet40_noJetID_v1)hist->jet40->Fill(jet_pt,HLT_PAJet40_noJetID_v1_Prescl);
                    if(HLT_PAJet60_noJetID_v1)hist->jet60->Fill(jet_pt,HLT_PAJet60_noJetID_v1_Prescl);
                    if(HLT_PAJet80_noJetID_v1)hist->jet80->Fill(jet_pt,HLT_PAJet80_noJetID_v1_Prescl);
                    if(HLT_PAJet100_noJetID_v1)hist->jet100->Fill(jet_pt,HLT_PAJet100_noJetID_v1_Prescl);
                }
            }
          else {
                 if(TMath::Abs(jet_eta-0.465)<=1.) {
              //   if(TMath::Abs(jet_eta)<=1.) {
               //     hist->jetpt[netabin]->Fill(jet_pt,weight);
                //    hist->jetpt[netabin]->Fill(rawpt,weight);
                    hist->jetpt[netabin]->Fill(jet_pt*jetweight,weight);
               //     hist->jetpt[netabin]->Fill(rawpt*jetweight,weight);
                    if(HLT_PAJet20_noJetID_v1)hist->jet20->Fill(jet_pt,HLT_PAJet20_noJetID_v1_Prescl);
                    if(HLT_PAJet40_noJetID_v1)hist->jet40->Fill(jet_pt,HLT_PAJet40_noJetID_v1_Prescl);
                    if(HLT_PAJet60_noJetID_v1)hist->jet60->Fill(jet_pt,HLT_PAJet60_noJetID_v1_Prescl);
                    if(HLT_PAJet80_noJetID_v1)hist->jet80->Fill(jet_pt,HLT_PAJet80_noJetID_v1_Prescl);
                    if(HLT_PAJet100_noJetID_v1)hist->jet100->Fill(jet_pt,HLT_PAJet100_noJetID_v1_Prescl);
                }
            }   
            for(int ieta =0; ieta<netabin;ieta++){
                if(coll=="PPb"){
                    if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
                  //  if((jet_eta)>deta[ieta]&&(jet_eta)<=deta[ieta+1]) dEtaBin=ieta;
                }
               else {
                 if((jet_eta-0.465)>deta[ieta]&&(jet_eta-0.465)<=deta[ieta+1]) dEtaBin=ieta;
               //     if((jet_eta)>deta[ieta]&&(jet_eta)<=deta[ieta+1]) dEtaBin=ieta;
                }
            }//assign the eta bin for jets
   //         if(dEtaBin!=-1) hist->jetpt[dEtaBin]->Fill(rawpt*jetweight,weight);
            if(dEtaBin!=-1) hist->jetpt[dEtaBin]->Fill(jet_pt*jetweight,weight);
        //    if(dEtaBin!=-1) hist->jetpt[dEtaBin]->Fill(jet_pt,weight);
         //   if(dEtaBin!=-1) hist->jetpt[dEtaBin]->Fill(rawpt,weight);
        }//jet loop
 /*   
             //Leading Jets seach
      double leadingJetPt = -1. ;
      double leadingJetPhi = -1. ;
      double leadingJetEta = -1. ;
      Int_t leadingJetIndex = -1 ;
      double subleadingJetPt = -1. ;
      double subleadingJetPhi = -1. ;
      double subleadingJetEta = -1. ;
      Int_t subleadingJetIndex = -1 ;
      double leadres = 1. ;
      double subleadres = 1. ;
      jet1_vec.SetPtEtaPhi(0, 0, 0);
      jet2_vec.SetPtEtaPhi(0, 0, 0);
     for(int j = 0; j < nref ; j++) {
         if (fabs(jteta[j])>3.0) continue;
         if (rawpt[j]<22) continue;
         if (jtpt[j]>leadingJetPt) {
            leadingJetPt = jtpt[j];
            leadingJetIndex = j;
            leadingJetPhi = jtphi[j];
            leadingJetEta = jteta[j];
         }
      }

        //  if(coll=="PPb" || coll=="PbP") residual*=(fUE->Eval(leadingJetPt))*C_rel->GetBinContent(C_rel->FindBin(leadingJetEta)); 
       if(!IsMC){
          if(coll=="PPb" || coll=="PbP") leadres*=C_rel->GetBinContent(C_rel->FindBin(leadingJetEta)); 
       //    leadingJetPt*=leadres;
     }
        //if leading jets above cuts found
    if(leadingJetIndex!=-1){ //search for subleading jet
           jet1_vec.SetPtEtaPhi(leadingJetPt, leadingJetEta, leadingJetPhi);
        //   cout<<"we found a leading jet:" <<leadingJetPt <<endl ;
          if(coll=="PbP"){
           if(TMath::Abs(leadingJetEta-0.465)<=1.) hist->leadjetpt->Fill(leadingJetPt, weight);
           }
          else {
           if(TMath::Abs(leadingJetEta+0.465)<=1.) hist->leadjetpt->Fill(leadingJetPt, weight);
         }
          for(int j = 0; j < nref ; j++) {

         if (fabs(jteta[j])>3.0) continue;
         if (rawpt[j]<22) continue;
         if(j==leadingJetIndex) continue ;
            double jet_pt= jtpt[j];
            double jet_phi= jtphi[j];
           double jet_eta = jteta[j];
          jet2_vec.SetPtEtaPhi(jet_pt, jet_eta, jet_phi);
      //    if(TMath::Abs(jet2_vec.DeltaPhi(jet1_vec))<dphicut) continue ;
         if (jtpt[j]>subleadingJetPt) {
            subleadingJetPt = jtpt[j];
            subleadingJetEta = jteta[j];
            subleadingJetPhi = jtphi[j];
            subleadingJetIndex = j;
         }
      } //jet loop for subleading search  
     }
       if(!IsMC){
          if(coll=="PPb" || coll=="PbP") subleadres*=C_rel->GetBinContent(C_rel->FindBin(subleadingJetEta)); 
        //   subleadingJetPt*=subleadres;
     }
     if(subleadingJetIndex!=-1){
         jet2_vec.SetPtEtaPhi(subleadingJetPt, subleadingJetEta, subleadingJetPhi);
          hist->leadjetptDPhi->Fill(leadingJetPt, TMath::Abs(jet2_vec.DeltaPhi(jet1_vec)), weight);
          hist->leadjetptratio->Fill(leadingJetPt, subleadingJetPt/leadingJetPt, weight);
         if(leadingJetPt>=subleadingJetPt && subleadingJetPt >subleadjetcut){
         double dijetM = TMath::Sqrt(2*leadingJetPt*subleadingJetPt*(TMath::CosH(leadingJetEta-subleadingJetEta)-TMath::Cos(leadingJetPhi-subleadingJetPhi)));
       double dijetEta = (leadingJetEta+subleadingJetEta)/2. ;
          hist->leadjetptDiEta->Fill(leadingJetPt, dijetEta, weight);
          hist->MjjDiEta->Fill(dijetM, dijetEta, weight);
         //  cout<<" subleading jet:" <<subleadingJetPt <<endl ;
          if(coll=="PbP"){
           if(TMath::Abs(subleadingJetEta-0.465)<=1.) hist->subjetpt->Fill(subleadingJetPt, weight);
           }
          else {
           if(TMath::Abs(subleadingJetEta+0.465)<=1.) hist->subjetpt->Fill(subleadingJetPt, weight);
         }
               int dEtaBin = -1. ;
          if(coll=="PPb"){
            if(TMath::Abs(dijetEta+0.465)<=1.) hist->dijetMass->Fill(dijetM, weight);
           }
          else if(coll=="PbP"){
             if(TMath::Abs(dijetEta-0.465)<=1.) hist->dijetMass->Fill(dijetM, weight);
          }
          else {
            if(TMath::Abs(dijetEta+0.465)<=1.) hist->dijetMass->Fill(dijetM, weight);
         }
         for(Int_t ieta = 0 ; ieta <netabin; ieta++){
            if(coll=="PPb"){
               if((dijetEta+0.465)>deta[ieta]&&(dijetEta+0.465)<=deta[ieta+1]) dEtaBin = ieta ;
             }
            else if(coll=="PbP"){
               if((dijetEta-0.465)>deta[ieta]&&(dijetEta-0.465)<=deta[ieta+1]) dEtaBin = ieta ;
             }
           else {
               if((dijetEta+0.465)>deta[ieta]&&(dijetEta+0.465)<=deta[ieta+1]) dEtaBin = ieta ;
            }
         } //assign the eta bin for jets
    //    if(dEtaBin!=-1){
    //    hist->dijetMassEtaBin[dEtaBin]->Fill(dijetM, weight);
    //    }
   //   } //if leading and subleading jet pt cuts satisfied

      //for data-MC comparison on the residual correction
       if(fabs(leadingJetEta)<1.3){
        double B =  2*(subleadingJetPt-leadingJetPt)/(subleadingJetPt+leadingJetPt);
        if((subleadingJetPt+leadingJetPt)/2.>100.)hist->ProbeEtaMeanB->Fill(subleadingJetEta, B);  
        if(fabs(subleadingJetEta)<1.3) 
           hist->AvgPtBalanceMid->Fill((subleadingJetPt+leadingJetPt)/2., 2*(subleadingJetPt-leadingJetPt)/(subleadingJetPt+leadingJetPt), weight);
      else 
           hist->AvgPtBalanceFor->Fill((subleadingJetPt+leadingJetPt)/2., 2*(subleadingJetPt-leadingJetPt)/(subleadingJetPt+leadingJetPt), weight);

      }
       if(fabs(subleadingJetEta)<1.3){
                double B =  2*(leadingJetPt-subleadingJetPt)/(subleadingJetPt+leadingJetPt);
        if((subleadingJetPt+leadingJetPt)/2.>100.)hist->ProbeEtaMeanB->Fill(leadingJetEta, B);
        if(fabs(leadingJetEta)<1.3)
           hist->AvgPtBalanceMid->Fill((subleadingJetPt+leadingJetPt)/2., 2*(leadingJetPt-subleadingJetPt)/(subleadingJetPt+leadingJetPt), weight);
      else
           hist->AvgPtBalanceFor->Fill((subleadingJetPt+leadingJetPt)/2., 2*(leadingJetPt-subleadingJetPt)/(subleadingJetPt+leadingJetPt), weight);

      }
    } //if found subleading jets 
 
        //for inclusive track analysis, without jet selection and requirement
       for(int itr = 0 ; itr < nTrk ; itr++){
                            double tr_pt = trkPt[itr];
                            double tr_eta = trkEta[itr];
                            if(TMath::Abs(tr_eta)>2.4) continue ;
                                     Int_t TrkEtaBin = -1 ;
                       for(Int_t ieta = 0 ; ieta <ntrketabin; ieta++){
                        if(coll=="PPb"){
                           if((tr_eta+0.465)>dtrketa[ieta]&&(tr_eta+0.465)<=dtrketa[ieta+1]) TrkEtaBin = ieta ;
                           }
                        else if(coll=="PbP"){
                            if((tr_eta-0.465)>dtrketa[ieta]&&(tr_eta-0.465)<=dtrketa[ieta+1]) TrkEtaBin = ieta ;
                           }
                        else
                      if((tr_eta+0.465)>dtrketa[ieta]&&(tr_eta+0.465)<=dtrketa[ieta+1]) TrkEtaBin = ieta ;
                     }
               double trkweight=1. ;
                if(doTrackCorrections){
                             trkweight = corr.getWeight(tr_pt,tr_eta,0.);
                            }
                        if(coll=="PPb"){
                           if(TMath::Abs(tr_eta+0.465)<=1.) hist->trkpt[ntrketabin]->Fill(tr_pt, weight*trkweight);
                         }
                       else if(coll=="PbP"){
                           if( TMath::Abs(tr_eta-0.465)<=1.) hist->trkpt[ntrketabin]->Fill(tr_pt, weight*trkweight);
                        }
                       else {
                         if(TMath::Abs(tr_eta+0.465)<=1.) hist->trkpt[ntrketabin]->Fill(tr_pt, weight*trkweight);
                       }
             if(TrkEtaBin!=-1) hist->trkpt[TrkEtaBin]->Fill(tr_pt, weight*trkweight);
            }  //! inclusive track loop

    // for jet-track analysis
   for(int j4i = 0; j4i < nref ; j4i++) {
    jet_vec.SetPtEtaPhi(0, 0, 0);
     track_vec.SetPtEtaPhi(0, 0, 0);
    double jet_pt= jtpt[j4i];
    double jet_eta = jteta[j4i];
    double jet_phi = jtphi[j4i];
    if (rawpt[j4i]<22) continue;
    if(TMath::Abs(jet_eta)<=3.){
     jetweight = 1;
       jet_vec.SetPtEtaPhi(jet_pt, jet_eta, jet_phi);
    // for track loop in each jet, do jet-track analysis
           for(int itr = 0 ; itr < nTrk ; itr++){
                            double tr_pt = trkPt[itr];
                            double tr_phi = trkPhi[itr];
                            double tr_eta = trkEta[itr];
                            if(TMath::Abs(tr_eta)>2.4) continue ;
                            Int_t TrkEtaBin = -1 ;
                              for(Int_t ieta = 0 ; ieta <ntrketabin; ieta++){
                        if(coll=="PPb"){
                           if((tr_eta+0.465)>dtrketa[ieta]&&(tr_eta+0.465)<=dtrketa[ieta+1]) TrkEtaBin = ieta ;
                           }
                        else if(coll=="PbP"){
                            if((tr_eta-0.465)>dtrketa[ieta]&&(tr_eta-0.465)<=dtrketa[ieta+1]) TrkEtaBin = ieta ;
                           }
                        else
                      if((tr_eta+0.465)>dtrketa[ieta]&&(tr_eta+0.465)<=dtrketa[ieta+1]) TrkEtaBin = ieta ;
                     }
            track_vec.SetPtEtaPhi(tr_pt, tr_eta, tr_phi);
                            double trkweight=1. ;
                   double  dr = jet_vec.DeltaR(track_vec);
                   if(dr>conesize) continue ;
                             if(doTrackCorrections){
                             trkweight = corr.getWeight(tr_pt,tr_eta,0.);
                            }
                         if(coll=="PPb"){
                            if( TMath::Abs(tr_eta+0.465)<=1.) hist->jettrkpt[ntrketabin]->Fill(tr_pt, weight*trkweight);
                         }
                        else if(coll=="PbP") {
                           if(TMath::Abs(tr_eta-0.465)<=1.) hist->jettrkpt[ntrketabin]->Fill(tr_pt, weight*trkweight);
                         }
                        else {
                           if(TMath::Abs(tr_eta+0.465)<=1.) hist->jettrkpt[ntrketabin]->Fill(tr_pt, weight*trkweight);
                        }
                        if(TrkEtaBin!=-1) hist->jettrkpt[TrkEtaBin]->Fill(tr_pt, weight*trkweight);
           } //end of track loop
        } //jet kinematics
    } //jet loop 
*/
        hist->JetAbove20->SetBinContent(1,Nevents*weight);
        hist->JetAbove40->SetBinContent(1,Nevt_40_60*weight);
        hist->JetAbove60->SetBinContent(1,Nevt_60_75*weight);
        hist->JetAbove75->SetBinContent(1,Nevt_75_95*weight);
        hist->JetAbove95->SetBinContent(1,Nevt_95_120*weight);
        hist->JetAbove120->SetBinContent(1,Nevt_120*weight);

    }
    
    hist->Write();
    std::cout << "working done\n";
    tFile->Close() ;
    hist->Delete();
    tJet = NULL;
}
