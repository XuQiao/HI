#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TCut.h"
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/Correction.C"

using namespace std;

// ******* GLOBAL DECLARATIONS **********
const int QCDpthatBins = 8; //MC samples
const int sigma = 70 ;  //units in mb
const double Ncoll = 6.9 ; //for pPb Ncoll
//***************************************

const double deta[]={-2.2, -1.2, -0.7, -0.3, 0.3, 0.7,1.2,2.2} ;
const int netabin = sizeof(deta)/sizeof(Double_t)-1 ;
const Double_t jetPtBin[]={3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;
const int pthatbin[10] = {15,30,50,80,120,170,220,280,370, 9999};
const double wght[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 1.088E-07, 0};
const Double_t TrkBin[]={0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 138., 155.2, 172.4, 189.6, 206.8};
const int nTrkBin = sizeof(TrkBin)/sizeof(Double_t)-1 ;
const int ntrketabin = 7 ;
const double dtrketa[]={-1.8, -1.3, -0.8, -0.3, 0.3, 0.8,1.3,1.8} ;

const float pPbLumi=21.1; //unit in nb, need to convert to mb for cross sections
const bool SavePlot=kFALSE;
const bool IsMC = kFALSE ;
const TString coll = "PPb";
const TString algo = "akPu3PF" ;
TString TrigName = "Jet20" ; //"JetMB","Jet20","Jet40","Jet60","Jet80","Jet100"
double weight;


class JetData
{
public:
    JetData(char *fileName, char *jetTree) {	
        cout <<"Open "<<fileName<<endl;
        tFile = new TFile(fileName,"read");
        tJet = (TTree*)tFile->Get(jetTree);
	nevents=tJet->GetEntries();
        tJet->SetBranchAddress("nref",&nref);
        tJet->SetBranchAddress("rawpt" , &rawpt );
        tJet->SetBranchAddress("jtpt" , &jtpt );
        tJet->SetBranchAddress("refpt", &refpt);
        tJet->SetBranchAddress("refeta", &refeta);
        tJet->SetBranchAddress("jteta", &jteta);

        tJet->SetBranchAddress("nTrk",&nTrk);
        tJet->SetBranchAddress("trkPt",&trkPt);
        tJet->SetBranchAddress("trkPhi",&trkPhi);
        tJet->SetBranchAddress("trkEta",&trkEta);
        
        tJet->SetBranchAddress("pthat",&pthat);
        tJet->SetBranchAddress("weight",&weight);
        tJet->SetBranchAddress("vz",&vz);
        tJet->SetBranchAddress("run",&run);
        tJet->SetBranchAddress("bin",&bin);
        tJet->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
        tJet->SetBranchAddress("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_noJetID_v1);
        tJet->SetBranchAddress("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_noJetID_v1);
        tJet->SetBranchAddress("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_noJetID_v1);
        tJet->SetBranchAddress("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_noJetID_v1);
        tJet->SetBranchAddress("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_noJetID_v1);
        tJet->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
    };
    
    TFile *tFile;
    TTree *tJet;
int nevents;
    double rawpt;
    double jtpt;
    double refpt;
    double refeta;
    double jteta;
    double weight;
    double pthat;
    int nref;
    int nTrk;
    double trkPt;
    double trkEta;
    double trkPhi;
    int bin;
    int run;
    double vz;
    double hiHFplusEta4;
    int HLT_PAJet20_noJetID_v1;
    int HLT_PAJet40_noJetID_v1;
    int HLT_PAJet60_noJetID_v1;
    int HLT_PAJet80_noJetID_v1;
    int HLT_PAJet100_noJetID_v1;
    int HLT_PAZeroBiasPixel_SingleTrack_v1;

};

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
    
    TH1F * jetpt[netabin+1];

//    TH1F * refjetpt;
//    TH2F * refptJES;
    
};
Histos::Histos()
{
    NEvents = new TH1D("Nevents_0-100%","Nevents_0-100%",100,0,2.);
    NEvents -> Sumw2();
    
    JetAbove20 = new TH1D("NeventsWithJetBelow40_0-100%","NeventsWithJetBelow40_0-100%",1,0,1);
    JetAbove20 -> Sumw2();
    
    JetAbove40 = new TH1D("NeventsWithJetAbove40_0-100%","NeventsWithJetAbove40_0-100%",1,0,1);
    JetAbove60 = new TH1D("NeventsWithJetAbove60_0-100%","NeventsWithJetAbove60_0-100%",1,0,1);
    JetAbove75 = new TH1D("NevntsWithJetAbove75_0-100%","NeventsWithJetAbove75_0-100%",1,0,1);
    JetAbove95 = new TH1D("NevntsWithJetAbove95_0-100%","NeventsWithJetAbove95_0-100%",1,0,1);
    JetAbove120 = new TH1D("NevntsWithJetAbove120_0-100%","NeventsWithJetAbove120_0-100%",1,0,1);
    
    for(int ieta=0; ieta <=netabin; ieta++){
        if(ieta==netabin)
            jetpt[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",-1.0*10, 1.0*10),Form("jetptEtaBin%.f_%.f",-1.0*10, 1.0*10),1000,0.,1000.);
            else
                jetpt[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
        
        jetpt[ieta]->Sumw2();
    }
}

void Histos::Write()
{
    TString dataType;
    TString out_name;
    if(IsMC) dataType = "MC";
    else dataType = "DATA";
    
    out_name=Form("%s%s%s%sJetSpectra.root",dataType.Data(),coll.Data(),TrigName.Data(),algo.Data());
    
    TFile *out_file = new TFile(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/NewTree/%s",out_name.Data()),"RECREATE");
    
    NEvents->Write();
    JetAbove20->Write();
    JetAbove40->Write();
    JetAbove60->Write();
    JetAbove75->Write();
    JetAbove95->Write();
    JetAbove120->Write();
   // jetpt->Write();
    for(int ieta=0; ieta<=netabin; ieta++) jetpt[ieta]->Write();
    
    out_file->Close();
    cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void MakeSpectFromTree(int isMC = 0, double recoJetPtCut = 30.,double trackMaxPtCut = 0)
{
  
  gStyle->SetErrorX(0.);
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetOptLogz(1);
  gStyle->SetPadRightMargin(0.13);	
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);


  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

    
  // input files
    std::cout<<"start working\n";
    cout<<"analyze trig = "<<TrigName<<endl;
    cout <<"MC = " << isMC <<endl;
    
    // Initialize Histograms
    Histos *hist = new Histos();
    
    if (IsMC) cout<<"This is a MC analysis"<<endl;
    else cout<< "This is a data analysis"<<endl;
    

      char *fileName= NULL;
    
    if(IsMC==kTRUE)
        fileName = (char*)Form("%sMCKurtForestNoBoost_ppReco_%s_QCDjetTrig_JetPt0noIPupperCut.root", coll.Data(), algo.Data());
    else {
        if(TrigName=="MB")
            fileName = (char*)Form("%sdata_ppReco_%s_AlljetTrigKurtTrCombFile0_400_JetPt0noIPupperCut.root", coll.Data(), algo.Data());
        else
            fileName = (char*)Form("%sdata_ppReco_%s_AlljetTrigKurtTrCombFile0_400_JetPt0noIPupperCut.root", coll.Data(), algo.Data());
    }
    

//    TTree *nt = (TTree*)f->Get("nt");
    JetData *tree = new JetData(fileName,(char*)"nt"); // pp data

    Int_t Nevents = 0;
    Int_t Nevt_40_60 = 0;
    Int_t Nevt_60_75 = 0;
    Int_t Nevt_75_95 = 0;
    Int_t Nevt_95_120 = 0;
    Int_t Nevt_120 = 0;

    TString CorrfileName;
    CorrfileName="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/pPbHijing_EffCorr_VRw.root";
    Corrector corr_Qiao(CorrfileName);
    

    Int_t nentries = tree->nevents;
    
    for(int i=0; i<nentries; i++){
        tree->GetEntry(i);
        
        weight =1.;
        
        double HFbin = hiHFplusEta4;
        
        if(IsMC==kFALSE){
            int jetTr2;
            if(TrigName=="Jet20") jetTr2 = tree->HLT_PAJet20_noJetID_v1;
            else if (TrigName=="Jet40") jetTr2 = tree->HLT_PAJet40_noJetID_v1;
            else if (TrigName=="Jet60") jetTr2 = tree->HLT_PAJet60_noJetID_v1;
            else if (TrigName=="Jet80") jetTr2 = tree->HLT_PAJet80_noJetID_v1;
            else if (TrigName=="Jet100")jetTr2 = tree->HLT_PAJet100_noJetID_v1;
            else jetTr2 = tree->HLT_PAZeroBiasPixel_SingleTrack_v1;
            if(!jetTr2) continue;
        }
        if(IsMC==kFALSE) {
            if(tree->run>211256) continue ; //only for pPb runs at moment
        }
        
        double leadingJet = -999.;
        for(int j4i = 0; j4i < tree->nref; j4i++) {
            if(tree->jtpt[j4i] > leadingJet ) leadingJet = jtpt[j4i]; //find leading jet
        }
  
        bool jetAbove = false;
        bool jetAbove40 = false;
        bool jetAbove60 = false;
        bool jetAbove75 = false;
        bool jetAbove95 = false;
        bool jetAbove120 = false;
        
        if(leadingJet>40.) jetAbove40 = true;
        if(leadingJet>60.) jetAbove60 = true;
        if(leadingJet>75.) jetAbove75 = true;
        if(leadingJet>95.) jetAbove95 = true;
        if(leadingJet>120.) jetAbove120 = true;
        
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
        
        if(!jetAbove) continue;//jet trigger cut
        
        if(i%1000==1)    cout<<"analyzing "<< i <<" th event"<<endl;
        
        if(IsMC==kFALSE)
            weight*=corr_Qiao.getEventWeightHFPlus4bak(HFbin,kTRUE);//efficiency weight
        else
          weight*=tree->weight; // weight for MC
        
        Histos->NEvents->Fill(1,weight);

        for(int j4i = 0; j4i < tree->nref; j4i++){
            double jet_pt = tree->jtpt[j4i];
            double jet_eta = tree->jteta[j4i];
            int dEtaBin = -1.;
            if(coll=="PPb"){
                if(TMath::Abs(jet_eta+0.465)<=1.) Histos->jetpt[netabin]->Fill(jet_pt,weight);
            }
            
            for(int ieta =0; ieta<netabin;ieta++){
                if(coll=="PPb"){
                    if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
                }
            }//assign the eta bin for jets
            if(dEtaBin!=-1) Histos->jetpt[dEtaBin]->Fill(jet_pt,weight);
        }//jet loop
        
        Histos->JetAbove20->SetBinContent(1,Nevents*weight);
        Histos->JetAbove40->SetBinContent(1,Nevt_40_60*weight);
        Histos->JetAbove60->SetBinContent(1,Nevt_60_75*weight);
        Histos->JetAbove75->SetBinContent(1,Nevt_75_95*weight);
        Histos->JetAbove95->SetBinContent(1,Nevt_95_120*weight);
        Histos->JetAbove120->SetBinContent(1,Nevt_120*weight);

    }
    
    Histos->Write();
    std::cout << "working done\n";
}
