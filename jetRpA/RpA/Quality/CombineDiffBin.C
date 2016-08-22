//
//  CombineDiffBin.C
//  
//
//  Created by Yaxian Mao on 3/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include "TLine.h"
#include <TF1.h>
#include <TMath.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include <TMath.h>
#include <TSystem.h>
#include "THStack.h"
#include <iostream>
using namespace std;



const Int_t kw = 800;
const Int_t kh = 800;
enum Data_t {kPP, kPPb, kPbP, kHI} ;
enum Method_t {kRC, kER} ;
enum Level_t {kReco, kGen} ;
const Double_t jetPtBin[]={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300, 500};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;

const int ntrkptbin = 4 ;

void CombineDiffBin(Data_t data = kPPb, Method_t met=kER, Level_t lev=kReco, double trackcut = 0.){

    TString bkg ;
    if(met==kER) {
        const double etalimit = 0.3 ; 
         bkg = "ER";
    }
    else {
        bkg = "RC";
        const double etalimit = 2.0 ; 
    }
    TString coll ;
    if(data==kPP) coll="PP";
    else if(data==kPPb) coll="PPb";
    else if(data==kPbP) coll="PbP";
    else coll = "HI" ;
    TString effTab = "2013HistIterTrkCorrtest" ;  //"HistIterTrkCorrtestFilterCut"; //  "Trk" ; "HistIterTrkCorrtest" ;
    TString trkSel ="Sim" ;
    
    TString fileName ;
    if(data==kPP){
        ////for pp MC sample
//        const int Nfile = 7 ;
//        const int pthat[]={50, 80,120, 170,200, 250, 300};
//        const double cs[]={1.021e-3, 9.913e-5, 1.128e-5, 1.470e-6, 5.310e-7, 1.192e-7,3.176e-8};
//////*******     //if use data sample under (no simTrk info):
/////        /mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged
// //       const double Nevt[]={131386, 128738, 59039, 64384, 56432, 37976,52929};
//       // //*******     //if use data sampe under (with simTrk info):
//        ///        /net/hisrv0001/home/zhukova/scratch/HIHighPt/forest/pthat*        
//        const double Nevt[]={131386, 128738, 59039, 64384, 56432, 37976,18189};
//        //2011 MC for pp with HI tracking
//        const int Nfile = 6 ;
//        const int pthat[]={80,120, 170,200, 250, 300};
//        const double cs[]={9.913e-5, 1.128e-5, 1.470e-6, 5.310e-7, 1.192e-7,3.176e-8};
//        const double Nevt[]={128738, 59039, 64384, 56432, 37976,18189};

//      //2013 pp @2.76 TeV MC sample for 2.76TeV
//        const int Nfile = 6 ;
//        const int pthat[]={50, 80,120, 170,220, 280};
//        const double cs[]={1.025E-03, 9.865E-05, 1.129E-05,1.465E-06,2.837E-07,5.323E-08};
//        //2013 pp MC sample with HI tracking
//      //  const double Nevt[]={161000, 149500, 159000, 182500, 193000,157000};
//
//        //2013 pp @2.76TeV MC sample with pp tracking
//        const double Nevt[]={161000, 149500, 159000, 182500, 193000,157000};
        
        //2013 pp @5.02 TeV MC sample for 5.02 TeV
        const int Nfile = 8 ;
        const int pthat[]={15, 30, 50, 80,120, 170,220, 280};
        const double cs[]={5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07};
        //2013 pp @ 5.02 TeV MC sample with pp tracking
        const double Nevt[]={197175, 199375, 182175, 199375, 164900, 200975, 180200,194375};

        const int nbin = 1 ;

    }
    else if(data==kPPb){
        //2013 ppb MC sample
        const int Nfile = 9 ;
        const int pthat[]={15, 30, 50, 80,120, 170,220, 280, 370};
        const double cs[]={5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07,1.088E-07};
        //2013 ppb MC sample with pp tracking
        const double Nevt[]={332075, 332075, 323475, 308275, 212875, 273275, 315275,208875, 135100};
        const int nbin = 1 ;
 
    }
    else if(data==kPbP){
        //2013 ppb MC sample
        const int Nfile = 7 ;
        const int pthat[]={50, 80,120, 170,220, 280, 370};
        const double cs[]={3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07,1.088E-07};
        //2013 pbp MC sample with pp tracking
        const double Nevt[]={362400, 358000, 306400, 360800, 352000, 315200, 338000};
        
        const int nbin = 1 ;
        
    }
    else {
        //for PbPb MC sample
        const int Nfile = 7 ;
        const int pthat[]={80,100, 120, 170,200, 250, 300};
        const double cs[]={9.913e-5, 3.069e-5, 1.128e-5, 1.470e-6, 5.310e-7, 1.192e-7,3.176e-8};
        const double Nevt[]={146162, 593463, 141263, 62944, 54952, 37856, 53009};

//        const int Nfile = 3 ;
//        const int pthat[]={80,100, 170};
//        const double cs[]={9.913e-5, 3.069e-5,1.470e-6};
//        const double Nevt[]={146162, 593463, 62944};
        const double lum = 1.5e5 ;// PbPb data luminosity=129mub^-1, convert to mb^-1
        const int nbin = 5 ;
    }
    TFile * f[Nfile] ;
    TString name[Nfile] ;
    TList * list[Nfile] ;
    double weight[Nfile] ;
//    const char * kDir = "/Users/ymao/group/CMS/trkeff" ;     
//    const char * kDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/FullStatPbPb/RebinJetPt/TrkCorr4CentBin" ;     
//    const char * kDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/pPb/GplusCut/CMEtaBin" ;
 //   const char * kDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/NoMatchGen/NoTrkAlgoCut/MyselfTrkEff/AfterCWR" ;
  //  const char * kDir = "/Users/ymao/group/CMS/anaOutputs/pPb/RpPb" ;
    const char * kDir = "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual" ;
    int N ;

    for(int ifile=0; ifile<Nfile; ifile++){ 
//       name[ifile]=Form("MCPP_ChargedJetPt80_40Trk4_BkgEta5Phi0_CenBin1_Nrbin6_Pthat%dmerged_hiForest.root",pthat[ifile]);
 //     name[ifile]=Form("MCPP_ChargedJetPt80_40Trk1Eta2_RCbkgEtaBound20_CenBin1_Nrbin6_Pthat%dpp276Dijet%d_merged.root",pthat[ifile], pthat[ifile]);
        if(data==kPP){
            if(nbin==1){
                if(lev==kGen)
                    name[ifile]=Form("MC%s_RefJetPt100_%sTrk%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1pthat%d_pt%d_JEC_ppHiIterativeTrack_P01_prod24_v84_merged_forest_0.root",coll.Data(),trkSel.Data(), trackcut, bkg.Data(),nbin,pthat[ifile],pthat[ifile]);
//                    name[ifile]=Form("MC%s_RefJetPt100_%sTrk%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1pthat%d_mergedFile.root",coll.Data(),trkSel.Data(), trackcut, bkg.Data(),nbin,pthat[ifile]);
                else  
          //          name[ifile]=Form("pp%d.root",pthat[ifile]);
        //            name[ifile]=Form("MC%s2011_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%dpthat%d_mergedFile.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin, pthat[ifile]);
     //               name[ifile]=Form("MC%s_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%dpthat%d_pt%d_JEC_ppHiIterativeTrack_P01_prod24_v84_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin, pthat[ifile],pthat[ifile]);
    //            name[ifile]=Form("MC%s_Ak3PFFliterPtGluonJetPt50_GenCharge%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%dpthat%d_pt%d_pp2013_P01_prod22_v81_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin, pthat[ifile], pthat[ifile]);
                    
                    //pp @5.02TeV MC
//                    name[ifile]=Form("MC%s_Ak3PFIncJetPt50_GenCharge%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%dpthat%d_pt%d_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin, pthat[ifile], pthat[ifile]);
//                name[ifile]=Form("MC%s_AkPu3CaloJetNoResidualGenCharge%s%.fEtaBin9_pt%d_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut, pthat[ifile]);
//                name[ifile]=Form("MC%s_AllAkPu3PFJetGenCharge%s%.fEtaCut10_pt%d_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut, pthat[ifile]);
//             name[ifile]=Form("MC%s_AkPu3PFJetNoResidualGenCharge%s%.fEtaBin7_HFsumEta4Bin1_pt%d_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut ,pthat[ifile]);
	//             name[ifile]=Form("MC%s_AkPu3PFJetNoResidualGenCharge%s%.fEtaCut10_pt%d_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut ,pthat[ifile]);
	             name[ifile]=Form("MC%s_ak3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt%d_HP04_hiforest77_hiSignal.root",coll.Data(),pthat[ifile]);
            }
            
            else {
                if(lev==kReco)  
                 //   name[ifile]=Form("MC%s_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin%dpthat%d_mergedFile.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin, pthat[ifile]);
                 name[ifile]=Form("MC%s_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%dpthat%d_pt%d_JEC_ppHiIterativeTrack_P01_prod24_v84_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin, pthat[ifile],pthat[ifile]);
            }
            
        }
        else if(data==kPPb){
 //           name[ifile]=Form("MC%s_Ak3PFIncJetPt50_GenCharge%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_pt%d_HP04_prod16_v77_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin, pthat[ifile]);
//            name[ifile]=Form("MC%s_AkPu3PFJetNoEtaShiftGenCharge%s%.fEtaCut10pthat%d_pt%d_HP04_prod16_v77_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut,pthat[ifile], pthat[ifile]);
           
//if(pthat[ifile]==280) 
//name[ifile]=Form("MC%s_AllAkPu3PFJetNoResidualGenCharge%s%.fEtaCut10_newbin_pt%d_HP04_prod16_v77_merged_forest_0_v2.root",coll.Data(),effTab.Data(), trackcut,pthat[ifile], pthat[ifile]);
//else 
name[ifile]=Form("MC%s_ak3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt%d_HP04_prod16_v77_merged_forest_0.root", coll.Data(),pthat[ifile]);
 
        }
        else if(data==kPbP){
            name[ifile]=Form("MC%s_Ak3PFIncJetPt50_GenCharge%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_pt%d_HP05_prod24_v84_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin, pthat[ifile]);
            
        }
        else{
            if(pthat[ifile]==50||pthat[ifile]==80||pthat[ifile]==100||pthat[ifile]==170){
                if(lev==kGen) 
                    name[ifile]=Form("MC%s_RefJetPt100_%sTrk%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1pthat%d_Dijet%d_HydjetDrum_v27_mergedV1.root",coll.Data(),trkSel.Data(), trackcut, bkg.Data(),nbin,pthat[ifile], pthat[ifile]);        
                else 
                    name[ifile]=Form("MC%s_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_Dijet%d_HydjetDrum_v27_mergedV1.root",coll.Data(), effTab.Data(), trackcut, bkg.Data(),nbin, ntrkptbin, pthat[ifile]);
            }
            else {
                if(lev==kGen)
                    name[ifile]=Form("MC%s_RefJetPt100_%sTrk%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1pthat%d_Dijet%d_HydjetDrum_v28_mergedV1.root",coll.Data(),trkSel.Data(), trackcut, bkg.Data(),nbin,pthat[ifile], pthat[ifile]);   
                else 
                    name[ifile]=Form("MC%s_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_Dijet%d_HydjetDrum_v28_mergedV1.root",coll.Data(), effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin, pthat[ifile]);
            }
    //        name[ifile]=Form("pbpb%d.root",pthat[ifile]);
        }
        
        f[ifile] = TFile::Open(Form("%s/%s", kDir, name[ifile].Data()), "readonly");
        
      list[ifile] =(TList*)f[ifile]->GetListOfKeys();
      if(ifile==0){
        N = list[ifile]->GetEntries();
cout<<"N = "<<N<<endl;
      }
        if(ifile==Nfile-1)weight[ifile]=cs[ifile]/Nevt[ifile];
        else weight[ifile]=(cs[ifile]-cs[ifile+1])/Nevt[ifile];
        
        cout <<"coll = " << coll << " ifile = " << ifile << " weight = " << weight[ifile] <<endl ;
    }
      const int Nhist = N ;
      TH1 * h[Nfile][Nhist];
      TH1 * htot[Nhist];
    for(int ifile=0; ifile<Nfile; ifile++){ 
      for(int ihist = 0 ; ihist < Nhist ; ihist++){
         h[ifile][ihist] = (TH1*)f[ifile]->Get(list[ifile]->At(ihist)->GetName()) ;
        //h[ifile][ihist]->Sumw2();
        h[ifile][ihist]->Scale(weight[ifile]);
      }
    }
    TString outname ;
//    TString outname ="merged_MCPP_ChargedJetPt80_40Trk4_BkgEta5Phi0_CenBin1_Nrbin6_merged_hiForest.root" ;
//    TString outname ="merged_MCPP_ChargedJetPt80_40Trk1Eta2_RCbkgEtaBound20_CenBin1_Nrbin6_pp276Dijet.root" ;
    switch(data){
        case kPP:
            if(nbin==1){
                if(lev==kGen)
                    outname=Form("mergedCSdiff_MC%s_RefJetPt100_%sTrk%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_ppHiIterativeTrack_P01_prod24_v84_merged_forest_0.root",coll.Data(),trkSel.Data(), trackcut, bkg.Data(),nbin);
//                    outname=Form("mergedCSdiff_MC%s_RefJetPt100_%sTrk%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_mergedFile.root",coll.Data(),trkSel.Data(), trackcut, bkg.Data(),nbin);
                else  
           //         outname =Form("merged_Pelin_PYTHIA_output.root") ;
            //        outname =Form("mergedCSdiff_MC%s2011_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin%d_pp276Dijet_merged.root", coll.Data(),effTab.Data(), trackcut,bkg.Data(), nbin, ntrkptbin) ;
          //          outname =Form("mergedCSdiff_MC%s2013_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin%d_pp276Dijet_ppHiIterativeTrack_merged.root", coll.Data(),effTab.Data(), trackcut,bkg.Data(), nbin, ntrkptbin) ;
       //             outname =Form("mergedCSdiff_MC%s_Ak3PFFliterPtGluonJetPt50_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin%d_pp2013_P01_prod22_v81_merged.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin);
                    //pp @5.02 TeV MC 
//                    outname =Form("mergedCSdiff_MC%s_Ak3PFIncJetPt50_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_HP04_hiforest77_hiSignal5TeV.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin);
 //               outname =Form("mergedCSdiff_MC%s_AkPu3CaloIncJetNoResidual_%s%.fEtaBin9_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut);
 //               outname=Form("mergedCSdiff_MC%s_AllAkPu3PFIncJetGenCharge%s%.fEtaCut10_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut);
 //           outname=Form("mergedCSdiff_MC%s_AkPu3PFIncJetNoResidual%s%.fEtaBin7_HFsumEta4Bin1_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut);
    // 	outname= Form("mergedCSdiff_MC%s_AkPu3PFIncJetNoResidual%s%.fEtaCut10_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut);
    	outname= Form("mergedCSdiff_MC%s_ak3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_hiforest77_hiSignal.root",coll.Data());
            }
            else {
              //  outname =Form("mergedCSdiff_MC%s_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin%d_pp276Dijet_merged.root", coll.Data(),effTab.Data(), trackcut,bkg.Data(), nbin, ntrkptbin) ;
               outname =Form("mergedCSdiff_MC%s_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin%d_pp276Dijet_ppHiIterativeTrack_merged.root", coll.Data(),effTab.Data(), trackcut,bkg.Data(), nbin, ntrkptbin) ;

            }
            break ;
        case kPPb:
   //         outname =Form("mergedCSdiff_MC%s_Ak3PFIncJetPt50_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_HP04_prod16_v77_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin);
      //      outname =Form("mergedCSdiff_MC%s_AkPu3PFIncJetNoEtaShift_%s%.fEtaCut10_HP04_prod16_v77_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut);
//            outname =Form("mergedCSdiff_MC%s_AllAkPu3PFIncJetGenCharge%s%.fEtaCut10_HP04_prod16_v77_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut);
            outname =Form("mergedCSdiff_MC%s_ak3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root",coll.Data());
            break ;
        case kPbP:
            outname =Form("mergedCSdiff_MC%s_Ak3PFIncJetPt50_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_HP05_prod24_v84_merged_forest_0.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin,ntrkptbin);
            break ;
        case kHI:
            if(lev==kGen)
                outname =Form("mergedCSdiff_MC%s_RefJetPt100_%sTrk%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_Dijet_HydjetDrum_v27_mergedV1.root",coll.Data(),trkSel.Data(), trackcut, bkg.Data(),nbin);
            else 
                outname =Form("mergedCSdiff_MC%s_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_%sbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_Dijet_HydjetDrum_v27v28.root",coll.Data(),effTab.Data(), trackcut, bkg.Data(),nbin, ntrkptbin) ;
        //        outname =Form("merged_Pelin_PYTHIAHydjet_output.root") ;
            break ;
    }

    TFile * fmerge = new TFile(Form("%s/%s", kDir,outname.Data()), "RECREATE");

   //for display
    //if wants to display
    TCanvas  * c1 = new TCanvas("c1", "c1", kw, kh);
    c1->SetFillColor(0);
    c1->SetBorderSize(0);
    c1->SetFrameBorderMode(0); 
    gStyle->SetOptStat(0);    
    
    TLegend *t1=new TLegend(0.65,0.58,0.85,0.88);
    t1->SetFillColor(0);
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextFont(63);
    t1->SetTextSize(17);
    
    TH1D * hFrame = new TH1D("hFrame","frame",1000,0.,1000.);    
 //   hFrame-> SetAxisRange(4.,500.,"X");
    hFrame-> SetAxisRange(30,1000.,"X");
    if(data==kHI) hFrame->SetAxisRange(1.e-10, 10., "Y") ;
  //  else hFrame->SetAxisRange(1.e-10, 1., "Y") ;
    else hFrame->SetAxisRange(1.e-12, 1.e-2, "Y") ;
    hFrame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    
    const Int_t dis = 9 ;

    for(int ihist = 0 ; ihist < Nhist ; ihist++){

            for(int ifile=0; ifile<Nfile; ifile++){ 
                if(ifile==0)
                    htot[ihist] = (TH1*)h[ifile][ihist]->Clone(Form("%s",list[ifile]->At(ihist)->GetName()));
                else 
                    htot[ihist]->Add(h[ifile][ihist]);
                
            }
        
        if(htot[ihist]->GetEntries()>0){
      //      cout <<"name =" << htot[ihist]->GetName()<<endl ;

//            if(Contains(htot[ihist]->GetName(),"jetpt_")){
//                cout <<"name =" << htot[ihist]->GetName()<<endl ;
//            }
            
   //         htot[ihist]->Scale(lum);
            htot[ihist]->Write();
        }

    

 //pick only one histogram for display
        if(ihist==dis){
            c1->cd() ;
            c1->SetLogy() ;
            cout <<list[0]->At(dis)->GetName()<<endl ;
            
            hFrame->DrawCopy();
//            for(Int_t i = 1 ; i < =nJetPtBin ; i++){
//                htot[dis]->SetBinContent(i,htot[dis]->GetBinContent(i)/htot[dis]->GetBinWidth(i));  
//                htot[dis]->SetBinError(i,htot[dis]->GetBinError(i)/htot[dis]->GetBinWidth(i)); 
//            }

//            htot[dis]->Rebin(10);
//            htot[dis]->Scale(1./10.);
            normalizeByBinWidth(htot[dis]);
            htot[dis]->SetLineColor(1);
            htot[dis]->SetMarkerColor(1);
            htot[dis]->SetMarkerStyle(20);
            htot[dis]->SetMarkerSize(1.5);
            htot[dis]->SetTitle(Form("%s",list[0]->At(dis)->GetName()));
            htot[dis]->DrawCopy("same PE");
            t1->AddEntry(htot[dis], "merged", "P");
            for(int ifile=0; ifile<Nfile; ifile++){ 
//                for(Int_t i = 1 ; i < =nJetPtBin ; i++){
//                    h[ifile][dis]->SetBinContent(i,h[ifile][dis]->GetBinContent(i)/h[ifile][dis]->GetBinWidth(i));  
//                    h[ifile][dis]->SetBinError(i,h[ifile][dis]->GetBinError(i)/h[ifile][dis]->GetBinWidth(i)); 
//                }

//                h[ifile][dis]->Rebin(10);
//                h[ifile][dis]->Scale(1./10.);
//                h[ifile][dis]->Scale(lum);
                normalizeByBinWidth(h[ifile][dis]);
                h[ifile][dis]->SetLineColor(ifile+2);
                h[ifile][dis]->SetMarkerColor(ifile+2);
                h[ifile][dis]->SetMarkerStyle(21+ifile);
                h[ifile][dis]->DrawCopy("same PE");
                t1->AddEntry(h[ifile][dis], Form("\#hat{p}_{T}=%dGeV/c",pthat[ifile]), "P");
            }
            
            t1->Draw("same");
            c1->Update();
 
        }   
    }
    fmerge->Close() ;
    
}
TH1 * normalizeByBinWidth(TH1 *histo) {
    for(int i = 1; i <= histo->GetNbinsX(); i++) {
        float content = histo->GetBinContent(i);
        float error = histo->GetBinError(i);
        histo->SetBinContent(i,content/histo->GetBinWidth(i));
        histo->SetBinError(i,error/histo->GetBinWidth(i));
    }
    return histo ;
}
