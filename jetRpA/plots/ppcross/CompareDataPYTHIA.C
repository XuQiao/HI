/*
 *  Perform data && MC (PYTHIA Z2) comparison for 7TeV pp collisions, data points reading from HepData below
 *  http://hepdata.cedar.ac.uk/view/ins902309
 *  http://hepdata.cedar.ac.uk/view/ins1208923
 *
 *  Created by Yaxian Mao on 04/12/14.
 *  Copyright 2014 Vanderbilt University (US). All rights reserved.
 *
 */

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
#include "TProfile.h"
#include "TGraphErrors.h"
#include <iostream>
#include <iostream.h>
#include <stdlib.h>
#include <fstream.h>
#include <iomanip>
#include "TGraphAsymmErrors.h"
using namespace std;

bool SavePlot = kFALSE ;
const bool SaveFile = kFALSE ;

const double ptboundary[]={846, 846, 846, 800, 570};
//const double ptboundary[]={600, 550, 500, 450, 300};
const Double_t jetPtBin[]={56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790, 846};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;

void CompareDataPYTHIA(TString sys="")
{
    const int nFile = 500 ;
  TString InputDir = "/Users/ymao/group/CMS/anaOutputs/pPb/RpPb";
    TFile *fMC = new TFile(Form("%s/AnaGENJetR357_7000GeV_Apr15_Z2Combined.root",InputDir.Data()));//7.00TeV QCD-11-004 binning
    TFile *fpowhegR5 = new TFile(Form("/Users/ymao/group/CMS/macros/ppRef/POWHEG/TuneZ2_POWHEG_CT10NLO_7TeVR5_supfactor250_ktmin20_%dfile.root", nFile));//7.00TeV QCD-11-004 binning
    TFile *fpowhegR7 = new TFile(Form("/Users/ymao/group/CMS/macros/ppRef/POWHEG/TuneZ2_POWHEG_CT10NLO_7TeVR7_supfactor250_ktmin20_%dfile.root", nFile));//7.00TeV QCD-11-004 binning

    TFile *fpowheg2R5 = new TFile(Form("/Users/ymao/group/CMS/macros/ppRef/POWHEG/TuneZ2_POWHEG_CT10NLO_5TeVR5_supfactor250_ktmin20_%dfile.root", nFile));//7.00TeV QCD-11-004 binning
    TFile *fpowheg2R7 = new TFile(Form("/Users/ymao/group/CMS/macros/ppRef/POWHEG/TuneZ2_POWHEG_CT10NLO_5TeVR7_supfactor250_ktmin20_%dfile.root", nFile));//7.00TeV QCD-11-004 binning
    
    const double deta[]={0., 0.5, 1.0, 1.5, 2.0, 2.5} ;
    int nptpoint[]={33, 30, 27, 24, 19};
//    const double deta[]={0., 0.5, 1.0, 1.5, 2.0} ;
//    int nptpoint[]={33, 30, 27, 24};
    
    TString plotsdir = "/Users/ymao/group/CMS/plots/pA/JetRpA/HIN-14-001/AfterQM14/NewOfficialMC";

    int netabin = sizeof(deta)/sizeof(Double_t)-1 ;
    const int N  = netabin ;
    
    char*R5dirname[] = {
        "ak5GenJetSpectrum_QCD10001_00_05",
        "ak5GenJetSpectrum_QCD10001_05_10",
        "ak5GenJetSpectrum_QCD10001_10_15",
        "ak5GenJetSpectrum_QCD10001_15_20",
        "ak5GenJetSpectrum_QCD10001_20_25"
         };
    char*R7dirname[] = {
        "ak7GenJetSpectrum_QCD10001_00_05",
        "ak7GenJetSpectrum_QCD10001_05_10",
        "ak7GenJetSpectrum_QCD10001_10_15",
        "ak7GenJetSpectrum_QCD10001_15_20",
        "ak7GenJetSpectrum_QCD10001_20_25"
    };

    
    int colorCode[N] = {
        1,2,4,6,1
    };
    int markerCode[N] = {
        20,29,34,33,27
    };
    int lineCode[N] = {
        1, 2, 3, 4,5
    };
    double scalepower[N]={4, 3, 2, 1,0};
    double scalef[N];
    
    for(int j=0; j<netabin; j++) scalef[j]=pow(10.0, scalepower[j]);

    char *R5dataname[N]={"R5Y0_0p5", "R5Y0p5_1p0", "R5Y1p0_1p5", "R5Y1p5_2p0", "R5Y2p0_2p5"
    };
    char *R7dataname[N]={"R7Y0_0p5", "R7Y0p5_1p0", "R7Y1p0_1p5", "R7Y1p5_2p0", "R7Y2p0_2p5"
    };

    char *R5powheg[N]={"d01-x01-y01", "d02-x01-y01", "d03-x01-y01", "d04-x01-y01", "d05-x01-y01"
    };
    char *R7powheg[N]={"d01-x01-y01", "d02-x01-y01", "d03-x01-y01", "d04-x01-y01", "d05-x01-y01"
    };

    cout <<"N =" << N <<endl ;
    
    TH1F * h7000_R5[N];
    TH1F * h7000_R7[N];
    TH1F * hratio_mc[N];
    TH1F * hratio_pow[N];
    TH1F * hratio_pow2[N];
    TH1F * hratio_data[N];
 //   TH1F * hratio_10011[N];
    TH1F * ConeDep[N];

    TFile * fR5DataFile;
    TFile * fR7DataFile;
    
    TH1F * hData2760[N];
    
    TGraphErrors * jetR05[N];
    TGraphErrors * jetR07[N];
    TGraphErrors * jet5TeV[N];
    
    TH1F * powhegR5[N];
    TH1F * powhegR7[N];
    TH1F * rebinpowhegR5[N];
    TH1F * rebinpowhegR7[N];

    TH1F * powheg2R5[N];
    TH1F * powheg2R7[N];
    TH1F * rebinpowheg2R5[N];
    TH1F * rebinpowheg2R7[N];

    TH1F * R5jet[N];
    TH1F * R7jet[N];

    TH1F * mcR5[N];
    TH1F * mcR7[N];

    
    TH1F *  hRatio[N];

    for(int ieta=0; ieta<N; ieta++)
    {
        R5jet[ieta] = new TH1F(Form("R5jetcsAbsEta%.f",deta[ieta+1]*10),Form("R5jetcsAbsEta%.f",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        R5jet[ieta]->Sumw2();
        R7jet[ieta] = new TH1F(Form("R7jetcsAbsEta%.f",deta[ieta+1]*10),Form("R7jetcsAbsEta%.f",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        R7jet[ieta]->Sumw2();

        mcR5[ieta] = new TH1F(Form("R5mcjetcsAbsEta%.f",deta[ieta+1]*10),Form("R5mcjetcsAbsEta%.f",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        mcR5[ieta]->Sumw2();
        mcR7[ieta] = new TH1F(Form("R7mcjetcsAbsEta%.f",deta[ieta+1]*10),Form("R7mcjetcsAbsEta%.f",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        mcR7[ieta]->Sumw2();

        powhegR5[ieta] = new TH1F(Form("powhegR5AbsEta%.f",deta[ieta+1]*10),Form("powhegR5AbsEta%.f",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        powhegR5[ieta]->Sumw2();

        powhegR7[ieta] = new TH1F(Form("powhegR7AbsEta%.f",deta[ieta+1]*10),Form("powhegR7AbsEta%.f",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        powhegR7[ieta]->Sumw2();

        powheg2R5[ieta] = new TH1F(Form("powheg5TeVR5AbsEta%.f",deta[ieta+1]*10),Form("powheg5TeVR5AbsEta%.f",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        powheg2R5[ieta]->Sumw2();
        
        powheg2R7[ieta] = new TH1F(Form("powheg5TeVR7AbsEta%.f",deta[ieta+1]*10),Form("powheg5TeVR7AbsEta%.f",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        powheg2R7[ieta]->Sumw2();

        h7000_R5[ieta] = (TH1F*) fMC->Get(Form("%s/JetSpectrum",R5dirname[ieta]));
        h7000_R7[ieta] = (TH1F*) fMC->Get(Form("%s/JetSpectrum",R7dirname[ieta]));
        
        rebinpowhegR5[ieta] = (TH1F*) fpowhegR5->Get(Form("%s",R5powheg[ieta]));
        rebinpowhegR5[ieta]->Scale(1./(nFile));
        rebinpowhegR7[ieta] = (TH1F*) fpowhegR7->Get(Form("%s",R7powheg[ieta]));
//        rebinpowhegR7[ieta]->Sumw2();
        rebinpowhegR7[ieta]->Scale(1./(nFile));

        rebinpowheg2R5[ieta] = (TH1F*) fpowheg2R5->Get(Form("%s",R5powheg[ieta]));
        rebinpowheg2R5[ieta]->Scale(1./(nFile));
        rebinpowheg2R7[ieta] = (TH1F*) fpowheg2R7->Get(Form("%s",R7powheg[ieta]));
        //        rebinpowhegR7[ieta]->Sumw2();
        rebinpowheg2R7[ieta]->Scale(1./(nFile));

   //     cout << "ieta =" << ieta <<endl ;
        
    }
    
    double R5ptCen[N][100], R7ptCen[N][100], xerrminus[N][100], xerrplus[N][100], R5yval[N][100],R5ystatup[N][100], R5ystatdown[N][100], R5ysysup[N][100], R5ysysdown[N][100],R7yval[N][100],R7ystatup[N][100], R7ystatdown[N][100], R7ysysup[N][100], R7ysysdown[N][100];
    ifstream infile;
    
    std::string inname ;
    std::string inname2 ;
    
    
    inname = Form("/Users/ymao/group/CMS/macros/ppRef/xTScaling/SMP13002/InclusiveJets_Table_AK5_new.txt");
    cout <<"input = " <<inname.c_str() <<endl ;
    //  FILE* infile = fopen(inname.c_str(), "r");
    infile.open(inname.c_str());
    std::string line;
    int ifile = 0;
    int np=0;
    
    //  infile.open(inname.c_str());
    // string line;
    Int_t row = 0;
    double v[1000][100];
    // cout << "going to reading files !!!" << endl ;
    //   cout << "infile = " << infile <<endl ;
    if(infile.fail()) {
        cerr << "unable to open file for reading" << endl;
        exit(1);
    }
    
    while(std::getline(infile, line) && row < 1000 )
    {
        int column = 0;
        //     cout << "inside loop now 2222222 !!!" << endl ;
        
        istringstream lineStream(line);
        while( lineStream >> v[row][column] && column < 100 ) column++;
        row++;
        //      cout << "row =" << row  << "column = " << column <<endl ;
    }
    infile.close();
    Int_t ieta=0;
    for(int irow = 0 ; irow < row ; irow++){
        if(v[irow][0]>2.0 || v[irow][3]> jetPtBin[nJetPtBin]) continue ;
        for(int ibin = 0 ; ibin < netabin ; ibin++){
            if(v[irow][0]==deta[ibin] && v[irow][1]==deta[ibin+1]) ieta=ibin ;
        }
        if(v[irow][2]==56) np=0 ;
        cout << " ieta = " << ieta << " np = " << np <<endl ;
        R5ptCen[ieta][np]=v[irow][2]+(v[irow][3]-v[irow][2])/2.;
        R5yval[ieta][np]=v[irow][4];
        R5ystatup[ieta][np]=(v[irow][5])*v[irow][4];
        R5ystatdown[ieta][np]=(v[irow][6])*v[irow][4];
        R5ysysup[ieta][np]=TMath::Sqrt(v[irow][7]*v[irow][7]+v[irow][10]*v[irow][10]+v[irow][12]*v[irow][12]+v[irow][14]*v[irow][14]+v[irow][16]*v[irow][16]+v[irow][18]*v[irow][18]+v[irow][20]*v[irow][20]+v[irow][22]*v[irow][22]+v[irow][24]*v[irow][24]+v[irow][26]*v[irow][26]+v[irow][28]*v[irow][28]+v[irow][30]*v[irow][30]+v[irow][32]*v[irow][32]+v[irow][34]*v[irow][34]+v[irow][36]*v[irow][36]+v[irow][38]*v[irow][38]+v[irow][40]*v[irow][40]+v[irow][42]*v[irow][42]+v[irow][44]*v[irow][44]);
        R5ysysdown[ieta][np]=TMath::Sqrt(v[irow][7]*v[irow][7]+v[irow][10]*v[irow][10]+v[irow][11]*v[irow][11]+v[irow][13]*v[irow][13]+v[irow][15]*v[irow][15]+v[irow][17]*v[irow][17]+v[irow][19]*v[irow][19]+v[irow][21]*v[irow][21]+v[irow][23]*v[irow][23]+v[irow][25]*v[irow][25]+v[irow][27]*v[irow][27]+v[irow][29]*v[irow][29]+v[irow][31]*v[irow][31]+v[irow][33]*v[irow][33]+v[irow][35]*v[irow][35]+v[irow][37]*v[irow][37]+v[irow][39]*v[irow][39]+v[irow][41]*v[irow][41]+v[irow][43]*v[irow][43]);
        
        
        int point = -1 ;
        for(int ibin = 0 ; ibin < nJetPtBin ; ibin++){
            if(v[irow][2]==jetPtBin[ibin] && v[irow][3]==jetPtBin[ibin+1]) point=ibin ;
        }
        if(point <0) continue ;
        int hbin = R5jet[ieta]->FindBin(R5ptCen[ieta][np]);
//        if(sys=="SysUp")
//            R5jet[ieta]->SetBinContent(hbin, R5yval[ieta][point]+R5ysysup[ieta][point]*R5yval[ieta][point]);
//        else if(sys=="SysDown")
//            R5jet[ieta]->SetBinContent(hbin, R5yval[ieta][point]-R5ysysdown[ieta][point]*R5yval[ieta][point]);
//        else
            R5jet[ieta]->SetBinContent(hbin, R5yval[ieta][point]);
        R5jet[ieta]->SetBinError(hbin, R5ystatup[ieta][point]);
        
        //      sys[np]=v2*TMath::Sqrt(2);
        np++;
        nptpoint[ieta] = np ;
        
    }

    
    inname = Form("/Users/ymao/group/CMS/macros/ppRef/xTScaling/SMP13002/InclusiveJets_Table_AK7_new.txt");
    cout <<"input = " <<inname.c_str() <<endl ;
    //  FILE* infile = fopen(inname.c_str(), "r");
    infile.open(inname.c_str());
    std::string line;
    int ifile = 0;
    int np=0;
    
    //  infile.open(inname.c_str());
    // string line;
    Int_t row = 0;
    double v[1000][100];
    // cout << "going to reading files !!!" << endl ;
    //   cout << "infile = " << infile <<endl ;
    if(infile.fail()) {
        cerr << "unable to open file for reading" << endl;
        exit(1);
    }
    
    while(std::getline(infile, line) && row < 1000 )
    {
        int column = 0;
        //     cout << "inside loop now 2222222 !!!" << endl ;
        
        istringstream lineStream(line);
        while( lineStream >> v[row][column] && column < 100 ) column++;
        row++;
        //      cout << "row =" << row  << "column = " << column <<endl ;
    }
    infile.close();
    Int_t ieta=0;
    for(int irow = 0 ; irow < row ; irow++){
        if(v[irow][0]>2.0 || v[irow][3]> jetPtBin[nJetPtBin]) continue ;
        for(int ibin = 0 ; ibin < netabin ; ibin++){
            if(v[irow][0]==deta[ibin] && v[irow][1]==deta[ibin+1]) ieta=ibin ;
        }
        if(v[irow][2]==56) np=0 ;
//        cout << " ieta = " << ieta << " np = " << np <<endl ;
        R7ptCen[ieta][np]=v[irow][2]+(v[irow][3]-v[irow][2])/2.;
        R7yval[ieta][np]=v[irow][4];
        R7ystatup[ieta][np]=(v[irow][5])*R7yval[ieta][np];
        R7ystatdown[ieta][np]=(v[irow][6])*R7yval[ieta][np];
        R7ysysup[ieta][np]=TMath::Sqrt(v[irow][7]*v[irow][7]+v[irow][10]*v[irow][10]+v[irow][12]*v[irow][12]+v[irow][14]*v[irow][14]+v[irow][16]*v[irow][16]+v[irow][18]*v[irow][18]+v[irow][20]*v[irow][20]+v[irow][22]*v[irow][22]+v[irow][24]*v[irow][24]+v[irow][26]*v[irow][26]+v[irow][28]*v[irow][28]+v[irow][30]*v[irow][30]+v[irow][32]*v[irow][32]+v[irow][34]*v[irow][34]+v[irow][36]*v[irow][36]+v[irow][38]*v[irow][38]+v[irow][40]*v[irow][40]+v[irow][42]*v[irow][42]+v[irow][44]*v[irow][44]);
        R7ysysdown[ieta][np]=TMath::Sqrt(v[irow][7]*v[irow][7]+v[irow][10]*v[irow][10]+v[irow][11]*v[irow][11]+v[irow][13]*v[irow][13]+v[irow][15]*v[irow][15]+v[irow][17]*v[irow][17]+v[irow][19]*v[irow][19]+v[irow][21]*v[irow][21]+v[irow][23]*v[irow][23]+v[irow][25]*v[irow][25]+v[irow][27]*v[irow][27]+v[irow][29]*v[irow][29]+v[irow][31]*v[irow][31]+v[irow][33]*v[irow][33]+v[irow][35]*v[irow][35]+v[irow][37]*v[irow][37]+v[irow][39]*v[irow][39]+v[irow][41]*v[irow][41]+v[irow][43]*v[irow][43]);
        
        int point = -1 ;
        for(int ibin = 0 ; ibin < nJetPtBin ; ibin++){
            if(v[irow][2]==jetPtBin[ibin] && v[irow][3]==jetPtBin[ibin+1]) point=ibin ;
        }
        if(point <0) continue ;
        int hbin = R7jet[ieta]->FindBin(R7ptCen[ieta][np]);
//        if(sys=="SysUp")
//            R7jet[ieta]->SetBinContent(hbin, R7yval[ieta][point]+R5ysysup[ieta][point]*R7yval[ieta][point]);
//        else if(sys=="SysDown")
//            R7jet[ieta]->SetBinContent(hbin, R7yval[ieta][point]-R5ysysdown[ieta][point]*R7yval[ieta][point]);
//        else
            R7jet[ieta]->SetBinContent(hbin, R7yval[ieta][point]);
        R7jet[ieta]->SetBinError(hbin, R5ystatup[ieta][point]);
        
        cout << " R7 data =" << R7ptCen[ieta][np] << " eta = "<< ieta << " value =" << R7yval[ieta][np] <<endl ;
        //      sys[np]=v2*TMath::Sqrt(2);
        np++;
        nptpoint[ieta] = np ;
        
    }

    if(SaveFile){
        ofstream myfile;
        std::string outName ;
        myfile << fixed << setprecision(5);
        for(int ieta = 0 ; ieta <netabin; ieta++){
            outName =Form("/Users/ymao/group/CMS/macros/ppRef/xTScaling/SMP13002/InclusiveJets_Table_AK7_Eta%d.txt",ieta);
            myfile.open(outName.c_str());
            for(int i=0;i<nptpoint[ieta];i++){
                myfile << R7ptCen[ieta][i]<<" \t "<< R7yval[ieta][i] <<" \t "<< R7ystatup[ieta][i]/R7yval[ieta][i] <<" \t "<< R7ystatdown[ieta][i]/R7yval[ieta][i] <<" \t "<< R7ysysup[ieta][i] <<" \t "<< R7ysysdown[ieta][i] <<" \t "<<"\n";
            }
            myfile.close();
        }
    }
    
    for(int ieta = 0 ; ieta <netabin; ieta++){
//        R5jet[ieta]->Scale(1./(deta[ieta+1]-deta[ieta]));
//        R7jet[ieta]->Scale(1./(deta[ieta+1]-deta[ieta]));

        for(int ibin = 0 ; ibin < mcR5[ieta]->GetNbinsX(); ibin++){
            double binCenter = mcR5[ieta]->GetBinCenter(ibin);
            int hbin = h7000_R5[ieta]->FindBin(binCenter);
            mcR5[ieta]->SetBinContent(ibin, h7000_R5[ieta]->GetBinContent(hbin));
            mcR5[ieta]->SetBinError(ibin, h7000_R5[ieta]->GetBinError(hbin));
        }
        
        for(int ibin = 0 ; ibin < mcR7[ieta]->GetNbinsX(); ibin++){
            double binCenter = mcR7[ieta]->GetBinCenter(ibin);
            int hbin = h7000_R7[ieta]->FindBin(binCenter);
            mcR7[ieta]->SetBinContent(ibin, h7000_R7[ieta]->GetBinContent(hbin));
            mcR7[ieta]->SetBinError(ibin, h7000_R7[ieta]->GetBinError(hbin));
        }
        for(int ibin = 0 ; ibin < powhegR5[ieta]->GetNbinsX(); ibin++){
            double binCenter = powhegR5[ieta]->GetBinCenter(ibin);
            int hbin = rebinpowhegR5[ieta]->FindBin(binCenter);
            powhegR5[ieta]->SetBinContent(ibin, rebinpowhegR5[ieta]->GetBinContent(hbin));
            powhegR5[ieta]->SetBinError(ibin, rebinpowhegR5[ieta]->GetBinError(hbin));
        }
        for(int ibin = 0 ; ibin < powhegR7[ieta]->GetNbinsX(); ibin++){
            double binCenter = powhegR7[ieta]->GetBinCenter(ibin);
            int hbin = rebinpowhegR7[ieta]->FindBin(binCenter);
            powhegR7[ieta]->SetBinContent(ibin, rebinpowhegR7[ieta]->GetBinContent(hbin));
            powhegR7[ieta]->SetBinError(ibin, rebinpowhegR7[ieta]->GetBinError(hbin));
            
        }
        for(int ibin = 0 ; ibin < powheg2R5[ieta]->GetNbinsX(); ibin++){
            double binCenter = powheg2R5[ieta]->GetBinCenter(ibin);
            int hbin = rebinpowheg2R5[ieta]->FindBin(binCenter);
            powheg2R5[ieta]->SetBinContent(ibin, rebinpowheg2R5[ieta]->GetBinContent(hbin));
            powheg2R5[ieta]->SetBinError(ibin, rebinpowheg2R5[ieta]->GetBinError(hbin));
        }
        for(int ibin = 0 ; ibin < powheg2R7[ieta]->GetNbinsX(); ibin++){
            double binCenter = powheg2R7[ieta]->GetBinCenter(ibin);
            int hbin = rebinpowheg2R7[ieta]->FindBin(binCenter);
            powheg2R7[ieta]->SetBinContent(ibin, rebinpowheg2R7[ieta]->GetBinContent(hbin));
            powheg2R7[ieta]->SetBinError(ibin, rebinpowheg2R7[ieta]->GetBinError(hbin));
            
        }
//        if(ieta>0){
//            powhegR5[ieta]->Scale(1./(deta[ieta+1]-deta[ieta]));
//            powhegR7[ieta]->Scale((1./deta[ieta+1]-deta[ieta]));
//        }
//        jetR05[ieta] = new TGraphErrors(nptpoint[ieta], R5ptCen[ieta], R5yval[ieta], 0., R5ystatup[ieta]);
//        jetR07[ieta] = new TGraphErrors(nptpoint[ieta], R7ptCen[ieta], R7yval[ieta], 0., R7ystatup[ieta]);
//
        hratio_data[ieta] = (TH1F*)R5jet[ieta]->Clone(Form("DataCSRatioR7_R5_EtaBin%.f_%.f", deta[ieta]*10, deta[ieta+1]*10));
        hratio_data[ieta]->Divide(R7jet[ieta]);

        hratio_mc[ieta] = (TH1F*) mcR5[ieta]->Clone(Form("hPYTHIAratio_R7_R5_%s",R7dataname[ieta]));
        hratio_mc[ieta]->Divide(mcR7[ieta]);
        
        hratio_pow[ieta] = (TH1F*) powhegR5[ieta]->Clone(Form("hPOWHEGratio_R7_R5_EtaBin%d",ieta));
        hratio_pow[ieta]->Divide(powhegR7[ieta]);

        hratio_pow2[ieta] = (TH1F*) powheg2R5[ieta]->Clone(Form("hPOWHEG5TeVratio_R7_R5_EtaBin%d",ieta));
        hratio_pow2[ieta]->Divide(powheg2R7[ieta]);

        for(int ibin = 0 ; ibin < hratio_data[ieta]->GetNbinsX(); ibin++){
            double binCenter = hratio_data[ieta]->GetBinCenter(ibin);
            if(binCenter>ptboundary[ieta]){
                hratio_data[ieta]->SetBinContent(ibin, -1.e9);
                hratio_data[ieta]->SetBinError(ibin,-1.e-9 );

            }
        }
        for(int ibin = 0 ; ibin < hratio_mc[ieta]->GetNbinsX(); ibin++){
            double binCenter = hratio_mc[ieta]->GetBinCenter(ibin);
            if(binCenter>ptboundary[ieta]){
                hratio_mc[ieta]->SetBinContent(ibin, -1.e9);
                hratio_mc[ieta]->SetBinError(ibin,-1.e-9 );
                
            }
        }

        for(int ibin = 0 ; ibin < hratio_pow[ieta]->GetNbinsX(); ibin++){
            double binCenter = hratio_pow[ieta]->GetBinCenter(ibin);
            if(binCenter>ptboundary[ieta]){
                hratio_pow[ieta]->SetBinContent(ibin, -1.e9);
                hratio_pow[ieta]->SetBinError(ibin,-1.e-9 );
                
            }
        }
        
        hRatio[ieta] = (TH1F*)hratio_data[ieta]->Clone(Form("DataMCDoubleRatioR7_R5_EtaBin%.f_%.f", deta[ieta]*10, deta[ieta+1]*10));
 //       hRatio[ieta]->Divide(hratio_mc[ieta]);
        hRatio[ieta]->Divide(hratio_pow[ieta]);
//        hRatio[ieta] = (TH1F*)powhegR5[ieta]->Clone(Form("DataMCDoubleRatioR7_R5_EtaBin%.f_%.f", deta[ieta]*10, deta[ieta+1]*10));
//        hRatio[ieta]->Divide(R5jet[ieta]);
        
        for(int ibin = 0 ; ibin < hRatio[ieta]->GetNbinsX(); ibin++){
            double binCenter = hRatio[ieta]->GetBinCenter(ibin);
            if(binCenter>ptboundary[ieta]){
                hRatio[ieta]->SetBinContent(ibin, -1.e9);
                hRatio[ieta]->SetBinError(ibin,-1.e-9 );
                
            }
        }
    }

    
    TCanvas *c1 = new TCanvas("c1a", "c1",0,0,959,389);
    makeMultiPanelCanvas(c1,4,1,0.0,0.0,0.2,0.2,0.02);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    gStyle->SetOptStat(0);
//    gStyle->SetPadBottomMargin(0.12);
//    gStyle->SetPadTopMargin   (0.025);
//    gStyle->SetPadLeftMargin  (0.15);
//    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
    
    TCanvas *c2 = new TCanvas("c2a", "c2",0,0,600,600);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    c2->Range(0,0,1,1);
    c2->SetFillColor(0);
    c2->SetBorderMode(0);
    c2->SetBorderSize(2);
    //c1->SetLogy();
    c2->SetTickx(1);
    c2->SetTicky(1);
    c2->SetLeftMargin(0.13);
    c2->SetRightMargin(0.06);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.16);
    c2->SetFrameFillStyle(0);
    c2->SetFrameBorderMode(0);
    
    //   TLegend *t1=new TLegend(0.45,0.6,0.8,0.92);
    TLegend *t1=new TLegend(0.15,0.78,0.35,0.92);
    t1->SetFillColor(0);
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextFont(43);
    t1->SetTextSize(20);
    
    //   TLegend *t2=new TLegend(0.55,0.65,0.8,0.9);
    TLegend *t2=new TLegend(0.15,0.18,0.35,0.35);
    t2->SetFillColor(0);
    t2->SetBorderSize(0);
    t2->SetFillStyle(0);
    t2->SetTextFont(63);
    t2->SetTextSize(24);
    
    TLegend *t3=new TLegend(0.55,0.7,0.7,0.92);
    //      TLegend *t3=new TLegend(0.25,0.7,0.8,0.85);
    t3->SetFillColor(0);
    t3->SetBorderSize(0);
    t3->SetFillStyle(0);
    t3->SetTextFont(63);
    t3->SetTextSize(24);
    
    TLegend *t4=new TLegend(0.48,0.62,0.8,0.85);
    t4->SetFillColor(0);
    t4->SetBorderSize(0);
    t4->SetFillStyle(0);
    t4->SetTextFont(63);
    t4->SetTextSize(24);
    
    TH1F * hDum = new TH1F("hhdum","",2000,0., 2000.);
        hDum->SetAxisRange(55., 650., "X");
    hDum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    //    hDum->GetXaxis()->SetRangeUser(20., 600.);
    hDum->GetXaxis()->SetTitleSize(0.035);
    hDum->GetXaxis()->SetTitleOffset(1.5);
    hDum->GetXaxis()->SetLabelSize(0.035);
    hDum->GetYaxis()->SetTitleSize(0.035);
    hDum->GetYaxis()->SetLabelSize(0.035);
    hDum->GetYaxis()->SetTitleOffset(1.5);
    hDum->GetXaxis()->CenterTitle();
    hDum->GetYaxis()->CenterTitle();
    hDum->GetYaxis()->SetNdivisions(510);
    hDum->GetXaxis()->SetNdivisions(510);

    c1->cd();
    c1->SetLogy();
//   c1->SetLogx();
    hDum->SetMinimum(0.6);
    hDum->SetMaximum(1.1);
//    hDum->SetMinimum(0.);
//    hDum->SetMaximum(2.);
    hhdum->GetYaxis()->SetTitle("R=0.5/R=0.7");

    //   drawCMSpp(0.18,0.9,ppLumi);

    for(int ieta = 0 ; ieta <netabin-1; ieta++){
        c1->cd(ieta+1);
//        c1->cd(ieta+1)->SetLogy();
//        hDum->SetMinimum(1.e-3);
//        hDum->SetMaximum(1.e8);
//        hhdum->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}d#eta [pb/GeV]");
        fixedFontHist(hDum, 1.5, 1.6);
        hDum->DrawCopy();
        if(ieta==0){
        drawText(Form("anti-k_{T} PF jet"),0.4,0.80,16);
        //    drawText(Form("PYTHIA Z2"),0.2,0.9,18);
        drawText(Form(" CMS, pp #sqrt{s} = 7 TeV"),0.4,0.88,16);
        }
        //    int ieta = current ;
        hratio_mc[ieta]->SetMarkerStyle(29);
        hratio_mc[ieta]->SetMarkerColor(6);
        hratio_mc[ieta]->SetMarkerSize(1.5);
        hratio_mc[ieta]->SetLineColor(6);
        hratio_mc[ieta]->SetLineStyle(1);
        hratio_mc[ieta]->Draw("same L");

        hratio_data[ieta]->SetMarkerStyle(20);
        hratio_data[ieta]->SetMarkerColor(1);
        hratio_data[ieta]->SetMarkerSize(1.5);
        hratio_data[ieta]->SetLineColor(1);
        hratio_data[ieta]->Draw("same PE");


        hratio_pow[ieta]->SetMarkerStyle(27);
        hratio_pow[ieta]->SetMarkerColor(2);
        hratio_pow[ieta]->SetMarkerSize(1.5);
        hratio_pow[ieta]->SetLineColor(2);
        hratio_pow[ieta]->SetLineStyle(2);
        hratio_pow[ieta]->Draw("same L");

        hratio_pow2[ieta]->SetMarkerStyle(28);
        hratio_pow2[ieta]->SetMarkerColor(4);
        hratio_pow2[ieta]->SetMarkerSize(1.5);
        hratio_pow2[ieta]->SetLineColor(4);
        hratio_pow2[ieta]->SetLineStyle(4);
        hratio_pow2[ieta]->Draw("same L");

//        R5jet[ieta]->SetMarkerStyle(20);
//        R5jet[ieta]->SetMarkerColor(1);
//        R5jet[ieta]->SetMarkerSize(1.5);
//        R5jet[ieta]->SetLineColor(1);
//        R5jet[ieta]->Draw("same PE");
//        
//        powhegR5[ieta]->SetMarkerStyle(27);
//        powhegR5[ieta]->SetMarkerColor(2);
//        powhegR5[ieta]->SetMarkerSize(1.5);
//        powhegR5[ieta]->SetLineColor(2);
//        powhegR5[ieta]->SetLineStyle(2);
//        powhegR5[ieta]->Draw("same L");

//        R7jet[ieta]->SetMarkerStyle(20);
//        R7jet[ieta]->SetMarkerColor(1);
//        R7jet[ieta]->SetMarkerSize(1.5);
//        R7jet[ieta]->SetLineColor(1);
//        R7jet[ieta]->Draw("same PE");
//        
//        powhegR7[ieta]->SetMarkerStyle(27);
//        powhegR7[ieta]->SetMarkerColor(2);
//        powhegR7[ieta]->SetMarkerSize(1.5);
//        powhegR7[ieta]->SetLineColor(2);
//        powhegR7[ieta]->SetLineStyle(2);
//        powhegR7[ieta]->Draw("same L");

        regSun(55.,1.,650.,1.,1, 1);

            if(ieta==0)
                drawText(Form("|y| < %.1f ", deta[ieta+1]),0.3, 0.25, 17);
              else
                  drawText(Form("%.1f < |y| < %.1f ", deta[ieta],deta[ieta+1]),0.2, 0.25, 17);
        if(ieta==2){
            t1->AddEntry(hratio_data[ieta], "DATA: SMP-13-002","P");
            t1->AddEntry(hratio_mc[ieta], "PYTHIA","P");
            t1->AddEntry(hratio_pow[ieta], "POWHEG:7TeV","P");
            t1->AddEntry(hratio_pow2[ieta], "POWHEG:5.02TeV","P");
//            t1->AddEntry(R5jet[ieta], "DATA: SMP-13-002","P");
//            t1->AddEntry(powhegR5[ieta], "POWHEG","P");
            t1->Draw("same");
  
        }
        
    }
    if(SavePlot)c1->Print(Form("%s/pp7TeVPFJetDataPOWHEGDiffConeSpectraRatio.gif", plotsdir.Data()));
    c1->Update();
    
    c2->cd();
    hDum->SetMinimum(0.75);
    hDum->SetMaximum(1.25);
//    hDum->SetMinimum(0.5);
//    hDum->SetMaximum(1.5);
    hhdum->GetYaxis()->SetTitle("Data/POWHEG");
    fixedFontHist(hDum, 1.3, 1.5);
    hDum->DrawCopy();
    drawText(Form("anti-k_{T} PF jet"),0.2,0.75,17);
    //    drawText(Form("PYTHIA Z2"),0.2,0.9,18);
    drawText(Form("CMS, pp #sqrt{s} = 7 TeV"),0.2,0.82,17);
    
     for(int ieta = 0 ; ieta <netabin-1; ieta++){
        //    int ieta = current ;
        hRatio[ieta]->SetMarkerStyle(markerCode[ieta]);
        hRatio[ieta]->SetMarkerColor(colorCode[ieta]);
        hRatio[ieta]->SetMarkerSize(1.5);
        hRatio[ieta]->SetLineColor(colorCode[ieta]);
        hRatio[ieta]->Draw("same PE");
         if(ieta==0)
             t3->AddEntry(hRatio[ieta], Form("|y| < %.1f ", deta[ieta+1]),"P");
         
         else
             t3->AddEntry(hRatio[ieta], Form("%.1f < |y| < %.1f ", deta[ieta],deta[ieta+1], scalepower[ieta]),"P");

        t3->Draw("same");
    }
    
    regSun(55.,1.,650.,1.,1, 1);
    
    if(SavePlot)c2->Print(Form("%s/pp7TeVakPFJetDataPOWHEGRatio.gif", plotsdir.Data()));
    c2->Update();
 

}

//---------------------------------------------------
void drawText(const char *text, float xp, float yp, int textSize=15){
    TLatex *tex = new TLatex(xp,yp,text);
    tex->SetTextFont(63);
    //tex->SetTextSize(20);
    tex->SetTextSize(textSize);
    //tex->SetTextSize(0.05);
    tex->SetTextColor(kBlack);
    tex->SetLineWidth(1);
    tex->SetNDC();
    tex->Draw();
}
void drawText2(const char *text, float xp, float yp, int textSize=15){
    TLatex *tex = new TLatex(xp,yp,text);
    tex->SetTextFont(63);
    tex->SetTextSize(textSize);
    //    tex->SetTextColor(kBlack);
    tex->SetTextColor(kRed);
    tex->SetLineWidth(1);
    tex->SetNDC();
    tex->Draw();
}

void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
    if (canv==0) {
        Error("makeMultiPanelCanvas","Got null canvas.");
        return;
    }
    canv->Clear();
    
    TPad* pad[columns][rows];
    
    Float_t Xlow[columns];
    Float_t Xup[columns];
    Float_t Ylow[rows];
    Float_t Yup[rows];
    Float_t PadWidth =
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
                      (1.0/(1.0-edge))+(Float_t)columns-2.0);
    Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
                        (1.0/(1.0-edge))+(Float_t)rows-2.0);
    Xlow[0] = leftOffset;
    Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
    Xup[columns-1] = 1;
    Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
    
    Yup[0] = 1;
    Ylow[0] = 1.0-PadHeight/(1.0-edge);
    Ylow[rows-1] = bottomOffset;
    Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
    
    for(Int_t i=1;i<columns-1;i++) {
        Xlow[i] = Xup[0] + (i-1)*PadWidth;
        Xup[i] = Xup[0] + (i)*PadWidth;
    }
    Int_t ct = 0;
    for(Int_t i=rows-2;i>0;i--) {
        Ylow[i] = Yup[rows-1] + ct*PadHeight;
        Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
        ct++;
    }
    
    TString padName;
    for(Int_t i=0;i<columns;i++) {
        for(Int_t j=0;j<rows;j++) {
            canv->cd();
            padName = Form("p_%d_%d",i,j);
            pad[i][j] = new TPad(padName.Data(),padName.Data(),
                                 Xlow[i],Ylow[j],Xup[i],Yup[j]);
            if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
            else pad[i][j]->SetLeftMargin(0);
            
            if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
            else pad[i][j]->SetRightMargin(0);
            
            if(j==0) pad[i][j]->SetTopMargin(edge);
            else pad[i][j]->SetTopMargin(0);
            
            if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
            else pad[i][j]->SetBottomMargin(0);
            
            pad[i][j]->Draw();
            pad[i][j]->cd();
            pad[i][j]->SetNumber(columns*j+i+1);
        }
    }
}
void drawCMSmc(float px, float py, TString coll) {
    TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    cms->SetTextFont(63);
    cms->SetTextSize(17);
    cms->SetNDC();
    cms->Draw();
    if(coll=="HI")TLatex *lumi = new TLatex(px,py-0.05,"PYTHIA+HYDJET");
    else TLatex *lumi = new TLatex(px,py-0.05,"PYTHIA");
    lumi->SetTextFont(63);
    lumi->SetTextSize(15);
    lumi->SetNDC();
    lumi->Draw();
}

void drawCMS(float px, float py, float nLumi) {
    TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    cms->SetTextFont(63);
    cms->SetTextSize(17);
    cms->SetNDC();
    cms->Draw();
    TLatex *lumi = new TLatex(px,py-0.10,Form("#intL dt = %.1f #mub^{-1}",nLumi));
    lumi->SetTextFont(63);
    lumi->SetTextSize(15);
    lumi->SetNDC();
    lumi->Draw();
}

void drawCMSpp(float px, float py, float nLumi) {
    TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    cms->SetTextFont(63);
    cms->SetTextSize(15);
    cms->SetNDC();
    cms->Draw();
    TLatex *lumi = new TLatex(px,py-0.05,Form("#intL dt = %.1f pb^{-1}",nLumi));
    //    TLatex *lumi = new TLatex(px,py-0.10,Form("#intL dt = %.1f nb^{-1}",nLumi));
    lumi->SetTextFont(63);
    lumi->SetTextSize(15);
    lumi->SetNDC();
    lumi->Draw();
}
void fixedFontHist(TH1 * h, Float_t xoffset=1.3, Float_t yoffset=1.2)
{
    h->SetLabelFont(43,"X");
    h->SetLabelFont(43,"Y");
    //h->SetLabelOffset(0.01);
    h->SetLabelSize(22);
    h->SetTitleFont(44);
    h->SetTitleSize(16);
    h->SetLabelSize(18,"Y");
    h->SetLabelSize(18,"X");
    h->SetTitleFont(43,"Y");
    h->SetTitleSize(22,"Y");
    h->SetTitleSize(20,"X");
    h->SetTitleOffset(xoffset,"X");
    h->SetTitleOffset(yoffset,"Y");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
}

void drawSys(TGraph *h, double *sys, double width=5, int theColor= kYellow, int fillStyle = -1, int lineStyle = -1)
{
    for (int i=0;i<h->GetN();i++)
    {
        double val;
        double theX;
        h->GetPoint(i,theX,val);
        double err = val * sys[i]/100;
        TBox *b = new TBox(theX-width,val-err,theX+width,val+err);
        
        b->SetLineColor(theColor);
        b->SetFillColor(theColor);
        if ( fillStyle > -1 ) b->SetFillStyle(fillStyle);
        if ( lineStyle > -1 ) b->SetLineStyle(lineStyle);
        
        b->Draw();
    }
}
void onSun(double x1=0,double y1=0,double x2=1,double y2=1,int color=1, double width=1)
{
    TLine* t1 = new TLine(x1,y1,x2,y2);
    t1->SetLineWidth(width);
    t1->SetLineStyle(1);
    t1->SetLineColor(color);
    t1->Draw();
}
void regSun(double x1=0,double y1=0,double x2=1,double y2=1,int color=1, double width=1)
{
    TLine* t1 = new TLine(x1,y1,x2,y2);
    t1->SetLineWidth(width);
    t1->SetLineStyle(3);
    t1->SetLineColor(color);
    t1->Draw();
}
TH1 * convertToInvYield(TH1 *hist) {
    for(int i = 1; i<=hist->GetNbinsX(); i++) {
        double content = hist->GetBinContent(i);
        double pt = hist->GetBinCenter(i);
        double error = hist->GetBinError(i);
        
        double new_content = content/(2.*TMath::Pi()*pt);
        double new_error = error/(2.*TMath::Pi()*pt);
        
        hist->SetBinContent(i,new_content);
        hist->SetBinError(i,new_error);
    }
    return hist ;
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
