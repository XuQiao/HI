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
const TString JetIDName[]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","eSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN","PPcut","PPcutTight","PPcutTighter","SumSumpt1"};
const int nJetID = sizeof(JetIDName)/sizeof(TString);
const double binbound_pt_coarse[]={30,50,80,100,600};
const int Nbin_pt_coarse=sizeof(binbound_pt_coarse)/sizeof(double)-1;
TString algo="akPu3PF";//"akPu3PF"
TString coll = "PPb";

class hist_class{
public:
    hist_class();
    void Write();
    TH1F * refjetpt;
    TH1F * genjetpt;
    TH1F * jetpt;
    TH1F * jetpt_jetidclosure;
    TH1F * jetpt_fake;
    TH1F * jetpt_real;
    TH1F * jetpt_fake1;
    TH1F * jetpt_real1;
    TH2F * jetptEta;
    TH2F * jetptphi;
    TH2F * jetEtaphi;
    TH2F * jetptjetid[nJetID];
    TH2F * jetptjetid_real[nJetID];
    TH2F * jetptjetid_fake[nJetID];
    TH2F * jetptjetid_real1[nJetID];
    TH2F * jetptjetid_fake1[nJetID];
    TH2F * jetid12jetid14[Nbin_pt_coarse];
    TH2F * jetid12jetid14_fake1[Nbin_pt_coarse];
    TH2F * jetid12jetid14_real1[Nbin_pt_coarse];
    TH2F * rawptJES;
    TH2F * refptJES;
    TH1F * Vz;
    TH1F * VzW;
    TH1F * Cent;
    TH1F * CentW;
    TH1F * Pthat;
    TH1F * PthatAfXw;
    TH2F* Spikeetaphi;
    TH2F* Spikeetapt;
    TH2F* Spikeetaphi_real;
    TH2F* Spikeetapt_real;
    TH2F* Spikeetaphi_fake;
    TH2F* Spikeetapt_fake;

    TH2F * rawptJESEtaBin[netabin];
    TH1F * refjetptEtaBin[netabin];
    TH1F * genjetptEtaBin[netabin];
    TH2F * jetptjetidEtaBin[netabin][nJetID];
    TH2F * jetptjetidEtaBin_fake[netabin][nJetID];
    TH2F * jetptjetidEtaBin_real[netabin][nJetID];
    TH2F * jetptjetidEtaBin_fake1[netabin][nJetID];
    TH2F * jetptjetidEtaBin_real1[netabin][nJetID];
    TH2F * refptJESEtaBin[netabin];
    TH1F * jetptEtaBin[netabin];
    TH1F * jetptEtaBin_jetidclosure[netabin];
    TH1F * jetptEtaBin_fake[netabin];
    TH1F * jetptEtaBin_real[netabin];
    TH1F * jetptEtaBin_fake1[netabin];
    TH1F * jetptEtaBin_real1[netabin];
};

hist_class::hist_class()
{
    refjetpt = new TH1F(Form("refjetpt"),Form("refjetpt"),1000,0.,1000.);
    refjetpt -> Sumw2(); 
    genjetpt = new TH1F(Form("genjetpt"),Form("genjetpt"),1000,0.,1000.);
    genjetpt -> Sumw2(); 
    jetpt = new TH1F(Form("jetpt"),Form("jetpt"),1000,0.,1000.);
    jetpt -> Sumw2(); 
    jetpt_jetidclosure = new TH1F(Form("jetpt_jetidclosure"),Form("jetpt_jetidclosure"),1000,0.,1000.);
    jetpt_jetidclosure -> Sumw2(); 
    jetpt_fake = new TH1F(Form("jetpt_fake"),Form("jetpt_fake"),1000,0.,1000.);
    jetpt_fake -> Sumw2(); 
    jetpt_real = new TH1F(Form("jetpt_real"),Form("jetpt_real"),1000,0.,1000.);
    jetpt_real -> Sumw2(); 
    jetpt_fake1 = new TH1F(Form("jetpt_fake1"),Form("jetpt_fake1"),1000,0.,1000.);
    jetpt_fake1 -> Sumw2(); 
    jetpt_real1 = new TH1F(Form("jetpt_real1"),Form("jetpt_real1"),1000,0.,1000.);
    jetpt_real1 -> Sumw2(); 
    jetptEta = new TH2F(Form("jetptEta"),Form("jetptEta"),1000,0.,1000.,100,-5.,5.);
    jetptEta -> Sumw2(); 
    jetptphi = new TH2F(Form("jetptphi"),Form("jetptphi"),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi -> Sumw2(); 
    jetEtaphi = new TH2F(Form("jetEtaphi"),Form("jetEtaphi"),100,-5.,5.,200, -TMath::Pi(), TMath::Pi());
    jetEtaphi -> Sumw2(); 
    for(int ipt=0;ipt<Nbin_pt_coarse;ipt++){ 
    jetid12jetid14[ipt] = new TH2F(Form("jetid12jetid14_%d",ipt),Form("jetid12jetid14_%d",ipt),200,0.,2.,200,0.,2.);//Added
    jetid12jetid14[ipt] -> Sumw2();
    jetid12jetid14_fake1[ipt] = new TH2F(Form("jetid12jetid14_%d_fake1",ipt),Form("jetid12jetid14_%d_fake1",ipt),200,0.,2.,200,0.,2.);//Added
    jetid12jetid14_fake1[ipt] -> Sumw2();
    jetid12jetid14_real1[ipt] = new TH2F(Form("jetid12jetid14_%d_real1",ipt),Form("jetid12jetid14_%d_real1",ipt),200,0.,2.,200,0.,2.);//Added
    jetid12jetid14_real1[ipt] -> Sumw2();
    }
        Spikeetaphi = new TH2F(Form("Spikeetaphi"),Form("Spikeetaphi"),200,-5.,5.,200, -TMath::Pi(), TMath::Pi());
        Spikeetaphi -> Sumw2();
        Spikeetapt = new TH2F(Form("Spikeetapt"),Form("Spikeetapt"),200,-5.,5.,1000,0.,1000.);
        Spikeetapt -> Sumw2();
        Spikeetaphi_real = new TH2F(Form("Spikeetaphi_real1"),Form("Spikeetaphi_real1"),200,-5.,5.,200, -TMath::Pi(), TMath::Pi());
        Spikeetaphi_real -> Sumw2();
        Spikeetapt_real = new TH2F(Form("Spikeetapt_real1"),Form("Spikeetapt_real1"),200,-5.,5.,1000,0.,1000.);
        Spikeetapt_real -> Sumw2();
        Spikeetaphi_fake = new TH2F(Form("Spikeetaphi_fake1"),Form("Spikeetaphi_fake1"),200,-5.,5.,200, -TMath::Pi(), TMath::Pi());
        Spikeetaphi_fake -> Sumw2();
        Spikeetapt_fake = new TH2F(Form("Spikeetapt_fake1"),Form("Spikeetapt_fake1"),200,-5.,5.,1000,0.,1000.);
        Spikeetapt_fake -> Sumw2();
	
        for(int ijetid=0;ijetid<nJetID;ijetid++){
        if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr")){
        jetptjetid[ijetid] = new TH2F(Form("jetpt%s",JetIDName[ijetid].Data()), Form("jetpt%s",JetIDName[ijetid].Data()), 1000, 0., 1000., 200, 0., 2.);      //Added
        jetptjetid_real[ijetid] = new TH2F(Form("jetpt%s_real",JetIDName[ijetid].Data()), Form("jetpt%s_real",JetIDName[ijetid].Data()), 1000, 0., 1000., 200, 0., 2.);       //Added
        jetptjetid_fake[ijetid] = new TH2F(Form("jetpt%s_fake",JetIDName[ijetid].Data()), Form("jetpt%s_fake",JetIDName[ijetid].Data()), 1000, 0., 1000., 200, 0., 2.);       //Added
        jetptjetid_real1[ijetid] = new TH2F(Form("jetpt%s_real1",JetIDName[ijetid].Data()), Form("jetpt%s_real1",JetIDName[ijetid].Data()), 1000, 0., 1000., 200, 0., 2.);       //Added
        jetptjetid_fake1[ijetid] = new TH2F(Form("jetpt%s_fake1",JetIDName[ijetid].Data()), Form("jetpt%s_fake1",JetIDName[ijetid].Data()), 1000, 0., 1000., 200, 0., 2.);       //Added
        }
        else if(JetIDName[ijetid].Contains("N") || JetIDName[ijetid].Contains("PP")){
        jetptjetid[ijetid] = new TH2F(Form("jetpt%s",JetIDName[ijetid].Data()), Form("jetpt%s",JetIDName[ijetid].Data()), 1000, 0., 1000., 100, 0., 100.);    //Added
        jetptjetid_real[ijetid] = new TH2F(Form("jetpt%s_real",JetIDName[ijetid].Data()), Form("jetpt%s_real",JetIDName[ijetid].Data()), 1000, 0., 1000., 100, 0., 100.);     //Added
        jetptjetid_fake[ijetid] = new TH2F(Form("jetpt%s_fake",JetIDName[ijetid].Data()), Form("jetpt%s_fake",JetIDName[ijetid].Data()), 1000, 0., 1000., 100, 0., 100.);     //Added
        jetptjetid_real1[ijetid] = new TH2F(Form("jetpt%s_real1",JetIDName[ijetid].Data()), Form("jetpt%s_real1",JetIDName[ijetid].Data()), 1000, 0., 1000., 100, 0., 100.);     //Added
        jetptjetid_fake1[ijetid] = new TH2F(Form("jetpt%s_fake1",JetIDName[ijetid].Data()), Form("jetpt%s_fake1",JetIDName[ijetid].Data()), 1000, 0., 1000., 100, 0., 100.);     //Added
        }
        else{
        jetptjetid[ijetid] = new TH2F(Form("jetpt%s",JetIDName[ijetid].Data()), Form("jetpt%s",JetIDName[ijetid].Data()), 1000, 0., 1000., 3000, 0., 300.);   //Added
        jetptjetid_real[ijetid] = new TH2F(Form("jetpt%s_real",JetIDName[ijetid].Data()), Form("jetpt%s_real",JetIDName[ijetid].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptjetid_fake[ijetid] = new TH2F(Form("jetpt%s_fake",JetIDName[ijetid].Data()), Form("jetpt%s_fake",JetIDName[ijetid].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptjetid_real1[ijetid] = new TH2F(Form("jetpt%s_real1",JetIDName[ijetid].Data()), Form("jetpt%s_real1",JetIDName[ijetid].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptjetid_fake1[ijetid] = new TH2F(Form("jetpt%s_fake1",JetIDName[ijetid].Data()), Form("jetpt%s_fake1",JetIDName[ijetid].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        }
        jetptjetid[ijetid]->Sumw2();      //Added
        jetptjetid_real[ijetid]->Sumw2(); //Added
        jetptjetid_fake[ijetid]->Sumw2(); //Added
        jetptjetid_real1[ijetid]->Sumw2(); //Added
        jetptjetid_fake1[ijetid]->Sumw2(); //Added
        }

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
    CentW = new TH1F(Form("CentW"),Form("CentW"),100,0.,100.);
    CentW-> Sumw2(); 
    Pthat = new TH1F(Form("Pthat"),Form("Pthat"),1000,0.,1000.);
    Pthat-> Sumw2(); 
    PthatAfXw = new TH1F(Form("PthatAfXw"),Form("PthatAfXw"),1000,0.,1000.);
    PthatAfXw-> Sumw2(); 

    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[ieta]->Sumw2();
      jetptEtaBin_jetidclosure[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_jetidclosure",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f_jetidclosure",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_jetidclosure[ieta]->Sumw2();
      jetptEtaBin_fake[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_fake",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin_fake%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_fake[ieta]->Sumw2();
      jetptEtaBin_real[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_real",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin_real%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_real[ieta]->Sumw2();
      jetptEtaBin_fake1[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_fake1",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin_fake1%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_fake1[ieta]->Sumw2();
      jetptEtaBin_real1[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_real1",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin_real1%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin_real1[ieta]->Sumw2();
       for(int ijetid=0;ijetid<nJetID;ijetid++){
        if(JetIDName[ijetid].Contains("pt") || JetIDName[ijetid].Contains("Maxr")){
        jetptjetidEtaBin[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptjetidEtaBin_fake[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_fake",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_fake",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptjetidEtaBin_real[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_real",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_real",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptjetidEtaBin_fake1[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_fake1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_fake1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptjetidEtaBin_real1[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_real1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_real1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 2.);    //Added
        }
        else if(JetIDName[ijetid].Contains("PP") || JetIDName[ijetid].Contains("N")){
        jetptjetidEtaBin[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 100, 0., 100.);    //Added
        jetptjetidEtaBin_real[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_real",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_real",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 100, 0., 100.);    //Added
        jetptjetidEtaBin_fake[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_fake",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_fake",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 100, 0., 100.);    //Added
        jetptjetidEtaBin_real1[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_real1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_real1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 100, 0., 100.);    //Added
        jetptjetidEtaBin_fake1[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_fake1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_fake1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 100, 0., 100.);    //Added
       }
        else{
        jetptjetidEtaBin[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 3000, 0., 300.);   //Added
        jetptjetidEtaBin_fake[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_fake",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_fake",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 3000, 0., 300.);   //Added
        jetptjetidEtaBin_real[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_real",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_real",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 3000, 0., 300.);   //Added
        jetptjetidEtaBin_fake1[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_fake1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_fake1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 3000, 0., 300.);   //Added
        jetptjetidEtaBin_real1[ieta][ijetid] = new TH2F(Form("jetpt%sEtaBin%.f_%.f_real1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), Form("jetpt%sEtaBin%.f_%.f_real1",JetIDName[ijetid].Data(),deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 3000, 0., 300.);   //Added
        }
        jetptjetidEtaBin[ieta][ijetid]->Sumw2();      //Added
        jetptjetidEtaBin_fake[ieta][ijetid]->Sumw2();      //Added
        jetptjetidEtaBin_real[ieta][ijetid]->Sumw2();      //Added
        jetptjetidEtaBin_fake1[ieta][ijetid]->Sumw2();      //Added
        jetptjetidEtaBin_real1[ieta][ijetid]->Sumw2();      //Added
       }
      refjetptEtaBin[ieta] = new TH1F(Form("refjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("refjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      refjetptEtaBin[ieta]->Sumw2();
      genjetptEtaBin[ieta] = new TH1F(Form("genjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("genjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      genjetptEtaBin[ieta]->Sumw2();
      refptJESEtaBin[ieta] = new TH2F(Form("refptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("refptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      refptJESEtaBin[ieta]->Sumw2();
      rawptJESEtaBin[ieta] = new TH2F(Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[ieta]->Sumw2();
      }
}

void hist_class::Write()
{
  TString dataType;
  TString out_name;
  dataType = "MC";
  out_name=Form("%s%s%s_useskim.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/tmp/xuq7/%s",out_name.Data()),"RECREATE");  
    refjetpt->Write();
    genjetpt->Write();
    jetpt->Write();
    jetpt_jetidclosure->Write();
    jetpt_fake->Write();
    jetpt_fake->Write();
    jetpt_real->Write();
    jetpt_fake1->Write();
    jetpt_real1->Write();
    jetptEta->Write();
    jetptphi->Write();
    jetEtaphi->Write();
    Spikeetaphi->Write();
    Spikeetapt->Write();
    Spikeetaphi_real->Write();
    Spikeetapt_real->Write();
    Spikeetaphi_fake->Write();
    Spikeetapt_fake->Write();
    for(int ipt=0;ipt<Nbin_pt_coarse;ipt++){ 
    jetid12jetid14[ipt] -> Write();
    jetid12jetid14_fake1[ipt] -> Write();
    jetid12jetid14_real1[ipt] -> Write();
    }
        for(int ijetid=0;ijetid<nJetID;ijetid++){
        jetptjetid[ijetid]->Write();
        jetptjetid_real[ijetid]->Write();
        jetptjetid_fake[ijetid]->Write();
        jetptjetid_real1[ijetid]->Write();
        jetptjetid_fake1[ijetid]->Write();
        }
    rawptJES->Write();
    refptJES->Write();
    Vz->Write();
    VzW->Write();
    Cent->Write();
    CentW->Write();
    Pthat->Write();
    PthatAfXw->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[ieta]->Write(); 
	 jetptEtaBin_jetidclosure[ieta]->Write(); 
         jetptEtaBin_fake[ieta]->Write();
         jetptEtaBin_real[ieta]->Write();
         jetptEtaBin_fake1[ieta]->Write();
         jetptEtaBin_real1[ieta]->Write();
	 for(int ijetid=0;ijetid<nJetID;ijetid++){
                 jetptjetidEtaBin[ieta][ijetid]->Write();
                 jetptjetidEtaBin_fake[ieta][ijetid]->Write();
                 jetptjetidEtaBin_real[ieta][ijetid]->Write();
                 jetptjetidEtaBin_fake1[ieta][ijetid]->Write();
                 jetptjetidEtaBin_real1[ieta][ijetid]->Write();
        }
	 genjetptEtaBin[ieta]->Write(); 
	 refjetptEtaBin[ieta]->Write(); 
	 refptJESEtaBin[ieta]->Write(); 
	 rawptJESEtaBin[ieta]->Write(); 
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void PPbanalyzeRpAskimTree()
{
  std::cout<<"start working\n";
  hist_class *my_hists = new hist_class();

  cout<<"Analyzing MC!"<<endl;
 TString user = getenv("USER");
 // TFile *f = new TFile(Form("/cms/store/user/ymao/pA5TEV/Mixing/STARTHI53V27/merged/%sMCOfficialForestNewVzWeightAddHLT_ppReco_akPu3PF_QCDjetTrigJECv8_JetPt0pthatLowerCut.root",coll.Data()));
  TFile *f = new TFile(Form("/cms/store/user/qixu/jetRpA/skimTree/MC%s%sskimfile0_10final.root",coll.Data(),algo.Data()));
 // TFile *f = new TFile(Form("/cms/store/user/qixu/jetRpA/skimTree/MC%s%sskimfile0_10.root",coll.Data(),algo.Data()));
 // TFile *f = new TFile(Form("/cms/store/user/ymao/pA5TEV/Mixing/STARTHI53V27/merged/MC%s%sskimFullInfoLowerpthatCutfile0_10.root",coll.Data(),algo.Data()));
  
  TTree *nt = (TTree*)f->Get("nt");
  TFile *fjetidcorr = TFile::Open(Form("../JetID/fjetidcorr%s_self.root",coll.Data()));
  TH1D* hjetidcorr = (TH1D*)fjetidcorr->Get("jetidcorr_-10_10");
    TH1D* hjetidcorrEtaBin[netabin];
  for(int ieta=0;ieta<netabin;ieta++){
    hjetidcorrEtaBin[ieta] = (TH1D*)fjetidcorr->Get(Form("jetidcorr_%.f_%.f",deta[ieta]*10,deta[ieta+1]*10));
  }

  Float_t jtpt[100],jtpu[100],jteta[100],geneta[100],genphi[100],jtphi[100],rawpt[100],refpt[100],refeta[100],genpt[100];
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
//  nt->SetBranchAddress("jtpu",jtpu);
  nt->SetBranchAddress("refpt",refpt);
  nt->SetBranchAddress("refeta",refeta);
  nt->SetBranchAddress("genpt",genpt);
  nt->SetBranchAddress("subid",subid);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("geneta",geneta);
//  nt->SetBranchAddress("genphi",genphi);
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
    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
if(TMath::Abs(vz)>15) continue;
 if(!HLT_PAJet20_noJetID_v1 && !HLT_PAJet40_noJetID_v1 && !HLT_PAJet60_noJetID_v1 && !HLT_PAJet80_noJetID_v1 && !HLT_PAJet100_noJetID_v1 ) continue;
        
	my_hists->Vz->Fill(vz,xSecWeight);
	my_hists->VzW->Fill(vz,xSecWeight*vzWeight);
	my_hists->Cent->Fill(hiBin,xSecWeight*vzWeight);
	my_hists->CentW->Fill(hiBin,weight);
	my_hists->Pthat->Fill(pthat);
	my_hists->PthatAfXw->Fill(pthat,xSecWeight);

for(int j4i = 0; j4i < nref; j4i++){
    double jet_pt = jtpt[j4i];
 //   double jet_pu = jtpu[j4i];
    double ref_pt = refpt[j4i];
    double ref_eta = refeta[j4i];
    double jet_eta = jteta[j4i];  
    double jet_phi = jtphi[j4i]; 
    double raw_pt = rawpt[j4i];
	if(raw_pt<22.0 || fabs(jet_eta)>5.0) continue;
        if(jet_pt>4*pthat) continue;
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
      //  if((chargedN == 0 || chargedSum == 0) && TMath::Abs(jet_eta)< 2.4) continue;  // jet id selection
	double PPTighter0 = (double)(neutralSum/jet_pt < 0.8 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 1.0 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4));
	double PPTighter1 = (double)(neutralSum/jet_pt < 0.8 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 1.0 && chargedSum/jet_pt < 0.95 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
	double PPTighter2 = (double)(neutralSum/jet_pt < 0.9 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 1.0 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
	double PPTighter3 = (double)(neutralSum/jet_pt < 0.9 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 1.0 && chargedSum/jet_pt < 0.95 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
	double PPTighter4 = (double)(neutralSum/jet_pt < 0.9 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 0.95 && chargedSum/jet_pt < 0.95 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
	double PPTighter5 = (double)(neutralSum/jet_pt < 0.8 && eSum/jet_pt < 1.0 && photonSum/jet_pt < 1.0 && (chargedSum+neutralSum+muSum+eSum)/jet_pt < 0.95 && chargedSum/jet_pt < 0.95 && ((chargedSum/jet_pt>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) );
        double PPTighter = PPTighter0*TMath::Power(2,5)+PPTighter1*TMath::Power(2,4)+PPTighter2*TMath::Power(2,3)+PPTighter3*TMath::Power(2,2)+PPTighter4*TMath::Power(2,1)+PPTighter5*TMath::Power(2,0)+0.5;
	double jetidv[nJetID]={chargedMax,chargedSum,neutralMax,neutralSum,photonMax,photonSum,chargedMax/jet_pt,chargedSum/jet_pt,neutralMax/jet_pt,neutralSum/jet_pt,photonMax/jet_pt,photonSum/jet_pt,eSum/jet_pt,(chargedSum+neutralSum+photonSum+muSum+eSum)/jet_pt,(chargedSum+neutralSum+photonSum+muSum+eSum)/raw_pt,neutralMax/TMath::Max(chargedSum,neutralSum),(double)chargedN,(double)neutralN,(double)photonN,(double)(neutralSum/jet_pt<1.0 && eSum/jet_pt<1.0 && photonSum/jet_pt<1.0 && ((chargedSum>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) ), (double)(neutralSum/jet_pt<0.9 && eSum/jet_pt<1.0 && photonSum/jet_pt<0.9 && ((chargedSum>0 && chargedN>0 && TMath::Abs(jet_eta)<2.4) || TMath::Abs(jet_eta) >=2.4) ),PPTighter,(chargedSum+neutralSum+muSum+eSum)/jet_pt};
        if(TMath::Abs(jet_eta)<=2.4){
                my_hists->jetptEta->Fill(jet_pt,jet_eta,weight);
                my_hists->jetptphi->Fill(jet_pt,jet_phi,weight);
                my_hists->jetEtaphi->Fill(jet_eta,jet_phi,weight);
        }

        int matchflag=0;
        int matchflag1=0;
        if(TMath::Abs(jet_eta)<=2){
	for(int j5i = 0; j5i < ngen ; j5i++) {
        double gen_pt=genpt[j5i];
//        double gen_eta=geneta[j5i];
//        double gen_phi=genphi[j5i];
        TVector3 jetVec, genVec;
        jetVec.SetPtEtaPhi(jet_pt,jet_eta,jet_phi);
        genVec.SetPtEtaPhi(jet_pt,jet_eta,jet_phi);
        double deltaR = jetVec.DeltaR(genVec);
        if(deltaR<0.3){ 
        matchflag++;
        }
        if(deltaR<0.3 && TMath::Abs(jet_pt-gen_pt)/jet_pt < 0.5){
           matchflag1++; 
        }
        }
        }
        
        if(neutralSum/jet_pt > 0.8) {
            my_hists->Spikeetaphi->Fill(jet_eta,jet_phi, weight);
            my_hists->Spikeetapt->Fill(jet_eta,jet_pt, weight);
	if(matchflag1==0){	//Added
            my_hists->Spikeetaphi_fake->Fill(jet_eta,jet_phi, weight);
            my_hists->Spikeetapt_fake->Fill(jet_eta,jet_pt, weight);
        }
        else{
            my_hists->Spikeetaphi_real->Fill(jet_eta,jet_phi, weight);
            my_hists->Spikeetapt_real->Fill(jet_eta,jet_pt, weight);
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

	if(TMath::Abs(ref_eta) <= 1.) {
        my_hists->refjetpt->Fill(ref_pt, weight);
       }
    int dEtaBin = -1. ;
        for(Int_t ieta = 0 ; ieta <netabin; ieta++){
               if(ref_eta>deta[ieta]&&ref_eta<=deta[ieta+1]) dEtaBin = ieta ;
         } //assign the eta bin for jets
      if(dEtaBin!=-1){
        my_hists->refjetptEtaBin[dEtaBin]->Fill(ref_pt, weight);
       	my_hists->refptJES->Fill(ref_pt,jet_pt/ref_pt,weight);
       }


      dEtaBin = -1.;
     if(TMath::Abs(jet_eta)<=1.) {
    for(int ipt=0;ipt<Nbin_pt_coarse;ipt++){ 
        if(jet_pt>binbound_pt_coarse[ipt] && jet_pt <= binbound_pt_coarse[ipt+1]){
	  my_hists->jetid12jetid14[ipt]->Fill(jetidv[12], jetidv[14], weight);	//Added
	if(matchflag1==0)	//Added
	  my_hists->jetid12jetid14_fake1[ipt]->Fill(jetidv[12], jetidv[14], weight);	//Added
        else
	  my_hists->jetid12jetid14_real1[ipt]->Fill(jetidv[12], jetidv[14], weight);	//Added
        }
    }
	for(int ijetid=0;ijetid<nJetID;ijetid++){
	my_hists->jetptjetid[ijetid]->Fill(jet_pt, jetidv[ijetid], weight);	//Added
	if(ref_pt<0)	//Added
	  my_hists->jetptjetid_fake[ijetid]->Fill(jet_pt, jetidv[ijetid], weight);	//Added
	else
	  my_hists->jetptjetid_real[ijetid]->Fill(jet_pt, jetidv[ijetid], weight);	//Added
	if(matchflag1==0)	//Added
	  my_hists->jetptjetid_fake1[ijetid]->Fill(jet_pt, jetidv[ijetid], weight);	//Added
	else
	  my_hists->jetptjetid_real1[ijetid]->Fill(jet_pt, jetidv[ijetid], weight);	//Added
        }
	if(ref_pt<0)	//Added
          my_hists->jetpt_fake->Fill(jet_pt,weight);
	else
          my_hists->jetpt_real->Fill(jet_pt,weight);
	if(matchflag1==0)	//Added
          my_hists->jetpt_fake1->Fill(jet_pt,weight);
	else
          my_hists->jetpt_real1->Fill(jet_pt,weight);
	my_hists->jetpt->Fill(jet_pt,weight);
        if((int)PPTighter2)
	my_hists->jetpt_jetidclosure->Fill(jet_pt,weight*hjetidcorr->GetBinContent(hjetidcorr->FindBin(jet_pt)));
       	my_hists->rawptJES->Fill(raw_pt,jet_pt/raw_pt,weight);
	}	

      for(int ieta =0; ieta<netabin;ieta++){
          if(jet_eta>deta[ieta]&&jet_eta<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 my_hists->jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
        if((int)PPTighter2)
	my_hists->jetptEtaBin_jetidclosure[dEtaBin]->Fill(jet_pt,weight*hjetidcorrEtaBin[dEtaBin]->GetBinContent(hjetidcorrEtaBin[dEtaBin]->FindBin(jet_pt)));
	for(int ijetid=0;ijetid<nJetID;ijetid++){
            my_hists->jetptjetidEtaBin[dEtaBin][ijetid]->Fill(jet_pt,jetidv[ijetid],weight);
            if(ref_pt<0)
                my_hists->jetptjetidEtaBin_fake[dEtaBin][ijetid]->Fill(jet_pt,jetidv[ijetid],weight);
            else
                my_hists->jetptjetidEtaBin_real[dEtaBin][ijetid]->Fill(jet_pt,jetidv[ijetid],weight);
            if(matchflag1==0)
                my_hists->jetptjetidEtaBin_fake1[dEtaBin][ijetid]->Fill(jet_pt,jetidv[ijetid],weight);
            else
                my_hists->jetptjetidEtaBin_real1[dEtaBin][ijetid]->Fill(jet_pt,jetidv[ijetid],weight);
        }
         if(ref_pt<0)
            my_hists->jetptEtaBin_fake[dEtaBin]->Fill(jet_pt,weight);
         else
            my_hists->jetptEtaBin_real[dEtaBin]->Fill(jet_pt,weight);
         if(matchflag1==0)
            my_hists->jetptEtaBin_fake1[dEtaBin]->Fill(jet_pt,weight);
         else
            my_hists->jetptEtaBin_real1[dEtaBin]->Fill(jet_pt,weight);
	 my_hists->rawptJESEtaBin[dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,weight);
       	 my_hists->refptJESEtaBin[dEtaBin]->Fill(ref_pt,jet_pt/ref_pt,weight);
	}
      } //loop over jet
	
	for(int j5i = 0; j5i < ngen ; j5i++) {
    double gen_pt=genpt[j5i];
    double gen_eta=geneta[j5i];
        if(coll =="PPb") 
            gen_eta = gen_eta+0.465;
        if(coll =="PbP") 
            gen_eta = gen_eta-0.465;
        if(TMath::Abs(gen_eta)<=1.) {
          my_hists->genjetpt->Fill(gen_pt, weight);
       }
      int  dEtaBin=-1;
        for(Int_t ieta = 0 ; ieta <netabin; ieta++){
               if(gen_eta>deta[ieta]&&gen_eta<=deta[ieta+1]) dEtaBin = ieta ;
         } //assign the eta bin for jets
      if(dEtaBin!=-1){
        my_hists->genjetptEtaBin[dEtaBin]->Fill(gen_pt, weight);
        }
	} //end of gen jet loop
}	//loop over tree

  my_hists->Write();
  cout<<"working done\n";
 
}  
