#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include "iomanip"

#define PI 3.14159265359
using namespace std; 

//Correlator();
/*int main(){
  Correlator();
}
*/
void Correlator(){

using namespace std;

   vector<TVector3> *pVect_trg;
   vector< vector<TVector3> > *pVectVect_trg;
   vector<TVector3> *pVect_ass;
   vector< vector<TVector3> > *pVectVect_ass;
   vector<double> *zvtxVect;

  TH1::SetDefaultSumw2();
  TH1D* hMult = new TH1D("mult","N",250,0,250);
  TH2D* hSignal = new TH2D("signal","#Delta#eta,#Delta#phi",32+1,-4.8-0.15,4.8+0.15,32-1,-PI/2+PI/32,3*PI/2-PI/32);
  TH2D* hBackground = new TH2D("background","#Delta#eta, #Delta#phi",32+1,-4.8-0.15,4.8+0.15,32-1,-PI/2+PI/32,3*PI/2-PI/32);
  TH2D* hCorrelation = new TH2D("correlation","#Delta#eta, #Delta#phi",32+1,-4.8-0.15,4.8+0.15,32-1,-PI/2+PI/32,3*PI/2-PI/32);

  pVectVect_trg = new vector< vector<TVector3> >;
  pVectVect_ass = new vector< vector<TVector3> >;
  zvtxVect = new vector<double>;


   double etaMin_trg_=-2.4;
   double etaMax_trg_=2.4;
   double etaMin_ass_=-2.4;
   double etaMax_ass_=2.4;
   double ptMin_trg_=1.0;
   double ptMax_trg_=10.0;
   double ptMin_ass_=1.0;
   double ptMax_ass_=10.0;
   int bkgFactor_=10;

  Int_t rmult=0;
  Float_t phir0;
  Float_t phi0[20000];
  Float_t pt0[20000];
  Float_t eta0[20000];

     char *file1 = "/store/user/tuos/tmp/STEGfixM/STEGv2onlyFixM100_50k.root";
     TFile *f1 = new TFile(file1);
     TTree *t1 = (TTree*)f1->Get("tree");
     t1->SetBranchAddress("npg", &rmult);
     t1->SetBranchAddress("ptg", &pt0);
     t1->SetBranchAddress("etag", &eta0);
     t1->SetBranchAddress("phig", &phi0);
     t1->SetBranchAddress("phirg", &phir0);
     long neve=t1->GetEntries();
     cout<<"Nevts = "<<neve<<endl;
     for(long ne=0; ne<neve; ne++){
       if(ne%1000==0)  cout<<"Have run "<<ne<<" events ,"<<" First pass"<<endl;
       t1->GetEntry(ne);

       pVect_trg = new vector<TVector3>;
       pVect_ass = new vector<TVector3>;

       for(int nh=0;nh<rmult;nh++){
         double eta = eta0[nh];
         double phi = phi0[nh];
         double pt  = pt0[nh];

         TVector3 pvector;
         pvector.SetPtEtaPhi(pt,eta,phi);
         if(eta0[nh]<=etaMax_trg_ && eta0[nh]>=etaMin_trg_ && pt0[nh]<=ptMax_trg_ && pt0[nh]>=ptMin_trg_) pVect_trg->push_back(pvector);
         if(eta0[nh]<=etaMax_ass_ && eta0[nh]>=etaMin_ass_ && pt0[nh]<=ptMax_ass_ && pt0[nh]>=ptMin_ass_) pVect_ass->push_back(pvector);
    }
//cout<<"work here"<<endl;

    int nMult_trg = (int)pVect_trg->size();
    int nMult_ass = (int)pVect_ass->size();
    hMult->Fill(nMult_trg);

    for(int ntrg=0;ntrg<nMult_trg;ntrg++){
      TVector3 pvector_trg = (*pVect_trg)[ntrg];
      double eta_trg = pvector_trg.Eta();
      double phi_trg = pvector_trg.Phi();

      for(int nass=0;nass<nMult_ass;nass++){
        TVector3 pvector_ass = (*pVect_ass)[nass];
        double eta_ass = pvector_ass.Eta();
        double phi_ass = pvector_ass.Phi();

        double deltaEta=eta_ass-eta_trg;
        double deltaPhi=phi_ass-phi_trg;
        if(deltaPhi>PI)
          deltaPhi=deltaPhi-2*PI;
        if(deltaPhi<-PI)
          deltaPhi=deltaPhi+2*PI;
        if(deltaPhi>-PI && deltaPhi<-PI/2.)
          deltaPhi=deltaPhi+2*PI;

        if(deltaEta==0 && deltaPhi==0) continue;
        hSignal->Fill(fabs(deltaEta),fabs(deltaPhi),1.0/4.0/nMult_trg);
        hSignal->Fill(-fabs(deltaEta),fabs(deltaPhi),1.0/4.0/nMult_trg);
        hSignal->Fill(fabs(deltaEta),-fabs(deltaPhi),1.0/4.0/nMult_trg);
        hSignal->Fill(-fabs(deltaEta),-fabs(deltaPhi),1.0/4.0/nMult_trg);
        hSignal->Fill(fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0/nMult_trg);
        hSignal->Fill(-fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0/nMult_trg);
      }
     } // end of loop over particles

    pVectVect_trg->push_back(*pVect_trg);
    pVectVect_ass->push_back(*pVect_ass);
    //zvtxVect->push_back(bestvz);

    delete pVect_trg;
    delete pVect_ass;

     } // end of loop over events


  int nevttotal_trg = (int)pVectVect_trg->size();
  int nevttotal_ass = (int)pVectVect_ass->size();

  for(int nround=0;nround<bkgFactor_;nround++)
  {
      cout<<"round = "<<nround<<endl;
    for(int nevt_ass=0; nevt_ass<nevttotal_ass; nevt_ass++)
    {
      int nevt_trg = gRandom->Integer(nevttotal_trg);
      if(nevt_trg == nevt_ass) { nevt_ass--; continue; }
      //if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) { nevt_ass--; continue; }

      vector<TVector3> pVectTmp_trg = (*pVectVect_trg)[nevt_trg];
      vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
      int nMult_trg = pVectTmp_trg.size();
      int nMult_ass = pVectTmp_ass.size();

      for(int ntrg=0;ntrg<nMult_trg;ntrg++)
      {
        TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
        double eta_trg = pvectorTmp_trg.Eta();
        double phi_trg = pvectorTmp_trg.Phi();
        for(int nass=0;nass<nMult_ass;nass++)
        {
          TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
          double eta_ass = pvectorTmp_ass.Eta();
          double phi_ass = pvectorTmp_ass.Phi();

          double deltaEta=eta_ass-eta_trg;
          double deltaPhi=phi_ass-phi_trg;
          if(deltaPhi>PI)
            deltaPhi=deltaPhi-2*PI;
          if(deltaPhi<-PI)
            deltaPhi=deltaPhi+2*PI;
          if(deltaPhi>-PI && deltaPhi<-PI/2.)
            deltaPhi=deltaPhi+2*PI;

          if(deltaEta==0 && deltaPhi==0) continue;
          hBackground->Fill(fabs(deltaEta),fabs(deltaPhi),1.0/nMult_trg/nMult_ass);
          hBackground->Fill(-fabs(deltaEta),fabs(deltaPhi),1.0/nMult_trg/nMult_ass);
          hBackground->Fill(fabs(deltaEta),-fabs(deltaPhi),1.0/nMult_trg/nMult_ass);
          hBackground->Fill(-fabs(deltaEta),-fabs(deltaPhi),1.0/nMult_trg/nMult_ass);
          hBackground->Fill(fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/nMult_trg/nMult_ass);
          hBackground->Fill(-fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/nMult_trg/nMult_ass);
        }
      }
    }
  }


//  int nEvent = hMult->Integral(3,10000);
//  hSignal->Scale(1.0/nEvent);
//  hBackground->Scale(hSignal->Integral()/hBackground->Integral());
  hCorrelation->Add(hSignal);
//  hCorrelation->Add(hBackground,-1);
  hCorrelation->Divide(hBackground);
//  hCorrelation->Scale(hSignal->Integral());



  TFile *file = new TFile("correlationSTEG.root", "RECREATE");
  cout << "writing to file" << endl;
  hMult->Write();
  hSignal->Write();
  hBackground->Write();
  hCorrelation->Write();
  file->Write();
  cout << "closing file" << endl;
  file->Close();
  cout << "THE END" << endl;


}
 
