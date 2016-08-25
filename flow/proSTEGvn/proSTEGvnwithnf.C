#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "SetSeedOwn.C"
using namespace std;

const Int_t maxnh = 20000;
Int_t b_npg, b_n, n;
Float_t b_phirg;
Float_t b_ptg[maxnh], b_etag[maxnh], b_phig[maxnh];

void proSTEGvnwithnf()
{
  
  int mult = atoi(getenv("MULT"));

  int tot_num=50000;  //50k events
  double MeanMult=(double)mult;
  double RMSMult=10;
  
  int ifile = atoi(getenv("IFILE")); 

  //simple toy event generator
  //TFile f(Form("/lio/lfs/cms/store/user/qixu/flow/NewSTEG/pPbDataV205m%d/vndata_50k_%d.root",mult,ifile), "RECREATE","ROOT file with histograms & tree");
  TFile f(Form("/tmp/xuq7/pPbDataV205m%d/vndata_50k_%d.root",mult,ifile), "RECREATE","ROOT file with histograms & tree");
  //TFile f(Form("vndata_50k_%d.root",mult,ifile), "RECREATE","ROOT file with histograms & tree");
  TTree *tree = new TTree("tree","Event tree with a few branches");
//  tree->Branch("npg", &b_npg, "npg/I");   // # of particles;
//  tree->Branch("phirg", &b_phirg, "phirg/F");  // RP angle;
  tree->Branch("n", &b_n, "n/I");          // same as npg;
  tree->Branch("ptg", &b_ptg, "ptg[n]/F");  // ;
  tree->Branch("etag", &b_etag, "etag[n]/F");
  tree->Branch("phig", &b_phig, "phig[n]/F");
  
  TF1 *EtaDistr = new TF1("EtaDistr","exp(-(x-2.1)^2/6.3)+exp(-(x+2.1)^2/6.3)",-4,4);
  //TF1 *PhiDistr = new TF1("PhiDistr","1+2*[0]*cos(x)+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x)+2*[4]*cos(5*x)+2*[5]*cos(6*x)",0,2.*TMath::Pi());
  TF1 *PhiDistr = new TF1("PhiDistr","1+2*[0]*cos(2*x)+2*[1]*cos(3*x)+2*[2]*cos(4*x)+2*[3]*cos(5*x)+2*[4]*cos(6*x)",0,2.*TMath::Pi());
  //TF1 *PtDistr  = new TF1("PtDistr","exp (-(x/.40))+0.0015*exp (-(x/1.5))", 0.2,10);	//V~0.12
  //TF1 *PtDistr  = new TF1("PtDistr","exp (-(x/0.90))+0.15*exp (-(x/15))", 0.1,10);	//V~=0.06
  TF1 *PtDistr  = new TF1("PtDistr","0.03*(exp (-(x/0.594540))+0.00499506*exp (-(x/1.89391)))", 0.3,6.0);	//Real Data
  //  TF1 *PtDistr = new TF1("PtDistr","[0]*x*TMath::Power(1+(sqrt(x*x+[1]*[1])-[1]/[2]),-[3])",0.2,10);
	//PtDistr->SetParameters(118.836,-0.335972,0.759243,118.836);	//Real data fit with Tsallis
  //TF1 *V1vsEta = new TF1("V1vsEta","-exp(-(x+1)^2)/20-x/30+exp(-(x-1)^2)/20",-2.4,2.4); 
  //TF1 *V2vsPt   = new TF1("V2vsPt","((x/3)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);
  //TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.1,10);
 //	V2vsPt->SetParameters(4.81159,1.80783,3.69272,3.11889,0.931485);	//Real data V~0.05
  //	V2vsPt->SetParameters(5,1.8,3,1.8,0.8); //V~0.06
  TF1 *V2vsPt = new TF1("V2vsPt","((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);
  TF1 *V3vsPt = new TF1("V3vsPt","((x/3.2)^2.3/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.3,6.0);
  TF1 *V4vsPt = new TF1("V4vsPt","((x/4.8)^2.1/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.3,6.0);
  TF1 *V5vsPt = new TF1("V5vsPt","((x/6.0)^3.2/(1+(x/11.4)^2.1))*(.00005+(1/x)^1.4)",0.3,6.0);
  TF1 *V6vsPt = new TF1("V6vsPt","((x/5.6)^2.4/(1+(x/4.7)^2.1))*(.00005+(1/x)^1.4)",0.3,6.0);
 
    UInt_t iniseed = SetSeedOwn();
    gRandom->SetSeed(iniseed);
  TRandom3 *rnd = new TRandom3(0);
  
  double v1, v2, v3, v4, v5, v6, ph, myphi, mypt, myeta, phi0, Psi;
  int nnf,k;
  double neta, nphi, npt;
  int slicept;
  
  for(int i=0; i<tot_num; i++){ 
    
    Psi = rnd->Uniform(0.0,2.*TMath::Pi());
    //Psi=0;
    b_phirg = Psi; 
    b_npg = rnd->Gaus(MeanMult,RMSMult); 
    int b_npgsame = (int)b_npg * 0.05;
    int b_npgback = (int)b_npg * 0.05;
    n = 0;
  
    for(int j=0; j<b_npg;j++ ){
      mypt = PtDistr->GetRandom();
      myeta =  EtaDistr->GetRandom();

      //v1=V1vsEta->Eval(myeta);
      v2=V2vsPt->Eval(mypt);
      v3=V3vsPt->Eval(mypt);
      v4=V4vsPt->Eval(mypt);
      v5=V5vsPt->Eval(mypt);
      v6=V6vsPt->Eval(mypt);

      b_etag[n] = myeta;
      b_ptg[n]  = mypt;
      //PhiDistr->SetParameters(v1,v2,v3,v4,v5,v6);
      PhiDistr->SetParameters(v2,v3,v4,v5,v6);
      
      myphi = PhiDistr->GetRandom(); // randon selection dn/dphi

      myphi = myphi+Psi; // angle in lab frame -- needed for correct cumulant v2
      if (myphi>2.*TMath::Pi()) myphi=myphi-2.*TMath::Pi(); // 0 - 2*Pi
      
      b_phig[n] = myphi; // save angle in lab frame
      n++;
    }
    for(int j=0;j<b_npgsame;j++){
      int k = rnd->Integer(n);
      b_ptg[n] = b_ptg[k];
      b_etag[n] = b_etag[k]+2*(k%2-0.5)*1e-3;
      b_phig[n] = b_phig[k]+2*(k%2-0.5)*1e-3;
      n++;
    }
    for(int j=0;j<b_npgback;j++){
      int k = rnd->Integer(n-b_npgsame);
      b_ptg[n] = b_ptg[k];
      b_etag[n] = -(b_etag[k]+2*(k%2-0.5)*1e-3);
      b_phig[n] = TMath::Pi()+(b_phig[k]+2*(k%2-0.5)*1e-3);
      if (b_phig[n]>2.*TMath::Pi()) b_phig[n]=b_phig[n]-2.*TMath::Pi(); // 0 - 2*Pi
      n++;
    } // End of loop over particles
 
    if (i%10000 == 0) cout << i << " " << "events processed" << endl;

    b_n = n;
    tree->Fill();
  } // End of loop over events
  
  cout << "writing tree" << endl;
  tree->Write();
  cout << "writing to file" << endl;
  f.Write();
  cout << "closing file" << endl;
  f.Close();
  cout << "THE END" << endl;
}

