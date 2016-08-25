#include <stdio.h>
#include <complex>
#include <iostream>
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include <TMath.h>
#define NB_SLICE 1
#define NB_SLICE_PT 15
#define sqr(x) ((x)*(x))
#define Pi 3.14159265359 
//#define PI TMath::Pi()

 
using namespace std; 


EPAnalyzer(){

  ofstream outv2;
  ofstream output;
  ofstream outv3;
  ofstream output3;

  int tt,k1,k2,nh,Nh,islc,islcpt,inttemp;
  int s=0;
  Int_t rmult=0;
  int mult=0;
  long Totevt=0;
  long neve=0; 
  long Neve;
  double rpmult=0.;
  double rmultmean=0.;
  double Trmultmean=0.;  
  Float_t phi0[20000];
  Float_t pt0[20000];
  Float_t eta0[20000];
  Float_t vz0;

  //double w2=1.;
  double temp=0.;

  double rpmult2[NB_SLICE][NB_SLICE_PT];
  double totmult[NB_SLICE][NB_SLICE_PT]; 
  double EPR=0.;  double EPR1=0.;  double EPR15=0.;  double EPR2=0.;
  double v22mean[NB_SLICE][NB_SLICE_PT];
  double err22mean[NB_SLICE][NB_SLICE_PT];
  double EPR3=0.; double EPR31=0.; double EPR315=0.;  double EPR32=0.;
  double v33mean[NB_SLICE][NB_SLICE_PT];
  double err33mean[NB_SLICE][NB_SLICE_PT];
  double v22mean1[NB_SLICE][NB_SLICE_PT];
  double v33mean1[NB_SLICE][NB_SLICE_PT];
  double v22mean15[NB_SLICE][NB_SLICE_PT];
  double v33mean15[NB_SLICE][NB_SLICE_PT];
  double v22mean2[NB_SLICE][NB_SLICE_PT];
  double v33mean2[NB_SLICE][NB_SLICE_PT];

  double minptcut[NB_SLICE_PT] = {0.2,0.4,0.6,0.8,1.1,1.4,1.75,2.2,2.75,3.45,4.2,5.05,6.0,7.2,8.5};
  double maxptcut[NB_SLICE_PT] = {0.4,0.6,0.8,1.1,1.4,1.75,2.2,2.75,3.45,4.2,5.05,6.0,7.2,8.5,10.0}; //15

  double v2EPtemp;
  double v2EPCtemp;
  double v3EPtemp;
  double v3EPCtemp;

  double MeanPt[NB_SLICE_PT];
  double MofPt[NB_SLICE_PT];

  for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) {   
    MeanPt[islcpt]=0.;
    MofPt[islcpt]=0.;

    for (islc=0; islc<NB_SLICE; islc++) {
      v22mean[islc][islcpt]=0.;
      v22mean1[islc][islcpt]=0.;
      v22mean15[islc][islcpt]=0.;
      v22mean2[islc][islcpt]=0.;
      err22mean[islc][islcpt]=0.;
      v33mean[islc][islcpt]=0.;
      v33mean1[islc][islcpt]=0.;
      v33mean15[islc][islcpt]=0.;
      v33mean2[islc][islcpt]=0.;
      err33mean[islc][islcpt]=0.;
    }
    }

  TH1F * pt = new TH1F("pt", "pt, Toy Model ",200,0,20);
  TH1F * eta = new TH1F("eta", "eta, Toy Model ",200,-4,4);
  TH1F * phi = new TH1F("phi", "phi, Toy Model ",200,-4,4);
  TH1F * psiPM = new TH1F("psiPM", "psiPM, Toy Model ",200,-6,6);
  TH1F * psiPM3 = new TH1F("psiPM3", "psiPM3, Toy Model ",200,-6,6);

  TH1F * PsiEP = new TH1F("PsiEP", "PsiEP, Toy Model, 50k, <M>=500",200,-4,4);
  TH1F * v2EP = new TH1F("v2EP", "v2EP, Toy Model, 10k, <M>=500",400,-0.5,0.5);
  TH1F * v2EPC = new TH1F("v2EPC", "v2EPC, Toy Model, 10k, <M>=500",400,-0.5,0.5);
  v2EP->Sumw2(); v2EPC->Sumw2();
  TH1F * PsiEP3 = new TH1F("PsiEP3", "PsiEP, Toy Model, 50k, <M>=500",200,-4,4);
  TH1F * v3EP = new TH1F("v3EP", "v2EP, Toy Model, 10k, <M>=500",400,-0.5,0.5);
  TH1F * v3EPC = new TH1F("v3EPC", "v2EPC, Toy Model, 10k, <M>=500",400,-0.5,0.5);
  v3EP->Sumw2(); v3EPC->Sumw2();


  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);

  PsiEP->GetXaxis()->SetTitle("#Psi_{EP}");
  PsiEP->GetXaxis()->CenterTitle();
  v2EP->GetXaxis()->SetTitle("v_{2}{EP}");
  v2EP->GetXaxis()->CenterTitle();
  v2EPC->GetXaxis()->SetTitle("v_{2}{EPC}");
  v2EPC->GetXaxis()->CenterTitle();
  PsiEP3->GetXaxis()->SetTitle("#Psi_{EP3}");
  PsiEP3->GetXaxis()->CenterTitle();
  v3EP->GetXaxis()->SetTitle("v_{3}{EP}");
  v3EP->GetXaxis()->CenterTitle();
  v3EPC->GetXaxis()->SetTitle("v_{3}{EPC}");
  v3EPC->GetXaxis()->CenterTitle();

Int_t b_nh;
Float_t b_ep;
Float_t b_v2ep;
Float_t b_v2epc;

Float_t b_ep3;
Float_t b_v3ep;
Float_t b_v3epc;

  TTree *tree = new TTree("tree","Event tree with a few branches");
  tree->Branch("nh", &b_nh, "nh/I");
  tree->Branch("ep", &b_ep, "ep/F");
  tree->Branch("v2ep", &b_v2ep, "v2ep/F");
  tree->Branch("v2epc", &b_v2epc, "v2epc/F");
  tree->Branch("ep3", &b_ep3, "ep3/F");
  tree->Branch("v3ep", &b_v3ep, "v3ep/F");
  tree->Branch("v3epc", &b_v3epc, "v3epc/F");


//double EP=0; double sumsin=0; double sumcos=0; double Ssin; double Scos;
double EPp2=0; double sumsinp2=0; double sumcosp2=0; double Ssinp2; double Scosp2;
double EPm2=0; double sumsinm2=0; double sumcosm2=0; double Ssinm2; double Scosm2;
double EPp15=0; double sumsinp15=0; double sumcosp15=0; double Ssinp15; double Scosp15;
double EPm15=0; double sumsinm15=0; double sumcosm15=0; double Ssinm15; double Scosm15;
double EPp1=0; double sumsinp1=0; double sumcosp1=0; double Ssinp1; double Scosp1;
double EPm1=0; double sumsinm1=0; double sumcosm1=0; double Ssinm1; double Scosm1;
double EPp=0; double sumsinp=0; double sumcosp=0; double Ssinp; double Scosp;
double EPm=0; double sumsinm=0; double sumcosm=0; double Ssinm; double Scosm;
double phimpsi=0;
double phimpsi2=0;
double EP3p=0; double sumsin3p=0; double sumcos3p=0; double Ssin3p; double Scos3p;
double EP3m=0; double sumsin3m=0; double sumcos3m=0; double Ssin3m; double Scos3m;
double EP3p1=0; double sumsin3p1=0; double sumcos3p1=0; double Ssin3p1; double Scos3p1;
double EP3m1=0; double sumsin3m1=0; double sumcos3m1=0; double Ssin3m1; double Scos3m1;
double EP3p15=0; double sumsin3p15=0; double sumcos3p15=0; double Ssin3p15; double Scos3p15;
double EP3m15=0; double sumsin3m15=0; double sumcos3m15=0; double Ssin3m15; double Scos3m15;
double EP3p2=0; double sumsin3p2=0; double sumcos3p2=0; double Ssin3p2; double Scos3p2;
double EP3m2=0; double sumsin3m2=0; double sumcos3m2=0; double Ssin3m2; double Scos3m2;
double phimpsi3=0;
double phimpsi32=0;

     TChain t1("t1");
     t1.Add("../../proSTEGvn/100_test_events/vndata.root/tree");
     t1.SetBranchAddress("npg", &rmult);
     t1.SetBranchAddress("ptg", &pt0);
     t1.SetBranchAddress("etag", &eta0);
     t1.SetBranchAddress("phig", &phi0);
     //t1.SetBranchAddress("vz", &vz0);
     neve=t1.GetEntries();
     Neve=0;

    TFile *f = new TFile("phipsi_dist.root","recreate");

     cout<<"Nevts = "<<neve<<endl;
     //for(long ne=0; ne<2000; ne++)
     for(long ne=0; ne<neve; ne++)
     {
       if(ne%1000==0)  cout<<"Have run "<<ne<<" of all "<<neve<<" events; "<<endl;
       t1.GetEntry(ne);

       //if(rmult<110) continue;
       Neve++;
       Nh=0;  
       sumsinp=0; sumcosp=0;  sumsinm=0; sumcosm=0;
       sumsinp1=0; sumcosp1=0;  sumsinm1=0; sumcosm1=0;
       sumsinp15=0; sumcosp15=0;  sumsinm15=0; sumcosm15=0;
       sumsinp2=0; sumcosp2=0;  sumsinm2=0; sumcosm2=0;
       v2EPtemp=0;
       v2EPCtemp=0;
       sumsin3p=0; sumcos3p=0; sumsin3m=0; sumcos3m=0;
       sumsin3p1=0; sumcos3p1=0;  sumsin3m1=0; sumcos3m1=0;
       sumsin3p15=0; sumcos3p15=0;  sumsin3m15=0; sumcos3m15=0;
       sumsin3p2=0; sumcos3p2=0;  sumsin3m2=0; sumcos3m2=0;
       v3EPtemp=0;
       v3EPCtemp=0;

       //rmultmean+=((float) rmult);
        for(nh=0;nh<rmult;nh++){
           //temp=pt0[nh];
            pt->Fill(pt0[nh]);
            phi->Fill(phi0[nh]);
            eta->Fill(eta0[nh]);
           if(eta0[nh]<-2.40||eta0[nh]>2.40) continue;
           if(pt0[nh]<0.2||pt0[nh]>10.0) continue;
           temp=1.0;
           if(eta0[nh]>0){
              sumsinp+=temp*sin(2*phi0[nh]);
              sumcosp+=temp*cos(2*phi0[nh]);
              sumsin3p+=temp*sin(3*phi0[nh]);
              sumcos3p+=temp*cos(3*phi0[nh]);
             }
           if(eta0[nh]<0){
              sumsinm+=temp*sin(2*phi0[nh]);
              sumcosm+=temp*cos(2*phi0[nh]);
              sumsin3m+=temp*sin(3*phi0[nh]);
              sumcos3m+=temp*cos(3*phi0[nh]);
             }
           if(eta0[nh]>1.0){
              sumsinp1+=temp*sin(2*phi0[nh]);
              sumcosp1+=temp*cos(2*phi0[nh]);
              sumsin3p1+=temp*sin(3*phi0[nh]);
              sumcos3p1+=temp*cos(3*phi0[nh]);
             }
           if(eta0[nh]<-1.0){
              sumsinm1+=temp*sin(2*phi0[nh]);
              sumcosm1+=temp*cos(2*phi0[nh]);
              sumsin3m1+=temp*sin(3*phi0[nh]);
              sumcos3m1+=temp*cos(3*phi0[nh]);
             }
           if(eta0[nh]>1.5){
              sumsinp15+=temp*sin(2*phi0[nh]);
              sumcosp15+=temp*cos(2*phi0[nh]);
              sumsin3p15+=temp*sin(3*phi0[nh]);
              sumcos3p15+=temp*cos(3*phi0[nh]);
             }
           if(eta0[nh]<-1.5){
              sumsinm15+=temp*sin(2*phi0[nh]);
              sumcosm15+=temp*cos(2*phi0[nh]);
              sumsin3m15+=temp*sin(3*phi0[nh]);
              sumcos3m15+=temp*cos(3*phi0[nh]);
             }
           if(eta0[nh]>2.0){
              sumsinp2+=temp*sin(2*phi0[nh]);
              sumcosp2+=temp*cos(2*phi0[nh]);
              sumsin3p2+=temp*sin(3*phi0[nh]);
              sumcos3p2+=temp*cos(3*phi0[nh]);
             }
           if(eta0[nh]<-2.0){
              sumsinm2+=temp*sin(2*phi0[nh]);
              sumcosm2+=temp*cos(2*phi0[nh]);
              sumsin3m2+=temp*sin(3*phi0[nh]);
              sumcos3m2+=temp*cos(3*phi0[nh]);
             }

          }
       //Ssinp1=sumsinp1; Scosm1=sumcosm1;
       EPp=1.0/2*atan2(sumsinp,sumcosp); // EP in (-pi/2, pi/2);
       EPm=1.0/2*atan2(sumsinm,sumcosm); // EP in (-pi/2, pi/2);
       EPp1=1.0/2*atan2(sumsinp1,sumcosp1); // EP in (-pi/2, pi/2);
       EPm1=1.0/2*atan2(sumsinm1,sumcosm1); // EP in (-pi/2, pi/2);
       EPp15=1.0/2*atan2(sumsinp15,sumcosp15); // EP in (-pi/2, pi/2);
       EPm15=1.0/2*atan2(sumsinm15,sumcosm15); // EP in (-pi/2, pi/2);
       EPp2=1.0/2*atan2(sumsinp2,sumcosp2); // EP in (-pi/2, pi/2);
       EPm2=1.0/2*atan2(sumsinm2,sumcosm2); // EP in (-pi/2, pi/2);
       EP3p=1.0/3*atan2(sumsin3p,sumcos3p); // EP3 in (-pi/3, pi/3);
       EP3m=1.0/3*atan2(sumsin3m,sumcos3m); // EP3 in (-pi/3, pi/3);
       EP3p1=1.0/3*atan2(sumsin3p1,sumcos3p1); // EP3 in (-pi/3, pi/3);
       EP3m1=1.0/3*atan2(sumsin3m1,sumcos3m1); // EP3 in (-pi/3, pi/3);
       EP3p15=1.0/3*atan2(sumsin3p15,sumcos3p15); // EP3 in (-pi/3, pi/3);
       EP3m15=1.0/3*atan2(sumsin3m15,sumcos3m15); // EP3 in (-pi/3, pi/3);
       EP3p2=1.0/3*atan2(sumsin3p2,sumcos3p2); // EP3 in (-pi/3, pi/3);
       EP3m2=1.0/3*atan2(sumsin3m2,sumcos3m2); // EP3 in (-pi/3, pi/3);

       PsiEP3->Fill(EP3p);
       EPR3+=cos(3*(EP3p-EP3m));
       EPR31+=cos(3*(EP3p1-EP3m1));
       EPR315+=cos(3*(EP3p15-EP3m15));
       EPR32+=cos(3*(EP3p2-EP3m2));
       b_ep3=EP3p; 

       psiPM->Fill(EPp-EPm);
       psiPM3->Fill(EP3p-EP3m);

       PsiEP->Fill(EPp);
       EPR+=cos(2*(EPp-EPm));
       EPR1+=cos(2*(EPp1-EPm1));
       EPR15+=cos(2*(EPp15-EPm15));
       EPR2+=cos(2*(EPp2-EPm2));
       b_ep=EPp;


       for(nh=0;nh<rmult;nh++)
	 {
           if(phi0[nh]>Pi) phi0[nh]=phi0[nh]-2*Pi; // phi in (-pi,pi)
           if(eta0[nh]<-2.40||eta0[nh]>2.40) continue;           
           if(pt0[nh]<0.2||pt0[nh]>10.0) continue;

           //EP=1.0/2*atan2(Ssin-temp*sin(2*phi0[nh]),Scos-temp*cos(2*phi0[nh]));     //Remove autocorrelations;  EP in (-pi/2,pi/2);
           //EP3=1.0/3*atan2(Ssin3-temp*sin(3*phi0[nh]),Scos3-temp*cos(3*phi0[nh]));     //Remove autocorrelations;  EP3 in (-pi/3,pi/3);

      for(islcpt=0;islcpt<NB_SLICE_PT;islcpt++)
        {
//          if((eta0[nh]>minptcut[islcpt]&&eta0[nh]<=maxptcut[islcpt])){
          if((pt0[nh]>minptcut[islcpt]&&pt0[nh]<=maxptcut[islcpt])){
             inttemp=islcpt;
          }
       }
      islcpt=inttemp;
      islc=0;
//      MeanPt[islcpt]+=eta0[nh];
      MeanPt[islcpt]+=pt0[nh];
      MofPt[islcpt]+=1.;
                  //phimpsi=phi0[nh]-EP;  // phimpsi in (-3*pi/2,3*pi/2);
                  //if(phimpsi<-Pi) phimpsi=phimpsi+2*Pi;
                  //if(phimpsi>Pi) phimpsi=phimpsi-2*Pi;    //phimpsi in (-pi,pi);
           if(eta0[nh]<0){           
             v2EPtemp+=cos(2.*(phi0[nh]-EPp));
             v2EPCtemp+=cos(2.*(phi0[nh]-EPp))/0.931942;
             v22mean[0][islcpt]+=cos(2.*(phi0[nh]-EPp));
             v22mean1[0][islcpt]+=cos(2.*(phi0[nh]-EPp1));
             v22mean15[0][islcpt]+=cos(2.*(phi0[nh]-EPp15));
             v22mean2[0][islcpt]+=cos(2.*(phi0[nh]-EPp2));
            }
           if(eta0[nh]>0){
             v2EPtemp+=cos(2.*(phi0[nh]-EPm));
             v2EPCtemp+=cos(2.*(phi0[nh]-EPm))/0.931942;
             v22mean[0][islcpt]+=cos(2.*(phi0[nh]-EPm));
             v22mean1[0][islcpt]+=cos(2.*(phi0[nh]-EPm1));
             v22mean15[0][islcpt]+=cos(2.*(phi0[nh]-EPm15));
             v22mean2[0][islcpt]+=cos(2.*(phi0[nh]-EPm2));
            }
           if(eta0[nh]<0){
             v3EPtemp+=cos(3.*(phi0[nh]-EP3p));
             v3EPCtemp+=cos(3.*(phi0[nh]-EP3p))/0.842135;
             v33mean[0][islcpt]+=cos(3.*(phi0[nh]-EP3p));
             v33mean1[0][islcpt]+=cos(3.*(phi0[nh]-EP3p1));
             v33mean15[0][islcpt]+=cos(3.*(phi0[nh]-EP3p15));
             v33mean2[0][islcpt]+=cos(3.*(phi0[nh]-EP3p2));
            }	     
           if(eta0[nh]>0){
             v3EPtemp+=cos(3.*(phi0[nh]-EP3m));
             v3EPCtemp+=cos(3.*(phi0[nh]-EP3m))/0.842135;
             v33mean[0][islcpt]+=cos(3.*(phi0[nh]-EP3m));
             v33mean1[0][islcpt]+=cos(3.*(phi0[nh]-EP3m1));
             v33mean15[0][islcpt]+=cos(3.*(phi0[nh]-EP3m15));
             v33mean2[0][islcpt]+=cos(3.*(phi0[nh]-EP3m2));
            }


	   Nh++;


	 } // end of loop over particles

       rmultmean+=((float) Nh);

       v2EPtemp=v2EPtemp/Nh;
       v2EPCtemp=v2EPCtemp/Nh;
       v2EP->Fill(v2EPtemp);
       v2EPC->Fill(v2EPCtemp);
       b_v2ep=v2EPtemp;
       b_v2epc=v2EPCtemp;

       v3EPtemp=v3EPtemp/Nh;
       v3EPCtemp=v3EPCtemp/Nh;
       v3EP->Fill(v3EPtemp);
       v3EPC->Fill(v3EPCtemp);
       b_v3ep=v3EPtemp;
       b_v3epc=v3EPCtemp;

       tree->Fill();


     } // end of loop over events
    output.open("output.dat");
    output3.open("output3.dat");

     Trmultmean=rmultmean;

     Totevt=Neve;

     cout<<"Totevt: "<<Totevt <<  " Trmult: " <<Trmultmean<<",  <Mult> = "<<Trmultmean/Totevt<<endl;
     //output<<"Totevt: "<<Totevt <<  " Trmult: " <<Trmultmean<<",  <Mult> = "<<Trmultmean/Totevt<<endl;

  rmultmean=Trmultmean/Totevt;

  //EPR=EPR/Totevt;  cout<<"EPR= "<<EPR<<endl;
  //EPR3=EPR3/Totevt;  cout<<"EPR3= "<<EPR3<<endl;
  if(EPR>0) EPR=sqrt(EPR/Totevt); else EPR=0.0001;  cout<<"EPR= "<<EPR<<endl;
  if(EPR3>0) EPR3=sqrt(EPR3/Totevt); else EPR3=0.0001;   cout<<"EPR3= "<<EPR3<<endl;
  if(EPR1>0) EPR1=sqrt(EPR1/Totevt); else EPR1=0.0001;   cout<<"EPR1= "<<EPR1<<endl;
  if(EPR31>0) EPR31=sqrt(EPR31/Totevt); else EPR31=0.0001;   cout<<"EPR31= "<<EPR31<<endl;
  if(EPR15>0) EPR15=sqrt(EPR15/Totevt); else EPR15=0.0001;   cout<<"EPR15= "<<EPR15<<endl;
  if(EPR315>0) EPR315=sqrt(EPR315/Totevt); else EPR315=0.0001;   cout<<"EPR315= "<<EPR315<<endl;
  if(EPR2>0) EPR2=sqrt(EPR2/Totevt); else EPR2=0.0001;   cout<<"EPR2= "<<EPR2<<endl;
  if(EPR32>0) EPR32=sqrt(EPR32/Totevt); else EPR32=0.0001;   cout<<"EPR32= "<<EPR32<<endl;

    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) {
      MeanPt[islcpt]/=MofPt[islcpt];
        //output<<"M of pt = "<<MofPt[islcpt]<<"   v22mean[0][islcpt] = "<<v22mean[0][islcpt]<<endl;
        output<<MofPt[islcpt]<<endl;
        output3<<MofPt[islcpt]<<endl;
      //cout<<"pt= "<<pt[islcpt]<<", Mean pt = "<<MeanPt[islcpt]<<endl;
      v22mean[0][islcpt]/=MofPt[islcpt];
      v22mean1[0][islcpt]/=MofPt[islcpt];
      v22mean15[0][islcpt]/=MofPt[islcpt];
      v22mean2[0][islcpt]/=MofPt[islcpt];
      err22mean[0][islcpt]=1.0/sqrt(2*MofPt[islcpt]);
      cout<<MeanPt[islcpt]<<" ep  "<<v22mean[0][islcpt]<<"   "<<err22mean[0][islcpt]<<endl;
      cout<<MeanPt[islcpt]<<" ep1  "<<v22mean1[0][islcpt]<<"   "<<err22mean[0][islcpt]<<endl;
      cout<<MeanPt[islcpt]<<" ep15  "<<v22mean15[0][islcpt]<<"   "<<err22mean[0][islcpt]<<endl;
      cout<<MeanPt[islcpt]<<" ep2  "<<v22mean2[0][islcpt]<<"   "<<err22mean[0][islcpt]<<endl;
        
      v33mean[0][islcpt]/=MofPt[islcpt];
      v33mean1[0][islcpt]/=MofPt[islcpt];
      v33mean15[0][islcpt]/=MofPt[islcpt];
      err33mean[0][islcpt]=1.0/sqrt(2*MofPt[islcpt]);
      cout<<MeanPt[islcpt]<<" ep3  "<<v33mean[0][islcpt]<<"   "<<err33mean[0][islcpt]<<endl;
      cout<<MeanPt[islcpt]<<" ep31  "<<v33mean1[0][islcpt]<<"   "<<err33mean[0][islcpt]<<endl;
      cout<<MeanPt[islcpt]<<" ep315  "<<v33mean15[0][islcpt]<<"   "<<err33mean[0][islcpt]<<endl;
      cout<<MeanPt[islcpt]<<" ep32  "<<v33mean2[0][islcpt]<<"   "<<err33mean[0][islcpt]<<endl;
    }
     output<<"Totevt: "<<Totevt <<  " Trmult: " <<Trmultmean<<",  <Mult> = "<<Trmultmean/Totevt<<endl;
     output<<"EPR= "<<EPR<<endl;
     output<<"EPR1= "<<EPR1<<endl;
     output<<"EPR15= "<<EPR15<<endl;
     output<<"EPR2= "<<EPR2<<endl;
      output.close();
     output3<<"Totevt: "<<Totevt <<  " Trmult: " <<Trmultmean<<",  <Mult> = "<<Trmultmean/Totevt<<endl;
     output3<<"EPR3= "<<EPR3<<endl;
     output3<<"EPR31= "<<EPR31<<endl;
     output3<<"EPR315= "<<EPR315<<endl;
     output3<<"EPR32= "<<EPR32<<endl;
      output3.close();

    outv2.open("v2.txt");
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv2<<MeanPt[islcpt]<<"   "<<v22mean[0][islcpt]/EPR<<"   "<<err22mean[0][islcpt]<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv2<<MeanPt[islcpt]<<"   "<<v22mean1[0][islcpt]/EPR1<<"   "<<err22mean[0][islcpt]<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv2<<MeanPt[islcpt]<<"   "<<v22mean15[0][islcpt]/EPR15<<"   "<<err22mean[0][islcpt]<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv2<<MeanPt[islcpt]<<"   "<<v22mean2[0][islcpt]/EPR2<<"   "<<err22mean[0][islcpt]<<endl;
      outv2.close();    

    outv3.open("v3.txt");
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv3<<MeanPt[islcpt]<<"   "<<v33mean[0][islcpt]/EPR3<<"   "<<err33mean[0][islcpt]<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv3<<MeanPt[islcpt]<<"   "<<v33mean1[0][islcpt]/EPR31<<"   "<<err33mean[0][islcpt]<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv3<<MeanPt[islcpt]<<"   "<<v33mean15[0][islcpt]/EPR315<<"   "<<err33mean[0][islcpt]<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv3<<MeanPt[islcpt]<<"   "<<v33mean2[0][islcpt]/EPR32<<"   "<<err33mean[0][islcpt]<<endl;
      outv3.close();    

  //TFile *f = new TFile("phipsi_dist.root","recreate");
  //f->cd();
  PsiEP->Write();
  v2EP->Write();
  v2EPC->Write();
  PsiEP3->Write();
  v3EP->Write();
  v3EPC->Write();

  pt->Write();
  phi->Write();
  eta->Write();
  psiPM->Write();
  psiPM3->Write();

  cout << "writing tree" << endl;
  tree->Write();
  cout << "writing to file" << endl;  
  f->Write();
  cout << "closing file" << endl;
  f->Close();
  cout << "THE END" << endl;

}
 
