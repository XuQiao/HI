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
#define Pi 3.14159265359 
#define NBin_PT 10

using namespace std; 
 
QCumulant(){

  hpt=new TH1D("pt","ptdist_allevents",120,0,12);

  const int NMultBin=1;  
  double MultCut[NMultBin+1];
  double MultCut[0]=4; double MultCut[NMultBin]=800;
  int nh, ipt, inttemp;
  double Nh;
  Int_t rmult=0;
  long nevt=0; 
  long Nevt[NMultBin]={0};
  //ref;
  double SumMult[NMultBin]={0.};
  double AveMult[NMultBin]={0.};
  //diff;
  double SumMultp[NMultBin][NBin_PT]={0.};
  double AveMultp[NMultBin][NBin_PT]={0.};
  Float_t phi0[20000];
  Float_t pt0[20000];
  Float_t eta0[20000];
  Float_t vz0;
  double minptcut[NBin_PT] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0};
  double maxptcut[NBin_PT] = {0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,10.0}; //10

  ///ebe;
  double w2, w4, w6, w8;
  double AziCorr22, AziCorr24, AziCorr26, AziCorr28;
  const std::complex<double> i(0.,1.);
  std::complex<double> Q2, Q4, Q422, Q6, Q8;
  //diff
  double meanPt[NMultBin][NBin_PT]={0.0};
  double Nhp[NBin_PT];
  double w2p[NBin_PT], w4p[NBin_PT], w6p[NBin_PT], w8p[NBin_PT];
  double AziCorr22p[NBin_PT], AziCorr24p[NBin_PT], AziCorr26p[NBin_PT], AziCorr28p[NBin_PT];
  std::complex<double> q2p[NBin_PT], q4p[NBin_PT], q6p[NBin_PT], q8p[NBin_PT];   
  /// all events;
  double AveAziCorr22[NMultBin]={0.0}; double AveAziCorr24[NMultBin]={0.0}; double AveAziCorr26[NMultBin]={0.0};  double AveAziCorr28[NMultBin]={0.0};
  double Cumulant22[NMultBin]={0.0}; double Cumulant24[NMultBin]={0.0};  double Cumulant26[NMultBin]={0.0}; double Cumulant28[NMultBin]={0.0};
  double V22[NMultBin]={0.0}; double V24[NMultBin]={0.0}; double V26[NMultBin]={0.0}; double V28[NMultBin]={0.0};
  //diff 
  double AveAziCorr22p[NMultBin][NBin_PT]={0.0}; double AveAziCorr24p[NMultBin][NBin_PT]={0.0}; double AveAziCorr26p[NMultBin][NBin_PT]={0.0}; double AveAziCorr28p[NMultBin][NBin_PT]={0.0};
  double Cumulant22p[NMultBin][NBin_PT]={0.0}; double Cumulant24p[NMultBin][NBin_PT]={0.0}; double Cumulant26p[NMultBin][NBin_PT]={0.0}; double Cumulant28p[NMultBin][NBin_PT]={0.0};
  double V22p[NMultBin][NBin_PT]={0.0}; double V24p[NMultBin][NBin_PT]={0.0}; double V26p[NMultBin][NBin_PT]={0.0}; double V28p[NMultBin][NBin_PT]={0.0};
  /// stat. uncertainty related; all events;
  double Sumw2[NMultBin]={0.0}; double Sumw4[NMultBin]={0.0};  double Sumw6[NMultBin]={0.0}; double Sumw8[NMultBin]={0.0};
  double Sumw2w4[NMultBin]={0.0}; double Sumw2sq[NMultBin]={0.0}; double Sumw4sq[NMultBin]={0.0}; 
  double Sumw2A2[NMultBin]={0.0}; double Sumw4A4[NMultBin]={0.0};  double Sumw6A6[NMultBin]={0.0}; double Sumw8A8[NMultBin]={0.0};
  double Sumw2A2sq[NMultBin]={0.0}; double Sumw4A4sq[NMultBin]={0.0};
  double Sumw2w4A2A4[NMultBin]={0.0}; double CovA2A4[NMultBin]={0.0};
  double AveAziCorr22sq[NMultBin]={0.0}; double AveAziCorr24sq[NMultBin]={0.0}; 
  double S2sq[NMultBin]={0.0}; double S4sq[NMultBin]={0.0};
  double SV22sq[NMultBin]={0.0}; double SV24sq[NMultBin]={0.0};
  //diff
  double Sumw2p[NMultBin][NBin_PT]={0.0}; double Sumw4p[NMultBin][NBin_PT]={0.0};  double Sumw6p[NMultBin][NBin_PT]={0.0}; double Sumw8p[NMultBin][NBin_PT]={0.0};
  double Sumw2A2p[NMultBin][NBin_PT]={0.0}; double Sumw4A4p[NMultBin][NBin_PT]={0.0}; double Sumw6A6p[NMultBin][NBin_PT]={0.0}; double Sumw8A8p[NMultBin][NBin_PT]={0.0};
  //diff error
  double Sumw2psq[NMultBin][NBin_PT]={0.0}; double S2psq[NMultBin][NBin_PT]={0.0};
  double AveAziCorr22psq[NMultBin][NBin_PT]={0.0}; double Sumw2pA2psq[NMultBin][NBin_PT]={0.0};
  double Sumw2w2p[NMultBin][NBin_PT]={0.0}; double CovA2A2p[NMultBin][NBin_PT]={0.0}; double Sumw2w2pA2A2p[NMultBin][NBin_PT]={0.0};
  double Sv22psq[NMultBin][NBin_PT]={0.0}; 
  double Sumw4psq[NMultBin][NBin_PT]={0.0}; double S4psq[NMultBin][NBin_PT]={0.0}; 
  double AveAziCorr24psq[NMultBin][NBin_PT]={0.0}; double Sumw4pA4psq[NMultBin][NBin_PT]={0.0};
  double Sumw2w4p[NMultBin][NBin_PT]={0.0}; double CovA2A4p[NMultBin][NBin_PT]={0.0}; double Sumw2w4pA2A4p[NMultBin][NBin_PT]={0.0};  
  double Sumw4w2p[NMultBin][NBin_PT]={0.0}; double CovA4A2p[NMultBin][NBin_PT]={0.0}; double Sumw4w2pA4A2p[NMultBin][NBin_PT]={0.0};
  double Sumw4w4p[NMultBin][NBin_PT]={0.0}; double CovA4A4p[NMultBin][NBin_PT]={0.0}; double Sumw4w4pA4A4p[NMultBin][NBin_PT]={0.0};
  double Sumw2pw4p[NMultBin][NBin_PT]={0.0}; double CovA2pA4p[NMultBin][NBin_PT]={0.0}; double Sumw2pw4pA2pA4p[NMultBin][NBin_PT]={0.0};
  double Sv24psq[NMultBin][NBin_PT]={0.0};

     char *file1 = "/home/xuq7/CMSSW_5_3_8_HI/src/UserCode/flow/proSTEGvn/vndata_10k.root";
     TFile *f1 = new TFile(file1);
     TTree *t1 = (TTree*)f1->Get("tree");
     t1->SetBranchAddress("npg", &rmult);
     t1->SetBranchAddress("ptg", &pt0);
     t1->SetBranchAddress("etag", &eta0);
     t1->SetBranchAddress("phig", &phi0);
     nevt=t1->GetEntries();

     cout<<"Nevts = "<<nevt<<endl;
     for(long ne=0; ne<nevt; ne++)
     //for(long ne=0; ne<2000; ne++)
     {
       if(ne%100==0)  cout<<"Have run "<<ne<<" of all "<<nevt<<" events;"<<endl;
       t1->GetEntry(ne);
 
       //event cut here
       //if(rmult<10) continue;
     for(int nm=0;nm<NMultBin;nm++)
      if(rmult>=MultCut[nm]&&rmult<MultCut[nm+1]){
       Nh=0; 
       w2=0.0; w4=0.0; w6=0.0; w8=0.0;
       AziCorr22=0.0; AziCorr24=0.0; AziCorr26=0.0; AziCorr28=0.0;
       Q2=0.0; Q4=0.0; Q422=0.0; Q6=0.0; Q8=0.0;
       for(ipt=0;ipt<NBin_PT;ipt++){
         Nhp[ipt]=0;        
         w2p[ipt]=0; w4p[ipt]=0; w6p[ipt]=0; w8p[ipt]=0;
         q2p[ipt]=0; q4p[ipt]=0; q6p[ipt]=0; q8p[ipt]=0; 
         AziCorr22p[ipt]=0;  AziCorr24p[ipt]=0;  AziCorr26p[ipt]=0;  AziCorr28p[ipt]=0;
       }
///ref;
       for(nh=0;nh<rmult;nh++){
           if(eta0[nh]<-2.40||eta0[nh]>2.40) continue;           
           if(pt0[nh]<0.2||pt0[nh]>10.0) continue;
	   hpt->Fill(pt0[nh]);
	   Nh++;
	 } // end of loop over particles
       if(Nh<4) continue;  // event cut; 
       Nevt[nm]++;
//diff;
     for(nh=0;nh<rmult;nh++){
       if(eta0[nh]<-2.40||eta0[nh]>2.40) continue;
       if(pt0[nh]<0.2||pt0[nh]>10.0) continue;
       Q2+=(cos(2.0*phi0[nh])+i*sin(2.0*phi0[nh]));
       Q4+=(cos(4.0*phi0[nh])+i*sin(4.0*phi0[nh]));
       //Q6+=(cos(6.0*phi0[nh])+i*sin(6.0*phi0[nh]));
       //Q8+=(cos(8.0*phi0[nh])+i*sin(8.0*phi0[nh]));
       for(ipt=0;ipt<NBin_PT;ipt++) {
          if((pt0[nh]>minptcut[ipt]&&pt0[nh]<=maxptcut[ipt])){
             inttemp=ipt;
          }
       }
       ipt=inttemp;
       Nhp[ipt]++; meanPt[nm][ipt]+=pt0[nh];
       q2p[ipt]+=(cos(2.0*phi0[nh])+i*sin(2.0*phi0[nh]));
       q4p[ipt]+=(cos(4.0*phi0[nh])+i*sin(4.0*phi0[nh]));
       //q6p[ipt]+=(cos(6.0*phi0[nh])+i*sin(6.0*phi0[nh]));
       //q8p[ipt]+=(cos(8.0*phi0[nh])+i*sin(8.0*phi0[nh]));
     } // end of loop over particles  
///ref;
       Q422=Q4*conj(Q2)*conj(Q2);
       SumMult[nm]+=Nh;
         w2=Nh*(Nh-1); 
         w4=Nh*(Nh-1)*(Nh-2)*(Nh-3);
         //w6=Nh*(Nh-1)*(Nh-2)*(Nh-3)*(Nh-4)*(Nh-5); 
         //w8=Nh*(Nh-1)*(Nh-2)*(Nh-3)*(Nh-4)*(Nh-5)*(Nh-6)*(Nh-7);
         AziCorr22=(abs(Q2)*abs(Q2)-Nh)/w2;
         double s24_1=pow(abs(Q2),4)+abs(Q4)*abs(Q4)-2.0*Q422.real();
         double s24_2=2.0*(2.0*(Nh-2)*abs(Q2)*abs(Q2)-Nh*(Nh-3));
         //AziCorr24=0; AziCorr26=0; AziCorr28=0;
         AziCorr24=(s24_1-s24_2)/w4;
         //AziCorr24=(pow(abs(Q2),4)+abs(Q4)*abs(Q4)-2.0*Q422.real())/w4 - 2.0*(2.0*(Nh-2)*abs(Q2)*abs(Q2)-Nh*(Nh-3))/w4;
//         AziCorr24=AziCorr24;   AziCorr24=AziCorr24;
//         cout<<"ne="<<ne<<", AziCorr22="<<AziCorr22<<", AziCorr24="<<AziCorr24<<endl;
/*         double s26_1=pow(abs(Q2),6)+9.0*pow(abs(Q4),2)*pow(abs(Q2),2)-6.0*real(Q4*Q2*conj(Q2)*conj(Q2)*conj(Q2));
         double s26_2=4.0*real(Q6*conj(Q2)*conj(Q2)*conj(Q2))-12.0*real(Q6*conj(Q4)*conj(Q2));
         double s26_3=18.0*(Nh-4)*real(Q4*conj(Q2)*conj(Q2))+4.0*pow(abs(Q6),2);
         double s26_4=9.0*(Nh-4)*(pow(abs(Q2),4)+pow(abs(Q4),2))+18.0*(Nh-2)*(Nh-5)*pow(abs(Q2),2);
         double s26_5=6.0*Nh*(Nh-4)*(Nh-5);
         AziCorr26=(s26_1+s26_2+s26_3-s26_4-s26_5)/w6;
         double s28_1=pow(abs(Q2),8)-12.0*Q4*Q2*Q2*conj(Q2)*conj(Q2)*conj(Q2)*conj(Q2);
         double s28_2=6.0*Q4*Q4*conj(Q2)*conj(Q2)*conj(Q2)*conj(Q2)+16.0*Q6*Q2*conj(Q2)*conj(Q2)*conj(Q2)*conj(Q2);
         double s28_3=96.0*Q6*Q2*conj(Q4)*conj(Q2)*conj(Q2)-12.0*Q8*conj(Q2)*conj(Q2)*conj(Q2)*conj(Q2);
         double s28_4=36.0*Q4*Q4*conj(Q4)*conj(Q2)*conj(Q2)+96.0*(Nh-6)*Q4*Q2*conj(Q2)*conj(Q2)*conj(Q2);
         double s28_5=72.0*Q8*conj(Q4)*conj(Q2)*conj(Q2)+48.0*Q6*Q2*conj(Q4)*conj(Q4);
         double s28_6=64.0*(Nh-6)*Q6*conj(Q2)*conj(Q2)*conj(Q2)+192.0*(Nh-6)*Q6*conj(Q4)*conj(Q2);
         double s28_7=96.0*Q8*conj(Q6)*conj(Q2)-36.0*Q8*conj(Q4)*conj(Q4);
         double s28_8=144.0*(Nh-7)*(Nh-4)*Q4*conj(Q2)*conj(Q2)+36.0*pow(abs(Q8),2)+64.0*pow(abs(Q6),2)*pow(abs(Q2),2);
         double s28_9=64.0*(Nh-6)*pow(abs(Q6),2)+9.0*pow(abs(Q4),4)+36.0*pow(abs(Q2),4)*pow(abs(Q4),2)-144.0*(Nh-6)*pow(abs(Q4),2)*pow(abs(Q2),2);
         double s28_10=72.0*(Nh-7)*(Nh-4)*(pow(abs(Q4),2)+pow(abs(Q2),4))-16.0*(Nh-6)*pow(abs(Q2),6);
         double s28_11=96.0*(Nh-7)*(Nh-6)*(Nh-2)*pow(abs(Q2),2);
         double s28_12=24.0*Nh*(Nh-7)*(Nh-6)*(Nh-5);
         AziCorr28=(s28_1+s28_2-s28_3-s28_4+s28_5-s28_6-s28_7-s28_8-s28_9+s28_10-s28_11+s28_12)/w8; */ AziCorr26=AziCorr22; AziCorr28=AziCorr24;
         //cout<<"ne="<<ne<<", AziCorr22="<<AziCorr22<<", AziCorr24="<<AziCorr24<<endl;

         //cout<<"ne="<<ne<<", AziCorr26="<<AziCorr26<<", AziCorr28="<<AziCorr28<<endl;
/*         AziCorr26=(pow(abs(Q2),6)+9.0*pow(abs(Q4),2)*pow(abs(Q2),2)-6.0*real(Q4*Q2*conj(Q2)*conj(Q2)*conj(Q2)) + 
                    4.0*real(Q6*conj(Q2)*conj(Q2)*conj(Q2))-12.0*real(Q6*conj(Q4)*conj(Q2)) +
                    18.0*(Nh-4)*real(Q4*conj(Q2)*conj(Q2))+4.0*pow(abs(Q6),2) -
                    9.0*(Nh-4)*(pow(abs(Q2),4)+pow(abs(Q4),2))+18.0*(Nh-2)*(Nh-5)*pow(abs(Q2),2) -
                    6.0*Nh*(Nh-4)*(Nh-5))/w6;
         AziCorr28=(pow(abs(Q2),8)-12.0*Q4*Q2*Q2*conj(Q2)*conj(Q2)*conj(Q2)*conj(Q2) +
                   6.0*Q4*Q4*conj(Q2)*conj(Q2)*conj(Q2)*conj(Q2)+16.0*Q6*Q2*conj(Q2)*conj(Q2)*conj(Q2)*conj(Q2) - 
                   96.0*Q6*Q2*conj(Q4)*conj(Q2)*conj(Q2)-12.0*Q8*conj(Q2)*conj(Q2)*conj(Q2)*conj(Q2) - 
                   36.0*Q4*Q4*conj(Q4)*conj(Q2)*conj(Q2)+96.0*(Nh-6)*Q4*Q2*conj(Q2)*conj(Q2)*conj(Q2) + 
                   72.0*Q8*conj(Q4)*conj(Q2)*conj(Q2)+48.0*Q6*Q2*conj(Q4)*conj(Q4) - 
                   64.0*(Nh-6)*Q6*conj(Q2)*conj(Q2)*conj(Q2)+192.0*(Nh-6)*Q6*conj(Q4)*conj(Q2) - 
                   96.0*Q8*conj(Q6)*conj(Q2)-36.0*Q8*conj(Q4)*conj(Q4) - 
                   144.0*(Nh-7)*(Nh-4)*Q4*conj(Q2)*conj(Q2)+36.0*pow(abs(Q8),2)+64.0*pow(abs(Q6),2)*pow(abs(Q2),2) - 
                   64.0*(Nh-6)*pow(abs(Q6),2)+9.0*pow(abs(Q4),4)+36.0*pow(abs(Q2),4)*pow(abs(Q4),2)-144.0*(Nh-6)*pow(abs(Q4),2)*pow(abs(Q2),2) + 
                   72.0*(Nh-7)*(Nh-4)*(pow(abs(Q4),2)+pow(abs(Q2),4))-16.0*(Nh-6)*pow(abs(Q2),6) - 
                   96.0*(Nh-7)*(Nh-6)*(Nh-2)*pow(abs(Q2),2) +
                   24.0*Nh*(Nh-7)*(Nh-6)*(Nh-5))/w8;  //  AziCorr26=AziCorr22; AziCorr28=AziCorr24; */
       Sumw2[nm]+=w2;  Sumw4[nm]+=w4; Sumw2w4[nm]+=w2*w4;
       //Sumw6[nm]+=w6;  Sumw8[nm]+=w8; 
       Sumw2sq[nm]+=w2*w2; Sumw4sq[nm]+=w4*w4;
       Sumw2A2[nm]+=w2*AziCorr22; Sumw4A4[nm]+=w4*AziCorr24;
       //Sumw6A6[nm]+=w6*AziCorr26; Sumw8A8[nm]+=w8*AziCorr28;
       Sumw2A2sq[nm]+=w2*AziCorr22*AziCorr22; Sumw4A4sq[nm]+=w4*AziCorr24*AziCorr24;
       Sumw2w4A2A4[nm]+=w2*w4*AziCorr22*AziCorr24;
       //cout<<"ne="<<ne<<", Nh="<<Nh<<", w2="<<w2<<", w4="<<w4<<", w6="<<w6<<", w8="<<w8<<endl;
       //cout<<"Sumw2A2="<<Sumw2A2[nm]<<", Sumw4A4="<<Sumw4A4[nm]<<", Sumw6A6="<<Sumw6A6[nm]<<", Sumw8A8="<<Sumw8A8[nm]<<endl;
//diff;
       for(ipt=0;ipt<NBin_PT;ipt++) {
         SumMultp[nm][ipt]+=((float) Nhp[ipt]);
         w2p[ipt]=Nhp[ipt]*(Nh-1);  w4p[ipt]=Nhp[ipt]*(Nh-1)*(Nh-2)*(Nh-3);
         if(w4p[ipt]>0){
           AziCorr22p[ipt]=(real(q2p[ipt]*conj(Q2))-Nhp[ipt])/w2p[ipt];
           AziCorr24p[ipt]=(real(q2p[ipt]*Q2*conj(Q2)*conj(Q2))-real(q4p[ipt]*conj(Q2)*conj(Q2))-real(q2p[ipt]*Q2*conj(Q4))+real(q4p[ipt]*conj(Q4)))/w4p[ipt] - 2.0*(Nhp[ipt]*abs(Q2)*abs(Q2)+(Nh-3)*real(q2p[ipt]*conj(Q2))-real(Q2*conj(q2p[ipt])))/w4p[ipt] + 2.0/(Nh-1)/(Nh-2);
         Sumw2p[nm][ipt]+=w2p[ipt];  Sumw4p[nm][ipt]+=w4p[ipt]; 
         Sumw2A2p[nm][ipt]+=w2p[ipt]*AziCorr22p[ipt]; Sumw4A4p[nm][ipt]+=w4p[ipt]*AziCorr24p[ipt];
         //diff error, 2
         Sumw2psq[nm][ipt]+=w2p[ipt]*w2p[ipt]; Sumw2pA2psq[nm][ipt]+=w2p[ipt]*AziCorr22p[ipt]*AziCorr22p[ipt];
         Sumw2w2p[nm][ipt]+=w2*w2p[ipt]; Sumw2w2pA2A2p[nm][ipt]+=w2*w2p[ipt]*AziCorr22*AziCorr22p[ipt];
         //diff error, 4
         Sumw4psq[nm][ipt]+=w4p[ipt]*w4p[ipt]; Sumw4pA4psq[nm][ipt]+=w4p[ipt]*AziCorr24p[ipt]*AziCorr24p[ipt];
         Sumw2w4p[nm][ipt]+=w2*w4p[ipt]; Sumw2w4pA2A4p[nm][ipt]+=w2*w4p[ipt]*AziCorr22*AziCorr24p[ipt];
         Sumw4w2p[nm][ipt]+=w4*w2p[ipt]; Sumw4w2pA4A2p[nm][ipt]+=w4*w2p[ipt]*AziCorr24*AziCorr22p[ipt];
         Sumw4w4p[nm][ipt]+=w4*w4p[ipt]; Sumw4w4pA4A4p[nm][ipt]+=w4*w4p[ipt]*AziCorr24*AziCorr24p[ipt];
         Sumw2pw4p[nm][ipt]+=w2p[ipt]*w4p[ipt]; Sumw2pw4pA2pA4p[nm][ipt]+=w2p[ipt]*w4p[ipt]*AziCorr22p[ipt]*AziCorr24p[ipt];
         }
       }
       
      }  // end of loop over mult cut;
     } // end of loop over events;

    for(int nm=0;nm<NMultBin;nm++){
     AveMult[nm]=SumMult[nm]*1.0/Nevt[nm];
     cout<<"Totevt: "<<Nevt[nm] <<  " Totmult: " <<SumMult[nm]<<",  <Mult> = "<<AveMult[nm]<<endl;
///ref;
     AveAziCorr22[nm]=Sumw2A2[nm]/Sumw2[nm];
     AveAziCorr24[nm]=Sumw4A4[nm]/Sumw4[nm];
     //AveAziCorr26[nm]=Sumw6A6[nm]/Sumw6[nm];
     //AveAziCorr28[nm]=Sumw8A8[nm]/Sumw8[nm];
     AveAziCorr22sq[nm]=Sumw2A2sq[nm]/Sumw2[nm];
     AveAziCorr24sq[nm]=Sumw4A4sq[nm]/Sumw4[nm];
     Cumulant22[nm]=AveAziCorr22[nm];    
     Cumulant24[nm]=AveAziCorr24[nm]-2.0*AveAziCorr22[nm]*AveAziCorr22[nm];
     //Cumulant26[nm]=AveAziCorr26[nm]-9.0*AveAziCorr22[nm]*AveAziCorr24[nm]+12.0*pow(AveAziCorr22[nm],3.0);
     //Cumulant28[nm]=AveAziCorr28[nm]-16.0*AveAziCorr26[nm]*AveAziCorr22[nm]-18.0*pow(AveAziCorr24[nm],2.0) + 
     //               144.0*AveAziCorr24[nm]*pow(AveAziCorr22[nm],2.0)-144.0*pow(AveAziCorr22[nm],4.0); 
     cout<<"Cumulant22="<<Cumulant22[nm]<<", Cumulant24="<<Cumulant24[nm]<<endl;
     //cout<<"Cumulant26="<<Cumulant26[nm]<<", Cumulant28="<<Cumulant28[nm]<<endl;
     if(Cumulant22[nm]>0.0)  V22[nm]=sqrt(Cumulant22[nm]);  //else V22[nm]=0.0;
     if(Cumulant24[nm]<0.0)  V24[nm]=pow(-Cumulant24[nm],0.25);  else V24[nm]=0.0; 
     //if(Cumulant26[nm]>0.0)  V26[nm]=pow(Cumulant26[nm]/4.0,1.0/6);  else V26[nm]=0.0;
     //if(Cumulant28[nm]<0.0)  V28[nm]=pow(-Cumulant28[nm]/33.0,1.0/8);  else V28[nm]=0.0; 
//diff;
    for(ipt=0;ipt<NBin_PT;ipt++) {
     AveMultp[nm][ipt]=SumMultp[nm][ipt]*1.0/Nevt[nm];
     meanPt[nm][ipt]=meanPt[nm][ipt]/SumMultp[nm][ipt];
     if(Sumw4p[nm][ipt]>0){
      AveAziCorr22p[nm][ipt]=Sumw2A2p[nm][ipt]/Sumw2p[nm][ipt];
      AveAziCorr24p[nm][ipt]=Sumw4A4p[nm][ipt]/Sumw4p[nm][ipt];
      //cout<<"AveAziCorr22p= "<<AveAziCorr22p[nm][ipt]<<", Sumw2p="<<Sumw2p[nm][ipt]<<endl;
      Cumulant22p[nm][ipt]=AveAziCorr22p[nm][ipt];
      Cumulant24p[nm][ipt]=AveAziCorr24p[nm][ipt]-2.0*AveAziCorr22p[nm][ipt]*AveAziCorr22[nm];
      //cout<<"Cumulant22p="<<Cumulant22p[nm][ipt]<<", Cumulant24p="<<Cumulant24p[nm][ipt]<<endl;
      if(V22[nm]>0) V22p[nm][ipt]=Cumulant22p[nm][ipt]/V22[nm]; else V22p[nm][ipt]=0.0;
      if(V24[nm]>0) V24p[nm][ipt]=-Cumulant24p[nm][ipt]/pow(V24[nm],3.0); else V24p[nm][ipt]=0.0;      
     }
     else {V22p[nm][ipt]=-1.0; V24p[nm][ipt]=-1.0; cout<<"m<1; so few POIs in all events;"<<endl; }
     //diff error
     AveAziCorr22psq[nm][ipt]=Sumw2pA2psq[nm][ipt]/Sumw2p[nm][ipt];
     AveAziCorr24psq[nm][ipt]=Sumw4pA4psq[nm][ipt]/Sumw4p[nm][ipt];
     //cout<<"AveAziCorr22psq="<<AveAziCorr22psq[nm][ipt]<<", AveAziCorr24psq="<<AveAziCorr24psq[nm][ipt]<<endl;
    }
// uncertainty calculation;
     S2sq[nm]=(AveAziCorr22sq[nm]-AveAziCorr22[nm]*AveAziCorr22[nm])*1.0/(1-Sumw2sq[nm]/Sumw2[nm]/Sumw2[nm]);
     S4sq[nm]=(AveAziCorr24sq[nm]-AveAziCorr24[nm]*AveAziCorr24[nm])*1.0/(1-Sumw4sq[nm]/Sumw4[nm]/Sumw4[nm]);
     CovA2A4[nm]=(Sumw2w4A2A4[nm]/Sumw2w4[nm]-AveAziCorr22[nm]*AveAziCorr24[nm])/(1-Sumw2w4[nm]/Sumw2[nm]/Sumw4[nm]);

     SV22sq[nm]=1.0/4.0/AveAziCorr22[nm]*Sumw2sq[nm]/Sumw2[nm]/Sumw2[nm]*S2sq[nm];
     if(Cumulant24[nm]<0.0) 
       SV24sq[nm]=1.0/pow(2.0*AveAziCorr22[nm]*AveAziCorr22[nm]-AveAziCorr24[nm],3.0/2)*(AveAziCorr22[nm]*AveAziCorr22[nm]*Sumw2sq[nm]/Sumw2[nm]/Sumw2[nm]*S2sq[nm]+1.0/16*Sumw4sq[nm]/Sumw4[nm]/Sumw4[nm]*S4sq[nm]-1.0/2*AveAziCorr22[nm]*Sumw2w4[nm]/Sumw2[nm]/Sumw4[nm]*CovA2A4[nm]);
     else SV24sq[nm]=0.0;
//diff;
    for(ipt=0;ipt<NBin_PT;ipt++) {
     if(Sumw4p[nm][ipt]>0){
      S2psq[nm][ipt]=(AveAziCorr22psq[nm][ipt]-AveAziCorr22p[nm][ipt]*AveAziCorr22p[nm][ipt])*1.0/(1-Sumw2psq[nm][ipt]/Sumw2p[nm][ipt]/Sumw2p[nm][ipt]);
      CovA2A2p[nm][ipt]=(Sumw2w2pA2A2p[nm][ipt]/Sumw2w2p[nm][ipt]-AveAziCorr22[nm]*AveAziCorr22p[nm][ipt])/(1-Sumw2w2p[nm][ipt]/Sumw2[nm]/Sumw2p[nm][ipt]);
      S4psq[nm][ipt]=(AveAziCorr24psq[nm][ipt]-AveAziCorr24p[nm][ipt]*AveAziCorr24p[nm][ipt])*1.0/(1-Sumw4psq[nm][ipt]/Sumw4p[nm][ipt]/Sumw4p[nm][ipt]);
      CovA2A4p[nm][ipt]=(Sumw2w4pA2A4p[nm][ipt]/Sumw2w4p[nm][ipt]-AveAziCorr22[nm]*AveAziCorr24p[nm][ipt])/(1-Sumw2w4p[nm][ipt]/Sumw2[nm]/Sumw4p[nm][ipt]);
      CovA4A2p[nm][ipt]=(Sumw4w2pA4A2p[nm][ipt]/Sumw4w2p[nm][ipt]-AveAziCorr24[nm]*AveAziCorr22p[nm][ipt])/(1-Sumw4w2p[nm][ipt]/Sumw4[nm]/Sumw2p[nm][ipt]);
      CovA4A4p[nm][ipt]=(Sumw4w4pA4A4p[nm][ipt]/Sumw4w4p[nm][ipt]-AveAziCorr24[nm]*AveAziCorr24p[nm][ipt])/(1-Sumw4w4p[nm][ipt]/Sumw4[nm]/Sumw4p[nm][ipt]);
      CovA2pA4p[nm][ipt]=(Sumw2pw4pA2pA4p[nm][ipt]/Sumw2pw4p[nm][ipt]-AveAziCorr22p[nm][ipt]*AveAziCorr24p[nm][ipt])/(1-Sumw2pw4p[nm][ipt]/Sumw2p[nm][ipt]/Sumw4p[nm][ipt]);

      Sv22psq[nm][ipt]=1.0/4.0/pow(AveAziCorr22[nm],3)*(pow(AveAziCorr22p[nm][ipt],2)*Sumw2sq[nm]/Sumw2[nm]/Sumw2[nm]*S2sq[nm]+4.0*pow(AveAziCorr22[nm],2)*Sumw2psq[nm][ipt]/Sumw2p[nm][ipt]/Sumw2p[nm][ipt]*S2psq[nm][ipt]-4.0*AveAziCorr22[nm]*AveAziCorr22p[nm][ipt]*Sumw2w2p[nm][ipt]/Sumw2[nm]/Sumw2p[nm][ipt]*CovA2A2p[nm][ipt]);
      if(V24[nm]>0)
        Sv24psq[nm][ipt]=1.0/pow(2*AveAziCorr22[nm]*AveAziCorr22[nm]-AveAziCorr24[nm],7.0/2)*(pow(2*AveAziCorr22[nm]*AveAziCorr22[nm]*AveAziCorr22p[nm][ipt]-3*AveAziCorr22[nm]*AveAziCorr24p[nm][ipt]+2*AveAziCorr24[nm]*AveAziCorr22p[nm][ipt],2)*Sumw2sq[nm]/Sumw2[nm]/Sumw2[nm]*S2sq[nm] + 9.0/16*pow(2*AveAziCorr22[nm]*AveAziCorr22p[nm][ipt]-AveAziCorr24p[nm][ipt],2)*Sumw4sq[nm]/Sumw4[nm]/Sumw4[nm]*S4sq[nm] + 4.0*AveAziCorr22[nm]*AveAziCorr22[nm]*pow(2*AveAziCorr22[nm]*AveAziCorr22[nm]-AveAziCorr24[nm],2)*Sumw2psq[nm][ipt]/Sumw2p[nm][ipt]/Sumw2p[nm][ipt]*S2psq[nm][ipt] + pow(2*AveAziCorr22[nm]*AveAziCorr22[nm]-AveAziCorr24[nm],2)*Sumw4psq[nm][ipt]/Sumw4p[nm][ipt]/Sumw4p[nm][ipt]*S4psq[nm][ipt] - 3.0/2*(2*AveAziCorr22[nm]*AveAziCorr22p[nm][ipt]-AveAziCorr24p[nm][ipt])*(2*AveAziCorr22[nm]*AveAziCorr22[nm]*AveAziCorr22p[nm][ipt]-3*AveAziCorr22[nm]*AveAziCorr24p[nm][ipt]+2*AveAziCorr24[nm]*AveAziCorr22p[nm][ipt])*Sumw2w4[nm]/Sumw2[nm]/Sumw4[nm]*CovA2A4[nm] - 4.0*AveAziCorr22[nm]*(2*AveAziCorr22[nm]*AveAziCorr22[nm]-AveAziCorr24[nm])*(2*AveAziCorr22[nm]*AveAziCorr22[nm]*AveAziCorr22p[nm][ipt]-3*AveAziCorr22[nm]*AveAziCorr24p[nm][ipt]+2*AveAziCorr24[nm]*AveAziCorr22p[nm][ipt])*Sumw2w2p[nm][ipt]/Sumw2[nm]/Sumw2p[nm][ipt]*CovA2A2p[nm][ipt] + 2.0*(2*AveAziCorr22[nm]*AveAziCorr22[nm]-AveAziCorr24[nm])*(2*AveAziCorr22[nm]*AveAziCorr22[nm]*AveAziCorr22p[nm][ipt]-3*AveAziCorr22[nm]*AveAziCorr24p[nm][ipt]+2*AveAziCorr24[nm]*AveAziCorr22p[nm][ipt])*Sumw2w4p[nm][ipt]/Sumw2[nm]/Sumw4p[nm][ipt]*CovA2A4p[nm][ipt] + 3.0*AveAziCorr22[nm]*(2*AveAziCorr22[nm]*AveAziCorr22[nm]-AveAziCorr24[nm])*(2*AveAziCorr22[nm]*AveAziCorr22p[nm][ipt]-AveAziCorr24p[nm][ipt])*Sumw4w2p[nm][ipt]/Sumw4[nm]/Sumw2p[nm][ipt]*CovA4A2p[nm][ipt] -3.0/2*(2*AveAziCorr22[nm]*AveAziCorr22[nm]-AveAziCorr24[nm])*(2*AveAziCorr22[nm]*AveAziCorr22p[nm][ipt]-AveAziCorr24p[nm][ipt])*Sumw4w4p[nm][ipt]/Sumw4[nm]/Sumw4p[nm][ipt]*CovA4A4p[nm][ipt] - 4.0*AveAziCorr22[nm]*pow(2*AveAziCorr22[nm]*AveAziCorr22[nm]-AveAziCorr24[nm],2)*Sumw2pw4p[nm][ipt]/Sumw2p[nm][ipt]/Sumw4p[nm][ipt]*CovA2pA4p[nm][ipt]);
      else Sv22psq[nm][ipt]=0.0;
     //cout<<"S2psq="<<S2psq[nm][ipt]<<", CovA2A2p="<<CovA2A2p[nm][ipt]<<", S4psq="<<S4psq[nm][ipt]<<", CovA2pA4p="<<CovA2pA4p[nm][ipt]<<endl;
     //cout<<"Sv22psq="<<Sv22psq[nm][ipt]<<", Sv24psq[nm][ipt]="<<Sv24psq[nm][ipt]<<endl;
     }
     else { Sv24psq[nm][ipt]=0.0; cout<<"In Error calculation; m<1; so few POIs in all events;"<<endl; }
    }

     cout<<"V22= "<<V22[nm]<<" +/- "<<sqrt(SV22sq[nm])<<endl;
     cout<<"V24= "<<V24[nm]<<" +/- "<<sqrt(SV24sq[nm])<<endl;
     //cout<<"V26= "<<V26[nm]<<" +/- "<<sqrt(SV24sq[nm])<<endl;
     //cout<<"V28= "<<V28[nm]<<" +/- "<<sqrt(SV24sq[nm])<<endl;
     cout<<"differential flow:"<<endl;
     for(ipt=0;ipt<NBin_PT;ipt++) 
        cout<<"v22= "<<V22p[nm][ipt]<<" +/- "<<sqrt(Sv22psq[nm][ipt])<<", v24= "<<V24p[nm][ipt]<<" +/- "<<sqrt(Sv24psq[nm][ipt])<<endl;

    } // end of loop over mult cut;
  ofstream output;
  ofstream outv2;
  ofstream outcumulant;
  output.open("output.txt");
    for(int nm=0;nm<NMultBin;nm++) for(ipt=0;ipt<NBin_PT;ipt++) output<<"<Multp>= "<<AveMultp[nm][ipt]<<", meanPt="<<meanPt[nm][ipt]<<", v22= "<<V22p[nm][ipt]<<", v24= "<<V24p[nm][ipt]<<endl;
    for(int nm=0;nm<NMultBin;nm++)  output<<"Totevt: "<<Nevt[nm] <<  " Totmult: " <<SumMult[nm]<<",  <Mult> = "<<AveMult[nm]<<endl;
    for(int nm=0;nm<NMultBin;nm++)  output<<"<Mult>= "<<AveMult[nm]<<",  V22= "<<V22[nm]<<" +/- "<<sqrt(SV22sq[nm])<<endl;
    for(int nm=0;nm<NMultBin;nm++)  output<<"<Mult>= "<<AveMult[nm]<<",  V24= "<<V24[nm]<<" +/- "<<sqrt(SV24sq[nm])<<endl;
    //for(int nm=0;nm<NMultBin;nm++)  output<<"<Mult>= "<<AveMult[nm]<<",  V26= "<<V26[nm]<<" +/- "<<sqrt(1.0/SumMult[nm]/2)<<endl;
    //for(int nm=0;nm<NMultBin;nm++)  output<<"<Mult>= "<<AveMult[nm]<<",  V28= "<<V28[nm]<<" +/- "<<sqrt(1.0/SumMult[nm]/2)<<endl;
  output.close();
  outv2.open("V2.txt");
    for(int nm=0;nm<NMultBin;nm++) for(ipt=0;ipt<NBin_PT;ipt++)  outv2<<meanPt[nm][ipt]<<"   "<<V22p[nm][ipt]<<"   "<<sqrt(Sv22psq[nm][ipt])<<"   "<<V24p[nm][ipt]<<"   "<<sqrt(Sv24psq[nm][ipt])<<endl;

    outv2<<endl<<endl;
    for(int nm=0;nm<NMultBin;nm++)  outv2<<AveMult[nm]<<"   "<<V22[nm]<<"   "<<sqrt(SV22sq[nm])<<"   "<<V24[nm]<<"   "<<sqrt(SV24sq[nm])<<endl;
    //for(int nm=0;nm<NMultBin;nm++)  outv2<<AveMult[nm]<<"   "<<V26[nm]<<"   "<<sqrt(1.0/SumMult[nm]/2)<<"   "<<V28[nm]<<"   "<<sqrt(1.0/SumMult[nm]/2)<<endl;
    //for(int nm=0;nm<NMultBin;nm++)  outv2<<AveMult[nm]<<"   "<<V24[nm]<<"   "<<sqrt(SV24sq[nm])<<endl;
  outv2.close();
  outcumulant.open("Cumulant.txt");
    for(int nm=0;nm<NMultBin;nm++) for(ipt=0;ipt<NBin_PT;ipt++) outcumulant<<AveMultp[nm][ipt]<<"   "<<Cumulant22p[nm][ipt]<<"   "<<Cumulant24p[nm][ipt]<<endl;
    for(int nm=0;nm<NMultBin;nm++)  outcumulant<<AveMult[nm]<<"   "<<Cumulant22[nm]<<"   "<<Cumulant24[nm]<<"   "<<Cumulant26[nm]<<"   "<<Cumulant28[nm]<<endl;
    //for(int nm=0;nm<NMultBin;nm++)  outcumulant<<AveMult[nm]<<"   "<<V24[nm]<<"   "<<sqrt(SV24sq[nm])<<endl;
  outcumulant.close();

  TFile *f = new TFile("dist_variables.root","recreate");
  f->cd();
 
  hpt->Write();

  f->Close();

 
}
 

