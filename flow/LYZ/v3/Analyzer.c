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
#define rootJ0 2.4048256
#define J1rootJ0 0.519147 // J1(rootJ0) 
#define J2rootJ0 0.431755 // J2(rootJ0) 
#define kr 200
#define kt 5
//#define Vmax 0.13*500
//#define Vstep 0.0006*500
#define Vmax 0.075*500
#define Vstep 0.00034*500
#define Vmax3 0.075*500
#define Vstep3 0.00034*500
#define Vmax5 0.035*500
#define Vstep5 0.00015*500
#define NB_SLICE 1
#define NB_SLICE_PT 13
#define sqr(x) ((x)*(x))
#define Pi 3.14159265359 
//#define PI TMath::Pi()

 
using namespace std; 

 
double BesselJ0(double x);
double BesselJ0(double x)
{
  double tmp=1., xn=1.;
  long n, Nmax;

  Nmax=int(2.*x)+4;
  for (n=1;n<Nmax;n++)
  {
    xn*=(-sqr(x/2./((double) n)));
    tmp+=xn;
  }
  return tmp;
}

Analyzer(){

  ofstream outR;
  ofstream outG;
  ofstream output;
  ofstream plotVint;
  ofstream plotvdiff;
  ofstream outv2;
  ofstream outv3;
  ofstream outv4;
  ofstream outv5;
  ofstream outv6;
  ofstream outv8;

    ofstream outevt;
    ofstream outevents;

  cent=new TH1D("cent_0-100%","Multdist",2000,0,10000);
  cent1=new TH1D("cent_0-10%","",2000,0,10000);
  cent2=new TH1D("cent_10-20%","",2000,0,10000);
  cent3=new TH1D("cent_20-30%","",2000,0,10000);
  cent4=new TH1D("cent_30-40%","",2000,0,10000);
  cent5=new TH1D("cent_40-50%","",2000,0,10000);
  cent6=new TH1D("cent_50-60%","",2000,0,10000);
  cent7=new TH1D("cent_60-70%","",2000,0,10000);
  cent8=new TH1D("cent_70-80%","",2000,0,10000);
  cent9=new TH1D("cent_80-90%","",2000,0,10000);
  cent10=new TH1D("cent_90-100%","",2000,0,10000);
  cent_used=new TH1D("cent_used","",2000,0,10000);
  npart_used=new TH1D("npart_used","npartdist",500,0,500);
  npartsig_used=new TH1D("npartsig_used","npartsigdist",500,0,500);
  ncoll_used=new TH1D("ncoll_used","ncolldist",1000,0,2000);
  ncollsig_used=new TH1D("ncollsig_used","ncollsigdist",1000,0,2000);

  hpt=new TH1D("pt","ptdist_allevents",200,0,14);
  hphi=new TH1D("phi","phidist_allevents",200,-4,4);
  heta=new TH1D("eta","etadist_allevents",200,-4,4);

  hptt=new TH1D("ptt","ptdistt_allevents",200,0,14);
  hphit=new TH1D("phit","phidistt_allevents",200,-1,7);
  hetat=new TH1D("etat","etadistt_allevents",200,-4,4);


  int tt,k1,k2,nh,Nh,islc,islcpt,inttemp;
  int s=0;
  Int_t rmult=0;
  int mult=0;
  Int_t cbin=0;
  Float_t cnpart=0;
  Float_t cnpartsig=0;
  Float_t cncoll=0;
  Float_t cncollsig=0;
  //long t;
  long Totevt=0;
  long neve=0; 
  long Neve;
  long evr;
  long nentries;
  double rpmult=0.;
  double rmultmean=0.;
  double Trmultmean=0.;  
  Float_t phi0[20000];
  Float_t pt0[20000];
  Float_t eta0[20000];
  //double w2=1.;
  double temp=0.;
  double temp2=0.; 
  double temp4=0.; 
  double temp6=0.; 
  double temp8=0.; 
  double tempG[kr];
  double tmpr=0.;
  double tmpi=0.;
  double tmpreal=0.;
  double V2mean=0.;
  double err2mean=0.;
  double realG2=0.;
  double imagG2=0.;
  double chi2=0.;
  double sigma2=0.;
  double forG2[kr]; 
  double V2[kt];
  double r02[kt];  
  double theta[kt];
  double Qtheta[kt];
  double Qx[kt];
  double Qy[kt];
  double Qxa[kt];
  double Qya[kt];
  double Q2[kt];
  double RZ[kr];
  double G[kr];
  double err2[kt];
  double G2mod2[kt][kr];
  double rpmult2[NB_SLICE][NB_SLICE_PT];
  double totmult[NB_SLICE][NB_SLICE_PT]; 
//16 bins - hiSelectedTracks

  double minptcut[NB_SLICE_PT] = {0.2,0.7,1.1,1.4,1.75,2.2,2.75,3.45,4.2,5.05,6.0,7.2,8.5};
  double maxptcut[NB_SLICE_PT] = {0.7,1.1,1.4,1.75,2.2,2.75,3.45,4.2,5.05,6.0,7.2,8.5,10.0}; //13

  double v22[kt][NB_SLICE][NB_SLICE_PT]; 
  double v22mean[NB_SLICE][NB_SLICE_PT]; 
  double err22[kt][NB_SLICE][NB_SLICE_PT];
  double err22mean[NB_SLICE][NB_SLICE_PT];

  double v24[kt][NB_SLICE][NB_SLICE_PT];
  double v24mean[NB_SLICE][NB_SLICE_PT];
  double err24[kt][NB_SLICE][NB_SLICE_PT];
  double err24mean[NB_SLICE][NB_SLICE_PT];

  double v26[kt][NB_SLICE][NB_SLICE_PT];
  double v26mean[NB_SLICE][NB_SLICE_PT];
  double err26[kt][NB_SLICE][NB_SLICE_PT];
  double err26mean[NB_SLICE][NB_SLICE_PT];

  double v28[kt][NB_SLICE][NB_SLICE_PT];
  double v28mean[NB_SLICE][NB_SLICE_PT];
  double err28[kt][NB_SLICE][NB_SLICE_PT];
  double err28mean[NB_SLICE][NB_SLICE_PT];

  double pt[NB_SLICE_PT];
  double MeanPt[NB_SLICE_PT];
  double MofPt[NB_SLICE_PT];

  const std::complex<double> i(0.,1.);
  std::complex<double> G2[kt][kr];
  std::complex<double> g2[kt][kr];
  std::complex<double> g02[kt];
  std::complex<double> dlng02_dz[kt],dG02_dz[kt];
  std::complex<double> dlng02_dw[kt][NB_SLICE][NB_SLICE_PT],dG02_dw[kt][NB_SLICE][NB_SLICE_PT],compratio[kt][NB_SLICE][NB_SLICE_PT];
  std::complex<double> dlng024_dw[kt][NB_SLICE][NB_SLICE_PT],dG024_dw[kt][NB_SLICE][NB_SLICE_PT],compratio4[kt][NB_SLICE][NB_SLICE_PT];
  std::complex<double> dlng026_dw[kt][NB_SLICE][NB_SLICE_PT],dG026_dw[kt][NB_SLICE][NB_SLICE_PT],compratio6[kt][NB_SLICE][NB_SLICE_PT];
  std::complex<double> dlng028_dw[kt][NB_SLICE][NB_SLICE_PT],dG028_dw[kt][NB_SLICE][NB_SLICE_PT],compratio8[kt][NB_SLICE][NB_SLICE_PT];
  //cout.precision(6);
  for(k1=0;k1<kt;k1++) dG02_dz[k1]=0.;

  for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) {   
    pt[islcpt]= minptcut[islcpt]/2. + maxptcut[islcpt]/2.;    
    MeanPt[islcpt]=0.;
    MofPt[islcpt]=0.;

    for (islc=0; islc<NB_SLICE; islc++) {
      //w22[islc][islcpt]=1.; 
      //w1[islc][islcpt]=1.;
      v22mean[islc][islcpt]=0.;
      err22mean[islc][islcpt]=0.;
      v24mean[islc][islcpt]=0.;
      err24mean[islc][islcpt]=0.;
      v26mean[islc][islcpt]=0.;
      err26mean[islc][islcpt]=0.;
      v28mean[islc][islcpt]=0.;
      err28mean[islc][islcpt]=0.;
      rpmult2[islc][islcpt]=0.;
      totmult[islc][islcpt]=0.;
    for(k1=0; k1<kt; k1++) 
     { 
       dG02_dw[k1][islc][islcpt]=0.;
       dlng02_dw[k1][islc][islcpt]=0.;
       v22[k1][islc][islcpt]=0.;  
       err22[k1][islc][islcpt]=0.;

       dG024_dw[k1][islc][islcpt]=0.;
       dlng024_dw[k1][islc][islcpt]=0.;
       v24[k1][islc][islcpt]=0.;  
       err24[k1][islc][islcpt]=0.;

       dG026_dw[k1][islc][islcpt]=0.;
       dlng026_dw[k1][islc][islcpt]=0.;
       v26[k1][islc][islcpt]=0.;  
       err26[k1][islc][islcpt]=0.;

       dG028_dw[k1][islc][islcpt]=0.;
       dlng028_dw[k1][islc][islcpt]=0.;
       v28[k1][islc][islcpt]=0.;  
       err28[k1][islc][islcpt]=0.;
	}
    }
    }
     for(k2=0; k2<kr; k2++){
       tempG[k2]=0.;
     }
     for(k1=0;k1<kt;k1++)  for(k2=0;k2<kr;k2++)
       {
         G2[k1][k2]=0.;
         g2[k1][k2]=0.;
         G2mod2[k1][k2]=0.;
       }
  for(k1=0; k1<kt; k1++){
    theta[k1]=((double) k1)*Pi/((double) kt); 
    Qxa[k1]=0.;
    Qya[k1]=0.;
    Q2[k1]=0.;
  }

  outR.open("R.dat");  
  for(k2=0; k2<kr; k2++){
    RZ[k2]=rootJ0/(Vmax-Vstep*k2);  
    outR<<RZ[k2]<<endl;
    forG2[k2]=0.;
  } 
  outR.close();

     rands = new TRandom();
     //TF1 *effpt = new TF1("effpt","exp(-(x-1.9)*(x-1.9)/1)/14+0.8-exp((-x+0.8))/4",0,10);
     //TF1 *effeta =new TF1 ("effeta","0.780831e*exp(-(x-0.00212541)*(x-0.00212541)/2/1.94182/1.94182)",-2.4,2.4);
     char *file1 = "../../vndata_5k.root";
     TFile *f1 = new TFile(file1);
     TTree *t1 = (TTree*)f1->Get("tree");
     t1->SetBranchAddress("npg", &rmult);
     t1->SetBranchAddress("ptg", &pt0);
     t1->SetBranchAddress("etag", &eta0);
     t1->SetBranchAddress("phig", &phi0);
     neve=t1->GetEntries();
     Neve=0;
     cout<<"Nevts = "<<neve<<endl;
     for(long ne=0; ne<neve; ne++)
     {
       if(ne%1000==0)  cout<<"Have run "<<ne<<" events ,"<<" First pass"<<endl;
       t1->GetEntry(ne);
       cent->Fill(rmult);

       //if(cbin<8||cbin>11) continue;
       Neve++;
       Nh=0;
       //rmultmean+=((float) rmult);
       for(k1=0; k1<kt; k1++)   for(k2=0; k2<kr; k2++) g2[k1][k2]=1.;
       for(k1=0; k1<kt; k1++)  {
	 Qtheta[k1]=0.;
	 Qx[k1]=0.;
	 Qy[k1]=0.;
       }
       for(nh=0;nh<rmult;nh++)

	 {
	   //if(eta0[nh]<-0.35||eta0[nh]>0.35) continue;
           if(eta0[nh]<-2.40||eta0[nh]>2.40) continue;           
           //if(eta0[nh]<-4.0||eta0[nh]>4.0) continue;           
           if(pt0[nh]<0.2||pt0[nh]>10.0) continue;
           hptt->Fill(pt0[nh]);
           hphit->Fill(phi0[nh]);
           hetat->Fill(eta0[nh]);
       //    if(rands->Uniform(0.0,1.0)>effeta->Eval(eta0[nh])) continue;
	   hpt->Fill(pt0[nh]);
	   hphi->Fill(phi0[nh]);
	   heta->Fill(eta0[nh]);
           for(k1=0; k1<kt; k1++)
	     {	 
               Qtheta[k1]+=(1.*cos(3.*phi0[nh]-theta[k1])); 
               //Qtheta[k1]+=(1.*cos(2.*phi0[Nh]-theta[k1]));
	       Qx[k1]+=cos(3.*phi0[nh]);
	       Qy[k1]+=sin(3.*phi0[nh]); 
	     }
	   Nh++;
	 } // end of loop over particles
       rmultmean+=((float) Nh);
       cent_used->Fill(Nh); 

       for(k1=0; k1<kt; k1++) {
	 Q2[k1]+=Qx[k1]*Qx[k1]+Qy[k1]*Qy[k1];
	 Qxa[k1]+=Qx[k1];
	 Qya[k1]+=Qy[k1];
       for(k2=0; k2<kr; k2++) {

	 G2[k1][k2]+=(cos(RZ[k2]*Qtheta[k1])+i*sin(RZ[k2]*Qtheta[k1])); 
       }    
       }
     } // end of loop over events

       Trmultmean=rmultmean;

     Totevt=Neve;

     cout<<"Totevt: "<<Totevt <<  " Trmult: " <<Trmultmean<<",  <Mult> = "<<Trmultmean/Totevt<<endl;

  rmultmean=Trmultmean/Totevt;
  for(k1=0; k1<kt; k1++) 
    { 
    sigma2+=(Q2[k1]/Totevt-(Qxa[k1]/Totevt)*(Qxa[k1]/Totevt)-(Qya[k1]/Totevt)*(Qya[k1]/Totevt));
      for(k2=0; k2<kr; k2++){  
      G2mod2[k1][k2]=sqr(abs(G2[k1][k2])/Totevt); 
      tempG[k2]+=abs(G2[k1][k2])/Totevt; 
      }
      k2=0;
      while (G2mod2[k1][k2]>G2mod2[k1][k2+1]){
	k2++; 
        //if(k2==199)
	//{ cout<<"G2mod2 Array index out of range! "<<endl; s++;   break; }
      }
      V2[k1]=Vmax-Vstep*(k2+(G2mod2[k1][k2-1]-G2mod2[k1][k2+1])/2./(G2mod2[k1][k2-1]-2.*G2mod2[k1][k2]+G2mod2[k1][k2+1]));
      r02[k1]=rootJ0/V2[k1]; 
      V2mean+=V2[k1]; 
    }      
    V2mean/=kt; 
    sigma2/=kt;
    sigma2-=V2mean*V2mean;
    //chi2=V2mean/sqrt(rmultmean);
    chi2=V2mean/sqrt(sigma2);
    cout<<"V2mean/sqrt(sigma2) = "<<chi2<<" ,  V2mean/sqrt(rmultmean) ="<<V2mean/sqrt(rmultmean)<<endl;


    outG.open("G.dat");
      for(k2=0; k2<kr; k2++){                                     
	outG<<tempG[k2]/kt<<endl;
      }
    outG.close();

   plotVint.open("plots_Vint.grf");
  for(k1=0; k1<kt; k1++) 
    {   
      err2[k1]=V2[k1]*sqrt((exp(-sqr(rootJ0/chi2)/2.)*BesselJ0(2.*rootJ0)+exp(sqr(rootJ0/chi2)/2.))/2./Totevt)/rootJ0/J1rootJ0;      

      plotVint<<theta[k1]/2.<<"   " <<V2[k1]/rmultmean<<"   " <<err2[k1]/rmultmean <<endl;
    }

     plotVint<<"theta      <V2[k2]>_evt      V2[k2]err" <<endl;

  temp2=0.; 
  for(k1=0; k1<kt; k1++) 

    temp2+=exp(sqr(rootJ0/chi2)*cos(theta[k1])/2.)*BesselJ0(2.*rootJ0*sin(theta[k1]/2.))
      +exp(-sqr(rootJ0/chi2)*cos(theta[k1])/2.)*BesselJ0(2.*rootJ0*cos(theta[k1]/2.));
  err2mean=V2mean*sqrt(temp2/2./Totevt/kt)/rootJ0/J1rootJ0;
    plotVint<<V2mean/rmultmean <<"   "<<err2mean/rmultmean<<endl;
    plotVint<<"V2" <<"   "<<"V2err"<<endl;
    plotVint.close();


     char *file1 = "../../vndata_5k.root";
     TFile *f1 = new TFile(file1);
     TTree *t1 = (TTree*)f1->Get("tree");
     t1->SetBranchAddress("npg", &rmult);
     t1->SetBranchAddress("ptg", &pt0);
     t1->SetBranchAddress("etag", &eta0);
     t1->SetBranchAddress("phig", &phi0);

     neve=t1->GetEntries();
     cout<<"Nevts = "<<neve<<endl;
   for(long ne=0; ne<neve; ne++)
     {
       if(ne%1000==0) cout<<"Have run "<<ne<<" events ,"<<" Second pass"<<endl;
       t1->GetEntry(ne);

       //if(cbin<8||cbin>11) continue;

    for(k1=0; k1<kt; k1++) {
	g02[k1]=1.;
        Qtheta[k1]=0.; 
        dlng02_dz[k1]=0.; 
  for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) {
    for (islc=0; islc<NB_SLICE; islc++) {
      dlng02_dw[k1][islc][islcpt]=0.;
      dlng024_dw[k1][islc][islcpt]=0.;
      dlng026_dw[k1][islc][islcpt]=0.;
      dlng028_dw[k1][islc][islcpt]=0.;
	  }
      }
    }
 

      Nh=0;
       for(nh=0;nh<rmult;nh++)
	 {
      //     if(rands->Uniform(0.0,1.0)>effeta->Eval(eta0[nh])) continue;
//	   if(eta0[nh]<-0.35||eta0[nh]>0.35) continue;
           //if(eta0[nh]<-4.0||eta0[nh]>4.0) continue;           
           if(eta0[nh]<-2.40||eta0[nh]>2.40) continue;
           if(pt0[nh]<0.2||pt0[nh]>10.0) continue;
      for(islcpt=0;islcpt<NB_SLICE_PT;islcpt++)  
	{
          //if((eta0[nh]>minptcut[islcpt]&&eta0[nh]<=maxptcut[islcpt]))
          if((pt0[nh]>minptcut[islcpt]&&pt0[nh]<=maxptcut[islcpt]))
	     inttemp=islcpt;
	}
      islcpt=inttemp;
      islc=0;
      rpmult2[islc][islcpt]+=1.; 
//      MeanPt[islcpt]+=eta0[nh];
      MeanPt[islcpt]+=pt0[nh];
      MofPt[islcpt]+=1.;
      temp2=0.;

      for(k1=0; k1<kt; k1++) 
	{
//	  //dlng02_dw[k1][islc][islcpt]+=(cos(2.*phi0[nh]-theta[k1])*(cos(r02[k1]*Qtheta[k1])+i*sin(r02[k1]*Qtheta[k1])));
//          dlng02_dw[k1][islc][islcpt]+=(cos(2.0*(2.*phi0[nh]-theta[k1]))*(cos(r02[k1]*Qtheta[k1])+i*sin(r02[k1]*Qtheta[k1])));
//          //dlng02_dw[k1][islc][islcpt]+=(cos(2.0*2.*phi0[nh]-theta[k1])*(cos(r02[k1]*Qtheta[k1])+i*sin(r02[k1]*Qtheta[k1])));
          temp2=1.*cos(3.*phi0[nh]-theta[k1]);
          g02[k1]*=(1.+i*r02[k1]*temp2);
          dlng02_dz[k1]+=(temp2/(1.+i*r02[k1]*temp2));
          //dlng02_dw[k1][islc][islcpt]+=(cos(2.*phi0[nh]-theta[k1])/(1.+i*r02[k1]*temp2));
          dlng02_dw[k1][islc][islcpt]+=(cos(1.0*(3.*phi0[nh]-theta[k1]))/(1.+i*r02[k1]*temp2));
          dlng024_dw[k1][islc][islcpt]+=(cos(2.0*(3.*phi0[nh]-theta[k1]))/(1.+i*r02[k1]*temp2));
          dlng026_dw[k1][islc][islcpt]+=(cos(3.0*(3.*phi0[nh]-theta[k1]))/(1.+i*r02[k1]*temp2));
          dlng028_dw[k1][islc][islcpt]+=(cos(4.0*(3.*phi0[nh]-theta[k1]))/(1.+i*r02[k1]*temp2));

	}  
      Nh++;
	 } // end of loop over particles

//       for(k1=0; k1<kt; k1++)   dlng02_dz[k1]=(Qtheta[k1]*(cos(r02[k1]*Qtheta[k1])+i*sin(r02[k1]*Qtheta[k1])));
    for(k1=0; k1<kt; k1++) {
//      dG02_dz[k1]+=dlng02_dz[k1]; 
      dG02_dz[k1]+=dlng02_dz[k1]*g02[k1];
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) {
    for (islc=0; islc<NB_SLICE; islc++) {
//      dG02_dw[k1][islc][islcpt]+=dlng02_dw[k1][islc][islcpt];
      dG02_dw[k1][islc][islcpt]+=(dlng02_dw[k1][islc][islcpt]*g02[k1]);
      dG024_dw[k1][islc][islcpt]+=(dlng024_dw[k1][islc][islcpt]*g02[k1]);
      dG026_dw[k1][islc][islcpt]+=(dlng026_dw[k1][islc][islcpt]*g02[k1]);
      dG028_dw[k1][islc][islcpt]+=(dlng028_dw[k1][islc][islcpt]*g02[k1]);
      }     
    } 
    }
     } // end of loop over events


    output.open("output.dat");   
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) {
      MeanPt[islcpt]/=MofPt[islcpt];
      //cout<<"pt= "<<pt[islcpt]<<", Mean pt = "<<MeanPt[islcpt]<<endl;
    for (islc=0; islc<NB_SLICE; islc++) {
      totmult[islc][islcpt]=rpmult2[islc][islcpt];
      output<<rpmult2[islc][islcpt]/Totevt<<endl;
      //cout<<"  rpmult2 = "<<rpmult2[islc][islcpt<<endl;
	  } 
    }

      output<<"chi2 = "<<chi2<<endl; 

      output.close();

    for(k1=0; k1<kt; k1++){
     dG02_dz[k1]/=Totevt;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) {
    for (islc=0; islc<NB_SLICE; islc++) {
	    dG02_dw[k1][islc][islcpt]/=rpmult2[islc][islcpt]; 
	    dG024_dw[k1][islc][islcpt]/=rpmult2[islc][islcpt]; 
	    dG026_dw[k1][islc][islcpt]/=rpmult2[islc][islcpt]; 
	    dG028_dw[k1][islc][islcpt]/=rpmult2[islc][islcpt]; 
	  }
    }
    }

  cout<<"Totevt = "<<Totevt<<"  Trmult = "<<Trmultmean<<endl;

  cout<<"For check#2, <Re>: "<<dG02_dz[3].real()<<" Im: "<<dG02_dz[3].imag()<<endl;

  //////// calcu of v2 .
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) {
    for (islc=0; islc<NB_SLICE; islc++) {
	    totmult[islc][islcpt]/=Totevt;

      for(k1=0; k1<kt; k1++){
//	compratio[k1][islc][islcpt]=dG02_dw[k1][islc][islcpt]/dG02_dz[k1];
        compratio[k1][islc][islcpt]=dG02_dw[k1][islc][islcpt]/dG02_dz[k1];
        compratio4[k1][islc][islcpt]=dG024_dw[k1][islc][islcpt]/dG02_dz[k1]/i;
        compratio6[k1][islc][islcpt]=dG026_dw[k1][islc][islcpt]/dG02_dz[k1]/i/i;
        compratio8[k1][islc][islcpt]=dG028_dw[k1][islc][islcpt]/dG02_dz[k1]/i/i/i;

//	    v22[k1][islc][islcpt]=V2[k1]*compratio[k1][islc][islcpt].real();
//            v22[k1][islc][islcpt]=rootJ0/2.0*V2[k1]*compratio[k1][islc][islcpt].real();
            v22[k1][islc][islcpt]=V2[k1]*compratio[k1][islc][islcpt].real();  // J10/J30 = 2.60878
	    err22[k1][islc][islcpt]=sqrt((exp(sqr(rootJ0/chi2)/2.) + exp(-sqr(rootJ0/chi2)/2.)*BesselJ0(2.*rootJ0))/totmult[islc][islcpt])/2./J1rootJ0; 
	    v22mean[islc][islcpt]+=v22[k1][islc][islcpt];

            v24[k1][islc][islcpt]=rootJ0/2.0*V2[k1]*compratio4[k1][islc][islcpt].real();  // J10/J30 = 2.60878
	    err24[k1][islc][islcpt]=sqrt((exp(sqr(rootJ0/chi2)/2.) + exp(-sqr(rootJ0/chi2)/2.)*BesselJ0(2.*rootJ0))/totmult[islc][islcpt])/2./J1rootJ0; 
	    v24mean[islc][islcpt]+=v24[k1][islc][islcpt];

            v26[k1][islc][islcpt]=2.60878*V2[k1]*compratio6[k1][islc][islcpt].real();  // J10/J30 = 2.60878
	    err26[k1][islc][islcpt]=sqrt((exp(sqr(rootJ0/chi2)/2.) + exp(-sqr(rootJ0/chi2)/2.)*BesselJ0(2.*rootJ0))/totmult[islc][islcpt])/2./J1rootJ0; 
	    v26mean[islc][islcpt]+=v26[k1][islc][islcpt];

            v28[k1][islc][islcpt]=8.01813*V2[k1]*compratio8[k1][islc][islcpt].real();  // J10/J40 = 8.01813
	    err28[k1][islc][islcpt]=sqrt((exp(sqr(rootJ0/chi2)/2.) + exp(-sqr(rootJ0/chi2)/2.)*BesselJ0(2.*rootJ0))/totmult[islc][islcpt])/2./J1rootJ0; 
	    v28mean[islc][islcpt]+=v28[k1][islc][islcpt];

	    if(k1==3)  cout<<"  real(num/den)= "<<compratio[k1][islc][islcpt].real()<<"  real(num/den)4= "<<compratio4[k1][islc][islcpt].real()<<"  real(num/den)6= "<<compratio6[k1][islc][islcpt].real()<<endl;
    }
    }
    }

   plotvdiff.open("plots_vdiff.grf");
    for(k1=0; k1<kt; k1++){
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) {
    for (islc=0; islc<NB_SLICE; islc++) {
      plotvdiff<<v22[k1][islc][islcpt]<<"   "<<err22[k1][islc][islcpt]<<endl;
    }
    }
    }

    plotvdiff<<"Totevt = "<<Totevt<<" Tottrack =  "<<Trmultmean<<" <Mult> = "<<Trmultmean/Totevt<<endl;
    plotvdiff.close();

    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) { 
    for (islc=0; islc<NB_SLICE; islc++) {
	v22mean[islc][islcpt]/=kt;
	v24mean[islc][islcpt]/=kt;
	v26mean[islc][islcpt]/=kt;
	v28mean[islc][islcpt]/=kt;
    }
    }

   temp2=0.;
   temp4=0.;
   temp6=0.;
   temp8=0.;
   for(k1=0; k1<kt; k1++) {
    temp2+=(exp(sqr(rootJ0/chi2)*cos(theta[k1])/2.)*
	    BesselJ0(2.*rootJ0*sin(theta[k1]/2.)) -
	    exp(-sqr(rootJ0/chi2)*cos(theta[k1])/2.)*
	    BesselJ0(2.*rootJ0*cos(theta[k1]/2.)))*cos(theta[k1]);
      temp4+=(exp(sqr(rootJ0/chi2)*cos(theta[k1])/2.)*BesselJ0(2.*rootJ0*sin(theta[k1]/2.)) +
      exp(-sqr(rootJ0/chi2)*cos(theta[k1])/2.)*BesselJ0(2.*rootJ0*cos(theta[k1]/2.)))*cos(2.0*theta[k1]);
      temp6+=(exp(sqr(rootJ0/chi2)*cos(theta[k1])/2.)*BesselJ0(2.*rootJ0*sin(theta[k1]/2.)) -
      exp(-sqr(rootJ0/chi2)*cos(theta[k1])/2.)*BesselJ0(2.*rootJ0*cos(theta[k1]/2.)))*cos(3.0*theta[k1]); 
      temp8+=(exp(sqr(rootJ0/chi2)*cos(theta[k1])/2.)*BesselJ0(2.*rootJ0*sin(theta[k1]/2.)) +
      exp(-sqr(rootJ0/chi2)*cos(theta[k1])/2.)*BesselJ0(2.*rootJ0*cos(theta[k1]/2.)))*cos(4.0*theta[k1]);
     }
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
    for (islc=0; islc<NB_SLICE; islc++){ 
       err22mean[islc][islcpt]=sqrt(temp2/totmult[islc][islcpt]/Totevt/kt)/2./J1rootJ0;
       err24mean[islc][islcpt]=sqrt(temp2/totmult[islc][islcpt]/Totevt/kt)/2./J2rootJ0;
       err26mean[islc][islcpt]=sqrt(temp6/totmult[islc][islcpt]/Totevt/kt)/2./0.199;  // J30 = 0.199 
       err28mean[islc][islcpt]=sqrt(temp8/totmult[islc][islcpt]/Totevt/kt)/2./0.0647467;  // J40 = 0.064746
     }
    //err22mean[bin]=sqrt(temp2/rpmult[bin]/kt)/2./J1rootJ0;
    outv2.open("v2.dat");
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv2<<MeanPt[islcpt]<<"   "<<v22mean[0][islcpt]<<"   "<<err22mean[0][islcpt]<<endl;
      outv2.close();

    outv4.open("v4.dat");
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv4<<MeanPt[islcpt]<<"   "<<v24mean[0][islcpt]<<"   "<<err24mean[0][islcpt]<<endl;
      outv4.close();

    outv6.open("v6.dat");
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv6<<MeanPt[islcpt]<<"   "<<v26mean[0][islcpt]<<"   "<<err26mean[0][islcpt]<<endl;
      outv6.close();

    outv8.open("v8.dat");
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
      outv8<<MeanPt[islcpt]<<"   "<<v28mean[0][islcpt]<<"   "<<err28mean[0][islcpt]<<endl;
      outv8.close();

    cout<<" v2 : "<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
        cout<<MeanPt[islcpt]<<"   "<<v22mean[0][islcpt]<<"   "<<err22mean[0][islcpt]<<endl; 

    cout<<" v4 : "<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
        cout<<MeanPt[islcpt]<<"   "<<v24mean[0][islcpt]<<"   "<<err24mean[0][islcpt]<<endl; 

    cout<<" v6 : "<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
        cout<<MeanPt[islcpt]<<"   "<<v26mean[0][islcpt]<<"   "<<err26mean[0][islcpt]<<endl; 

    cout<<" v8 : "<<endl;
    for (islcpt=0; islcpt<NB_SLICE_PT; islcpt++) 
        cout<<MeanPt[islcpt]<<"   "<<v28mean[0][islcpt]<<"   "<<err28mean[0][islcpt]<<endl; 
 


  TFile *f = new TFile("dist_mult.root","recreate");
  f->cd();
 
  cent->Write();
  cent1->Write();
  cent2->Write();
  cent3->Write();
  cent4->Write();
  cent5->Write();
  cent6->Write();
  cent7->Write();
  cent8->Write();
  cent9->Write();
  cent10->Write();
  cent_used->Write();
  npart_used->Write();
  npartsig_used->Write();
  ncoll_used->Write();
  ncollsig_used->Write();
  hpt->Write();
  hphi->Write();
  heta->Write();
  hptt->Write();
  hphit->Write();
  hetat->Write();

  f->Close();

 
}
 
