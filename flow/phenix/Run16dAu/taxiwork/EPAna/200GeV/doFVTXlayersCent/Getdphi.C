#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "TString.h"
#include "TProfile.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "RpPar.h"

using namespace std;

int GetTotalRun();
void FillGoodRun(int);
float GetReso(int , int, int , int , int, bool);
TString GetRun(int);
TString pro = "";
int taxi = 8760;
float eps = 1e-4;
float GoodRunFit[ncent][nhar][nsub][10000];

TString choosesub(int isub){
    TString str;
//        if(isub==0)
//         str = "FVTX1LS";
//        else if(isub==1)
//         str = "FVTX2LS";
//        else if(isub==2)
//         str = "FVTX3LS";
//        else if(isub==3)
//         str = "FVTX4LS";
//        else if(isub==5)
        if(isub==5)
         str = "FVTX1S";
        else if(isub==4)
         str = "BBCS";
        else if(isub==6)
         str = "FVTX1p2p3LS";
        else
         str = "ABORT";
    return str;
}

TString choosesub1(int isub){
    TString str;
    if(choosesub(isub).Contains("FVTX"))
        str = "BBCS";
    else if(choosesub(isub).Contains("BBC"))
        str = "FVTX1p2p3LS";
    else 
        str = "ABORT";
        return str;
}


void Getdphi(int ihar=1, int iangle1=0, int iangle2=0, bool usingCNTEP=0){
    TString str;
    int nrun = GetTotalRun();
    std::cout<<"Totally we have "<<nrun<<" runs/segments!"<<std::endl;
    FillGoodRun(ihar);
    std::cout<<"Filling Good run finished!"<<std::endl;

    if(nrun<0) exit(1);

     int iharE=0;
     if(nhar==1 || nhar ==2) iharE=1;
     TFile *fin;

    int n = ihar+1.0+iharE;
    cout<<"iangle1 = "<<iangle1<<" iangle2 = "<<iangle2<<endl;
    TFile *fout = new TFile(Form("dphiv%d.root",n),"recreate");

     for(int icent=0;icent<ncent;icent++){
//      for(int ihar=0;ihar<nhar;ihar++){
//          if(ihar!=1) continue;
       for(int isub=0;isub<nsub;isub++){
        str = choosesub(isub);
        TString UseCNTEP;
        if(str=="ABORT") continue;
        if(usingCNTEP)
         UseCNTEP = "UseCNTEP";
        else
         UseCNTEP = "NoUseCNTEP";
        std::cout<<UseCNTEP<<std::endl;
        std::cout<<"starting doing "<<str<<" v"<<n<<" analysis!"<<std::endl;
//         float reso = GetReso(iangle1, iangle2, icent,ihar,isub,usingCNTEP);
//         if(reso<=0) {std::cout<<"resolution is wrong!"<<std::endl; reso = 1.0;}
         TH2F* hvall = new TH2F(Form("hdphiall_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),Form("hdphiall_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),60,0,6,200,-4,4);
         TH2F* hvnall = new TH2F(Form("hdphinall_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),Form("hdphinall_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),60,0,6,100,-2,2);

        for(int iphi=0;iphi<nphi+1;iphi++){
         TH2F* hv = new TH2F(Form("hdphi_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),Form("hdphi_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),60,0,6,200,-4,4);
         TH2F* hvn = new TH2F(Form("hdphin_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),Form("hdphin_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),60,0,6,100,-2,2);
         string phistr = (iphi==0)?"_east":"_west";
         if(iphi==nphi) phistr = "";
         if(iphi<nphi){
        for(int irun=0;irun<nrun;irun++){
         fin = TFile::Open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/work/200GeV/output/%s",GetRun(irun).Data()));
         if(!(GoodRunFit[icent][ihar][isub][irun]>0.2 && GoodRunFit[icent][ihar][isub][irun]<3.0)){
         std::cout<<"cent = "<<icent<<"; n = "<<n<<" ;isub = "<<str<<" ;run = "<<GetRun(irun)<<" is bad run!"<<std::endl;
         fin->Close();
        continue;
         }
         TH2F *hvtemp;
         TH2F *hvsqtemp;
         hvtemp = (TH2F*)fin->Get(Form("v%s_%d_%d_%d_%d_%d",str.Data(),iangle1,iangle2,icent,ihar,iphi));
         hvntemp = (TH2F*)fin->Get(Form("vn%s_%d_%d_%d_%d_%d",str.Data(),iangle1,iangle2,icent,ihar,iphi));
         hv->Add(hvtemp);
         hvn->Add(hvntemp);
         fin->Close();
        }
         }
        hvall->Add(hv);
        hvnall->Add(hvn);
        if(iphi==nphi){
        hv = hvall;
        hvn = hvnall;
        }

        fout->cd();
        hv->Write();
        hvn->Write();

         for(int ipt=0;ipt<npt-1;ipt++){
             int xbinmin = hv->GetXaxis()->FindBin(ptbin[ipt]+eps);
             int xbinmax = hv->GetXaxis()->FindBin(ptbin[ipt+1]-eps);
           //  std::cout<<xbinmin<<" "<<xbinmax<<std::endl;
           //  std::cout<<ptbin[ipt]<<" "<<ptbin[ipt+1]<<std::endl;
            TH1F* hvProj = (TH1F*)hv->ProjectionY(Form("hvProj_%d",ipt),xbinmin,xbinmax);
            TH1F* hvnProj = (TH1F*)hvn->ProjectionY(Form("hvnProj_%d",ipt),xbinmin,xbinmax);
         }
         }
        }
       // }
     }
     fout->Close();

}

void FillGoodRun(int ihar){
    float pi = acos(-1.0);
    TString str;
    TFile *fin;
    int nrun = GetTotalRun();
    if(nrun<0) exit(1);
     for(int icent=0;icent<ncent;icent++){
 //     for(int ihar=0;ihar<nhar;ihar++){
 //         if(ihar!=1) continue;
       for(int isub=0;isub<nsub;isub++){
            str = choosesub(isub);
            if(str=="ABORT") continue;
            for(int irun=0;irun<nrun;irun++){
      //std::cout<<icent<<" "<<ihar<<" "<<isub<<" "<<irun<<std::endl;
        //fin = TFile::Open(Form("/phenix/plhf/xuq/taxi/%s%s/%d/data/%s.root",dataset.Data(),pro.Data(),taxi,GetRun(irun).Data()));
         fin = TFile::Open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/work/200GeV/output/%s",GetRun(irun).Data()));
        TH1F* hpsi = new TH1F("psi","psi",100,-pi,pi);
        for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
          TH1F* hpsitemp = (TH1F*)fin->Get(Form("psiFla_0_0_%d_%d_%d_%d",icent,ibbcz,ihar,isub));
          hpsi->Add(hpsitemp);
        }
        TF1 *fun = new TF1("fun","pol0",-pi,pi);
      if(hpsi->GetEntries()>1000){
	//hpsi->SetMarkerStyle(20);
	//hpsi->SetMarkerSize(0.6);
	//hpsi->SetMarkerColor(4);
	hpsi->SetMinimum(10);
	hpsi->Fit("fun","QR0");
	//float par=fun->GetParameter(0);
	//hpsi->SetMaximum(1.5*par);
	//hpsi->Draw();
	GoodRunFit[icent][ihar][isub][irun] = fun->GetChisquare()/fun->GetNDF();
        fin->Close();
      }
      else{
        GoodRunFit[icent][ihar][isub][irun] = -9999;
        fin->Close();
     }
     // GoodRunFit[icent][ihar][isub][irun] = 1.;
    }
       }
   //   }
     }
}


int GetTotalRun(){
    ifstream frun(dataset+".Lst");
    string RunNumber;
    int nrun = -1;
    while(!frun.eof()){
        frun>>RunNumber;
        nrun++;
    }
    return nrun;
}


TString GetRun(int irun){
    ifstream frun(dataset+".Lst");
    string RunNumbertemp;
    int jrun = 0;
    while(!frun.eof() && jrun<=irun){
        frun>>RunNumbertemp;
        jrun++;
    }
    TString RunNumber(RunNumbertemp);
    return RunNumber;
}
