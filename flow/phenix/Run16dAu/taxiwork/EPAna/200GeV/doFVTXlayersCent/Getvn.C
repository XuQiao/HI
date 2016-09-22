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


void Getvn(int ihar, int iangle1, int iangle2, bool usingCNTEP=0){
    TString str;
    int nrun = GetTotalRun();
    std::cout<<"Totally we have "<<nrun<<" runs/segments!"<<std::endl;
    FillGoodRun(ihar);
    std::cout<<"Filling Good run finished!"<<std::endl;

    if(nrun<0) exit(1);

     ofstream fout, foutraw, fout1, fout2;
     int iharE=0;
     if(nhar==1) iharE=1;
     TFile *fin;

    cout<<"iangle1 = "<<iangle1<<" iangle2 = "<<iangle2<<endl;

     for(int icent=0;icent<ncent;icent++){
//      for(int ihar=0;ihar<nhar;ihar++){
//          if(ihar!=1) continue;
       for(int isub=0;isub<nsub;isub++){
        int n = ihar+1.0+iharE;
        str = choosesub(isub);
        TString UseCNTEP;
        if(str=="ABORT") continue;
        if(usingCNTEP)
         UseCNTEP = "UseCNTEP";
        else
         UseCNTEP = "NoUseCNTEP";
        std::cout<<UseCNTEP<<std::endl;
        std::cout<<"starting doing "<<str<<" v"<<n<<" analysis!"<<std::endl;
         fout1.open(Form("Result/%s/res_%d%d_%d_%d_%s.dat",UseCNTEP.Data(),iangle1,iangle2,n,icent,str.Data())); //using str as event plane detector
         fout2.open(Form("Result/%s/psi_%d%d_%d_%d_%s.dat",UseCNTEP.Data(),iangle1,iangle2,n,icent,str.Data())); //using str as event plane detector
         float reso = GetReso(iangle1, iangle2, icent,ihar,isub,usingCNTEP);
         fout1<<reso<<std::endl;
         if(reso<=0) {std::cout<<"resolution is wrong!"<<std::endl; reso = 1.0;}
        for(int irun=0;irun<nrun;irun++){
         fout2<<GetRun(irun)<<" "<<GoodRunFit[icent][ihar][isub][irun]<<std::endl;
        }
         TH1D* hvobsall = new TH1D(Form("hvobsall_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),Form("hvobsall_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),60,0,6);
         TH1D* hvobsallsq = new TH1D(Form("hvobsallsq_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),Form("hvobsallsq_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),60,0,6);
         TH1D* hvobs2all = new TH1D(Form("hvobs2all_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),Form("hvobs2all_%d%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub),60,0,6);

        for(int iphi=0;iphi<nphi+1;iphi++){
         TH1D* hvobs = new TH1D(Form("hvobs_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),Form("hvobs_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),60,0,6);
         TH1D* hvobssq = new TH1D(Form("hvobssq_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),Form("hvobssq_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),60,0,6);
         TH1D* hvobs2 = new TH1D(Form("hvobs2_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),Form("hvobs2_%d%d_%d_%d_%d_%d",iangle1,iangle2,icent,ihar,isub,iphi),60,0,6);
         string phistr = (iphi==0)?"_east":"_west";
         if(iphi==nphi) phistr = "";
         cout<<"open v"<<n<<" file"<<endl;
         fout.open(Form("Result/%s/v%d_%d%d_%d%s_%s.dat",UseCNTEP.Data(),n,iangle1,iangle2,icent,phistr.c_str(),str.Data())); //using str as event plane detector
         cout<<"open v"<<n<<" raw file"<<endl;
         foutraw.open(Form("Result/%s/v%draw_%d%d_%d%s_%s.dat",UseCNTEP.Data(),n,iangle1,iangle2,icent,phistr.c_str(),str.Data())); //using str as event plane detector
         if(iphi<nphi){
        for(int irun=0;irun<nrun;irun++){
         //std::cout<<"cent = "<<icent<<"; n = "<<n<<" ;isub = "<<str<<" ;run = "<<irun<<" "<<phistr<<std::endl;
         //fin = TFile::Open(Form("/phenix/plhf/xuq/taxi/%s%s/%d/data/%s",dataset.Data(),pro.Data(),taxi,GetRun(irun).Data()));
         fin = TFile::Open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/work/200GeV/output/%s",GetRun(irun).Data()));
         if(!(GoodRunFit[icent][ihar][isub][irun]>0.2 && GoodRunFit[icent][ihar][isub][irun]<3.0)){
         std::cout<<"cent = "<<icent<<"; n = "<<n<<" ;isub = "<<str<<" ;run = "<<GetRun(irun)<<" is bad run!"<<std::endl;
         fin->Close();
        continue;
         }
         TProfile *hvobstemp;
         TProfile *hvobssqtemp;
         hvobstemp = (TProfile*)fin->Get(Form("vobs%s_%d_%d_%d_%d_%d",str.Data(),iangle1,iangle2,icent,ihar,iphi));
         hvobssqtemp = (TProfile*)fin->Get(Form("vobs%ssq_%d_%d_%d_%d_%d",str.Data(),iangle1,iangle2,icent,ihar,iphi));
         TH1D* hvobssumtemp = (TH1D*)hvobstemp->ProjectionX(Form("vobssum%s_%d_%d_%d",str.Data(),icent,ihar,iphi),"W");
         TH1D* hvobssqsumtemp = (TH1D*)hvobssqtemp->ProjectionX(Form("vobssqsum%s_%d_%d_%d",str.Data(),icent,ihar,iphi),"W"); //Add weighted v2
         TH1D* hvobssum2temp = (TH1D*)hvobstemp->ProjectionX(Form("vobssum2%s_%d_%d_%d",str.Data(),icent,ihar,iphi),"B");//Add Entries
         hvobs->Add(hvobssumtemp);
         hvobssq->Add(hvobssqsumtemp);
         hvobs2->Add(hvobssum2temp);
         fin->Close();
        }
         }
        hvobsall->Add(hvobs);
        hvobsallsq->Add(hvobssq);
        hvobs2all->Add(hvobs2);
        if(iphi==nphi){
        hvobs = hvobsall;
        hvobssq = hvobsallsq;
        hvobs2 = hvobs2all;
        }
         for(int ipt=0;ipt<npt-1;ipt++){
             int xbinmin = hvobs->GetXaxis()->FindBin(ptbin[ipt]+eps);
             int xbinmax = hvobs->GetXaxis()->FindBin(ptbin[ipt+1]-eps);
           //  std::cout<<xbinmin<<" "<<xbinmax<<std::endl;
           //  std::cout<<ptbin[ipt]<<" "<<ptbin[ipt+1]<<std::endl;
           // TH1F* hvobsProj = (TH1F*)hvobs->ProjectionY(Form("hvobsProj_%d",ipt),xbinmin,xbinmax);
           // TH1F* hvobssqProj = (TH1F*)hvobssq->ProjectionY(Form("hvobssqProj_%d",ipt),xbinmin,xbinmax);
            float Ntracks = hvobs2->Integral(xbinmin,xbinmax);
            float vobs = hvobs->Integral(xbinmin,xbinmax)/Ntracks;
            float vobssq = hvobssq->Integral(xbinmin,xbinmax)/Ntracks;
           // float vobssq = hvobssqProj->GetMean();
            float v = vobs/reso;
//            float verr = hvobsProj->GetRMS()/reso/sqrt(Ntracks);
            float verr = sqrt(vobssq-vobs*vobs)/reso/sqrt(Ntracks);
            TH1D* hvobsclone = (TH1D*)hvobs2->Clone("hvobsclone");
            hvobsclone->GetXaxis()->SetRangeUser(ptbin[ipt]+eps,ptbin[ipt+1]-eps);
            float pt = hvobsclone->GetMean();
            fout<<pt<<" "<<v<<" "<<verr<<" "<<std::endl;
            foutraw<<pt<<" "<<vobs<<" "<<verr*reso<<" "<<std::endl;
         }
        fout.close();
        foutraw.close();
         }
        fout1.close();
        fout2.close();
        }
       // }
     }

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



float GetReso(int iangle1,int  iangle2,int icent, int ihar, int isub, bool usingCNTEP){
 //   return 1;
    TString str1, str2;
    TFile *fin;
    int nrun = GetTotalRun();
    if(nrun<0) exit(1);
    int iharE=0;
    if(nhar==1) iharE=1;
    int n = ihar+1.0+iharE;
     ofstream fout;
     TH1F *hEPR1temp;
     TH1F *hEPR2temp;
     TH1F *hEPR1temp2;
     TH1F *hEPR2temp2;
     TH1F* hEPR1;
     TH1F* hEPR2;
     TH1F* hEPR12;
     TH1F* hEPR22;
    str1 = choosesub(isub);
    str2 = choosesub1(isub);
    if(str1=="ABORT" || str2 == "ABORT") return -9999;
         if(usingCNTEP){
        hEPR1 = new TH1F(Form("hEPR1_%d_%d_%d",icent,ihar,isub),Form("hEPR1_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
        hEPR2 = new TH1F(Form("hEPR2_%d_%d_%d",icent,ihar,isub),Form("hEPR2_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
         
        hEPR12 = new TH1F(Form("hEPR12_%d_%d_%d",icent,ihar,isub),Form("hEPR12_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
         hEPR22 = new TH1F(Form("hEPR22_%d_%d_%d",icent,ihar,isub),Form("hEPR22_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
         }
         else{
         hEPR1 = new TH1F(Form("hEPR1_%d_%d_%d",icent,ihar,isub),Form("hEPR1_%d_%d_%d",icent,ihar,isub),60,0,6);
         hEPR2 = new TH1F(Form("hEPR2_%d_%d_%d",icent,ihar,isub),Form("hEPR2_%d_%d_%d",icent,ihar,isub),60,0,6);
         
         hEPR12 = new TH1F(Form("hEPR12_%d_%d_%d",icent,ihar,isub),Form("hEPR12_%d_%d_%d",icent,ihar,isub),60,0,6);
         hEPR22 = new TH1F(Form("hEPR22_%d_%d_%d",icent,ihar,isub),Form("hEPR22_%d_%d_%d",icent,ihar,isub),60,0,6);
         }
        TH1F* hEPR3 = new TH1F(Form("hEPR3_%d_%d_%d",icent,ihar,isub),Form("hEPR3_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
        std::cout<<"Calculating Resolution..."<<std::endl;
        for(int irun=0;irun<nrun;irun++){
         //std::cout<<"cent = "<<icent<<"; n = "<<n<<" ;run = "<<irun<<" of total "<<nrun<<" runs"<<std::endl;
         //fin = TFile::Open(Form("/phenix/plhf/xuq/taxi/%s%s/%d/data/%s.root",dataset.Data(),pro.Data(),taxi,GetRun(irun).Data()));
         fin = TFile::Open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/work/200GeV/output/%s",GetRun(irun).Data()));
         if(!(GoodRunFit[icent][ihar][isub][irun]>0.2 && GoodRunFit[icent][ihar][isub][irun]<3.0)){
         std::cout<<"cent = "<<icent<<"; n = "<<n<<" ;isub = "<<str1<<" ;run = "<<GetRun(irun)<<" is bad run!"<<std::endl;
         fin->Close();
        continue;
         }
         if(usingCNTEP){
         hEPR1temp = (TH1F*)fin->Get(Form("EPR%s%s_%d_%d_%d_%d","CNT",str1.Data(),iangle1,iangle2,icent,ihar));
         hEPR2temp = (TH1F*)fin->Get(Form("EPR%s%s_%d_%d_%d_%d","CNT",str2.Data(),iangle1,iangle2,icent,ihar));
         hEPR1temp2 = new TH1F(Form("hEPRtemp12_%d_%d_%d",icent,ihar,isub),Form("hEPRtemp12_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
         hEPR2temp2 = new TH1F(Form("hEPRtemp22_%d_%d_%d",icent,ihar,isub),Form("hEPRtemp22_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
         }
         else{
            hEPR1temp = new TH1F(Form("hEPRtemp1_%d_%d_%d",icent,ihar,isub),Form("hEPRtemp1_%d_%d_%d",icent,ihar,isub),60,0,6);
            hEPR2temp = new TH1F(Form("hEPRtemp2_%d_%d_%d",icent,ihar,isub),Form("hEPRtemp2_%d_%d_%d",icent,ihar,isub),60,0,6);
            
            hEPR1temp2 = new TH1F(Form("hEPRtemp12_%d_%d_%d",icent,ihar,isub),Form("hEPRtemp12_%d_%d_%d",icent,ihar,isub),60,0,6);
            hEPR2temp2 = new TH1F(Form("hEPRtemp22_%d_%d_%d",icent,ihar,isub),Form("hEPRtemp22_%d_%d_%d",icent,ihar,isub),60,0,6);
             for(int iphi=0;iphi<nphi;iphi++){
         TProfile *hvobstemp1;
         TProfile *hvobstemp2;
         hvobstemp1 = (TProfile*)fin->Get(Form("vobs%s_%d_%d_%d_%d_%d",str1.Data(),iangle1,iangle2,icent,ihar,iphi));
         hvobstemp2 = (TProfile*)fin->Get(Form("vobs%s_%d_%d_%d_%d_%d",str2.Data(),iangle1,iangle2,icent,ihar,iphi));

         TH1D* hEPR1temptemp = hvobstemp1->ProjectionX(Form("hEPR1%s%s_%d_%d_%d","CNT",str1.Data(),icent,ihar,iphi),"W");
         TH1D* hEPR2temptemp = hvobstemp2->ProjectionX(Form("hEPR2%s%s_%d_%d_%d","CNT",str1.Data(),icent,ihar,iphi),"W");
         
         TH1D* hEPR1temptemp2 = hvobstemp1->ProjectionX(Form("hEPR12%s%s_%d_%d_%d","CNT",str1.Data(),icent,ihar,iphi),"B");
         TH1D* hEPR2temptemp2 = hvobstemp2->ProjectionX(Form("hEPR22%s%s_%d_%d_%d","CNT",str1.Data(),icent,ihar,iphi),"B");

            hEPR1temp->Add(hEPR1temptemp);
            hEPR2temp->Add(hEPR2temptemp);
            
            hEPR1temp2->Add(hEPR1temptemp2);
            hEPR2temp2->Add(hEPR2temptemp2);
            }
         }
         TH1F* hEPR3temp = (TH1F*)fin->Get(Form("EPR%s%s_%d_%d_%d_%d",str1.Data(),str2.Data(),iangle1,iangle2,icent,ihar));
         if(!hEPR3temp){
         hEPR3temp = (TH1F*)fin->Get(Form("EPR%s%s_%d_%d_%d_%d",str2.Data(),str1.Data(),iangle1,iangle2,icent,ihar));
         }
         hEPR1->Add(hEPR1temp);
         hEPR2->Add(hEPR2temp);
         if(!usingCNTEP){
         hEPR12->Add(hEPR1temp2);
         hEPR22->Add(hEPR2temp2);
         }
         hEPR3->Add(hEPR3temp);
         if(!usingCNTEP){
         delete hEPR1temp;
         delete hEPR2temp;
         }
         fin->Close();
        }
         if(usingCNTEP){
        if(hEPR1->GetMean()*hEPR3->GetMean()/hEPR2->GetMean()>0){
            float reso = sqrt(hEPR1->GetMean()*hEPR3->GetMean()/hEPR2->GetMean());
            return reso;
        }
        else return -9999;
         }
         else{
             int xbinmin = hEPR1->GetXaxis()->FindBin(0.4+eps);
             int xbinmax = hEPR1->GetXaxis()->FindBin(2.0-eps);
        if(hEPR1->Integral(xbinmin,xbinmax)/hEPR12->Integral(xbinmin,xbinmax)*hEPR3->GetMean()/(hEPR2->Integral(xbinmin,xbinmax)/hEPR22->Integral(xbinmin,xbinmax))>0){
            float reso = sqrt(hEPR1->Integral(xbinmin,xbinmax)/hEPR12->Integral(xbinmin,xbinmax)*hEPR3->GetMean()/(hEPR2->Integral(xbinmin,xbinmax)/hEPR22->Integral(xbinmin,xbinmax)));
            return reso;
        }
        else return -9999;
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
