#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "RpPar.h"

using namespace std;

int GetTotalRun();
void FillGoodRun();
float GetReso(int , int , int, bool);
TString GetRun(int);
float eps = 1e-4;
float GoodRunFit[ncent][nhar][nsub][10000];

TString choosesub(int isub){
    TString str;
        if(isub==0)
         str = "FVTX1LS";
        else if(isub==1)
         str = "FVTX2LS";
        else if(isub==2)
         str = "FVTX3LS";
        else if(isub==3)
         str = "FVTX4LS";
        else if(isub==5)
         str = "FVTX1S";
        else if(isub==4)
         str = "BBCS";
        else if(isub==6)
         str = "FVTX1p2p3LS";
        else if(isub==7)
         str = "FVTX1p2p4LS";
        else if(isub==8)
         str = "FVTX1N";
        else if(isub==9)
         str = "FVTXtrkS";
        else
         str = "ABORT";
    return str;
}

TString choosesub1(int isub){
    TString str;
    if(choosesub(isub).Contains("BBC") || choosesub(isub).Contains("N")) //BBC or FVTX1N
        str = "FVTX1S";
    else if(choosesub(isub).Contains("FVTX"))
        str = "BBCS";
    else 
        str = "ABORT";
        return str;
}


void Getvn2D(int icent, int ihar, bool usingCNTEP=0){
    TString str;
    TFile *fin;
    int nrun = GetTotalRun();
    std::cout<<"Totally we have "<<nrun<<" runs/segments!"<<std::endl;
    FillGoodRun();
    std::cout<<"Filling Good run finished!"<<std::endl;

    if(nrun<0) exit(1);

     ofstream fout, foutraw, fout1, fout2;
     int iharE=0;
     if(nhar==1||nhar==2) iharE=1;
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
        cout<<Form("Result/%s/res%d_%d_%s.dat",UseCNTEP.Data(),n,icent,str.Data())<<endl;
         fout1.open(Form("Result/%s/res%d_%d_%s.dat",UseCNTEP.Data(),n,icent,str.Data())); //using str as event plane detector
         fout2.open(Form("Result/%s/psi%d_%d_%s.dat",UseCNTEP.Data(),n,icent,str.Data())); //using str as event plane detector
         float reso = GetReso(icent,ihar,isub,usingCNTEP);
         fout1<<reso<<std::endl;
         if(reso<0) {std::cout<<"resolution is wrong!"<<std::endl; reso = -9999;}
        for(int irun=0;irun<nrun;irun++){
         fout2<<GetRun(irun)<<" "<<GoodRunFit[icent][ihar][isub][irun]<<std::endl;
        }
         TH2F* hvobsall = new TH2F(Form("hvobs_%d_%d_%d",icent,ihar,isub),Form("hvobs_%d_%d_%d",icent,ihar,isub),60,0,6,220,-1.1,1.1);
        // TH2F* hvobsallsq = new TH2F(Form("hvobssq_%d_%d_%d",icent,ihar,isub),Form("hvobssq_%d_%d_%d",icent,ihar,isub),60,0,6,220,-1.1,1.1);
         
        for(int iphi=0;iphi<nphi+1;iphi++){
         TH2F* hvobs = new TH2F(Form("hvobs_%d_%d_%d_%d",icent,ihar,isub,iphi),Form("hvobs_%d_%d_%d_%d",icent,ihar,isub,iphi),60,0,6,220,-1.1,1.1);
        // TH2F* hvobssq = new TH2F(Form("hvobssq_%d_%d_%d_%d",icent,ihar,isub,iphi),Form("hvobssq_%d_%d_%d_%d",icent,ihar,isub,iphi),60,0,6,220,-1.1,1.1);
         string phistr = (iphi==0)?"_east":"_west";
         if(iphi==nphi) phistr = "";
         cout<<"open v2 file"<<endl;
         fout.open(Form("Result/%s/v%d_%d%s_%s.dat",UseCNTEP.Data(),n,icent,phistr.c_str(),str.Data())); //using str as event plane detector
         cout<<"open v2raw file"<<endl;
         foutraw.open(Form("Result/%s/v%draw_%d%s_%s.dat",UseCNTEP.Data(),n,icent,phistr.c_str(),str.Data())); //using str as event plane detector
         if(iphi<nphi){
        for(int irun=0;irun<nrun;irun++){
         //std::cout<<"cent = "<<icent<<"; n = "<<n<<" ;isub = "<<str<<" ;run = "<<irun<<" "<<phistr<<std::endl;
         fin = TFile::Open(Form("/store/user/qixu/flow/Run16dAu/62GeV/%s",GetRun(irun).Data()));
         if(!(GoodRunFit[icent][ihar][isub][irun]>0.2 && GoodRunFit[icent][ihar][isub][irun]<3.0)){
         std::cout<<"cent = "<<icent<<"; n = "<<n<<" ;isub = "<<str<<" ;run = "<<GetRun(irun)<<" is bad run!"<<std::endl;
         fin->Close();
        continue;
         }
         TH2F* hvobstemp = (TH2F*)fin->Get(Form("vobs%s_%d_%d_%d",str.Data(),icent,ihar,iphi));
         //TH2F* hvobssqtemp = (TH2F*)fin->Get(Form("vobs%ssq_%d_%d_%d",str.Data(),icent,ihar,iphi));
         hvobs->Add(hvobstemp);
         //hvobssq->Add(hvobssqtemp);
         fin->Close();
        }
         }
        hvobsall->Add(hvobs);
        //hvobsallsq->Add(hvobssq);
        if(iphi==nphi){
        hvobs = hvobsall;
        //hvobssq = hvobsallsq;
        }
            TH1F* ptProj = (TH1F*)hvobs->ProjectionX(Form("hptProj"),0,-1);
         for(int ipt=0;ipt<npt-1;ipt++){
             int xbinmin = hvobs->GetXaxis()->FindBin(ptbin[ipt]+eps);
             int xbinmax = hvobs->GetXaxis()->FindBin(ptbin[ipt+1]-eps);
           //  std::cout<<xbinmin<<" "<<xbinmax<<std::endl;
           //  std::cout<<ptbin[ipt]<<" "<<ptbin[ipt+1]<<std::endl;
            TH1F* hvobsProj = (TH1F*)hvobs->ProjectionY(Form("hvobsProj_%d",ipt),xbinmin,xbinmax);
           // TH1F* hvobssqProj = (TH1F*)hvobssq->ProjectionY(Form("hvobssqProj_%d",ipt),xbinmin,xbinmax);
            float vobs = hvobsProj->GetMean();
            float Ntracks = hvobsProj->Integral();
           // float vobssq = hvobssqProj->GetMean();
            float v = vobs/reso;
            float verr=-9999;
            if(Ntracks>0)
            verr = hvobsProj->GetRMS()/reso/sqrt(Ntracks);
           // float verr = sqrt(vobssq/reso/reso-(v*v))/sqrt(Ntracks);
            ptProj->GetXaxis()->SetRangeUser(ptbin[ipt],ptbin[ipt+1]);
            float pt = ptProj->GetMean();
            fout<<pt<<" "<<v<<" "<<verr<<" "<<std::endl;
            foutraw<<pt<<" "<<vobs<<" "<<verr*reso<<" "<<std::endl;
         }
        fout.close();
        foutraw.close();
         }
        fout1.close();
        fout2.close();
        }
}

void FillGoodRun(){
    float pi = acos(-1.0);
    TString str;
    TFile *fin;
    int nrun = GetTotalRun();
    if(nrun<0) exit(1);
     for(int icent=0;icent<ncent;icent++){
      for(int ihar=0;ihar<nhar;ihar++){
        //  if(ihar!=0) continue;
       for(int isub=0;isub<nsub;isub++){
            str = choosesub(isub);
            if(str=="ABORT") continue;
            for(int irun=0;irun<nrun;irun++){
        //std::cout<<icent<<" "<<ihar<<" "<<isub<<" "<<irun<<std::endl;
         fin = TFile::Open(Form("/store/user/qixu/flow/Run16dAu/62GeV/%s",GetRun(irun).Data()));
        TH1F* hpsi = new TH1F("psi","psi",100,-pi,pi);
        for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
          TH1F* hpsitemp = (TH1F*)fin->Get(Form("psiFla_%d_%d_%d_%d",icent,ibbcz,ihar,isub));
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
    }
       }
      }
     }
}



float GetReso(int icent, int ihar, int isub, bool usingCNTEP){
    TString str1, str2;
    TFile *fin;
    int nrun = GetTotalRun();
    if(nrun<0) exit(1);
    int iharE=0;
    if(nhar==1||nhar==2) iharE=1;
    int n = ihar+1.0+iharE;
     ofstream fout;
     TH1F *hEPR1temp;
     TH1F *hEPR2temp;
    str1 = choosesub(isub);
    str2 = choosesub1(isub);
    if(str1=="ABORT" || str2 == "ABORT") return -9999;
        TH1F* hEPR1 = new TH1F(Form("hEPR1_%d_%d_%d",icent,ihar,isub),Form("hEPR1_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
        TH1F* hEPR2 = new TH1F(Form("hEPR2_%d_%d_%d",icent,ihar,isub),Form("hEPR2_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
        TH1F* hEPR3 = new TH1F(Form("hEPR3_%d_%d_%d",icent,ihar,isub),Form("hEPR3_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
        std::cout<<"Calculating Resolution..."<<std::endl;
        for(int irun=0;irun<nrun;irun++){
         //std::cout<<"cent = "<<icent<<"; n = "<<n<<" ;run = "<<irun<<" of total "<<nrun<<" runs"<<std::endl;
         fin = TFile::Open(Form("/store/user/qixu/flow/Run16dAu/62GeV/%s",GetRun(irun).Data()));
         if(!(GoodRunFit[icent][ihar][isub][irun]>0.2 && GoodRunFit[icent][ihar][isub][irun]<3.0)){
         std::cout<<"cent = "<<icent<<"; n = "<<n<<" ;isub = "<<str1<<" ;run = "<<GetRun(irun)<<" is bad run!"<<std::endl;
         fin->Close();
        continue;
         }
         if(usingCNTEP){
            if(isub>=9 && isub<12) //do not store CNT-FVTX1-4LS correlation histograms, only use NOCNTEP
                return 0.2;
         hEPR1temp = (TH1F*)fin->Get(Form("EPR%s%s_%d_%d","CNT",str1.Data(),icent,ihar));
         hEPR2temp = (TH1F*)fin->Get(Form("EPR%s%s_%d_%d","CNT",str2.Data(),icent,ihar));
         }
         else{
            hEPR1temp = new TH1F(Form("hEPRtemp1_%d_%d_%d",icent,ihar,isub),Form("hEPRtemp1_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
            hEPR2temp = new TH1F(Form("hEPRtemp2_%d_%d_%d",icent,ihar,isub),Form("hEPRtemp2_%d_%d_%d",icent,ihar,isub),220,-1.1,1.1);
             for(int iphi=0;iphi<nphi;iphi++){
         TH2F* hvobstemp1 = (TH2F*)fin->Get(Form("vobs%s_%d_%d_%d",str1.Data(),icent,ihar,iphi));
         TH2F* hvobstemp2 = (TH2F*)fin->Get(Form("vobs%s_%d_%d_%d",str2.Data(),icent,ihar,iphi));
         TH1D* hEPR1temptemp = hvobstemp1->ProjectionY(Form("hEPR1%s%s_%d_%d","CNT",str1.Data(),icent,ihar),hvobstemp1->GetXaxis()->FindBin(0.2+1e-4),hvobstemp1->GetXaxis()->FindBin(5.0-1e-4));
         TH1D* hEPR2temptemp = hvobstemp2->ProjectionY(Form("hEPR2%s%s_%d_%d","CNT",str1.Data(),icent,ihar),hvobstemp1->GetXaxis()->FindBin(0.2+1e-4),hvobstemp1->GetXaxis()->FindBin(5.0-1e-4));
            hEPR1temp->Add(hEPR1temptemp);
            hEPR2temp->Add(hEPR2temptemp);
            }
         }
         TH1F* hEPR3temp = (TH1F*)fin->Get(Form("EPR%s%s_%d_%d",str1.Data(),str2.Data(),icent,ihar));
         if(!hEPR3temp)
         hEPR3temp = (TH1F*)fin->Get(Form("EPR%s%s_%d_%d",str2.Data(),str1.Data(),icent,ihar));

         hEPR1->Add(hEPR1temp);
         hEPR2->Add(hEPR2temp);
         hEPR3->Add(hEPR3temp);
         if(!usingCNTEP){
         delete hEPR1temp;
         delete hEPR2temp;
         }
         fin->Close();
        }
        cout<<"CNT-"<<str1<<" "<<hEPR1->GetMean()<<endl;
        cout<<str2<<"-"<<str1<<" "<<hEPR3->GetMean()<<endl;
        cout<<"CNT-"<<str2<<" "<<hEPR2->GetMean()<<endl;
        if(hEPR1->GetMean()*hEPR3->GetMean()/hEPR2->GetMean()>0){
            float reso = sqrt(hEPR1->GetMean()*hEPR3->GetMean()/hEPR2->GetMean());
            return reso;
        }
        else return -9999;
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
