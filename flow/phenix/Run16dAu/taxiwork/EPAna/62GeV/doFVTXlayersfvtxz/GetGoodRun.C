#include "RpPar.h"
int GetTotalRun();
TString pro = "Pro104";
int taxi1 = 8570;
int taxi2 = 8559;

void GetGoodRun(){
    TString str;
    int nrun = GetTotalRun();
    if(nrun<0) exit("Empty run list file!");

     ofstream fout, fout1, fout2;
     int iharE=0;
     if(nhar==1) iharE=1;
     for(int icent=0;icent<ncent;icent++){
      for(int ihar=0;ihar<nhar;ihar++){
       for(int isub=0;isub<nsub;isub++){
        int n = ihar+1.0+iharE;
        if(isub==2)
         str = "BBCS";
        else if(isub==6)
         str = "FVTX1S";
        else if(isub==7)
         str = "FVTX2S";
        else continue;
         fout2.open(Form("Result/psi%d_%d_%s.dat",n,icent,str.Data())); //using str as event plane detector
        for(int irun=0;irun<nrun;irun++){
         cout<<"cent = "<<icent<<"; n = "<<ihar+1.0+iharE<<" ;isub = "<<str<<" ;run = "<<irun<<" of total "<<nrun<<" runs"<<endl;
         fout2<<GetRun(irun)<<" "<<GoodRun(icent,ihar,isub,irun)<<endl;
         if(!(GoodRun(icent,ihar,isub,irun)>0.2 && GoodRun(icent,ihar,isub,irun)<3.0)){
         cout<<"cent = "<<icent<<"; n = "<<ihar+1.0+iharE<<" ;isub = "<<str<<" ;run = "<<GetRun(irun)<<" is bad run!"<<endl;
         }
         }
        fout2.close();
         }
        }
     }
}

float GoodRun(int icent, int ihar, int isub, int irun){
    float pi = acos(-1);
    TF1 *fun = new TF1("fun","pol0",-pi,pi);
    TString str="";
    TFile *fin;
    int iharE;
    if(nhar==1) iharE=1;
    int n = ihar+1.0+iharE;

     ofstream fout;
        if(isub==2){
         str = "BBCS";
        }
        else if(isub==6){
         str = "FVTX1S";
        }
        else if(isub==7){
         str = "FVTX2S";
        }
        else return -9999;
        fin = TFile::Open(Form("/phenix/plhf/xuq/taxi/%s%s/%d/data/%d.root",dataset.Data(),pro.Data(),taxi2,GetRun(irun)));
        TH1F* hpsi = new TH1F("psi","psi",100,-pi,pi);
        for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
          hpsitemp = (TH1F*)fin->Get(Form("psi_%d_%d_%d_%d",icent,ibbcz,ihar,isub));
          hpsi->Add(hpsitemp);
        }
      if(hpsi->GetEntries()>10000){
	hpsi->SetMarkerStyle(20);
	hpsi->SetMarkerSize(0.6);
	hpsi->SetMarkerColor(4);
	hpsi->SetMinimum(10);
	hpsi->Fit("fun","QR0");
	float par=fun->GetParameter(0);
	hpsi->SetMaximum(1.5*par);
        fin->Close();
	//hpsi->Draw();
	return fun->GetChisquare()/fun->GetNDF();
      }
      else{
        fin->Close();
        return -9999;
      }

    }


int GetTotalRun(){
    ifstream frun(dataset+".Lst");
    int RunNumber;
    int nrun = -1;
    while(!frun.eof()){
        frun>>RunNumber;
        nrun++;
    }
    return nrun;
}


int GetRun(int irun){
    ifstream frun(dataset+".Lst");
    int RunNumber;
    int jrun = 0;
    while(!frun.eof() && jrun<=irun){
        frun>>RunNumber;
        jrun++;
    }
    return RunNumber;
}
