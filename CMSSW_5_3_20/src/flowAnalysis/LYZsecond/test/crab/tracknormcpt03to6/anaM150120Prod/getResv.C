#include <iomanip>
#include <fstream>
#include "par.h"

void getResv(){

	double theta[ntheta];
	int nFileAll=getnumber("filelist.dat");
	TFile *f;
		
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;


	ofstream fstrv;
	fstrv.open("v.txt");
	TVectorD totmult[nbin], totpt[nbin];	TVectorD Nevent, totmultall;
	TVectorD* V[nbin]; TVectorD* chi[nbin];
	TVectorD* dDRe[nbin]; TVectorD* dDIm[nbin];
	TVectorD* dNRe[nbin][ntheta]; TVectorD* dNIm[nbin][ntheta];
	TComplex dD[nbin][ntheta], dN[nbin][ntheta][npt];
	TVectorD avgmult[nbin];	TVectorD avgmultall;
	TVectorD avgpt[nbin];
	TVectorD v[nbin][ntheta],vmean[nbin],deltav[nbin][ntheta],deltavmean[nbin];
	TFile *infile = TFile::Open("mergedV_Prod_M150120.root");
	Nevent.ResizeTo(nbin);	Nevent.Zero();
	totmultall.ResizeTo(nbin);	totmultall.Zero();
	avgmultall.ResizeTo(nbin);
	for(int ibin=0;ibin<nbin;ibin++){
		totmult[ibin].ResizeTo(npt);	totmult[ibin].Zero();
		avgmult[ibin].ResizeTo(npt);	avgpt[ibin].ResizeTo(npt);
		totpt[ibin].ResizeTo(npt);	totpt[ibin].Zero();
		vmean[ibin].ResizeTo(npt);	deltavmean[ibin].ResizeTo(npt);
		V[ibin] = (TVectorD*) infile->Get(Form("D_%d/V",ibin));
		chi[ibin] = (TVectorD*) infile->Get(Form("D_%d/chi",ibin));
		for(int itheta=0;itheta<ntheta;itheta++){
			v[ibin][itheta].ResizeTo(npt);	deltav[ibin][itheta].ResizeTo(npt);
			dD[ibin][itheta]=0;
			for(int iptbin=0;iptbin<npt;iptbin++){
				dN[ibin][itheta][iptbin]=0;
			}
		}
	}

        for(int ifile=0; ifile<nFileAll; ifile++){
		f=TFile::Open(readline("filelist.dat",ifile).c_str());
		TVectorD* Nevent_t = (TVectorD*)f->Get("demo/Nevent");	
		TVectorD* totmultall_t = (TVectorD*)f->Get("demo/totmultall");
		for(int ibin=0;ibin<nbin;ibin++){
			TVectorD* totmult_t = (TVectorD*)f->Get(Form("demo/totmult_%d",ibin));
			TVectorD* totpt_t = (TVectorD*)f->Get(Form("demo/totpt_%d",ibin));
			dDRe[ibin] = (TVectorD*)f->Get(Form("demo/dDRe_%d",ibin));
			dDIm[ibin] = (TVectorD*)f->Get(Form("demo/dDIm_%d",ibin));
			Nevent[ibin]+=(*Nevent_t)[ibin];
			totmultall[ibin]+=(*totmultall_t)[ibin];
			for(int iptbin=0;iptbin<npt;iptbin++){
				totmult[ibin][iptbin]+=(*totmult_t)[iptbin];
				totpt[ibin][iptbin]+=(*totpt_t)[iptbin];
			}
			for(int itheta=0;itheta<ntheta;itheta++){
				dNRe[ibin][itheta] = (TVectorD*)f->Get(Form("demo/dNRe_%d_%d",ibin,itheta));
				dNIm[ibin][itheta] = (TVectorD*)f->Get(Form("demo/dNIm_%d_%d",ibin,itheta));
				dD[ibin][itheta]+=TComplex((*dDRe[ibin])[itheta],(*dDIm[ibin])[itheta]);
				for(int iptbin=0;iptbin<npt;iptbin++)
					dN[ibin][itheta][iptbin]+=TComplex((*dNRe[ibin][itheta])[iptbin],(*dNIm[ibin][itheta])[iptbin]);
			}
		}
	f->Close();
	}

	for(int ibin=0;ibin<nbin;ibin++)
		for(int itheta=0;itheta<ntheta;itheta++)
			dD[ibin][itheta]*=TComplex::Power(TComplex::I(),mm-1)/(Double_t)Nevent[ibin];
	
	for(int ibin=0;ibin<nbin;ibin++){
		avgmultall[ibin]=1.0*totmultall[ibin]/Nevent[ibin];
		fstrv<<trkbin[ibin]<<"-"<<trkbin[ibin+1]<<endl;
		for(int iptbin=0;iptbin<npt; iptbin++){
			vmean[ibin][iptbin]=0;
			deltavmean[ibin][iptbin]=0;
			avgmult[ibin][iptbin]=1.0*totmult[ibin][iptbin]/Nevent[ibin];
			avgpt[ibin][iptbin]=1.0*totpt[ibin][iptbin]/totmult[ibin][iptbin];
			for(int itheta=0;itheta<ntheta;itheta++){
				dN[ibin][itheta][iptbin]/=totmult[ibin][iptbin];
				TComplex Res=dN[ibin][itheta][iptbin]/dD[ibin][itheta];
				v[ibin][itheta][iptbin]=(*V[ibin])[itheta]*avgmultall[ibin]*TMath::BesselJ1(j01)/Besselj01(mm)*Res.Re();
				vmean[ibin][iptbin]+=v[ibin][itheta][iptbin];
				deltav[ibin][itheta][iptbin]=TMath::Cos(mm*nn*theta[itheta])/totmult[ibin][iptbin]*(TMath::Exp(j01*j01/2./(*chi[ibin])[itheta]/(*chi[ibin])[itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Power(-1,mm)*TMath::Exp(-j01*j01/2./(*chi[ibin])[itheta]/(*chi[ibin])[itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.)));
				deltavmean[ibin][iptbin]+=deltav[ibin][itheta][iptbin];
				//fstrv<<itheta<<"\t"<<v[ibin][itheta][iptbin]<<"\t"<<deltav[ibin][itheta][iptbin]<<endl;
			}
			deltavmean[ibin][iptbin]=TMath::Sqrt(deltavmean[ibin][iptbin])/2./Besselj01(mm);
			vmean[ibin][iptbin]/=ntheta;
			deltavmean[ibin][iptbin]/=ntheta;
			fstrv<<ptbin[iptbin]<<"-"<<ptbin[iptbin+1]<<"\t"<<vmean[ibin][iptbin]<<"\t"<<deltavmean[ibin][iptbin]<<endl;
		}
		fstrv<<endl;
	}
	
	TFile *fout = new TFile("mergedv_Prod.root","Recreate");
	for(ibin=0;ibin<nbin;ibin++){
	Nevent.Write("Nevent");
	TDirectory *dir = fout->mkdir(Form("D_%d",ibin));
	dir->cd();
	avgpt[ibin].Write("avgpt");
	vmean[ibin].Write("vmean");
	deltavmean[ibin].Write("deltavmean");
	}
	
	infile->Close();
	fout->Close();
	
}

double Besselj01(int m0){
        if(m0==1)
                return TMath::BesselJ1(j01);
        else if(m0==2)
                return 0.431755;
        else if(m0==3)
                return 0.199000;
        else if(m0==4)
                return 0.064747;
        else return 0;
}

