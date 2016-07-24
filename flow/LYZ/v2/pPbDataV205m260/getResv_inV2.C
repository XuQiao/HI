#include <iomanip>
#include <fstream>
#include "par.h"

void getResv_inV2(){

	double theta[ntheta];
	TFile *f[nFileAll];
		
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;


	ofstream fstrv;
	if(isSum)fstrv.open("v_inV2.txt");
	else fstrv.open("v_inV2_2.txt");
	TVectorD totmult[nbin], totpt[nbin];	TVectorD Nevent, totmultall;
	TVectorD V_int;
	TVectorD* V_mean;
	TVectorD* V[nbin]; TVectorD* chi[nbin];
	TVectorD* dDRe[nbin]; TVectorD* dDIm[nbin];
	TVectorD* dNRe[nbin][ntheta]; TVectorD* dNIm[nbin][ntheta];
	TComplex dD[nbin][ntheta], dN[nbin][ntheta][nptv];
	TVectorD avgmult[nbin];	TVectorD avgmultall;
	TVectorD avgpt[nbin];
	TVectorD v[nbin][ntheta],vmean[nbin],deltav[nbin][ntheta],deltavmean[nbin];
	if(isSum)	TFile *infile = TFile::Open("mergedV_Sum.root");
	else	TFile *infile = TFile::Open("mergedV_Prod.root");
	Nevent.ResizeTo(nbin);	Nevent.Zero();
	totmultall.ResizeTo(nbin);	totmultall.Zero();
	avgmultall.ResizeTo(nbin);
	V_int.ResizeTo(nbin);	V_int.Zero();
	for(int ibin=0;ibin<nbin;ibin++){
		totmult[ibin].ResizeTo(nptv);	totmult[ibin].Zero();
		avgmult[ibin].ResizeTo(nptv);	avgpt[ibin].ResizeTo(nptv);
		totpt[ibin].ResizeTo(nptv);	totpt[ibin].Zero();
		vmean[ibin].ResizeTo(nptv);	deltavmean[ibin].ResizeTo(nptv);
		V[ibin] = (TVectorD*) infile->Get(Form("D_%d/D_0/V",ibin));
		chi[ibin] = (TVectorD*) infile->Get(Form("D_%d/chi",ibin));
		V_mean = (TVectorD*) infile->Get(Form("D_%d/Vmean",ibin));
		for(int itheta=0;itheta<ntheta;itheta++){
			v[ibin][itheta].ResizeTo(nptv);	deltav[ibin][itheta].ResizeTo(nptv);
			dD[ibin][itheta]=0;
			for(int iptbin=0;iptbin<nptv;iptbin++){
				dN[ibin][itheta][iptbin]=0;
			}
		}
	}

        for(int ifile=0; ifile<nFileAll; ifile++){
	        if(isSum) f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/pPbDataV205m260/Anav_Prod_inV2_%d.root",ifile));
	        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/pPbDataV205m260/Anav_Prod2_inV2_%d.root",ifile));
		TVectorD* Nevent_t = (TVectorD*)f[ifile]->Get("Nevent");	
		TVectorD* totmultall_t = (TVectorD*)f[ifile]->Get("totmultall");
		for(int ibin=0;ibin<nbin;ibin++){
			TVectorD* totmult_t = (TVectorD*)f[ifile]->Get(Form("totmult_%d",ibin));
			TVectorD* totpt_t = (TVectorD*)f[ifile]->Get(Form("totpt_%d",ibin));
			dDRe[ibin] = (TVectorD*)f[ifile]->Get(Form("dDRe_%d",ibin));
			dDIm[ibin] = (TVectorD*)f[ifile]->Get(Form("dDIm_%d",ibin));
			Nevent[ibin]+=(*Nevent_t)[ibin];
			totmultall[ibin]+=(*totmultall_t)[ibin];
			for(int iptbin=0;iptbin<nptv;iptbin++){
				totmult[ibin][iptbin]+=(*totmult_t)[iptbin];
				totpt[ibin][iptbin]+=(*totpt_t)[iptbin];
			}
			for(int itheta=0;itheta<ntheta;itheta++){
				dNRe[ibin][itheta] = (TVectorD*)f[ifile]->Get(Form("dNRe_%d_%d",ibin,itheta));
				dNIm[ibin][itheta] = (TVectorD*)f[ifile]->Get(Form("dNIm_%d_%d",ibin,itheta));
				dD[ibin][itheta]+=TComplex((*dDRe[ibin])[itheta],(*dDIm[ibin])[itheta]);
				for(int iptbin=0;iptbin<nptv;iptbin++)
					dN[ibin][itheta][iptbin]+=TComplex((*dNRe[ibin][itheta])[iptbin],(*dNIm[ibin][itheta])[iptbin]);
			}
		}
	f[ifile]->Close();
	}

	for(int ibin=0;ibin<nbin;ibin++)
		for(int itheta=0;itheta<ntheta;itheta++)
			dD[ibin][itheta]*=TComplex::Power(TComplex::I(),mm-1)/(Double_t)Nevent[ibin];
	
	for(int ibin=0;ibin<nbin;ibin++){
			avgmultall[ibin]=1.0*totmultall[ibin]/Nevent[ibin];
		for(int iptbin=0;iptbin<nptv; iptbin++){
			vmean[ibin][iptbin]=0;
			deltavmean[ibin][iptbin]=0;
			avgmult[ibin][iptbin]=1.0*totmult[ibin][iptbin]/Nevent[ibin];
			avgpt[ibin][iptbin]=1.0*totpt[ibin][iptbin]/totmult[ibin][iptbin];
			for(int itheta=0;itheta<ntheta;itheta++){
				dN[ibin][itheta][iptbin]/=totmult[ibin][iptbin];
				TComplex Res=dN[ibin][itheta][iptbin]/dD[ibin][itheta];
				v[ibin][itheta][iptbin]=inV2*avgmultall[ibin]*TMath::BesselJ1(j01)/Besselj01(mm)*Res.Re();
				vmean[ibin][iptbin]+=v[ibin][itheta][iptbin];
				deltav[ibin][itheta][iptbin]=TMath::Cos(mm*nn*theta[itheta])/totmult[ibin][iptbin]*(TMath::Exp(j01*j01/2./(*chi[ibin])[0]/(*chi[ibin])[0]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Power(-1,mm)*TMath::Exp(-j01*j01/2./(*chi[ibin])[0]/(*chi[ibin])[0]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.)));
				deltavmean[ibin][iptbin]+=deltav[ibin][itheta][iptbin];
			//	fstrv<<itheta<<"\t"<<v[ibin][itheta][iptbin]<<"\t"<<deltav[ibin][itheta][iptbin]<<endl;
			}
		deltavmean[ibin][iptbin]=TMath::Sqrt(deltavmean[ibin][iptbin])/2./Besselj01(mm);
		//fstrv<<endl;
		vmean[ibin][iptbin]/=ntheta;
		deltavmean[ibin][iptbin]/=ntheta;
		fstrv<<ptbinv[iptbin]<<"-"<<ptbinv[iptbin+1]<<"\t"<<vmean[ibin][iptbin]<<"\t"<<deltavmean[ibin][iptbin]<<endl;
		V_int[ibin]+=vmean[ibin][iptbin]*totmult[ibin][iptbin];
		}
		V_int[ibin]/=totmultall[ibin];
		fstrv<<endl<<"V ref="<<inV2<<"\t"<<"V int="<<V_int[ibin]<<endl;
	}
	
	if(isSum)TFile *fout = new TFile("mergedv_Prod_inV2.root","Recreate");
	else TFile *fout = new TFile("mergedv_Prod2_inV2.root","Recreate");
	for(ibin=0;ibin<nbin;ibin++){
	Nevent.Write("Nevent");
	V_int.Write("V_int");
	TDirectory *dir = fout->mkdir(Form("D_%d",ibin));
	dir->cd();
	avgpt[ibin].Write("avgpt");
	totmult[ibin].Write("totmult");
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

