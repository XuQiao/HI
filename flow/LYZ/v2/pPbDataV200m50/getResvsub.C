#include <iomanip>
#include <fstream>
#include "par.h"

void getResvsub(){
	double theta[ntheta];
	TFile *f[nFileAll];
		
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;


	ofstream fstrv;
	if(isSum)fstrv.open("v_sub.txt");
	else fstrv.open("v_2_sub.txt");
	TVectorD totmult[nbin], totpt[nbin];	TVectorD Nevent, totmultall;
	TVectorD V_int,V_mean;
	TVectorD* V_mean_;
	TVectorD* V[nbin]; TVectorD* chi[nbin];
	TVectorD* dDRe[nbin]; TVectorD* dDIm[nbin];
	TVectorD* dNRe[nbin][ntheta]; TVectorD* dNIm[nbin][ntheta];
	TComplex dD[nbin][ntheta], dN[nbin][ntheta][nptv];
	TVectorD avgmult[nbin];	TVectorD avgmultall;
	TVectorD avgpt[nbin];
	TVectorD v[nbin][ntheta],vmean[nbin],deltav[nbin][ntheta],deltavmean[nbin];
	TVectorD vmeanmean[nbin],deltavmeanmean[nbin],sigmavmeanmean[nbin];
	for(int ibin=0;ibin<nbin;ibin++){
		vmeanmean[ibin].ResizeTo(nptv);	vmeanmean[ibin].Zero();
		deltavmeanmean[ibin].ResizeTo(nptv);	deltavmeanmean[ibin].Zero();
		sigmavmeanmean[ibin].ResizeTo(nptv);	sigmavmeanmean[ibin].Zero();
	}
	if(isSum)	TFile *infile = TFile::Open("mergedV_Sum_sub.root");
	else	TFile *infile = TFile::Open("mergedV_Prod_sub.root");
        if(isSum)TFile *fout = new TFile("mergedv_Prod_sub.root","Recreate");
        else TFile *fout = new TFile("mergedv_Prod2_sub.root","Recreate");
        int subsample = (int)(nFileAll/nsamples);

	for(int isample=0;isample<nsamples;isample++){
        int start=isample*subsample;
        int end=(isample+1)*subsample;
	Nevent.ResizeTo(nbin);	Nevent.Zero();
	totmultall.ResizeTo(nbin);	totmultall.Zero();
	avgmultall.ResizeTo(nbin);
	V_mean.ResizeTo(nbin);	V_mean.Zero();
	V_int.ResizeTo(nbin);	V_int.Zero();
	for(int ibin=0;ibin<nbin;ibin++){
		totmult[ibin].ResizeTo(nptv);	totmult[ibin].Zero();
		avgmult[ibin].ResizeTo(nptv);	avgpt[ibin].ResizeTo(nptv);
		totpt[ibin].ResizeTo(nptv);	totpt[ibin].Zero();
		V[ibin] = (TVectorD*) infile->Get(Form("s_%d/D_%d/D_0/V",isample,ibin));
		chi[ibin] = (TVectorD*) infile->Get(Form("s_%d/D_%d/chi",isample,ibin));
		V_mean_ = (TVectorD*) infile->Get(Form("s_%d/D_%d/Vmean",isample,ibin));
		vmean[ibin].ResizeTo(nptv);	deltavmean[ibin].ResizeTo(nptv);
		vmean[ibin].Zero();	deltavmean[ibin].Zero();
		for(int itheta=0;itheta<ntheta;itheta++){
			v[ibin][itheta].ResizeTo(nptv);	deltav[ibin][itheta].ResizeTo(nptv);
			dD[ibin][itheta]=0;
			for(int iptbin=0;iptbin<nptv;iptbin++){
				dN[ibin][itheta][iptbin]=0;
			}
		}
	}

        for(int ifile=start; ifile<end; ifile++){
	        if(isSum) f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/pPbDataV205m185/Anav_Prod_sub_%d.root",ifile));
	        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/pPbDataV205m185/Anav_Prod2_sub_%d.root",ifile));
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
				v[ibin][itheta][iptbin]=(*V[ibin])[itheta]*avgmultall[ibin]*TMath::BesselJ1(j01)/Besselj01(mm)*Res.Re();
				vmean[ibin][iptbin]+=v[ibin][itheta][iptbin];
				deltav[ibin][itheta][iptbin]=TMath::Cos(mm*nn*theta[itheta])/totmult[ibin][iptbin]*(TMath::Exp(j01*j01/2./(*chi[ibin])[0]/(*chi[ibin])[0]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Power(-1,mm)*TMath::Exp(-j01*j01/2./(*chi[ibin])[0]/(*chi[ibin])[0]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.)));
				deltavmean[ibin][iptbin]+=deltav[ibin][itheta][iptbin];
			//	fstrv<<itheta<<"\t"<<v[ibin][itheta][iptbin]<<"\t"<<deltav[ibin][itheta][iptbin]<<endl;
			}
		deltavmean[ibin][iptbin]=TMath::Sqrt(deltavmean[ibin][iptbin])/2./Besselj01(mm);
		//fstrv<<endl;
		vmean[ibin][iptbin]/=ntheta;
		deltavmean[ibin][iptbin]/=TMath::Sqrt(ntheta);
		fstrv<<ptbinv[iptbin]<<"-"<<ptbinv[iptbin+1]<<"\t"<<vmean[ibin][iptbin]<<"\t"<<deltavmean[ibin][iptbin]<<endl;
		V_int[ibin]+=vmean[ibin][iptbin]*totmult[ibin][iptbin];
		vmeanmean[ibin][iptbin]+=vmean[ibin][iptbin]/nsamples;
		deltavmeanmean[ibin][iptbin]+=deltavmean[ibin][iptbin]/nsamples;
		sigmavmeanmean[ibin][iptbin]+=TMath::Power(vmean[ibin][iptbin]/nsamples,2);
		}
		V_int[ibin]/=totmultall[ibin];
		V_mean[ibin]=(*V_mean_)[ibin];
		fstrv<<endl<<"V ref="<<V_mean[ibin]<<"\t"<<"V int="<<V_int[ibin]<<endl;
	}
	
        TDirectory *dirsample = fout->mkdir(Form("s_%d",isample));
        dirsample->cd();
	for(ibin=0;ibin<nbin;ibin++){
	Nevent.Write("Nevent");
	V_int.Write("V_int");
	TDirectory *dir = dirsample->mkdir(Form("D_%d",ibin));
	dir->cd();
	avgpt[ibin].Write("avgpt");
	totmult[ibin].Write("totmult");
	vmean[ibin].Write("vmean");
	deltavmean[ibin].Write("deltavmean");
	}
	}	
	fout->cd();
	for(ibin=0;ibin<nbin;ibin++){
	TDirectory *dir_ = fout->mkdir(Form("D_%d",ibin));
	dir_->cd();
	fstrv<<endl;
		for(int iptbin=0;iptbin<nptv; iptbin++){
			sigmavmeanmean[ibin][iptbin]=TMath::Sqrt(sigmavmeanmean[ibin][iptbin]*nsamples-vmeanmean[ibin][iptbin]*vmeanmean[ibin][iptbin]);
			fstrv<<vmeanmean[ibin][iptbin]<<"\t"<<deltavmeanmean[ibin][iptbin]<<"\t"<<sigmavmeanmean[ibin][iptbin]<<endl;
		}
		avgpt[ibin].Write("avgpt");	vmeanmean[ibin].Write("vmeanmean");deltavmeanmean[ibin].Write("deltavmeanmean");sigmavmeanmean[ibin].Write("sigmavmeanmean");
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

