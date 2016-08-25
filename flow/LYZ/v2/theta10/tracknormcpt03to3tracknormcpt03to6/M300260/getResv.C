#include <iomanip>
#include <fstream>
#include "par.h"

void getResv(){

	double theta[ntheta];
	TFile *f[nFileAll];
		
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;


	string SumorProd = getenv("SUMORPROD");
	string mdir = getenv("DIR");
	ofstream fstrv;
	if(SumorProd=="Sum")fstrv.open("v.txt");
	else fstrv.open("v_2.txt");
	TVectorD totmult[nbin], totpt[nbin];	TVectorD Nevent, totmultall, totmultall_;
	TVectorD totmulthisto[nbin], totmulthistocorr[nbin];
	TVectorD V_int, V_interr;
	TVectorD V_intcorr, V_intcorrerr;
	TVectorD* V_mean;
	TVectorD* deltaV_mean;
	TVectorD* V[nbin]; TVectorD* chi[nbin];
	TVectorD* dDRe[nbin]; TVectorD* dDIm[nbin];
	TVectorD* dNRe[nbin][ntheta]; TVectorD* dNIm[nbin][ntheta];
	TComplex dD[nbin][ntheta], dN[nbin][ntheta][nptv];
	TVectorD avgmult[nbin];	TVectorD avgmultall;
	TVectorD avgpt[nbin];
	TVectorD v[nbin][ntheta],vmean[nbin],deltav[nbin][ntheta],deltavmean[nbin];
	if(SumorProd=="Sum")	TFile *infile = TFile::Open("mergedV_Sum.root");
	else	TFile *infile = TFile::Open("mergedV_Prod.root");
	Nevent.ResizeTo(nbin);	Nevent.Zero();
	totmultall.ResizeTo(nbin);	totmultall.Zero();
	totmultall_.ResizeTo(nbin);	totmultall_.Zero();
	avgmultall.ResizeTo(nbin);
	V_int.ResizeTo(nbin);	V_int.Zero();
	V_interr.ResizeTo(nbin);	V_interr.Zero();
	V_intcorr.ResizeTo(nbin);	V_intcorr.Zero();
	V_intcorrerr.ResizeTo(nbin);	V_intcorrerr.Zero();
	for(int ibin=0;ibin<nbin;ibin++){
		totmult[ibin].ResizeTo(nptv);	totmult[ibin].Zero();
		totmulthisto[ibin].ResizeTo(nptv);	
		totmulthistocorr[ibin].ResizeTo(nptv);
		avgmult[ibin].ResizeTo(nptv);	avgpt[ibin].ResizeTo(nptv);
		totpt[ibin].ResizeTo(nptv);	totpt[ibin].Zero();
		vmean[ibin].ResizeTo(nptv);	deltavmean[ibin].ResizeTo(nptv);
		V[ibin] = (TVectorD*) infile->Get(Form("D_%d/D_0/V",ibin));
		chi[ibin] = (TVectorD*) infile->Get(Form("D_%d/chi",ibin));
		V_mean = (TVectorD*) infile->Get(Form("D_%d/Vmean",ibin));
		deltaV_mean = (TVectorD*) infile->Get(Form("D_%d/deltaVmean",ibin));
		for(int itheta=0;itheta<ntheta;itheta++){
			v[ibin][itheta].ResizeTo(nptv);	deltav[ibin][itheta].ResizeTo(nptv);
			dD[ibin][itheta]=0;
			for(int iptbin=0;iptbin<nptv;iptbin++){
				dN[ibin][itheta][iptbin]=0;
			}
		}
	}

        for(int ifile=0; ifile<nFileAll; ifile++){
	        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/theta10/tracknormcpt03to3tracknormcpt03to6/%s/Anav_Prod_%d.root",mdir.c_str(),ifile));
	        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/theta10/tracknormcpt03to3tracknormcpt03to6/%s/Anav_Prod2_%d.root",mdir.c_str(),ifile));
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
	
	TH1D* hpt[nbin];
	TH1D* hpteffcorr[nbin];
	TFile *fhisto = TFile::Open("histomerged.root");
	TFile *feff = TFile::Open("/home/xuq7/HI/flow/LYZ/v2/TrackCorrections_HIJING_538_OFFICIAL_Mar24.root");
	TH2F* heff = (TH2F*)feff->Get("rTotalEff3D");
      	TH1D* hpteff = (TH1D*)heff->ProjectionY("hpteff",heff->GetXaxis()->FindBin(-2.4),heff->GetXaxis()->FindBin(2.4)-1,"o");
	TArrayD *ptarr = (TArrayD*)heff->GetYaxis()->GetXbins();
	double *ptbinhisto = ptarr->GetArray();
	int NbinX = heff->GetXaxis()->GetNbins();
	int NbinY = heff->GetYaxis()->GetNbins();
	hpteff->Scale(1.0/NbinX);
        
	for(int ibin=0;ibin<nbin;ibin++)
		for(int itheta=0;itheta<ntheta;itheta++)
			dD[ibin][itheta]*=TComplex::Power(TComplex::I(),mm-1)/(Double_t)Nevent[ibin];
	
	for(int ibin=0;ibin<nbin;ibin++){
		avgmultall[ibin]=1.0*totmultall[ibin]/Nevent[ibin];
                hpt[ibin] = (TH1D*)fhisto->Get(Form("D_%d/hpt",ibin));
                TH1D* hptre = (TH1D*)hpt[ibin]->Rebin(NbinY,"hptre",ptbinhisto);
                hpteffcorr[ibin] = (TH1D*)hptre->Clone(Form("hpteffcorr_%d",ibin));
                hpteffcorr[ibin]->Divide(hpteff);
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
		if(ptbinv[iptbin+1]>3.0) continue;
                totmulthisto[ibin][iptbin]=hpt[ibin]->Integral(hpt[ibin]->GetXaxis()->FindBin(ptbinv[iptbin]),hpt[ibin]->GetXaxis()->FindBin(ptbinv[iptbin+1])-1);
                totmulthistocorr[ibin][iptbin]=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[iptbin]),hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[iptbin+1])-1);
		V_int[ibin]+=vmean[ibin][iptbin]*totmult[ibin][iptbin];
		V_interr[ibin]+=deltavmean[ibin][iptbin]*totmult[ibin][iptbin];
		V_intcorr[ibin]+=vmean[ibin][iptbin]*totmulthistocorr[ibin][iptbin];
		V_intcorrerr[ibin]+=deltavmean[ibin][iptbin]*totmulthistocorr[ibin][iptbin];
		totmultall_[ibin]+=totmult[ibin][iptbin];
		}
		V_int[ibin]/=totmultall_[ibin];
		V_interr[ibin]/=totmultall_[ibin];
		V_intcorr[ibin]/=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1);
		V_intcorrerr[ibin]/=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1);
	}

        fstrv<<endl<<"pt range\t\t"<<"totmult"<<"\t\t"<<"totmult from histo"<<"\t"<<"totmult corrected"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		for(int iptbin=0;iptbin<nptv; iptbin++){
                fstrv<<ptbinv[iptbin]<<"-"<<ptbinv[iptbin+1]<<"\t\t"<<totmult[ibin][iptbin]<<"\t"<<totmulthisto[ibin][iptbin]<<"\t"<<totmulthistocorr[ibin][iptbin]<<endl;
		}
		fstrv<<"Integral\t\t"<<totmultall_[ibin]<<"\t"<<hpt[ibin]->Integral(hpt[ibin]->GetXaxis()->FindBin(ptbinv[0]),hpt[ibin]->GetXaxis()->FindBin(3.0)-1)<<"\t"<<hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1)<<endl;
		fstrv<<"V ref="<<(*V_mean)[ibin]<<"\t"<<"V int="<<V_int[ibin]<<"\t"<<"V int corr="<<V_intcorr[ibin]<<endl;
		fstrv<<"V ref err="<<(*deltaV_mean)[ibin]<<"\t"<<"V int err="<<V_interr[ibin]<<"\t"<<"V int corr err="<<V_intcorrerr[ibin]<<endl;
	}

	if(SumorProd=="Sum")TFile *fout = new TFile("mergedv_Prod.root","Recreate");
	else TFile *fout = new TFile("mergedv_Prod2.root","Recreate");
	for(ibin=0;ibin<nbin;ibin++){
	Nevent.Write("Nevent");
	V_int.Write("V_int");
	V_interr.Write("V_interr");
	V_intcorr.Write("V_intcorr");
	V_intcorrerr.Write("V_intcorrerr");
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

