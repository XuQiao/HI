#include <iomanip>
#include <fstream>
#include "par.h"

void getResvsub(){

	double theta[ntheta];
	TFile *f[nFileAll];
		
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;


	string SumorProd = getenv("SUMORPROD");
	string mdir = getenv("DIR");
        int iloop = atoi(getenv("ILOOP"));
	ofstream fstrv;
//	if(SumorProd=="Sum")fstrv.open("v_sub.txt");
//	else fstrv.open("v_2_sub.txt");
	fstrv.open(Form("v_2_sub_loop%d.txt",iloop));
	
	TVectorD Nevent[nbin], totmultall[nbin], totmultall_[nbin], avgavgpt[nbin];
	TVectorD totmult[nbin], totpt[nbin];
	TVectorD totmulthisto[nbin], totmulthistocorr[nbin];
	TVectorD V_int[nbin], V_interr[nbin];
	TVectorD V_intcorr[nbin], V_intcorrerr[nbin];
	TVectorD* V_mean;
	TVectorD* deltaV_mean;
	TVectorD vmeanmean[nbin], deltavmeanmean[nbin], sigmavmeanmean[nbin];
	TVectorD* V[nbin]; TVectorD* chi[nbin];
	TVectorD* dDRe[nbin]; TVectorD* dDIm[nbin];
	TVectorD* dNRe[nbin][ntheta]; TVectorD* dNIm[nbin][ntheta];
	TComplex dD[nbin][ntheta], dN[nbin][ntheta][nptv];
	TVectorD avgmult[nbin];	TVectorD avgmultall;
	TVectorD avgpt[nbin];
	TVectorD v[nbin][ntheta],vmean[nbin],deltav[nbin][ntheta],deltavmean[nbin];
	TVectorD IFILE[nbin];
//	if(SumorProd=="Sum")	TFile *infile = TFile::Open("mSum_sub.root");
//	else	TFile *infile = TFile::Open("mergedV_Prod_sub.root");
	TFile *infile = TFile::Open(Form("parin_%d.root",iloop));
//	if(SumorProd=="Sum")TFile *fout = new TFile("mergedv_Prod_sub.root","Recreate");
//       else TFile *fout = new TFile("mergedv_Prod2_sub.root","Recreate");
        TFile *fout = new TFile(Form("mergedv_Prod2_sub_loop%d.root",iloop),"Recreate");

	 TVectorD Nevent_;       Nevent_.ResizeTo(nbin);  Nevent_.Zero();
	TH1D* hpt[nbin];
	TH1D* hpteffcorr[nbin];
	for(int ibin=0;ibin<nbin;ibin++){
        for(int ifile=0; ifile<nFileAll; ifile++){
        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/v2sensdebug/tracknormcpt03to6/loop%d/%s/Anav_Prod_%d.root",iloop,mdir.c_str(),ifile));
        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/v2sensdebug/tracknormcpt03to6/loop%d/%s/Anav_Prod2_%d.root",iloop,mdir.c_str(),ifile));
        TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
        Nevent_[ibin]+=(*Nevent_t)[ibin];
        f[ifile]->Close();
        }
	IFILE[ibin].ResizeTo(nsamples);
        Nevent[ibin].ResizeTo(nsamples);  Nevent[ibin].Zero();
        totmultall[ibin].ResizeTo(nsamples);      totmultall[ibin].Zero();
        totmultall_[ibin].ResizeTo(nsamples);      totmultall_[ibin].Zero();
        TDirectory *dir0 = fout->mkdir(Form("D_%d",ibin));


      	vmeanmean[ibin].ResizeTo(nptv);  vmeanmean[ibin].Zero();
        deltavmeanmean[ibin].ResizeTo(nptv);  deltavmeanmean[ibin].Zero();
       	sigmavmeanmean[ibin].ResizeTo(nptv);  sigmavmeanmean[ibin].Zero();
	avgavgpt[ibin].ResizeTo(nptv);	avgavgpt[ibin].Zero();
	V_interr[ibin].ResizeTo(nsamples);	V_interr[ibin].Zero();
	V_int[ibin].ResizeTo(nsamples);	V_int[ibin].Zero();
	V_intcorr[ibin].ResizeTo(nsamples);	V_intcorr[ibin].Zero();
	V_intcorrerr[ibin].ResizeTo(nsamples);	V_intcorrerr[ibin].Zero();
        for(int isample=0;isample<nsamples;isample++){

	TVectorD Nevent0; Nevent0.ResizeTo(nbin);	Nevent0.Zero();
	avgmultall.ResizeTo(nbin);
		totmult[ibin].ResizeTo(nptv);	totmult[ibin].Zero();
		totmulthisto[ibin].ResizeTo(nptv);	
		totmulthistocorr[ibin].ResizeTo(nptv);
		avgmult[ibin].ResizeTo(nptv);	avgpt[ibin].ResizeTo(nptv);
		totpt[ibin].ResizeTo(nptv);	totpt[ibin].Zero();
		vmean[ibin].ResizeTo(nptv);	deltavmean[ibin].ResizeTo(nptv);
		V[ibin] = (TVectorD*) infile->Get(Form("D_%d/D_0/Vcalc",ibin));
	//	chi[ibin] = (TVectorD*) infile->Get(Form("D_%d/s_%d/chi",ibin,isample));
		V_mean = (TVectorD*) infile->Get(Form("D_%d/Vcalcmean",ibin));
	//	deltaV_mean = (TVectorD*) infile->Get(Form("D_%d/s_%d/deltaVmean",ibin,isample));
		for(int itheta=0;itheta<ntheta;itheta++){
			v[ibin][itheta].ResizeTo(nptv);	deltav[ibin][itheta].ResizeTo(nptv);
			dD[ibin][itheta]=0;
			for(int iptbin=0;iptbin<nptv;iptbin++){
				dN[ibin][itheta][iptbin]=0;
			}
		}

        for(int ifile=0; ifile<nFileAll; ifile++){
                if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/v2sensdebug/tracknormcpt03to6/loop%d/%s/Anav_Prod_%d.root",iloop,mdir.c_str(),ifile));
                else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/v2sensdebug/tracknormcpt03to6/loop%d/%s/Anav_Prod2_%d.root",iloop,mdir.c_str(),ifile));
		TVectorD* Nevent_t = (TVectorD*)f[ifile]->Get("Nevent");	
		TVectorD* totmultall_t = (TVectorD*)f[ifile]->Get("totmultall");
		Nevent0[ibin] += (*Nevent_t)[ibin];
                double start=isample*Nevent_[ibin]/nsamples;
                double end=(isample+1)*Nevent_[ibin]/nsamples;
                if(Nevent0[ibin]>start && Nevent0[ibin]<=end){
                        IFILE[ibin][isample]=(double)ifile;
                        Nevent[ibin][isample] += (*Nevent_t)[ibin];
			TVectorD* totmult_t = (TVectorD*)f[ifile]->Get(Form("totmult_%d",ibin));
			TVectorD* totpt_t = (TVectorD*)f[ifile]->Get(Form("totpt_%d",ibin));
			dDRe[ibin] = (TVectorD*)f[ifile]->Get(Form("dDRe_%d",ibin));
			dDIm[ibin] = (TVectorD*)f[ifile]->Get(Form("dDIm_%d",ibin));
			totmultall[ibin][isample]+=(*totmultall_t)[ibin];
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
	
	TFile *fhisto = TFile::Open("histomerged.root");
	TFile *feff = TFile::Open("/home/xuq7/HI/flow/LYZ/v2/TrackCorrections_HIJING_538_OFFICIAL_Mar24.root");
	TH2F* heff = (TH2F*)feff->Get("rTotalEff3D");
      	TH1D* hpteff = (TH1D*)heff->ProjectionY("hpteff",heff->GetXaxis()->FindBin(-2.4),heff->GetXaxis()->FindBin(2.4)-1,"o");
	TArrayD *ptarr = (TArrayD*)heff->GetYaxis()->GetXbins();
	double *ptbinhisto = ptarr->GetArray();
	int NbinX = heff->GetXaxis()->GetNbins();
	int NbinY = heff->GetYaxis()->GetNbins();
	hpteff->Scale(1.0/NbinX);
        
		for(int itheta=0;itheta<ntheta;itheta++)
			dD[ibin][itheta]*=TComplex::Power(TComplex::I(),mm-1)/(Double_t)Nevent[ibin][isample];
	
		avgmultall[ibin]=1.0*totmultall[ibin][isample]/Nevent[ibin][isample];
                hpt[ibin] = (TH1D*)fhisto->Get(Form("D_%d/hpt",ibin));
                TH1D* hptre = (TH1D*)hpt[ibin]->Rebin(NbinY,"hptre",ptbinhisto);
                hpteffcorr[ibin] = (TH1D*)hptre->Clone(Form("hpteffcorr_%d",ibin));
                hpteffcorr[ibin]->Divide(hpteff);
		for(int iptbin=0;iptbin<nptv; iptbin++){
			vmean[ibin][iptbin]=0;
	//		deltavmean[ibin][iptbin]=0;
			avgmult[ibin][iptbin]=1.0*totmult[ibin][iptbin]/Nevent[ibin][isample];
			avgpt[ibin][iptbin]=1.0*totpt[ibin][iptbin]/totmult[ibin][iptbin];
			for(int itheta=0;itheta<ntheta;itheta++){
				dN[ibin][itheta][iptbin]/=totmult[ibin][iptbin];
				TComplex Res=dN[ibin][itheta][iptbin]/dD[ibin][itheta];
				v[ibin][itheta][iptbin]=(*V[ibin])[itheta]*avgmultall[ibin]*TMath::BesselJ1(j01)/Besselj01(mm)*Res.Re();
				vmean[ibin][iptbin]+=v[ibin][itheta][iptbin];
	//			deltav[ibin][itheta][iptbin]=TMath::Cos(mm*nn*theta[itheta])/totmult[ibin][iptbin]*(TMath::Exp(j01*j01/2./(*chi[ibin])[0]/(*chi[ibin])[0]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Power(-1,mm)*TMath::Exp(-j01*j01/2./(*chi[ibin])[0]/(*chi[ibin])[0]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.)));
	//			deltavmean[ibin][iptbin]+=deltav[ibin][itheta][iptbin];
			//	fstrv<<itheta<<"\t"<<v[ibin][itheta][iptbin]<<"\t"<<deltav[ibin][itheta][iptbin]<<endl;
			}
	//	deltavmean[ibin][iptbin]=TMath::Sqrt(deltavmean[ibin][iptbin])/2./Besselj01(mm);
		//fstrv<<endl;
		vmean[ibin][iptbin]/=ntheta;
	//	deltavmean[ibin][iptbin]/=TMath::Sqrt(ntheta);
		fstrv<<ptbinv[iptbin]<<"-"<<ptbinv[iptbin+1]<<"\t"<<vmean[ibin][iptbin]<<"\t"<<deltavmean[ibin][iptbin]<<endl;
		if(ptbinv[iptbin+1]>3.0) continue;
                totmulthisto[ibin][iptbin]=hpt[ibin]->Integral(hpt[ibin]->GetXaxis()->FindBin(ptbinv[iptbin]),hpt[ibin]->GetXaxis()->FindBin(ptbinv[iptbin+1])-1);
                totmulthistocorr[ibin][iptbin]=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[iptbin]),hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[iptbin+1])-1);
		V_int[ibin][isample]+=vmean[ibin][iptbin]*totmult[ibin][iptbin];
	//	V_interr[ibin][isample]+=deltavmean[ibin][iptbin]*totmult[ibin][iptbin];
		V_intcorr[ibin][isample]+=vmean[ibin][iptbin]*totmulthistocorr[ibin][iptbin];
	//	V_intcorrerr[ibin][isample]+=deltavmean[ibin][iptbin]*totmulthistocorr[ibin][iptbin];
		totmultall_[ibin][isample]+=totmult[ibin][iptbin];
		}
		V_int[ibin][isample]/=totmultall_[ibin][isample];
	//	V_interr[ibin][isample]/=totmultall_[ibin][isample];
		V_intcorr[ibin][isample]/=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1);
	//	V_intcorrerr[ibin][isample]/=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1);

        fstrv<<endl<<"pt range\t\t"<<"totmult"<<"\t\t"<<"totmult from histo"<<"\t"<<"totmult corrected"<<endl;
		for(int iptbin=0;iptbin<nptv; iptbin++){
                fstrv<<ptbinv[iptbin]<<"-"<<ptbinv[iptbin+1]<<"\t\t"<<totmult[ibin][iptbin]<<"\t"<<totmulthisto[ibin][iptbin]<<"\t"<<totmulthistocorr[ibin][iptbin]<<endl;
		}
		fstrv<<"Integral\t\t"<<totmultall[ibin][isample]<<"\t"<<hpt[ibin]->Integral(hpt[ibin]->GetXaxis()->FindBin(ptbinv[0]),hpt[ibin]->GetXaxis()->FindBin(3.0)-1)<<"\t"<<hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(ptbinv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1)<<endl;
		fstrv<<"V ref="<<(*V_mean)[0]<<"\t"<<"V int="<<V_int[ibin][isample]<<"\t"<<"V int corr="<<V_intcorr[ibin][isample]<<endl;
//		fstrv<<"V ref err="<<(*deltaV_mean)[0]<<"\t"<<"V int err="<<V_interr[ibin][isample]<<"\t"<<"V int corr err="<<V_intcorrerr[ibin][isample]<<endl;
	
	TDirectory *dirsample = dir0->mkdir(Form("s_%d",isample));
        dirsample->cd();
	avgpt[ibin].Write("avgpt");
	totmult[ibin].Write("totmult");
	vmean[ibin].Write("vmean");
//	deltavmean[ibin].Write("deltavmean");
	for(int iptbin=0;iptbin<nptv; iptbin++){
		avgavgpt[ibin][iptbin]+=avgpt[ibin][iptbin]*Nevent[ibin][isample];
		vmeanmean[ibin][iptbin]+=vmean[ibin][iptbin]/nsamples;
 //       	deltavmeanmean[ibin][iptbin]+=deltavmean[ibin][iptbin]/nsamples;
	        sigmavmeanmean[ibin][iptbin]+=TMath::Power(vmean[ibin][iptbin]/nsamples,2);
	}
	}//subsample loop
	dir0->cd();
        IFILE[ibin].Write("IFILE");
	V_int[ibin].Write("V_int");
//	V_interr[ibin].Write("V_interr");
	V_intcorr[ibin].Write("V_intcorr");
//	V_intcorrerr[ibin].Write("V_intcorrerr");
	Nevent[ibin].Write("Nevent");
	totmultall[ibin].Write("totmultall");
	for(int iptbin=0;iptbin<nptv; iptbin++){
		avgavgpt[ibin][iptbin]/=Nevent_[ibin];
        	sigmavmeanmean[ibin][iptbin]=TMath::Sqrt(sigmavmeanmean[ibin][iptbin]*nsamples-vmeanmean[ibin][iptbin]*vmeanmean[ibin][iptbin])/TMath::Sqrt(nsamples);
		}
	avgavgpt[ibin].Write("avgavgpt");
        vmeanmean[ibin].Write("vmeanmean");
//	deltavmeanmean[ibin].Write("deltavmeanmean");
        sigmavmeanmean[ibin].Write("sigmavmeanmean");
	}//ntrk bin loop
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

