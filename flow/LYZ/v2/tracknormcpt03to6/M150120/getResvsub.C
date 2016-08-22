#include <iomanip>
#include <fstream>
#include "par.h"
#include "../par.h"

void getResvsub(int ispt=1){

	double theta[ntheta];
	TFile *f[nFileAll];
		
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;

        int nvv;
        const double *binv;

        if(ispt==1){             nvv = nptv;       binv = ptbinv;}
            else{             nvv = netav;          binv = etabinv;}

	string SumorProd = getenv("SUMORPROD");
	string mdir = getenv("DIR");
	ofstream fstrv;
        if(ispt==1){
	if(SumorProd=="Sum")fstrv.open("v_sub.txt");
	else fstrv.open("v_2_sub.txt");
        }
        else{
	if(SumorProd=="Sum")fstrv.open("v_eta_sub.txt");
	else fstrv.open("v_2_coarse_eta_sub.txt");
        }
	
	TVectorD Nevent[nbin], totmultall[nbin], totmultall_[nbin];
	TVectorD totmult[nbin], totpt[nbin],toteta[nbin];
	TVectorD totmulthisto[nbin], totmulthistocorr[nbin];
	TVectorD V_int[nbin], V_interr[nbin];
	TVectorD V_intcorr[nbin], V_intcorrerr[nbin];
	TVectorD* V_mean;
	TVectorD* deltaV_mean;
	TVectorD vmeanmean[nbin], deltavmeanmean[nbin], sigmavmeanmean[nbin], avgavgpt[nbin], avgavgeta[nbin];
	TVectorD* V[nbin]; TVectorD* chi[nbin];
	TVectorD* dDRe[nbin]; TVectorD* dDIm[nbin];
	TVectorD* dNRe[nbin][ntheta]; TVectorD* dNIm[nbin][ntheta];
	TComplex dD[nbin][ntheta], dN[nbin][ntheta][nvv];
	TVectorD avgmult[nbin];	TVectorD avgmultall;
	TVectorD avgpt[nbin],avgeta[nbin];
	TVectorD v[nbin][ntheta],vmean[nbin],deltav[nbin][ntheta],deltavmean[nbin];
	TVectorD IFILE[nbin];
	if(SumorProd=="Sum")	TFile *infile = TFile::Open("mergedV_Sum.root");
	else	TFile *infile = TFile::Open("mergedV_Prod.root");
        if(ispt==1){
	if(SumorProd=="Sum")TFile *fout = new TFile("mergedv_Prod_sub.root","Recreate");
	else TFile *fout = new TFile("mergedv_Prod2_sub.root","Recreate");
        }
        else if(ispt==-1){
	if(SumorProd=="Sum")TFile *fout = new TFile("mergedv_Prod_coarse_eta_sub.root","Recreate");
	else TFile *fout = new TFile("mergedv_Prod2_coarse_eta_sub.root","Recreate");
        }
        else{
	if(SumorProd=="Sum")TFile *fout = new TFile("mergedv_Prod_eta_sub.root","Recreate");
	else TFile *fout = new TFile("mergedv_Prod2_eta_sub.root","Recreate");
        }

	TVectorD Nevent_;       Nevent_.ResizeTo(nbin);  Nevent_.Zero();
	TH1D* hpt[nbin];
	TH1D* hpteffcorr[nbin];
	for(int ibin=0;ibin<nbin;ibin++){
        for(int ifile=0; ifile<nFileAll; ifile++){
            if(ispt==1){
        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anav_Prod_%d.root",mdir.c_str(),ifile));
        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anav_Prod2_%d.root",mdir.c_str(),ifile));
            }
            else if(ispt==-1){
        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/coarsebins/%s/Anaveta_Prod_%d.root",mdir.c_str(),ifile));
        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/coarsebins/%s/Anaveta_Prod2_%d.root",mdir.c_str(),ifile));
            }
            else{
        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anaveta_Prod_%d.root",mdir.c_str(),ifile));
        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anaveta_Prod2_%d.root",mdir.c_str(),ifile));
            }
        TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
        Nevent_[ibin]+=(*Nevent_t)[ibin];
        f[ifile]->Close();
        }
	IFILE[ibin].ResizeTo(nsamples);
        Nevent[ibin].ResizeTo(nsamples);  Nevent[ibin].Zero();
        totmultall[ibin].ResizeTo(nsamples);      totmultall[ibin].Zero();
        totmultall_[ibin].ResizeTo(nsamples);      totmultall_[ibin].Zero();
        TDirectory *dir0 = fout->mkdir(Form("D_%d",ibin));


      	vmeanmean[ibin].ResizeTo(nvv);  vmeanmean[ibin].Zero();
        deltavmeanmean[ibin].ResizeTo(nvv);  deltavmeanmean[ibin].Zero();
       	sigmavmeanmean[ibin].ResizeTo(nvv);  sigmavmeanmean[ibin].Zero();
	avgavgpt[ibin].ResizeTo(nvv);	avgavgpt[ibin].Zero();
	avgavgeta[ibin].ResizeTo(nvv);	avgavgeta[ibin].Zero();
	V_interr[ibin].ResizeTo(nsamples);	V_interr[ibin].Zero();
	V_int[ibin].ResizeTo(nsamples);	V_int[ibin].Zero();
        if(ispt==1){
	V_intcorr[ibin].ResizeTo(nsamples);	V_intcorr[ibin].Zero();
	V_intcorrerr[ibin].ResizeTo(nsamples);	V_intcorrerr[ibin].Zero();
        }
        for(int isample=0;isample<nsamples;isample++){

	TVectorD Nevent0; Nevent0.ResizeTo(nbin);	Nevent0.Zero();
	avgmultall.ResizeTo(nbin);
		totmult[ibin].ResizeTo(nvv);	totmult[ibin].Zero();
                if(ispt==1){
		totmulthisto[ibin].ResizeTo(nvv); totmulthisto[ibin].Zero();
		totmulthistocorr[ibin].ResizeTo(nvv);totmulthistocorr[ibin].Zero();
                }
		avgmult[ibin].ResizeTo(nvv);	avgpt[ibin].ResizeTo(nvv);  avgeta[ibin].ResizeTo(nvv);
		totpt[ibin].ResizeTo(nvv);	totpt[ibin].Zero();
		toteta[ibin].ResizeTo(nvv);	toteta[ibin].Zero();
		vmean[ibin].ResizeTo(nvv);	deltavmean[ibin].ResizeTo(nvv);
		V[ibin] = (TVectorD*) infile->Get(Form("D_%d/V",ibin,isample));
		chi[ibin] = (TVectorD*) infile->Get(Form("chi",ibin,isample));
		V_mean = (TVectorD*) infile->Get(Form("Vmean",ibin,isample));
		deltaV_mean = (TVectorD*) infile->Get(Form("deltaVmean",ibin,isample));
		for(int itheta=0;itheta<ntheta;itheta++){
			v[ibin][itheta].ResizeTo(nvv);	deltav[ibin][itheta].ResizeTo(nvv);
			dD[ibin][itheta]=0;
			for(int ivbin=0;ivbin<nvv;ivbin++){
				dN[ibin][itheta][ivbin]=0;
			}
		}

        for(int ifile=0; ifile<nFileAll; ifile++){
            if(ispt==1){
	        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anav_Prod_%d.root",mdir.c_str(),ifile));
	        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anav_Prod2_%d.root",mdir.c_str(),ifile));
            }
            else if(ispt==-1){
	        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/coarsebins/%s/Anaveta_Prod_%d.root",mdir.c_str(),ifile));
	        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/coarsebins/%s/Anaveta_Prod2_%d.root",mdir.c_str(),ifile));
            }
            else{
	        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anaveta_Prod_%d.root",mdir.c_str(),ifile));
	        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anaveta_Prod2_%d.root",mdir.c_str(),ifile));
            }
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
			TVectorD* toteta_t = (TVectorD*)f[ifile]->Get(Form("toteta_%d",ibin));
			dDRe[ibin] = (TVectorD*)f[ifile]->Get(Form("dDRe_%d",ibin));
			dDIm[ibin] = (TVectorD*)f[ifile]->Get(Form("dDIm_%d",ibin));
			totmultall[ibin][isample]+=(*totmultall_t)[ibin];
			for(int ivbin=0;ivbin<nvv;ivbin++){
				totmult[ibin][ivbin]+=(*totmult_t)[ivbin];
				totpt[ibin][ivbin]+=(*totpt_t)[ivbin];
				toteta[ibin][ivbin]+=(*toteta_t)[ivbin];
			}
			for(int itheta=0;itheta<ntheta;itheta++){
				dNRe[ibin][itheta] = (TVectorD*)f[ifile]->Get(Form("dNRe_%d_%d",ibin,itheta));
				dNIm[ibin][itheta] = (TVectorD*)f[ifile]->Get(Form("dNIm_%d_%d",ibin,itheta));
				dD[ibin][itheta]+=TComplex((*dDRe[ibin])[itheta],(*dDIm[ibin])[itheta]);
				for(int ivbin=0;ivbin<nvv;ivbin++){
                                    if(ispt==1)
					dN[ibin][itheta][ivbin]+=TComplex((*dNRe[ibin][itheta])[ivbin],(*dNIm[ibin][itheta])[ivbin]);
                                    else{
					if(ifile<188) dN[ibin][itheta][ivbin]+=TComplex((*dNRe[ibin][itheta])[ivbin],(*dNIm[ibin][itheta])[ivbin]);
                                        else dN[ibin][itheta][nvv-ivbin-1]+=TComplex((*dNRe[ibin][itheta])[ivbin],(*dNIm[ibin][itheta])[ivbin]);
                                }
                                }
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
                //TH1D* hptre = (TH1D*)hpt[ibin]->Rebin(NbinY,"hptre",ptbinhisto);
                //hpteffcorr[ibin] = (TH1D*)hptre->Clone(Form("hpteffcorr_%d",ibin));
                //hpteffcorr[ibin]->Divide(hpteff);
		for(int ivbin=0;ivbin<nvv; ivbin++){
			vmean[ibin][ivbin]=0;
			deltavmean[ibin][ivbin]=0;
			avgmult[ibin][ivbin]=1.0*totmult[ibin][ivbin]/Nevent[ibin][isample];
			avgpt[ibin][ivbin]=1.0*totpt[ibin][ivbin]/totmult[ibin][ivbin];
			avgeta[ibin][ivbin]=1.0*toteta[ibin][ivbin]/totmult[ibin][ivbin];
			for(int itheta=0;itheta<ntheta;itheta++){
				dN[ibin][itheta][ivbin]/=totmult[ibin][ivbin];
				TComplex Res=dN[ibin][itheta][ivbin]/dD[ibin][itheta];
				v[ibin][itheta][ivbin]=(*V[ibin])[itheta]*avgmultall[ibin]*TMath::BesselJ1(j01)/Besselj01(mm)*Res.Re();
				vmean[ibin][ivbin]+=v[ibin][itheta][ivbin];
				deltav[ibin][itheta][ivbin]=TMath::Cos(mm*nn*theta[itheta])/totmult[ibin][ivbin]*(TMath::Exp(j01*j01/2./(*chi[ibin])[0]/(*chi[ibin])[0]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Power(-1,mm)*TMath::Exp(-j01*j01/2./(*chi[ibin])[0]/(*chi[ibin])[0]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.)));
				deltavmean[ibin][ivbin]+=deltav[ibin][itheta][ivbin];
			//	fstrv<<itheta<<"\t"<<v[ibin][itheta][ivbin]<<"\t"<<deltav[ibin][itheta][ivbin]<<endl;
			}
		deltavmean[ibin][ivbin]=TMath::Sqrt(deltavmean[ibin][ivbin])/2./Besselj01(mm);
		//fstrv<<endl;
		vmean[ibin][ivbin]/=ntheta;
		deltavmean[ibin][ivbin]/=TMath::Sqrt(ntheta);
		fstrv<<binv[ivbin]<<"-"<<binv[ivbin+1]<<"\t"<<vmean[ibin][ivbin]<<"\t"<<deltavmean[ibin][ivbin]<<endl;
                if(ispt==1){
		if(binv[ivbin+1]>3.0) continue;
                totmulthisto[ibin][ivbin]=hpt[ibin]->Integral(hpt[ibin]->GetXaxis()->FindBin(binv[ivbin]),hpt[ibin]->GetXaxis()->FindBin(binv[ivbin+1])-1);
//                totmulthistocorr[ibin][ivbin]=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(binv[ivbin]),hpteffcorr[ibin]->GetXaxis()->FindBin(binv[ivbin+1])-1);
                totmulthistocorr[ibin][ivbin]=totmulthisto[ibin][ivbin]/hpteff->GetBinContent(hpteff->FindBin((binv[ivbin]+binv[ivbin+1])/2));
		V_int[ibin][isample]+=vmean[ibin][ivbin]*totmult[ibin][ivbin];
		V_interr[ibin][isample]+=deltavmean[ibin][ivbin]*totmult[ibin][ivbin];
		V_intcorr[ibin][isample]+=vmean[ibin][ivbin]*totmulthistocorr[ibin][ivbin];
		V_intcorrerr[ibin][isample]+=deltavmean[ibin][ivbin]*totmulthistocorr[ibin][ivbin];
                }
		totmultall_[ibin][isample]+=totmult[ibin][ivbin];
		}
		V_int[ibin][isample]/=totmultall_[ibin][isample];
		V_interr[ibin][isample]/=totmultall_[ibin][isample];
		//V_intcorr[ibin][isample]/=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(binv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1);
		//V_intcorr[ibin][isample]/=totmulthistocorr[ibin];
		//V_intcorrerr[ibin][isample]/=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(binv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1);

                if(ispt==1){
        fstrv<<endl<<"pt range\t\t"<<"totmult"<<"\t\t"<<"totmult from histo"<<"\t"<<"totmult corrected"<<endl;
		for(int ivbin=0;ivbin<nvv; ivbin++){
                fstrv<<binv[ivbin]<<"-"<<binv[ivbin+1]<<"\t\t"<<totmult[ibin][ivbin]<<"\t"<<totmulthisto[ibin][ivbin]<<"\t"<<totmulthistocorr[ibin][ivbin]<<endl;
		}
		//fstrv<<"Integral\t\t"<<totmultall[ibin][isample]<<"\t"<<hpt[ibin]->Integral(hpt[ibin]->GetXaxis()->FindBin(binv[0]),hpt[ibin]->GetXaxis()->FindBin(3.0)-1)<<"\t"<<hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(binv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1)<<endl;
		fstrv<<"V ref="<<(*V_mean)[ibin]<<"\t"<<"V int="<<V_int[ibin][isample]<<"\t"<<"V int corr="<<V_intcorr[ibin][isample]<<endl;
		fstrv<<"V ref err="<<(*deltaV_mean)[ibin]<<"\t"<<"V int err="<<V_interr[ibin][isample]<<"\t"<<"V int corr err="<<V_intcorrerr[ibin][isample]<<endl;
                }
	TDirectory *dirsample = dir0->mkdir(Form("s_%d",isample));
        dirsample->cd();
	avgpt[ibin].Write("avgpt");
	avgeta[ibin].Write("avgeta");
	totmult[ibin].Write("totmult");
	vmean[ibin].Write("vmean");
	deltavmean[ibin].Write("deltavmean");
	for(int ivbin=0;ivbin<nvv; ivbin++){
		avgavgpt[ibin][ivbin]+=avgpt[ibin][ivbin]*Nevent[ibin][isample];
		avgavgeta[ibin][ivbin]+=avgeta[ibin][ivbin]*Nevent[ibin][isample];
		vmeanmean[ibin][ivbin]+=vmean[ibin][ivbin]/nsamples;
        	deltavmeanmean[ibin][ivbin]+=deltavmean[ibin][ivbin]/nsamples;
	        sigmavmeanmean[ibin][ivbin]+=TMath::Power(vmean[ibin][ivbin]/nsamples,2);
	}
	}//subsample loop
	dir0->cd();
        IFILE[ibin].Write("IFILE");
	V_int[ibin].Write("V_int");
	V_interr[ibin].Write("V_interr");
	V_intcorr[ibin].Write("V_intcorr");
	V_intcorrerr[ibin].Write("V_intcorrerr");
	Nevent[ibin].Write("Nevent");
	totmultall[ibin].Write("totmultall");
	for(int ivbin=0;ivbin<nvv; ivbin++){
		avgavgpt[ibin][ivbin]/=Nevent_[ibin];
		avgavgeta[ibin][ivbin]/=Nevent_[ibin];
        	sigmavmeanmean[ibin][ivbin]=TMath::Sqrt(sigmavmeanmean[ibin][ivbin]*nsamples-vmeanmean[ibin][ivbin]*vmeanmean[ibin][ivbin])/TMath::Sqrt(nsamples);
		}
	avgavgpt[ibin].Write("avgavgpt");
	avgavgeta[ibin].Write("avgavgeta");
        vmeanmean[ibin].Write("vmeanmean");
	deltavmeanmean[ibin].Write("deltavmeanmean");
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

