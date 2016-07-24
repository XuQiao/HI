#include <iomanip>
#include <fstream>
#include "par.h"

void getResv(int ispt){

	double theta[ntheta];
	TFile *f[nFileAll];
		
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
        
        int nvv;
        const double *binv;

        if(ispt){             nvv = nptv;       binv = ptbinv;}
            else{             nvv = netav;          binv = etabinv;}

	string SumorProd = getenv("SUMORPROD");
	string mdir = getenv("DIR");
	ofstream fstrv;
        if(ispt){
	if(SumorProd=="Sum")fstrv.open("v.txt");
	else fstrv.open("v_2.txt");
        }
        else{
	if(SumorProd=="Sum")fstrv.open("v_eta.txt");
	else fstrv.open("v_2_eta.txt");
        }
	TVectorD totmult[nbin], totpt[nbin], toteta[nbin];	TVectorD Nevent, totmultall, totmultall_;
	TVectorD V_int, V_interr;
	TVectorD* V_mean;
	TVectorD* deltaV_mean;
	TVectorD* V[nbin]; TVectorD* chi;
	TVectorD V_intcorr, V_intcorrerr;
	TVectorD totmulthisto[nbin], totmulthistocorr[nbin];
	TVectorD* dDRe[nbin]; TVectorD* dDIm[nbin];
	TVectorD* dNRe[nbin][ntheta]; TVectorD* dNIm[nbin][ntheta];
	TComplex dD[nbin][ntheta], dN[nbin][ntheta][nvv];
	TVectorD avgmult[nbin];	TVectorD avgmultall;
	TVectorD avgpt[nbin], avgeta[nbin];
	TVectorD v[nbin][ntheta],vmean[nbin],deltav[nbin][ntheta],deltavmean[nbin];
	if(SumorProd=="Sum")	TFile *infile = TFile::Open("mergedV_Sum.root");
	else	TFile *infile = TFile::Open("mergedV_Prod.root");
	Nevent.ResizeTo(nbin);	Nevent.Zero();
	totmultall.ResizeTo(nbin);	totmultall.Zero();
	totmultall_.ResizeTo(nbin);	totmultall_.Zero();
	avgmultall.ResizeTo(nbin);
	V_int.ResizeTo(nbin);	V_int.Zero();
	V_interr.ResizeTo(nbin);	V_interr.Zero();
        if(ispt){
	V_intcorr.ResizeTo(nbin);	V_intcorr.Zero();
	V_intcorrerr.ResizeTo(nbin);	V_intcorrerr.Zero();
        }
	for(int ibin=0;ibin<nbin;ibin++){
		totmult[ibin].ResizeTo(nvv);	totmult[ibin].Zero();
                if(ispt){
		totmulthisto[ibin].ResizeTo(nvv);	
		totmulthistocorr[ibin].ResizeTo(nvv);
                }
		avgmult[ibin].ResizeTo(nvv);	avgpt[ibin].ResizeTo(nvv);  avgeta[ibin].ResizeTo(nvv);
		totpt[ibin].ResizeTo(nvv);	totpt[ibin].Zero(); 
                toteta[ibin].ResizeTo(nvv);     toteta[ibin].Zero();
		vmean[ibin].ResizeTo(nvv);	deltavmean[ibin].ResizeTo(nvv);
		V[ibin] = (TVectorD*) infile->Get(Form("D_%d/V",ibin));
		chi = (TVectorD*) infile->Get(Form("chi"));
		V_mean = (TVectorD*) infile->Get(Form("Vmean"));
		deltaV_mean = (TVectorD*) infile->Get(Form("deltaVmean"));
		for(int itheta=0;itheta<ntheta;itheta++){
			v[ibin][itheta].ResizeTo(nvv);	deltav[ibin][itheta].ResizeTo(nvv);
			dD[ibin][itheta]=0;
			for(int ivbin=0;ivbin<nvv;ivbin++){
				dN[ibin][itheta][ivbin]=0;
			}
		}
	}
        
        for(int ifile=0; ifile<nFileAll; ifile++){
            if(ispt){
	        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anav_Prod_%d.root",mdir.c_str(),ifile));
	        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anav_Prod2_%d.root",mdir.c_str(),ifile));
            }
            else{
	        if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anaveta_Prod_%d.root",mdir.c_str(),ifile));
	        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to6/%s/Anaveta_Prod2_%d.root",mdir.c_str(),ifile));
            }
		TVectorD* Nevent_t = (TVectorD*)f[ifile]->Get("Nevent");	
		TVectorD* totmultall_t = (TVectorD*)f[ifile]->Get("totmultall");
		for(int ibin=0;ibin<nbin;ibin++){
			TVectorD* totmult_t = (TVectorD*)f[ifile]->Get(Form("totmult_%d",ibin));
			TVectorD* totpt_t = (TVectorD*)f[ifile]->Get(Form("totpt_%d",ibin));
			TVectorD* toteta_t = (TVectorD*)f[ifile]->Get(Form("toteta_%d",ibin));
			dDRe[ibin] = (TVectorD*)f[ifile]->Get(Form("dDRe_%d",ibin));
			dDIm[ibin] = (TVectorD*)f[ifile]->Get(Form("dDIm_%d",ibin));
			Nevent[ibin]+=(*Nevent_t)[ibin];
			totmultall[ibin]+=(*totmultall_t)[ibin];
			for(int ivbin=0;ivbin<nvv;ivbin++){
				totmult[ibin][ivbin]+=(*totmult_t)[ivbin];
				totpt[ibin][ivbin]+=(*totpt_t)[ivbin];
				toteta[ibin][ivbin]+=(*toteta_t)[ivbin];
			}
			for(int itheta=0;itheta<ntheta;itheta++){
				dNRe[ibin][itheta] = (TVectorD*)f[ifile]->Get(Form("dNRe_%d_%d",ibin,itheta));
				dNIm[ibin][itheta] = (TVectorD*)f[ifile]->Get(Form("dNIm_%d_%d",ibin,itheta));
				dD[ibin][itheta]+=TComplex((*dDRe[ibin])[itheta],(*dDIm[ibin])[itheta]);
				for(int ivbin=0;ivbin<nvv;ivbin++)
					dN[ibin][itheta][ivbin]+=TComplex((*dNRe[ibin][itheta])[ivbin],(*dNIm[ibin][itheta])[ivbin]);
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
	        for(int ivbin=0;ivbin<nvv; ivbin++){
			vmean[ibin][ivbin]=0;
			deltavmean[ibin][ivbin]=0;
			avgmult[ibin][ivbin]=1.0*totmult[ibin][ivbin]/Nevent[ibin];
			avgpt[ibin][ivbin]=1.0*totpt[ibin][ivbin]/totmult[ibin][ivbin];
			avgeta[ibin][ivbin]=1.0*toteta[ibin][ivbin]/totmult[ibin][ivbin];
			for(int itheta=0;itheta<ntheta;itheta++){
				dN[ibin][itheta][ivbin]/=totmult[ibin][ivbin];
				TComplex Res=dN[ibin][itheta][ivbin]/dD[ibin][itheta];
				v[ibin][itheta][ivbin]=(*V[ibin])[itheta]*avgmultall[ibin]*TMath::BesselJ1(j01)/Besselj01(mm)*Res.Re();
				vmean[ibin][ivbin]+=v[ibin][itheta][ivbin];
				deltav[ibin][itheta][ivbin]=TMath::Cos(mm*nn*theta[itheta])/totmult[ibin][ivbin]*(TMath::Exp(j01*j01/2./(*chi)[ibin]/(*chi)[ibin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Power(-1,mm)*TMath::Exp(-j01*j01/2./(*chi)[ibin]/(*chi)[ibin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.)));
				deltavmean[ibin][ivbin]+=deltav[ibin][itheta][ivbin];
			//	fstrv<<itheta<<"\t"<<v[ibin][itheta][ivbin]<<"\t"<<deltav[ibin][itheta][ivbin]<<endl;
			}
		deltavmean[ibin][ivbin]=TMath::Sqrt(deltavmean[ibin][ivbin])/2./Besselj01(mm);
		//fstrv<<endl;
		vmean[ibin][ivbin]/=ntheta;
		deltavmean[ibin][ivbin]/=TMath::Sqrt(ntheta);
		fstrv<<binv[ivbin]<<"-"<<binv[ivbin+1]<<"\t"<<vmean[ibin][ivbin]<<"\t"<<deltavmean[ibin][ivbin]<<endl;
                if(ispt)
                totmulthisto[ibin][ivbin]=hpt[ibin]->Integral(hpt[ibin]->GetXaxis()->FindBin(binv[ivbin]),hpt[ibin]->GetXaxis()->FindBin(binv[ivbin+1])-1);
                
		V_int[ibin]+=vmean[ibin][ivbin]*totmult[ibin][ivbin];
		V_interr[ibin]+=deltavmean[ibin][ivbin]*totmult[ibin][ivbin];
		totmultall_[ibin]+=totmult[ibin][ivbin];
                if(ispt){
		if(binv[ivbin+1]>3.0) continue;
                totmulthistocorr[ibin][ivbin]=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(binv[ivbin]),hpteffcorr[ibin]->GetXaxis()->FindBin(binv[ivbin+1])-1);
		V_intcorr[ibin]+=vmean[ibin][ivbin]*totmulthistocorr[ibin][ivbin];
		V_intcorrerr[ibin]+=deltavmean[ibin][ivbin]*totmulthistocorr[ibin][ivbin];
		}
                }
		V_int[ibin]/=totmultall_[ibin];
		V_interr[ibin]/=totmultall_[ibin];
                if(ispt){
		V_intcorr[ibin]/=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(binv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1);
		V_intcorrerr[ibin]/=hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(binv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1);
                }
	}
                if(ispt){
                fstrv<<endl<<"pt range\t\t"<<"totmult"<<"\t\t"<<"totmult from histo"<<"\t"<<"totmult corrected"<<endl;
	        for(int ibin=0;ibin<nbin;ibin++){
		for(int ivbin=0;ivbin<nvv; ivbin++){
                fstrv<<binv[ivbin]<<"-"<<binv[ivbin+1]<<"\t\t"<<totmult[ibin][ivbin]<<"\t"<<totmulthisto[ibin][ivbin]<<"\t"<<totmulthistocorr[ibin][ivbin]<<endl;
		}
		fstrv<<"Integral\t\t"<<totmultall_[ibin]<<"\t"<<hpt[ibin]->Integral(hpt[ibin]->GetXaxis()->FindBin(binv[0]),hpt[ibin]->GetXaxis()->FindBin(3.0)-1)<<"\t"<<hpteffcorr[ibin]->Integral(hpteffcorr[ibin]->GetXaxis()->FindBin(binv[0]),hpteffcorr[ibin]->GetXaxis()->FindBin(3.0)-1)<<endl;
		fstrv<<"V ref="<<(*V_mean)[ibin]<<"\t"<<"V int="<<V_int[ibin]<<"\t"<<"V int corr="<<V_intcorr[ibin]<<endl;
	        fstrv<<"V ref err="<<(*deltaV_mean)[ibin]<<"\t"<<"V int err="<<V_interr[ibin]<<"\t"<<"V int corr err="<<V_intcorrerr[ibin]<<endl;
	    }
                }
            else 
            for(int ibin=0;ibin<nbin;ibin++){
                for(int ivbin=0;ivbin<nvv; ivbin++){
                    fstrv<<binv[ivbin]<<"-"<<binv[ivbin+1]<<"\t\t"<<totmult[ibin][ivbin]<<endl;
                }
		fstrv<<"V ref="<<(*V_mean)[ibin]<<"\t"<<"V int="<<V_int[ibin]<<endl;
		fstrv<<"V ref err="<<(*deltaV_mean)[ibin]<<"\t"<<"V int err="<<V_interr[ibin]<<endl;
            }
        if(ispt){
	if(SumorProd=="Sum")TFile *fout = new TFile("mergedv_Prod.root","Recreate");
	else TFile *fout = new TFile("mergedv_Prod2.root","Recreate");
        }
        else{
	if(SumorProd=="Sum")TFile *fout = new TFile("mergedv_Prod_eta.root","Recreate");
	else TFile *fout = new TFile("mergedv_Prod2_eta.root","Recreate");
        }
	for(ibin=0;ibin<nbin;ibin++){
	Nevent.Write("Nevent");
	V_int.Write("V_int");
	V_interr.Write("V_interr");
        if(ispt){
	V_intcorr.Write("V_intcorr");
	V_intcorrerr.Write("V_intcorrerr");
        }
	TDirectory *dir0 = fout->mkdir(Form("D_%d",ibin));
	dir0->cd();
	avgpt[ibin].Write("avgpt");
	avgeta[ibin].Write("avgeta");
	totmult[ibin].Write("totmult");
	vmean[ibin].Write("vmean");
	deltavmean[ibin].Write("deltavmean");
         for(int itheta=0;itheta<ntheta;itheta++){
            TDirectory *dir1 = dir0->mkdir(Form("D_%d",itheta));dir1->cd();
            v[ibin][itheta].Write("v");
            deltav[ibin][itheta].Write("deltav");
        }
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

