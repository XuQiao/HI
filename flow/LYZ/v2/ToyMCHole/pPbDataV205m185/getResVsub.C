#include <iomanip>
#include <fstream>
#include "par.h"

//------------------------Deal with the number and get the output stored in txt files and root files
void getResVsub(){
	double Vmax[nbin], eps[nbin];
	for(int ibin=0; ibin<nbin ;ibin++){
        	Vmax[ibin]=0.065*(trkbin[ibin]+30);
        	eps[ibin]=0.00025*(trkbin[ibin]+30);
	}
	ofstream  fstrV;
	double theta[ntheta];
	TVectorD avgmult, avgpt, totptall;
	TVectorD Qx1, Qy1, Q2;
	TVectorD Gmod2[nbin][ntheta];
	TVectorD sigma2[nbin],deltaV[nbin];
	TVectorD sigma2_,chi_;
	TVectorD deltaVmean, Vmean;
	TVectorD deltaVmeanmean, Vmeanmean, sigmaVmeanmean;
	TVectorD r[nbin];
	TVectorD r0[nbin], r01[nbin], V[nbin], chi[nbin];
	TVectorD GRe[nbin][ntheta]; TVectorD* GRe_t[nbin][ntheta];
	TVectorD GIm[nbin][ntheta]; TVectorD* GIm_t[nbin][ntheta];
	TComplex G[nbin][ntheta][nstepr];
	if(isSum)	fstrV.open("V_Sum_sub.txt");
	else	fstrV.open("V_Prod_sub.txt");
	if(isSum)       TString outname = "mergedV_Sum_sub.root";
        else    TString outname="mergedV_Prod_sub.root";
        TFile *outf = new TFile(outname,"Recreate");
        int subsample = (int)(nFileAll/nsamples);
	TFile *f[nFileAll];

        for(int ibin=0; ibin<nbin; ibin++){
		r[ibin].ResizeTo(nstepr);
                for(int ir=0; ir<nstepr; ir++){
        if(isSimple==0)  r[ibin][ir]=j01/(Vmax[ibin]-eps[ibin]*ir);	
	else		 r[ibin][ir]=0.00025*20*(ir+1);
		}
	}
               
	Vmeanmean.ResizeTo(nbin);  Vmeanmean.Zero();
	deltaVmeanmean.ResizeTo(nbin);  deltaVmeanmean.Zero();
	sigmaVmeanmean.ResizeTo(nbin);  sigmaVmeanmean.Zero();
	for(int isample=0;isample<nsamples;isample++){
	cout<<"start with "<<isample<<"\t";
        int start=isample*subsample;
        int end=(isample+1)*subsample;
	TVectorD Nevent;	Nevent.ResizeTo(nbin);  Nevent.Zero();
        TVectorD totmultall;	totmultall.ResizeTo(nbin);      totmultall.Zero();
	
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;

	for(int ibin=0;ibin<nbin;ibin++){
			r0[ibin].ResizeTo(ntheta);
			r01[ibin].ResizeTo(ntheta);
			sigma2[ibin].ResizeTo(ntheta);
			V[ibin].ResizeTo(ntheta);
			deltaV[ibin].ResizeTo(ntheta);
			chi[ibin].ResizeTo(ntheta);
                        for(int itheta=0;itheta<ntheta;itheta++){
				Gmod2[ibin][itheta].ResizeTo(nstepr);
				GRe[ibin][itheta].ResizeTo(nstepr);
				GRe[ibin][itheta].Zero();
				GIm[ibin][itheta].ResizeTo(nstepr);
				GIm[ibin][itheta].Zero();
			}
		}

		avgmult.ResizeTo(nbin);	deltaVmean.ResizeTo(nbin);	Vmean.ResizeTo(nbin);
		deltaVmean.Zero();	Vmean.Zero();	deltaVmean.Zero();
		totptall.ResizeTo(nbin);	totptall.Zero();
		avgpt.ResizeTo(nbin);
                Qx1.ResizeTo(nbin);     Qx1.Zero();
                Qy1.ResizeTo(nbin);     Qy1.Zero();
                Q2.ResizeTo(nbin);     Q2.Zero();
		sigma2_.ResizeTo(nbin);sigma2_.Zero();
		chi_.ResizeTo(nbin);chi_.Zero();

	for(int ifile=start; ifile<end; ifile++){
		if(isSum) f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/ToyMCHole/pPbDataV205m185/AnaV_Sum_%d.root",ifile));
		else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/ToyMCHole/pPbDataV205m185/AnaV_Prod_%d.root",ifile));
		TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
		TVectorD* totmultall_t =  (TVectorD*)f[ifile]->Get(Form("totmultall"));
			TVectorD* Qx1_t =  (TVectorD*)f[ifile]->Get(Form("Qx1"));
			TVectorD* Qy1_t =  (TVectorD*)f[ifile]->Get(Form("Qy1"));
			TVectorD* Q2_t =  (TVectorD*)f[ifile]->Get(Form("Q2"));
			TVectorD* totptall_t =  (TVectorD*)f[ifile]->Get(Form("totptall"));
		for(int ibin=0;ibin<nbin;ibin++){
					for(int itheta=0;itheta<ntheta;itheta++){
					GRe_t[ibin][itheta] = (TVectorD*)f[ifile]->Get(Form("GRe_%d_%d",ibin,itheta));
					GIm_t[ibin][itheta] = (TVectorD*)f[ifile]->Get(Form("GIm_%d_%d",ibin,itheta));
					for(ir=0; ir<nstepr; ir++){
						GRe[ibin][itheta][ir] += (*GRe_t[ibin][itheta])[ir];
						GIm[ibin][itheta][ir] += (*GIm_t[ibin][itheta])[ir];
					}
				}
				Qx1[ibin] += (*Qx1_t)[ibin];
				Qy1[ibin] += (*Qy1_t)[ibin];
				Q2[ibin] += (*Q2_t)[ibin];
				totptall[ibin] += (*totptall_t)[ibin];
			Nevent[ibin] += (*Nevent_t)[ibin];
			totmultall[ibin] += (*totmultall_t)[ibin];	
		}
		f[ifile]->Close();
	}
		
	fstrV<<setprecision(4)<<fixed;
	fstrV<<"ibin"<<"\t"<<"itheta"<<"\t"<<"r0"<<"\t"<<"V"<<"\t"<<"sigma2"<<"\t"<<"chi"<<"\t"<<"Vn Errors"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
			for(int itheta=0;itheta<ntheta;itheta++){
				for(ir=0; ir<nstepr; ir++){
					G[ibin][itheta][ir]=TComplex(GRe[ibin][itheta][ir],GIm[ibin][itheta][ir]);	
					G[ibin][itheta][ir]/=Nevent[ibin];
					Gmod2[ibin][itheta][ir]=TMath::Power(TComplex::Abs(G[ibin][itheta][ir]),2);
				}
				for(ir=0; ir<nstepr-1; ir++)
					if(ir!=0 && Gmod2[ibin][itheta][ir]<=Gmod2[ibin][itheta][ir-1] && Gmod2[ibin][itheta][ir]<=Gmod2[ibin][itheta][ir+1]) break;
				if(ir!=0 && ir<nstepr-1)	r01[ibin][itheta]=r[ibin][ir];
				else if(ir==0)	{cout<<"ibin="<<ibin<<"\t"<<"itheta="<<itheta<<"\tminimum lies on ir = 0, please select proper range!"<<endl;	continue;}
				else 	{cout<<"ibin="<<ibin<<"\t"<<"itheta="<<itheta<<"\tminimum lies on ir = maximum "<<nstepr-1<<", please select proper range!"<<endl;	continue;}
				avgmult[ibin]=1.0*totmultall[ibin]/Nevent[ibin];
				avgpt[ibin]=1.0*totptall[ibin]/totmultall[ibin];
				if(isSimple==0)	V[ibin][itheta]=Vmax[ibin]-ir*eps[ibin]+eps[ibin]*(Gmod2[ibin][itheta][ir+1]-Gmod2[ibin][itheta][ir-1])/2./(Gmod2[ibin][itheta][ir-1]-2*Gmod2[ibin][itheta][ir]+Gmod2[ibin][itheta][ir+1]);
				else V[ibin][itheta]=j01/r0[ibin][itheta]; //simple method
				r0[ibin][itheta]=j01/V[ibin][itheta];
				V[ibin][itheta]/=avgmult[ibin];
				//sigma2[ibin][itheta]=Q2[ibin]/Nevent[ibin]-(Qx1[ibin]/Nevent[ibin])*(Qx1[ibin]/Nevent[ibin])-(Qy1[ibin]/Nevent[ibin])*(Qy1[ibin]/Nevent[ibin])-(V[ibin][itheta]*avgmult[ibin])*(V[ibin][itheta]*avgmult[ibin]);
				sigma2[ibin][itheta]=Q2[ibin]/Nevent[ibin]-(Qx1[ibin]/Nevent[ibin])*(Qx1[ibin]/Nevent[ibin])-(Qy1[ibin]/Nevent[ibin])*(Qy1[ibin]/Nevent[ibin]);
				sigma2_[ibin]+=sigma2[ibin][itheta];
				Vmean[ibin]+=V[ibin][itheta];
				chi[ibin][itheta]=V[ibin][itheta]*avgmult[ibin]/TMath::Sqrt(sigma2[ibin][itheta]);
				//deltaV[ibin][itheta]=V[ibin][itheta]/j01/TMath::BesselJ1(j01)*TMath::Sqrt((TMath::Exp(j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta])+TMath::Exp(-j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta])*TMath::BesselJ0(2*j01))/2./Nevent[ibin]);
				}
			sigma2_[ibin]/=ntheta;
			Vmean[ibin]/=ntheta;
			sigma2_[ibin]-=TMath::Power(Vmean[ibin]*avgmult[ibin],2);
			chi_[ibin]=Vmean[ibin]*avgmult[ibin]/TMath::Sqrt(sigma2_[ibin]);
			//deltaVmean[ibin]+=TMath::Exp(j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.));
			for(int itheta=0;itheta<ntheta;itheta++){
			deltaV[ibin][itheta]=V[ibin][itheta]/j01/TMath::BesselJ1(j01)*TMath::Sqrt((TMath::Exp(j01*j01/2./chi_[ibin]/chi_[ibin])+TMath::Exp(-j01*j01/2./chi_[ibin]/chi_[ibin])*TMath::BesselJ0(2*j01))/2./Nevent[ibin]);
			deltaVmean[ibin]+=TMath::Exp(j01*j01/2./chi_[ibin]/chi_[ibin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi_[ibin]/chi_[ibin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.));
			fstrV<<ibin<<"\t"<<itheta<<"\t"<<r0[ibin][itheta]<<"\t"<<V[ibin][itheta]<<"\t"<<sigma2[ibin][itheta]<<"\t"<<chi[ibin][itheta]<<"\t"<<deltaV[ibin][itheta]<<endl;
			}
			fstrV<<endl;
			deltaVmean[ibin]=Vmean[ibin]/j01/TMath::BesselJ1(j01)*TMath::Sqrt(deltaVmean[ibin]/ntheta/2./Nevent[ibin]);
		}
		fstrV<<endl;

	fstrV<<"ibin"<<"\t"<<"avgmult"<<"\t"<<"avgpt"<<"\t"<<"Vn mean"<<"\t"<<"Vn mean Error"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
			fstrV<<ibin<<"\t"<<avgmult[ibin]<<"\t"<<avgpt[ibin]<<"\t"<<Vmean[ibin]<<"\t"<<deltaVmean[ibin]<<endl;
		}
		fstrV<<endl;
	
        TDirectory *dirsample = outf->mkdir(Form("s_%d",isample));
        dirsample->cd();
        Nevent.Write("Nevent");
       	totmultall.Write("totmultall");
		Qx1.Write("Qx1");
                Qy1.Write("Qy1");
                Q2.Write("Q2");
		totptall.Write("totptall");
		avgmult.Write("avgmult");
		avgpt.Write("avgpt");
		Vmean.Write("Vmean");
		deltaVmean.Write("deltaVmean");
		chi_.Write("chi");
        
	for(int ibin=0;ibin<nbin;ibin++){
                TDirectory *dir0 = dirsample->mkdir(Form("D_%d",ibin));dir0->cd();
		r[ibin].Write("r");
		Vmeanmean[ibin]+=Vmean[ibin]/nsamples;
		deltaVmeanmean[ibin]+=deltaVmean[ibin]/nsamples;
		sigmaVmeanmean[ibin]+=TMath::Power(Vmean[ibin]/nsamples,2);
			sigma2[ibin].Write("sigma2");	chi[ibin].Write("chi0");	deltaV[ibin].Write("deltaV");
			r0[ibin].Write("r0");	r01[ibin].Write("r01");	V[ibin].Write("V");

        		for(int itheta=0;itheta<ntheta;itheta++){
                        	TDirectory *dir1 = dir0->mkdir(Form("D_%d",itheta));dir1->cd();
                        	GRe[ibin][itheta].Write(Form("GRe"));
                        	GIm[ibin][itheta].Write(Form("GIm"));
                        	Gmod2[ibin][itheta].Write(Form("G2"));
                	}
        	}
        }
	fstrV<<"ibin"<<"\t"<<"Vn mean mean"<<"\t"<<"Vn mean mean Error"<<"\t"<<"Vn mean deviation"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		sigmaVmeanmean[ibin]=TMath::Sqrt(sigmaVmeanmean[ibin]*nsamples-Vmeanmean[ibin]*Vmeanmean[ibin]);
		fstrV<<endl<<ibin<<"\t"<<Vmeanmean[ibin]<<"\t"<<deltaVmeanmean[ibin]<<"\t"<<sigmaVmeanmean[ibin]<<endl;
	}
	outf->cd();
	Vmeanmean.Write("Vmeanmean");
	deltaVmeanmean.Write("deltaVmeanmean");
	sigmaVmeanmean.Write("sigmaVmeanmean");	
	
	cout<<endl;
	
	outf->Close();
}
