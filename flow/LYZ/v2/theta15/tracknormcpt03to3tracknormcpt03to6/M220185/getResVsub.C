#include <iomanip>
#include <fstream>
#include "par.h"

void getResVsub(){
	string SumorProd = getenv("SUMORPROD");

//------------------------Deal with the number and get the output stored in txt files and root files
	double Vmax[nbin], eps[nbin];
	string dir = getenv("DIR");
	for(int ibin=0; ibin<nbin ;ibin++){
        	Vmax[ibin]=0.065*(trkbin[ibin]+30);
        	eps[ibin]=0.00025*(trkbin[ibin]+30);
	}
	ofstream  fstrV;
	double theta[ntheta];
	TVectorD avgmult[nbin], avgpt[nbin], totmult[nbin], totpt[nbin];
	TVectorD Qx1[nbin], Qy1[nbin], Q2[nbin];
	TVectorD Gmod2[nbin][nptV][ntheta];
	TVectorD sigma2[nbin][nptV],deltaV[nbin][nptV];
	TVectorD sigma2_[nbin],chi_[nbin];
	TVectorD deltaVmean[nbin], Vmean[nbin];
	TVectorD Vmeanmean, deltaVmeanmean, sigmaVmeanmean;
	TVectorD r[nbin];
	TVectorD r0[nbin][nptV], r01[nbin][nptV], V[nbin][nptV], chi[nbin][nptV];
	TVectorD GRe[nbin][nptV][ntheta]; TVectorD* GRe_t[nbin][nptV][ntheta];
	TVectorD GIm[nbin][nptV][ntheta]; TVectorD* GIm_t[nbin][nptV][ntheta];
	TComplex G[nbin][nptV][ntheta][nstepr];
	if(SumorProd=="Sum")	fstrV.open("V_Sum_sub.txt");
	else	fstrV.open("V_Prod_sub.txt");
	if(SumorProd=="Sum")       TString outname = "mergedV_Sum_sub.root";
        else    TString outname="mergedV_Prod_sub.root";
      	TFile *outf = new TFile(outname,"Recreate");
	TFile *f[nFileAll];

        for(int ibin=0; ibin<nbin; ibin++){
		r[ibin].ResizeTo(nstepr);
                for(int ir=0; ir<nstepr; ir++){
        if(isSimple==0)  r[ibin][ir]=j01/(Vmax[ibin]-eps[ibin]*ir);	
	else		 r[ibin][ir]=0.00025*20*(ir+1);
		}
	}
         
	TVectorD Nevent_;	Nevent_.ResizeTo(nbin);  Nevent_.Zero();
	for(int ifile=0; ifile<nFileAll; ifile++){
	if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/%s/AnaV_Sum_%d.root",dir.c_str(),ifile));
        else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/%s/AnaV_Prod_%d.root",dir.c_str(),ifile));
        TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
	for(int ibin=0;ibin<nbin;ibin++){
		Nevent_[ibin]+=(*Nevent_t)[ibin];
	}
	f[ifile]->Close();
	}
	for(int ibin=0;ibin<nbin;ibin++)
	cout<<Nevent_[ibin]<<endl;	
      
	Vmeanmean.ResizeTo(nbin);  Vmeanmean.Zero();
	deltaVmeanmean.ResizeTo(nbin);  deltaVmeanmean.Zero();
	sigmaVmeanmean.ResizeTo(nbin);  sigmaVmeanmean.Zero();

	for(int isample=0;isample<nsamples;isample++){
	TVectorD Nevent;	Nevent.ResizeTo(nbin);  Nevent.Zero();
	TVectorD Nevent0;	Nevent0.ResizeTo(nbin);  Nevent0.Zero();
        TVectorD totmultall;	totmultall.ResizeTo(nbin);      totmultall.Zero();
	
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;

	for(int ibin=0;ibin<nbin;ibin++){
		for(int iptbin=0;iptbin<nptV;iptbin++){
			r0[ibin][iptbin].ResizeTo(ntheta);
			r01[ibin][iptbin].ResizeTo(ntheta);
			sigma2[ibin][iptbin].ResizeTo(ntheta);
			V[ibin][iptbin].ResizeTo(ntheta);
			deltaV[ibin][iptbin].ResizeTo(ntheta);
			chi[ibin][iptbin].ResizeTo(ntheta);
                        for(int itheta=0;itheta<ntheta;itheta++){
				Gmod2[ibin][iptbin][itheta].ResizeTo(nstepr);
				GRe[ibin][iptbin][itheta].ResizeTo(nstepr);
				GRe[ibin][iptbin][itheta].Zero();
				GIm[ibin][iptbin][itheta].ResizeTo(nstepr);
				GIm[ibin][iptbin][itheta].Zero();
			}
		}

		avgmult[ibin].ResizeTo(nptV);	deltaVmean[ibin].ResizeTo(nptV);	Vmean[ibin].ResizeTo(nptV);
		deltaVmean[ibin].Zero();	Vmean[ibin].Zero();
		Vmean[ibin].Zero();	deltaVmean[ibin].Zero();
		totpt[ibin].ResizeTo(nptV);	totpt[ibin].Zero();
		totmult[ibin].ResizeTo(nptV);	totmult[ibin].Zero();
		avgpt[ibin].ResizeTo(nptV);
                Qx1[ibin].ResizeTo(nptV);     Qx1[ibin].Zero();
                Qy1[ibin].ResizeTo(nptV);     Qy1[ibin].Zero();
                Q2[ibin].ResizeTo(nptV);     Q2[ibin].Zero();
		sigma2_[ibin].ResizeTo(nptV);sigma2_[ibin].Zero();
		chi_[ibin].ResizeTo(nptV);chi_[ibin].Zero();
	}

	for(int ifile=0; ifile<nFileAll; ifile++){
		if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/%s/AnaV_Sum_%d.root",dir.c_str(),ifile));
		else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/%s/AnaV_Prod_%d.root",dir.c_str(),ifile));
		TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
		TVectorD* totmultall_t =  (TVectorD*)f[ifile]->Get(Form("totmultall"));
		for(int ibin=0;ibin<nbin;ibin++){
			Nevent0[ibin] += (*Nevent_t)[ibin];
			double start=isample*Nevent_[ibin]/nsamples;
			double end=(isample+1)*Nevent_[ibin]/nsamples;
			if(Nevent0[ibin]>start && Nevent0[ibin]<=end){
			Nevent[ibin] += (*Nevent_t)[ibin];
			for(int iptbin=0;iptbin<nptV;iptbin++){
					for(int itheta=0;itheta<ntheta;itheta++){
					GRe_t[ibin][iptbin][itheta] = (TVectorD*)f[ifile]->Get(Form("GRe_%d_%d_%d",ibin,iptbin,itheta));
					GIm_t[ibin][iptbin][itheta] = (TVectorD*)f[ifile]->Get(Form("GIm_%d_%d_%d",ibin,iptbin,itheta));
					for(ir=0; ir<nstepr; ir++){
						GRe[ibin][iptbin][itheta][ir] += (*GRe_t[ibin][iptbin][itheta])[ir];
						GIm[ibin][iptbin][itheta][ir] += (*GIm_t[ibin][iptbin][itheta])[ir];
					}
				}
			}
			TVectorD* Qx1_t =  (TVectorD*)f[ifile]->Get(Form("Qx1_%d",ibin));
			TVectorD* Qy1_t =  (TVectorD*)f[ifile]->Get(Form("Qy1_%d",ibin));
			TVectorD* Q2_t =  (TVectorD*)f[ifile]->Get(Form("Q2_%d",ibin));
			TVectorD* totmult_t =  (TVectorD*)f[ifile]->Get(Form("totmult_%d",ibin));
			TVectorD* totpt_t =  (TVectorD*)f[ifile]->Get(Form("totpt_%d",ibin));
			for(int iptbin=0;iptbin<nptV;iptbin++){
				Qx1[ibin][iptbin] += (*Qx1_t)[iptbin];
				Qy1[ibin][iptbin] += (*Qy1_t)[iptbin];
				Q2[ibin][iptbin] += (*Q2_t)[iptbin];
				totmult[ibin][iptbin] += (*totmult_t)[iptbin];
				totpt[ibin][iptbin] += (*totpt_t)[iptbin];
			}
			totmultall[ibin] += (*totmultall_t)[ibin];	
		}
		}
		f[ifile]->Close();
	}
		
	fstrV<<setprecision(4)<<fixed;
	fstrV<<"ibin"<<"\t"<<"iptbin"<<"\t"<<"itheta"<<"\t"<<"r0"<<"\t"<<"V"<<"\t"<<"sigma2"<<"\t"<<"chi"<<"\t"<<"Vn Errors"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		for(int iptbin=0;iptbin<nptV;iptbin++){
			for(int itheta=0;itheta<ntheta;itheta++){
				for(ir=0; ir<nstepr; ir++){
					G[ibin][iptbin][itheta][ir]=TComplex(GRe[ibin][iptbin][itheta][ir],GIm[ibin][iptbin][itheta][ir]);	
					G[ibin][iptbin][itheta][ir]/=Nevent[ibin];
					Gmod2[ibin][iptbin][itheta][ir]=TMath::Power(TComplex::Abs(G[ibin][iptbin][itheta][ir]),2);
				}
				for(ir=0; ir<nstepr-1; ir++)
					if(ir!=0 && Gmod2[ibin][iptbin][itheta][ir]<=Gmod2[ibin][iptbin][itheta][ir-1] && Gmod2[ibin][iptbin][itheta][ir]<=Gmod2[ibin][iptbin][itheta][ir+1]) break;
				if(ir!=0 && ir<nstepr-1)	r01[ibin][iptbin][itheta]=r[ibin][ir];
				else if(ir==0)	{cout<<"ibin="<<ibin<<"\t"<<"iptbin="<<iptbin<<"\t"<<"itheta="<<itheta<<"\tminimum lies on ir = 0, please select proper range!"<<endl;	continue;}
				else 	{cout<<"ibin="<<ibin<<"\t"<<"iptbin="<<iptbin<<"\t"<<"itheta="<<itheta<<"\tminimum lies on ir = maximum "<<nstepr-1<<", please select proper range!"<<endl;	continue;}
				avgmult[ibin][iptbin]=1.0*totmult[ibin][iptbin]/Nevent[ibin];
				avgpt[ibin][iptbin]=1.0*totpt[ibin][iptbin]/totmult[ibin][iptbin];
				if(isSimple==0)	V[ibin][iptbin][itheta]=Vmax[ibin]-ir*eps[ibin]+eps[ibin]*(Gmod2[ibin][iptbin][itheta][ir+1]-Gmod2[ibin][iptbin][itheta][ir-1])/2./(Gmod2[ibin][iptbin][itheta][ir-1]-2*Gmod2[ibin][iptbin][itheta][ir]+Gmod2[ibin][iptbin][itheta][ir+1]);
				else V[ibin][iptbin][itheta]=j01/r0[ibin][iptbin][itheta]; //simple method
				r0[ibin][iptbin][itheta]=j01/V[ibin][iptbin][itheta];
				V[ibin][iptbin][itheta]/=avgmult[ibin][iptbin];
				//sigma2[ibin][iptbin][itheta]=Q2[ibin][iptbin]/Nevent[ibin]-(Qx1[ibin][iptbin]/Nevent[ibin])*(Qx1[ibin][iptbin]/Nevent[ibin])-(Qy1[ibin][iptbin]/Nevent[ibin])*(Qy1[ibin][iptbin]/Nevent[ibin])-(V[ibin][iptbin][itheta]*avgmult[ibin][iptbin])*(V[ibin][iptbin][itheta]*avgmult[ibin][iptbin]);
				sigma2[ibin][iptbin][itheta]=Q2[ibin][iptbin]/Nevent[ibin]-(Qx1[ibin][iptbin]/Nevent[ibin])*(Qx1[ibin][iptbin]/Nevent[ibin])-(Qy1[ibin][iptbin]/Nevent[ibin])*(Qy1[ibin][iptbin]/Nevent[ibin]);
				sigma2_[ibin][iptbin]+=sigma2[ibin][iptbin][itheta];
				Vmean[ibin][iptbin]+=V[ibin][iptbin][itheta];
				chi[ibin][iptbin][itheta]=V[ibin][iptbin][itheta]*avgmult[ibin][iptbin]/TMath::Sqrt(sigma2[ibin][iptbin][itheta]);
				//deltaV[ibin][iptbin][itheta]=V[ibin][iptbin][itheta]/j01/TMath::BesselJ1(j01)*TMath::Sqrt((TMath::Exp(j01*j01/2./chi[ibin][iptbin][itheta]/chi[ibin][iptbin][itheta])+TMath::Exp(-j01*j01/2./chi[ibin][iptbin][itheta]/chi[ibin][iptbin][itheta])*TMath::BesselJ0(2*j01))/2./Nevent[ibin]);
				}
			sigma2_[ibin][iptbin]/=ntheta;
			Vmean[ibin][iptbin]/=ntheta;
			sigma2_[ibin][iptbin]-=TMath::Power(Vmean[ibin][iptbin]*avgmult[ibin][iptbin],2);
			chi_[ibin][iptbin]=Vmean[ibin][iptbin]*avgmult[ibin][iptbin]/TMath::Sqrt(sigma2_[ibin][iptbin]);
			//deltaVmean[ibin][iptbin]+=TMath::Exp(j01*j01/2./chi[ibin][iptbin][itheta]/chi[ibin][iptbin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi[ibin][iptbin][itheta]/chi[ibin][iptbin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.));
			for(int itheta=0;itheta<ntheta;itheta++){
			deltaV[ibin][iptbin][itheta]=V[ibin][iptbin][itheta]/j01/TMath::BesselJ1(j01)*TMath::Sqrt((TMath::Exp(j01*j01/2./chi_[ibin][iptbin]/chi_[ibin][iptbin])+TMath::Exp(-j01*j01/2./chi_[ibin][iptbin]/chi_[ibin][iptbin])*TMath::BesselJ0(2*j01))/2./Nevent[ibin]);
			deltaVmean[ibin][iptbin]+=TMath::Exp(j01*j01/2./chi_[ibin][iptbin]/chi_[ibin][iptbin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi_[ibin][iptbin]/chi_[ibin][iptbin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.));
			fstrV<<ibin<<"\t"<<iptbin<<"\t"<<itheta<<"\t"<<r0[ibin][iptbin][itheta]<<"\t"<<V[ibin][iptbin][itheta]<<"\t"<<sigma2[ibin][iptbin][itheta]<<"\t"<<chi[ibin][iptbin][itheta]<<"\t"<<deltaV[ibin][iptbin][itheta]<<endl;
			}
			fstrV<<endl;
			deltaVmean[ibin][iptbin]=Vmean[ibin][iptbin]/j01/TMath::BesselJ1(j01)*TMath::Sqrt(deltaVmean[ibin][iptbin]/ntheta/2./Nevent[ibin]);
		}
		fstrV<<endl;
	}

	fstrV<<"ibin"<<"\t"<<"iptbin"<<"\t"<<"avgmult"<<"\t"<<"avgpt"<<"\t"<<"Vn mean"<<"\t"<<"Vn mean Error"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		for(int iptbin=0;iptbin<nptV;iptbin++){
			fstrV<<ibin<<"\t"<<iptbin<<"\t"<<avgmult[ibin][iptbin]<<"\t"<<avgpt[ibin][iptbin]<<"\t"<<Vmean[ibin][iptbin]<<"\t"<<deltaVmean[ibin][iptbin]<<endl;
		}
		fstrV<<endl;
	}
	
        TDirectory *dirsample = outf->mkdir(Form("s_%d",isample));
        dirsample->cd();
        Nevent.Write("Nevent");
       	totmultall.Write("totmultall");
        
	for(int ibin=0;ibin<nbin;ibin++){
                TDirectory *dir0 = dirsample->mkdir(Form("D_%d",ibin,itheta));dir0->cd();
		r[ibin].Write("r");
		Qx1[ibin].Write("Qx1");
                Qy1[ibin].Write("Qy1");
                Q2[ibin].Write("Q2");
		totmult[ibin].Write("totmult");
		totpt[ibin].Write("totpt");
		avgmult[ibin].Write("avgmult");
		avgpt[ibin].Write("avgpt");
		Vmean[ibin].Write("Vmean");
		deltaVmean[ibin].Write("deltaVmean");
		Vmeanmean[ibin]+=Vmean[ibin][0]/nsamples;
		deltaVmeanmean[ibin]+=deltaVmean[ibin][0]/nsamples;
		sigmaVmeanmean[ibin]+=TMath::Power(Vmean[ibin][0]/nsamples,2);
		chi_[ibin].Write("chi");

		for(int iptbin=0;iptbin<nptV;iptbin++){
                        TDirectory *dir1 = dir0->mkdir(Form("D_%d",iptbin));dir1->cd();
			sigma2[ibin][iptbin].Write("sigma2");	chi[ibin][iptbin].Write("chi0");	deltaV[ibin][iptbin].Write("deltaV");
			r0[ibin][iptbin].Write("r0");	r01[ibin][iptbin].Write("r01");	V[ibin][iptbin].Write("V");

        		for(int itheta=0;itheta<ntheta;itheta++){
                        	TDirectory *dir2 = dir1->mkdir(Form("D_%d",itheta));dir2->cd();
                        	GRe[ibin][iptbin][itheta].Write(Form("GRe"));
                        	GIm[ibin][iptbin][itheta].Write(Form("GIm"));
                        	Gmod2[ibin][iptbin][itheta].Write(Form("G2"));
                	}
        	}
	}
	
	}
	
	for(int ibin=0;ibin<nbin;ibin++){
                sigmaVmeanmean[ibin]=TMath::Sqrt(sigmaVmeanmean[ibin]*nsamples-Vmeanmean[ibin]*Vmeanmean[ibin])/TMath::Sqrt(nsamples);
        }
	outf->cd();
	Vmeanmean.Write("Vmeanmean");
	deltaVmeanmean.Write("deltaVmeanmean");
	sigmaVmeanmean.Write("sigmaVmeanmean");
	outf->Close();
}
