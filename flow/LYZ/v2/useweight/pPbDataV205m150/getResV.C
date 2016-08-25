#include <iomanip>
#include <iostream>
#include <fstream>
#include "par.h"

void getResV(){

//------------------------Deal with the number and get the output stored in txt files and root files
	
        double Vmax[nbin], eps[nbin];
        for(int ibin=0; ibin<nbin ;ibin++){
            Vmax[ibin]=0.065*(trkbin[ibin]+30);
            eps[ibin]=0.00025*(trkbin[ibin]+30);
        }
	ofstream  fstrV;
	double theta[ntheta];
	TVectorD Gmod2[nbin][ntheta];
	TVectorD sigma2[nbin],deltaV[nbin];
	TVectorD r[nbin];
	TVectorD r0[nbin], r01[nbin], V[nbin], chi[nbin];
	TVectorD GRe[nbin][ntheta]; TVectorD* GRe_t[nbin][ntheta];
	TVectorD GIm[nbin][ntheta]; TVectorD* GIm_t[nbin][ntheta];
	TComplex G[nbin][ntheta][nstepr];
	if(isSum)	fstrV.open("V_Sum.txt");
	else	fstrV.open("V_Prod.txt"); 
	TFile *f[nFileAll];

        for(int ibin=0; ibin<nbin; ibin++){
		    r[ibin].ResizeTo(nstepr);
                    for(int ir=0; ir<nstepr; ir++){
                        if(isSimple==0)  r[ibin][ir]=j01/(Vmax[ibin]-eps[ibin]*ir);	
	                else		 r[ibin][ir]=0.00025*20*(ir+1);
		    }
        }
                
	TVectorD Nevent;	Nevent.ResizeTo(nbin);  Nevent.Zero();
        TVectorD totmultall;	totmultall.ResizeTo(nbin);      totmultall.Zero();
        TVectorD totwall;	totwall.ResizeTo(nbin);      totwall.Zero();
        TVectorD avgmultall;	avgmultall.ResizeTo(nbin);      avgmultall.Zero();
        TVectorD tottrk;	tottrk.ResizeTo(nbin);      tottrk.Zero();
	TVectorD totptall;      totptall.ResizeTo(nbin);    totptall.Zero();
	TVectorD totetaall;      totetaall.ResizeTo(nbin);    totetaall.Zero();
        TVectorD avgtrk;	avgtrk.ResizeTo(nbin);      avgtrk.Zero();
	TVectorD avgmult;       avgmult.ResizeTo(nbin);
	TVectorD avgw;       avgw.ResizeTo(nbin);
        TVectorD deltaVmean;    deltaVmean.ResizeTo(nbin);	
        TVectorD Vmean;         Vmean.ResizeTo(nbin);
        TVectorD avgpt;         avgpt.ResizeTo(nbin);
        TVectorD avgeta;         avgeta.ResizeTo(nbin);
        TVectorD Qx1;           Qx1.ResizeTo(nbin); Qx1.Zero();
        TVectorD Qy1;           Qy1.ResizeTo(nbin); Qy1.Zero();
        TVectorD Q2;           Q2.ResizeTo(nbin); Q2.Zero();
        TVectorD sigma2_;       sigma2_.ResizeTo(nbin); sigma2_.Zero();
        TVectorD chi_;           chi_.ResizeTo(nbin); chi_.Zero();
	
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

        for(int ifile=0; ifile<nFileAll; ifile++){
		if(isSum) f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/useweight/pPbDataV205m150/AnaV_Sum_%d.root",ifile));
		else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/useweight/pPbDataV205m150/AnaV_Prod_%d.root",ifile));
		TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
		TVectorD* totmultall_t =  (TVectorD*)f[ifile]->Get(Form("totmultall"));
		TVectorD* totwall_t =  (TVectorD*)f[ifile]->Get(Form("totwall"));
		TVectorD* tottrk_t =  (TVectorD*)f[ifile]->Get(Form("tottrk"));
		TVectorD* totptall_t =  (TVectorD*)f[ifile]->Get(Form("totptall"));
		TVectorD* totetaall_t =  (TVectorD*)f[ifile]->Get(Form("totetaall"));
		TVectorD* Qx1_t =  (TVectorD*)f[ifile]->Get(Form("Qx1"));
	        TVectorD* Qy1_t =  (TVectorD*)f[ifile]->Get(Form("Qy1"));
		TVectorD* Q2_t =  (TVectorD*)f[ifile]->Get(Form("Q2"));
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
				totetaall[ibin] += (*totetaall_t)[ibin];
			        Nevent[ibin] += (*Nevent_t)[ibin];
			        totmultall[ibin] += (*totmultall_t)[ibin];	
			        totwall[ibin] += (*totwall_t)[ibin];	
			        tottrk[ibin] += (*tottrk_t)[ibin];	
		}
		f[ifile]->Close();
	}
	
	fstrV<<setprecision(4)<<fixed;
	fstrV<<"ibin"<<"\t"<<"itheta"<<"\t"<<"r0"<<"\t"<<"V"<<"\t"<<"sigma2"<<"\t"<<"chi"<<"\t"<<"Vn Errors"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		avgmultall[ibin]=totmultall[ibin]/Nevent[ibin];
		avgtrk[ibin]=tottrk[ibin]/Nevent[ibin];
		avgw[ibin]=totwall[ibin]/Nevent[ibin];
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
				avgeta[ibin]=1.0*totetaall[ibin]/totmultall[ibin];
				if(isSimple==0)	V[ibin][itheta]=Vmax[ibin]-ir*eps[ibin]+eps[ibin]*(Gmod2[ibin][itheta][ir+1]-Gmod2[ibin][itheta][ir-1])/2./(Gmod2[ibin][itheta][ir-1]-2*Gmod2[ibin][itheta][ir]+Gmod2[ibin][itheta][ir+1]);
				else V[ibin][itheta]=j01/r0[ibin][itheta]; //simple method
				r0[ibin][itheta]=j01/V[ibin][itheta];
				V[ibin][itheta]/=avgw[ibin];
				//V[ibin][itheta]/=avgmult[ibin]/(avgw[ibin]/avgmult[ibin]);
				sigma2[ibin][itheta]=Q2[ibin]/Nevent[ibin]-(Qx1[ibin]/Nevent[ibin])*(Qx1[ibin]/Nevent[ibin])-(Qy1[ibin]/Nevent[ibin])*(Qy1[ibin]/Nevent[ibin]);
				sigma2_[ibin]+=sigma2[ibin][itheta];
				Vmean[ibin]+=V[ibin][itheta];
				chi[ibin][itheta]=V[ibin][itheta]*avgmult[ibin]/TMath::Sqrt(sigma2[ibin][itheta]);
				}
			//deltaVmean[ibin]+=TMath::Exp(j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.));
			sigma2_[ibin]/=ntheta;
                        Vmean[ibin]/=ntheta;
                        sigma2_[ibin]-=TMath::Power(Vmean[ibin]*avgmult[ibin],2);
                        chi_[ibin]=Vmean[ibin]*avgmult[ibin]/TMath::Sqrt(sigma2_[ibin]);
                        for(int itheta=0;itheta<ntheta;itheta++){
				deltaV[ibin][itheta]=V[ibin][itheta]/j01/TMath::BesselJ1(j01)*TMath::Sqrt((TMath::Exp(j01*j01/2./chi_[ibin]/chi_[ibin])+TMath::Exp(-j01*j01/2./chi_[ibin]/chi_[ibin])*TMath::BesselJ0(2*j01))/2./Nevent[ibin]);
                        	deltaVmean[ibin]+=TMath::Exp(j01*j01/2./chi_[ibin]/chi_[ibin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi_[ibin]/chi_[ibin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.));
	        	  fstrV<<ibin<<"\t"<<itheta<<"\t"<<r0[ibin][itheta]<<"\t"<<V[ibin][itheta]<<"\t"<<sigma2[ibin][itheta]<<"\t"<<chi[ibin][itheta]<<"\t"<<deltaV[ibin][itheta]<<endl;
			}
                        fstrV<<"\t\t\t\t\t"<<sigma2_[ibin]<<"\t"<<chi_[ibin]<<endl;
			deltaVmean[ibin]=Vmean[ibin]/j01/TMath::BesselJ1(j01)*TMath::Sqrt(deltaVmean[ibin]/ntheta/2./Nevent[ibin]);
		fstrV<<endl;
	}

	fstrV<<"ibin"<<"\t"<<"avgmult"<<"\t"<<"avgpt"<<"\t"<<"avgeta"<<"\t"<<"Vn mean"<<"\t"<<"Vn mean Error"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
			fstrV<<ibin<<"\t"<<avgmult[ibin]<<"\t"<<avgpt[ibin]<<"\t"<<avgeta[ibin]<<"\t"<<Vmean[ibin]<<"\t"<<deltaVmean[ibin]<<endl;
		fstrV<<endl;
	}
	
	fstrV<<"ibin"<<"\t"<<"Nevent"<<"\t"<<"totmultall"<<"\t"<<"avgmultall"<<"\t"<<"avgtrk"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		fstrV<<ibin<<"\t"<<Nevent[ibin]<<"\t"<<totmultall[ibin]<<"\t"<<avgmultall[ibin]<<"\t"<<avgtrk[ibin]<<endl;		
		fstrV<<endl;
	}

	if(isSum)	TString outname = "mergedV_Sum.root";
	else 	TString outname="mergedV_Prod.root";
	TFile *outf = new TFile(outname,"Recreate");
        Nevent.Write("Nevent");
       	totmultall.Write("totmultall");
       	avgmultall.Write("avgmultall");
       	tottrk.Write("tottrk");
       	avgtrk.Write("avgtrk");
	Qx1.Write("Qx1");
        Qy1.Write("Qy1");
        Q2.Write("Q2");
        avgmult.Write("avgmult");
        avgpt.Write("avgpt");
        avgeta.Write("avgeta");
        Vmean.Write("Vmean");
        deltaVmean.Write("deltaVmean");
        chi_.Write("chi");
        sigma2_.Write("sigma2");
        
	for(int ibin=0;ibin<nbin;ibin++){
                TDirectory *dir0 = outf->mkdir(Form("D_%d",ibin));dir0->cd();
		r[ibin].Write("r");
		sigma2[ibin].Write("sigma2");	chi[ibin].Write("chi0");	deltaV[ibin].Write("deltaV");
		r0[ibin].Write("r0");	r01[ibin].Write("r01");	V[ibin].Write("V");

        		for(int itheta=0;itheta<ntheta;itheta++){
                        	TDirectory *dir1 = dir0->mkdir(Form("D_%d",itheta));dir1->cd();
                        	GRe[ibin][itheta].Write(Form("GRe"));
                        	GIm[ibin][itheta].Write(Form("GIm"));
                        	Gmod2[ibin][itheta].Write(Form("G2"));
                	}
        	}
	outf->Close();
}
