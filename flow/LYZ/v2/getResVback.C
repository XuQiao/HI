#include <iomanip>
#include <fstream>
#include "par.h"

void getResV(){

//------------------------Deal with the number and get the output stored in txt files and root files
	
	string SumorProd = getenv("SUMORPROD");
	string dir = getenv("DIR");
        double Vmax[nbin][netaV+1], eps[nbin][netaV+1];
        for(int ibin=0; ibin<nbin ;ibin++){
            for(int ietabin=0;ietabin<netaV; ietabin++){
                Vmax[ibin][ietabin]=0.065*(trkbin[ibin]+30)/netaV;
                eps[ibin][ietabin]=0.00025*(trkbin[ibin]+30)/netaV;
                }
            Vmax[ibin][netaV]=0.065*(trkbin[ibin]+30);
            eps[ibin][netaV]=0.00025*(trkbin[ibin]+30);
        }
	ofstream  fstrV;
	double theta[ntheta];
	TVectorD avgmult[nbin], avgpt[nbin], totmult[nbin], totpt[nbin];
	TVectorD Qx1[nbin], Qy1[nbin], Q2[nbin];
	TVectorD Gmod2[nbin][netaV+1][ntheta];
	TVectorD sigma2[nbin][netaV+1],deltaV[nbin][netaV+1];
	TVectorD sigma2_[nbin],chi_[nbin];
	TVectorD deltaVmean[nbin], Vmean[nbin];
	TVectorD r[nbin][netaV+1];
	TVectorD r0[nbin][netaV+1], r01[nbin][netaV+1], V[nbin][netaV+1], chi[nbin][netaV+1];
	TVectorD GRe[nbin][netaV+1][ntheta]; TVectorD* GRe_t[nbin][netaV+1][ntheta];
	TVectorD GIm[nbin][netaV+1][ntheta]; TVectorD* GIm_t[nbin][netaV+1][ntheta];
	TComplex G[nbin][netaV+1][ntheta][nstepr];
	if(SumorProd=="Sum")	fstrV.open("V_Sum.txt");
	else	fstrV.open("V_Prod.txt"); 
	TFile *f[nFileAll];

        for(int ibin=0; ibin<nbin; ibin++){
		for(int ietabin=0;ietabin<=netaV;ietabin++){
		    r[ibin][ietabin].ResizeTo(nstepr);
		    r[ibin][ietabin].Zero();
                    for(int ir=0; ir<nstepr; ir++){
                        if(isSimple==0)  r[ibin][ietabin][ir]=j01/(Vmax[ibin][ietabin]-eps[ibin][ietabin]*ir);	
	                else		 r[ibin][ietabin][ir]=0.00025*20*(ir+1);
		    }
	    }
        }
                
	TVectorD Nevent;	Nevent.ResizeTo(nbin);  Nevent.Zero();
        TVectorD totmultall;	totmultall.ResizeTo(nbin);      totmultall.Zero();
        TVectorD avgmultall;	avgmultall.ResizeTo(nbin);      avgmultall.Zero();
        TVectorD tottrk;	tottrk.ResizeTo(nbin);      tottrk.Zero();
        TVectorD avgtrk;	avgtrk.ResizeTo(nbin);      avgtrk.Zero();
	
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;

	for(int ibin=0;ibin<nbin;ibin++){
		for(int ietabin=0;ietabin<=netaV;ietabin++){
			r0[ibin][ietabin].ResizeTo(ntheta);
			r01[ibin][ietabin].ResizeTo(ntheta);
			sigma2[ibin][ietabin].ResizeTo(ntheta);
			V[ibin][ietabin].ResizeTo(ntheta);
			deltaV[ibin][ietabin].ResizeTo(ntheta);
			chi[ibin][ietabin].ResizeTo(ntheta);
                        for(int itheta=0;itheta<ntheta;itheta++){
				Gmod2[ibin][ietabin][itheta].ResizeTo(nstepr);
				GRe[ibin][ietabin][itheta].ResizeTo(nstepr);
				GRe[ibin][ietabin][itheta].Zero();
				GIm[ibin][ietabin][itheta].ResizeTo(nstepr);
				GIm[ibin][ietabin][itheta].Zero();
			}
		}

		avgmult[ibin].ResizeTo(netaV+1);	deltaVmean[ibin].ResizeTo(netaV+1);	Vmean[ibin].ResizeTo(netaV+1);
		totpt[ibin].ResizeTo(netaV+1);	totpt[ibin].Zero();
		totmult[ibin].ResizeTo(netaV+1);	totmult[ibin].Zero();
		avgpt[ibin].ResizeTo(netaV+1);
                Qx1[ibin].ResizeTo(netaV+1);     Qx1[ibin].Zero();
                Qy1[ibin].ResizeTo(netaV+1);     Qy1[ibin].Zero();
                Q2[ibin].ResizeTo(netaV+1);     Q2[ibin].Zero();
		sigma2_[ibin].ResizeTo(netaV+1);sigma2_[ibin].Zero();
		chi_[ibin].ResizeTo(netaV+1);chi_[ibin].Zero();
	}

        for(int ifile=0; ifile<nFileAll; ifile++){
		if(SumorProd=="Sum") f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/etatracknormcpt03to6/%s/AnaV_Sum_%d.root",dir.c_str(),ifile));
		else f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/etatracknormcpt03to6/%s/AnaV_Prod_%d.root",dir.c_str(),ifile));
		TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
		TVectorD* totmultall_t =  (TVectorD*)f[ifile]->Get(Form("totmultall"));
		TVectorD* tottrk_t =  (TVectorD*)f[ifile]->Get(Form("tottrk"));
		for(int ibin=0;ibin<nbin;ibin++){
			for(int ietabin=0;ietabin<=netaV;ietabin++){
				for(int itheta=0;itheta<ntheta;itheta++){
				    GRe_t[ibin][ietabin][itheta] = (TVectorD*)f[ifile]->Get(Form("GRe_%d_%d_%d",ibin,ietabin,itheta));
				    GIm_t[ibin][ietabin][itheta] = (TVectorD*)f[ifile]->Get(Form("GIm_%d_%d_%d",ibin,ietabin,itheta));
				    for(ir=0; ir<nstepr; ir++){
				        GRe[ibin][ietabin][itheta][ir] += (*GRe_t[ibin][ietabin][itheta])[ir];
					GIm[ibin][ietabin][itheta][ir] += (*GIm_t[ibin][ietabin][itheta])[ir];
				    }
				}
			}
			TVectorD* Qx1_t =  (TVectorD*)f[ifile]->Get(Form("Qx1_%d",ibin));
			TVectorD* Qy1_t =  (TVectorD*)f[ifile]->Get(Form("Qy1_%d",ibin));
			TVectorD* Q2_t =  (TVectorD*)f[ifile]->Get(Form("Q2_%d",ibin));
			TVectorD* totmult_t =  (TVectorD*)f[ifile]->Get(Form("totmult_%d",ibin));
			TVectorD* totpt_t =  (TVectorD*)f[ifile]->Get(Form("totpt_%d",ibin));
			for(int ietabin=0;ietabin<=netaV;ietabin++){
				Qx1[ibin][ietabin] += (*Qx1_t)[ietabin];
				Qy1[ibin][ietabin] += (*Qy1_t)[ietabin];
				Q2[ibin][ietabin] += (*Q2_t)[ietabin];
				totmult[ibin][ietabin] += (*totmult_t)[ietabin];
				totpt[ibin][ietabin] += (*totpt_t)[ietabin];
			}
			Nevent[ibin] += (*Nevent_t)[ibin];
			totmultall[ibin] += (*totmultall_t)[ibin];	
			tottrk[ibin] += (*tottrk_t)[ibin];	
		}
		f[ifile]->Close();
	}
	
	fstrV<<setprecision(4)<<fixed;
	fstrV<<"ibin"<<"\t"<<"ietabin"<<"\t"<<"itheta"<<"\t"<<"r0"<<"\t"<<"V"<<"\t"<<"sigma2"<<"\t"<<"chi"<<"\t"<<"Vn Errors"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		avgmultall[ibin]=totmultall[ibin]/Nevent[ibin];
		avgtrk[ibin]=tottrk[ibin]/Nevent[ibin];
		for(int ietabin=0;ietabin<=netaV;ietabin++){
			for(int itheta=0;itheta<ntheta;itheta++){
				for(ir=0; ir<nstepr; ir++){
					G[ibin][ietabin][itheta][ir]=TComplex(GRe[ibin][ietabin][itheta][ir],GIm[ibin][ietabin][itheta][ir]);	
					G[ibin][ietabin][itheta][ir]/=Nevent[ibin];
					Gmod2[ibin][ietabin][itheta][ir]=TMath::Power(TComplex::Abs(G[ibin][ietabin][itheta][ir]),2);
				}
				for(ir=0; ir<nstepr-1; ir++)
					if(ir!=0 && Gmod2[ibin][ietabin][itheta][ir]<=Gmod2[ibin][ietabin][itheta][ir-1] && Gmod2[ibin][ietabin][itheta][ir]<=Gmod2[ibin][ietabin][itheta][ir+1]) break;
				if(ir!=0 && ir<nstepr-1)	r01[ibin][ietabin][itheta]=r[ibin][ietabin][ir];
				else if(ir==0)	{cout<<"ibin="<<ibin<<"\t"<<"ietabin="<<ietabin<<"\t"<<"itheta="<<itheta<<"\tminimum lies on ir = 0, please select proper range!"<<endl;	continue;}
				else 	{cout<<"ibin="<<ibin<<"\t"<<"ietabin="<<ietabin<<"\t"<<"itheta="<<itheta<<"\tminimum lies on ir = maximum "<<nstepr-1<<", please select proper range!"<<endl;	continue;}
				avgmult[ibin][ietabin]=1.0*totmult[ibin][ietabin]/Nevent[ibin];
				avgpt[ibin][ietabin]=1.0*totpt[ibin][ietabin]/totmult[ibin][ietabin];
				if(isSimple==0)	V[ibin][ietabin][itheta]=Vmax[ibin][ietabin]-ir*eps[ibin][ietabin]+eps[ibin][ietabin]*(Gmod2[ibin][ietabin][itheta][ir+1]-Gmod2[ibin][ietabin][itheta][ir-1])/2./(Gmod2[ibin][ietabin][itheta][ir-1]-2*Gmod2[ibin][ietabin][itheta][ir]+Gmod2[ibin][ietabin][itheta][ir+1]);
				else V[ibin][ietabin][itheta]=j01/r01[ibin][ietabin][itheta]; //simple method
				r0[ibin][ietabin][itheta]=j01/V[ibin][ietabin][itheta];
				V[ibin][ietabin][itheta]/=avgmult[ibin][ietabin];
				sigma2[ibin][ietabin][itheta]=Q2[ibin][ietabin]/Nevent[ibin]-(Qx1[ibin][ietabin]/Nevent[ibin])*(Qx1[ibin][ietabin]/Nevent[ibin])-(Qy1[ibin][ietabin]/Nevent[ibin])*(Qy1[ibin][ietabin]/Nevent[ibin]);
				sigma2_[ibin][ietabin]+=sigma2[ibin][ietabin][itheta];
				Vmean[ibin][ietabin]+=V[ibin][ietabin][itheta];
				chi[ibin][ietabin][itheta]=V[ibin][ietabin][itheta]*avgmult[ibin][ietabin]/TMath::Sqrt(sigma2[ibin][ietabin][itheta]);
				}
			//deltaVmean[ibin][ietabin]+=TMath::Exp(j01*j01/2./chi[ibin][ietabin][itheta]/chi[ibin][ietabin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi[ibin][ietabin][itheta]/chi[ibin][ietabin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.));
			sigma2_[ibin][ietabin]/=ntheta;
                        Vmean[ibin][ietabin]/=ntheta;
                        sigma2_[ibin][ietabin]-=TMath::Power(Vmean[ibin][ietabin]*avgmult[ibin][ietabin],2);
                        chi_[ibin][ietabin]=Vmean[ibin][ietabin]*avgmult[ibin][ietabin]/TMath::Sqrt(sigma2_[ibin][ietabin]);
                        for(int itheta=0;itheta<ntheta;itheta++){
				deltaV[ibin][ietabin][itheta]=V[ibin][ietabin][itheta]/j01/TMath::BesselJ1(j01)*TMath::Sqrt((TMath::Exp(j01*j01/2./chi_[ibin][ietabin]/chi_[ibin][ietabin])+TMath::Exp(-j01*j01/2./chi_[ibin][ietabin]/chi_[ibin][ietabin])*TMath::BesselJ0(2*j01))/2./Nevent[ibin]);
                        	deltaVmean[ibin][ietabin]+=TMath::Exp(j01*j01/2./chi_[ibin][ietabin]/chi_[ibin][ietabin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi_[ibin][ietabin]/chi_[ibin][ietabin]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.));
				if(ietabin==netaV) fstrV<<ibin<<"\t"<<ietabin<<"\t"<<itheta<<"\t"<<r0[ibin][ietabin][itheta]<<"\t"<<V[ibin][ietabin][itheta]<<"\t"<<sigma2[ibin][ietabin][itheta]<<"\t"<<chi[ibin][ietabin][itheta]<<"\t"<<deltaV[ibin][ietabin][itheta]<<endl;
			}
			deltaVmean[ibin][ietabin]=Vmean[ibin][ietabin]/j01/TMath::BesselJ1(j01)*TMath::Sqrt(deltaVmean[ibin][ietabin]/ntheta/2./Nevent[ibin]);
		}
		fstrV<<endl;
	}

	fstrV<<"ibin"<<"\t"<<"ietabin"<<"\t"<<"avgmult"<<"\t"<<"avgpt"<<"\t"<<"Vn mean"<<"\t"<<"Vn mean Error"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		for(int ietabin=0;ietabin<=netaV;ietabin++){
			fstrV<<ibin<<"\t"<<ietabin<<"\t"<<avgmult[ibin][ietabin]<<"\t"<<avgpt[ibin][ietabin]<<"\t"<<Vmean[ibin][ietabin]<<"\t"<<deltaVmean[ibin][ietabin]<<endl;
		}
		fstrV<<endl;
	}
	
	fstrV<<"ibin"<<"\t"<<"Nevent"<<"\t"<<"totmultall"<<"\t"<<"avgmultall"<<"\t"<<"avgtrk"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		fstrV<<ibin<<"\t"<<Nevent[ibin]<<"\t"<<totmultall[ibin]<<"\t"<<avgmultall[ibin]<<"\t"<<avgtrk[ibin]<<endl;		
		fstrV<<endl;
	}

	if(SumorProd=="Sum")	TString outname = "mergedV_Sum.root";
	else 	TString outname="mergedV_Prod.root";
	TFile *outf = new TFile(outname,"Recreate");
        Nevent.Write("Nevent");
       	totmultall.Write("totmultall");
       	avgmultall.Write("avgmultall");
       	tottrk.Write("tottrk");
       	avgtrk.Write("avgtrk");
        
	for(int ibin=0;ibin<nbin;ibin++){
                TDirectory *dir0 = outf->mkdir(Form("D_%d",ibin,itheta));dir0->cd();
		Qx1[ibin].Write("Qx1");
                Qy1[ibin].Write("Qy1");
                Q2[ibin].Write("Q2");
		totmult[ibin].Write("totmult");
		totpt[ibin].Write("totpt");
		avgmult[ibin].Write("avgmult");
		avgpt[ibin].Write("avgpt");
		Vmean[ibin].Write("Vmean");
		deltaVmean[ibin].Write("deltaVmean");
		chi_[ibin].Write("chi");

		for(int ietabin=0;ietabin<=netaV;ietabin++){
                        TDirectory *dir1 = dir0->mkdir(Form("D_%d",ietabin));dir1->cd();
		        r[ibin][ietabin].Write("r");
			sigma2[ibin][ietabin].Write("sigma2");	chi[ibin][ietabin].Write("chi0");	deltaV[ibin][ietabin].Write("deltaV");
			r0[ibin][ietabin].Write("r0");	r01[ibin][ietabin].Write("r01");	V[ibin][ietabin].Write("V");

        		for(int itheta=0;itheta<ntheta;itheta++){
                        	TDirectory *dir2 = dir1->mkdir(Form("D_%d",itheta));dir2->cd();
                        	GRe[ibin][ietabin][itheta].Write(Form("GRe"));
                        	GIm[ibin][ietabin][itheta].Write(Form("GIm"));
                        	Gmod2[ibin][ietabin][itheta].Write(Form("G2"));
                	}
        	}
	}
	outf->Close();
}
