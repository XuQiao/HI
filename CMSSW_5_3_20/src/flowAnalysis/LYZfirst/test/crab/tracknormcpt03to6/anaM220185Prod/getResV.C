#include <iomanip>
#include <fstream>
#include <vector>
#include "par.h"

void getResV(){

//------------------------Deal with the number and get the output stored in txt files and root files
	int nFileAll=getnumber("filelist.dat");
	ofstream  fstrV;
	double theta[ntheta];
	TVectorD avgmult;
	TVectorD Gmod2[nbin][ntheta];
	TVectorD sigma2[nbin],deltaV[nbin];
	TVectorD deltaVmean, Vmean;
	TVectorD r[nbin];
	TVectorD r0[nbin], V[nbin], chi[nbin];
	TVectorD GRe[nbin][ntheta]; TVectorD* GRe_t[nbin][ntheta];
	TVectorD GIm[nbin][ntheta]; TVectorD* GIm_t[nbin][ntheta];
	TComplex G[nbin][ntheta][nstepr];
	fstrV.open("V.txt");
	TFile *f;
		
	double Vmax[nbin], eps[nbin];
        for(int ibin=0; ibin<nbin ;ibin++){
                Vmax[ibin]=0.065*(trkbin[ibin]+30);
                eps[ibin]=0.00025*(trkbin[ibin]+30);
        }
	
        for(int ibin=0; ibin<nbin; ibin++){
		r[ibin].ResizeTo(nstepr);
                for(int ir=0; ir<nstepr; ir++){
                        r[ibin][ir]=j01/(Vmax[ibin]-eps[ibin]*ir);	
	//		  r[ibin][ir]=0.00025*20*(ir+1);
		}
	}
	
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;

	for(int ibin=0;ibin<nbin;ibin++){
			r0[ibin].ResizeTo(ntheta);
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

                TVectorD Nevent;	Nevent.ResizeTo(nbin);  Nevent.Zero();
                TVectorD totmultall;	totmultall.ResizeTo(nbin);      totmultall.Zero();
		avgmult.ResizeTo(nbin);	deltaVmean.ResizeTo(nbin);	Vmean.ResizeTo(nbin);
                TVectorD Qx1; Qx1.ResizeTo(nbin);     Qx1.Zero();
                TVectorD Qy1; Qy1.ResizeTo(nbin);     Qy1.Zero();
                TVectorD Q2; Q2.ResizeTo(nbin);     Q2.Zero();

		for(int ifile=0; ifile<nFileAll; ifile++){
			f=TFile::Open(readline("filelist.dat",ifile).c_str());
			for(int ibin=0;ibin<nbin;ibin++){
				for(int itheta=0;itheta<ntheta;itheta++){
					GRe_t[ibin][itheta] = (TVectorD*)f->Get(Form("demo/GRe_%d_%d",ibin,itheta));
					GIm_t[ibin][itheta] = (TVectorD*)f->Get(Form("demo/GIm_%d_%d",ibin,itheta));
					for(ir=0; ir<nstepr; ir++){
						GRe[ibin][itheta][ir] += (*GRe_t[ibin][itheta])[ir];
						GIm[ibin][itheta][ir] += (*GIm_t[ibin][itheta])[ir];
						}
				}
			}
			TVectorD* Qx1_t =  (TVectorD*)f->Get(Form("demo/Qx1"));
			TVectorD* Qy1_t =  (TVectorD*)f->Get(Form("demo/Qy1"));
			TVectorD* Q2_t =  (TVectorD*)f->Get(Form("demo/Q2"));
			TVectorD* Nevent_t =  (TVectorD*)f->Get(Form("demo/Nevent"));
			TVectorD* totmultall_t =  (TVectorD*)f->Get(Form("demo/totmultall"));
			for(int ibin=0;ibin<nbin;ibin++){
				Qx1[ibin] += (*Qx1_t)[ibin];
				Qy1[ibin] += (*Qy1_t)[ibin];
				Q2[ibin] += (*Q2_t)[ibin];
				Nevent[ibin] += (*Nevent_t)[ibin];
				totmultall[ibin] += (*totmultall_t)[ibin];
			}
		}
		
	TFile *outf = new TFile("mergedV_Prod.root","Recreate");
	Qx1.Write("Qx1");
	Qy1.Write("Qy1");
	Q2.Write("Q2");
	Nevent.Write("Nevent");
	totmultall.Write("totmultall");

	fstrV<<setprecision(4)<<fixed;
	fstrV<<"ibin"<<"\t"<<"itheta"<<"\t"<<"r0"<<"\t"<<"V"<<"\t"<<"sigma"<<"\t"<<"chi"<<"\t"<<"Vn Errors"<<endl;
		for(int ibin=0;ibin<nbin;ibin++){
			for(int itheta=0;itheta<ntheta;itheta++){
				for(ir=0; ir<nstepr; ir++){
					G[ibin][itheta][ir]=TComplex(GRe[ibin][itheta][ir],GIm[ibin][itheta][ir]);	
					G[ibin][itheta][ir]/=Nevent[ibin];
					Gmod2[ibin][itheta][ir]=TMath::Power(TComplex::Abs(G[ibin][itheta][ir]),2);
				}
			for(ir=0; ir<nstepr-1; ir++)
				if(ir!=0 && Gmod2[ibin][itheta][ir]<=Gmod2[ibin][itheta][ir-1] && Gmod2[ibin][itheta][ir]<=Gmod2[ibin][itheta][ir+1]) break;
			if(ir!=0 && ir<nstepr-1)	r0[ibin][itheta]=r[ibin][ir];
			else if(ir==0)	{cout<<"ibin="<<ibin<<"\titheta="<<itheta<<"\tminimum lies on ir = 0, please select proper range!"<<endl;	continue;}
			else 	{cout<<"ibin="<<ibin<<"\titheta="<<itheta<<"\tminimum lies on ir = maximum "<<nstepr-1<<", please select proper range!"<<endl;	continue;}
			//V[ibin][itheta]=j01/r0[itheta]/avgmult; //simple method
			avgmult[ibin]=1.0*totmultall[ibin]/Nevent[ibin];
			V[ibin][itheta]=Vmax[ibin]-ir*eps[ibin]+eps[ibin]*(Gmod2[ibin][itheta][ir+1]-Gmod2[ibin][itheta][ir-1])/2./(Gmod2[ibin][itheta][ir-1]-2*Gmod2[ibin][itheta][ir]+Gmod2[ibin][itheta][ir+1]);
			V[ibin][itheta]/=avgmult[ibin];
			sigma2[ibin][itheta]=Q2[ibin]/Nevent[ibin]-(Qx1[ibin]/Nevent[ibin])*(Qx1[ibin]/Nevent[ibin])-(Qy1[ibin]/Nevent[ibin])*(Qy1[ibin]/Nevent[ibin])-(V[ibin][itheta]*avgmult[ibin])*(V[ibin][itheta]*avgmult[ibin]);
			Vmean[ibin]+=V[ibin][itheta];
			chi[ibin][itheta]=V[ibin][itheta]*avgmult[ibin]/TMath::Sqrt(sigma2[ibin][itheta]);
			deltaV[ibin][itheta]=V[ibin][itheta]/j01/TMath::BesselJ1(j01)*TMath::Sqrt((TMath::Exp(j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta])+TMath::Exp(-j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta])*TMath::BesselJ0(2*j01))/2./Nevent[ibin]);
			fstrV<<ibin<<"\t"<<itheta<<"\t"<<r0[ibin][itheta]<<"\t"<<V[ibin][itheta]<<"\t"<<sigma2[ibin][itheta]<<"\t"<<chi[ibin][itheta]<<"\t"<<deltaV[ibin][itheta]<<endl;
			deltaVmean[ibin]+=TMath::Exp(j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(nn*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi[ibin][itheta]/chi[ibin][itheta]*TMath::Cos(nn*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(nn*theta[itheta]/2.));
		}
		fstrV<<endl;
		Vmean[ibin]/=ntheta;
		deltaVmean[ibin]=Vmean[ibin]/j01/TMath::BesselJ1(j01)*TMath::Sqrt(deltaVmean[ibin]/ntheta/2./Nevent[ibin]);
	}
	
	fstrV<<"ibin"<<"\t"<<"avgmult"<<"\t"<<"Vn mean"<<"\t"<<"Vn mean Error"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		fstrV<<ibin<<"\t"<<avgmult[ibin]<<"\t"<<Vmean[ibin]<<"\t"<<deltaVmean[ibin]<<endl;
	}
	
        for(int ibin=0;ibin<nbin;ibin++){
        	for(int itheta=0;itheta<ntheta;itheta++){
                        TDirectory *dir = outf->mkdir(Form("D_%d_%d",ibin,itheta));dir->cd();
                        GRe[ibin][itheta].Write(Form("GRe",ibin,itheta));
                        GIm[ibin][itheta].Write(Form("GIm",ibin,itheta));
                        Gmod2[ibin][itheta].Write(Form("G2",ibin,itheta));
                }
        }

	for(int ibin=0;ibin<nbin;ibin++){
		TDirectory *dir = outf->mkdir(Form("D_%d",ibin));dir->cd();
		r[ibin].Write("r");
		r0[ibin].Write("r0");	V[ibin].Write("V");
		sigma2[ibin].Write("sigma2");	chi[ibin].Write("chi");	deltaV[ibin].Write("deltaV");
	}

	TDirectory *dir = outf;	dir->cd();
	
	avgmult.Write("avgmult");	Vmean.Write("Vmean");	deltaVmean.Write("deltaVmean");
}

int getnumber(char* name){
        std::ifstream backstory (name);
        std::string line;
        int nlines=0;
        if (backstory.is_open())
        {
                while (std::getline(backstory, line))
                        nlines++;
                backstory.close();
                return nlines;
        }
        else
        {
                std::cout << "Unable to open file" << std::endl << std::endl;
                return 0;
        }

}

std::string readline(char* name, int iline){
        std::ifstream backstory(name);
        std::string line;
        if (backstory.is_open())
                if(backstory.good()){
                        for(int i = 0; i < iline+1; ++i)
                           getline(backstory, line);
                }
        return line;
}
