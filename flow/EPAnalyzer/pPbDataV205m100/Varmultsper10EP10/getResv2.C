#include <iomanip>
#include <fstream>
#include "par.h"

void getResv2(){
	
	ifstream f_("merged.root");
	if(!f_.good()){
	system("hadd merged.root /scratch/xuq7/flow/EPAnalyzer/pPbDataV205m100/Varmultper10EP10/AnaEP_*.root");	
	}
	TFile *fmerged = TFile::Open("merged.root");
	TH1D* hEPPhip[nbin][neta];
	TH1D* hEPPhin[nbin][neta];
	TH1D* hEPPhic[nbin][neta];
	for(int ibin=0;ibin<nbin;ibin++){
		for(int ieta=0;ieta<neta;ieta++){
			hEPPhip[ibin][ieta] = (TH1D*)fmerged->Get(Form("hEPPhip_%d_%d",ibin,ieta));
			hEPPhin[ibin][ieta] = (TH1D*)fmerged->Get(Form("hEPPhin_%d_%d",ibin,ieta));
			hEPPhic[ibin][ieta] = (TH1D*)fmerged->Get(Form("hEPPhic_%d_%d",ibin,ieta));
		}
	}
	ofstream  fstrV;
	fstrV.open("v2.txt");
	TFile *f[nFileAll];

       TVectorD Nevent;	Nevent.ResizeTo(nbin);  Nevent.Zero();
       TVectorD totmultall;	totmultall.ResizeTo(nbin);      totmultall.Zero();
       TVectorD EPR[nbin], EPRpc[nbin], EPRpm[nbin], EPRcm[nbin];	
	for(int ibin=0;ibin<nbin;ibin++){
	EPR[ibin].ResizeTo(neta);	EPR[ibin].Zero();
	EPRpc[ibin].ResizeTo(neta);	EPRpc[ibin].Zero();
	EPRpm[ibin].ResizeTo(neta);	EPRpm[ibin].Zero();
	EPRcm[ibin].ResizeTo(neta);	EPRcm[ibin].Zero();
	}
       TVectorD avgmult;	avgmult.ResizeTo(nbin);
       TVectorD totmult[nbin], totpt[nbin], avgpt[nbin];
       TVectorD Vobs[nbin][neta], Vcorr[nbin][neta], Vcorrerr[nbin][neta], Vcorr3[nbin][neta], Vcorr3err[nbin][neta], Vobserr[nbin][neta];
       TVectorD totmultp[nbin][neta], totmultn[nbin][neta], totmultc[nbin][neta];
       for(int ibin=0;ibin<nbin;ibin++){
       	totmult[ibin].ResizeTo(npt);  totmult[ibin].Zero();
       	totpt[ibin].ResizeTo(npt);  totpt[ibin].Zero();
       	avgpt[ibin].ResizeTo(npt);  avgpt[ibin].Zero();
		for(int ieta=0;ieta<neta;ieta++){
		       	Vobs[ibin][ieta].ResizeTo(npt);  Vobs[ibin][ieta].Zero();
		       	Vobserr[ibin][ieta].ResizeTo(npt);  Vobserr[ibin][ieta].Zero();
		       	totmultp[ibin][ieta].ResizeTo(npt);  totmultp[ibin][ieta].Zero();
		       	totmultn[ibin][ieta].ResizeTo(npt);  totmultn[ibin][ieta].Zero();
		       	totmultc[ibin][ieta].ResizeTo(npt);  totmultc[ibin][ieta].Zero();
		       	Vcorr[ibin][ieta].ResizeTo(npt);  Vcorr[ibin][ieta].Zero();
		       	Vcorrerr[ibin][ieta].ResizeTo(npt);  Vcorrerr[ibin][ieta].Zero();
		       	Vcorr3[ibin][ieta].ResizeTo(npt);  Vcorr3[ibin][ieta].Zero();
		       	Vcorr3err[ibin][ieta].ResizeTo(npt);  Vcorr3err[ibin][ieta].Zero();
		}
	}

	for(int ifile=0; ifile<nFileAll; ifile++){
		f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/EPAnalyzer/pPbDataV205m100/Varmultper10EP10/AnaEP_%d.root",ifile));
		TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
		TVectorD* totmultall_t =  (TVectorD*)f[ifile]->Get(Form("totmultall"));
		for(int ibin=0;ibin<nbin;ibin++){
			TVectorD* EPR_t =  (TVectorD*)f[ifile]->Get(Form("EPR_%d",ibin));
			TVectorD* EPRpc_t =  (TVectorD*)f[ifile]->Get(Form("EPRpc_%d",ibin));
			TVectorD* EPRpm_t =  (TVectorD*)f[ifile]->Get(Form("EPRpm_%d",ibin));
			TVectorD* EPRcm_t =  (TVectorD*)f[ifile]->Get(Form("EPRcm_%d",ibin));
			TVectorD* totmult_t = (TVectorD*)f[ifile]->Get(Form("totmult_%d",ibin));
			TVectorD* totpt_t = (TVectorD*)f[ifile]->Get(Form("totpt_%d",ibin));
			Nevent[ibin] += (*Nevent_t)[ibin];
			totmultall[ibin] += (*totmultall_t)[ibin];
			for(int iptbin=0;iptbin<npt;iptbin++){
				totmult[ibin][iptbin] += (*totmult_t)[iptbin];
				totpt[ibin][iptbin] += (*totpt_t)[iptbin];
			}
			for(int ieta=0;ieta<neta;ieta++){
				EPR[ibin][ieta] += (*EPR_t)[ieta];
				EPRpc[ibin][ieta] += (*EPRpc_t)[ieta];
				EPRpm[ibin][ieta] += (*EPRpm_t)[ieta];
				EPRcm[ibin][ieta] += (*EPRcm_t)[ieta];
				TVectorD* Vobs_t = (TVectorD*)f[ifile]->Get(Form("v2obs_%d_%d",ibin,ieta));
				TVectorD* totmultp_t = (TVectorD*)f[ifile]->Get(Form("totmultp_%d_%d",ibin,ieta));
				TVectorD* totmultn_t = (TVectorD*)f[ifile]->Get(Form("totmultn_%d_%d",ibin,ieta));
				TVectorD* totmultc_t = (TVectorD*)f[ifile]->Get(Form("totmultc_%d_%d",ibin,ieta));
				for(int iptbin=0;iptbin<npt;iptbin++){
					Vobs[ibin][ieta][iptbin] += (*Vobs_t)[iptbin];
					totmultp[ibin][ieta][iptbin] += (*totmultp_t)[iptbin];
					totmultn[ibin][ieta][iptbin] += (*totmultn_t)[iptbin];
					totmultc[ibin][ieta][iptbin] += (*totmultc_t)[iptbin];
				}
			}
		}	
		f[ifile]->Close();
	}
		
	TFile *outf = new TFile("mergedVobs.root","Recreate");
	Nevent.Write("Nevent");
	totmultall.Write("totmultall");
	for(int ibin=0;ibin<nbin;ibin++){
		for(int ieta=0;ieta<neta;ieta++){
			hEPPhip[ibin][ieta]->Write();
			hEPPhin[ibin][ieta]->Write();
			hEPPhic[ibin][ieta]->Write();
		}
	}

	fstrV<<setprecision(4)<<fixed;
	fstrV<<"ibin"<<"\t"<<"ieta"<<"\t"<<"iptbin"<<"\t"<<"avgpt"<<"\t"<<"v2obs"<<"\t"<<"v2obserr"<<"v2corr"<<"\t"<<"v2correrr"<<"\t"<<"v2corr3"<<"\t"<<"v2corr3err"<<"\t"<<"R"<<"\t"<<"R3"<<endl;
	for(int ibin=0;ibin<nbin;ibin++){
		avgmult[ibin]=(double)totmultall[ibin]/Nevent[ibin];
		for(int iptbin=0;iptbin<npt;iptbin++){
			avgpt[ibin][iptbin]=totpt[ibin][iptbin]/totmult[ibin][iptbin];
		}
		for(int ieta=0;ieta<neta;ieta++){
			//EPR[ibin][ieta]=EPR[ibin][ieta]/Nevent[ibin];
			EPR[ibin][ieta]=sqrt(EPR[ibin][ieta]/Nevent[ibin]);
			EPRpc[ibin][ieta]=EPRpc[ibin][ieta]/Nevent[ibin];
			EPRcm[ibin][ieta]=EPRcm[ibin][ieta]/Nevent[ibin];
			EPRpm[ibin][ieta]=EPRpm[ibin][ieta]/Nevent[ibin];
			for(int iptbin=0;iptbin<npt;iptbin++){
				Vobs[ibin][ieta][iptbin]/=totmultc[ibin][ieta][iptbin];
				Vobserr[ibin][ieta][iptbin]=sqrt(1.0/totmultc[ibin][ieta][iptbin]);
				Vcorr[ibin][ieta][iptbin]=Vobs[ibin][ieta][iptbin]/EPR[ibin][ieta];
				Vcorrerr[ibin][ieta][iptbin]=Vobserr[ibin][ieta][iptbin]/EPR[ibin][ieta];
				Vcorr3[ibin][ieta][iptbin]=Vobs[ibin][ieta][iptbin]/sqrt(EPRpc[ibin][ieta]*EPRpm[ibin][ieta]/EPRcm[ibin][ieta]);
				Vcorr3err[ibin][ieta][iptbin]=Vobserr[ibin][ieta][iptbin]/sqrt(EPRpc[ibin][ieta]*EPRpm[ibin][ieta]/EPRcm[ibin][ieta]);
			fstrV<<trkbin[ibin+1]<<"-"<<trkbin[ibin]<<"\t"<<etap[ieta]<<"-"<<etan[ieta]<<"\t"<<ptbin[iptbin]<<"-"<<ptbin[iptbin+1]<<"\t"<<avgpt[ibin][iptbin]<<"\t"<<Vobs[ibin][ieta][iptbin]<<"\t"<<Vobserr[ibin][ieta][iptbin]<<"\t"<<Vcorr[ibin][ieta][iptbin]<<"\t"<<Vcorrerr[ibin][ieta][iptbin]<<"\t"<<Vcorr3[ibin][ieta][iptbin]<<"\t"<<Vcorr3err[ibin][ieta][iptbin]<<"\t"<<EPR[ibin][ieta]<<"\t"<<sqrt(EPRpc[ibin][ieta]*EPRpm[ibin][ieta]/EPRcm[ibin][ieta])<<endl;
			}
			fstrV<<endl;
		}
		TDirectory *dir0 = outf->mkdir(Form("D_%d",ibin));	dir0->cd();
		totmult[ibin].Write("totmult");
		avgpt[ibin].Write("avgpt");
		EPR[ibin].Write("EPR");
		EPRpc[ibin].Write("EPRpc");
		EPRpm[ibin].Write("EPRpm");
		EPRcm[ibin].Write("EPRcm");
		for(int ieta=0;ieta<neta;ieta++){
			TDirectory *dir1 = dir0->mkdir(Form("E_%d",ieta));	dir1->cd();
			Vobs[ibin][ieta].Write("v2obs");
			Vobserr[ibin][ieta].Write("v2obserr");
			totmultp[ibin][ieta].Write("totmultp");
			totmultn[ibin][ieta].Write("totmultn");
			totmultc[ibin][ieta].Write("totmultc");
			Vcorr[ibin][ieta].Write("v2");
			Vcorrerr[ibin][ieta].Write("v2err");
			Vcorr3[ibin][ieta].Write("v23sub");
			Vcorr3err[ibin][ieta].Write("v23suberr");
		}
	}
	outf->cd();
	avgmult.Write("avgmult");
	outf->Close();
	
	remove("merged.root");
}
