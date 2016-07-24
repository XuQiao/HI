#include <iomanip>
#include <iostream>
#include <fstream>
#include "par.h"

void getResV(){

//------------------------Deal with the number and get the output stored in txt files and root files
	
        TF1 *ffit[nbin];
        TF1 *ffit0[nbin][ntheta];
	TFile *f[nFileAll];

	double theta[ntheta];
	TVectorD Nevent;	Nevent.ResizeTo(nbin);  Nevent.Zero();
        TVectorD totmultall;	totmultall.ResizeTo(nbin);      totmultall.Zero();
        TVectorD avgmultall;	avgmultall.ResizeTo(nbin);      avgmultall.Zero();
        TVectorD tottrk;	tottrk.ResizeTo(nbin);      tottrk.Zero();
	TVectorD totptall;      totptall.ResizeTo(nbin);    totptall.Zero();
	TVectorD totetaall;      totetaall.ResizeTo(nbin);    totetaall.Zero();
        TVectorD avgtrk;	avgtrk.ResizeTo(nbin);      avgtrk.Zero();
        TVectorD avgpt;         avgpt.ResizeTo(nbin);
        TVectorD avgeta;         avgeta.ResizeTo(nbin);
        TH1D* hq[nbin][ntheta];
        TH1D* hqx[nbin];
        TH1D* hqy[nbin];
        TH1D* hq2[nbin];
        TH1D* hq2nonf[nbin];
	
	string dir = getenv("DIR");
	for(int itheta=0;itheta<ntheta;itheta++)
        	theta[itheta]=itheta*TMath::Pi()/ntheta/nn;

	for(int ibin=0;ibin<nbin;ibin++){
                        for(int itheta=0;itheta<ntheta;itheta++){
                                hq[ibin][itheta] = new TH1D(Form("hq"),Form("hq"),1000,0,10);
                        }
                        hqx[ibin] = new TH1D(Form("hqx"),Form("hqx"),1000,0,10);
                        hqy[ibin] = new TH1D(Form("hqy"),Form("hqy"),1000,0,10);
                        hq2[ibin] = new TH1D(Form("hq2"),Form("hq2"),1000,0,10);
                        hq2nonf[ibin] = new TH1D(Form("hq2nonf"),Form("hq2nonf"),1000,0,10);
	}

        for(int ifile=0; ifile<nFileAll; ifile++){
		f[ifile] = TFile::Open(Form("/scratch/xuq7/flow/pbsjoboutput/qfit/tracknormcpt03to6/%s/AnaV_Sum_%d.root",dir.c_str(),ifile));
		TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
		TVectorD* totmultall_t =  (TVectorD*)f[ifile]->Get(Form("totmultall"));
		TVectorD* tottrk_t =  (TVectorD*)f[ifile]->Get(Form("tottrk"));
		TVectorD* totptall_t =  (TVectorD*)f[ifile]->Get(Form("totptall"));
		TVectorD* totetaall_t =  (TVectorD*)f[ifile]->Get(Form("totetaall"));
		for(int ibin=0;ibin<nbin;ibin++){
                                TH1D* hqx_t = (TH1D*)f[ifile]->Get(Form("hqx_%d",ibin));
                                hqx[ibin]->Add(hqx_t);
                                TH1D* hqy_t = (TH1D*)f[ifile]->Get(Form("hqy_%d",ibin));
                                hqy[ibin]->Add(hqy_t);
                                TH1D* hq2_t = (TH1D*)f[ifile]->Get(Form("hq2_%d",ibin));
                                hq2[ibin]->Add(hq2_t);
                                TH1D* hq2nonf_t = (TH1D*)f[ifile]->Get(Form("hq2nonf_%d",ibin));
                                hq2nonf[ibin]->Add(hq2nonf_t);
				for(int itheta=0;itheta<ntheta;itheta++){
                                    TH1D* hq_t = (TH1D*)f[ifile]->Get(Form("hq_%d_%d",ibin,itheta));
                                    hq[ibin][itheta]->Add(hq_t);
			        }
				totptall[ibin] += (*totptall_t)[ibin];
				totetaall[ibin] += (*totetaall_t)[ibin];
			        Nevent[ibin] += (*Nevent_t)[ibin];
			        totmultall[ibin] += (*totmultall_t)[ibin];	
			        tottrk[ibin] += (*tottrk_t)[ibin];	
		}
		f[ifile]->Close();
	}
	
	for(int ibin=0;ibin<nbin;ibin++){
		avgmultall[ibin]=totmultall[ibin]/Nevent[ibin];
		avgtrk[ibin]=tottrk[ibin]/Nevent[ibin];
			for(int itheta=0;itheta<ntheta;itheta++){
				avgpt[ibin]=1.0*totptall[ibin]/totmultall[ibin];
				avgeta[ibin]=1.0*totetaall[ibin]/totmultall[ibin];
				}
	}

	TString outname = "mergedV_Sum.root";
	TFile *outf = new TFile(outname,"Recreate");
        Nevent.Write("Nevent");
       	totmultall.Write("totmultall");
       	avgmultall.Write("avgmultall");
       	tottrk.Write("tottrk");
       	avgtrk.Write("avgtrk");
        avgpt.Write("avgpt");
        avgeta.Write("avgeta");
        
	for(int ibin=0;ibin<nbin;ibin++){
                TDirectory *dir0 = outf->mkdir(Form("D_%d",ibin));dir0->cd();
                hqx[ibin]->Write(); hqy[ibin]->Write();hq2[ibin]->Write();hq2nonf[ibin]->Write();

        		for(int itheta=0;itheta<ntheta;itheta++){
                        	TDirectory *dir1 = dir0->mkdir(Form("D_%d",itheta));dir1->cd();
                                hq[ibin][itheta]->Write();
                        }
        	}
	outf->Close();
}
