#include "TComplex.h"
#include "TString.h"
#include "TVector.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include <fstream>
#include <iostream>
#include "par.h"

//
// class declaration
//

class qfit {
   public:
      qfit(TString);
      ~qfit();

	void beginJob(int ispt_=1) ;
	void calcV();
	void calcv();
     	void endJobV(TString) ;
     	void endJobv(TString) ;

      // ----------member data ---------------------------
   private:
            
        int nvv, ispt;
        const double *binv;
	TVectorD Nevent, totmultall, tottrk, totptall,totetaall;
	TString filename;
	double theta[ntheta];
	TVectorD totmultv[nbin];
	TVectorD totptv[nbin];
	TVectorD totetav[nbin];
        TVectorD q22;
        TVectorD q24;
        TH1D* hq[nbin][ntheta];
        TH1D* hqx[nbin];
        TH1D* hqy[nbin];
        TH1D* hq2[nbin];
        TH1D* hq2nonf[nbin];

};

//
// constants, enums and typedefs

//double Besselj01(int);



qfit::qfit(TString filenamev)
{	
	filename = filenamev;
}

qfit::~qfit(){
    for(int ibin=0; ibin<nbin; ibin++){
     	for(int itheta=0;itheta<ntheta;itheta++){
                delete   hq[ibin][itheta];
        }
        delete   hqx[ibin];
        delete   hqy[ibin];
        delete   hq2[ibin];
        delete   hq2nonf[ibin];
    }
}


// ------------ method called for each event  ------------
void
qfit::calcV()
{
   using namespace std;
	
        double Q[ntheta];	double Qx;	double Qy;
	int mult, ntrk;
        Float_t phi[10000], eta[10000], pt[10000];
        TFile *infile = TFile::Open(filename);
        TTree* t = (TTree*)infile->Get("tree");
        t->SetBranchAddress("n",&mult);
        t->SetBranchAddress("ptg",pt);
        t->SetBranchAddress("phig",phi);
        t->SetBranchAddress("etag",eta);
        int nevent = t->GetEntries();
        for(int ievt=0; ievt<nevent; ievt++){
                t->GetEntry(ievt);
	//	if(ievt%10000==0) cout<<"has processed "<<ievt<<" events"<<endl;
			for(int itheta=0;itheta<ntheta;itheta++){
       			Q[itheta]=0;
	        	}
			Qx=0,Qy=0;
		ntrk = mult; int xbin=-1;

		for(int j=0;j<nbin;j++)
                if(ntrk<trkbin[j]&&ntrk>=trkbin[j+1])
                        xbin=j;
        	if(xbin<0 || xbin==nbin) continue;	
                tottrk[xbin]+=ntrk;
                int xmult=0;
		for(int imult=0;imult<mult;imult++){
			if(eta[imult]<-2.40||eta[imult]>2.40) continue;
                        if(pt[imult]<ptmin||pt[imult]>ptmax) continue; //event selection
			Qx+=1.*cos(nn*phi[imult]);
			Qy+=1.*sin(nn*phi[imult]);
                        totptall[xbin]+=pt[imult];
                        totetaall[xbin]+=eta[imult];
			totmultall[xbin]++;
                        xmult++;
		}
                if(xmult==0) continue;
			for(int itheta=0;itheta<ntheta;itheta++){
				Q[itheta]=Qx*TMath::Cos(nn*theta[itheta])+Qy*TMath::Sin(nn*theta[itheta]);
                                hq[xbin][itheta]->Fill(TMath::Abs(Q[itheta])/TMath::Sqrt(xmult));
		}
                hqx[xbin]->Fill(TMath::Abs(Qx)/TMath::Sqrt(xmult));
                hqy[xbin]->Fill(TMath::Abs(Qy)/TMath::Sqrt(xmult));
                double q2 = TMath::Sqrt(Qx*Qx+Qy*Qy)/TMath::Sqrt(xmult);
                hq2[xbin]->Fill(q2);
                hq2nonf[xbin]->Fill(TMath::Sqrt(2)*q2);
                q22[xbin]+=q2*q2;
                q24[xbin]+=q2*q2*q2*q2;
		Nevent[xbin]++;
	}
	infile->Close();
}


void qfit::calcv(){
		
        int mult, ntrk;
        Float_t phi[10000], eta[10000], pt[10000];
        TFile *infile = TFile::Open(filename);
        TTree* t = (TTree*)infile->Get("tree");
        t->SetBranchAddress("n",&mult);
        t->SetBranchAddress("ptg",pt);
        t->SetBranchAddress("phig",phi);
        t->SetBranchAddress("etag",eta);
        int nevent = t->GetEntries();

	for(int ievt=0; ievt<nevent; ievt++){
		t->GetEntry(ievt);
	//	if(ievt%100==0) cout<<"has processed "<<ievt<<" events"<<endl;
		
		ntrk = mult; int xbin=-1;
		for(int j=0;j<nbin;j++)
                	if(ntrk<trkbin[j]&&ntrk>=trkbin[j+1])
                	        xbin=j;
        	if(xbin<0 || xbin == nbin) continue;
		tottrk[xbin]+=ntrk;
		for(int imult=0;imult<mult;imult++){
			if(eta[imult]<-2.40||eta[imult]>2.40) continue;
			if(pt[imult]<ptmin||pt[imult]>ptmax) continue; //event selection	
			int xv=-1;
			for(int ivbin=0;ivbin<nvv; ivbin++){
                            if(ispt && pt[imult]>binv[ivbin] && pt[imult]<=binv[ivbin+1])
                                xv = ivbin;
                            if((!ispt) && eta[imult]>binv[ivbin] && eta[imult]<=binv[ivbin+1])
                                xv = ivbin;
                        }
			if(xv<0 || xv==nvv)	continue;
			totptv[xbin][xv]+=pt[imult];
			totetav[xbin][xv]+=eta[imult];
                        totmultv[xbin][xv]++;
		        totmultall[xbin]++;
		}

	Nevent[xbin]++;
	}
	infile->Close();
}


// ------------ method called once each job just before starting event loop  ------------
void 
qfit::beginJob(int ispt_)
{
    ispt = ispt_;
    if(ispt_){             nvv = nptv;       binv = ptbinv;}
    else{             nvv = netav;          binv = etabinv;}

	Nevent.ResizeTo(nbin);	totmultall.ResizeTo(nbin), tottrk.ResizeTo(nbin), totptall.ResizeTo(nbin), totetaall.ResizeTo(nbin);
	Nevent.Zero();	totmultall.Zero(),	tottrk.Zero(); totptall.Zero(); totetaall.Zero();
        q22.ResizeTo(nbin);q24.ResizeTo(nbin);
        q22.Zero();q24.Zero();

        	for(int itheta=0;itheta<ntheta;itheta++){
               		theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
		}
	
	for(int ibin=0; ibin<nbin; ibin++){
                totmultv[ibin].ResizeTo(nvv);	totmultv[ibin].Zero();  
                totptv[ibin].ResizeTo(nvv);	totptv[ibin].Zero();  
                totetav[ibin].ResizeTo(nvv);	totetav[ibin].Zero();  
                for(int itheta=0;itheta<ntheta;itheta++){
                    hq[ibin][itheta] = new TH1D(Form("hq_%d_%d",ibin,itheta),Form("hq_%d_%d",ibin,itheta),1000,0,10);
                    hq[ibin][itheta]->Sumw2();
                }
                hqx[ibin] = new TH1D(Form("hqx_%d",ibin),Form("hqx_%d",ibin),1000,0,10);
                hqy[ibin] = new TH1D(Form("hqy_%d",ibin),Form("hqy_%d",ibin),1000,0,10);
                hq2[ibin] = new TH1D(Form("hq2_%d",ibin),Form("hq2_%d",ibin),1000,0,10);
                hq2nonf[ibin] = new TH1D(Form("hq2nonf_%d",ibin),Form("hq2nonf_%d",ibin),1000,0,10);
                hqx[ibin]->Sumw2();
                hqy[ibin]->Sumw2();
                hq2[ibin]->Sumw2();
                hq2nonf[ibin]->Sumw2();
        }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
qfit::endJobV(TString outstr) 
{
        
	TFile *fs = new TFile(outstr,"Recreate");
	fs->cd();
	Nevent.Write("Nevent");
	totmultall.Write("totmultall");
	totptall.Write("totptall");
	totetaall.Write("totetaall");
	tottrk.Write("tottrk");
        q22.Write("q22");
        q24.Write("q24");
	for(int ibin=0; ibin<nbin; ibin++){
        		for(int itheta=0;itheta<ntheta;itheta++){
                                hq[ibin][itheta]->Write();
		}
                hqx[ibin]->Write();
                hqy[ibin]->Write();
                hq2[ibin]->Write();
                hq2nonf[ibin]->Write();
	}
	fs->Close();
}

void qfit::endJobv(TString outstr){
	
   using namespace std;

	TFile *fs = new TFile(outstr,"Recreate");
	fs->cd();
        Nevent.Write("Nevent");
        totmultall.Write("totmultall");
	tottrk.Write("tottrk");

	for(int ibin=0; ibin<nbin; ibin++){
                for(int itheta=0;itheta<ntheta;itheta++){
                }
		totmultv[ibin].Write(Form("totmult_%d",ibin));
		totptv[ibin].Write(Form("totpt_%d",ibin));
		totetav[ibin].Write(Form("toteta_%d",ibin));
        }
	fs->Close();


}
