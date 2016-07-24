#include "TString.h"
#include "TVector.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TComplex.h"
#include "TH1.h"
#include "TH2.h"
#include <fstream>
#include <iostream>
#include "par.h"

void DrawQDis(){
    int ifile=0;
    TFile *f = TFile::Open(Form("/cms/store/user/qixu/flow/NewSTEG/pPbDataV205m150/vndata_50k_%d.root",ifile));
    TTree *t = (TTree*)f->Get("tree");
    TFile *fout = new TFile("Qdis.root","Recreate");
    double theta[ntheta];
    Float_t eta[1000];
    Float_t phi[1000];
    Float_t pt[1000];

    TComplex g[nstepr];
    TVectorD GRe;
    TVectorD GIm;
    TVectorD r;
    GRe.ResizeTo(nstepr);
    GIm.ResizeTo(nstepr);
    GRe.Zero();
    GIm.Zero();
    r.ResizeTo(nstepr);
    double Q,Qx,Qy;
    int mult;
    TH1D* hQ = new TH1D("hQ","hQ",100,-20,20);
    TH1D* hQx = new TH1D("hQx","hQx",100,-20,20);
    TH1D* hQy = new TH1D("hQy","hQy",100,-20,20);
    TH1D* hq = new TH1D("hq","hq",100,-5,5);
    TH1D* hqx = new TH1D("hqx","hqx",100,-5,5);
    TH1D* hqy = new TH1D("hqy","hqy",100,-5,5);
    t->SetBranchAddress("etag",eta);
    t->SetBranchAddress("phig",phi);
    t->SetBranchAddress("ptg",pt);
    t->SetBranchAddress("n",&mult);
    for(int ir=0; ir<nstepr; ir++)
        r[ir]=0.00025*20*(ir+1);
    for(int itheta=0;itheta<ntheta;itheta++)
        theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
    Int_t N = t->GetEntries();
    for(int ievt=0; ievt<N; ievt++){
   if(ievt%5000==0)        cout<<"Processing "<<ievt<<" events"<<endl;
                t->GetEntry(ievt);
		Qx=0,Qy=0;
		for(int itheta=0;itheta<ntheta;itheta++)
                for(int ir=0; ir<nstepr; ir++)
                g[ir]=0.;
		int nTracks = mult;	int xbin=-1;
		for(int j=0;j<nbin;j++)
                if(nTracks<trkbin[j]&&nTracks>=trkbin[j+1])
                        xbin=j;
        	if(xbin<0 || xbin==nbin) continue;	
		for(int imult=0;imult<mult;imult++){
			if(eta[imult]<-2.40||eta[imult]>2.40) continue;
                        if(pt[imult]<ptmin||pt[imult]>ptmax) continue; //event selection
			Qx+=1.*cos(nn*phi[imult]);
			Qy+=1.*sin(nn*phi[imult]);
	//	    for(int itheta=0;itheta<ntheta;itheta++){
        //                double temp=TMath::Cos(nn*(phi[imult]-theta[itheta]));
          //          }
                }
		for(int itheta=0;itheta<ntheta;itheta++){
		    Q=Qx*TMath::Cos(nn*theta[itheta])+Qy*TMath::Sin(nn*theta[itheta]);
                }

	    for(int ir=0; ir<nstepr; ir++){
                
		g[ir]=TComplex::Exp(TComplex::I()*r[ir]*Qx);
	        GRe[ir]+=1./N*g[ir].Re();
                GIm[ir]+=1./N*g[ir].Im();
             }
           double q = Q/TMath::Sqrt(mult);
           double qx = Qx/TMath::Sqrt(mult);
           double qy = Qy/TMath::Sqrt(mult);
           hQ->Fill(Q);
           hQx->Fill(Qx);
           hQy->Fill(Qy);
           hq->Fill(q);
           hqx->Fill(qx);
           hqy->Fill(qy);
       }
    fout->cd();
    hQ->Write();
    hQx->Write();
    hQy->Write();
    hq->Write();
    hqx->Write();
    hqy->Write();
    r.Write(Form("r"));
    GRe.Write(Form("GRe"));
    GIm.Write(Form("GIm"));
    fout->Close();
}

