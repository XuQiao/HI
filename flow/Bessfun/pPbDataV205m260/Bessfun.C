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

class Bessfun {
   public:
      Bessfun(TString);
      ~Bessfun();

	void beginJob(int ispt_=1) ;
	void calcV(int);
	void calcv(TString, int, int);
     	void endJobV(TString) ;
     	void endJobv(TString) ;

      // ----------member data ---------------------------
   private:
            
        int nvv, ispt;
        const double *binv;
	TVectorD Nevent, totmultall, tottrk, totptall,totetaall, Qx1, Qy1, Q2;
	TString filename;
	double Qx2[nbin],Qy2[nbin];
	double theta[ntheta];
        double r[nbin][nstepr];
	TVectorD GRe[nbin][ntheta];
	TVectorD GIm[nbin][ntheta];
	TVectorD dDRe[nbin], dNRe[nbin][ntheta];
	TVectorD dDIm[nbin], dNIm[nbin][ntheta];
	TVectorD totmultv[nbin];
	TVectorD totptv[nbin];
	TVectorD totetav[nbin];
        TH1D* hq[nbin][ntheta];
        TH1D* hqx[nbin];
        TH1D* hqy[nbin];
        TH1D* hq2[nbin];
        TH1D* hq2nonf[nbin];

};

//
// constants, enums and typedefs

//double Besselj01(int);



Bessfun::Bessfun(TString filenamev)
{	
	filename = filenamev;
}

Bessfun::~Bessfun(){
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
Bessfun::calcV(int way)	//way=0: Prod way=1: Sum
{
   using namespace std;


	TComplex g[ntheta][nstepr];
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
                	for(int ir=0; ir<nstepr; ir++)
                		if(way==0)
                        		g[itheta][ir]=1.;
                        	else
                        		g[itheta][ir]=0.;
	        	}
			Qx=0,Qy=0;

		ntrk = mult; int nTracks = ntrk;	int xbin=-1;
		for(int j=0;j<nbin;j++)
                if(nTracks<trkbin[j]&&nTracks>=trkbin[j+1])
                        xbin=j;
        	if(xbin<0 || xbin==nbin) continue;	
                tottrk[xbin]+=ntrk;
		for(int imult=0;imult<mult;imult++){
		//	if(eta[imult]<-2.40||eta[imult]>2.40) continue;
                //       if(pt[imult]<ptmin||pt[imult]>ptmax) continue; //event selection
			Qx+=1.*cos(nn*phi[imult]);
			Qy+=1.*sin(nn*phi[imult]);
			for(int itheta=0;itheta<ntheta;itheta++){
                        	double temp=TMath::Cos(nn*(phi[imult]-theta[itheta]));
				for(int ir=0; ir<nstepr; ir++)
	                           if(way==0){
					g[itheta][ir]*=TComplex(1.,r[xbin][ir]*temp);
                                   }
                	        //   else 
				//	g[itheta][0]+=temp;
                	}
                        totptall[xbin]+=pt[imult];
                        totetaall[xbin]+=eta[imult];
			totmultall[xbin]++;
		}
			Qx1[xbin]+=Qx;	Qy1[xbin]+=Qy;
			Qx2[xbin]+=Qx*Qx;	Qy2[xbin]+=Qy*Qy;
			for(int itheta=0;itheta<ntheta;itheta++){
				Q[itheta]=Qx*TMath::Cos(nn*theta[itheta])+Qy*TMath::Sin(nn*theta[itheta]);
				for(int ir=0; ir<nstepr; ir++){
					if(way!=0)	g[itheta][ir]=TComplex::Exp(TComplex(0,r[xbin][ir]*Qx));//g[itheta][0]
					GRe[xbin][itheta][ir]+=g[itheta][ir].Re();
					GIm[xbin][itheta][ir]+=g[itheta][ir].Im();
				}
                    hq[xbin][itheta]->Fill(TMath::Abs(Q[itheta])/TMath::Sqrt(mult));
		}
                hqx[xbin]->Fill(TMath::Abs(Qx)/TMath::Sqrt(mult));
                hqy[xbin]->Fill(TMath::Abs(Qy)/TMath::Sqrt(mult));
                hq2[xbin]->Fill(TMath::Sqrt(Qx*Qx+Qy*Qy)/TMath::Sqrt(mult));
                hq2nonf[xbin]->Fill(TMath::Sqrt(2)*TMath::Sqrt(Qx*Qx+Qy*Qy)/TMath::Sqrt(mult));
		Nevent[xbin]++;
	}
	infile->Close();
}


void Bessfun::calcv(TString res, int way, int isample){	//way=0: product way=1: sum
	TComplex g0[ntheta], dDsum[ntheta], dNsum[ntheta][nvv];
	TVectorD* r0res[nbin]; 
	//TVectorD* Vres[nbin]; TVectorD* chires[nbin];
	TFile *fres = TFile::Open(res);
	for(int ibin=0;ibin<nbin;ibin++){
		if(isample>=0)
                	r0res[ibin] = (TVectorD*)fres->Get(Form("s_%d/D_%d/r0",isample,ibin));
                else if(isample==-1)
			r0res[ibin] = (TVectorD*)fres->Get(Form("D_%d/r0",ibin));
                else
		        r0res[ibin] = (TVectorD*)fres->Get(Form("D_%d/r02",ibin));
	}
		
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
		
		for(int ivbin=0; ivbin<nvv;ivbin++)
                        for(int itheta=0;itheta<ntheta;itheta++)
                        dNsum[itheta][ivbin]=0.;
                for(int itheta=0;itheta<ntheta;itheta++){
                        if(way==0)
				g0[itheta]=1.;
			else
                        	g0[itheta]=0.;
                        dDsum[itheta]=0.;
                }
		ntrk = mult; int nTracks = ntrk;	int xbin=-1;
		for(int j=0;j<nbin;j++)
                	if(nTracks<trkbin[j]&&nTracks>=trkbin[j+1])
                	        xbin=j;
        	if(xbin<0 || xbin == nbin) continue;
		tottrk[xbin]+=ntrk;
		for(int imult=0;imult<mult;imult++){
	//		if(eta[imult]<-2.40||eta[imult]>2.40) continue;
	//		if(pt[imult]<ptmin||pt[imult]>ptmax) continue; //event selection	
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
			for(int itheta=0;itheta<ntheta;itheta++){
                                Double_t temp=TMath::Cos(nn*(phi[imult]-theta[itheta]));
				TComplex temp1(1.,(*r0res[xbin])[itheta]*temp);
                                if(way==0)
					g0[itheta]*=temp1;
                                else 
					g0[itheta]+=temp;
                                dDsum[itheta]+=temp/temp1;
                                dNsum[itheta][xv]+=TMath::Cos(mm*nn*(phi[imult]-theta[itheta]))/temp1;
			}
		totmultall[xbin]++;
		}

		for(int itheta=0;itheta<ntheta;itheta++){
			dDRe[xbin][itheta]+=(g0[itheta]*dDsum[itheta]).Re();
			dDIm[xbin][itheta]+=(g0[itheta]*dDsum[itheta]).Im();
			for(int ivbin=0;ivbin<nvv; ivbin++){
				dNRe[xbin][itheta][ivbin]+=(g0[itheta]*dNsum[itheta][ivbin]).Re();
				dNIm[xbin][itheta][ivbin]+=(g0[itheta]*dNsum[itheta][ivbin]).Im();
                	}
		}
	Nevent[xbin]++;
	}
	fres->Close();
	infile->Close();
}


// ------------ method called once each job just before starting event loop  ------------
void 
Bessfun::beginJob(int ispt_)
{
  double Vmax[nbin], eps[nbin];
  for(int ibin=0; ibin<nbin ;ibin++){
        Vmax[ibin]=0.065*(trkbin[ibin]+30);
        eps[ibin]=0.00025*(trkbin[ibin]+30);
  }
    ispt = ispt_;
    if(ispt_){             nvv = nptv;       binv = ptbinv;}
    else{             nvv = netav;          binv = etabinv;}

	Nevent.ResizeTo(nbin);	totmultall.ResizeTo(nbin), tottrk.ResizeTo(nbin), totptall.ResizeTo(nbin), totetaall.ResizeTo(nbin);
		Qx1.ResizeTo(nbin);	Qy1.ResizeTo(nbin);	Q2.ResizeTo(nbin);
	Nevent.Zero();	totmultall.Zero(),	tottrk.Zero(); totptall.Zero(); totetaall.Zero();
		Qx1.Zero();	Qy1.Zero();	Q2.Zero();
	for(int ibin=0; ibin<nbin; ibin++){
        		Qx2[ibin]=0;  Qy2[ibin]=0;
	}

        	for(int itheta=0;itheta<ntheta;itheta++){
               		theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
			for(int ibin=0; ibin<nbin; ibin++){
					GRe[ibin][itheta].ResizeTo(nstepr); 	GRe[ibin][itheta].Zero();
					GIm[ibin][itheta].ResizeTo(nstepr);	GIm[ibin][itheta].Zero();
				}
		}
	
	for(int ibin=0; ibin<nbin; ibin++){
        	    for(int ir=0; ir<nstepr; ir++)
                        if(isSimple==0)	r[ibin][ir]=j01/(Vmax[ibin]-eps[ibin]*ir);
                	else	r[ibin][ir]=0.00025*20*(ir+1);
                totmultv[ibin].ResizeTo(nvv);	totmultv[ibin].Zero();  
                totptv[ibin].ResizeTo(nvv);	totptv[ibin].Zero();  
                totetav[ibin].ResizeTo(nvv);	totetav[ibin].Zero();  
		dDRe[ibin].ResizeTo(ntheta);	dDRe[ibin].Zero();
		dDIm[ibin].ResizeTo(ntheta);	dDIm[ibin].Zero();
                for(int itheta=0;itheta<ntheta;itheta++){
                    dNRe[ibin][itheta].ResizeTo(nvv);	dNRe[ibin][itheta].Zero();
                    dNIm[ibin][itheta].ResizeTo(nvv);	dNIm[ibin][itheta].Zero();
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
Bessfun::endJobV(TString outstr) 
{
        
//fstrGRe, fstrGIm
//	ofstream fstrQx1, fstrQy1, fstrQ2, fstrmult;

//        fstrGRe.open("GRe.txt",std::ios::app);
//        fstrGIm.open("GIm.txt",std::ios::app);
//        fstrQx1.open("Qx1.txt",std::ios::app);
//        fstrQy1.open("Qy1.txt",std::ios::app);
 //       fstrQ2.open("Q2.txt",std::ios::app);
 //       fstrmult.open("mult.txt",std::ios::app);

	for(int ibin=0; ibin<nbin; ibin++){
//			fstrmult<<Nevent[ibin]<<"\t";
//        		fstrmult<<totmultall[ibin]<<"\t";
//			fstrmult<<tottrk[ibin]<<std::endl;
			Q2[ibin]+=Qx2[ibin]+Qy2[ibin];
                	for(int ir=0; ir<nstepr; ir++){
                        	for(int itheta=0;itheta<ntheta;itheta++){
 //                               	fstrGRe<<GRe[ibin][itheta][ir]<<"\t";
  //                              	fstrGIm<<GIm[ibin][itheta][ir]<<"\t";
				}
  //                      	fstrGRe<<std::endl;
  //                      	fstrGIm<<std::endl;
                	}
 //               	fstrGRe<<std::endl<<std::endl;
 //               	fstrGIm<<std::endl<<std::endl;
 //               	fstrQx1<<Qx1[ibin]<<"\t";
 //               	fstrQy1<<Qy1[ibin]<<"\t";
 //               	fstrQ2<<Qx2[ibin]+Qy2[ibin]<<"\t";
        }
//	fstrmult.close();
//        fstrGRe.close();
//        fstrGIm.close();
 //       fstrQx1.close();
//        fstrQy1.close();
//        fstrQ2.close();
	TFile *fs = new TFile(outstr,"Recreate");
	fs->cd();
	Nevent.Write("Nevent");
	totmultall.Write("totmultall");
	totptall.Write("totptall");
	totetaall.Write("totetaall");
	tottrk.Write("tottrk");
	Qx1.Write("Qx1");
	Qy1.Write("Qy1");
	Q2.Write("Q2");
	for(int ibin=0; ibin<nbin; ibin++){
        		for(int itheta=0;itheta<ntheta;itheta++){
				GRe[ibin][itheta].Write(Form("GRe_%d_%d",ibin,itheta));
				GIm[ibin][itheta].Write(Form("GIm_%d_%d",ibin,itheta));
                                hq[ibin][itheta]->Write();
		}
                hqx[ibin]->Write();
                hqy[ibin]->Write();
                hq2[ibin]->Write();
                hq2nonf[ibin]->Write();
	}
	fs->Close();
}

void Bessfun::endJobv(TString outstr){
	
   using namespace std;

        ofstream fstrdNRe, fstrdNIm, fstrdDRe, fstrdDIm, fstrmult;
//        fstrdNRe.open("dNRe.txt",ios::app);
//        fstrdNIm.open("dNIm.txt",ios::app);
 //       fstrdDRe.open("dDRe.txt",ios::app);
//        fstrdDIm.open("dDIm.txt",ios::app);
//        fstrmult.open("mult_.txt",ios::app);


	for(int ibin=0;ibin<nbin;ibin++){
  //      	fstrmult<<Nevent[ibin]<<endl;
        	for(int ivbin=0;ivbin<nvv; ivbin++){
   //             	fstrmult<<(double)totmultv[ibin][ivbin]/Nevent[ibin]<<"\t";
    //            	fstrmult<<(double)totptv[ibin][ivbin]/totmultv[ibin][ivbin]<<"\t";
		}
        	for(int itheta=0;itheta<ntheta;itheta++){
                	for(int ivbin=0;ivbin<nvv; ivbin++){
//                        	fstrdNRe<<dNRe[ibin][itheta][ivbin]<<"\t";
 //                       	fstrdNIm<<dNIm[ibin][itheta][ivbin]<<"\t";
			}
  //              	fstrdNRe<<endl;
   //             	fstrdNIm<<endl;
//                	fstrdDRe<<dDRe[ibin][itheta]<<"\t";
 //               	fstrdDIm<<dDIm[ibin][itheta]<<"\t";
        	}
  //      	fstrmult<<tottrk[ibin]<<endl;
  //      	fstrmult<<endl<<totmultall[ibin]<<endl;
	}	
   //     fstrmult<<endl;
    //    fstrdNRe<<endl;
 //       fstrdNIm<<endl;
 //       fstrdDRe<<endl;
 //       fstrdDIm<<endl;
	
	TFile *fs = new TFile(outstr,"Recreate");
	fs->cd();
        Nevent.Write("Nevent");
        totmultall.Write("totmultall");
	tottrk.Write("tottrk");

	for(int ibin=0; ibin<nbin; ibin++){
		dDRe[ibin].Write(Form("dDRe_%d",ibin));
		dDIm[ibin].Write(Form("dDIm_%d",ibin));
                for(int itheta=0;itheta<ntheta;itheta++){
                        dNRe[ibin][itheta].Write(Form("dNRe_%d_%d",ibin,itheta));
                        dNIm[ibin][itheta].Write(Form("dNIm_%d_%d",ibin,itheta));
                }
		totmultv[ibin].Write(Form("totmult_%d",ibin));
		totptv[ibin].Write(Form("totpt_%d",ibin));
		totetav[ibin].Write(Form("toteta_%d",ibin));
        }
	fs->Close();


}
