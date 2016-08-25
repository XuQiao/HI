#include "TComplex.h"
#include "TString.h"
#include "TVector.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include "par.h"

//
// class declaration
//

class LYZ {
   public:
      LYZ(TString);
      ~LYZ(){};

	void beginJob() ;
	void calcV(int);
	void calcv(TString, int, int);
     	void endJobV(TString) ;
     	void endJobv(TString) ;

      // ----------member data ---------------------------
   private:

	TVectorD Nevent, totmultall, tottrk;
	TString filename;
        TVectorD Qx1[nbin];	
	TVectorD Qy1[nbin];
	TVectorD Q2[nbin];	
	double Qx2[nbin][netaV+1],Qy2[nbin][netaV+1];
	double theta[ntheta];
        double r[nbin][netaV+1][nstepr];
	TVectorD GRe[nbin][netaV+1][ntheta];
	TVectorD GIm[nbin][netaV+1][ntheta];
	TVectorD dDRe[nbin], dNRe[nbin][ntheta];
	TVectorD dDIm[nbin], dNIm[nbin][ntheta];
	TVectorD totmultV[nbin];
	TVectorD totmultv[nbin];
	TVectorD totptV[nbin];
	TVectorD totptv[nbin];

};

//
// constants, enums and typedefs

//double Besselj01(int);



LYZ::LYZ(TString filenamev)
{	
	filename = filenamev;
}

// ------------ method called for each event  ------------
void
LYZ::calcV(int way)	//way=0: Prod way=1: Sum
{
   using namespace std;


	TComplex g[netaV+1][ntheta][nstepr];
	double Q[netaV+1][ntheta];	double Qx[netaV+1];	double Qy[netaV+1];
	int mult, ntrk;
        double phi[10000], eta[10000], pt[10000];
        TFile *infile = TFile::Open(filename);
        TTree* t = (TTree*)infile->Get("demo/TrackTree");
        t->SetBranchAddress("mult",&mult);
        t->SetBranchAddress("ntrk",&ntrk);
        t->SetBranchAddress("pt",pt);
        t->SetBranchAddress("phi",phi);
        t->SetBranchAddress("eta",eta);
        int nevent = t->GetEntries();
        for(int ievt=0; ievt<nevent; ievt++){
                t->GetEntry(ievt);
	//	if(ievt%10000==0) cout<<"has processed "<<ievt<<" events"<<endl;
		for(int ietabin=0; ietabin<=netaV;ietabin++){
			for(int itheta=0;itheta<ntheta;itheta++){
       			Q[ietabin][itheta]=0;
                	for(int ir=0; ir<nstepr; ir++)
                		if(way==0)
                        		g[ietabin][itheta][ir]=1.;
                        	else
                        		g[ietabin][itheta][ir]=0.;
	        	}
			Qx[ietabin]=0,Qy[ietabin]=0;
		}

		int nTracks = ntrk;	int xbin=-1;
		for(int j=0;j<nbin;j++)
                if(nTracks<trkbin[j]&&nTracks>=trkbin[j+1])
                        xbin=j;
        	if(xbin<0 || xbin==nbin) continue;	
                tottrk[xbin]+=ntrk;
		for(int imult=0;imult<mult;imult++){
			if(eta[imult]<-2.40||eta[imult]>2.40) continue;
                        if(pt[imult]<ptmin||pt[imult]>ptmax) continue; //event selection
			int xeta=-1;
			for(int ietabin=0;ietabin<netaV; ietabin++){
                        	if(eta[imult]>etabinV[ietabin] && eta[imult]<=etabinV[ietabin+1])
                        		xeta = ietabin;
                        }
			if(xeta<0 || xeta==netaV)	continue;
                        totptV[xbin][xeta]+=pt[imult];
                        totptV[xbin][netaV]+=pt[imult];
                        totmultV[xbin][xeta]++;
                        totmultV[xbin][netaV]++;
			Qx[xeta]+=1.*cos(nn*phi[imult]);
			Qy[xeta]+=1.*sin(nn*phi[imult]);
                        Qx[netaV]+=1.*cos(nn*phi[imult]);
			Qy[netaV]+=1.*sin(nn*phi[imult]);
			for(int itheta=0;itheta<ntheta;itheta++){
                        	double temp=TMath::Cos(nn*(phi[imult]-theta[itheta]));
				for(int ir=0; ir<nstepr; ir++)
	                           if(way==0){
					g[xeta][itheta][ir]*=TComplex(1.,r[xbin][xeta][ir]*temp);
					g[netaV][itheta][ir]*=TComplex(1.,r[xbin][netaV][ir]*temp);
                                   }
                	        //   else 
				//	g[itheta][0]+=temp;
                	}
			totmultall[xbin]++;
		}
		for(int ietabin=0; ietabin<=netaV;ietabin++){
			Qx1[xbin][ietabin]+=Qx[ietabin];	Qy1[xbin][ietabin]+=Qy[ietabin];
			Qx2[xbin][ietabin]+=Qx[ietabin]*Qx[ietabin];	Qy2[xbin][ietabin]+=Qy[ietabin]*Qy[ietabin];
			for(int itheta=0;itheta<ntheta;itheta++){
				Q[ietabin][itheta]=Qx[ietabin]*TMath::Cos(nn*theta[itheta])+Qy[ietabin]*TMath::Sin(nn*theta[itheta]);
				for(int ir=0; ir<nstepr; ir++){
					if(way!=0)	g[ietabin][itheta][ir]=TComplex::Exp(TComplex(0,r[xbin][ietabin][ir]*Q[ietabin][itheta]));//g[itheta][0]
					GRe[xbin][ietabin][itheta][ir]+=g[ietabin][itheta][ir].Re();
					GIm[xbin][ietabin][itheta][ir]+=g[ietabin][itheta][ir].Im();
				}
			}
		}

		Nevent[xbin]++;
	}
	infile->Close();
}


void LYZ::calcv(TString res, int way, int isample){	//way=0: product way=1: sum
	TComplex g0[ntheta], dDsum[ntheta], dNsum[ntheta][nptv];
	TVectorD* r0res[nbin]; 
	//TVectorD* Vres[nbin]; TVectorD* chires[nbin];
	TFile *fres = TFile::Open(res);
	int xeta=netaV;
	for(int ibin=0;ibin<nbin;ibin++){
		if(isample>=0)
                	r0res[ibin] = (TVectorD*)fres->Get(Form("D_%d/s_%d/D_%d/r0",ibin,isample,xeta));
                else
			r0res[ibin] = (TVectorD*)fres->Get(Form("D_%d/D_%d/r0",ibin,xeta));
	}
		
        int mult, ntrk;
        double phi[10000], eta[10000], pt[10000];
        TFile *infile = TFile::Open(filename);
        TTree* t = (TTree*)infile->Get("demo/TrackTree");
        t->SetBranchAddress("mult",&mult);
        t->SetBranchAddress("ntrk",&ntrk);
        t->SetBranchAddress("pt",pt);
        t->SetBranchAddress("phi",phi);
        t->SetBranchAddress("eta",eta); 
        int nevent = t->GetEntries();
	
	for(int ievt=0; ievt<nevent; ievt++){
		t->GetEntry(ievt);
	//	if(ievt%100==0) cout<<"has processed "<<ievt<<" events"<<endl;
		
		for(int iptbin=0; iptbin<nptv;iptbin++)
                        for(int itheta=0;itheta<ntheta;itheta++)
                        dNsum[itheta][iptbin]=0.;
                for(int itheta=0;itheta<ntheta;itheta++){
                        if(way==0)
				g0[itheta]=1.;
			else
                        	g0[itheta]=0.;
                        dDsum[itheta]=0.;
                }
		int nTracks = ntrk;	int xbin=-1;
		for(int j=0;j<nbin;j++)
                	if(nTracks<trkbin[j]&&nTracks>=trkbin[j+1])
                	        xbin=j;
        	if(xbin<0 || xbin == nbin) continue;
		tottrk[xbin]+=ntrk;
		for(int imult=0;imult<mult;imult++){
			if(eta[imult]<-2.40||eta[imult]>2.40) continue;
			if(pt[imult]<ptmin||pt[imult]>ptmax) continue; //event selection	
			int xpt=-1;
			for(int iptbin=0;iptbin<nptv; iptbin++){
                                if(pt[imult]>ptbinv[iptbin] && pt[imult]<=ptbinv[iptbin+1])
                                xpt = iptbin;
                        }
			if(xpt<0 || xpt==nptv)	continue;
			totptv[xbin][xpt]+=pt[imult];
                        totmultv[xbin][xpt]++;
			for(int itheta=0;itheta<ntheta;itheta++){
                                Double_t temp=TMath::Cos(nn*(phi[imult]-theta[itheta]));
				TComplex temp1(1.,(*r0res[xbin])[itheta]*temp);
                                if(way==0)
					g0[itheta]*=temp1;
                                else 
					g0[itheta]+=temp;
                                dDsum[itheta]+=temp/temp1;
                                dNsum[itheta][xpt]+=TMath::Cos(mm*nn*(phi[imult]-theta[itheta]))/temp1;
			}
		totmultall[xbin]++;
		}

		for(int itheta=0;itheta<ntheta;itheta++){
			dDRe[xbin][itheta]+=(g0[itheta]*dDsum[itheta]).Re();
			dDIm[xbin][itheta]+=(g0[itheta]*dDsum[itheta]).Im();
			for(int iptbin=0;iptbin<nptv; iptbin++){
				dNRe[xbin][itheta][iptbin]+=(g0[itheta]*dNsum[itheta][iptbin]).Re();
				dNIm[xbin][itheta][iptbin]+=(g0[itheta]*dNsum[itheta][iptbin]).Im();
                	}
		}
	Nevent[xbin]++;
	}
	fres->Close();
	infile->Close();
}


// ------------ method called once each job just before starting event loop  ------------
void 
LYZ::beginJob()
{
  double Vmax[nbin][netaV+1], eps[nbin][netaV+1];
  for(int ibin=0; ibin<nbin ;ibin++){
	for(int ietabin=0;ietabin<netaV; ietabin++){
            Vmax[ibin][ietabin]=0.065*(trkbin[ibin]+30)/netaV;
            eps[ibin][ietabin]=0.00025*(trkbin[ibin]+30)/netaV;
        }
        Vmax[ibin][netaV]=0.065*(trkbin[ibin]+30);
        eps[ibin][netaV]=0.00025*(trkbin[ibin]+30);
  }

	Nevent.ResizeTo(nbin);	totmultall.ResizeTo(nbin), tottrk.ResizeTo(nbin);
	Nevent.Zero();	totmultall.Zero(),	tottrk.Zero();
	for(int ibin=0; ibin<nbin; ibin++){
		for(int ietabin=0;ietabin<=netaV; ietabin++){
        		Qx2[ibin][ietabin]=0;  Qy2[ibin][ietabin]=0;
		}
	}

        	for(int itheta=0;itheta<ntheta;itheta++){
               		theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
			for(int ibin=0; ibin<nbin; ibin++){
				for(int ietabin=0;ietabin<=netaV; ietabin++){
					GRe[ibin][ietabin][itheta].ResizeTo(nstepr); 	GRe[ibin][ietabin][itheta].Zero();
					GIm[ibin][ietabin][itheta].ResizeTo(nstepr);	GIm[ibin][ietabin][itheta].Zero();
				}
        		}
		}
	
	for(int ibin=0; ibin<nbin; ibin++){
		for(int ietabin=0;ietabin<=netaV; ietabin++)
        	    for(int ir=0; ir<nstepr; ir++)
                        if(isSimple==0)	r[ibin][ietabin][ir]=j01/(Vmax[ibin][ietabin]-eps[ibin][ietabin]*ir);
                	else	r[ibin][ietabin][ir]=0.00025*20*(ir+1);
		Qx1[ibin].ResizeTo(netaV+1);	Qy1[ibin].ResizeTo(netaV+1);	Q2[ibin].ResizeTo(netaV+1);
		Qx1[ibin].Zero();	Qy1[ibin].Zero();	Q2[ibin].Zero();
                totmultV[ibin].ResizeTo(netaV+1);	totmultV[ibin].Zero();  
                totmultv[ibin].ResizeTo(nptv);	totmultv[ibin].Zero();  
                totptV[ibin].ResizeTo(netaV+1);	totptV[ibin].Zero();  
                totptv[ibin].ResizeTo(nptv);	totptv[ibin].Zero();  
		dDRe[ibin].ResizeTo(ntheta);	dDRe[ibin].Zero();
		dDIm[ibin].ResizeTo(ntheta);	dDIm[ibin].Zero();
                for(int itheta=0;itheta<ntheta;itheta++){
                    dNRe[ibin][itheta].ResizeTo(nptv);	dNRe[ibin][itheta].Zero();
                    dNIm[ibin][itheta].ResizeTo(nptv);	dNIm[ibin][itheta].Zero();
                }
        }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LYZ::endJobV(TString outstr) 
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
		for(int ietabin=0;ietabin<=netaV; ietabin++){
			Q2[ibin][ietabin]+=Qx2[ibin][ietabin]+Qy2[ibin][ietabin];
                	for(int ir=0; ir<nstepr; ir++){
                        	for(int itheta=0;itheta<ntheta;itheta++){
 //                               	fstrGRe<<GRe[ibin][ietabin][itheta][ir]<<"\t";
  //                              	fstrGIm<<GIm[ibin][ietabin][itheta][ir]<<"\t";
				}
  //                      	fstrGRe<<std::endl;
  //                      	fstrGIm<<std::endl;
                	}
 //               	fstrGRe<<std::endl<<std::endl;
 //               	fstrGIm<<std::endl<<std::endl;
 //               	fstrQx1<<Qx1[ibin][ietabin]<<"\t";
 //               	fstrQy1<<Qy1[ibin][ietabin]<<"\t";
 //               	fstrQ2<<Qx2[ibin][ietabin]+Qy2[ibin][ietabin]<<"\t";
		}       
//                fstrGRe<<std::endl<<std::endl;
//               fstrGIm<<std::endl<<std::endl;
 //               fstrQx1<<std::endl;
 //               fstrQy1<<std::endl;
 //               fstrQ2<<std::endl;
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
	tottrk.Write("tottrk");
	for(int ibin=0; ibin<nbin; ibin++){
		for(int ietabin=0;ietabin<=netaV; ietabin++){
        		for(int itheta=0;itheta<ntheta;itheta++){
				GRe[ibin][ietabin][itheta].Write(Form("GRe_%d_%d_%d",ibin,ietabin,itheta));
				GIm[ibin][ietabin][itheta].Write(Form("GIm_%d_%d_%d",ibin,ietabin,itheta));
			}
		}
	Qx1[ibin].Write(Form("Qx1_%d",ibin));
	Qy1[ibin].Write(Form("Qy1_%d",ibin));
	Q2[ibin].Write(Form("Q2_%d",ibin));
	totptV[ibin].Write(Form("totpt_%d",ibin));
	totmultV[ibin].Write(Form("totmult_%d",ibin));
	}
	fs->Close();
}

void LYZ::endJobv(TString outstr){
	
   using namespace std;

        ofstream fstrdNRe, fstrdNIm, fstrdDRe, fstrdDIm, fstrmult;
//        fstrdNRe.open("dNRe.txt",ios::app);
//        fstrdNIm.open("dNIm.txt",ios::app);
 //       fstrdDRe.open("dDRe.txt",ios::app);
//        fstrdDIm.open("dDIm.txt",ios::app);
//        fstrmult.open("mult_.txt",ios::app);


	for(int ibin=0;ibin<nbin;ibin++){
  //      	fstrmult<<Nevent[ibin]<<endl;
        	for(int iptbin=0;iptbin<nptv; iptbin++){
   //             	fstrmult<<(double)totmultv[ibin][iptbin]/Nevent[ibin]<<"\t";
    //            	fstrmult<<(double)totptv[ibin][iptbin]/totmultv[ibin][iptbin]<<"\t";
		}
        	for(int itheta=0;itheta<ntheta;itheta++){
                	for(int iptbin=0;iptbin<nptv; iptbin++){
//                        	fstrdNRe<<dNRe[ibin][itheta][iptbin]<<"\t";
 //                       	fstrdNIm<<dNIm[ibin][itheta][iptbin]<<"\t";
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
        }
	fs->Close();


}
