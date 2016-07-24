#include <fstream>
#include <iostream>
#include "TComplex.h"
#include "TNamed.h"
#include "TObject.h"
#include "TString.h"
#include "TKey.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TIterator.h"
#include "TDirectory.h"

//Usage:

const int npt=10;
const int ntheta=5;
const int nstepr=500;

double Besselj01(int);

class LYZ: public TNamed{
	private:
		TString instr;
		ios_base::openmode mode;
		int p0, m, n;
		int mult; Float_t phi[1000], eta[1000], pt[1000];
		Long_t nevent;
		double Qx1,Qx2,Qy1,Qy2,theta[ntheta];
		double r[nstepr];
		TComplex G[ntheta][nstepr];
		double r0[ntheta], V[ntheta], chi[ntheta];
		TComplex dD[ntheta];	TComplex dN[ntheta][npt];
		int Nevent,totmultall;
		int totmult[npt];
	public:
		LYZ(){ p0=5; n=2; m=1; instr = "/home/xuq7/CMSSW_6_2_3_patch1/src/flow/proSTEGvn/vndataSTEG_50k.root";};
		LYZ(TString strname){ p0=5; n=2; m=1; instr = strname;};
		LYZ(int p0value, int nvalue, int mvalue, TString strname){ p0=p0value; n=nvalue; m=mvalue; instr = strname;};
		~LYZ(){};
		void init();
		void init(int, int);
		void calcV(int way=0);
		void calcv(int way=0);
		void getResV();
		void getResv();
		void Add(const LYZ *l,const int isV);
		void OwnCopy(const LYZ *l);
    
	ClassDef(LYZ, 1); //The class title
};

ClassImp(LYZ);

const double Vmax=0.15;
const double eps=0.0002;
const double j01=2.404826;
const double ptbin[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,10.0};
//const double ptbin[]={0.2,0.3,0.5,0.8,1.0,1.2,1.5,1.8,2.0,2.2,2.5,2.8,3.0,3.3,3.8,4.0,4.5,5.0,5.4,6.0,7.0,8.0,9.0,10.0};

void LYZ::init(){
        Nevent=0;
	totmultall=0;
	Qx1=0;	Qy1=0;	Qx2=0;	Qy2=0;
        for(int itheta=0;itheta<ntheta;itheta++){
                theta[itheta]=itheta*TMath::Pi()/p0/n;
        	for(int ir=0; ir<nstepr; ir++)
                        G[itheta][ir]=0;
        }

        for(int ir=0; ir<nstepr; ir++)
		r[ir]=j01/((Vmax-eps*ir)*nstepr);
	//r[ir]=eps*(ir+1);

        for(int ipt=0; ipt<npt;ipt++){
                totmult[ipt]=0;
        	for(int itheta=0;itheta<ntheta;itheta++){
                dN[itheta][ipt]=0.;
                dD[itheta]=0.;
		}
	}
}

void LYZ::init(int number, int isapp){
	if(number>0)
		nevent = number;
	else nevent = 0;	
	if(isapp==0)	
		mode = ios::trunc;
	else	
		mode = ios::app;
	this->init();
}

void LYZ::calcV(int way){	//way=0: product way=1: sum
	TComplex g[ntheta][nstepr];
	double Q[ntheta];
	ofstream fstrG, fstrQx1, fstrQy1, fstrQ2, fstrmult;
        fstrG.open("G.txt",mode);
	fstrQx1.open("Qx1.txt",mode);
	fstrQy1.open("Qy1.txt",mode);
	fstrQ2.open("Q2.txt",mode);
	TFile *infile = TFile::Open(instr);
	TTree* t = (TTree*)infile->Get("tree");
        t->SetBranchAddress("n",&mult);
        t->SetBranchAddress("ptg",pt);
        t->SetBranchAddress("phig",phi);
        t->SetBranchAddress("etag",eta);
	if(nevent==0) nevent = t->GetEntries();
	for(int ievt=0; ievt<nevent; ievt++){
		t->GetEntry(ievt);
		if(ievt%100==0) cout<<"has processed "<<ievt<<" events"<<endl;
		
			double Qx=0,Qy=0;

		for(int itheta=0;itheta<ntheta;itheta++){
			Q[itheta]=0;
			for(int ir=0; ir<nstepr; ir++)
                		if(way==0)
					g[itheta][ir]=1.;
				else
					g[itheta][ir]=0.;
                	}

		for(int ipart=0; ipart<mult; ipart++){
			if(eta[ipart]<-2.40||eta[ipart]>2.40) continue;
			if(pt[ipart]<0.2||pt[ipart]>10.0) continue; //event selection
				Qx+=1.*cos(n*phi[ipart]);
				Qy+=1.*sin(n*phi[ipart]);
			/*for(int itheta=0;itheta<ntheta;itheta++){
                                double temp=TMath::Cos(n*(phi[ipart]-theta[itheta]));
				for(int ir=0; ir<nstepr; ir++)
	                              if(way==0)
						g[itheta][ir]*=TComplex(1.,r[ir]*temp);
                	                else 
						g[itheta][0]+=temp;
                        }*/
			totmultall++;
		}

			Qx1+=Qx;	Qy1+=Qy;
			Qx2+=Qx*Qx;	Qy2+=Qy*Qy;
		for(int itheta=0;itheta<ntheta;itheta++){
			Q[itheta]=Qx*TMath::Cos(n*theta[itheta])+Qy*TMath::Sin(n*theta[itheta]);
			for(int ir=0; ir<nstepr; ir++){
				if(way!=0)	g[itheta][ir]=TComplex::Exp(TComplex(0,r[ir]*Q[itheta]));//g[itheta][0]
				G[itheta][ir]+=g[itheta][ir];
			}
		}

	Nevent++;
	}
	for(int ir=0; ir<nstepr; ir++){
		for(int itheta=0;itheta<ntheta;itheta++)
        		fstrG<<G[itheta][ir]<<"\t";
              	fstrG<<endl;
	}
		fstrQx1<<Qx1<<endl;
		fstrQy1<<Qy1<<endl;
		fstrQ2<<Qx2+Qy2<<endl;
}

void LYZ::calcv(int way){	//way=0: product way=1: sum
	TComplex g0[ntheta], dDsum[ntheta], dNsum[ntheta][npt];
	int ipt=0;
 	ofstream fstrdN, fstrdD, fstrmult;
        fstrdN.open("dN.txt",mode);
        fstrdD.open("dD.txt",mode);
        fstrmult.open("mult.txt",mode);
	TFile *infile = TFile::Open(instr);
	TTree* t = (TTree*)infile->Get("tree");
        t->SetBranchAddress("n",&mult);
        t->SetBranchAddress("ptg",pt);
        t->SetBranchAddress("phig",phi);
        t->SetBranchAddress("etag",eta);
	if(nevent==0) nevent = t->GetEntries();
	for(int ievt=0; ievt<nevent; ievt++){
		t->GetEntry(ievt);
		if(ievt%100==0) cout<<"has processed "<<ievt<<" events"<<endl;
		
		for(int iptbin=0; iptbin<npt;iptbin++)
                        for(int itheta=0;itheta<ntheta;itheta++)
                        dNsum[itheta][iptbin]=0.;
                for(int itheta=0;itheta<ntheta;itheta++){
                        if(way==0)
				g0[itheta]=1.;
			else
                        	g0[itheta]=0.;
                        dDsum[itheta]=0.;
                }
		for(int ipart=0; ipart<mult; ipart++){
			if(eta[ipart]<-2.40||eta[ipart]>2.40) continue;
			if(pt[ipart]<0.2||pt[ipart]>10.0) continue; //event selection	
			for(int iptbin=0;iptbin<npt; iptbin++){
                                if(pt[ipart]>ptbin[iptbin] && pt[ipart]<=ptbin[iptbin+1])
                                ipt = iptbin;
                        }
                        totmult[ipt]++;
			for(int itheta=0;itheta<ntheta;itheta++){
                                Double_t temp=TMath::Cos(n*(phi[ipart]-theta[itheta]));
				TComplex temp1(1.,r0[itheta]*temp);
                                if(way==0)
					g0[itheta]*=temp1;
                                else 
					g0[itheta]+=temp;
                                dDsum[itheta]+=temp/temp1;
                                dNsum[itheta][ipt]+=TMath::Cos(m*n*(phi[ipart]-theta[itheta]))/temp1;
			}
		totmultall++;
		}

		for(int itheta=0;itheta<ntheta;itheta++){
			dD[itheta]+=g0[itheta]*dDsum[itheta];
			for(int iptbin=0;iptbin<npt; iptbin++){
				dN[itheta][iptbin]+=g0[itheta]*dNsum[itheta][iptbin];
                	}
		}
	Nevent++;
	}

        fstrmult<<Nevent<<endl;
        for(int iptbin=0;iptbin<npt; iptbin++)
        	fstrmult<<totmult[iptbin]<<"\t";
	for(int itheta=0;itheta<ntheta;itheta++){
		for(int iptbin=0;iptbin<npt; iptbin++)
                      	fstrdN<<dN[itheta][iptbin]<<"\t";
                fstrdN<<endl;
                fstrdD<<dD[itheta]<<"\t";
	}
	fstrmult<<endl<<totmultall<<endl;
	fstrmult<<endl;
        fstrdN<<endl;
	fstrdD<<endl;
}

void LYZ::getResV(){
	ofstream fstrG2, fstrV;
	int ir;
	double avgmult=1.0*totmultall/Nevent;
	double Gmod2[ntheta][nstepr];
	double sigma2[ntheta],deltaV[ntheta];
	double deltaVmean=0, Vmean=0;
	fstrV.open("V.txt");
        fstrG2.open("G2.txt",mode);

        for(int itheta=0;itheta<ntheta;itheta++){
                for(ir=0; ir<nstepr; ir++){
                        G[itheta][ir]/=Nevent;
                        Gmod2[itheta][ir]=TMath::Power(TComplex::Abs(G[itheta][ir]),2);
			fstrG2<<r[ir]<<"\t"<<Gmod2[itheta][ir]<<"\t";
		}
		fstrG2<<endl;
		for(ir=0; ir<nstepr-1; ir++)
			if(Gmod2[itheta][ir]<=Gmod2[itheta][ir+1]) break;
		if(ir!=0)	r0[itheta]=r[ir];
		//V[itheta]=j01/r0[itheta]/avgmult; //simple method
		V[itheta]=Vmax-ir*eps+eps*(Gmod2[itheta][ir+1]-Gmod2[itheta][ir-1])/2./(Gmod2[itheta][ir-1]-2*Gmod2[itheta][ir]+Gmod2[itheta][ir+1]);
		fstrV<<"r0 "<<itheta<<endl;
		fstrV<<r0[itheta]<<endl;
		fstrV<<"Vn "<<itheta<<endl;
		fstrV<<V[itheta]<<endl;
		sigma2[itheta]=(Qx2+Qy2)/Nevent-(Qx1/Nevent)*(Qx1/Nevent)-(Qy1/Nevent)*(Qy1/Nevent)-(V[itheta]*avgmult)*(V[itheta]*avgmult);
		Vmean+=V[itheta];
		chi[itheta]=V[itheta]*avgmult/TMath::Sqrt(sigma2[itheta]);
		fstrV<<chi[itheta]<<endl;
		deltaV[itheta]=V[itheta]/j01/TMath::BesselJ1(j01)*TMath::Sqrt((TMath::Exp(j01*j01/2./chi[itheta]/chi[itheta])+TMath::Exp(-j01*j01/2./chi[itheta]/chi[itheta])*TMath::BesselJ0(2*j01))/2./Nevent);
		fstrV<<"Vn Errors "<<itheta<<endl;
		fstrV<<deltaV[itheta]<<endl<<endl;
		deltaVmean+=TMath::Exp(j01*j01/2./chi[itheta]/chi[itheta]*TMath::Cos(n*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(n*theta[itheta]/2.))+TMath::Exp(-j01*j01/2./chi[itheta]/chi[itheta]*TMath::Cos(n*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(n*theta[itheta]/2.));
	}
	Vmean/=ntheta;
	deltaVmean=Vmean/j01/TMath::BesselJ1(j01)*TMath::Sqrt(deltaVmean/ntheta/2./Nevent);
	fstrV<<"Vn mean "<<endl;
	fstrV<<Vmean<<endl;
	fstrV<<"Vn Errors mean "<<endl;
	fstrV<<deltaVmean<<endl;
}

void LYZ::getResv(){
	ofstream fstrv;
	fstrv.open("v.txt");
	double avgmult=1.0*totmultall/Nevent;
	double v[ntheta][npt],vmean[npt],deltav[ntheta][npt],deltavmean[npt];
	for(int itheta=0;itheta<ntheta;itheta++){
		dD[itheta]*=TComplex::Power(TComplex::I(),m-1)/(Double_t)Nevent;
	}
	for(int iptbin=0;iptbin<npt; iptbin++){
		vmean[iptbin]=0;
		deltavmean[iptbin]=0;
		for(int itheta=0;itheta<ntheta;itheta++){
			dN[itheta][iptbin]/=totmult[iptbin];
			TComplex Res=dN[itheta][iptbin]/dD[itheta];
			v[itheta][iptbin]=V[itheta]*avgmult*TMath::BesselJ1(j01)/Besselj01(m)*Res.Re();
			vmean[iptbin]+=v[itheta][iptbin];
			deltav[itheta][iptbin]=TMath::Sqrt(TMath::Cos(m*n*theta[itheta])/totmult[iptbin]*(TMath::Exp(j01*j01/2./chi[itheta]/chi[itheta]*TMath::Cos(n*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Sin(n*theta[itheta]/2.))+TMath::Power(-1,m)*TMath::Exp(-j01*j01/2./chi[itheta]/chi[itheta]*TMath::Cos(n*theta[itheta]))*TMath::BesselJ0(2*j01*TMath::Cos(n*theta[itheta]/2.))))/2./Besselj01(m);
			deltavmean[iptbin]+=deltav[itheta][iptbin];
		}
		vmean[iptbin]/=ntheta;
		deltavmean[iptbin]/=ntheta;
		fstrv<<vmean[iptbin]<<"\t"<<deltavmean[iptbin]<<endl;
	}
}

void LYZ::Add(const LYZ *l,const int isV){
	if(isV==1){
		Nevent+=l->Nevent;
		totmultall+=l->totmultall;
		Qx1+=l->Qx1;	Qy1+=l->Qy1;	Qx2+=l->Qx2;	Qy2+=l->Qy2;
		for(int itheta=0;itheta<ntheta;itheta++){
                	for(int ir=0; ir<nstepr; ir++)
				G[itheta][ir]+=l->G[itheta][ir];
		}
	}
	else{
		Nevent+=l->Nevent;
		totmultall+=l->totmultall;
		for(int itheta=0;itheta<ntheta;itheta++){
			dD[itheta]+=l->dD[itheta];
		}
		for(int iptbin=0;iptbin<npt; iptbin++){
			totmult[iptbin]+=l->totmult[iptbin];
			for(int itheta=0;itheta<ntheta;itheta++){
				dN[itheta][iptbin]+=l->dN[itheta][iptbin];
			}
		}
	}
}

void LYZ::OwnCopy(const LYZ *l){
	for(int itheta=0;itheta<ntheta;itheta++){
		r0[itheta]=l->r0[itheta];
		V[itheta]=l->V[itheta];
		chi[itheta]=l->chi[itheta];
	}
}

double Besselj01(int m0){
	if(m0==1)
		return TMath::BesselJ1(j01);
	else if(m0==2)
		return 0.431755;
	else if(m0==3)
		return 0.199000;
	else if(m0==4)
		return 0.064747;
	else return 0;
}

void storeV(){
       	TFile *f = new TFile("LYZV.root","Update");
	TString fstr(getenv("FILENAME"));
	LYZ *def;
	if(fstr.IsNull()){
                def=new LYZ();
                def->SetName(Form("default"));
        }
        else{
                def=new LYZ(fstr);
	}
	cout<<"Integrated Vn"<<endl;
       	def->init(0,0);
       	def->calcV(1);
	f->cd();
       	def->Write();
	f->Close();
}

void getResultV(){
        TFile *f = TFile::Open("LYZV.root","ReadOnly");
        TFile *fo = TFile::Open("LYZres.root","Recreate");
        LYZ *sum = new LYZ();
	sum->init();
        TIter nextkey(f->GetListOfKeys());
   	TKey *key;
   	while ((key = (TKey*)nextkey())) {
      		LYZ *def = (LYZ*)key->ReadObj();
        	sum->Add(def,1);
        }
	sum->getResV();
	fo->cd();
	sum->Write("sumV");
	fo->Close();
	f->Close();
}

void storev(){
	TFile *fo = new TFile("LYZres.root","ReadOnly");
	TFile *f = new TFile("LYZv.root","Recreate");
	TString fstr = (TString)getenv("FILENAME");
	LYZ *def;
	fo->cd();
	LYZ *sum = (LYZ*)fo->Get("sumV");
	if(fstr.IsNull()){
                def=new LYZ();
        	def->SetName(Form("default"));
	}
        else
                def=new LYZ(fstr);
	def->init(0,0);
        def->OwnCopy(sum);
        def->calcv(0);
	f->cd();
	def->Write();
	fo->Close();
	f->Close();
	}

void getResultv(){
        TFile *f = TFile::Open("LYZv.root","ReadOnly");
        TFile *fo = TFile::Open("LYZres.root","Update");
	LYZ *sum = new LYZ();
        sum->init();
        TIter nextkey(f->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)nextkey())) {
                LYZ *def = (LYZ*)key->ReadObj();
                sum->Add(def,0);
		sum->OwnCopy(def);
        }
        sum->getResv();
        fo->cd();
        sum->Write("sumv");
	fo->Close();
        f->Close();
}
