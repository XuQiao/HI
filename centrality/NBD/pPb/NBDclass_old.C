#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <TH1.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TRandom.h>
#include <TFile.h>
#include <TString.h>
#include <TList.h>
#include <THStack.h>
#include <TVectorT.h>
#include <TMath.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include "/home/xuq7/HI/centrality/NBD/sort.h"

using namespace std;
//const int nGlau=9;
class NBD{
private:	
	double mumin, mumax, mustep;
	double kmin, kmax, kstep;
	int N;
public:
	TVectorD method;
	TVectorD xmin,xmax;
	TObjString dataname;	//real data root file
	TObjString Glaubername;	//Glauber file root file
	TObjString histoname;	//real data histogram name
	TVectorD chis, Ndf;
	TVectorD mubest, kbest;
	TVectorD NcollAver, NpartAver, BAver;
	TVectorD centbin;	//centrality bins(%)
	TVectorD kpoint;	//bin boundary
	TVectorD centbin_;	//centrality bins using data(%)
	TVectorD kpoint_;	//bin boundary using data(%)
	THStack* Npartdis;
	TGraph2D* Grgrid;

	NBD(){};
	NBD(TString datastr, TString Glauberstr,TString histonamev){
		dataname.SetString(datastr.Data());
		Glaubername.SetString(Glauberstr.Data());
		histoname.SetString(histonamev.Data());
	        mubest.ResizeTo(1);  kbest.ResizeTo(1);
	        chis.ResizeTo(1);  Ndf.ResizeTo(1);
	};
	~NBD(){};
	void initmu(double muminv, double mumaxv, double mustepv){
		mumin=muminv;	mumax=mumaxv;	mustep=mustepv;
	};
	void initk(double kminv, double kmaxv, double kstepv){
		kmin=kminv;	kmax=kmaxv;	kstep=kstepv;
	};
	void initx(double xminv, double xmaxv){
		xmin.ResizeTo(1); xmax.ResizeTo(1);
		xmin[0]=xminv;	xmax[0]=xmaxv;
	};
        void initN(double *binv, int Nv, double methodv){
                centbin.ResizeTo(Nv);   kpoint.ResizeTo(Nv);
                centbin_.ResizeTo(Nv);  kpoint_.ResizeTo(Nv);
		method.ResizeTo(1);
                if(methodv==0){
                        centbin=TVectorD(Nv,binv);
                        centbin_=TVectorD(Nv,binv);
                }
                else{
                        kpoint=TVectorD(Nv,binv);
                        kpoint_=TVectorD(Nv,binv);
                }
                N=Nv-1;
                method[0]=methodv;
                NpartAver.ResizeTo(N); NcollAver.ResizeTo(N);BAver.ResizeTo(N);
                Npartdis = new THStack("Npartdis","");
        };
        void assign(double muv, double kv){
                mubest[0]=muv;  kbest[0]=kv;
        };
        void fit();
        void calcvar();

};

void NBD::fit(){
	TFile *fdata = TFile::Open(dataname.GetName());
	TH1D *histo_obs = (TH1D*)fdata->Get(histoname.GetName());
	TFile *fGlauber = TFile::Open(Glaubername.GetName());
	int binnum = histo_obs->GetNbinsX();
	double Minx = histo_obs->GetXaxis()->GetXmin();
	double Maxx = histo_obs->GetXaxis()->GetXmax();
	double binsize = (Double_t)(Maxx-Minx)/binnum;
	int xbinmin=(int)((xmin[0]-Minx)/binsize);
	int xbinmax=(int)((xmax[0]-Minx)/binsize);
	TH1D *histo_exp = new TH1D("histo_exp","Simulated distribution;Multiplicity;# of events",binnum,Minx,Maxx);
	TF1 *NBD_fun = new TF1("NBD_fun","[0]*TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1]))*TMath::Power([2]/[1],x)/TMath::Power([2]/[1]+1,x+[1])",0,100);
	int ibin;
	TH1D *histo_obs_norm = (TH1D*)histo_obs->Clone();
	histo_obs_norm->Scale(1/histo_obs->Integral(xbinmin,xbinmax));
	double Bin_obs[binnum],Bin_obs_norm[binnum];	
	for (ibin=0; ibin<binnum; ibin++){
		Bin_obs[ibin] = histo_obs->GetBinContent(ibin);
		Bin_obs_norm[ibin]= histo_obs_norm->GetBinContent(ibin);
	}
	int npar=0;
        Double_t chi_square_min[4], chi_square_max[4], kfit[4], mufit[4];
        Double_t chi_add[4], chi_square[4];
	int ncounts=0, ncounts_=0;
	double mu,k;
	std::vector<double> muvector, kvector, chisvector;
	UInt_t iniseed = gRandom->GetSeed();	// reproduce the "random" numbers
	//ofstream fout("chis3.txt");
	for(mu=mumin;mu<=mumax;mu+=mustep){
		for(k=kmin;k<=kmax;k+=kstep){
		//	if(npar%50==0)	cout<<"Have run "<<npar<<" parameter sets"<<endl; 
			NBD_fun->SetParameter(0,1);	//[0]: Normalized constant
			NBD_fun->SetParameter(1,k);	//[1]: k value
			NBD_fun->SetParameter(2,mu);	//[2]: mu value
			TTree *t = (TTree*)fGlauber->Get("nt_p_Pb");

			Float_t Ncoll, Npart, B;	Long_t Nevent;

			t->SetBranchAddress("Ncoll",&Ncoll);
			t->SetBranchAddress("Npart",&Npart);
			t->SetBranchAddress("B",&B);

			Nevent = (Long_t) t->GetEntries();
			Long_t Ev;	Int_t Bino;	Double_t Para, Bi_Para;
			gRandom->SetSeed(iniseed);
		/*	TUnuran unr;
			TUnuranDiscrDist dist(NBD_bin, NBD_bin+Maxfun);
			unr.Init(dist, "method=dgt");	*/

			for (Ev=0; Ev<Nevent; Ev++){
				//if(Ev%100000==0)	cout<<"\t"<<"Have run "<<Ev<<" events"<<endl;
				t->GetEntry(Ev);
				Para = 0; //make sure that Para doesn't accumulate through loops
				for(Bino=0; Bino<Ncoll; Bino++){
				//	Bi_Para = unr.SampleDiscr();
					Bi_Para = NBD_fun->GetRandom();
					Para += Bi_Para;
				}
				histo_exp->Fill(Para);
			}

			Double_t SumEvent, scale;
			SumEvent = histo_exp->Integral(xbinmin,xbinmax);
			if(!SumEvent){
				continue;
			}
			scale = 1/SumEvent;
			TH1D *histo_exp_norm = (TH1D*)histo_exp->Clone();
			histo_exp_norm->Scale(scale);
			Double_t Bin_exp_norm[binnum], Bin_exp[binnum];
		
			for (ibin=0; ibin<binnum; ibin++){
				Bin_exp_norm[ibin] = histo_exp_norm->GetBinContent(ibin);
				Bin_exp[ibin] = histo_exp->GetBinContent(ibin);
			}
			for(int i=0;i<4;i++)	chi_square[i]=0;
			ncounts_=0;
			for(ibin=xbinmin;ibin<xbinmax;ibin++){	
				if(Bin_exp_norm[ibin]==0);
				else {	chi_add[0] = TMath::Power((Bin_obs[ibin]-histo_obs->Integral(xbinmin,xbinmax)*Bin_exp_norm[ibin]),2)/(histo_obs->Integral(xbinmin,xbinmax)*Bin_exp_norm[ibin]);	chi_square[0]+=chi_add[0];}
				if(Bin_obs[ibin]==0);
				else{ chi_add[1] = TMath::Power((Bin_obs[ibin]-histo_obs->Integral(xbinmin,xbinmax)*Bin_exp_norm[ibin]),2)/Bin_obs[ibin];	chi_square[1]+=chi_add[1];}
				if(Bin_obs_norm[ibin]==0);
				else{ ncounts_++; chi_add[2] = TMath::Power((histo_exp->Integral(xbinmin,xbinmax)*Bin_obs_norm[ibin]-Bin_exp[ibin]),2)/(histo_exp->Integral(xbinmin,xbinmax)*Bin_obs_norm[ibin]); chi_square[2]+=chi_add[2];}
				if(Bin_exp[ibin]==0);
				else{ chi_add[3] = TMath::Power((histo_exp->Integral(xbinmin,xbinmax)*Bin_obs_norm[ibin] - Bin_exp[ibin]),2)/Bin_exp[ibin]; 	chi_square[3]+=chi_add[3];}
			}
			cout<<mu<<"\t"<<k<<"\t"<<chi_square[2]<<endl;
				muvector.push_back(mu);
				kvector.push_back(k);
				chisvector.push_back(chi_square[2]);
			for(int i=0;i<4;i++){
				if(npar==0){
					chi_square_min[i]=chi_square[i];
					chi_square_max[i]=chi_square[i];
				}
				if(chi_square_min[i]<chi_square[i]);
				else{
					chi_square_min[i]=chi_square[i];
					mufit[i]=mu;
					kfit[i]=k;
					ncounts=ncounts_;
				}
				if(chi_square_max[i]>chi_square[i]);
				else{
					chi_square_max[i]=chi_square[i];
				}
			}
			npar++;
			histo_exp->Reset("M");
			histo_exp_norm->Reset("M");
		}
	}
	double *amu = &muvector[0];
	double *ak = &kvector[0];
	double *achis = &chisvector[0];
	if(chi_square_max[2]>=chi_square_min[2]){
	cout<<"{"<<mufit[2]<<","<<kfit[2]<<"}"<<endl;
	cout<<chi_square_min[2]<<"\t"<<ncounts-3<<endl;
	mubest[0]=mufit[2];kbest[0]=kfit[2];
	chis[0]=chi_square_min[2];
	Ndf[0]=ncounts-3;
	}
	Grgrid = new TGraph2D("Grgrid","",chisvector.size(),amu,ak,achis);
}

void NBD::calcvar(){
	double mumin_temp = mumin;	double mumax_temp = mumax;	
	double kmin_temp = kmin;	double kmax_temp = kmax;	
	this->initmu(mubest[0],mubest[0],mustep);
	this->initk(kbest[0],kbest[0],kstep);
	this->fit();
	this->initmu(mumin_temp,mumax_temp,mustep);
	this->initk(kmin_temp,kmax_temp,kstep);
	TFile *fdata = TFile::Open(dataname.GetName());
	TH1D *histo_obs = (TH1D*)fdata->Get(histoname.GetName());
	TFile *fGlauber = TFile::Open(Glaubername.GetName());
	int binnum = histo_obs->GetNbinsX();
        double Minx = histo_obs->GetXaxis()->GetXmin();
        double Maxx = histo_obs->GetXaxis()->GetXmax();
	TH1D *histo_exp = new TH1D("histo_exp","Simulated distribution;Multiplicity;# of events",binnum,Minx,Maxx);
	TF1 *NBD_fun = new TF1("NBD_fun","[0]*TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1]))*TMath::Power([2]/[1],x)/TMath::Power([2]/[1]+1,x+[1])",0,100);
     	NBD_fun->SetParameter(0,1);
        NBD_fun->SetParameter(1,kbest[0]);
        NBD_fun->SetParameter(2,mubest[0]);
	TTree *t = (TTree*) fGlauber ->Get("nt_p_Pb");
	Float_t Ncoll, Npart, B;
	Long_t Nevent;

	t->SetBranchAddress("Ncoll",&Ncoll);
	t->SetBranchAddress("Npart",&Npart);
	t->SetBranchAddress("B",&B);

	Nevent = (Long_t) t->GetEntries();
	Long_t Ev;	Int_t Bino;	Double_t Para, Bi_Para, Mult;
	
	TString treestr=Form("tree%.f.root",kbest[0]*mubest[0]*1000);
	TFile *treefile = new TFile(treestr,"Recreate");
	TTree *GlauEvent = new TTree("GlauEvent","Glauber Events");
	GlauEvent -> Branch("Mult",&Para,"Mult/D");
	GlauEvent -> Branch("Ncoll",&Ncoll,"Ncoll/F");
	GlauEvent -> Branch("Npart",&Npart,"Npart/F");
	GlauEvent -> Branch("B",&B,"B/F");

	for(Ev=0; Ev<Nevent; Ev++){
		//if(Ev%100000==0)	 cout<<"Have run "<<Ev<<" events"<<endl;
		t->GetEntry(Ev);
		Para = 0; //make sure that Para doesn't accumulate through loops
		for(Bino=0; Bino<Ncoll; Bino++){
		//	Bi_Para = unr.SampleDiscr();
			Bi_Para = NBD_fun->GetRandom();
			Para += Bi_Para;
		}
		histo_exp->Fill(Para);
		GlauEvent->Fill();
	}
        
	GlauEvent->SetBranchAddress("Ncoll",&Ncoll);
        GlauEvent->SetBranchAddress("Npart",&Npart);
        GlauEvent->SetBranchAddress("B",&B);
        GlauEvent->SetBranchAddress("Mult",&Mult);
	
	Int_t i;
	if(method[0]==0){
/*	double Glau_array[1000000];
	Long_t count_Ev=0;
	
	for(Ev=0;Ev<Nevent;Ev++){ 
		GlauEvent->GetEntry(Ev);
		if(Mult>=Minx&&Mult<=Maxx){
			count_Ev++;
			Glau_array[count_Ev-1]=Mult;
		}
	}
	int low=0;
	int high=count_Ev-1;
	int k0[N+1];
	for(i=0;i<N+1;i++){
		k0[i]=count_Ev*(1-centbin[i]);
		kpoint[i]=findkth(Glau_array,low,high,k0[i]);
		count_Ev=0;
		for(Ev=0;Ev<Nevent;Ev++){
			GlauEvent->GetEntry(Ev);
			if(Mult>=Minx&&Mult<=Maxx){
				count_Ev++;
				Glau_array[count_Ev-1]=Mult;
			}
		}
	}
*/
	
	for(i=0;i<N+1;i++){
	kpoint[i] = findpoint(histo_exp,centbin[i]);
	kpoint_[i] = findpoint(histo_obs,centbin[i]);
	}

	}
        
	TH1D* hNpart[N+1];
	for(i=0; i<N; i++){     //Initialization
	hNpart[i] = new TH1D(Form("Npart_%d-%d",i,i+1),Form("Npart distribution for %dth bin",i),80,0,80);
	}
	hNpart[N] = new TH1D(Form("Npart_0-%d",N),Form("Npart distribution for all bin"),80,0,80);

        std::vector<double> PartEvent(N);std::vector<double> PartEvent_(N+1);

        double TotalEvent=0;

	for(Ev=0; Ev<Nevent; Ev++){
		//if(Ev%100000==0) cout<<"Have run "<<Ev<<" events"<<endl;
            	GlauEvent->GetEntry(Ev);
              	for(i=0; i<N; i++){
	        	if(Mult>kpoint[i]&&Mult<=kpoint[i+1]){
        	        	NcollAver[i] += Ncoll;
                	        NpartAver[i] += Npart;
                        	BAver[i] += B;
			   	PartEvent[i]++;
                     	 	TotalEvent++;
				hNpart[i]->Fill(Npart);
			}
		}
		hNpart[N]->Fill(Npart);
                for(i=0; i<N+1; i++){
			if(method[0]!=0 && Mult>kpoint[i])
				PartEvent_[i]++;
               	}
	}

        if(method[0]!=0){
                for(i=0; i<N+1; i++){
		centbin[i]=(double)PartEvent_[i]/TotalEvent;
		centbin_[i]=(double)histo_obs->Integral(histo_obs->GetXaxis()->FindBin(kpoint[i]),histo_obs->GetXaxis()->GetLast())/histo_obs->Integral();
		}
        }

        for(i=0; i<N; i++){
        	if(PartEvent[i]){
                	NcollAver[i]/=PartEvent[i];
                        NpartAver[i]/=PartEvent[i];
                        BAver[i]/=PartEvent[i];
			hNpart[i]->Scale(1.0/TotalEvent);
			Npartdis->Add(hNpart[i]);
         	}
	}
                hNpart[N]->Scale(1.0/TotalEvent);
		Npartdis->Add(hNpart[N]);
		treefile->Close();
		remove(treestr);
}
