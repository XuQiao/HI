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
const int maxNpart = 40;
class GammaF{
private:	
	double thetamin, thetamax, thetastep;
	double kmin, kmax, kstep;
	int N;
public:
	TVectorD method;
	TVectorD xmin,xmax;
	TVectorD k0, theta0;
	TObjString dataname;	//real data root file
	TObjString Glaubername;	//Glauber file root file
	TObjString histoname;	//real data histogram name
	TObjString MCfitname;
	TVectorD chis, Ndf;
	TVectorD thetabest, kbest;
	TVectorD NcollAver, NpartAver, BAver;
	TVectorD centbin;	//centrality bins(%)
	TVectorD kpoint;	//bin boundary
	TVectorD centbin_;	//centrality bins using data(%)
	TVectorD kpoint_;	//bin boundary using data(%)
	THStack* Npartdis;
	TGraph2D* Grgrid;

	GammaF(){};
	GammaF(TString datastr, TString Glauberstr, TString histonamev, TString MCfitstr){
		dataname.SetString(datastr.Data());
		Glaubername.SetString(Glauberstr.Data());
		histoname.SetString(histonamev.Data());
		MCfitname.SetString(MCfitstr.Data());
		ifstream fin(MCfitname.GetName());
	        double k0err,theta0err;
		k0.ResizeTo(1);theta0.ResizeTo(1);
	     	fin>>k0[0]>>k0err;
        	fin>>theta0[0]>>theta0err;
        	//k0[0]=1.39;theta0[0]=3.41;
	        thetabest.ResizeTo(1);  kbest.ResizeTo(1);
	        chis.ResizeTo(1);  Ndf.ResizeTo(1);
	};
	~GammaF(){};
	void inittheta(double thetaminv, double thetamaxv, double thetastepv){
		thetamin=thetaminv;	thetamax=thetamaxv;	thetastep=thetastepv;
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
        void assign(double thetav, double kv){
                thetabest[0]=thetav;  kbest[0]=kv;
        };
        void fit();
        void calcvar();

};

void GammaF::fit(){
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
        TF1 *gammafun[maxNpart];
//      TF1 *gaus = new TF1("gaus","[0]*TMath::Gaus(x,[1]-[2],[3])",0,100);
	int ibin;
	TH1D *histo_obs_norm = (TH1D*)histo_obs->Clone();
	histo_obs_norm->Scale(1/histo_obs->Integral(xbinmin,xbinmax));
	double Bin_obs[binnum],Bin_obs_norm[binnum];	
	for (ibin=0; ibin<binnum; ibin++){
		Bin_obs[ibin] = histo_obs->GetBinContent(ibin);
		Bin_obs_norm[ibin]= histo_obs_norm->GetBinContent(ibin);
	}
	int npar=0;
        Double_t chi_square_min[4], chi_square_max[4], kfit[4], thetafit[4];
        Double_t chi_add[4], chi_square[4];
	int ncounts=0, ncounts_=0;
	double theta,k;
	std::vector<double> thetavector, kvector, chisvector;
	UInt_t iniseed = gRandom->GetSeed();	// reproduce the "random" numbers
	//ofstream fout("chis3.txt");
	for(theta=thetamin;theta<=thetamax;theta+=thetastep){
		for(k=kmin;k<=kmax;k+=kstep){
		        for(int iNpart=0;iNpart<maxNpart;iNpart++){
				gammafun[iNpart] = new TF1("gammafun","TMath::GammaDist(x,[0],0,[1])",0,100);
   				double k_=k0[0]+k*(iNpart-2);
        			double theta_=theta0[0]+theta*TMath::Log(iNpart-1);
        			gammafun[iNpart]->SetParameter(0,k_);   //[1]: k value
        			gammafun[iNpart]->SetParameter(1,theta_);       //[2]: theta value
        		}

        	//	gammafun->SetParameters(k,theta);
		//	if(npar%50==0)	cout<<"Have run "<<npar<<" parameter sets"<<endl; 
			TTree *t = (TTree*)fGlauber->Get("nt_p_Pb");

			Float_t Ncoll, Npart, B;	Long_t Nevent;

			t->SetBranchAddress("Ncoll",&Ncoll);
			t->SetBranchAddress("Npart",&Npart);
			t->SetBranchAddress("B",&B);

			Nevent = (Long_t) t->GetEntries();
			Long_t Ev;	Int_t Bino;	Double_t Para, Bi_Para;
			gRandom->SetSeed(iniseed);
		/*	TUnuran unr;
			TUnuranDiscrDist dist(GammaF_bin, GammaF_bin+Maxfun);
			unr.Init(dist, "method=dgt");	*/

			for (Ev=0; Ev<Nevent; Ev++){
				//if(Ev%10000==0)	cout<<"\t"<<"Have run "<<Ev<<" events"<<endl;
				t->GetEntry(Ev);
				//Para = 0; //make sure that Para doesn't accuthetalate through loops
				//for(Bino=0; Bino<Npart; Bino++){
				//	Bi_Para = unr.SampleDiscr();
				//	Bi_Para =gammafun->Eval(10);
				Para = gammafun[(int)Npart]->GetRandom();
				//	Para += Bi_Para;
				//}
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
			cout<<theta<<"\t"<<k<<"\t"<<chi_square[2]<<endl;
				thetavector.push_back(theta);
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
					thetafit[i]=theta;
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
	double *atheta = &thetavector[0];
	double *ak = &kvector[0];
	double *achis = &chisvector[0];
	if(chi_square_max[2]>=chi_square_min[2]){
	cout<<"{"<<thetafit[2]<<","<<kfit[2]<<"}"<<endl;
	cout<<chi_square_min[2]<<"\t"<<ncounts-3<<endl;
	thetabest[0]=thetafit[2];kbest[0]=kfit[2];
	chis[0]=chi_square_min[2];
	Ndf[0]=ncounts-3;
	}
	Grgrid = new TGraph2D("Grgrid","",chisvector.size(),atheta,ak,achis);
}

void GammaF::calcvar(){
	double thetamin_temp = thetamin;	double thetamax_temp = thetamax;	
	double kmin_temp = kmin;	double kmax_temp = kmax;	
	this->inittheta(thetabest[0],thetabest[0],thetastep);
	this->initk(kbest[0],kbest[0],kstep);
	this->fit();
	this->inittheta(thetamin_temp,thetamax_temp,thetastep);
	this->initk(kmin_temp,kmax_temp,kstep);
	TFile *fdata = TFile::Open(dataname.GetName());
	TH1D *histo_obs = (TH1D*)fdata->Get(histoname.GetName());
	TFile *fGlauber = TFile::Open(Glaubername.GetName());
	int binnum = histo_obs->GetNbinsX();
        double Minx = histo_obs->GetXaxis()->GetXmin();
        double Maxx = histo_obs->GetXaxis()->GetXmax();
	TH1D *histo_exp = new TH1D("histo_exp","Simulated distribution;Multiplicity;# of events",binnum,Minx,Maxx);
        TF1 *gammafun[maxNpart];
        //TF1 *gammafun = new TF1("gammafun","TMath::GammaDist(x,[0],0,[1])",0,100);
        for(int iNpart=0;iNpart<maxNpart;iNpart++){
	        gammafun[iNpart] = new TF1("gammafun","TMath::GammaDist(x,[0],0,[1])",0,100);
                double k_=k0[0]+kbest[0]*(iNpart-2);
        	double theta_=theta0[0]+thetabest[0]*TMath::Log(iNpart-1);
                gammafun[iNpart]->SetParameter(0,k_);   //[1]: k value
                gammafun[iNpart]->SetParameter(1,theta_);       //[2]: theta value
        }
	TTree *t = (TTree*) fGlauber ->Get("nt_p_Pb");
	Float_t Ncoll, Npart, B;
	Long_t Nevent;

	t->SetBranchAddress("Ncoll",&Ncoll);
	t->SetBranchAddress("Npart",&Npart);
	t->SetBranchAddress("B",&B);

	Nevent = (Long_t) t->GetEntries();
	Long_t Ev;	Int_t Bino;	Double_t Para, Bi_Para, Mult;
	
	TString treestr=Form("tree%.f.root",kbest[0]*thetabest[0]*1000);
	TFile *treefile = new TFile(treestr,"Recreate");
	TTree *GlauEvent = new TTree("GlauEvent","Glauber Events");
	GlauEvent -> Branch("Mult",&Para,"Mult/D");
	GlauEvent -> Branch("Ncoll",&Ncoll,"Ncoll/F");
	GlauEvent -> Branch("Npart",&Npart,"Npart/F");
	GlauEvent -> Branch("B",&B,"B/F");

	for(Ev=0; Ev<Nevent; Ev++){
		//if(Ev%1000==0)	 cout<<"Have run "<<Ev<<" events"<<endl;
		t->GetEntry(Ev);
		//Para = 0; //make sure that Para doesn't accuthetalate through loops
		//for(Bino=0; Bino<Npart; Bino++){
		//	Bi_Para = unr.SampleDiscr();
		Para = gammafun[(int)Npart]->GetRandom();
		//Para += Bi_Para;
		//}
		histo_exp->Fill(Para);
		GlauEvent->Fill();
	}
        GlauEvent->Write();
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
