#include "parameter.h"
#include "par.h"
void hPYphocalc(){
	gStyle->SetOptStat(kFALSE);
	const int maxNpart = 100;
	int sth=0, Gth=0;
        TFile *f = TFile::Open(outG);
        if(sth==0){TString dirname = "std";}
        else if(sth==1){TString dirname ="Gri055";}
        else {TString dirname ="Gri101";}
        TObjString* dataname = (TObjString*)f->Get(Form("dataname"));
        TObjString* histoname = (TObjString*)f->Get(Form("histoname"));
	TFile *fdata = TFile::Open(dataname->GetName());
                TString name;
        if(Gth==0)
                name = "G0";
        else if(Gth<nGlau)
                name = Form("Glau_%d",Gth);
        else
                name = Form("bin_%d",Gth-nGlau+1);
        TObjString* Glaubername = (TObjString*)f->Get(Form("%s/%s/Glaubername",dirname.Data(),name.Data()));
        TVectorD* k0 = (TVectorD*)f->Get(Form("%s/%s/k0",dirname.Data(),name.Data()));
        TVectorD* theta0 = (TVectorD*)f->Get(Form("%s/%s/theta0",dirname.Data(),name.Data()));
        TVectorD* xmin = (TVectorD*)f->Get(Form("%s/%s/xmin",dirname.Data(),name.Data()));
        TVectorD* xmax = (TVectorD*)f->Get(Form("%s/%s/xmax",dirname.Data(),name.Data()));
        TVectorD* thetabest = (TVectorD*)f->Get(Form("%s/%s/thetabest",dirname.Data(),name.Data()));
        TVectorD* kbest = (TVectorD*)f->Get(Form("%s/%s/kbest",dirname.Data(),name.Data()));
	TVectorD* kpoint = (TVectorD*)f->Get(Form("%s/%s/kpoint",dirname.Data(),name.Data()));
        TVectorD* NcollAver = (TVectorD*)f->Get(Form("%s/%s/NcollAver",dirname.Data(),name.Data()));
	TVectorD* centbin = (TVectorD*)f->Get(Form("%s/%s/centbin",dirname.Data(),name.Data()));

        TFile *fGlauber = TFile::Open(Glaubername->GetName());
        //(*k0)[0]=1.39;  (*kbest)[0]=0.425;
        //(*theta0)[0]=3.41;      (*thetabest)[0]=1.30;
        TF1 *gammafun[maxNpart];
        TF1 *gammafunevt[maxNpart];
        TF1 *gammafunnucl[maxNpart];
        TF1 *gammafunnuclNcoll[maxNpart];
        double kevt = (*k0)[0]-(*kbest)[0];
        for(int iNpart=0;iNpart<maxNpart;iNpart++){
        	gammafun[iNpart] = new TF1("gammafun","TMath::GammaDist(x,[0],0,[1])",0,200);
		gammafunevt[iNpart] = new TF1("gammafunevt","TMath::GammaDist(x,[0],0,[1])",0,200);
	        gammafunnucl[iNpart] = new TF1("gammafunnucl","TMath::GammaDist(x,[0],0,[1])",0,200);
        	gammafunnuclNcoll[iNpart] = new TF1("gammafunnuclNcoll","TMath::GammaDist(x,[0],0,[1])",0,200);
                double k_=(*k0)[0]+(*kbest)[0]*(iNpart-2);
                double theta_=(*theta0)[0]+(*thetabest)[0]*TMath::Log(iNpart-1);
                gammafun[iNpart]->SetParameter(0,k_);   //[1]: k value
                gammafun[iNpart]->SetParameter(1,theta_);       //[2]: theta value
		gammafunevt[iNpart]->SetParameter(0,kevt);
		gammafunevt[iNpart]->SetParameter(1,theta_);
		gammafunnucl[iNpart]->SetParameter(0,(*kbest)[0]);
		gammafunnucl[iNpart]->SetParameter(1,theta_);
		gammafunnuclNcoll[iNpart]->SetParameter(0,(*kbest)[0]*(iNpart-1));
		gammafunnuclNcoll[iNpart]->SetParameter(1,theta_);
		if(iNpart==2){
			gammafunnuclNcoll[iNpart]->SetNpx(1e4);
			gammafunnuclNcoll[iNpart]->SetRange(1e-11,200);
		}
       }

        TTree *t = (TTree*)fGlauber->Get("nt_p_Pb");
        Float_t Ncoll, Npart, B;        Long_t Nevent;
        t->SetBranchAddress("Ncoll",&Ncoll);
        t->SetBranchAddress("Npart",&Npart);
        t->SetBranchAddress("B",&B);

        Nevent = (Long_t) t->GetEntries();
        Long_t Ev;      Int_t Bino;     Double_t Para, Para_nucl, Para_p, Para_evt, Bi_Para_nucl, Bi_Para_evt;			

	double yUCM[8]={};
	double yPCM[8]={};
	double yVCM[8]={};
	double yUCM_[200]={};
	double yPCM_[200]={};
	double yVCM_[200]={};
	double C=1e-4;
	double PNcoll[maxNpart]={};
	TH1D *histo_obs = (TH1D*)fdata->Get(histoname->GetName());
        TH1D *histo_obs_norm = (TH1D*)histo_obs->Clone();
        histo_obs_norm->Scale(1/histo_obs->Integral());
	TH1D* hUCM = new TH1D("hUCM","hUCM",200,0,200);
	TH1D* hPCM = new TH1D("hPCM","hPCM",200,0,200);
	TH1D* hVCM = new TH1D("hVCM","hVCM",200,0,200);
	TH2D* NcollvsET = new TH2D("NcollvsET","NcollvsET",100,0,100,2000,0,400);
	for(Ev=0;Ev<Nevent;Ev++){
		t->GetEntry(Ev);
		PNcoll[(int)Ncoll]++;
	}
	for(int i=0;i<maxNpart;i++){
		PNcoll[i]/=Nevent;
		cout<<PNcoll[i]<<"\t";
	}
		cout<<endl;
	for(Ev=0;Ev<Nevent;Ev++){
		if(Ev%100000==0)	cout<<"\t"<<"Have run "<<Ev<<" events"<<endl;
		t->GetEntry(Ev);
		Para = gammafun[(int)Npart]->GetRandom();
		Para_nucl = gammafunnuclNcoll[(int)Npart]->GetRandom();
		Para_p = gammafunnuclNcoll[(int)Npart]->GetRandom();
		Para_evt = 0;
		for(int i=0;i<N-1;i++)
			if(Para>=(*kpoint)[i] && Para<(*kpoint)[i+1])
				int ibin = i;
		for(int Bino=0;Bino<Ncoll;Bino++){
			Bi_Para_evt = gammafunevt[(int)Npart]->GetRandom();
			Para_evt += Bi_Para_evt;
		}	
		double PNcollET = gammafun[(int)Npart]->Eval(Para);
//		double k = gammafun[(int)Npart]->GetParameter(0);
                double theta=(*theta0)[0]+(*thetabest)[0]*TMath::Log(Npart-1);
		double YNcollUCM = C*Ncoll;
		double YNcollPCM = C/1.0/(*kbest)[0]/theta*(Para_nucl);
		double YNcollVCM = C/2.0*(Para_nucl/(*kbest)[0]/theta+Ncoll);
		yUCM[ibin] += PNcoll[(int)Ncoll]*PNcollET*YNcollUCM;
		yPCM[ibin] += PNcoll[(int)Ncoll]*PNcollET*YNcollPCM;
		yVCM[ibin] += PNcoll[(int)Ncoll]*PNcollET*YNcollVCM;
		yUCM_[(int)Para] += PNcoll[(int)Ncoll]*PNcollET*YNcollUCM;
		yPCM_[(int)Para] += PNcoll[(int)Ncoll]*PNcollET*YNcollPCM;
		yVCM_[(int)Para] += PNcoll[(int)Ncoll]*PNcollET*YNcollVCM;

		NcollvsET->Fill(Ncoll,Para);
	}
	for(int ibin=1;ibin<hUCM->GetNbinsX();ibin++){
		hUCM->SetBinContent(ibin,yUCM_[ibin-1]);			
		hPCM->SetBinContent(ibin,yPCM_[ibin-1]);			
		hVCM->SetBinContent(ibin,yVCM_[ibin-1]);
	}
	TCanvas *c1 = new TCanvas();
	TCanvas *c2 = new TCanvas();
	c1->SetLogy();
	c2->SetLogx();
	c2->SetLogy();
	c2->SetLogz();
	c1->cd();
	TH1D* hFrame = new TH1D("","",200,0,200);
	hFrame->SetTitle("");
	hFrame->GetXaxis()->SetTitle("HF #Sigma E_{T} |#eta|>4");
	hFrame->GetYaxis()->SetTitle("Yield no units");
	hFrame->GetXaxis()->SetRangeUser(0,150);
	hFrame->GetYaxis()->SetRangeUser(1e-6,1);
	hFrame->Draw();
	histo_obs_norm->SetMarkerStyle(20);
	histo_obs_norm->SetMarkerSize(1.0);
	histo_obs_norm->SetMarkerColor(1);
	histo_obs_norm->Draw("Psame");
	hUCM->SetMarkerStyle(24);
	hUCM->SetMarkerSize(1.0);
	hUCM->SetMarkerColor(2);
        hPCM->SetMarkerStyle(29);
        hPCM->SetMarkerSize(1.0);
        hPCM->SetMarkerColor(4);
        hVCM->SetMarkerStyle(34);
        hVCM->SetMarkerSize(1.0);
        hVCM->SetMarkerColor(5);
	hUCM->Draw("Psame");
	hPCM->Draw("Psame");
	hVCM->Draw("Psame");
	TLegend *leg = new TLegend(0.1,0.2,0.5,0.45);	
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.03);
	leg->AddEntry(histo_obs_norm,"minimum bias events","lp");
	leg->AddEntry(hUCM,"hard scattering events(UCM)","lp");
	leg->AddEntry(hPCM,"hard scattering events(PCM)","lp");
	leg->AddEntry(hVCM,"hard scattering events(VCM)","lp");
	leg->Draw("same");
	c1->Print("paperfig3_CMS.png");
	c2->cd();
	gStyle->SetOptStat("nemr");
	NcollvsET->GetXaxis()->SetTitle("Ncoll");
	NcollvsET->GetYaxis()->SetTitle("HF #Sigma E_{T} |#eta|>4");
	NcollvsET->Draw("colz");
	c2->Print("NcollvsET2D.png");
	ofstream fstr("result_CMS.dat");
	fstr<<"i"<<"\t"<<"centbin"<<"\t"<<"kpoint"<<"\t"<<"NcollAver"<<"\t"<<"UCM"<<"\t"<<"PCM"<<"\t"<<"VCM"<<"\t"<<"pho1"<<"\t"<<"pho2"<<"\t"<<"MB"<<endl;
	for(int i=0;i<N-1;i++){
		fstr<<i<<"\t"<<(*centbin)[i]*100<<"% to "<<(*centbin)[i+1]*100<<"% \t"<<(*kpoint)[i]<<" to "<<(*kpoint)[i+1]<<"\t"<<(*NcollAver)[i]<<"\t"<<yUCM[i]<<"\t"<<yPCM[i]<<"\t"<<yVCM[i]<<"\t"<<yPCM[i]/yUCM[i]<<"\t"<<yVCM[i]/yUCM[i]<<"\t"<<"undetermined"<<endl;
	}
}
