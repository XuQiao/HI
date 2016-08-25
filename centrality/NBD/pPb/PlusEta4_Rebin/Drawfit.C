#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include "par.h"
#include <vector>
void Drawfit(){
	int sth=0, Gth=0;
	TFile *f = TFile::Open(outG);
	if(sth==0){TString dirname = "std";}
        else if(sth==1){TString dirname ="Gri055";}
        else {TString dirname ="Gri101";}
        TObjString* dataname = (TObjString*)f->Get(Form("dataname"));
        TObjString* histoname = (TObjString*)f->Get(Form("histoname"));
	        TString name;
        if(Gth==0)
                name = "G0";
        else if(Gth<nGlau)
                name = Form("Glau_%d",Gth);
        else
                name = Form("bin_%d",Gth-nGlau+1);
        TObjString* Glaubername = (TObjString*)f->Get(Form("%s/%s/Glaubername",dirname.Data(),name.Data()));
        TVectorD* xmin = (TVectorD*)f->Get(Form("%s/%s/xmin",dirname.Data(),name.Data()));
        TVectorD* xmax = (TVectorD*)f->Get(Form("%s/%s/xmax",dirname.Data(),name.Data()));
        TVectorD* mubest = (TVectorD*)f->Get(Form("%s/%s/mubest",dirname.Data(),name.Data()));
        TVectorD* kbest = (TVectorD*)f->Get(Form("%s/%s/kbest",dirname.Data(),name.Data()));
        TVectorD* Ndf = (TVectorD*)f->Get(Form("%s/%s/Ndf",dirname.Data(),name.Data()));
        TVectorD* chis = (TVectorD*)f->Get(Form("%s/%s/chis",dirname.Data(),name.Data()));
	TVectorD *kpoint = (TVectorD*)f->Get(Form("%s/%s/kpoint",dirname.Data(),name.Data()));
	TFile *fdata = TFile::Open(dataname->GetString());
	TH1D *histo_obs = (TH1D*)fdata->Get(histoname->GetString());
	TFile *fGlauber = TFile::Open(Glaubername->GetString());
        int binnum = histo_obs->GetNbinsX();
        double Minx = histo_obs->GetXaxis()->GetXmin();
        double Maxx = histo_obs->GetXaxis()->GetXmax();
	double binsize = (Double_t)(Maxx-Minx)/binnum;
        int xbinmin=(int)(((*xmin)[0]-Minx)/binsize);
        int xbinmax=(int)(((*xmax)[0]-Minx)/binsize);
	TH1D *histo_exp = new TH1D("histo_exp","Simulated distribution;Multiplicity;Event Fraction",binnum,Minx,Maxx);
	//histo_obs->Sumw2();
	//histo_exp->Sumw2();
	Int_t ibin;
		
	TH1D *histo_obs_norm = (TH1D*)histo_obs->Clone();
        histo_obs_norm->Scale(1/histo_obs->Integral(xbinmin,xbinmax));
	TF1 *NBD_fun = new 
	TF1("NBD_fun","[0]*TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1]))*TMath::Power([2]/[1],x)/TMath::Power([2]/[1]+1,x+[1])",0,100);
	NBD_fun->SetParameter(0,1);	//[0]: Normalized constant
	NBD_fun->SetParameter(1,(*kbest)[0]);	//[1]: k value
	NBD_fun->SetParameter(2,(*mubest)[0]);	//[2]: mu value
		
	TTree *t = (TTree*) fGlauber->Get("nt_p_Pb");
	Long_t Nevent;

	Nevent = (Long_t) t->GetEntries();

	Long_t Ev;	Int_t Bino;	Double_t Para, Bi_Para, Mult;
	Float_t Ncoll;
	t->SetBranchAddress("Ncoll",&Ncoll);

	for(Ev=0; Ev<Nevent; Ev++){
		if(Ev%100000==0)	 cout<<"Have run "<<Ev<<" events"<<endl;
		t->GetEntry(Ev);
		Para = 0; //make sure that Para doesn't accumulate through loops
		for(Bino=0; Bino<Ncoll; Bino++){
                         Bi_Para = NBD_fun->GetRandom();
                         Para += Bi_Para;
		}	
		histo_exp->Fill(Para);
	}
	Double_t SumEvent, scale;
	SumEvent = histo_exp->Integral(xbinmin,xbinmax);
	scale = 1/SumEvent;
	TH1D *histo_exp_norm = (TH1D*) histo_exp->Clone();
	histo_exp_norm->Scale(scale);

        gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas();
	c1->SetLogy();

        histo_obs_norm->GetXaxis()->SetRangeUser(0,100);
        histo_obs_norm->SetMaximum(1.0);
        histo_obs_norm->GetXaxis()->SetTitle("HF #Sigma E_{T} #eta>4");
        histo_obs_norm->GetYaxis()->SetTitle("Event Fraction");
        histo_obs_norm->SetTitle("real data vs. Best fit for #mu and k");

	histo_obs_norm->SetLineColor(1);
	histo_obs_norm->SetMarkerStyle(24);
	histo_obs_norm->SetMarkerColor(1);
        histo_obs_norm->SetMarkerSize(1.5);
	histo_obs_norm->Draw("P");

        TLegend *leg = new TLegend(0.56, 0.7, 0.8, 0.9);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0.035);
        leg->SetTextFont(42);
        leg->SetTextSize(0.045);
        leg->AddEntry(histo_obs_norm,"CMS pPb real data","p");
        leg->AddEntry(histo_exp_norm,Form("NBD Fit, #mu=%.3f, k=%.3f",(*mubest)[0],(*kbest)[0]),"l");
	leg->Draw("same");
        TLatex *tex1= new TLatex(0.2,0.8,Form("#chi^{2}/NDF=%.1f/%.f",(*chis)[0],(*Ndf)[0]));
	tex1->SetNDC();
        tex1->SetTextColor(1);
        tex1->SetTextSize(0.045);
        tex1->SetTextFont(42);
	tex1->Draw("same");

        for(int i=0;i<N-1;i++){
                TH1D *h2_Clone = (TH1D*)histo_exp_norm->Clone();
                h2_Clone->SetFillColor(color_[i]);
                h2_Clone->GetXaxis()->SetRangeUser((*kpoint)[i],(*kpoint)[i+1]);
		h2_Clone->Draw("same");
                }
	c1->SaveAs("stdfit.png");	
}
