#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include "par.h"
#include <vector>
void Drawfit(){

        if (!TClass::GetDict("NBD")) {
                gROOT->ProcessLine(".L /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/NBDclass.C+");
        }
/*
	TFile *f = TFile::Open("G.root");
        dir = f->GetDirectory("Gri101");
        def = new TList();
        tl = dir->GetListOfKeys();
        for(int iGlau=0; iGlau<nGlau; iGlau++){
                key = ((TKey*) tl->At(iGlau));
                NBD* temp = (NBD*)key->ReadObj();
                def->Add(temp);
        }
*/	
	int sth=0, Gth=0;
	TFile *f = TFile::Open(Form("G%d%d.root",sth,Gth));
        TList *tl = f->GetListOfKeys();
        TKey *key = ((TKey*) tl->At(0));
        NBD *n0 = (NBD*)key->ReadObj();
	TFile *fdata = TFile::Open(n0->dataname);
	TFile *fGlauber = TFile::Open(n0->Glaubername);
	TH1D *histo_obs = (TH1D*)fdata->Get(n0->histoname);
        int binnum = histo_obs->GetNbinsX();
        double Minx = histo_obs->GetXaxis()->GetXmin();
        double Maxx = histo_obs->GetXaxis()->GetXmax();
	double binsize = (Double_t)(Maxx-Minx)/binnum;
        int xbinmin=(int)((n0->xmin-Minx)/binsize);
        int xbinmax=(int)((n0->xmax-Minx)/binsize);
	TH1D *histo_exp = new TH1D("histo_exp","Simulated distribution;Multiplicity;# of events",binnum,Minx,Maxx);
	Int_t ibin;
		
	TH1D *histo_obs_norm = (TH1D*)histo_obs->Clone();
        histo_obs_norm->Scale(1/histo_obs->Integral(xbinmin,xbinmax));
	cout<<n0->xmin<<endl<<n0->xmax<<endl;
	cout<<xbinmin<<endl<<xbinmax<<endl;
	TF1 *NBD_fun = new 
	TF1("NBD_fun","[0]*TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1]))*TMath::Power([2]/[1],x)/TMath::Power([2]/[1]+1,x+[1])",0,100);
	NBD_fun->SetParameter(0,1);	//[0]: Normalized constant
	NBD_fun->SetParameter(1,n0->kbest);	//[1]: k value
	NBD_fun->SetParameter(2,n0->mubest);	//[2]: mu value
	cout<<n0->kbest<<endl<<n0->mubest<<endl;	
	cout<<fdata<<endl<<fGlauber<<endl<<n0->histoname<<endl;	
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

	TCanvas *c1 = new TCanvas();
	c1->SetLogy();
        gStyle->SetOptStat(kFALSE);

        histo_obs_norm->GetXaxis()->SetRangeUser(0,Maxx);
        histo_obs_norm->GetXaxis()->SetTitle("Mult");
        histo_obs_norm->GetYaxis()->SetTitle("# of events");
        histo_obs_norm->SetTitle("real data vs. Best fit for #mu and k");

	histo_obs_norm->SetLineColor(1);
	histo_obs_norm->SetMarkerStyle(24);
        histo_obs_norm->SetMarkerSize(1.2);
	histo_obs_norm->Draw();

        TLegend *leg = new TLegend(0.56, 0.7, 0.8, 0.9);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0.035);
        leg->SetTextFont(42);
        leg->SetTextSize(0.045);
        leg->AddEntry(histo_obs_norm,"CMS pPb real data","p");
        leg->AddEntry(histo_exp_norm,Form("NBD Fit, #mu=%.3f, k=%.3f",n0->mubest,n0->kbest),"l");
	leg->Draw("same");
        TLatex *tex1= new TLatex(0.2,0.8,Form("#chi^{2}/NDF=%.3f/%d",n0->chis,n0->Ndf));
	tex1->SetNDC();
        tex1->SetTextColor(1);
        tex1->SetTextSize(0.045);
        tex1->SetTextFont(42);
	tex1->Draw("same");

        for(int i=0;i<N-1;i++){
                TH1D *h2_Clone = (TH1D*)histo_exp_norm->Clone();
                h2_Clone->SetFillColor(color_[i]);
//                h2_Clone->GetXaxis()->SetRangeUser(n0->kpoint[i],n0->kpoint[i+1]);
		h2_Clone->Draw("same");
                }

}
