#include "../parameter.h"
#include "par.h"
#include <vector>
void Drawfit(){
        gStyle->SetErrorX(0);
	int sth=0, Gth=0;
	TFile *f = TFile::Open(outG);
        TString dirname;
	if(sth==0){dirname = "std";}
        else if(sth==1){dirname ="Gri055";}
        else { dirname ="Gri101";}
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
        TVectorD* fbest = (TVectorD*)f->Get(Form("%s/%s/fbest",dirname.Data(),name.Data()));
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
		
	TTree *t = (TTree*) fGlauber->Get("nt_Pb_Pb");
	Long_t Nevent;

	Nevent = (Long_t) t->GetEntries();

	Long_t Ev;	Int_t Bino;	Double_t Para, Bi_Para, Mult;
	Float_t Ncoll, Npart;
	t->SetBranchAddress("Ncoll",&Ncoll);
	t->SetBranchAddress("Npart",&Npart);

	for(Ev=0; Ev<Nevent; Ev++){
		if(Ev%100000==0)	 cout<<"Have run "<<Ev<<" events"<<endl;
		t->GetEntry(Ev);
		Para = 0; //make sure that Para doesn't accumulate through loops

                int Nanc = (int) ((*fbest)[0]*Npart + (1-(*fbest)[0]) * Ncoll);
		for(Bino=0; Bino<Nanc; Bino++){
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

        histo_obs_norm->GetXaxis()->SetRangeUser(0,5000);
        histo_obs_norm->SetMaximum(1.0);
        histo_obs_norm->GetXaxis()->SetTitle("HF #Sigma E_{T} |#eta|>3");
        histo_obs_norm->GetYaxis()->SetTitle("Event Fraction");
        histo_obs_norm->SetTitle("");

	histo_obs_norm->SetLineColor(1);
	histo_obs_norm->SetMarkerStyle(24);
	histo_obs_norm->SetMarkerColor(1);
        histo_obs_norm->SetMarkerSize(1.0);
	//histo_obs_norm->Rebin(5);
	//histo_obs_norm->Scale(1./5.);
	histo_obs_norm->DrawCopy("P");

        TLegend *leg = new TLegend(0.46, 0.75, 0.8, 0.9);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.045);
        leg->AddEntry(histo_obs_norm,"CMS PbPb real data","p");
        leg->AddEntry(histo_exp_norm,Form("NBD Fit"));
	leg->Draw("same");
        TLatex *tex = new TLatex(0.5,0.70, Form("#mu=%.2f, k=%.2f, f=%.2f",(*mubest)[0],(*kbest)[0],(*fbest)[0]));
        TLatex *tex1= new TLatex(0.2,0.8,Form("#chi^{2}/NDF=%.1f/%.f",(*chis)[0],(*Ndf)[0]));
        TLatex *tex2= new TLatex(0.6,0.65,Form("fit range %.f<x<%.f",(*xmin)[0],(*xmax)[0]));
	tex->SetNDC();
	tex1->SetNDC();
	tex2->SetNDC();
        tex->SetTextColor(1);
        tex1->SetTextColor(1);
        tex2->SetTextColor(1);
        tex->SetTextSize(0.045);
        tex1->SetTextSize(0.045);
        tex2->SetTextSize(0.045);
        tex->SetTextFont(42);
        tex1->SetTextFont(42);
        tex2->SetTextFont(42);
	tex->Draw("same");
	tex1->Draw("same");
	tex2->Draw("same");

        for(int i=0;i<N-1;i++){
                TH1D *h2_Clone = (TH1D*)histo_exp_norm->Clone();
                h2_Clone->SetFillColor(color_[i]);
                h2_Clone->GetXaxis()->SetRangeUser((*kpoint)[i],(*kpoint)[i+1]);
		h2_Clone->DrawCopy("same");
                }
	histo_obs_norm->Draw("Psame");
        gPad->RedrawAxis();
	
        TCanvas *c2 = new TCanvas();
	c2->SetLogy();
        histo_obs_norm->GetXaxis()->SetRangeUser(0,500);
        histo_obs_norm->SetMaximum(1.0);
        histo_obs_norm->GetXaxis()->SetTitle("HF #Sigma E_{T} |#eta|>3");
        histo_obs_norm->GetYaxis()->SetTitle("Event Fraction");
        histo_obs_norm->SetTitle("");

	histo_obs_norm->SetLineColor(1);
	histo_obs_norm->SetMarkerStyle(24);
	histo_obs_norm->SetMarkerColor(1);
        histo_obs_norm->SetMarkerSize(0.8);
	histo_obs_norm->Draw("P");

        for(int i=0;i<N-1;i++){
                TH1D *h2_Clone = (TH1D*)histo_exp_norm->Clone();
                h2_Clone->SetFillColor(color_[i]);
                h2_Clone->GetXaxis()->SetRangeUser((*kpoint)[i],(*kpoint)[i+1]);
		h2_Clone->Draw("same");
                }
	histo_obs_norm->Draw("Psame");
        gPad->RedrawAxis();
	c1->SaveAs(Form("%sfit_%s.png",dirname.Data(),name.Data()));	
	c2->SaveAs(Form("%sfit_%s_zoom.png",dirname.Data(),name.Data()));	
}
