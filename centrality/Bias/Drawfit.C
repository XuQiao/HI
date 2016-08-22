#include "parameter.h"
#include "par.h"
#include <vector>
#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
const int maxNpart = 40;
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
        TVectorD* k0 = (TVectorD*)f->Get(Form("%s/%s/k0",dirname.Data(),name.Data()));
        TVectorD* theta0 = (TVectorD*)f->Get(Form("%s/%s/theta0",dirname.Data(),name.Data()));
        TVectorD* xmin = (TVectorD*)f->Get(Form("%s/%s/xmin",dirname.Data(),name.Data()));
        TVectorD* xmax = (TVectorD*)f->Get(Form("%s/%s/xmax",dirname.Data(),name.Data()));
        TVectorD* thetabest = (TVectorD*)f->Get(Form("%s/%s/thetabest",dirname.Data(),name.Data()));
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
	TH1D *histo_exp = new TH1D("histo_exp","Sithetalated distribution;Multiplicity;Event Fraction",binnum,Minx,Maxx);
	//histo_obs->Sumw2();
	//histo_exp->Sumw2();
	Int_t ibin;
		
	TH1D *histo_obs_norm = (TH1D*)histo_obs->Clone();
        histo_obs_norm->Scale(1/histo_obs->Integral(xbinmin,xbinmax));
	TF1* gammafun[maxNpart];
	for(int iNpart=0;iNpart<maxNpart;iNpart++){
                gammafun[iNpart] = new TF1("gammafun","TMath::GammaDist(x,[0],0,[1])",0,100);
                double k_=(*k0)[0]+(*kbest)[0]*(iNpart-2);
                double theta_=(*theta0)[0]+(*thetabest)[0]*TMath::Log(iNpart-1);
                gammafun[iNpart]->SetParameter(0,k_);   //[1]: k value
                gammafun[iNpart]->SetParameter(1,theta_);       //[2]: theta value
        }
		
	TTree *t = (TTree*) fGlauber->Get("nt_p_Pb");
	Long_t Nevent;

	Nevent = (Long_t) t->GetEntries();

	Long_t Ev;	Int_t Bino;	Double_t Para, Bi_Para, Mult;
	Float_t Ncoll, Npart;
	t->SetBranchAddress("Ncoll",&Ncoll);
	t->SetBranchAddress("Npart",&Npart);
	for(Ev=0; Ev<100000; Ev++){
		if(Ev%100000==0)	 cout<<"Have run "<<Ev<<" events"<<endl;
		t->GetEntry(Ev);
		//Para = 0; //make sure that Para doesn't accumulate through loops
		//for(Bino=0; Bino<Npart; Bino++){
                Para = gammafun[(int)Npart]->GetRandom();
                //         Para += Bi_Para;
		//}	
		histo_exp->Fill(Para);
	}
	Double_t SumEvent, scale;
	SumEvent = histo_exp->Integral(xbinmin,xbinmax);
	scale = 1/SumEvent;
	TH1D *histo_exp_norm = (TH1D*) histo_exp->Clone();
	histo_exp_norm->Scale(scale);

        gStyle->SetOptStat(kFALSE);
	c1 = new TCanvas("c1","",500,800);
	makeMultiPanelCanvas(c1,1,2,-0.1,0.04,0.12,0.12,0.01);
	c1->cd(1)->SetLogy();

        histo_obs_norm->GetXaxis()->SetRangeUser(Minx,100);
        histo_obs_norm->SetMaximum(1.0);
        histo_obs_norm->GetXaxis()->SetTitle("");
        histo_obs_norm->GetYaxis()->SetTitle("Event Fraction");
        histo_obs_norm->GetYaxis()->SetTitleSize(0.06);
        histo_obs_norm->GetYaxis()->SetLabelSize(0.04);
        histo_obs_norm->GetYaxis()->CenterTitle();
        histo_obs_norm->SetTitle("");

	histo_obs_norm->SetLineColor(1);
	histo_obs_norm->SetMarkerStyle(24);
	histo_obs_norm->SetMarkerColor(1);
        histo_obs_norm->SetMarkerSize(1.5);
	histo_obs_norm->Draw("P");

        TLegend *leg = new TLegend(0.46, 0.7, 0.7, 0.9);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0.035);
        leg->SetTextFont(42);
        leg->SetTextSize(0.045);
        leg->AddEntry(histo_obs_norm,"CMS pPb real data","p");
        leg->AddEntry(histo_exp_norm,Form("Gamma Fit, #theta=%.3f, k=%.3f",(*thetabest)[0],(*kbest)[0]),"l");
	leg->Draw("same");
        TLatex *tex1= new TLatex(0.7,0.55,Form("#chi^{2}/NDF=%.1f/%.f",(*chis)[0],(*Ndf)[0]));
        TLatex *tex2= new TLatex(0.7,0.65,Form("fit range %.f<x<%.f",(*xmin)[0],(*xmax)[0]));
	tex1->SetNDC();
	tex2->SetNDC();
        tex1->SetTextColor(1);
        tex2->SetTextColor(1);
        tex1->SetTextSize(0.045);
        tex2->SetTextSize(0.045);
        tex1->SetTextFont(42);
        tex2->SetTextFont(42);
	tex1->Draw("same");
	tex2->Draw("same");

        for(int i=0;i<N-1;i++){
                TH1D *h2_Clone = (TH1D*)histo_exp_norm->Clone();
                h2_Clone->SetFillColor(color_[i]);
                h2_Clone->GetXaxis()->SetRangeUser((*kpoint)[i],(*kpoint)[i+1]);
		h2_Clone->Draw("same");
                }
	gPad->RedrawAxis();
	c1->cd(2);
	
	TH1D* ratio = (TH1D*)histo_obs_norm->Clone();
	ratio->Divide(histo_exp_norm);
	ratio->GetYaxis()->SetTitle("Data/fit");
	ratio->GetXaxis()->SetTitle("HF #Sigma E_{T} |#eta|>4");
	ratio->GetYaxis()->CenterTitle();
	ratio->GetXaxis()->CenterTitle();
	ratio->GetXaxis()->SetTitleSize(0.06);
	ratio->GetXaxis()->SetLabelSize(0.04);
	ratio->GetYaxis()->SetRangeUser(0,1.95);
	ratio->SetMarkerSize(1);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(1);
	ratio->Draw("P");
	TLine *l = new TLine(0,1,100,1);
	l->SetLineStyle(2);
	l->Draw("same");
	c1->SaveAs(Form("%sfit.png",dirname.Data()));	
}
