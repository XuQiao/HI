#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
#include "parameter.h"
#include <vector>
void Draw3panelfit(){
        gStyle->SetOptStat(kFALSE);
	TString outG="G0.root";
	TString dirname;
	c1 = new TCanvas("c1"," ",1200,800);
	makeMultiPanelCanvas(c1,3,2,0,0,0.2,0.1,0.03);
	for(int sth=0;sth<3;sth++){
	int cenvar=0, Gth=0;
	//TFile *f = TFile::Open(Form("%s/%s",cenvardir[cenvar].Data(),outG.Data()));
	TFile *f = TFile::Open(Form("%s",outG.Data()));
	if(sth==0)	dirname = "std";
        else if(sth==1)	dirname ="Gri055";
        else dirname ="Gri101";
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
	TH1D *histo_exp = new TH1D("histo_exp","Simulated distribution;Multiplicity;Event Fraction",binnum,Minx,Maxx);
	//histo_obs->Sumw2();
	//histo_exp->Sumw2();
	Int_t ibin;
		
	TH1D *histo_obs_norm = (TH1D*)histo_obs->Clone();
        histo_obs_norm->Scale(1/histo_obs->Integral(xbinmin,xbinmax));
        TF1 *gammafun = new TF1("gammafun","TMath::GammaDist(x,[0],0,[1])",0,100);
	gammafun->SetParameter(0,(*kbest)[0]);	//[1]: k value
	gammafun->SetParameter(1,(*thetabest)[0]);	//[2]: theta value
		
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
                         Bi_Para = gammafun->GetRandom();
                         Para += Bi_Para;
		}	
		histo_exp->Fill(Para);
	}
	Double_t SumEvent, scale;
	SumEvent = histo_exp->Integral(xbinmin,xbinmax);
	scale = 1/SumEvent;
	TH1D *histo_exp_norm = (TH1D*) histo_exp->Clone();
	histo_exp_norm->Scale(scale);
	c1->cd(sth+1)->SetLogy();
	/*c1->cd(sth+1)->SetFrameFillColor(0);
        c1->cd(sth+1)->SetFrameBorderSize(0);
        c1->cd(sth+1)->SetFrameBorderMode(0);
        if(sth==0)c1->cd(sth+1)->SetLeftMargin(0.12);
        else c1->cd(sth+1)->SetLeftMargin(0);
        c1->cd(sth+1)->SetBottomMargin(0.1);
        c1->cd(sth+1)->SetTopMargin(0.02);
        if(sth==2)c1->cd(sth+1)->SetRightMargin(0.05);
        else c1->cd(sth+1)->SetRightMargin(0);
	*/
        c1->cd(sth+1)->SetTicks(-1);

	
	TH1D* hFrame = (TH1D*) histo_exp->Clone();
	hFrame->Reset();
        if(cenvar==0)hFrame->GetXaxis()->SetRangeUser(0,100);
        else if(cenvar==1) hFrame->GetXaxis()->SetRangeUser(0,70);
	else hFrame->GetXaxis()->SetRangeUser(0,200);
        if(cenvar==2)hFrame->SetMaximum(0.5);
        else hFrame->SetMaximum(1.0);
        hFrame->GetXaxis()->SetTitle("");
        if(sth==0)	hFrame->GetYaxis()->SetTitle("Event Fraction");
	else hFrame->GetYaxis()->SetTitle("");
        hFrame->SetTitle("");
	hFrame->GetXaxis()->CenterTitle(0);
        hFrame->GetYaxis()->CenterTitle(1);
        if(sth==0)hFrame->GetYaxis()->SetTitleOffset(1.2);
        else hFrame->GetYaxis()->SetTitleOffset(0);
        if(sth==0){
		hFrame->GetYaxis()->SetLabelSize(0.05);
		hFrame->GetYaxis()->SetTitleSize(0.06);
	}
        else{
		hFrame->GetYaxis()->SetLabelSize(0);
		hFrame->GetYaxis()->SetTitleSize(0);
	}
        hFrame->GetXaxis()->SetLabelSize(0.05);
        hFrame->GetXaxis()->SetTitleSize(0.05);
        hFrame->GetXaxis()->SetTitleOffset(1.0);

	hFrame->DrawCopy();

	histo_obs_norm->SetLineColor(1);
	histo_obs_norm->SetLineWidth(0.5);
	histo_obs_norm->SetMarkerStyle(24);
	histo_obs_norm->SetMarkerColor(1);
        histo_obs_norm->SetMarkerSize(1.1);

        TLegend *leg = new TLegend(0.50, 0.7, 0.75, 0.80);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0.035);
        leg->SetTextFont(42);
        leg->SetTextSize(0.045);
        leg->AddEntry(histo_obs_norm,"CMS pPb real data","p");
        leg->AddEntry(histo_exp_norm,Form("Gamma Fit","l"));
	if(sth==2) leg->Draw("same");
        TLatex *tex1= new TLatex(0.18,0.89,Form("#chi^{2}/NDF=%.1f/%.f",(*chis)[0],(*Ndf)[0]));
        TLatex *tex2= new TLatex(0.18,0.85,Form("fit range (%.f , %.f)",(*xmin)[0],(*xmax)[0]));
        TLatex *tex3= new TLatex(0.60,0.83,Form("#theta=%.3f, k=%.3f",(*thetabest)[0],(*kbest)[0]));
        TLatex *tex4= new TLatex(0.60,0.87,Grilist[sth]);
	tex1->SetNDC();
	tex2->SetNDC();
	tex3->SetNDC();
	tex4->SetNDC();
        tex1->SetTextColor(1);
        tex2->SetTextColor(1);
        tex3->SetTextColor(1);
        tex4->SetTextColor(1);
        tex1->SetTextSize(0.05);
        tex2->SetTextSize(0.05);
        tex3->SetTextSize(0.05);
        tex4->SetTextSize(0.05);
        tex1->SetTextFont(42);
        tex2->SetTextFont(42);
        tex3->SetTextFont(42);
        tex4->SetTextFont(42);
	tex1->Draw("same");
	tex2->Draw("same");
	tex3->Draw("same");
	tex4->Draw("same");

        for(int i=0;i<kpoint->GetUpb();i++){
                TH1D *h2_Clone = (TH1D*)histo_exp_norm->Clone();
                h2_Clone->SetFillColor(color_[i]);
                h2_Clone->GetXaxis()->SetRangeUser((*kpoint)[i],(*kpoint)[i+1]);
		h2_Clone->Draw("same");
                }
	histo_obs_norm->Draw("Psame");
	gPad->RedrawAxis();
	c1->cd(sth+4);
	TH1D* histo_ratio = (TH1D*)histo_exp_norm->Clone();
	histo_ratio->Divide(histo_obs_norm);
	histo_ratio->SetMarkerStyle(20);
	histo_ratio->SetMarkerSize(1.1);
	hFrame->GetYaxis()->SetRangeUser(0,4.5);
        hFrame->GetXaxis()->SetTitle(cenvariable[cenvar]);
        if(sth==0)	hFrame->GetYaxis()->SetTitle("fit/data");
	hFrame->DrawCopy();
	histo_ratio->Draw("Psame");
	TLine *l = new TLine(0,1,histo_ratio->GetXaxis()->GetXmax(),1);
	l->SetLineStyle(2);
	l->SetLineWidth(1.2);
	l->Draw("same");
	}
	c1->Print(Form("3panelfit.png"));	
	c1->Print(Form("3panelfit.pdf"));	
}
