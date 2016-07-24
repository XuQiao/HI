#include "/home/xuq7/CMSSW_6_2_3_patch1/src/Centrality/NBD/parameter.h"
void Draw3paneleff(){
        gStyle->SetOptStat(kFALSE);
	TString outG="G1.root";
	TCanvas *c1  = new TCanvas("","",1000,400);
	c1->Divide(3);
	gStyle->SetErrorX(0);
	double hfbin[]={0,1,2,3,4,6,8,10,13,16,20,25,30,40,55,70,90};
	int nhfbin = 16;
	TString name, dirname;
	TFile *fDSeff = TFile::Open("/scratch/xuq7/Centrality/pPbHijing_EffCorr_forNBD.root");
        TFile *ftreff = TFile::Open("/scratch/xuq7/Centrality/pPbHist_Hijing_TrandEs.root");
        TH1D* hbef = (TH1D*)ftreff->Get("hHFEnergy4");
        TH1D* rehbef = (TH1D*)hbef->Rebin(nhfbin,"rehHFEnergy4",hfbin);
        TH1D* haft = (TH1D*)ftreff->Get("hHFEnergy4_tr");
        TH1D* rehaft = (TH1D*)haft->Rebin(nhfbin,"rehHFEnergy4_tr",hfbin);
        TGraphAsymmErrors *gtreff = new TGraphAsymmErrors();
        gtreff->BayesDivide(rehaft,rehbef);

        TGraphAsymmErrors *geff = (TGraphAsymmErrors*)fDSeff->Get("regEffHF4");
        for(int i=0;i<geff->GetN();i++){
                geff->SetPointEXlow(i,0);
                geff->SetPointEXhigh(i,0);
                gtreff->SetPointEXlow(i,0);
                gtreff->SetPointEXhigh(i,0);
        }
	
	geff->SetTitle("");
        geff->SetMarkerStyle(33);
        geff->SetMarkerColor(5);
        geff->SetMarkerSize(1.3);

	gtreff->SetTitle("");
        gtreff->SetMarkerStyle(21);
        gtreff->SetMarkerColor(6);
        gtreff->SetMarkerSize(1.3);

	for(int cenvar=0;cenvar<3;cenvar++){
	TFile *f = TFile::Open(Form("%s/%s",cenvardir[cenvar].Data(),outG.Data()));
	c1->cd(cenvar+1);
        TObjString* dataname = (TObjString*)f->Get(Form("dataname"));
        TObjString* histoname = (TObjString*)f->Get(Form("histoname"));
	TFile *fdata = TFile::Open(dataname->GetString());
	TH1D *histo_obs = (TH1D*)fdata->Get(histoname->GetString());
	histo_obs->Sumw2();
	TH1D* hFrame = histo_obs->Clone();
	hFrame->Reset();
	hFrame->SetTitle("");
        hFrame->SetMaximum(1.2);
        hFrame->SetMinimum(0);
        if(cenvar==0)hFrame->GetYaxis()->SetTitle("ratio");
	else hFrame->GetYaxis()->SetTitle("");
	hFrame->GetXaxis()->SetRangeUser(0,100);
        hFrame->GetXaxis()->CenterTitle(0);
        hFrame->GetYaxis()->CenterTitle(1);
        hFrame->GetYaxis()->SetTitleOffset(1.1);
        hFrame->GetXaxis()->SetTitleOffset(1.1);
        hFrame->GetXaxis()->SetTitleSize(0.056);
        hFrame->GetYaxis()->SetTitleSize(0.056);
        hFrame->GetXaxis()->SetLabelSize(0.05);
        hFrame->GetYaxis()->SetLabelSize(0.05);

        int binnum = histo_obs->GetNbinsX();
        double Minx = histo_obs->GetXaxis()->GetXmin();
        double Maxx = histo_obs->GetXaxis()->GetXmax();
        double binsize = (Double_t)(Maxx-Minx)/binnum;
	hFrame->Draw();
        TLegend *leg = new TLegend(0.40, 0.2, 0.65, 0.4);
	
	for(int sth=0;sth<3;sth++){
        int Gth=0;
        if(sth==0){dirname = "std";}
        else if(sth==1){dirname ="Gri055";}
        else {dirname ="Gri101";}
        hFrame->GetXaxis()->SetTitle(cenvariable[cenvar]);

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
        int xbinmin=(int)(((*xmin)[0]-Minx)/binsize);
        int xbinmax=(int)(((*xmax)[0]-Minx)/binsize);
	TFile *fGlauber = TFile::Open(Glaubername->GetString());
	TH1D *histo_exp = new TH1D("histo_exp","Simulated distribution;Multiplicity;Event Fraction",binnum,Minx,Maxx);
	histo_exp->Sumw2();
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


	rehisto_obs_norm = (TH1D*)histo_obs_norm->Rebin(nhfbin,"rehisto_obs_norm",hfbin);
	normalizeByBinWidth(rehisto_obs_norm);
	rehisto_exp_norm = (TH1D*)histo_exp_norm->Rebin(nhfbin,"rehisto_exp_norm",hfbin);
	normalizeByBinWidth(rehisto_exp_norm);
	TH1D* ratio = (TH1D*)rehisto_obs_norm->Clone("ratio");
	ratio->Divide(rehisto_exp_norm);

	ratio->SetTitle("");
	ratio->SetLineColor(1);
	ratio->SetMarkerStyle(20);
	ratio->SetMarkerColor(color_[sth]);
        ratio->SetMarkerSize(1.3);
	ratio->Draw("Psame");
	geff->Draw("Psame");
	gtreff->Draw("Psame");
        leg->AddEntry(ratio,Form("data/fit %s",Grilist[sth].Data()),"p");
	}
	
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0.035);
        leg->SetTextFont(42);
        leg->SetTextSize(0.045);
        leg->AddEntry(geff,"DS efficiency","p");
        leg->AddEntry(gtreff,"Event selection efficiency","p");
	if(cenvar==2) leg->Draw("same");
	
	TLine *l = new TLine(0,1,90,1);
	l->SetLineStyle(2);
	l->Draw("same");
	}
	c1->SaveAs(Form("ratiovseff.png"));	
	c1->SaveAs(Form("ratiovseff.pdf"));	
}

TH1 * normalizeByBinWidth(TH1 *histo) {
for(int i = 1; i <= histo->GetNbinsX(); i++) {
float content = histo->GetBinContent(i);
float error = histo->GetBinError(i);
histo->SetBinContent(i,content/histo->GetBinWidth(i));
histo->SetBinError(i,error/histo->GetBinWidth(i));
}
return histo ;
}
