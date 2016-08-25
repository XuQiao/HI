void Drawdatadiff(){
	TFile *f1 = TFile::Open("double_side/G0.root");
	TFile *f2 = TFile::Open("PlusEta4_Rebin/G0.root");
        TObjString* dataname1 = (TObjString*)f1->Get(Form("dataname"));
        TObjString* dataname2 = (TObjString*)f2->Get(Form("dataname"));
        TObjString* histoname1 = (TObjString*)f1->Get(Form("histoname"));
        TObjString* histoname2 = (TObjString*)f2->Get(Form("histoname"));
	TFile *fdata1 = TFile::Open(dataname1->GetString());
	TFile *fdata2 = TFile::Open(dataname2->GetString());
	TH1D *histo_obs1 = (TH1D*)fdata1->Get(histoname1->GetString());
	histo_obs1->Sumw2();
	TH1D *histo_obs2 = (TH1D*)fdata2->Get(histoname2->GetString());
	histo_obs2->Sumw2();
	TString dirname = "std";
	TString name =	"G0";
        TVectorD* xmin = (TVectorD*)f1->Get(Form("%s/%s/xmin",dirname.Data(),name.Data()));
        TVectorD* xmax = (TVectorD*)f1->Get(Form("%s/%s/xmax",dirname.Data(),name.Data()));
        int binnum = histo_obs1->GetNbinsX();
        double Minx = histo_obs1->GetXaxis()->GetXmin();
        double Maxx = histo_obs1->GetXaxis()->GetXmax();
        double binsize = (Double_t)(Maxx-Minx)/binnum;
        int xbinmin=(int)(((*xmin)[0]-Minx)/binsize);
        int xbinmax=(int)(((*xmax)[0]-Minx)/binsize);

        TH1D *histo_obs_norm1 = (TH1D*)histo_obs1->Clone();
        TH1D *histo_obs_norm2 = (TH1D*)histo_obs2->Clone();
        histo_obs_norm1->Scale(1/histo_obs1->Integral(xbinmin,xbinmax));
        histo_obs_norm2->Scale(1/histo_obs2->Integral(xbinmin,xbinmax));

/*	double hfbin[]={0,1,2,3,4,6,8,10,13,16,20,25,30,40,55,70,90};
	int nhfbin = 16;
	rehisto_obs_norm = (TH1D*)histo_obs_norm->Rebin(nhfbin,"rehisto_obs_norm",hfbin);
	normalizeByBinWidth(rehisto_obs_norm);
	rehisto_exp_norm = (TH1D*)histo_exp_norm->Rebin(nhfbin,"rehisto_exp_norm",hfbin);
	normalizeByBinWidth(rehisto_exp_norm);
*/
        TCanvas *c1 = new TCanvas();
	c1->SetLogy();
        gStyle->SetOptStat(kFALSE);
	histo_obs_norm1->SetLineColor(1);
	histo_obs_norm2->SetLineColor(2);
	histo_obs_norm1->SetMarkerStyle(20);
	histo_obs_norm2->SetMarkerStyle(33);
	histo_obs_norm1->SetMarkerColor(1);
	histo_obs_norm2->SetMarkerColor(2);
        histo_obs_norm1->SetMarkerSize(1.3);
        histo_obs_norm2->SetMarkerSize(1.3);
	histo_obs_norm2->SetTitle("");
	histo_obs_norm2->GetYaxis()->SetTitle("Event Fraction");
	histo_obs_norm2->GetXaxis()->SetTitle("HF deposit Energy");
	histo_obs_norm2->Draw("P");	
	histo_obs_norm1->Draw("Psame");	

        TLegend *leg = new TLegend(0.66, 0.4, 0.85, 0.6);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0.035);
        leg->SetTextFont(42);
        leg->SetTextSize(0.045);
        leg->AddEntry(histo_obs_norm1,"HF #Sigma E_{T} |#eta|>4","p");
        leg->AddEntry(histo_obs_norm2,"HF #Sigma E_{T} #eta>4","p");
	leg->Draw("same");
	
        TLine *l1 = new TLine(12.61,0,12.61,0.3);
        TLine *l3 = new TLine(18.99,0,18.99,0.3);
        TLine *l2 = new TLine(6.63,0,6.63,0.3);
        TLine *l4 = new TLine(11.27,0,11.27,0.3);
        l1->SetLineStyle(2);
        l1->SetLineWidth(2);
        l1->SetLineColor(1);
        l3->SetLineStyle(2);
        l3->SetLineWidth(2);
        l3->SetLineColor(1);
        l2->SetLineStyle(2);
        l2->SetLineWidth(2);
        l2->SetLineColor(2);
        l4->SetLineStyle(2);
        l4->SetLineWidth(2);
        l4->SetLineColor(2);
        l1->Draw("same");
        l3->Draw("same");
        l2->Draw("same");
        l4->Draw("same");
	c1->SaveAs("datacompare.png");	

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
