#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
void DrawPFQuality(){
	TString mult = "M150120";
        TString dir10 = "PFcandpt01to10tracknormcpt03to6";
        TString dir1 = "PFcandpt03to6tracknormcpt03to6";
        TString dir2 = "tracknormcpt03to6";
        const int nID=6;
        const int marker[nID] = {20,24,27,30,29,32};
        const int color[nID] = {1,2,3,4,6,3};
        const TString particleleg[nID]={"all PF particles","charged hardons","electron","muon","photon","neutral hadron"};
        TH1D* hpt_PF[nID];
        TH1D* heta_PF[nID];
        TH1D* hphi_PF[nID];
        TH1D* hid_PF;
        for(int iID=0;iID<nID;iID++){
            hpt_PF[iID]=plot(dir1,mult,iID,"pt",color[iID],marker[iID]);
            heta_PF[iID]=plot(dir1,mult,iID,"eta",color[iID],marker[iID]);
            hphi_PF[iID]=plot(dir1,mult,iID,"phi",color[iID],marker[iID]);
        }
            hid_PF=plot(dir1,mult,-1,"id",color[0],marker[0]);
        TH1D* hpt_tr03to6 = plot(dir2,mult,-1,"pt",2,20);
        TH1D* heta_tr03to6 = plot(dir2,mult,-1,"eta",2,20);
        TH1D* hphi_tr03to6 = plot(dir2,mult,-1,"phi",2,20);

	c1 = new TCanvas("c1"," ",800,800);
	c2 = new TCanvas("c2"," ",800,800);
	c3 = new TCanvas("c3"," ",800,800);
	c4 = new TCanvas("c4"," ",600,600);
        TLatex *tl = new TLatex();
        tl->SetTextSize(0.04);
        tl->SetNDC();
//        c1->Divide(3,1);
	c1->cd()->SetLogy();
//	makeMultiPanelCanvas(c1,3,1,0,0,0.25,0.2,0.03);
        gStyle->SetOptFit(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetErrorX(0);
	TH1D* hFrame_pt = new TH1D("hFrame_pt","",100,0,10);
        hFrame_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hFrame_pt->GetXaxis()->SetRangeUser(0,10);
        hFrame_pt->GetYaxis()->SetRangeUser(1,3e9);
        hFrame_pt->GetYaxis()->SetTitle("# of particles");
        hFrame_pt->GetXaxis()->SetTitleOffset(1);
	hFrame_pt->Draw();
        for(int iID=0;iID<nID;iID++){
        hpt_PF[iID]->Draw("Psame");
        }
        hpt_tr03to6->Draw("Psame");
	TLegend *tl1 = new TLegend(0.3,0.12,0.45,0.28);
	TLegend *tl2 = new TLegend(0.55,0.12,0.64,0.28);
	tl1->SetTextSize(0.022);
	tl1->SetBorderSize(0);
	tl1->SetFillColor(0);
	tl2->SetTextSize(0.022);
	tl2->SetBorderSize(0);
	tl2->SetFillColor(0);
	//tl->AddEntry(hpt_PFon,"PFcandidate charged","p");
        for(int iID=0;iID<4;iID++){
	tl1->AddEntry(hpt_PF[iID],particleleg[iID],"p");
        }
        for(int iID=4;iID<nID;iID++){
	tl2->AddEntry(hpt_PF[iID],particleleg[iID],"p");
        }
	tl2->AddEntry(hpt_tr03to6,"tracks 0.3<p_{T}<6.0 (GeV/c)","p");
	tl1->Draw("same");
	tl2->Draw("same");
//	TH1D* ratio = (TH1D*)hpt_PFon->Clone();
//	ratio->Divide(hpt_t03to6);
//	TCanvas *c2 = new TCanvas();
//	ratio->Draw();
        c2->cd()->SetLogy();
	TH1D* hFrame_eta = new TH1D("hFrame_eta","",100,-5,5);
        hFrame_eta->GetXaxis()->SetTitle("#eta");
        hFrame_eta->GetXaxis()->SetRangeUser(-5,5);
        hFrame_eta->GetYaxis()->SetRangeUser(1,3e9);
        hFrame_eta->GetYaxis()->SetTitle("# of particles");
        hFrame_eta->GetXaxis()->SetTitleOffset(1);
        hFrame_eta->Draw();
        for(int iID=0;iID<nID;iID++){
        heta_PF[iID]->Draw("Psame");
        }
        heta_tr03to6->Draw("Psame");
	tl1->Draw("same");
	tl2->Draw("same");
        tl->DrawLatex(0.4,0.8,"0.3<p_{T}<6.0(GeV/c)");
        c3->cd()->SetLogy();
	TH1D* hFrame_phi = new TH1D("hFrame_phi","",100,-4,4);
	hFrame_phi->GetXaxis()->SetTitle("#phi");
        hFrame_phi->GetXaxis()->SetRangeUser(-4,4);
        hFrame_phi->GetYaxis()->SetRangeUser(1e1,3e9);
        hFrame_phi->GetYaxis()->SetTitle("# of particles");
        hFrame_phi->GetXaxis()->SetTitleOffset(1);
        hFrame_phi->Draw();
        for(int iID=0;iID<nID;iID++){
        hphi_PF[iID]->Draw("Psame");
        }
        hphi_tr03to6->Draw("Psame");
	tl1->Draw("same");
	tl2->Draw("same");
        tl->DrawLatex(0.4,0.8,"0.3<p_{T}<6.0(GeV/c)");
        c4->cd()->SetLogy();
	hid_PF->GetXaxis()->SetTitle("particle id");
        hid_PF->GetYaxis()->SetRangeUser(1e1,3e9);
        hid_PF->GetYaxis()->SetTitle("# of particles");
        hid_PF->GetXaxis()->SetTitleOffset(1);
        hid_PF->Draw();
	c1->SaveAs(Form("PFQualify_particlespt_%s.png",dir1.Data()));
	c2->SaveAs(Form("PFQualify_particleseta_%s.png",dir1.Data()));
	c3->SaveAs(Form("PFQualify_particlesphi_%s.png",dir1.Data()));
	c4->SaveAs("PFQualify_particlesid.png");

}
	
TH1D* plot(TString dir, TString Tmult, int id, TString var,int color, int marker){
	TString pf;
	int xbin=0;
	TFile *fhisto = TFile::Open(Form("../%s/%s/histomerged.root",dir.Data(),Tmult.Data()));
        if(id>=0){
	TH1D* histo = (TH1D*)fhisto->Get(Form("D_%d/h%s_%d",xbin,var.Data(),id));
        }
        else{
	TH1D* histo = (TH1D*)fhisto->Get(Form("D_%d/h%s",xbin,var.Data()));
        }
	int rebin=histo->GetNbinsX()/100;
        if(var!="id"){
	histo->Rebin(rebin);
        }
        histo->SetMarkerColor(color);
        histo->SetLineColor(color);
        histo->SetMarkerStyle(marker);
        histo->SetMarkerSize(1);
        cout<<histo->GetEntries()<<endl;
	return histo;
}
