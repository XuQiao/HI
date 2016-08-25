void Draweventeff(){
    gStyle->SetOptStat(kFALSE);
    gStyle->SetErrorX(0);
    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    TH1D* heff = Draw(1,20);
    TH1D* hFrame = new TH1D("","",200,0,200);
    hFrame->SetTitle("");
    c1->SetLogy();
    c1->SetLogx();
    hFrame->GetXaxis()->SetTitle("Multiplicity");
    hFrame->GetXaxis()->SetRangeUser(0,200);
    hFrame->GetYaxis()->SetRangeUser(1e-6,10);
    hFrame->GetYaxis()->SetTitle("event efficiency");
    hFrame->GetYaxis()->SetTitleOffset(1.1);
    TLegend *leg = new TLegend(0.3,0.75,0.85,0.88);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.025);
    TLine *l = new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
    l->SetLineStyle(2);
    l->SetLineWidth(1.2);
    hFrame->Draw();
    heff->Draw("Psame");
    l->Draw("same");
    //leg->Draw("same");
    
    c1->Print(Form("eventeff.png"));
}
TH1D* Draw(int color, int marker){

    TFile *f = new TFile(Form("CorrHijing_recoMC.root"),"ReadOnly");
    TH1D* heff_= (TH1D*)f->Get("hwevent");
    heff = inverse(heff_);
    heff->SetMarkerStyle(marker);
    heff->SetMarkerColor(color);
    heff->SetLineColor(color);
    heff->SetMarkerSize(1.3);
    //f->Close();
    return heff;
}

TH1D* inverse(TH1D* h){
    TH1D* h_ = (TH1D*)h->Clone(h->GetName());
    for(int ibin=0;ibin<h_->GetNbinsX();ibin++){
        if(h->GetBinContent(ibin)==0)
            h_->SetBinContent(ibin,0);
        else{
            h_->SetBinContent(ibin,1.0/h->GetBinContent(ibin));
            h_->SetBinError(ibin,h->GetBinError(ibin)*TMath::Power(h->GetBinContent(ibin),-2));
        }
    }
    return h_;
}
