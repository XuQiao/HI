void Drawtrackeff(){
    gStyle->SetOptStat(kFALSE);
    TString type="eta";
    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    TH1D* heff = Draw(type,1,20);
    TH1D* hFrame = new TH1D("","",2030,-3,200);
    hFrame->SetTitle("");
    if(type=="pt"){
        c1->SetLogx();
   //   c1->SetLogy();
        hFrame->GetXaxis()->SetTitle("p_{T}");
        hFrame->GetXaxis()->SetRangeUser(0,200);
    }
    else if(type=="eta"){
        hFrame->GetXaxis()->SetTitle("#eta");
        hFrame->GetXaxis()->SetRangeUser(-3,3);
    }
    hFrame->GetYaxis()->SetRangeUser(0,1);
    hFrame->GetYaxis()->SetTitle("tracking efficiency");
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
    
    c1->Print(Form("trackeff_%s.png",type.Data()));
}
TH1D* Draw(TString type, int color, int marker){

    TFile *f = new TFile(Form("trackCorr.root"),"ReadOnly");
    TH2F* heptrack= (TH2F*)f->Get("heptrack");
    if(type=="pt"){
        TH1D* heff = (TH1D*)heptrack->ProjectionX("heff",heptrack->GetYaxis()->FindBin(-2.4),heptrack->GetYaxis()->FindBin(2.4)-1,"i");
        heff->Scale(1.0/(heptrack->GetYaxis()->FindBin(2.4)-heptrack->GetYaxis()->FindBin(-2.4)));
        heff->Draw();
        cout<<(heptrack->GetYaxis()->FindBin(2.4)-heptrack->GetYaxis()->FindBin(-2.4))<<endl;
    }
    else if(type=="eta"){
        TH1D* heff = (TH1D*)heptrack->ProjectionY("heff",heptrack->GetXaxis()->FindBin(0.1),-1,"i");
        heff->Scale(1.0/(heptrack->GetNbinsX()-heptrack->GetXaxis()->FindBin(0.1)+1));
    }
    heff->SetMarkerStyle(marker);
    heff->SetMarkerColor(color);
    heff->SetLineColor(color);
    heff->SetMarkerSize(1.5);
  //  f->Close();
    return heff;
}



