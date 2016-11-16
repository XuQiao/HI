void DrawFVTXclus(){
    TFile *fin = TFile::Open("output_perform.root");
    TH2F* hfvtxclusxy[8];
    TH2F* hfvtxclusrphi[8];
    TH2F* hfvtxclusetaphi[8];
    gStyle->SetOptStat(kFALSE);
    for(int ilayer=0;ilayer<8;ilayer++){
        hfvtxclusxy[ilayer] = (TH2F*)fin->Get(Form("fvtxclusxy_%d",ilayer));
        hfvtxclusrphi[ilayer] = (TH2F*)fin->Get(Form("fvtxclusrphi_%d",ilayer));
        hfvtxclusetaphi[ilayer] = (TH2F*)fin->Get(Form("fvtxclusetaphi_%d",ilayer));

        TCanvas *c1 = new TCanvas("c1","c1",600,600);
        hfvtxclusxy[ilayer]->GetXaxis()->SetTitle("fvtx X(cm)");
        hfvtxclusxy[ilayer]->GetYaxis()->SetTitle("fvtx Y(cm)");
        hfvtxclusxy[ilayer]->Draw("colz");
        c1->Print(Form("fig/fvtxclusxy_%d.png",ilayer));
        TCanvas *c2 = new TCanvas("c2","c2",600,600);
        hfvtxclusrphi[ilayer]->GetXaxis()->SetTitle("fvtx r(cm)");
        hfvtxclusrphi[ilayer]->GetXaxis()->SetRangeUser(2,20);
        hfvtxclusrphi[ilayer]->GetYaxis()->SetTitle("fvtx phi(rad)");
        hfvtxclusrphi[ilayer]->Draw("colz");
        c2->Print(Form("fig/fvtxclusrphi_%d.png",ilayer));
        TCanvas *c3 = new TCanvas("c3","c3",600,600);
        hfvtxclusetaphi[ilayer]->GetXaxis()->SetTitle("fvtx eta");
        hfvtxclusetaphi[ilayer]->GetYaxis()->SetTitle("fvtx phi(rad)");
        hfvtxclusetaphi[ilayer]->Draw("colz");
        c3->Print(Form("fig/fvtxclusetaphi_%d.png",ilayer));
    }
}
