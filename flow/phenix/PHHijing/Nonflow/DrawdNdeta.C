void DrawdNdeta(){
    gStyle->SetOptStat(kTRUE);
    TFile *f = TFile::Open("mergedouthisto_pAu.root");
    //TFile *f = TFile::Open("outhisto.root");
    TString coll = "pAu";
    TH1F* hKF = (TH1F*)f->Get("hKF");
    TH1F* hy = (TH1F*)f->Get("hy");
    TH1F* hy1 = (TH1F*)f->Get("hy1");
    TH1F* heta = (TH1F*)f->Get("heta");
    TH1F* hMass = (TH1F*)f->Get("hMass");
    TH1F* hn = (TH1F*)f->Get("hn");
    TH1F* hnch = (TH1F*)f->Get("hnch");
    TH1F* hphi = (TH1F*)f->Get("hphi");
    TH2F* hpteta = (TH2F*)f->Get("hpteta");
    TH1F* hpt = (TH1F*)hpteta->ProjectionX();
    TCanvas *c1 = new TCanvas("c1");
    hy->GetXaxis()->SetTitle("y");
    hy->GetYaxis()->SetTitle("# of particles");
    hy->Draw();
    TLatex t;
    t.SetNDC();
 //   t.DrawLatex(0.1,0.4,Form("#pi_{0}=%d, #pi^{#pm}=%d,K^{#pm}=%d",18308,32043,3714));
 //   t.DrawLatex(0.1,0.3,Form("#pi_{0}/#pi^{#pm}=%.2f,K^{#pm}/#pi^{#pm}=%.2f",18308./32043,3714./32043));

    c1->Print(Form("fig/%s_hy.png",coll.Data()));
    TCanvas *c2 = new TCanvas("c2");
    c2->cd();
    hy1->GetXaxis()->SetTitle("y");
    hy1->GetYaxis()->SetTitle("# of particles");
    hy1->Draw();
 //   t.DrawLatex(0.1,0.4,Form("#pi_{0}=%d, #pi_{#pm}=%d,K^{#pm}=%d",11676,20281,2530));
 //   t.DrawLatex(0.1,0.3,Form("#pi_{0}/#pi_{#pm}=%.2f,K^{#pm}/#pi^{#pm}=%.2f",11676./20281,2530./20281));
    c2->Print(Form("fig/%s_hy1.png",coll.Data()));

    TCanvas *c3 = new TCanvas("c3");
    c3->SetLogy();
    hKF->GetXaxis()->SetTitle("particle id");
    hKF->GetYaxis()->SetTitle("# of particles");
    hKF->Draw();
    c3->Print(Form("fig/%s_pid.png",coll.Data()));

    TCanvas *c4 = new TCanvas("c4");
    heta->GetXaxis()->SetTitle("#eta");
    heta->GetYaxis()->SetTitle("# of particles");
    heta->Draw();
    c4->Print(Form("fig/%s_eta.png",coll.Data()));

    TCanvas *c5 = new TCanvas("c5");
    hMass->GetXaxis()->SetTitle("mass(GeV)");
    hMass->GetYaxis()->SetTitle("# of particles");
    hMass->Draw();
    c5->Print(Form("fig/%s_mass.png",coll.Data()));

    TCanvas *c6 = new TCanvas("c6");
    c6->SetLogy();
    hn->GetXaxis()->SetTitle("N");
    hn->GetYaxis()->SetTitle("# of events");
    hn->Draw();
    c6->Print(Form("fig/%s_N.png",coll.Data()));

    TCanvas *c7 = new TCanvas("c7");
    c7->SetLogy();
    hnch->GetXaxis()->SetTitle("N_{ch}");
    hnch->GetYaxis()->SetTitle("# of events");
    hnch->Draw();
    c7->Print(Form("fig/%s_Nch.png",coll.Data()));

    TCanvas *c8 = new TCanvas("c8");
    c8->SetLogy();
    hpt->GetXaxis()->SetTitle("particle p_{T}(GeV/c)");
    hpt->GetYaxis()->SetTitle("# of particles");
    hpt->Draw();
    c8->Print(Form("fig/%s_pt.png",coll.Data()));

    TCanvas *c9 = new TCanvas("c9");
    hphi->GetXaxis()->SetRangeUser(0,6.2);
    hphi->GetYaxis()->SetRangeUser(0,6.2);
    hphi->GetXaxis()->SetTitle("#phi");
    hphi->GetYaxis()->SetTitle("# of particles");
    hphi->Draw();
    c9->Print(Form("fig/%s_phi.png",coll.Data()));
    
}
