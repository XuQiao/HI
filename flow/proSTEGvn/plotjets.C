void plotjets(){
    TFile *f = TFile::Open("/store/user/qixu/flow/STEGwithnf/pPbDataV205m100/Varmults10b10/vndata_50k_0.root");
    int ievent = 1000;
    int n;
    float etag[10000];
    float phig[10000];
    float ptg[10000];
    TH2F* hetaphi = new TH2F("hetaphi","",100,-4,4,200,0,6.5);
    TTree *t = f->Get("tree");
    t->SetBranchAddress("n",&n);
    t->SetBranchAddress("etag",&etag);
    t->SetBranchAddress("phig",&phig);
    t->GetEntry(ievent);
    for(int ip = 0;ip<n;ip++){
        hetaphi->Fill(etag[ip],phig[ip]);
    }
    hetaphi->GetXaxis()->SetTitle("#eta");
    hetaphi->GetYaxis()->SetTitle("#phi");
    TCanvas *c1 = new TCanvas();
    //c1->SetLogy();
    hetaphi->Draw("colz");
    c1->Print("hetaphi.png");
}


