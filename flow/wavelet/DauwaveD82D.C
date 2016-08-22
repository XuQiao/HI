void DauwaveD82D(){
    TFile *f = TFile::Open("pSTEG.root");
    TTree *t = (TTree*)f->Get("tree");
    const int nbinx=200;
    const int nbiny=200;
    TH2F* h = new TH2F("h","h",nbinx,-5,5,nbiny,-4,4);
    double eta,phi;
    t->SetBranchAddress("eta",&eta);
    t->SetBranchAddress("phi",&phi);
    Int_t N = t->GetEntries();
    for(int imult = 0;imult<N; imult++){
        t->GetEntry(imult);
        h->Fill(eta,phi);
    }
    double binc[nbinx][nbiny];
    for(int i=0;i<nbinx;i++){
        for(int j=0;j<nbiny;j++){
            binc[i][j] = h->GetBinContent(i,j);
        }
    }

}
