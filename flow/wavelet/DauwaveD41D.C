void DauwaveD41D(){
    int ifile=0;
    TFile *f = TFile::Open(Form("/cms/store/user/qixu/flow/STEG/pPbDataV200m100/vndata_50k_%d.root",ifile));
    TTree *t = (TTree*)f->Get("tree");
    const int nbinx=200;
    const int nbiny=200;
    TH2F* h = new TH2F("h","h",nbinx,-5,5,nbiny,0,8);
    Float_t eta[1000],phi[1000];
    int n;
    t->SetBranchAddress("etag",eta);
    t->SetBranchAddress("phig",phi);
    t->SetBranchAddress("n",&n);
    Int_t N = t->GetEntries();
    for(int ievt=0;ievt<N;ievt++){
        t->GetEntry(ievt);
       for(int imult = 0;imult<n; imult++){
        h->Fill(eta[imult],phi[imult]);
       }
    }
    h->Draw("surf1");
    double binc[nbinx][nbiny];
    for(int i=0;i<nbinx;i++){
        for(int j=0;j<nbiny;j++){
            binc[i][j] = h->GetBinContent(i,j);
        }
    }

}
