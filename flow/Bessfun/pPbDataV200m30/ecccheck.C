void ecccheck(){
    int ifile=0;
    TFile *f = TFile::Open(Form("/cms/store/user/qixu/flow/NewSTEG/pPbDataV205m300/vndata_50k_%d.root",ifile));
    TTree *t = (TTree*)f->Get("tree");
    Float_t eta[1000];
    Float_t phi[1000];
    int n;
    t->SetBranchAddress("etag",eta);
    t->SetBranchAddress("phig",phi);
    t->SetBranchAddress("n",&n);
    //Int_t N = t->GetEntries();
    Int_t N = 10000;
    double sum=0;
    for(int i=0;i<N;i++){
        t->GetEntry(i);
        for(int imult = 0;imult<n; imult++){
            sum+=TMath::Cos(2*phi[imult]);
    }
    }
    cout<<sum/N/300<<endl;
}

