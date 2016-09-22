void getmult(){
    const int ncent = 6;
    int centmin[ncent] = {0,6,11,21,41,61};
    int centmax[ncent] = {5,10,20,40,60,80};
    TFile *f=TFile::Open("../../../work/200GeV/output_3corr.root");
    TFile *f1=TFile::Open("../../../work/200GeV/output_3corr_cntbbcsfvtxs.root");
//    for(int icent=0;icent<ncent;icent++){
        //TH2F* hbbcsnfvtxs = (TH2F*)f->Get(Form("hbbcsnfvtxs_%d",icent));
        TH2F* hcentnfvtxs = (TH2F*)f->Get(Form("hcentnfvtxs"));
        TH2F* hcentbbcs = (TH2F*)f1->Get(Form("hcentbbcs"));
    for(int icent=0;icent<ncent;icent++){
        cout<<icent<<endl;
        h1 = hcentnfvtxs->ProjectionY("h1",centmin[icent]+1,centmax[icent]+1);
        h2 = hcentbbcs->ProjectionY("h2",centmin[icent]+1,centmax[icent]+1);
        //hbbcsnfvtxs->GetYaxis()->SetRangeUser(0,50);
        //cout<<"mean bbc mult = "<<hbbcsnfvtxs->GetMean(1)<<endl;
        //cout<<"mean fvtx mult = "<<hbbcsnfvtxs->GetMean(2)<<endl;
        cout<<"mean bbc mult = "<<h2->GetMean(1)<<endl;
        cout<<"mean fvtx mult = "<<h1->GetMean(1)<<endl;
    }
}
