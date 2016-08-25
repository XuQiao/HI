void handmerge(){
    TString dir = "/scratch/xuq7/phenix/pAu/Perform/AnapAumbst_";
    TH2F* hcentcnteta; 
    TH2F* hcentcnteta1; 
    TH2F* hcentcnteta11; 
    TH2F* hcentcnteta12; 
    TH2F* hcentcnteta21; 
    TH2F* hcentcnteta22; 
    TH2F* hcentcntpt; 
    TH2F* hcentntrk; 
    TH2F* hcentntrk1; 
    TFile *f[1513];
    for(int j=5;j<10;j++){
    for(int i=0;i<153;i++){
        cout<<i<<endl;
        if(j*153+i>=1513) continue;
        f[j*153+i] = TFile::Open(Form("%s%d.root",dir.Data(),j*153+i));
        //try{
        /*
            TH2F* hcentcnteta_tmp = (TH2F*)f[i]->Get("hcentcnteta");
            TH2F* hcentcntpt_tmp = (TH2F*)f[i]->Get("hcentcntpt");
            TH2F* hcentntrk_tmp = (TH2F*)f[i]->Get("hcentntrk");
        */
            TH2F* hcentcnteta_tmp = (TH2F*)f[j*153+i]->Get("hcentbbct0fracsouth_4");
            TH2F* hcentcnteta_tmp1 = (TH2F*)f[j*153+i]->Get("hcentbbcmixt0fracsouth_4");
            TH2F* hcentcnteta_tmp11 = (TH2F*)f[j*153+i]->Get("hcentbbct0fracsouth_7");
            TH2F* hcentcnteta_tmp12 = (TH2F*)f[j*153+i]->Get("hcentbbcmixt0fracsouth_7");
            TH2F* hcentcnteta_tmp21 = (TH2F*)f[j*153+i]->Get("hcentbbct0fracsouth_0");
            TH2F* hcentcnteta_tmp22 = (TH2F*)f[j*153+i]->Get("hcentbbcmixt0fracsouth_0");
            TH2F* hcentcntpt_tmp = (TH2F*)f[j*153+i]->Get("hcent1cent2_0");
            TH2F* hcentntrk_tmp = (TH2F*)f[j*153+i]->Get("hcent1cent2_4");
            TH2F* hcentntrk_tmp1 = (TH2F*)f[j*153+i]->Get("hcent1cent2_9");
       // }
       // catch( const std::exception& e){
       //     cout<<e.what()<<endl;
       //     continue;
       // }
        if(i==0){
            hcentcnteta = (TH2F*)hcentcnteta_tmp->Clone("hcentcnteta_add");
            hcentcnteta1 = (TH2F*)hcentcnteta_tmp1->Clone("hcentcnteta_add1");
            hcentcnteta11 = (TH2F*)hcentcnteta_tmp11->Clone("hcentcnteta_add11");
            hcentcnteta12 = (TH2F*)hcentcnteta_tmp12->Clone("hcentcnteta_add12");
            hcentcnteta21 = (TH2F*)hcentcnteta_tmp11->Clone("hcentcnteta_add21");
            hcentcnteta22 = (TH2F*)hcentcnteta_tmp12->Clone("hcentcnteta_add22");
            hcentcntpt = (TH2F*)hcentcntpt_tmp->Clone("hcentcntpt_add");
            hcentntrk = (TH2F*)hcentntrk_tmp->Clone("hcentntrk_add");
            hcentntrk1 = (TH2F*)hcentntrk_tmp1->Clone("hcentntrk_add1");
        }
        else{
            hcentcnteta->Add(hcentcnteta_tmp);
            hcentcnteta1->Add(hcentcnteta_tmp1);
            hcentcnteta11->Add(hcentcnteta_tmp11);
            hcentcnteta12->Add(hcentcnteta_tmp12);
            hcentcnteta21->Add(hcentcnteta_tmp21);
            hcentcnteta22->Add(hcentcnteta_tmp22);
            hcentcntpt->Add(hcentcntpt_tmp);
            hcentntrk->Add(hcentntrk_tmp);
            hcentntrk1->Add(hcentntrk_tmp1);
        }
       // f[i]->Close();
    }
    TFile *fout = new TFile(Form("merged_AnapAumbst%d.root",j),"Update");
    fout->cd();
    hcentcnteta->Write("hcentbbct0fracsouth_4",TObject::kOverwrite);
    hcentcnteta1->Write("hcentbbcmixt0fracsouth_4",TObject::kOverwrite);
    hcentcnteta11->Write("hcentbbct0fracsouth_7",TObject::kOverwrite);
    hcentcnteta12->Write("hcentbbcmixt0fracsouth_7",TObject::kOverwrite);
    hcentcnteta21->Write("hcentbbct0fracsouth_0",TObject::kOverwrite);
    hcentcnteta22->Write("hcentbbcmixt0fracsouth_0",TObject::kOverwrite);
    hcentcntpt->Write("hcent1cent2_0",TObject::kOverwrite);
    hcentntrk->Write("hcent1cent2_4",TObject::kOverwrite);
    hcentntrk1->Write("hcent1cent2_9",TObject::kOverwrite);
    fout->Close();
    }
}

