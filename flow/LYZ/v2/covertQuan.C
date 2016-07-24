void covertQuan(){
    TFile *fin = TFile::Open("pPb.root");
//    TFile *fout = new TFile("pPb_xq.root","Recreate");
    int mult[4] = {120,150,185,220};
    int n[4] = {2,4,6,8};
    for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
        TGraphErrors *g = (TGraphErrors*)fin->Get(Form("v2Pt_%d_%d",n[j],mult[i]));
        if(n[j]==4){
        cout<<mult[i]<<endl;
        for(int ibin=0;ibin<g->GetN();ibin++){
            cout<<g->GetX()[ibin]<<",";
        }
        cout<<endl;
        for(int ibin=0;ibin<g->GetN();ibin++){
            cout<<g->GetY()[ibin]<<",";
        }
        cout<<endl;
        for(int ibin=0;ibin<g->GetN();ibin++){
            cout<<g->GetEY()[ibin]<<",";
        }
        cout<<endl;
        }
       // fout->cd();
       // g->Write("AP");
       // g1->Write();
    }
    }
}

