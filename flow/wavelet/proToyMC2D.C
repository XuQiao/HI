void proToyMC2D(){
    double eta1 = 2 ;
    double eta2 = -3;
    double phi1 = -2 ;
    double phi2 = 2;
    double A = 40;
    double sigma = 0.4;
    Double_t eta, phi;
    TF2 *f= new TF2("f","1./TMath::Sqrt((x-[0])*(x-[0])+(y-[1])*(y-[1]))+[2]*TMath::Exp(-0.5*((x-[3])*(x-[3])+(y-[4])*(y-[4]))/([5]*[5]))",-5,5,-4,4);
    f->SetNpx(3000);
    f->SetNpy(3000);
    f->SetParameters(eta1,eta2,A,phi1,phi2,sigma);
    TTree* tree = new TTree("tree","tree with #eta and #phi");
    TFile *fout = new TFile("pSTEG.root","Recreate");
    tree->Branch("eta",&eta,"eta/D");
    tree->Branch("phi",&phi,"phi/D");
    const int N = 100000;
    for(int imult = 0 ;imult<N; imult++){
      if(imult%1000==0)  cout<<"processing "<<imult<<" particles"<<endl;
    f->GetRandom2(eta,phi);
    h->Fill(eta,phi);
    tree->Fill();
    }
 /*   h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle("#phi");
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->Draw("surf1");
    */
    fout->cd();
    tree->Write();
}
