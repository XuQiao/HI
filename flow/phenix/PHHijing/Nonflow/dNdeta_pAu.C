#define MAXP 100000
void dNdeta_pAu(int i){
    //TFile *f = TFile::Open(Form("/phenix/plhf3/xuq/MCgen/Hijing/phhijing_pAu_%d.root",i));
    TFile *f = TFile::Open(Form("/phenix/plhf3/xuq2/MCgen/phhijing_dAu_%d.root",i));
    TTree *tree = (TTree*)f->Get("T");
    TBranchElement *array = (TBranchElement*)T->GetBranch("part_array");
    tree->SetMakeClass(1);
    TH1::SetDefaultSumw2(kTRUE);
    TH1F* hKF = new TH1F("hKF","particle id",10000,-5000,5000);
    TH1F* hn = new TH1F("hn","Ntotal",1000,0,1000);
    TH1F* hnch = new TH1F("hnch","Ncharge",500,0,500);
    TH1F* heta = new TH1F("heta","particle eta distr",120,-6,6);
    TH1F* hy = new TH1F("hy","particle y distr",120,-6,6);
    TH1F* hy1 = new TH1F("hy1","particle y distr in |y|<1",150,-1.5,1.5);
    TH1F* hphi = new TH1F("hphi","",200,0,8);
    TH1F* hMass= new TH1F("hMass","mass distribution",500,0,5);
    TH2F* hpteta = new TH2F("hpteta","",100,0,10,120,-6,6);
    TH2F* hpty = new TH2F("hpty","",100,0,10,120,-6,6);
    
    float px[MAXP];
    float py[MAXP];
    float pz[MAXP];
    float E[MAXP];
    int KS[MAXP];
    int KF[MAXP];
    float M[MAXP];
    int n;
    int nPion = 0;
    int nPrim = 0;
    int nSecond = 0;
    int nPionC = 0;
    int nKaonC = 0;
    int nout=0;
    tree->SetBranchAddress("part_array.fKS",&KS);
    tree->SetBranchAddress("part_array.fKF",&KF);
    tree->SetBranchAddress("part_array.fPx",&px);
    tree->SetBranchAddress("part_array.fPy",&py);
    tree->SetBranchAddress("part_array.fPz",&pz);
    tree->SetBranchAddress("part_array.fEnergy",&E);
    tree->SetBranchAddress("part_array.fMass",&M);
    array->SetAddress(&n);
    for(int ievent = 0;ievent < tree->GetEntries(); ievent++){
        int nch=0;
        tree->GetEntry(ievent);
 //       array->GetEntry(ievent);
 //       int n = array->GetNdata();
        if(n>MAXP)  {cout<<n<<"!"<<endl;  continue;}
        if(ievent%10000==0) cout<<ievent<<endl;
        for(int iparticle = 0;iparticle < n; iparticle++){
     //     if(KS[iparticle] < 1 || KS[iparticle] > 10) continue;
     //     if(iparticle>10000) cout<<KF[iparticle]<<" ";
            hKF->Fill(KF[iparticle]);
            if(KF[iparticle] == 111) nPion++;
            if(fabs(KF[iparticle]) == 211) nPionC++;
            if(fabs(KF[iparticle]) == 321) nKaonC++;
            if(fabs(KF[iparticle]) > 5000) nout++;
         //   if(KF[iparticle]!=313 && KF[iparticle]!=421 && KF[iparticle]!=311 && KF[iparticle]!d221 && KF[iparticle]!=111 && KF[iparticle]!=223 && KF[iparticle]!=333 && KF[iparticle]!=443 && KF[iparticle]!=2112 && KF[iparticle]!=3122) nch++;
            if(!(fabs(KF[iparticle]) == 211 || fabs(KF[iparticle]) == 213 || fabs(KF[iparticle]) == 321 || fabs(KF[iparticle]) == 323 ||fabs(KF[iparticle]) == 2212)) continue;
            if(KS[iparticle]==1) nPrim++;
            if(KS[iparticle]==11) nSecond++;
            hMass->Fill(M[iparticle]);
            if(KS[iparticle]!=1) continue;
            float pt = sqrt(px[iparticle]**2+py[iparticle]**2);
            float theta = atan(pt/pz[iparticle]);
            if(theta<0) theta = TMath::Pi()+theta;
            float eta = -log(tan(theta/2));
            if(px[iparticle]>0){
                float phi = atan(py[iparticle]/px[iparticle]);
            if(py[iparticle]<0)
                float phi = 2*TMath::Pi()+atan(py[iparticle]/px[iparticle]);
            }
            else
                float phi = TMath::Pi()+atan(py[iparticle]/px[iparticle]);
            if((E[iparticle]+pz[iparticle])/(E[iparticle]-pz[iparticle])<=0) continue;
            float y = 0.5*log((E[iparticle]+pz[iparticle])/(E[iparticle]-pz[iparticle]));
            hy->Fill(y);
            heta->Fill(eta);
            if(fabs(y)<=1)
            hy1->Fill(y);

            hphi->Fill(phi);
            hpteta->Fill(pt,eta);
            hpty->Fill(pt,y);
        }
        hn->Fill(n);
        hnch->Fill(nch);
    }
    TFile *fout = new TFile(Form("output/outhisto_dAu_%d.root",i),"Recreate");
    fout->cd();
    hKF->Write();
    heta->Write();
    hy->Write();
    hy1->Write();
    hphi->Write();
    hpteta->Write();
    hpty->Write();
    hMass->Write();
    hn->Write();
    hnch->Write();
    fout->Close();
    cout<<nPion<<"\t"<<nPionC<<"\t"<<nKaonC<<"\t"<<nout<<"\t"<<nPrim<<"\t"<<nSecond<<endl;
}
