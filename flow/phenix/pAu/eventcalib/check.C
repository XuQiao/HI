void check(){
     //  TString name = "../filelistfvtx.dat";
     //  TString name1 = "../filelistfvtx1.dat";
     //  TString name2 = "../filelistfvtx2.dat";
     TString   name = "../filelistmb.dat";
     TString   name1 = "../filelistmb1.dat";
     TString   name2 = "../filelistmb2.dat";
     TString   name3 = "../filelistmb3.dat";
       // for(int i=0;i<1376;i++){
       int nbbc;
       int nbbc1;
       int nbbc2;
       int nbbc3;
       float bbcv;
       float bbcv1;
       float bbcv2;
       float bbcv3;
        for(int i=0;i<1376;i++){
            TString ifile(readline(Form("%s",name.Data()),i));
            TString ifile1(readline(Form("%s",name1.Data()),i));
            TString ifile2(readline(Form("%s",name2.Data()),i));
            TString ifile3(readline(Form("%s",name3.Data()),i));
            TFile *f = TFile::Open(ifile);
            TFile *f1 = TFile::Open(ifile1);
            TFile *f2 = TFile::Open(ifile2);
            TFile *f3 = TFile::Open(ifile3);
            TTree *t = (TTree*)f->Get("tree");
            TTree *t1 = (TTree*)f1->Get("tree");
            TTree *t2 = (TTree*)f2->Get("tree");
            TTree *t3 = (TTree*)f3->Get("tree");
            t->SetBranchAddress("nbbc",&nbbc);
            t->SetBranchAddress("bbcv",&bbcv);
            t1->SetBranchAddress("nbbc",&nbbc1);
            t2->SetBranchAddress("nbbc",&nbbc2);
            t3->SetBranchAddress("bbcv",&bbcv3);
            cout<<t->GetEntries()<<"\t"<<t1->GetEntries()<<"\t"<<t2->GetEntries()<<"\t"<<t3->GetEntries()<<"\t";
            t->GetEntry(t->GetEntries()-134);
            t1->GetEntry(t1->GetEntries()-134);
            t2->GetEntry(t2->GetEntries()-134);
            t3->GetEntry(t3->GetEntries()-134);
            if(t->GetEntries()!=t2->GetEntries() || bbcv!=bbcv3) cout<<"BAD!";
            cout<<endl;
            f->Close();
            f1->Close();
            f2->Close();
            f3->Close();
        }
    }

std::string readline(char* name, int iline){
        std::ifstream backstory(name);
        std::string line;
        if (backstory.is_open())
                if(backstory.good()){
                        for(int i = 0; i < iline+1; ++i)
                           getline(backstory, line);
                }
        return line;
}
