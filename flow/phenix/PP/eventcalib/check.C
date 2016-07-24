void check(){
    /*
     TString   name = "../filelistMWGfvtx.dat";
     TString   namec = "../filelistMWGfvtxcalib.dat";
     TString   name1 = "../filelistMWGfvtx1.dat";
     TString   name2 = "../filelistMWGfvtx2.dat";
     TString   name1c = "../filelistMWGfvtx1calib.dat";
     TString   name2c = "../filelistMWGfvtx2calib.dat";
    
    */
     TString   name = "../filelistMWGmb.dat";
     TString   namec = "../filelistMWGmbcalib.dat";
     TString   name1 = "../filelistMWGmb1.dat";
     TString   name2 = "../filelistMWGmb2.dat";
     TString   name1c = "../filelistMWGmb1calib.dat";
     TString   name2c = "../filelistMWGmb2calib.dat";
    
   //  ofstream fout("checkfvtx.out");
     //ofstream fout("checkfvtx1.out");
    ofstream fout("checkmb.out");
       // for(int i=0;i<1376;i++){
       float bbcv;
       float bbcvc;
       float bbcv1;
       float bbcv1c;
       float bbcv2;
       int nfvtxtrack1;
       int nfvtxtrack1c;
       int nfvtxtrack2;
       int nfvtxtrack2c;
        for(int i=0;i<754;i++){
            if(i!=693 && i!=44)continue;
    //    for(int i=256;i<258;i++){
          //  if(i!=170 && i!=171 && i!=174 && i!=178 && i!=257 && i!=300 && i!=460 && i!=491 && i!=535 && i!=621)continue; 
        //    if(i!=43)continue; 
            TString ifile(readline(Form("%s",name.Data()),i));
            TString ifilec(readline(Form("%s",namec.Data()),i));
            TString ifile1(readline(Form("%s",name1.Data()),i));
            TString ifile2(readline(Form("%s",name2.Data()),i));
            TString ifile1c(readline(Form("%s",name1c.Data()),i));
            TString ifile2c(readline(Form("%s",name2c.Data()),i));
            TFile *f = TFile::Open(ifile);
            TFile *fc = TFile::Open(ifilec);
            TFile *f1 = TFile::Open(ifile1);
            TFile *f2 = TFile::Open(ifile2);
            TFile *f1c = TFile::Open(ifile1c);
            TFile *f2c = TFile::Open(ifile2c);
            TTree *t = (TTree*)f->Get("tree");
            TTree *tc = (TTree*)fc->Get("tree");
            TTree *t1 = (TTree*)f1->Get("tree");
            TTree *t2 = (TTree*)f2->Get("tree");
            TTree *t1c = (TTree*)f1c->Get("tree");
            TTree *t2c = (TTree*)f2c->Get("tree");
            t->SetBranchAddress("bbcv",&bbcv);
            tc->SetBranchAddress("bbcv",&bbcvc);
            t1->SetBranchAddress("bbcv",&bbcv1);
            t1c->SetBranchAddress("bbcv",&bbcv1c);
            t2->SetBranchAddress("bbcv",&bbcv2);
            t1->SetBranchAddress("nfvtxtrack",&nfvtxtrack1);
            t1c->SetBranchAddress("nfvtxtrack",&nfvtxtrack1c);
            t2->SetBranchAddress("nfvtxtrack",&nfvtxtrack2);
            t2c->SetBranchAddress("nfvtxtrack",&nfvtxtrack2c);
            cout<<t->GetEntries()<<"\t"<<tc->GetEntries()<<"\t"<<t1->GetEntries()<<"\t"<<t1c->GetEntries()<<"\t"<<t2->GetEntries()<<"\t"<<t2c->GetEntries()<<"\t";
            fout<<t->GetEntries()<<"\t"<<tc->GetEntries()<<"\t"<<t1->GetEntries()<<"\t"<<t1c->GetEntries()<<"\t"<<t2->GetEntries()<<"\t"<<t2c->GetEntries()<<"\t";
           // cout<<ifile<<"\t";//<<t2->GetEntries();
            t->GetEntry(t->GetEntries()-100);
            tc->GetEntry(tc->GetEntries()-100);
            t1->GetEntry(t1->GetEntries()-100);
            t1c->GetEntry(t1c->GetEntries()-100);
            t2->GetEntry(t1->GetEntries()-100);
            t2c->GetEntry(t1c->GetEntries()-100);
//            t->GetEntry(0);
//            t1->GetEntry(0);
  //          cout<<bbcv<<"\t"<<bbcvc<<"\t"<<bbcv1<<"\t"<<bbcv1c<<"\t"<<bbcv2<<endl;
  //          cout<<nfvtxtrack1<<"\t"<<nfvtxtrack1c<<"\t"<<nfvtxtrack2<<"\t"<<nfvtxtrack2c<<endl;
 //           fout<<t->GetEntries()<<"\t"<<tc->GetEntries()<<"\t"<<t1->GetEntries()<<"\t"<<t1c->GetEntries()<<"\t"<<t2->GetEntries()<<"\t"<<t2c->GetEntries();
            if(tc->GetEntries()!=t2c->GetEntries() || nfvtxtrack1c != nfvtxtrack2c) 
          //  if(tc->GetEntries()!=t1c->GetEntries() || bbcvc!=bbcv1c) 
            {
                cout<<"BAD";
                fout<<"BAD";
            }
            if(tc->GetEntries()!=0)
            if(1.0*t->GetEntries() / tc->GetEntries()>1.01)
            {
                cout<<"Number BAD!";
                fout<<"Number BAD!";
            }
            if(t1c->GetEntries()!=0)
            if(1.0*t1->GetEntries() / t1c->GetEntries()>1.01)
            {
                cout<<"Number1 BAD!";
                fout<<"Number1 BAD!";
            }
            if(t2c->GetEntries()!=0)
            if(1.0*t2->GetEntries() / t2c->GetEntries()>1.01)
            {
                cout<<"Number2 BAD!";
                fout<<"Number2 BAD!";
            }
//            if(bbcv!=bbcv1) cout<<"BAD!";
            cout<<endl;
            fout<<endl;
            f->Close();
            fc->Close();
            f1->Close();
            f1c->Close();
            f2->Close();
            f2c->Close();
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
