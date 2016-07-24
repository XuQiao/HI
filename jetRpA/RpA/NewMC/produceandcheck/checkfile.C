void checkfile(){
    int n;
    for(int i=96;i<98;i++){
    //TFile *f1 = new TFile(Form("/cms/store/user/qixu/jetRpA/skimTree/DATAPPbakPu3PFskimfile%d_%d.root",i*50,(i+1)*50),"ReadOnly");
    //TFile *f2 = new TFile(Form("/cms/store/user/ymao/pA5TEV/PAHighPt/HIRun2013-PromptReco-v1/skimfile/ForwardPb/DATAPPbakPu3PFskimFullInfofile%d_%d.root",i*50,(i+1)*50),"ReadOnly");
    TFile *f2 = new TFile(Form("/cms/store/user/ymao/pA5TEV/PAHighPt/HIRun2013-PromptReco-v1/skimfile/DATAPbPakPu3PFskimFullInfofile%d_%d.root",i*50,(i+1)*50),"ReadOnly");
    //TTree *nt1 = (TTree*)f1->Get("nt");
    TTree *nt2 = (TTree*)f2->Get("nt");
    //int nevent1 = nt1->GetEntries();
    int nevent2 = nt2->GetEntries();
    cout<<i<<"\t"<<nevent2<<endl;
        //"\t"<<nevent2<<endl;
//    n+=nevents;
//    nt1= NULL;
    nt2= NULL;
 //   f1->Close();
    f2->Close();
    }
    //cout<<n<<endl;
}
