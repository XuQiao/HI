void checkEvent(){
	int i;
	for(i=0;i<462;i++){
		string oldname=readline("/home/xuq7/tmp/flow/CMSSW_5_3_8_HI_patch2/src/Analysis/LYZfirst/test/crab/tracknormcpt03to6/anaM185120/filelist.dat",i);
		TFile *fsold=TFile::Open(oldname.c_str());
		TFile *fs=TFile::Open(Form("jobout/AnaV_Sum_%d.root",i));
//		TFile *fp=TFile::Open(Form("jobout/AnaV_Prod_%d.root",i));
		TVectorD *neve1 = (TVectorD*)fs->Get("totmultall");
		TVectorD *neve1old = (TVectorD*)fsold->Get("demo/totmultall");
//		TVectorD *neve2 = (TVectorD*)fp->Get("Nevent");
		int nevent1 = (*neve1)[0];
		int nevent1old = (*neve1old)[0];
//		int nevent2 = (*neve2)[0];
//		if(nevent1!=nevent1old)
		cout<<i<<"\t"<<nevent1<<"\t"<<nevent1old<<endl;//"\t"<<nevent2<<endl;
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
