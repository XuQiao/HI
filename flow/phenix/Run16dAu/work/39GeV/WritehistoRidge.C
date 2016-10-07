void Writehisto(){
    int i=atoi(getenv("I"));
    string trig="dAu39all";
    TString trigtype(trig);
    TString name, name1;
    name = "tree.lst";
    name1 = "tree1.lst";
        std::ifstream corrs("Run16dAu39GeV.lst");
        int index=0; int run=0;
        for(int irun=0;irun<i+1;irun++){
        corrs>>index>>run;
        }

        //RidgedAuRun16 *pl = new RidgedAuRun16(readline(Form("%s",name.Data()),i),Form("/scratch/xuq7/phenix/Run16dAu/39GeV/EPAnaFull%s_%d.root",trigtype.Data(),i));
        RidgedAuRun16 *pl = new RidgedAuRun16(readline(Form("%s",name.Data()),i), readline(Form("%s",name1.Data()),i), Form("/store/user/qixu/flow/Run16dAu/39GeV/Ridge%s_%d.root",trigtype.Data(),i));
        pl->Init();
        pl->Inittree();
        pl->process_event();
        pl->End();
        delete pl;
}

std::vector<TString> readline(char* name, int i){
        std::ifstream backstory(name);
        std::vector<TString> vecline;
        std::string line;
        std::ifstream corrs("Run16dAu39GeV.lst");
        int index=0; int run=0;
        for(int irun=0;irun<i+1;irun++){
        corrs>>index>>run;
        }
        if (backstory.is_open())
                if(backstory.good()){
                    while(!backstory.eof()){
                            getline(backstory, line);
                            if(i>=0){
                            TString tline(line);
                            if(!tline.Contains(Form("%d",run)))continue;
                            vecline.push_back(tline);
                            }
                        }
                }
        return vecline;
}
