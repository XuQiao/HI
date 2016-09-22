void Writehisto(){
    int i=atoi(getenv("I"));
    string trig="dAu62all";
    TString trigtype(trig);
    TString name, name1;
    name = "tree.lst";
    name1 = "tree1.lst";
        std::ifstream corrs("Run16dAu62GeV.lst");
        int index=0; int run=0;
        for(int irun=0;irun<i+1;irun++){
        corrs>>index>>run;
        }

        //EPAnaRun16alltree *pl = new EPAnaRun16alltree(readline(Form("%s",name.Data()),i),Form("/scratch/xuq7/phenix/Run16dAu/39GeV/EPAnaFull%s_%d.root",trigtype.Data(),i));
        EPAnaRun16alltree *pl = new EPAnaRun16alltree(readline(Form("%s",name.Data()),i), readline(Form("%s",name1.Data()),i), Form("output/EPAnaFull%s_%d.root",trigtype.Data(),i));
        pl->Init();
        pl->Inittree();
        pl->SetcalFlag(0);
        pl->process_event();
        pl->SetcalFlag(1);
        pl->process_event();
        pl->SetcalFlag(2);
        pl->process_event();
        pl->SetcalFlag(3);
        pl->process_event();
        pl->End();
        delete pl;
}

std::vector<TString> readline(char* name, int i){
        std::ifstream backstory(name);
        std::vector<TString> vecline;
        std::string line;
        std::ifstream corrs("Run16dAu62GeV.lst");
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
