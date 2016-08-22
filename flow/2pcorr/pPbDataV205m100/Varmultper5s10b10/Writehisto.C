void Writehisto(){
    int start=atoi(getenv("BEGIN"));
    int end=atoi(getenv("END"));
    TString name = "filelist.dat";
    for(int i=start;i<end;i++){
        cout << "ifile:" <<i << endl;
        Ridgecntbbc *pl = new Ridgecntbbc(readline(Form("%s",name.Data()),i), Form("/scratch/xuq7/flow/2pcorr/pPbDataV205m100/Varmultper5s10b10/pPbAnacntbbc_%d.root",i));
        //Ridgecntcnt *pl = new Ridgecntcnt(readline(Form("%s",name.Data()),i), Form("AnaPPcntcnt%d.root",i));
        //Ridgefull *pl = new Ridgefull(readline(Form("%s",name.Data()),i), Form("AnaPPfull%d.root",i));
        pl->Init();
        pl->process_event();
        pl->End();
        delete pl;
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
