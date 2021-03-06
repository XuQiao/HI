void Writehisto(){
    int start=atoi(getenv("BEGIN"));
    int end=atoi(getenv("END"));
    string trig=getenv("TRIG");
    TString trigtype(trig);
    TString name,name1,name2;
    if(trigtype.Contains("fvtx"))  {
        name = "../filelistfvtx.dat";
        name1 = "../filelistfvtx2.dat"; //event calibrated
        name2 = "../filelistfvtx3.dat";
    }
    else if(trigtype.Contains("mb")){
        name = "../filelistmb.dat";
        name1 = "../filelistmb2.dat"; //event calibrated
        name2 = "../filelistmb3.dat";
    }
    else exit();
    for(int i=start;i<end;i++){
        cout << i << endl;
        Ridgecntbbc *pl = new Ridgecntbbc(readline(Form("%s",name.Data()),i), readline(Form("%s",name1.Data()),i),readline(Form("%s",name2.Data()),i),Form("/scratch/xuq7/phenix/pAu/Ridgecntbbc/Ana%s_%d.root",trigtype.Data(),i));
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
