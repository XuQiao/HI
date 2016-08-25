void Writehisto(){
    int start=atoi(getenv("BEGIN"));
    int end=atoi(getenv("END"));
    string trig=getenv("TRIG");
    TString trigtype(trig);
    TString name,name1,name2;
//    if(trigtype.Contains("fvtx"))  name = "../filelistfvtx.dat";
//    else if(trigtype.Contains("mb"))  name = "../filelistmb.dat";
    if(trigtype.Contains("fvtx")) { 
        name = "../filelistMWGfvtxcalib.dat";
        name1 = "../filelistMWGfvtx1calib.dat";
        name2 = "../filelistMWGfvtx2calib.dat";
    }
    else if(trigtype.Contains("mb")){
        name = "../filelistMWGmbcalib.dat";
        name1 = "../filelistMWGmb1calib.dat";
        name2 = "../filelistMWGmb2calib.dat";
    }
    else exit();
    for(int i=start;i<end;i++){
        cout << i << endl;
        //Perform *pl = new Perform(readline(Form("%s",name.Data()),i), Form("output/Ana%s_%d.root",trigtype.Data(),i));
        //PerformMWG *pl = new PerformMWG(readline(Form("%s",name.Data()),i), readline(Form("%s",name1.Data()),i), readline(Form("%s",name2.Data()),i), Form("/store/user/qixu/.phenix/PP/Perform/output/AnaMWG%s_%d.root",trigtype.Data(),i));
        PerformMWG *pl = new PerformMWG(readline(Form("%s",name.Data()),i), readline(Form("%s",name1.Data()),i), readline(Form("%s",name2.Data()),i), Form("output/AnaMWG%s_%d.root",trigtype.Data(),i));
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
