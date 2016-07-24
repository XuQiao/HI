void calihisto(){
    int start=atoi(getenv("BEGIN"));
    int end=atoi(getenv("END"));
    string trig=getenv("TRIG");
    TString trigtype(trig);
    TString name,namec,name1,name1c;
    if(trigtype.Contains("fvtx"))  {
        name = "../filelistMWGfvtx.dat";
        namec = "../filelistMWGfvtxcalib.dat";
       // name1 = "../filelistMWGfvtx1.dat";
        name1 = "../filelistMWGfvtx2.dat";
       // name1c = "../filelistMWGfvtx1calib.dat";
        name1c = "../filelistMWGfvtx2calib.dat";
    }
    else if(trigtype.Contains("mb")){
        name = "../filelistMWGmb1calib.dat";
        namec = "../filelistMWGmbcalib.dat";
       // name1 = "../filelistMWGmb1.dat";
        name1 = "../filelistMWGmb2.dat";
       // name1c = "../filelistMWGmb1calib.dat";
        name1c = "../filelistMWGmb2calib.dat";
    }
    else exit();
    for(int i=start;i<end;i++){
        cout << i << endl;
//        m *pl = new m(readline(Form("%s",name.Data()),i),readline(Form("%s",namec.Data()),i), readline(Form("%s",name1.Data()),i));
        //mm *pl = new mm(readline(Form("%s",name.Data()),i),readline(Form("%s",name1.Data()),i), readline(Form("%s",namec.Data()),i), readline(Form("%s",name1c.Data()),i));
        mm *pl = new mm(readline(Form("%s",name.Data()),i),readline(Form("%s",name1.Data()),i), readline(Form("%s",name1c.Data()),i));
       // m1 *pl = new m1(readline(Form("%s",name.Data()),i),readline(Form("%s",name2.Data()),i));
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
