void test(){
int nfile=0;
    std::string infile = "PAHighPtPurdueForest_4th.txt";
//    std::string infile = "test.txt";
std::ifstream instr(infile.c_str(),std::ifstream::in);
    std::string fname;
    std::string line;
        int iline;
    for (iline = 0; std::getline(instr, line); ++iline)
    ;
    int nline=iline;
cout<<nline<<endl;
    vector<string> listVector;
instr.clear();
instr.seekg(0, ios::beg);
    for(int ifile=0; ifile<5; ifile++){
      instr>>fname;
      listVector.push_back(fname);
}
TString inname=Form("%s",listVector.at(nfile).c_str());
	cout<<inname<<endl;
}           
