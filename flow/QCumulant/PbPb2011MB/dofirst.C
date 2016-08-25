#include "par.h"
void dofirst(){
        gInterpreter->AddIncludePath("/home/xuq7/HI/flow/QCumulant/correlations-0.9_test");
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	gROOT->ProcessLine(".L Analyze.C+");
	for(int i=start;i<end;i++){
	    cout<<"start "<<i<<" th job"<<endl;
            string name = readline("filelist.dat",i);
    	    gROOT->ProcessLine(Form("Analyze(\"%s\",%d,0,\"%s\",%d)",method.Data(),maxN,name.c_str(),i));
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
