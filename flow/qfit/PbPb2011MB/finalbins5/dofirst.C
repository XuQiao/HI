#include "par.h"
void dofirst(){
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	TString name;
	for(int i=start;i<end;i++){
		name=Form("/scratch/xuq7/flow/pbsjoboutput/qfit/PbPb2011MB/finalbins5/AnaV_Sum_%d.root",i);
		remove(name.Data());
		qfit *l = new qfit(readline("./filelist.dat",i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		l->calcV();
		l->endJobV(name);
                delete l;
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
