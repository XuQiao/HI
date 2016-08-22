#include "par.h"
void dofirst(){
        gInterpreter->AddIncludePath("/home/xuq7/HI/flow/QCumulant/correlations-0.9_test");
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	gROOT->ProcessLine(".L Analyze.C+");
	for(int i=start;i<end;i++){
	    cout<<"start "<<i<<" th job"<<endl;
    	    gROOT->ProcessLine(Form("Analyze(\"%s\",%d,%d,%d)",method.Data(),maxN,useWeight,i));
	}
}
