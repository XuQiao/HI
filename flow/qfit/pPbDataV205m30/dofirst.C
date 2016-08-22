#include "par.h"
void dofirst(){
	gROOT->ProcessLine(".L qfit.C+");
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	TString name;
	for(int i=start;i<end;i++){
		name=Form("/scratch/xuq7/flow/pbsjoboutput/qfit/pPbDataV205m30/AnaV_Sum_%d.root",i);
		remove(name.Data());
		qfit *l = new qfit(Form("%s/vndata_50k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		l->calcV();
		l->endJobV(name);
                delete l;
	}
}
