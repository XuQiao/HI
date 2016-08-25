#include "par.h"
void dofirst(){
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	TString name;
	for(int i=start;i<end;i++){
		name=Form("/scratch/xuq7/flow/pbsjoboutput/qfit/pPbDataV205m50/AnaV_Sum_%d.root",i);
		remove(name.Data());
		qfit *l = new qfit(Form("%s/vndata_50k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
	        l->calcV();
		l->endJobV(name);
                delete l;
	}
}
