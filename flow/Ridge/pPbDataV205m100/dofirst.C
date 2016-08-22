#include "par.h"
void dofirst(){
	gROOT->ProcessLine(".L ridgeAn.C+");
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	TString name;
	for(int i=start;i<end;i++){
		name=Form("%s/Ana_%d.root",outdir.Data(),i);
		remove(name.Data());
		ridge *l = new ridge(Form("%s/vndata_50k_%d.root",indir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
                l->beginJob();
		l->calcS();
		l->calcB();
		l->endJob(name);
                delete l;
	}
}
