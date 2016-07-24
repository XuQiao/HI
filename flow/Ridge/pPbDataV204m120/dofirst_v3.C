#include "par.h"
void dofirst_v3(){
	gROOT->ProcessLine(".L ridgeAn_v3.C+");
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	TString name;
	for(int i=start;i<end;i++){
		name=Form("%s/Anav3_%d.root",outdir.Data(),i);
		remove(name.Data());
		ridge *l = new ridge(Form("%s/STEGv2onlyFixM120_50k.root",indir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
                l->beginJob();
		l->calc();
		l->endJob(name);
                delete l;
	}
}
