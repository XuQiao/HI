#include "par.h"
void dofirst(){
	gROOT->ProcessLine(".L Bessfun.C+");
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	TString name;
	for(int i=start;i<end;i++){
		if(isSum)	name=Form("/scratch/xuq7/flow/pbsjoboutput/Bessfun/pPbDataV200m120/AnaV_Sum_%d.root",i);
		else	        name=Form("/scratch/xuq7/flow/pbsjoboutput/Bessfun/pPbDataV200m120/AnaV_Prod_%d.root",i);
		remove(name.Data());
		Bessfun *l = new Bessfun(Form("%s/vndata_50k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		if(isSum)	l->calcV(1);
		else l->calcV(0);
		l->endJobV(name);
                delete l;
	}
}
