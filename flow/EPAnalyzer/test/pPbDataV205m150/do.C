#include "par.h"
void do(){
	gROOT->ProcessLine(".L ../EP.C+");
	for(int i=0;i<nFileAll;i++){
		remove(Form("/scratch/xuq7/flow/EPAnalyzer/pPbDataV205m150/AnaEP_%d.root",i));
		EP *l = new EP(Form("%s/vndata_50k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		l->calcEP();
		l->endJobEP(Form("/scratch/xuq7/flow/EPAnalyzer/pPbDataV205m150/AnaEP_%d.root",i));
	}
}
