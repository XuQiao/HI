#include "par.h"
void do(){
	gROOT->ProcessLine(".L ../EP.C+");
	for(int i=0;i<nFileAll;i++){
		remove(Form("./jobout/AnaEP_%d.root",i));
		EP *l = new EP(Form("%s/vndata_10k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		l->calcEP();
		l->endJobEP(Form("./jobout/AnaEP_%d.root",i));
	}
}
