#include "par.h"
void dofirst(){
	remove("GRe.txt");
	remove("GIm.txt");
	remove("Q2.txt");
	remove("Qx1.txt");
	remove("Qy1.txt");
	remove("mult.txt");
	gROOT->ProcessLine(".L ../LYZ.C+");
	for(int i=0;i<nFileAll;i++){
		remove(Form("./jobout/AnaV_Sum_%d.root",i));
		LYZ *l = new LYZ(Form("%s/vndata_10k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		l->calcV(1);
		l->endJobV(Form("./jobout/AnaV_Sum_%d.root",i));
	}
}
