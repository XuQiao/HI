#include "par.h"
void dosecond(){
	remove("dDRe.txt");
	remove("dDIm.txt");
	remove("dNRe.txt");
	remove("dNIm.txt");
	remove("mult_.txt");
	gROOT->ProcessLine(".L ../LYZ.C+");
	for(int i=0;i<nFileAll;i++){
		remove(Form("./jobout/Anav_Prod_%d.root",i));
		LYZ *l = new LYZ(Form("%s/vndata_10k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		l->calcv("mergedV_Sum.root",0);
		l->endJobv(Form("./jobout/Anav_Prod_%d.root",i));
	}
}
