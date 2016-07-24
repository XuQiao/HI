#include "par.h"
void dofirst(){
	//gSystem->AddIncludePath("-D__USE_XOPEN2K8");
	remove("GRe.txt");
	remove("GIm.txt");
	remove("Q2.txt");
	remove("Qx1.txt");
	remove("Qy1.txt");
	remove("mult.txt");
	gROOT->ProcessLine(".L LYZ.C+");
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	TString name;
	for(int i=start;i<end;i++){
		if(isSum)	name=Form("/scratch/xuq7/flow/pbsjoboutput/pPbDataV205m220/AnaV_Sum_%d.root",i);
		else	        name=Form("/scratch/xuq7/flow/pbsjoboutput/pPbDataV205m220/AnaV_Prod_%d.root",i);
		remove(name.Data());
		LYZ *l = new LYZ(Form("%s/vndata_50k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		if(isSum)	l->calcV(1);
		else l->calcV(0);
		l->endJobV(name);
	}
}
