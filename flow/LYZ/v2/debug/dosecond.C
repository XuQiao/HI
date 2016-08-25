#include "par.h"
void dosecond(){
	//gSystem->AddIncludePath("-D__USE_XOPEN2K8");
	remove("dDRe.txt");
	remove("dDIm.txt");
	remove("dNRe.txt");
	remove("dNIm.txt");
	remove("mult_.txt");
	gROOT->ProcessLine(".L LYZ.C+");
        int start=0;
        int end=1;
        TString name;
	for(int i=start;i<end;i++){	
		if(isSum)	name=Form("Anav_Prod_%d.root",i);
		else		name=Form("Anav_Prod2_%d.root",i);
		remove(name.Data());
		LYZ *l = new LYZ(Form("%s/vndata_50k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		if(isSum)	l->calcv("mergedV_Sum.root",0,-1);
		else	l->calcv("mergedV_Prod.root",0,-1);
		l->endJobv(name);
	}
}
