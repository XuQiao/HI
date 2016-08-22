#include "par.h"
void dofirst(){
	remove("GRe.txt");
	remove("GIm.txt");
	remove("Q2.txt");
	remove("Qx1.txt");
	remove("Qy1.txt");
	remove("mult.txt");
	int start=atoi(getenv("START"));
	int end=atoi(getenv("END"));
	TString name;
	string SumorProd=getenv("SUMORPROD");
	for(int i=start;i<end;i++){
		if(SumorProd=="Sum")	name=Form("./jobout/AnaV_Sum_%d.root",i);
		else	name=Form("./jobout/AnaV_Prod_%d.root",i);
		remove(name.Data());
		LYZ *l = new LYZ(Form("%s/vndata_50k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		if(SumorProd=="Sum")	l->calcV(1);
		else l->calcV(0);
		l->endJobV(name);
	}
}

void dosecond(){
        gSystem->AddIncludePath("-D__USE_XOPEN2K8");
        remove("dDRe.txt");
        remove("dDIm.txt");
        remove("dNRe.txt");
        remove("dNIm.txt");
        remove("mult_.txt");
	string SumorProd=getenv("SUMORPROD");
        int start=atoi(getenv("START"));
        int end=atoi(getenv("END"));
        TString name;
        for(int i=start;i<end;i++){
                if(SumorProd=="Sum")       name=Form("./jobout/Anav_Prod_%d.root",i);
                else            name=Form("./jobout/Anav_Prod2_%d.root",i);
                remove(name.Data());
                LYZ *l = new LYZ(Form("%s/vndata_50k_%d.root",dir.Data(),i));
                cout<<"start "<<i<<" th job"<<endl;
                l->beginJob();
                if(SumorProd=="Sum")       l->calcv("mergedV_Sum.root",0);
                else    l->calcv("mergedV_Prod.root",0);
                l->endJobv(name);
        }
}

