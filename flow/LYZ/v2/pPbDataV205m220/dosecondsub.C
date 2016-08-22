#include "par.h"
void dosecondsub(){
	//gSystem->AddIncludePath("-D__USE_XOPEN2K8");
	remove("dDRe.txt");
	remove("dDIm.txt");
	remove("dNRe.txt");
	remove("dNIm.txt");
	remove("mult_.txt");
	gROOT->ProcessLine(".L ../LYZ.C+");
        int start=atoi(getenv("START"));
        int end=atoi(getenv("END"));
        int subsample = (int)(nFileAll/nsamples);
        TString name;
	for(int i=start;i<end;i++){	
		if(isSum)	name=Form("/lio/lfs/cms/store/user/qixu/flow/pbsjoboutput/pPbDataV205m220/Anav_Prod_sub_%d.root",i);
		else		name=Form("/lio/lfs/cms/store/user/qixu/flow/pbsjoboutput/pPbDataV205m220/Anav_Prod2_sub_%d.root",i);
		remove(name.Data());
		LYZ *l = new LYZ(Form("%s/vndata_50k_%d.root",dir.Data(),i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
        	for(int isample=0;isample<nsamples;isample++){
			if(i>=isample*subsample&&i<(isample+1)*subsample){
				if(isSum)	l->calcv("mergedV_Sum_sub.root",0,isample);
				else	l->calcv("mergedV_Prod_sub.root",0,isample);
				l->endJobv(name);
			}
		}
	}
}
