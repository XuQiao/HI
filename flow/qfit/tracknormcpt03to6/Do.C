#include "TSystem.h"
#include "TRint.h"
#include "TROOT.h"
#include <fstream>
#include <iostream>
#include "par.h"

void dofirst(){
	TString name;
	int start=atoi(getenv("BEGIN"));
	int end=atoi(getenv("END"));
	string dir=getenv("DIR");
	for(int i=start;i<end;i++){
		name=Form("/scratch/xuq7/flow/pbsjoboutput/qfit/tracknormcpt03to6/%s/AnaV_Sum_%d.root",dir.c_str(),i);
		remove(name.Data());
		qfit *l = new qfit(readline("./filelist.dat",i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		l->calcV();
		l->endJobV(name);
                delete l;
	}
}

void dosecond(){
        int start=atoi(getenv("BEGIN"));
        int end=atoi(getenv("END"));
	string dir=getenv("DIR");
        TString name;
        for(int i=start;i<end;i++){
                name=Form("/scratch/xuq7/flow/pbsjoboutput/qfit/tracknormcpt03to6/%s/Anav_Sum_%d.root",dir.c_str(),i);
                remove(name.Data());
                qfit *l = new qfit(readline("filelist.dat",i));
                cout<<"start "<<i<<" th job"<<endl;
                l->beginJob(1);
                l->calcv();
                l->endJobv(name);
        }
}

void dosecondsub(){
        int start=atoi(getenv("BEGIN"));
        int end=atoi(getenv("END"));
        string dir=getenv("DIR");
        TString name, infname;	
	TFile *inf;
   	infname="mergedV_Sum_sub.root";
        else    infname="mergedV_Prod_sub.root";
        for(int i=start;i<end;i++){
                name=Form("/scratch/xuq7/flow/pbsjoboutput/qfit/tracknormcpt03to6/%s/Anav_Sum_sub_%d.root",dir.c_str(),i);
                remove(name.Data());
                qfit *l = new qfit(readline("filelist.dat",i));
                cout<<"start "<<i<<" th job"<<endl;
                l->beginJob(1);
		for(int isample=0;isample<nsamples;isample++){
                        if(i>getfilenumber(infname,isample-1)&&i<=getfilenumber(infname,isample)){
			l->calcv(infname,0,isample);
		}
                l->endJobv(name);
        	}
	}
}

void dothird(){
        int start=atoi(getenv("BEGIN"));
        int end=atoi(getenv("END"));
        string dir=getenv("DIR");
        TString name;
        for(int i=start;i<end;i++){
                name=Form("/scratch/xuq7/flow/pbsjoboutput/qfit/tracknormcpt03to6/coarsebins/%s/Anaveta_Sum_%d.root",dir.c_str(),i);
                remove(name.Data());
                qfit *l = new qfit(readline("filelist.dat",i));
                cout<<"start "<<i<<" th job"<<endl;
                l->beginJob(0);
                l->calcv("mergedV_Prod.root",0,-1);
                l->endJobv(name);
        }
}


int getfilenumber(TString inf, int isample){
	if(isample==-1) return -1;
	TFile *f = TFile::Open(inf);
	int xbin=0;
	TVectorD* IFILE = (TVectorD*)f->Get(Form("D_%d/IFILE",xbin));
	return (int)(*IFILE)[isample];
}

int getnumber(char* name){
        std::ifstream backstory (name);
        std::string line;
        int nlines=0;
        if (backstory.is_open())
        {
                while (std::getline(backstory, line))
                        nlines++;
                backstory.close();
                return nlines;
        }
        else
        {
                std::cout << "Unable to open file" << std::endl << std::endl;
                return 0;
        }

}

std::string readline(char* name, int iline){
        std::ifstream backstory(name);
        std::string line;
        if (backstory.is_open())
                if(backstory.good()){
                        for(int i = 0; i < iline+1; ++i)
                           getline(backstory, line);
                }
        return line;
}
