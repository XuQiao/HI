#include "TSystem.h"
#include "TRint.h"
#include "TROOT.h"
#include <fstream>
#include <iostream>
#include "par.h"

void dofirst(){
	remove("GRe.txt");
	remove("GIm.txt");
	remove("Q2.txt");
	remove("Qx1.txt");
	remove("Qy1.txt");
	remove("mult.txt");
	TString name;
	string SumorProd=getenv("SUMORPROD");
	int start=atoi(getenv("BEGIN"));
	int end=atoi(getenv("END"));
	string dir=getenv("DIR");
	for(int i=start;i<end;i++){
		if(SumorProd=="Sum")	name=Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/tracktc/%s/AnaV_Sum_%d.root",dir.c_str(),i);
		else	name=Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/tracktc/%s/AnaV_Prod_%d.root",dir.c_str(),i);
		remove(name.Data());
		LYZ *l = new LYZ(readline("./filelist.dat",i));
		cout<<"start "<<i<<" th job"<<endl;
		l->beginJob();
		if(SumorProd=="Sum")	l->calcV(1);
		else l->calcV(0);
		l->endJobV(name);
	}
}

void dosecond(){
        remove("dDRe.txt");
        remove("dDIm.txt");
        remove("dNRe.txt");
        remove("dNIm.txt");
        remove("mult_.txt");
	string SumorProd=getenv("SUMORPROD");
        int start=atoi(getenv("BEGIN"));
        int end=atoi(getenv("END"));
	string dir=getenv("DIR");
        TString name;
        for(int i=start;i<end;i++){
                if(SumorProd=="Sum")       name=Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/tracktc/%s/Anav_Prod_%d.root",dir.c_str(),i);
                else            name=Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/tracktc/%s/Anav_Prod2_%d.root",dir.c_str(),i);
                remove(name.Data());
                LYZ *l = new LYZ(readline("filelist.dat",i));
                cout<<"start "<<i<<" th job"<<endl;
                l->beginJob(1);
                if(SumorProd=="Sum")       l->calcv("mergedV_Sum.root",0,-1);
                else    l->calcv("mergedV_Prod.root",0,-1);
                l->endJobv(name);
        }
}

void dosecondsub(){
        remove("dDRe.txt");
        remove("dDIm.txt");
        remove("dNRe.txt");
        remove("dNIm.txt");
        remove("mult_.txt");
        string SumorProd=getenv("SUMORPROD");
        int start=atoi(getenv("BEGIN"));
        int end=atoi(getenv("END"));
        string dir=getenv("DIR");
        TString name, infname;	
	TFile *inf;
   	if(SumorProd=="Sum")       infname="mergedV_Sum_sub.root";
        else    infname="mergedV_Prod_sub.root";
        for(int i=start;i<end;i++){
                if(SumorProd=="Sum")       name=Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/tracktc/%s/Anav_Prod_sub_%d.root",dir.c_str(),i);
                else            name=Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/tracktc/%s/Anav_Prod2_sub_%d.root",dir.c_str(),i);
                remove(name.Data());
                LYZ *l = new LYZ(readline("filelist.dat",i));
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
        string SumorProd=getenv("SUMORPROD");
        int start=atoi(getenv("BEGIN"));
        int end=atoi(getenv("END"));
        string dir=getenv("DIR");
        TString name;
        for(int i=start;i<end;i++){
                if(SumorProd=="Sum")       name=Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/tracktc/%s/Anaveta_Prod_%d.root",dir.c_str(),i);
                else            name=Form("/scratch/xuq7/flow/pbsjoboutput/tracknormcpt03to3tracknormcpt03to6/tracktc/%s/Anaveta_Prod2_%d.root",dir.c_str(),i);
                remove(name.Data());
                LYZ *l = new LYZ(readline("filelist.dat",i));
                cout<<"start "<<i<<" th job"<<endl;
                l->beginJob(0);
                if(SumorProd=="Sum")       l->calcv("mergedV_Sum.root",0,-1);
                else    l->calcv("mergedV_Prod.root",0,-1);
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
