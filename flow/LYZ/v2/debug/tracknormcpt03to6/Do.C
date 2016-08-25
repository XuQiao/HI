#include "TSystem.h"
#include "TRint.h"
#include "TROOT.h"
#include <fstream>
#include <iostream>

void dofirst(){
	remove("GRe.txt");
	remove("GIm.txt");
	remove("Q2.txt");
	remove("Qx1.txt");
	remove("Qy1.txt");
	remove("mult.txt");
	TString name;
	string SumorProd="Prod";
	int start=0;
	int end=1;
	for(int i=start;i<end;i++){
		if(SumorProd=="Sum")	name=Form("AnaV_Sum_%d.root",i);
		else	name=Form("AnaV_Prod_%d.root",i);
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
	string SumorProd="Prod";
        int start=0;
        int end=1;
        TString name;
        for(int i=start;i<end;i++){
                if(SumorProd=="Sum")       name=Form("Anav_Prod_%d.root",i);
                else            name=Form("Anav_Prod2_%d.root",i);
                remove(name.Data());
                LYZ *l = new LYZ(readline("filelist.dat",i));
                cout<<"start "<<i<<" th job"<<endl;
                l->beginJob();
                if(SumorProd=="Sum")       l->calcv("mergedV_Sum.root",0);
                else    l->calcv("mergedV_Prod.root",0);
                l->endJobv(name);
        }
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
