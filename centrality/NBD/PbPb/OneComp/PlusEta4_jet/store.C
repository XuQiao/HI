#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include "par.h"
void store(){
	if (!TClass::GetDict("NBD")) {
		gROOT->ProcessLine(".L /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/NBDclass.C+");
	}
	int Gth=atoi(getenv("GTH"));
	int sth=atoi(getenv("STH"));
	NBD *l;
	struct para1 var;
        if(sth==0){TString dirname = "std"; var = var1; l = new NBD(datafile,stdGlaulist[Gth],Gth,sth,histoname);}
        else if(sth==1){TString dirname ="Gri055"; var = var2; l = new NBD(datafile,Gri055Glaulist[Gth],Gth,sth,histoname);}
        else {TString dirname ="Gri101"; var = var3; l = new NBD(datafile,Gri101Glaulist[Gth],Gth,sth,histoname);}
        l->initmu(var.mumin,var.mumax,var.mustep);
        l->initk(var.kmin,var.kmax,var.kstep);
        l->initx(var.xmin,var.xmax);
//	l->fit();
	if(sth==0) l->assign(bestlist1[Gth].mubest,bestlist1[Gth].kbest);
	else if(sth==1) l->assign(bestlist2[Gth].mubest,bestlist2[Gth].kbest);
	else l->assign(bestlist3[Gth].mubest,bestlist3[Gth].kbest);

        l->initN(bin,N,method);
        l->calcvar();
/*	TFile *outfile = new TFile("G.root","Update");
	TDirectory *dir = outfile->GetDirectory(dirname);
	if (dir) dir->cd();else {outfile->mkdir(dirname);TDirectory *dir = outfile->GetDirectory(dirname);
//outfile->cd(dirname);
}
	dir->WriteObject(l,Form("Glau%d",Gth));
	//l->SetName(Form("Glau%d",Gth));
*/
        TFile *outfile = new TFile(Form("G%d%d.root",sth,Gth),"Recreate");
        //l->SetName(Form("Glau"));

  	l->Write();
        outfile->Close();
}
