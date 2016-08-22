#include "../parameter.h"
#include "par.h"
void store(){
	 if (!TClass::GetDict("NBD")) {
                gROOT->ProcessLine(".L /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/NBDclass.C+");
        }
	TH1::AddDirectory(kFALSE); 
	int Gth=atoi(getenv("GTH"));
	int sth=atoi(getenv("STH"));
	NBD *l;
	struct para1 var;
        if(sth==0){TString dirname = "std"; var = var1; if(Gth<nGlau) l = new NBD(datafile,stdGlaulist[Gth],histoname);else l = new NBD(datafile,stdGlaulist[0],histoname);}
        else if(sth==1){TString dirname ="Gri055"; var = var2; if(Gth<nGlau) l = new NBD(datafile,Gri055Glaulist[Gth],histoname);else l = new NBD(datafile,Gri055Glaulist[0],histoname);}
        else {TString dirname ="Gri101"; var = var3; if(Gth<nGlau) l = new NBD(datafile,Gri101Glaulist[Gth],histoname);else l = new NBD(datafile,Gri101Glaulist[0],histoname);}
	l->initmu(var.mumin,var.mumax,var.mustep);
     	l->initk(var.kmin,var.kmax,var.kstep);
	if(Gth<nGlau)
        	l->initx(var.xmin,var.xmax);
	else if(Gth-nGlau==0)
		l->initx(var.xmin-binshift,var.xmax);
	else
		l->initx(var.xmin+binshift,var.xmax);
//	l->fit();

 	if(sth==0) l->assign(bestlist1[Gth].mubest,bestlist1[Gth].kbest);
	else if(sth==1) l->assign(bestlist2[Gth].mubest,bestlist2[Gth].kbest);
	else l->assign(bestlist3[Gth].mubest,bestlist3[Gth].kbest);

        l->initN(bin,N,method);
        l->calcvar();
	TFile *outfile = new TFile(outG,"Update");
	if(Gth==0 && sth==0){
	l->dataname.Write("dataname",TObject::kOverwrite);
	l->histoname.Write("histoname",TObject::kOverwrite);
	}
	TDirectory *dir = outfile->GetDirectory(dirname);
	if (!dir) {outfile->mkdir(dirname);	TDirectory *dir = outfile->GetDirectory(dirname);}
	dir->cd();
	TString name;
	if(Gth==0)
		name = "G0";
	else if(Gth<nGlau)
		name = Form("Glau_%d",Gth);
	else
		name = Form("bin_%d",Gth-nGlau+1);
	TDirectory *subdir = dir->GetDirectory(name);
	if(!subdir) {dir->mkdir(name);	TDirectory *subdir = dir->GetDirectory(name);}
	subdir->cd();
	l->method.Write("method",TObject::kOverwrite);
	l->Glaubername.Write("Glaubername",TObject::kOverwrite);
	l->xmin.Write("xmin",TObject::kOverwrite);l->xmax.Write("xmax",TObject::kOverwrite);
	l->mubest.Write("mubest",TObject::kOverwrite);l->kbest.Write("kbest",TObject::kOverwrite);
	l->chis.Write("chis",TObject::kOverwrite);l->Ndf.Write("Ndf",TObject::kOverwrite);
	l->NcollAver.Write("NcollAver",TObject::kOverwrite);l->NpartAver.Write("NpartAver",TObject::kOverwrite);l->BAver.Write("BAver",TObject::kOverwrite);
	l->centbin.Write("centbin",TObject::kOverwrite);	l->kpoint.Write("kpoint",TObject::kOverwrite);
	l->centbin_.Write("centbin_",TObject::kOverwrite);	l->kpoint_.Write("kpoint_",TObject::kOverwrite);
	l->Npartdis->Write();
        l->Grgrid->Write("Grgrid",TObject::kOverwrite);
        outfile->Close();

}
