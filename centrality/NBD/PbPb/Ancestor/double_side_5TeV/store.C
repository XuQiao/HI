#include "../parameter.h"
#include "par.h"
void store(int type, int sth=0, int Gth=0){
	TH1::AddDirectory(kFALSE); 
	NBD *l;
	struct para1 var;
        TString dirname;
        double sss=0.1;
	if(sth==0){
	dirname = "std";
	if(Gth<nGlau)
            l = new NBD(datafile,stdGlaulist[Gth],histoname);
	else
            l = new NBD(datafile,stdGlaulist[0],histoname);
	var = var1;
	//l->initmu(bestlist1[Gth].mubest-sss,bestlist1[Gth].mubest+sss,sss/5);
	//l->initk(bestlist1[Gth].kbest-sss,bestlist1[Gth].kbest+sss,sss/5);
	//l->initf(bestlist1[Gth].fbest-sss,bestlist1[Gth].fbest+sss,sss/5);
	}
        
	else if(sth==1){
	dirname ="Gri055"; 
	if(Gth<nGlau)
            l = new NBD(datafile,Gri055Glaulist[Gth],histoname);
	else
            l = new NBD(datafile,Gri055Glaulist[0],histoname);
	var = var2;
	//l->initmu(bestlist2[Gth].mubest-sss,bestlist2[Gth].mubest+sss,sss/5);
	//l->initk(bestlist2[Gth].kbest-sss,bestlist2[Gth].kbest+sss,sss/5);
	//l->initf(bestlist2[Gth].fbest-sss,bestlist2[Gth].fbest+sss,sss/5);
	}
        
	else {
	dirname ="Gri101";
	if(Gth<nGlau)
            l = new NBD(datafile,Gri101Glaulist[Gth],histoname);
	else
            l = new NBD(datafile,Gri101Glaulist[0],histoname);
        var = var3;
	//l->initmu(bestlist3[Gth].mubest-sss,bestlist3[Gth].mubest+sss,sss/5);
	//l->initk(bestlist3[Gth].kbest-sss,bestlist3[Gth].kbest+sss,sss/5);
	//l->initf(bestlist3[Gth].fbest-sss,bestlist3[Gth].fbest+sss,sss/5);
	}

	l->initmu(var.mumin,var.mumax,var.mustep);
    	l->initk(var.kmin,var.kmax,var.kstep);
    	l->initf(var.fmin,var.fmax,var.fstep);
	if(Gth<nGlau)
        	l->initx(var.xmin,var.xmax);
	else if(Gth-nGlau==0)
		l->initx(var.xmin-binshift,var.xmax);
	else
		l->initx(var.xmin+binshift,var.xmax);
if(type==0)
	l->fit();	
if(type==1){
 	if(sth==0) l->assign(bestlist1[Gth].mubest,bestlist1[Gth].kbest,bestlist1[Gth].fbest);
	else if(sth==1) l->assign(bestlist2[Gth].mubest,bestlist2[Gth].kbest,bestlist2[Gth].fbest);
	else l->assign(bestlist3[Gth].mubest,bestlist3[Gth].kbest,bestlist3[Gth].fbest);

        l->initN(bin,N,method);
        l->calcvar();
	TFile *outfile = new TFile(outG,"Update");
	if(Gth==0 && sth==0){
	l->dataname.Write("dataname",TObject::kOverwrite);
	l->histoname.Write("histoname",TObject::kOverwrite);
	}
	TDirectory *dir = (TDirectory*)outfile->GetDirectory(dirname);
	if (!dir) {outfile->mkdir(dirname);	dir = (TDirectory*)outfile->GetDirectory(dirname);}
	dir->cd();
	TString name;
	if(Gth==0)
		name = "G0";
	else if(Gth<nGlau)
		name = Form("Glau_%d",Gth);
	else
		name = Form("bin_%d",Gth-nGlau+1);
	TDirectory *subdir = (TDirectory*)dir->GetDirectory(name);
	if(!subdir) {dir->mkdir(name);	subdir = (TDirectory*)dir->GetDirectory(name);}
	subdir->cd();
	l->method.Write("method",TObject::kOverwrite);
	l->Glaubername.Write("Glaubername",TObject::kOverwrite);
	l->xmin.Write("xmin",TObject::kOverwrite);l->xmax.Write("xmax",TObject::kOverwrite);
	l->mubest.Write("mubest",TObject::kOverwrite);l->kbest.Write("kbest",TObject::kOverwrite);l->fbest.Write("fbest",TObject::kOverwrite);
	l->chis.Write("chis",TObject::kOverwrite);l->Ndf.Write("Ndf",TObject::kOverwrite);
	l->NcollAver.Write("NcollAver",TObject::kOverwrite);l->NpartAver.Write("NpartAver",TObject::kOverwrite);l->BAver.Write("BAver",TObject::kOverwrite);
	l->NcollRMS.Write("NcollRMS",TObject::kOverwrite);l->NpartRMS.Write("NpartRMS",TObject::kOverwrite);l->BRMS.Write("BRMS",TObject::kOverwrite);
	l->centbin.Write("centbin",TObject::kOverwrite);	l->kpoint.Write("kpoint",TObject::kOverwrite);
	l->centbin_.Write("centbin_",TObject::kOverwrite);	l->kpoint_.Write("kpoint_",TObject::kOverwrite);
	l->Npartdis->Write("Npartdis",TObject::kSingleKey | TObject::kOverwrite);
	l->Grgrid->Write("Grgrid",TObject::kOverwrite);
        outfile->Close();
}
}
