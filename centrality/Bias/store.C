#include "parameter.h"
#include "par.h"
void store(){
	 if (!TClass::GetDict("GammaF")) {
                gROOT->ProcessLine(".L /home/xuq7/HI/centrality/Bias/GammaFit.C+");
        }
	TH1::AddDirectory(kFALSE); 
	int Gth=atoi(getenv("GTH"));
	int sth=atoi(getenv("STH"));
	GammaF *l;
	struct para1 var;
	TString MCfitdat = "/home/xuq7/HI/centrality/Bias/fitresultMC.dat";

	if(sth==0){
	TString dirname = "std";
	if(Gth<nGlau) l = new GammaF(datafile,stdGlaulist[Gth],histoname,MCfitdat);
	else l = new GammaF(datafile,stdGlaulist[0],histoname,MCfitdat);
	var=var1;
	l->inittheta(bestlist1[Gth].thetabest-0.01,bestlist1[Gth].thetabest+0.01,0.001);
	l->initk(bestlist1[Gth].kbest-0.01,bestlist1[Gth].kbest+0.01,0.001);
	}
        
	else if(sth==1){
	TString dirname ="Gri055"; 
	if(Gth<nGlau) l = new GammaF(datafile,Gri055Glaulist[Gth],histoname,MCfitdat);
	else l = new GammaF(datafile,Gri055Glaulist[0],histoname,MCfitdat);
	var=var2;
//	l->inittheta(bestlist2[Gth].thetabest-0.01,bestlist2[Gth].thetabest+0.01,0.001);
//	l->initk(bestlist2[Gth].kbest-0.01,bestlist2[Gth].kbest+0.01,0.001);
	}
        
	else {
	TString dirname ="Gri101";
	if(Gth<nGlau) l = new GammaF(datafile,Gri101Glaulist[Gth],histoname,MCfitdat);
	else l = new GammaF(datafile,Gri101Glaulist[0],histoname,MCfitdat);
	var=var3;
//	l->inittheta(bestlist3[Gth].thetabest-0.01,bestlist3[Gth].thetabest+0.01,0.001);
//	l->initk(bestlist3[Gth].kbest-0.01,bestlist3[Gth].kbest+0.01,0.001);
	}

//	l->inittheta(var.thetamin,var.thetamax,var.thetastep);
//   	l->initk(var.kmin,var.kmax,var.kstep);

	if(Gth<nGlau)
        	l->initx(var.xmin,var.xmax);
	else if(Gth-nGlau==0)
		l->initx(var.xmin-binshift,var.xmax);
	else
		l->initx(var.xmin+binshift,var.xmax);
//	l->fit();	

 	if(sth==0) l->assign(bestlist1[Gth].thetabest,bestlist1[Gth].kbest);
	else if(sth==1) l->assign(bestlist2[Gth].thetabest,bestlist2[Gth].kbest);
	else l->assign(bestlist3[Gth].thetabest,bestlist3[Gth].kbest);

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
	l->k0.Write("k0",TObject::kOverwrite);
	l->theta0.Write("theta0",TObject::kOverwrite);
	l->Glaubername.Write("Glaubername",TObject::kOverwrite);
	l->xmin.Write("xmin",TObject::kOverwrite);l->xmax.Write("xmax",TObject::kOverwrite);
	l->thetabest.Write("thetabest",TObject::kOverwrite);l->kbest.Write("kbest",TObject::kOverwrite);
	l->chis.Write("chis",TObject::kOverwrite);l->Ndf.Write("Ndf",TObject::kOverwrite);
	l->NcollAver.Write("NcollAver",TObject::kOverwrite);l->NpartAver.Write("NpartAver",TObject::kOverwrite);l->BAver.Write("BAver",TObject::kOverwrite);
	l->centbin.Write("centbin",TObject::kOverwrite);	l->kpoint.Write("kpoint",TObject::kOverwrite);
	l->centbin_.Write("centbin_",TObject::kOverwrite);	l->kpoint_.Write("kpoint_",TObject::kOverwrite);
	l->Npartdis->Write("Npartdis",TObject::kOverwrite);
	l->Grgrid->Write("Grgrid",TObject::kOverwrite);
        outfile->Close();

}
