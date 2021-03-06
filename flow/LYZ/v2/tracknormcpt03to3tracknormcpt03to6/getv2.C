#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
#include "par.h"
void getv2(){

const int ntotbin=5;
const int trkpointmin[ntotbin] = {120,150,185,220,260};
const int trkpointmax[ntotbin] = {150,185,220,260,300};
int ibin=0;
TFile *fout = new TFile(Form("v2diffdataLYZN%d.root",nsamples),"Recreate");
for(int i=0;i<ntotbin;i++){
	TFile *fProd = TFile::Open(Form("M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	TVectorD *vecDv2_Prod = (TVectorD*)fProd->Get(Form("D_%d/vmeanmean",ibin));
	TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/sigmavmeanmean",ibin));
	TVectorD *vecDavgpt_Prod = (TVectorD*)fProd->Get(Form("D_%d/avgavgpt",ibin));

	TFile *fProderr = TFile::Open(Form("M%d%d/mergedv_Prod2_eta_sub.root",trkpointmax[i],trkpointmin[i]));
	TVectorD *vecDv2_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/vmeanmean",ibin));
	TVectorD *vecDv2err_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/sigmavmeanmean",ibin));
	TVectorD *vecDavgeta_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/avgavgeta",ibin));

	double *avgeta_Proderr = vecDavgeta_Proderr->GetMatrixArray();
	double *v2_Proderr = vecDv2_Proderr->GetMatrixArray();
	double *v2err_Proderr = vecDv2err_Proderr->GetMatrixArray();
	int neta = vecDavgeta_Proderr->GetNrows();

	double *avgpt_Prod = vecDavgpt_Prod->GetMatrixArray();
	double *v2_Prod = vecDv2_Prod->GetMatrixArray();
	double *v2err_Prod = vecDv2err_Prod->GetMatrixArray();
	int npt = vecDavgpt_Prod->GetNrows();
	TGraphErrors *grProd=new TGraphErrors(npt,avgpt_Prod,v2_Prod,0,v2err_Prod);
	TGraphErrors *grProderr=new TGraphErrors(neta,avgeta_Proderr,v2_Proderr,0,v2err_Proderr);
        fout->cd();
	grProd->Write(Form("v2vspt%d_%d",trkpointmin[i],trkpointmax[i]));
	grProderr->Write(Form("v2vseta%d_%d",trkpointmin[i],trkpointmax[i]));
}
}

