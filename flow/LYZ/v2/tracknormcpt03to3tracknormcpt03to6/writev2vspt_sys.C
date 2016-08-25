#include "par.h"
#include <iomanip.h>
void writev2vspt_sys(){

const int ntotbin=5;
const int xbin=0;
ofstream fout("pPbsys.dat");
const int trkpointmin[ntotbin] = {120,150,185,220,260};
const int trkpointmax[ntotbin] = {150,185,220,260,300};

for(int i=0;i<ntotbin;i++){
	TFile *fProd = TFile::Open(Form("M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i],trkpointmax[i],trkpointmin[i]));
	TFile *fProdtc = TFile::Open(Form("../tracktc/tracknormcpt03to3tracknormcpt03to6/M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProdlc = TFile::Open(Form("../tracklc/tracknormcpt03to3tracknormcpt03to6/M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProdzl = TFile::Open(Form("../trackzvtxl/tracknormcpt03to3tracknormcpt03to6/M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProdzs = TFile::Open(Form("../trackzvtxs/tracknormcpt03to3tracknormcpt03to6/M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	//TVectorD *vecDv2_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/vmeanmean",ibin));
	//TVectorD *vecDv2err_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/sigmavmeanmean",ibin));
	//TVectorD *vecDavgpt_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/avgavgavgpt",ibin));
	
	TVectorD *vecDv2_Prod = (TVectorD*)fProd->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prod = (TVectorD*)fProd->Get(Form("D_%d/avgavgpt",xbin));

	TVectorD *vecDv2_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/avgavgpt",xbin));

	TVectorD *vecDv2_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/avgavgpt",xbin));
	
        TVectorD *vecDv2_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/avgavgpt",xbin));

	TVectorD *vecDv2_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/avgavgpt",xbin));

	//double *avgpt_Proderr = vecDavgpt_Proderr->GetMatrixArray();
	//double *v2_Proderr = vecDv2_Proderr->GetMatrixArray();
	//double *v2err_Proderr = vecDv2err_Proderr->GetMatrixArray();

	double *avgpt_Prod = vecDavgpt_Prod->GetMatrixArray();
	double *v2_Prod = vecDv2_Prod->GetMatrixArray();
	double *v2err_Prod = vecDv2err_Prod->GetMatrixArray();
	double *v2_Prodtc = vecDv2_Prodtc->GetMatrixArray();
	double *v2err_Prodtc = vecDv2err_Prodtc->GetMatrixArray();
	double *v2_Prodlc = vecDv2_Prodlc->GetMatrixArray();
	double *v2err_Prodlc = vecDv2err_Prodlc->GetMatrixArray();
	double *v2_Prodzl = vecDv2_Prodzl->GetMatrixArray();
	double *v2err_Prodzl = vecDv2err_Prodzl->GetMatrixArray();
	double *v2_Prodzs = vecDv2_Prodzs->GetMatrixArray();
	double *v2err_Prodzs = vecDv2err_Prodzs->GetMatrixArray();
	int npt = vecDavgpt_Prod->GetNrows();
        fout<<fixed<<setprecision(4)<<endl;
        fout<<"avgpt\tv2\tv2 err stat.\tv2 tight\tv2 loose\t3<|z|<15\t|z|<3"<<endl;
        for(int iptbin=0;iptbin<npt;iptbin++){
            fout<<avgpt_Prod[iptbin]<<"\t";
            fout<<v2_Prod[iptbin]<<"\t";
            fout<<v2err_Prod[iptbin]<<"\t";
            fout<<v2_Prodtc[iptbin]<<"\t";
            fout<<v2_Prodlc[iptbin]<<"\t";
            fout<<v2_Prodzl[iptbin]<<"\t";
            fout<<v2_Prodzs[iptbin]<<"\t";
            fout<<endl;
        }
	
        }
}

