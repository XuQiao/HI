#include "par.h"
#include <iomanip.h>
void writev2vspt_sys(){

const int ntotbin=4;
ofstream fout("PbPbsys.dat");

int shift=2;
for(int i=shift;i<shift+ntotbin;i++){
        int t = ntotbin-(i-shift);
	//TFile *fProderr = TFile::Open(Form("M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProd = TFile::Open(Form("mergedv_Prod2.root"));
	TFile *fProdtc = TFile::Open(Form("../tracktc/PbPb2011MB/mergedv_Prod2.root"));
	TFile *fProdlc = TFile::Open(Form("../tracklc/PbPb2011MB/mergedv_Prod2.root"));
	TFile *fProdzl = TFile::Open(Form("../trackzvtxl/PbPb2011MB/mergedv_Prod2.root"));
	TFile *fProdzs = TFile::Open(Form("../trackzvtxs/PbPb2011MB/mergedv_Prod2.root"));
	//TVectorD *vecDv2_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/vmeanmean",ibin));
	//TVectorD *vecDv2err_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/sigmavmeanmean",ibin));
	//TVectorD *vecDavgpt_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/avgavgpt",ibin));
	
	TVectorD *vecDv2_Prod = (TVectorD*)fProd->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prod = (TVectorD*)fProd->Get(Form("D_%d/avgpt",i));

	TVectorD *vecDv2_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/avgpt",i));

	TVectorD *vecDv2_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/avgpt",i));
	
        TVectorD *vecDv2_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/avgpt",i));

	TVectorD *vecDv2_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/avgpt",i));

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

