const int ntotbin=5;
const int trkpointmin[ntotbin] = {120,150,185,220,260};
const int trkpointmax[ntotbin] = {150,185,220,260,300};
const int ntheta=5;
const int ndir=6;
const int nbin = 1;
const double j01=2.404826;
TString dir[ndir] ={ "tracknormcpt03to6","tracknormcpt03to3","tracklcpt01to10","PFcandpt03to6tracknormcpt03to6","PFcandpt03to3tracknormcpt03to6","PFcandpt01to10tracknormcpt03to6"};
double *avgmult;
void calcinpar(){
    TVectorD r00[nbin];   TVectorD V[nbin];    TVectorD Vmean[nbin]; TVector Vcalc[nbin]; TVectorD Vcalcmean[nbin];
    double V2[ndir][ntheta], V2int[ndir][ntheta];
    double r0[ndir][ntheta], r0int[ndir][ntheta];
    double V2x[ndir];
    double r0x[ndir];
 //   TVectorD *avgmult[nbin];
    for(int i = 0; i<ntotbin ;i++){
 //       TFile *fin = TFile::Open(Form("M%d%d/mergedV_Prod.root",trkpointmax[i],trkpointmin[i]));
        
    for(int ibin=0;ibin<nbin;ibin++){
 //       avgmult[ibin] = (TVectorD*)fin->Get(Form("D_%d/avgmult",ibin));
        r00[ibin].ResizeTo(ntheta);
        Vmean[ibin].ResizeTo(1);
        Vcalcmean[ibin].ResizeTo(1);
        Vcalc[ibin].ResizeTo(ntheta);
        V[ibin].ResizeTo(ntheta);
        for(int idir=0;idir<ndir;idir++){
            double *V2ref = read(dir[idir],0,i,1);
            double *r0ref = read(dir[idir],0,i,0);
            for(int itheta=0; itheta<ntheta; itheta++){
                V2[idir][itheta]=V2ref[itheta];
                r0[idir][itheta]=r0ref[itheta];
            }
        }
        for(int itheta=0; itheta<ntheta; itheta++){
            for(int idir=0;idir<ndir;idir++){
               V2x[idir] = V2[idir][itheta];
               r0x[idir] = r0[idir][itheta];
            }
               double stepV2 = (TMath::MaxElement(ndir,V2x)-TMath::MinElement(ndir,V2x))/(ndir-3);
               double stepr0 = (TMath::MaxElement(ndir,r0x)-TMath::MinElement(ndir,r0x))/(ndir-3);
            for(int idir=0;idir<ndir;idir++){
               V2int[idir][itheta] = TMath::MinElement(ndir,V2x) + (idir-1) * stepV2;
               r0int[idir][itheta] = TMath::MinElement(ndir,r0x) + (idir-1) * stepr0;
            }
        }
    for(int iloop=0; iloop<ndir; iloop++){
        TFile *fout = new TFile(Form("M%d%d/parin_%d.root",trkpointmax[i],trkpointmin[i],iloop),"recreate");
        for(int itheta=0; itheta<ntheta; itheta++){
            V[ibin][itheta] = V2int[iloop][itheta];
            r00[ibin][itheta] = r0int[iloop][itheta];
            Vcalc[ibin][itheta] = j01/avgmult[ibin]/r0int[iloop][itheta];
 //           r0[ibin][itheta] = j01/V2int[iloop][itheta]/(*avgmult[ibin])[0];
        }
        Vmean[ibin][0] = TMath::Mean(ntheta,V2int[iloop]);
        Vcalcmean[ibin][0] = Vcalc[ibin].Sum()/ntheta;
        TDirectory *dir0 = fout->mkdir(Form("D_%d",ibin));
        dir0->cd();
        Vmean[ibin].Write("Vmean");
        Vcalcmean[ibin].Write("Vcalcmean");
        TDirectory *dir1 = dir0->mkdir("D_0");
        dir1->cd();
        r00[ibin].Write("r0");
        V[ibin].Write("V");
        Vcalc[ibin].Write("Vcalc");
        fout->Close();
    }
    }
    }

}

double *read(TString dir, int isSum, int i, int isV2){

    if(isSum==0) TFile *f = TFile::Open(Form("../../%s/M%d%d/mergedV_Prod.root",dir.Data(),trkpointmax[i],trkpointmin[i]));
    else TFile *f = TFile::Open(Form("../../%s/M%d%d/mergedV_Sum.root",dir.Data(),trkpointmax[i],trkpointmin[i]));
    int xbin=0;
    int xpt=0;
    const int nn = 2;
        TVectorD *vecNtrk = (TVectorD*)f->Get("avgtrk");
        TVectorD *vecavgmult = (TVectorD*)f->Get("avgmultall");
        if(dir=="tracknormcpt03to6"){
        TVectorD *vecV2mean_0=(TVectorD*)f->Get(Form("Vmean",xbin));
        TVectorD *vecV2_0=(TVectorD*)f->Get(Form("D_%d/V",xbin));
        TVectorD *vecr0_0=(TVectorD*)f->Get(Form("D_%d/r0",xbin));
        }
        else{
        TVectorD *vecV2mean_0=(TVectorD*)f->Get(Form("D_%d/Vmean",xbin));
        TVectorD *vecV2_0=(TVectorD*)f->Get(Form("D_%d/D_%d/V",xbin,xpt));
        TVectorD *vecr0_0=(TVectorD*)f->Get(Form("D_%d/D_%d/r0",xbin,xpt));
        }
      //  TVectorD *vecV2err_0=(TVectorD*)f->Get(Form("D_%d/D_%d/deltaV",xbin,xpt));
        double *V2_0=vecV2_0->GetMatrixArray();
        double *r0_0=vecr0_0->GetMatrixArray();
        if(dir=="tracknormcpt03to6")
        avgmult = vecavgmult->GetMatrixArray();
      //  double *V2err_0=vecV2err_0->GetMatrixArray();
        double theta[ntheta];
        for(int itheta=0;itheta<ntheta;itheta++){
            theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
        }
        
        f->Close();
        if(isV2)
            return V2_0;
        else
            return r0_0;
            
}
