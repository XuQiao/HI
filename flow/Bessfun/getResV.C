#include "par.h"
void getResV(){
	TVectorD r;
	//TVectorD r0, r01, V, chi;
	TVectorD Gmod2;
	TVectorD GRe; TVectorD* GRe_t;
	TVectorD GIm; TVectorD* GIm_t;
	TComplex G[nstepr];
            r.ResizeTo(nstepr);
	    Gmod2.ResizeTo(nstepr);
            GRe.ResizeTo(nstepr);
	    GRe.Zero();
	    GIm.ResizeTo(nstepr);
	    GIm.Zero();
	ofstream  fstrV;
	fstrV.open("V.txt"); 

    TFile *f = TFile::Open("Qdis.root");
    TFile *fout = new TFile("mergedV.root","Recreate");
    GRe_t = (TVectorD*)f->Get(Form("GRe"));
    GIm_t = (TVectorD*)f->Get(Form("GIm"));
    for(int ir=0; ir<nstepr; ir++){
        GRe[ir] += (*GRe_t)[ir];
        GIm[ir] += (*GIm_t)[ir];
    }
    for(int ir=0; ir<nstepr; ir++)
        r[ir]=0.00025*20*(ir+1);
    for(ir=0; ir<nstepr; ir++){
        G[ir]=TComplex(GRe[ir],GIm[ir]);
	//G[ir]/=Nevent[ibin];
        Gmod2[ir]=TMath::Power(TComplex::Abs(G[ir]),2);
        //Gmod2[ir]=TMath::Power(GRe[ir],2);
    }
    for(ir=0; ir<nstepr-1; ir++)
	if(ir!=0 && Gmod2[ir]<=Gmod2[ir-1] && Gmod2[ir]<=Gmod2[ir+1]) break;
    if(ir!=0 && ir<nstepr-1)	double r01=r[ir];
    else if(ir==0)	{cout<<"ibin="<<ibin<<"\t"<<"ietabin="<<ietabin<<"\t"<<"itheta="<<itheta<<"\tminimum lies on ir = 0, please select proper range!"<<endl;	continue;}
    else 	{cout<<"ibin="<<ibin<<"\t"<<"ietabin="<<ietabin<<"\t"<<"itheta="<<itheta<<"\tminimum lies on ir = maximum "<<nstepr-1<<", please select proper range!"<<endl;	continue;}
        double V=j01/r01;
        V/=150;
	fstrV<<r01<<"\t"<<V<<endl;
        fout->cd();
        r.Write("r");
        GRe.Write("GRe");
        GIm.Write("GIm");
        Gmod2.Write("G2");
        fout->Close();
}
  
