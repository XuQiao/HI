#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
void checkmatrix(){
        static const double boundaries_truth[] = {
            3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000
            //
            //3, 4, 5, 7, 9, 12, 15, 18,21,24,32,43,56,74,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000
            //
            };
        static const int nbins_truth = sizeof(boundaries_truth)/sizeof(Double_t)-1;
            
        TF1 * ftrue = new TF1("ftrue", "[0]*exp([1]/x)*pow(x,[2])*pow(1-x*cosh([4])/4000.,[3])", 0, 800) ;
        ftrue->SetParameters(7.12152e+07, -9.06678e+01,-6.31860e+00,3.51433e+00,0);
        TF1 *fres = new TF1("fres", "sqrt([0]*[0]+pow(([1]/sqrt(x)),2)+pow([2]/x, 2))", 0, 800);
        int Nevents =1.e6;
        fres->SetParNames("C", "S", "N");
        fres->SetParameters(0.0631308,  1.01648 ,  4.72812e-05);
        TH2F *mt = new TH2F(Form("mt_jet_0"),"",nbins_truth,boundaries_truth,nbins_truth,boundaries_truth);
 //       TH1D *mx = new TH1D(Form("mgen_jet_0"),"",nbins_truth,boundaries_truth);
        TH1D *mx = new TH1D(Form("mgen_jet_0"),"",1000,0,1000);
        TH1D *my = new TH1D(Form("mreco_jet_0"),"",nbins_truth,boundaries_truth);
        TFile *f = new TFile("matrix.root","Recreate");
        ftrue->SetNpx(10000);
        for(int i = 0 ; i < Nevents ; i++){
            if(i%1000000==0) cout<<"Processed "<<i <<" events"<<endl;
            double ptgen = ftrue->GetRandom(20, 800);
            double res = fres->Eval(ptgen)*ptgen;
            double ptreco = gRandom->Gaus(ptgen, res);
         //   if(ptreco<=20) continue;
           // int ibin = mt->GetXaxis()->FindBin(ptgen);
           // int jbin = mt->GetXaxis()->FindBin(ptreco);
            mt->Fill(ptgen, ptreco);
        }
         for(int i = 0 ; i < Nevents ; i++){
            if(i%1000000==0) cout<<"Processed "<<i <<" events"<<endl;
            double ptgen = ftrue->GetRandom(0, 800);
            double res = fres->Eval(ptgen)*ptgen;
            double ptreco = gRandom->Gaus(ptgen, res);
         //   if(ptreco<=20) continue;
            mx->Fill(ptgen);
            my->Fill(ptreco);
        }
        mt->Write();
        mx->Write();
        my->Write();
        f->Close();
}
