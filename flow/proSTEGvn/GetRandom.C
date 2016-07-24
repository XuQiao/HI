#include "TMath.h"
#include <stdlib.h>
#include <time.h>
Double_t GetRandom1(TF1 *f){
Double_t*	fAlpha;	//!Array alpha. for each bin in x the deconvolution r of fIntegral
Double_t*	fBeta;	//!Array beta. is approximated by x = alpha +beta*r *gamma*r**2
Double_t*	fGamma;	//!Array gamma.
Double_t*	fIntegral;	//![fNpx] Integral of function binned on fNpx bins
Int_t	fNpx;	//Number of points used for the graphical representation
Double_t	fXmax;	//Upper bounds for the range
Double_t	fXmin;	//Lower bounds for the range

 fXmin      = f->GetXmin();
 fXmax      = f->GetXmax();
   fNpx       = 100;
   fIntegral   = 0;
   fAlpha      = 0;
   fBeta       = 0;
   fGamma      = 0;

	

 if (fIntegral == 0) {
      fIntegral = new Double_t[fNpx+1];
      fAlpha    = new Double_t[fNpx+1];
      fBeta     = new Double_t[fNpx];
      fGamma    = new Double_t[fNpx];
      fIntegral[0] = 0;
      fAlpha[fNpx] = 0;
      Double_t integ;
      Int_t intNegative = 0;
      Int_t i;
      Bool_t logbin = kFALSE;
      Double_t dx;
      Double_t xmin = fXmin;
      Double_t xmax = fXmax;
      if (xmin > 0 && xmax/xmin> fNpx) {
         logbin =  kTRUE;
         fAlpha[fNpx] = 1;
         xmin = TMath::Log10(fXmin);
         xmax = TMath::Log10(fXmax);
      }
      dx = (xmax-xmin)/fNpx;
         
      Double_t *xx = new Double_t[fNpx+1];
      for (i=0;i<fNpx;i++) {
            xx[i] = xmin +i*dx;
      }
      xx[fNpx] = xmax;
      for (i=0;i<fNpx;i++) {
         if (logbin) {
            integ = f->Integral(TMath::Power(10,xx[i]), TMath::Power(10,xx[i+1]));
         } else {
            integ = f->Integral(xx[i],xx[i+1]);
         }
         if (integ < 0) {intNegative++; integ = -integ;}
         fIntegral[i+1] = fIntegral[i] + integ;
      }
      if (intNegative > 0) {
         Warning("GetRandom","function:%s has %d negative values: abs assumed",f->GetName(),intNegative);
      }
      if (fIntegral[fNpx] == 0) {
         delete [] xx;
         Error("GetRandom","Integral of function is zero");
         return 0;
      }
      Double_t total = fIntegral[fNpx];
      for (i=1;i<=fNpx;i++) {  // normalize integral to 1
         fIntegral[i] /= total;
      }
      //the integral r for each bin is approximated by a parabola
      //  x = alpha + beta*r +gamma*r**2
      // compute the coefficients alpha, beta, gamma for each bin
      Double_t x0,r1,r2,r3;
      for (i=0;i<fNpx;i++) {
         x0 = xx[i];
         r2 = fIntegral[i+1] - fIntegral[i];
         if (logbin) r1 = f->Integral(TMath::Power(10,x0),TMath::Power(10,x0+0.5*dx))/total;
         else        r1 = f->Integral(x0,x0+0.5*dx)/total;
         r3 = 2*r2 - 4*r1;
         if (TMath::Abs(r3) > 1e-8) fGamma[i] = r3/(dx*dx);
         else           fGamma[i] = 0;
         fBeta[i]  = r2/dx - fGamma[i]*dx;
         fAlpha[i] = x0;
         fGamma[i] *= 2;
      }
      delete [] xx;
   }

   // return random number
   //Double_t r  = gRandom->Rndm();
   TRandom3 *rnd0 = new TRandom3(0);
   //Double_t r = rnd0->Uniform(0.0,1.0);
   Double_t r = rnd0->Rndm();
   //srand(time(0));
   //Double_t r = (double)rand()/(double)(RAND_MAX+1);
   Int_t bin  = TMath::BinarySearch(fNpx,fIntegral,r);
   Double_t rr = r - fIntegral[bin];

   Double_t yy;
   if(fGamma[bin] != 0)
      yy = (-fBeta[bin] + TMath::Sqrt(fBeta[bin]*fBeta[bin]+2*fGamma[bin]*rr))/fGamma[bin];
   else
      yy = rr/fBeta[bin];
   Double_t x = fAlpha[bin] + yy;
   if (fAlpha[fNpx] > 0) return TMath::Power(10,x);
   return x;
}
