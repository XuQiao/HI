void polar(){
Double_t rmin=0;
   Double_t rmax=TMath::Pi()*2;
      Double_t r[1000];
         Double_t theta[1000];

            TF1 * fp1 = new TF1("fplot","1+0.5*sin(2*x)",rmin,rmax);
               for (Int_t ipt = 0; ipt < 1000; ipt++) {
                         r[ipt] = ipt*(rmax-rmin)/1000+rmin;
                               theta[ipt] = fp1->Eval(r[ipt]);
                                  }
   TGraphPolar * grP1 = new TGraphPolar(1000,r,theta);
  // TGraphPolargram *g = grP1->GetPolargram();
  // g->SetRangeRadial(0,2);
  // grP1->SetMaxRadial(0);
      grP1->SetTitle("");
         grP1->SetLineColor(2);
            grP1->Draw("AOP");
            }
