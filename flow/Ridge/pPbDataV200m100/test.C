void test(){
TH1D* h = new TH1D("","",1000,-5,5);
for(int i=0;i<1e5;i++){

   if(i%10000==0) cout<<i<<endl;
  // TF1 * f = new TF1("f","x",0,1);
  TF1 *f = new TF1("f","exp(-(x-2.1)^2/6.3)+exp(-(x+2.1)^2/6.3)",-2.4,2.4);
//TRandom3 *r1 = new TRandom3(0);
//double x = r1->Rndm();
double x = f->GetRandom();
//TRandom3 *r2 = new TRandom3(0);
//double y = r1->Rndm();
double y = f->GetRandom();
h->Fill(x-y);
}
h->Draw();
}
