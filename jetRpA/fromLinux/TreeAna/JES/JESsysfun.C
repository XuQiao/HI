void JESsysfun(){
TF1 *f=new TF1("JESsysfun","[0] + [1]/x + [2]/x/x + [3]/x/x/x",1,1000);
f->SetParameters(1,-0.893019,-23.1071,4217.43);
//f->Draw();
}

