void test(){
	TCanvas *c1 = new TCanvas;
	c1->Divide(3);
	c1->cd(1);
	
  	TF1 *PtDistr  = new TF1("PtDistr","exp (-(x/0.90))+0.15*exp (-(x/15))", 0.2,10);	//V~=0.06
//	TF1 *PtDistr = new TF1("PtDistr","[0]*x*TMath::Power(1+(sqrt(x*x+[1]*[1])-[1]/[2]),-[3])",0.2,10);   //Real data fit with Tsallis
     //   PtDistr->SetParameters(118.836,-0.335972,0.759243,7.5747);
        TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.2,10);      //Real data V~0.06
      //    V2vsPt->SetParameters(4.81159,1.80783,3.69272,3.11889,0.931485);
       V2vsPt->SetParameters(5,1.8,3,1.8,0.8);
  
	PtDistr->Draw();
	
        double ptmin=0.2;
        double ptmax=10.0;
       //TF1 *fint = new TF1("fint","[0]*x*TMath::Power(1+(sqrt(x*x+[1]*[1])-[1]/[2]),-[3])*((x/[4])^[5]/(1+(x/[6])^[7]))*(.00005+(1/x)^[8])",0.2,10);
        //fint->SetParameters(118.836,-0.335972,0.759243,7.5747,4.81159,1.80783,3.69272,3.11889,0.931485);
       // fint->SetParameters(118.836,-0.335972,0.759243,7.5747,4.7,1.8,3,1.8,0.8);
	TF1 *fint = new TF1("fint","(exp (-(x/0.90))+0.15*exp (-(x/15)))*((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.2,10);
         // V2vsPt->SetParameters(4.81159,1.80783,3.69272,3.11889,0.931485);
	fint->SetParameters(5,1.8,3,1.8,0.8);

	c1->cd(2);
	V2vsPt->Draw();
	c1->cd(3);
        fint->Draw();
	cout<<fint->Integral(ptmin,ptmax)/PtDistr->Integral(ptmin,ptmax)<<endl;
}
