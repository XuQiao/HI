#include "parameter.h"
#include <iostream>
#include <iomanip>
#include "par.h"
void merge(){
	TFile *f=TFile::Open(outG,"Update");
	TString dirname;
        TVectorD Ncoll,B,Npart;
        TVectorD Ncollerr1,Berr1,Nparterr1;
        TVectorD Ncollerr2,Berr2,Nparterr2;
        TVectorD Ncollerr,Berr,Nparterr;
	TVectorD *NcollAver[nGlau+2], *NpartAver[nGlau+2], *BAver[nGlau+2];
        ofstream output("sys.txt");   output<<setprecision(2)<<fixed;
	Ncoll.ResizeTo(nGlau+2);Npart.ResizeTo(nGlau+2);B.ResizeTo(nGlau+2);
	Ncollerr1.ResizeTo(N-1);Nparterr1.ResizeTo(N-1);Berr1.ResizeTo(N-1);
	Ncollerr2.ResizeTo(N-1);Nparterr2.ResizeTo(N-1);Berr2.ResizeTo(N-1);
	Ncollerr.ResizeTo(N-1);Nparterr.ResizeTo(N-1);Berr.ResizeTo(N-1);
	if(method==0)
        output<<"kpoint_simu\t"<<"kpoint_data\t"<<"centbin\t"<<"<Ncoll>"<<"\t"<<"<Ncoll> err1"<<"\t"<<"<Ncoll> err2"<<"\t"<<"<Ncoll> err"<<"\t"<<endl;
	else
        output<<"kpoint\t"<<"centbin_data\t"<<"centbin_simu\t"<<"<Ncoll>"<<"\t"<<"<Ncoll> err1"<<"\t"<<"<Ncoll> err2"<<"\t"<<"<Ncoll> err"<<"\t"<<endl;
	//<<"<Npart>"<<"\t"<<"<Npart> err"<<"\t\t"<<"<B>"<<"\t"<<"<B> err"<<endl;
	TString name;

	for(int sth=0;sth<3;sth++){
		Ncollerr1.Zero();Nparterr1.Zero();Berr1.Zero();
		Ncollerr2.Zero();Nparterr2.Zero();Berr2.Zero();
		Ncollerr.Zero();Nparterr.Zero();Berr.Zero();
		if(sth==0){dirname = "std";}
        	else if(sth==1){dirname ="Gri055";}
	        else {dirname ="Gri101";}
        	output<<dirname<<endl;
		TVectorD *centbin = (TVectorD*)f->Get(Form("%s/G0/centbin",dirname.Data()));
		TVectorD *kpoint = (TVectorD*)f->Get(Form("%s/G0/kpoint",dirname.Data()));
		if(method==0)
		TVectorD *kpoint_ = (TVectorD*)f->Get(Form("%s/G0/kpoint_",dirname.Data()));
		else
		TVectorD *centbin_ = (TVectorD*)f->Get(Form("%s/G0/centbin_",dirname.Data()));
        	for(int i=0;i<N-1;i++){
                	for(int iGlau=0;iGlau<nGlau+2; iGlau++){
        			if(iGlau==0)
                			name = "G0";
        			else if(iGlau<nGlau && iGlau>0)
                			name = Form("Glau_%d",iGlau);
       				else
                			name = Form("bin_%d",iGlau-nGlau+1);
				NcollAver[iGlau] = (TVectorD*)f->Get(Form("%s/%s/NcollAver",dirname.Data(),name.Data()));
				NpartAver[iGlau] = (TVectorD*)f->Get(Form("%s/%s/NpartAver",dirname.Data(),name.Data()));
				BAver[iGlau] = (TVectorD*)f->Get(Form("%s/%s/BAver",dirname.Data(),name.Data()));
                       		Ncoll[iGlau]=(*NcollAver[iGlau])[i]/(*NcollAver[0])[i];
                       		Npart[iGlau]=(*NpartAver[iGlau])[i]/(*NpartAver[0])[i];
                       		B[iGlau]=(*BAver[iGlau])[i]/(*BAver[0])[i];
			}
                	for(int iGlau=0;iGlau<nGlau+2; iGlau++){
                        	if(iGlau%2==1){
                                	double Ncoll_err_add=(TMath::Abs(Ncoll[iGlau]-1)>TMath::Abs(Ncoll[iGlau+1]-1))?Ncoll[iGlau]-1:Ncoll[iGlau+1]-1;
                                	double Npart_err_add=(TMath::Abs(Npart[iGlau]-1)>TMath::Abs(Npart[iGlau+1]-1))?Npart[iGlau]-1:Npart[iGlau+1]-1;
                                	double B_err_add=(TMath::Abs(B[iGlau]-1)>TMath::Abs(B[iGlau+1]-1))?B[iGlau]-1:B[iGlau+1]-1;
                                	if(iGlau<nGlau){
						Ncollerr1[i]+=TMath::Power(Ncoll_err_add,2);
                               			Nparterr1[i]+=TMath::Power(Npart_err_add,2);
                                		Berr1[i]+=TMath::Power(B_err_add,2);
					}
					else{	
						Ncollerr2[i]+=TMath::Power(Ncoll_err_add,2);
                               			Nparterr2[i]+=TMath::Power(Npart_err_add,2);
                                		Berr2[i]+=TMath::Power(B_err_add,2);
					}
						
						Ncollerr[i]+=TMath::Power(Ncoll_err_add,2);
                               			Nparterr[i]+=TMath::Power(Npart_err_add,2);
                                		Berr[i]+=TMath::Power(B_err_add,2);
				}
			}
        		Ncollerr1[i]=TMath::Sqrt(Ncollerr1[i])*100;//(*NcollAver[0])[i];
        		Ncollerr2[i]=TMath::Sqrt(Ncollerr2[i])*100;//(*NcollAver[0])[i];
        		Nparterr1[i]=TMath::Sqrt(Nparterr1[i])*100;//(*NpartAver[0])[i];
        		Nparterr2[i]=TMath::Sqrt(Nparterr2[i])*100;//(*NpartAver[0])[i];
        		Berr1[i]=TMath::Sqrt(Berr1[i])*100;//(*BAver[0])[i];
        		Berr2[i]=TMath::Sqrt(Berr2[i])*100;//(*BAver[0])[i];
			if(method==0)
        		output<<(*kpoint)[i]<<" to "<<(*kpoint)[i+1]<<"\t"<<(*kpoint_)[i]<<" to "<<(*kpoint_)[i+1]<<"\t"<<(*centbin)[i]*100<<"% to "<<(*centbin)[i+1]*100<<"%:"<<"\t"<<(*NcollAver[0])[i]<<"\t"<<Ncollerr1[i]<<"%\t"<<Ncollerr2[i]<<"%\t"<<TMath::Sqrt(Ncollerr[i])*100<<"%"<<endl;
			else
        		output<<(*kpoint)[i]<<" to "<<(*kpoint)[i+1]<<"\t"<<(*centbin_)[i]*100<<"% to "<<(*centbin_)[i+1]*100<<"%\t"<<(*centbin)[i]*100<<"% to "<<(*centbin)[i+1]*100<<"%:"<<"\t"<<(*NcollAver[0])[i]<<"\t"<<Ncollerr1[i]<<"%\t"<<Ncollerr2[i]<<"%\t"<<TMath::Sqrt(Ncollerr[i])*100<<"%"<<endl;
        		//output<<(*kpoint)[i]<<" to "<<(*kpoint)[i+1]<<"\t"<<(*centbin)[i]*100<<"% to "<<(*centbin)[i+1]*100<<"%:"<<"\t"<<(*NcollAver[0])[i]<<"\t"<<Ncollerr1[i]<<"\t"<<Ncollerr2[i]<<"\t"<<Ncollerr[i]<<endl;
//<<(*NpartAver[0])[i]<<"\t"<<Nparterr[i]<<"\t\t"<<(*BAver[0])[i]<<"\t"<<Berr[i]<<endl;
        	}
	output<<endl;
	for(i=0;i<N-1;i++){
        	Ncollerr[i]=TMath::Sqrt(Ncollerr[i])*(*NcollAver[0])[i];
        	Nparterr[i]=TMath::Sqrt(Nparterr[i])*(*NpartAver[0])[i];
        	Berr[i]=TMath::Sqrt(Berr[i])*(*BAver[0])[i];
	}
	f->cd(dirname);
	Ncollerr.Write("Ncollerr",TObject::kOverwrite);
	Nparterr.Write("Nparterr",TObject::kOverwrite);
	Berr.Write("Berr",TObject::kOverwrite);	

	std:::vector<double> xbin(N-1);
        for(int i=0;i<N-1;i++)
		//xbin[i]=(centbin[N-i]+centbin[N-i-1])*100/2;
        	xbin[i]=0.29+i;
		double *axbin = &xbin[0];
        TGraphErrors *gr;
       		double* aNcollAver  = NcollAver[0]->GetMatrixArray();
        	double* aNcollerr  = Ncollerr.GetMatrixArray();
        	double* aNpartAver  = NpartAver[0]->GetMatrixArray();
        	double* aNparterr  = Nparterr.GetMatrixArray();
        	double* aBAver  = BAver[0]->GetMatrixArray();
        	double* aBerr = Berr.GetMatrixArray();
	TString op;
	for(int option=0; option<3; option++){
        	if(option==0){   gr = new TGraphErrors(N-1,axbin,aNcollAver,0,aNcollerr); op="Ncoll_graph";}
        	else if(option==1){      gr = new TGraphErrors(N-1,axbin,aNpartAver,0,aNparterr);	op="Npart_graph";}
        	else    {gr = new TGraphErrors(N-1,axbin,aBAver,0,aBerr);	op="B_graph";}
		gr->Write(op,TObject::kOverwrite);
	}

	}

	f->Close();
}
