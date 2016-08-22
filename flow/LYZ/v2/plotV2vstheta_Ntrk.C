#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
const int ntotbin=5;
const int trkpointmin[ntotbin] = {120,150,185,220,260};
const int trkpointmax[ntotbin] = {150,185,220,260,300};
const int ntheta=5;
void plotV2vstheta_Ntrk(){

const int nbin24 = 12;
const double avgtrkbin[nbin24]={44.36,54.33,68.63,88.39,108.2,131.3,162.1,196.6,227.5,247.2,269.2,301.2};
const double V24[nbin24]={0.02965,0.03913,0.04832,0.04941,0.04822,0.04955,0.049,0.04805,0.04709,0.04665,0.04772,0.04797};
const double V24err[nbin24]={0.005967,0.004262,0.001496,0.001351,0.001599,0.0003922,0.0003065,0.0002939,0.0004568,0.0008684,0.001471,0.004329};
const int startpoint = 5;


c1 = new TCanvas("c1"," ",1200,800);
makeMultiPanelCanvas(c1,3,2,0,0,0.25,0.2,0.03);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    TH1D *hFrame = new TH1D("","",300,-1,2);
    hFrame->GetXaxis()->SetTitle("#theta");
    hFrame->GetYaxis()->SetTitle("referenceV_{2}");	
    hFrame->GetXaxis()->SetTitleSize(0.04);
    hFrame->GetYaxis()->SetTitleSize(0.04);
    hFrame->GetXaxis()->SetRangeUser(-0.1,1.5);
    hFrame->SetMaximum(0.080);
    hFrame->SetMinimum(0.020);

for(int i=0;i<ntotbin;i++){
	c1->cd(i+1);
	hFrame->Draw();
//	TGraphErrors* gr24=new TGraphErrors(npt24,pt,v24[i],0,v24err[i]);
//    	TGraphErrors* grPFfull=plot("PFcandpt01to10",0,4,33,i);
//    	TGraphErrors* grPFon=plot("PFcandpt01to10/OnlyCharged",0,4,29,i);
//    	TGraphErrors* grtracklc01to10=plot("tracklcpt01to10",0,2,34,i);
    	TGraphErrors* grtracknormcpt03to6=plot("tracknormcpt03to6",0,1,21,i);
    //	TGraphErrors* grtracknormc03to6Sum=plot("tracknormcpt03to6",1,4,29,i);
    	TGraphErrors* grtracknormcpt03to3=plot("tracknormcpt03to3",0,2,34,i);
    	TGraphErrors* grtracklcpt01to10=plot("tracklcpt01to10",0,3,29,i);
    	TGraphErrors* grPFcandpt03to6=plot("PFcandpt03to6tracknormcpt03to6",0,4,24,i);
    	TGraphErrors* grPFcandpt03to3=plot("PFcandpt03to3tracknormcpt03to6",0,6,33,i);
    	TGraphErrors* grPFcandpt01to10=plot("PFcandpt01to10tracknormcpt03to6",0,7,34,i);
	//grtracklc01to10->Draw("Psame");
	//grtracknormc03to6->Draw("Psame");
	//grPFfull->Draw("Psame");
	//grPFon->Draw("Psame");
   	grtracknormcpt03to6->Draw("Psame");
   	grtracknormcpt03to3->Draw("Psame");
   	grtracklcpt01to10->Draw("Psame");
   	//grtracknormc03to3->Draw("Psame");
   	grPFcandpt03to6->Draw("Psame");
   	grPFcandpt03to3->Draw("Psame");
   	grPFcandpt01to10->Draw("Psame");
//        if(i!=ntotbin-1) gr24->Draw("Psame");
    if(i<3)       TLine *l = new TLine(hFrame->GetXaxis()->GetBinLowEdge(hFrame->GetXaxis()->GetFirst()),V24[startpoint+i],hFrame->GetXaxis()->GetBinUpEdge(hFrame->GetXaxis()->GetLast()),V24[startpoint+i]);
    else if(i==3) {  
        TLine *l = new TLine(hFrame->GetXaxis()->GetBinLowEdge(hFrame->GetXaxis()->GetFirst()),V24[startpoint+i],hFrame->GetXaxis()->GetBinUpEdge(hFrame->GetXaxis()->GetLast()),V24[startpoint+i]);
        TLine *l_ = new TLine(hFrame->GetXaxis()->GetBinLowEdge(hFrame->GetXaxis()->GetFirst()),V24[startpoint+i+1],hFrame->GetXaxis()->GetBinUpEdge(hFrame->GetXaxis()->GetLast()),V24[startpoint+i+1]);
        l_->SetLineStyle(2);
        l_->Draw("same");
    }
    else      TLine *l = new TLine(hFrame->GetXaxis()->GetBinLowEdge(hFrame->GetXaxis()->GetFirst()),V24[startpoint+i+1],hFrame->GetXaxis()->GetBinUpEdge(hFrame->GetXaxis()->GetLast()),V24[startpoint+i+1]);
    l->SetLineStyle(2);
    l->Draw("same");
        if(i==0 || i==3)
                TLatex *tlx2 = new TLatex(0.3,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
        else
                TLatex *tlx2 = new TLatex(0.1,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
        tlx2->SetNDC();
        tlx2->Draw("same");
	}
        c1->cd(ntotbin+1);
        //TLatex *tlx1 = new TLatex(0.12,0.60,Form("%.1f<p_{T}<%.1f (GeV/c)",0.1,10.0));
        TLegend *tl = new TLegend(0.2,0.40,0.5,0.85);
        tl->SetFillColor(0);
        tl->SetBorderSize(0);
        tl->SetTextSize(0.065);
      //  tl->AddEntry(gr24,"v2 4-particle cum","lp");
      //  tl->AddEntry(grPFfull,"LYZ PF candidate full","lp");
      //  tl->AddEntry(grPFon,"LYZ PF candidate charged","lp");
      //  tl->AddEntry(grtracklc01to10,"LYZ track loose cut","lp");
        // tl->AddEntry(hFrame,Form("%.1f<p_{T}<%.1f (GeV/c)",0.1,10.0),"");
        tl->AddEntry(grtracknormcpt03to6,Form("track, %.1f < p_{T} < %.1f",0.3,6.0),"lp");
        tl->AddEntry(grtracknormcpt03to3,Form("track, %.1f < p_{T} < %.1f",0.3,3.0),"lp");
        tl->AddEntry(grtracklcpt01to10,Form("track, %.1f < p_{T} < %.1f",0.1,10.0),"lp");
       //tl->AddEntry(grtracknormc03to6Sum,"LYZ general track Sum","lp");
       tl->AddEntry(grPFcandpt03to6,Form("PF cand, %.1f < p_{T} < %.1f",0.3,6.0),"lp");
       tl->AddEntry(grPFcandpt03to3,Form("PF cand, %.1f < p_{T} < %.1f",0.3,3.0),"lp");
       tl->AddEntry(grPFcandpt01to10,Form("PF cand, %.1f < p_{T} < %.1f",0.1,10.0),"lp");
  //      tl->AddEntry(hFrame,Form("%.1f<p_{T}<%.1f (GeV/c)",0.3,6.0),"");
        //tl->AddEntry(grtracknormc03to3tracknormcpt03to6,"LYZ general track 0.3 to 6","lp");
        //tl->AddEntry(grtracknormc03to3,"LYZ general track ","lp");
        //tl->AddEntry(hFrame,Form("%.1f<p_{T}<%.1f (GeV/c)",0.3,3.0),"");
        //tlx1->SetNDC();
        //tlx1->SetTextSize(0.065);
        tlx2->SetTextSize(0.065);
        hFrame->Draw();
        //tlx1->Draw("same");
        tl->Draw("same");
        c1->Print("V2vstheta_Ntrk.png");
}

TGraphErrors* plot(TString dir, int isSum, int color, int marker, int i){

    if(isSum==0) TFile *f = TFile::Open(Form("%s/M%d%d/mergedV_Prod.root",dir.Data(),trkpointmax[i],trkpointmin[i]));
    else TFile *f = TFile::Open(Form("%s/M%d%d/mergedV_Sum.root",dir.Data(),trkpointmax[i],trkpointmin[i]));
    int xbin=0;
    int xpt=0;
    const int nn = 2;
	TVectorD *vecNtrk = (TVectorD*)f->Get("avgtrk");
        TVectorD *vecavgmult = (TVectorD*)f->Get("avgmultall");
        TVectorD *vecV2mean_0=(TVectorD*)f->Get(Form("D_%d/Vmean",xbin));
        TVectorD *vecV2_0=(TVectorD*)f->Get(Form("D_%d/D_%d/V",xbin,xpt));
        TVectorD *vecV2err_0=(TVectorD*)f->Get(Form("D_%d/D_%d/deltaV",xbin,xpt));
        double *V2_0=vecV2_0->GetMatrixArray();
        double *V2err_0=vecV2err_0->GetMatrixArray();
        double theta[ntheta];
        for(int itheta=0;itheta<ntheta;itheta++){
            theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
        }

        gr=new TGraphErrors(ntheta,theta,V2_0,0,V2err_0);
	
	gr->SetMarkerStyle(marker);
	gr->SetMarkerSize(1.6);
	gr->SetMarkerColor(color);
	gr->SetLineColor(color);
	f->Close();
	return gr;
}

