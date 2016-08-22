#include "data_cms.C"
#include "data_atlas.C"

void drawSys(TGraph *h, double *sys_minus, double *sys_plus, double *width_minus, double *width_plus, int theColor, int fillStyle = -1, int lineStyle = -1, int binmin, int binnum)
{
    for (int i=binmin;i<binnum+binmin;i++)
    {
        double val;
        double theX;
        h->GetPoint(i,theX,val);
//        double err = val * sys[i]/100;
        TBox *b = new TBox(theX-width_minus[i],val-sys_minus[i],theX+width_plus[i],val+sys_plus[i]);
        
        b->SetLineColor(theColor);
        b->SetFillColor(theColor);
        if ( fillStyle > -1 ) b->SetFillStyle(fillStyle);
        if ( lineStyle > -1 ) b->SetLineStyle(lineStyle);
        
        b->Draw();
    }
}

void PlotAll(double *p8387_yval,double *p8387_xval,double *p8387_xerrminus,double *p8387_xerrplus,double *p8387_yerrminus,double *p8387_yerrplus,double *p8387_ystatminus,double *p8387_ystatplus,int p8387_numpoints,double *p8075_yval,double *p8075_xval,double *p8075_xerrminus,double *p8075_xerrplus,double *p8075_yerrminus,double *p8075_yerrplus,double *p8075_ystatminus,double *p8075_ystatplus,int p8075_numpoints,const int intersect,const int offset_1, const int offset_2, int xmin,int xmax,double ymin,int ymax,char *yrange,char *dirname, ofstream &fout){

double ratio_yval[intersect];
double ratio_xval[intersect];
double ratio_xerrzero[intersect];
double ratio_xerrminus[intersect];
double ratio_xerrplus[intersect];
double ratio_ystatminus[intersect];
double ratio_ystatplus[intersect];
double ratio_yerrminus[intersect];
double ratio_yerrplus[intersect];
double ratio_ystatminus_corr[intersect];
double ratio_ystatplus_corr[intersect];
double ratio_yerrminus_corr[intersect];
double ratio_yerrplus_corr[intersect];
double ratio_yval_corr[intersect];

Char_t buf[1024],textcontent[1024];
sprintf(textcontent,"|y| range is %s",yrange);
sprintf(buf,"%s.root",dirname);
TFile *f = new TFile(buf,"Update");

TCanvas *c1 = new TCanvas("c1");
c1->Divide(2);
c1->cd(1)->SetLogy();
c1->cd(1)->SetLogx();

// Plot: p8387
TGraphAsymmErrors *p8387 = new TGraphAsymmErrors(p8387_numpoints, p8387_xval, p8387_yval, p8387_xerrminus, p8387_xerrplus, p8387_ystatminus, p8387_ystatplus);
p8387->SetMinimum(ymin);
f->cd();
sprintf(buf,"R2_%s",yrange);
p8387->Write(buf,TObject::kOverwrite);
p8387->SetMaximum(ymax);
p8387->SetMarkerStyle(24);
p8387->SetMarkerSize(0.8); 
p8387->SetMarkerColor(kBlack);
//p8387->SetTitle("Jet pt distribution, R=0.7");
p8387->Draw("AP");
drawSys(p8387,p8387_yerrminus,p8387_yerrplus,p8387_xerrminus, p8387_xerrplus,kBlue,3001,-1,0,p8387->GetN());
sprintf(buf,"plots/%s/R2_%s.png",dirname,yrange);
c1->cd(1)->Print(buf);

c1->cd(2)->SetLogy();
c1->cd(2)->SetLogx();
// Plot: p8387
TGraphAsymmErrors *p8075 = new TGraphAsymmErrors(p8075_numpoints, p8075_xval, p8075_yval, p8075_xerrminus, p8075_xerrplus, p8075_ystatminus, p8075_ystatplus);
f->cd();
sprintf(buf,"R1_%s",yrange);
p8075->Write(buf,TObject::kOverwrite);
p8075->SetMinimum(ymin);
p8075->SetMaximum(ymax);
p8075->SetMarkerStyle(20);
p8075->SetMarkerSize(0.8);
p8075->SetMarkerColor(kRed);
//p8075->SetTitle("Jet pt distribution, R=0.5");
p8075->Draw("AP");
drawSys(p8075,p8075_yerrminus,p8075_yerrplus,p8075_xerrminus, p8075_xerrplus,kRed,3004,-1,0,p8075->GetN());
sprintf(buf,"plots/%s/R1_%s.png",dirname,yrange);
c1->cd(2)->Print(buf);

TCanvas *c2= new TCanvas("c2");
c2->cd();
c2->SetLogy();
c2->SetLogx();
p8387->SetMinimum(ymin);
p8387->SetMaximum(ymax);
p8075->SetMinimum(ymin);
p8075->SetMaximum(ymax);
p8387->GetXaxis()->SetLimits(xmin,xmax);
p8387->GetXaxis()->SetMoreLogLabels(kTRUE);
//p8387->SetTitle("comparsion between CMS different jet cone size yield");
p8387->GetXaxis()->SetTitle("p_{T} in GeV");
p8387->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}dy (pb/GeV)");
p8387->SetMarkerSize(1.3);
p8075->SetMarkerStyle(21);
p8387->Draw("AP");
p8075->Draw("P");
//TLegend *leg= new TLegend(0.6,0.6,0.8,0.8);
//leg->AddEntry(p8387,"Cone size R=0.7","lp");
//leg->AddEntry(p8075,"Cone size R=0.5","lp");
//leg->Draw("same");

drawSys(p8387,p8387_yerrminus,p8387_yerrplus,p8387_xerrminus, p8387_xerrplus,kBlue,3001,-1,offset_2,intersect);
drawSys(p8075,p8075_yerrminus,p8075_yerrplus,p8075_xerrminus, p8075_xerrplus,kRed,3004,-1,offset_1,intersect);
sprintf(buf,"plots/%s/R_%s.png",dirname,yrange);
c2->Print(buf);

//Ratio plot
for(int i=0;i<intersect;i++){
ratio_yval[i]=p8075_yval[i+offset_1]/p8387_yval[i+offset_2];
ratio_xval[i]=(p8075_xval[i+offset_1]+p8387_xval[i+offset_2])/2;
ratio_xerrzero[i]=0;
ratio_xerrminus[i]=(p8387_xerrminus[i+offset_2]+p8075_xerrminus[i+offset_1])/2;
ratio_xerrplus[i]=(p8387_xerrplus[i+offset_2]+p8075_xerrplus[i+offset_1])/2;
ratio_ystatminus[i]=TMath::Sqrt(TMath::Power(p8387_ystatminus[i+offset_2]/p8387_yval[i+offset_2],2)+TMath::Power(p8075_ystatminus[i+offset_1]/p8075_yval[i+offset_1],2))*ratio_yval[i];
ratio_ystatplus[i]=TMath::Sqrt(TMath::Power(p8387_ystatplus[i+offset_2]/p8387_yval[i+offset_2],2)+TMath::Power(p8075_ystatplus[i+offset_1]/p8075_yval[i+offset_1],2))*ratio_yval[i];
ratio_yerrminus[i]=TMath::Sqrt(TMath::Power(p8387_yerrminus[i+offset_2]/p8387_yval[i+offset_2],2)+TMath::Power(p8075_yerrminus[i+offset_1]/p8075_yval[i+offset_1],2))*ratio_yval[i];
ratio_yerrplus[i]=TMath::Sqrt(TMath::Power(p8387_yerrplus[i+offset_2]/p8387_yval[i+offset_2],2)+TMath::Power(p8075_yerrplus[i+offset_1]/p8075_yval[i+offset_1],2))*ratio_yval[i];
}

TGraphAsymmErrors *Ratio = new TGraphAsymmErrors(intersect,ratio_xval,ratio_yval,ratio_xerrminus,ratio_xerrplus,ratio_ystatminus,ratio_ystatplus);
TCanvas *c3 = new  TCanvas("c3");
c3->cd();
Ratio->SetMaximum(2.5);
Ratio->SetMinimum(0);
Ratio->SetMarkerStyle(20);
Ratio->SetMarkerSize(0.8);
Ratio->SetTitle("Ratio plot");
Ratio->GetXaxis()->SetTitle("p_{T} in GeV");
//Ratio->GetYaxis()->SetTitle("#frac{R=0.5}{R=0.7}");
Ratio->Draw("AP");

//fit and do Ratio_corr plot
Ratio->Fit("pol0","","",xmin,xmax);
f->cd();
sprintf(buf,"Ratio_%s",yrange);
Ratio->Write(buf,TObject::kOverwrite);
TF1 *fConstant = Ratio->GetFunction("pol0");
double scalefactor=fConstant->GetParameter(0);
drawSys(Ratio,ratio_yerrminus,ratio_yerrplus,ratio_xerrminus,ratio_xerrplus,kBlue,3003,-1,0,intersect);
sprintf(buf,"plots/%s/ratio_%s.png",dirname,yrange);
c3->Print(buf);

fout<<dirname<<": "<<"|y| range is: "<<yrange<<", scale factor is:"<<scalefactor<<endl;

for(int i=0;i<intersect;i++){
ratio_yval_corr[i]=ratio_yval[i]/scalefactor;
ratio_ystatminus_corr[i]=ratio_ystatminus[i]/scalefactor;
ratio_ystatplus_corr[i]=ratio_ystatplus[i]/scalefactor;
ratio_yerrminus_corr[i]=ratio_yerrminus[i]/scalefactor;
ratio_yerrplus_corr[i]=ratio_yerrplus[i]/scalefactor;
}

TGraphAsymmErrors *Ratio_corr = new TGraphAsymmErrors(intersect,ratio_xval,ratio_yval_corr,ratio_xerrminus,ratio_xerrplus,ratio_ystatminus_corr,ratio_ystatplus_corr);
f->cd();
sprintf(buf,"Ratio_corr_%s",yrange);
Ratio_corr->Write(buf,TObject::kOverwrite);
TCanvas *c4 = new TCanvas("c4");
c4->cd();
Ratio_corr->SetMaximum(3);
Ratio_corr->SetMinimum(0);
Ratio_corr->SetMarkerStyle(20);
Ratio_corr->SetMarkerSize(0.8);
Ratio_corr->SetTitle("Ratio plot Scaled");
Ratio_corr->GetXaxis()->SetTitle("p_{T} in GeV");
//Ratio_corr->GetYaxis()->SetTitle("#frac{R=0.5}{R=0.7*scalefactor}");
Ratio_corr->Draw("AP");
drawSys(Ratio_corr,ratio_yerrminus_corr,ratio_yerrplus_corr,ratio_xerrminus, ratio_xerrplus,kBlue,3003,-1,0,intersect);
sprintf(buf,"plots/%s/ratio_corr_%s.png",dirname,yrange);
c4->Print(buf);

//compare after scale factor correction
const int p8387_N=32;
double p8387_yval_corr[p8387_N];
double p8387_ystatminus_corr[p8387_N];
double p8387_ystatplus_corr[p8387_N];
double p8387_yerrminus_corr[p8387_N];
double p8387_yerrplus_corr[p8387_N];

for(i=0;i<p8387_N;i++){
p8387_yval_corr[i]=p8387_yval[i]*scalefactor;
p8387_ystatminus_corr[i]=p8387_ystatminus[i]*scalefactor;
p8387_ystatplus_corr[i]=p8387_ystatplus[i]*scalefactor;
p8387_ystatminus_corr[i]=p8387_ystatminus[i]*scalefactor;
p8387_ystatplus_corr[i]=p8387_ystatplus[i]*scalefactor;
}

TGraphAsymmErrors *p8387_corr = new TGraphAsymmErrors(intersect, p8387_xval, p8387_yval_corr, p8387_xerrminus, p8387_xerrplus, p8387_ystatminus_corr, p8387_ystatplus_corr);
f->cd();
sprintf(buf,"R2_corr_%s",yrange);
p8387_corr->Write(buf,TObject::kOverwrite);
TCanvas *c5= new TCanvas("c5");
c5->cd();
c5->SetLogy();
c5->SetLogx();
p8387_corr->SetMinimum(ymin);
p8387_corr->SetMaximum(ymax);
p8387_corr->GetXaxis()->SetLimits(xmin,xmax);
p8387_corr->GetXaxis()->SetMoreLogLabels(kTRUE);
//p8387_corr->SetTitle("comparsion between CMS different jet cone size yield AFTER SCALE");
p8387_corr->GetXaxis()->SetTitle("p_{T} in GeV");
p8387_corr->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}dy (pb/GeV)");
p8387_corr->SetMarkerSize(1.3);
p8387_corr->SetMarkerStyle(24);
p8387_corr->Draw("AP");
p8075->Draw("P");
//TLegend *leg= new TLegend(0.6,0.6,0.8,0.8);
//leg->AddEntry(p8387_corr,"Cone size R=0.7, L=5.0 fb^{-1}","lp");
//leg->AddEntry(p8075,"Cone size R=0.5 L=34 pb^{-1}","lp");
//leg->Draw("same");

/*
TLatex T;
T.SetTextAlign(12);
T.DrawLatex(0.8,0.8,textcontent);
*/
drawSys(p8387_corr,p8387_yerrminus_corr,p8387_yerrplus_corr,p8387_xerrminus, p8387_xerrplus,kBlue,3001,-1,offset_2,intersect);
drawSys(p8075,p8075_yerrminus,p8075_yerrplus,p8075_xerrminus, p8075_xerrplus,kRed, 3004,-1,offset_1,intersect);
sprintf(buf,"plots/%s/R_corr_%s.png",dirname,yrange);
c5->Print(buf);

}

void jet_pt_check(){

Char_t *outtxt="scale_results.txt";
ofstream fout;
fout.open(outtxt);

//-----------------------------------CMS Jet Spectrum------------------------------------------------

PlotAll(p8387_d1x1y1_yval,p8387_d1x1y1_xval,p8387_d1x1y1_xerrminus,p8387_d1x1y1_xerrplus,p8387_d1x1y1_yerrminus,p8387_d1x1y1_yerrplus,p8387_d1x1y1_ystatminus,p8387_d1x1y1_ystatplus,p8387_d1x1y1_numpoints,p8075_d1x1y1_yval,p8075_d1x1y1_xval,p8075_d1x1y1_xerrminus,p8075_d1x1y1_xerrplus,p8075_d1x1y1_yerrminus,p8075_d1x1y1_yerrplus,p8075_d1x1y1_ystatminus,p8075_d1x1y1_ystatplus,p8075_d1x1y1_numpoints, intersect_05,offset_1_05,offset_2_05,xmin_05,xmax_05,ymin_05,ymax_05,yrange_05,"cms",fout);

//PlotAll(p8387_d2x1y1_yval,p8387_d2x1y1_xval,p8387_d2x1y1_xerrminus,p8387_d2x1y1_xerrplus,p8387_d2x1y1_yerrminus,p8387_d2x1y1_yerrplus,p8387_d2x1y1_ystatminus,p8387_d2x1y1_ystatplus, p8387_d2x1y1_numpoints,p8075_d2x1y1_yval,p8075_d2x1y1_xval,p8075_d2x1y1_xerrminus,p8075_d2x1y1_xerrplus,p8075_d2x1y1_yerrminus,p8075_d2x1y1_yerrplus,p8075_d2x1y1_ystatminus,p8075_d2x1y1_ystatplus, p8075_d2x1y1_numpoints, intersect_0510,offset_0510,xmin_0510,xmax_0510,ymin_0510,ymax_0510,yrange_0510,"cms",fout);

//PlotAll(p8387_d3x1y1_yval,p8387_d3x1y1_xval,p8387_d3x1y1_xerrminus,p8387_d3x1y1_xerrplus,p8387_d3x1y1_yerrminus,p8387_d3x1y1_yerrplus,p8387_d3x1y1_ystatminus,p8387_d3x1y1_ystatplus,p8387_d3x1y1_numpoints,p8075_d3x1y1_yval,p8075_d3x1y1_xval,p8075_d3x1y1_xerrminus,p8075_d3x1y1_xerrplus,p8075_d3x1y1_yerrminus,p8075_d3x1y1_yerrplus,p8075_d3x1y1_ystatminus,p8075_d3x1y1_ystatplus,p8075_d3x1y1_numpoints,intersect_1015,offset_1015,xmin_1015,xmax_1015,ymin_1015,ymax_1015,yrange_1015,"cms",fout);

//PlotAll(p8387_d4x1y1_yval,p8387_d4x1y1_xval,p8387_d4x1y1_xerrminus,p8387_d4x1y1_xerrplus,p8387_d4x1y1_yerrminus,p8387_d4x1y1_yerrplus,p8387_d4x1y1_ystatminus,p8387_d4x1y1_ystatplus,p8387_d4x1y1_numpoints,p8075_d4x1y1_yval,p8075_d4x1y1_xval,p8075_d4x1y1_xerrminus,p8075_d4x1y1_xerrplus,p8075_d4x1y1_yerrminus,p8075_d4x1y1_yerrplus,p8075_d4x1y1_ystatminus,p8075_d4x1y1_ystatplus,p8075_d4x1y1_numpoints,intersect_1520,offset_1520,xmin_1520,xmax_1520,ymin_1520,ymax_1520,yrange_1520,"cms",fout);

//PlotAll(p8387_d5x1y1_yval,p8387_d5x1y1_xval,p8387_d5x1y1_xerrminus,p8387_d5x1y1_xerrplus,p8387_d5x1y1_yerrminus,p8387_d5x1y1_yerrplus,p8387_d5x1y1_ystatminus,p8387_d5x1y1_ystatplus,p8387_d5x1y1_numpoints,p8075_d5x1y1_yval,p8075_d5x1y1_xval,p8075_d5x1y1_xerrminus,p8075_d5x1y1_xerrplus,p8075_d5x1y1_yerrminus,p8075_d5x1y1_yerrplus,p8075_d5x1y1_ystatminus,p8075_d5x1y1_ystatplus,p8075_d5x1y1_numpoints,intersect_2025,offset_2025,xmin_2025,xmax_2025,ymin_2025,ymax_2025,yrange_2025,"cms",fout);

//-----------------------------------ATLAS Full Jet Spectrum-------------------------------------------------

PlotAll(p8137_d8x1y1_yval,p8137_d8x1y1_xval,p8137_d8x1y1_xerrminus,p8137_d8x1y1_xerrplus,p8137_d8x1y1_yerrminus,p8137_d8x1y1_yerrplus,p8137_d8x1y1_ystatminus,p8137_d8x1y1_ystatplus,p8137_d8x1y1_numpoints,p8137_d1x1y1_yval,p8137_d1x1y1_xval,p8137_d1x1y1_xerrminus,p8137_d1x1y1_xerrplus,p8137_d1x1y1_yerrminus,p8137_d1x1y1_yerrplus,p8137_d1x1y1_ystatminus,p8137_d1x1y1_ystatplus,p8137_d1x1y1_numpoints, intersect_03,offset_1_03,offset_2_03,xmin_03,xmax_03,ymin_03,ymax_03,yrange_03,"atlas_fulljet",fout);

//PlotAll(p8137_d9x1y1_yval,p8137_d8x1y1_xval,p8137_d8x1y1_xerrminus,p8137_d8x1y1_xerrplus,p8137_d8x1y1_yerrminus,p8137_d8x1y1_yerrplus,p8137_d8x1y1_ystatminus,p8137_d8x1y1_ystatplus,p8137_d8x1y1_numpoints,p8137_d2x1y1_yval,p8137_d2x1y1_xval,p8137_d2x1y1_xerrminus,p8137_d2x1y1_xerrplus,p8137_d2x1y1_yerrminus,p8137_d2x1y1_yerrplus,p8137_d2x1y1_ystatminus,p8137_d2x1y1_ystatplus,p8137_d2x1y1_numpoints, intersect_0308,offset_1_0308,offset_2_0308,xmin_0308,xmax_0308,ymin_0308,ymax_0308,yrange_0308,"atlas_fulljet",fout);

//PlotAll(p8137_d10x1y1_yval,p8137_d10x1y1_xval,p8137_d10x1y1_xerrminus,p8137_d10x1y1_xerrplus,p8137_d10x1y1_yerrminus,p8137_d10x1y1_yerrplus,p8137_d10x1y1_ystatminus,p8137_d10x1y1_ystatplus,p8137_d10x1y1_numpoints,p8137_d3x1y1_yval,p8137_d3x1y1_xval,p8137_d3x1y1_xerrminus,p8137_d3x1y1_xerrplus,p8137_d3x1y1_yerrminus,p8137_d3x1y1_yerrplus,p8137_d3x1y1_ystatminus,p8137_d3x1y1_ystatplus,p8137_d3x1y1_numpoints, intersect,offset,xmin_0812,xmax_0812,ymin_0812,ymax_0812,yrange_0812,"atlas_fulljet",fout);

//PlotAll(p8137_d11x1y1_yval,p8137_d11x1y1_xval,p8137_d11x1y1_xerrminus,p8137_d11x1y1_xerrplus,p8137_d11x1y1_yerrminus,p8137_d11x1y1_yerrplus,p8137_d11x1y1_ystatminus,p8137_d11x1y1_ystatplus,p8137_d11x1y1_numpoints,p8137_d4x1y1_yval,p8137_d4x1y1_xval,p8137_d4x1y1_xerrminus,p8137_d4x1y1_xerrplus,p8137_d4x1y1_yerrminus,p8137_d4x1y1_yerrplus,p8137_d4x1y1_ystatminus,p8137_d4x1y1_ystatplus,p8137_d4x1y1_numpoints, intersect,offset,xmin_1221,xmax_1221,ymin_1221,ymax_1221,yrange_1221,"atlas_fulljet",fout);

//PlotAll(p8137_d12x1y1_yval,p8137_d12x1y1_xval,p8137_d12x1y1_xerrminus,p8137_d12x1y1_xerrplus,p8137_d12x1y1_yerrminus,p8137_d12x1y1_yerrplus,p8137_d12x1y1_ystatminus,p8137_d12x1y1_ystatplus,p8137_d12x1y1_numpoints,p8137_d5x1y1_yval,p8137_d5x1y1_xval,p8137_d5x1y1_xerrminus,p8137_d5x1y1_xerrplus,p8137_d5x1y1_yerrminus,p8137_d5x1y1_yerrplus,p8137_d5x1y1_ystatminus,p8137_d5x1y1_ystatplus,p8137_d5x1y1_numpoints, intersect,offset,xmin_2128,xmax_2128,ymin_2128,ymax_2128,yrange_2128,"atlas_fulljet",fout);

//PlotAll(p8137_d13x1y1_yval,p8137_d13x1y1_xval,p8137_d13x1y1_xerrminus,p8137_d13x1y1_xerrplus,p8137_d13x1y1_yerrminus,p8137_d13x1y1_yerrplus,p8137_d13x1y1_ystatminus,p8137_d13x1y1_ystatplus,p8137_d13x1y1_numpoints,p8137_d6x1y1_yval,p8137_d6x1y1_xval,p8137_d6x1y1_xerrminus,p8137_d6x1y1_xerrplus,p8137_d6x1y1_yerrminus,p8137_d6x1y1_yerrplus,p8137_d6x1y1_ystatminus,p8137_d6x1y1_ystatplus,p8137_d6x1y1_numpoints, intersect,offset,xmin_2836,xmax_2836,ymin_2836,ymax_2836,yrange_2836,"atlas_fulljet",fout);

//PlotAll(p8137_d14x1y1_yval,p8137_d14x1y1_xval,p8137_d14x1y1_xerrminus,p8137_d14x1y1_xerrplus,p8137_d14x1y1_yerrminus,p8137_d14x1y1_yerrplus,p8137_d14x1y1_ystatminus,p8137_d14x1y1_ystatplus,p8137_d14x1y1_numpoints,p8137_d7x1y1_yval,p8137_d7x1y1_xval,p8137_d7x1y1_xerrminus,p8137_d7x1y1_xerrplus,p8137_d7x1y1_yerrminus,p8137_d7x1y1_yerrplus,p8137_d7x1y1_ystatminus,p8137_d7x1y1_ystatplus,p8137_d7x1y1_numpoints, intersect,offset,xmin_3644,xmax_3644,ymin_3644,ymax_3644,yrange_3644,"atlas_fulljet",fout);

//----------------------------------ATLAS Charged Spetrum----------------------------------------------------------------------

PlotAll(p8062_d1x1y1_yval,p8062_d1x1y1_xval,p8062_d1x1y1_xerrminus,p8062_d1x1y1_xerrplus,p8062_d1x1y1_yerrminus,p8062_d1x1y1_yerrplus,p8062_d1x1y1_ystatminus,p8062_d1x1y1_ystatplus,p8062_d1x1y1_numpoints,p8062_d1x1y2_yval,p8062_d1x1y2_xval,p8062_d1x1y2_xerrminus,p8062_d1x1y2_xerrplus,p8062_d1x1y2_yerrminus,p8062_d1x1y2_yerrplus,p8062_d1x1y2_ystatminus,p8062_d1x1y2_ystatplus,p8062_d1x1y2_numpoints, intersect_05_at,offset_1_05_at,offset_2_05_at,xmin_05_at,xmax_05_at,ymin_05_at,ymax_05_at,yrange_05,"atlas_chargedjet",fout);
}
