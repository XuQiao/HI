
//#include <memory>
//void ellip()
{
  int i,j;
  double pt22[30];
  double v22[30];
  double err22[30];
  double test;
  double pt24[28];
  double v24[28];
  double err24[28];
  double pt32[30];
  double v32[30];
  double err32[30];
  double pt34[30];
  double v34[30];
  double err34[30];
  double pt43[30];
  double v43[30];
  double err43[30];
  double pt45[30];
  double v45[30];
  double err45[30];
  double pt52[30];
  double v52[30];
  double err52[30];
  double pt54[30];
  double v54[30];
  double err54[30];
  double pt63[30];
  double v63[30];
  double err63[30];
  double pt65[30];
  double v65[30];
  double err65[30];


  double ratio[30];
  double ratioerr[30];
  double M2[12];
  double X2[12];
  double temp;

  const int Nbin=13;
  ifstream inR22;
  ifstream inR24;
  ifstream inR32;
  ifstream inR34;
  ifstream inR43;
  ifstream inR45;
  ifstream inR52;
  ifstream inR54;
  ifstream inR63;
  ifstream inR65;
  ifstream inG;
  ifstream inRc;
  ifstream inM;
  ofstream outX;
  ifstream tests;
  TF1 *V1vsEta = new TF1("V1vsEta","-exp(-(x+1)^2)/20-x/30+exp(-(x-1)^2)/20",-2.4,2.4);
  TF1 *V2vsPt = new TF1("V2vsPt","((x/3.4)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);
  TF1 *V3vsPt = new TF1("V3vsPt","((x/3.2)^2.3/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.2,10);
  TF1 *V4vsPt = new TF1("V4vsPt","((x/4.8)^2.1/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.2,10);
  TF1 *V5vsPt = new TF1("V5vsPt","((x/6.0)^3.2/(1+(x/11.4)^2.1))*(.00005+(1/x)^1.4)",0.2,10);
  TF1 *V6vsPt = new TF1("V6vsPt","((x/5.6)^2.4/(1+(x/4.7)^2.1))*(.00005+(1/x)^1.4)",0.2,10);

  inR22.open("../Sum/v2/v2.dat");  if(!inR22.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR22>>pt22[i];  //if(R[i]>0.5||R[i]<0.)R[i]=0.9;
    inR22>>v22[i];
    inR22>>err22[i];
    cout<<pt22[i]<<"  "<<v22[i]<<"  "<<err22[i]<<endl;
  } 
  inR24.open("../Prod/v2/v2.dat");  if(!inR24.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR24>>pt24[i];  
    inR24>>v24[i];
    inR24>>err24[i];
  }
  inR32.open("../Sum/v3/v2.dat");  if(!inR32.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR32>>pt32[i];
    inR32>>v32[i];
    inR32>>err32[i];
  }
  inR34.open("../Prod/v3/v2.dat");  if(!inR34.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR34>>pt34[i];
    inR34>>v34[i];
    inR34>>err34[i];
  }
  inR43.open("../Sum/v2/v4.dat");  if(!inR43.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR43>>pt43[i];
    inR43>>v43[i];
    inR43>>err43[i];
  }
  inR45.open("../Prod/v2/v4.dat");  if(!inR45.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR45>>pt45[i];
    inR45>>v45[i];
    inR45>>err45[i];
  }
  inR52.open("../Sum/v2/v6.dat");  if(!inR52.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR52>>pt52[i];
    inR52>>v52[i];
    inR52>>err52[i];
  }
  inR54.open("../Prod/v2/v6.dat");  if(!inR54.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR54>>pt54[i];
    inR54>>v54[i];
    inR54>>err54[i];
  }
  inR63.open("../Sum/v3/v4.dat");  if(!inR63.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR63>>pt63[i];
    inR63>>v63[i];
    inR63>>err63[i];
  }
  inR65.open("../Prod/v3/v4.dat");  if(!inR65.good())   cout<<" input fail"<<endl;
  else cout<<" have put in"<<endl;
  for(i=0;i<Nbin;i++){
    inR65>>pt65[i];
    inR65>>v65[i];
    inR65>>err65[i];
  }


  TCanvas *c1 = new TCanvas("c1","c1",1,1,800,600);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderSize(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);
  c1->SetTopMargin(0.03);
  c1->SetRightMargin(0.03);
  c1->SetTicks(-1);
   //c1->SetGrid();
  TGraphErrors *graph = new TGraphErrors(Nbin,pt22,v22,0,err22);
  graph->SetTitle("");
  graph->GetXaxis()->SetRangeUser(0.2,10);
  graph->GetYaxis()->SetRangeUser(0.0,0.219);
  graph->GetYaxis()->SetTitle("v_{mn}{LYZ}");
  graph->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graph->GetYaxis()->SetTitleOffset(.92);
    graph->GetXaxis()->SetTitleOffset(.8);
    graph->GetYaxis()->SetTitleSize(0.06);
    graph->GetXaxis()->SetTitleSize(0.06);
    graph->GetYaxis()->SetTitleFont(42);
    graph->GetXaxis()->SetTitleFont(42);
    graph->GetYaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->CenterTitle(kTRUE);
    graph->GetXaxis()->CenterTitle(kTRUE);
  graph->SetMarkerStyle(24);
  graph->SetMarkerColor(2);
  graph->SetMarkerSize(1.5);
  graph->Draw("Ap");
  V2vsPt->SetLineStyle(2);
  V3vsPt->SetLineStyle(3);
  V4vsPt->SetLineStyle(4);
  V5vsPt->SetLineStyle(6);
  V6vsPt->SetLineStyle(8);
  V2vsPt->Draw("same");
  V3vsPt->Draw("same");
  V4vsPt->Draw("same");
  //V5vsPt->Draw("same");
  V6vsPt->Draw("same");

  TGraphErrors *g22 = new TGraphErrors(Nbin,pt22,v22,0,err22);
  g22->SetMarkerStyle(21);
  g22->SetMarkerColor(2);
  g22->SetMarkerSize(1.5);
  //g22->Draw("psame");
  TGraphErrors *g24 = new TGraphErrors(Nbin,pt24,v24,0,err24);
  g24->SetMarkerStyle(20);
  g24->SetMarkerColor(1);
  g24->SetMarkerSize(1.5);
  g24->Draw("psame");
  TGraphErrors *g32 = new TGraphErrors(Nbin,pt32,v32,0,err32);
  g32->SetMarkerStyle(25);
  g32->SetMarkerColor(2);
  g32->SetMarkerSize(1.5);
  g32->Draw("psame");
  TGraphErrors *g34 = new TGraphErrors(Nbin,pt34,v34,0,err34);
  g34->SetMarkerStyle(21);
  g34->SetMarkerColor(1);
  g34->SetMarkerSize(1.5);
  g34->Draw("psame");
  TGraphErrors *g43 = new TGraphErrors(Nbin,pt43,v43,0,err43);
  g43->SetMarkerStyle(26);
  g43->SetMarkerColor(2);
  g43->SetMarkerSize(1.5);
  g43->Draw("psame");
  TGraphErrors *g45 = new TGraphErrors(Nbin,pt45,v45,0,err45);
  g45->SetMarkerStyle(22);
  g45->SetMarkerColor(1);
  g45->SetMarkerSize(1.5);
  g45->Draw("psame");
  TGraphErrors *g52 = new TGraphErrors(Nbin,pt52,v52,0,err52);
  g52->SetMarkerStyle(27);
  g52->SetMarkerColor(2);
  g52->SetMarkerSize(1.5);
  g52->Draw("psame");
  TGraphErrors *g54 = new TGraphErrors(Nbin,pt54,v54,0,err54);
  g54->SetMarkerStyle(23);
  g54->SetMarkerColor(1);
  g54->SetMarkerSize(1.5);
  g54->Draw("psame");
  TGraphErrors *g63 = new TGraphErrors(Nbin,pt63,v63,0,err63);
  g63->SetMarkerStyle(30);
  g63->SetMarkerColor(2);
  g63->SetMarkerSize(1.5);
  g63->Draw("psame");
  TGraphErrors *g65 = new TGraphErrors(Nbin,pt65,v65,0,err65);
  g65->SetMarkerStyle(29);
  g65->SetMarkerColor(1);
  g65->SetMarkerSize(1.5);
  g65->Draw("psame");
double ymin=0.75;
double ymax=0.95;
    TLegend *leg0 = new TLegend(0.14,ymin-0.02,0.3,ymax);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    leg0->SetTextColor(1);
    leg0->SetTextSize(0.048);
    //leg0->SetLineStyle(0.06);
    leg0->AddEntry(V2vsPt,"input v_{2}","l");
    leg0->AddEntry(V3vsPt,"input v_{3}","l");
    leg0->AddEntry(V4vsPt,"input v_{4}","l");
    //leg0->AddEntry(V5vsPt,"input v_{6}","l");
    leg0->AddEntry(V6vsPt,"input v_{6}","l");
    leg0->Draw();
    TLegend *leg0 = new TLegend(0.53,ymin,0.63,ymax);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    leg0->SetTextColor(1);
    leg0->SetTextSize(0.048);
    //leg0->SetLineStyle(0.06);
    leg0->AddEntry(graph," v_{2=12}{Sum}","p");
    leg0->AddEntry(g32," v_{3=13}{Sum}","p");
    leg0->AddEntry(g43," v_{4=22}{Sum}","p");
    leg0->AddEntry(g52," v_{6=32}{Sum}","p");
    leg0->AddEntry(g63," v_{6=23}{Sum}","p");
    leg0->Draw();
    TLegend *leg0 = new TLegend(0.74,ymin,0.84,ymax);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    leg0->SetTextColor(1);
    leg0->SetTextSize(0.048);
    //leg0->SetLineStyle(0.06);
    leg0->AddEntry(g24," v_{2=12}{Prod}","p");
    leg0->AddEntry(g34," v_{3=13}{Prod}","p");
    leg0->AddEntry(g45," v_{4=22}{Prod}","p");
    leg0->AddEntry(g54," v_{6=32}{Prod}","p");
    leg0->AddEntry(g65," v_{6=23}{Prod}","p");
    leg0->Draw();


   TLatex *tex = new TLatex(2.5.,0.2,"STEG,  50k");
   tex->SetTextColor(1);
   tex->SetTextFont(42);
   tex->SetTextSize(0.045);
   tex->Draw();
   TLatex *tex = new TLatex(7.3,0.24,"Cumulant method");
   tex->SetTextColor(1);
   tex->SetTextFont(42);
   tex->SetTextSize(0.05);
   //tex->Draw();


  //c1->Print("grvnnLYZ_pt.pdf");
  //c1->Print("grvnnLYZ_pt.png");
  

}

