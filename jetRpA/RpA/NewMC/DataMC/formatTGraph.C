void formatTGraph(TGraphErrors& t1, double Ymax, double Ymin, double& col, double& Mstyle, double& fill, double& style, const char* titx, const char* tity ){
  //void format1Dhisto(TH1& h1, string& xTitle, double Ymax, double Ymin){

  //h1.SetTitle(";XXXX;XXXX");
  if(Ymax!=-1 && Ymin!=-1) t1.GetYaxis()->SetRangeUser(Ymax, Ymin);
  //if(Ymax==-1 && Ymin!=-1) t1.GetYaxis()->SetMinimum(Ymin);
  t1.SetMarkerColor(col);
  t1.SetMarkerStyle(Mstyle);
  t1.SetLineColor(col);
  t1.SetFillColor(fill);
  t1.SetFillStyle(style);
  t1.SetMarkerSize(1.1);
  t1.GetXaxis()->SetTitle(titx);
  t1.GetYaxis()->SetTitle(tity);
  t1.GetXaxis()->CenterTitle();
  t1.GetYaxis()->CenterTitle();
  //cout<<"The title is : "<<tit<<endl;

  return;
}

void format1Dhisto(TH2& h2, double Ymax, double Ymin, double& col, double& Mstyle, double& fill, double& style, const char* titx, const char* tity ){
  //void format1Dhisto(TH1& h1, string& xTitle, double Ymax, double Ymin){

  //h2.SetTitle(";XXXX;XXXX");
  if(Ymax!=-1 && Ymin!=-1) h2.GetYaxis()->SetRangeUser(Ymax, Ymin);
  //if(Ymax==-1 && Ymin!=-1) h2.GetYaxis()->SetMinimum(Ymin);
  //h2.SetMarkerColor(col);
  //h2.SetMarkerStyle(Mstyle);
  //h2.SetLineColor(col);
  //h2.SetFillColor(fill);
  //h2.SetFillStyle(style);
  h2.GetXaxis()->SetTitle(titx);
  h2.GetYaxis()->SetTitle(tity);
  h2.GetXaxis()->CenterTitle();
  h2.GetYaxis()->CenterTitle();
  //cout<<"The title is : "<<tit<<endl;

  return;
}







