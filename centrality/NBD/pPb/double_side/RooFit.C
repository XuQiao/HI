{
using namespace RooFit;
TH1F h("h", "h", 3, 0.5, 3.5);
for (int i=0;i<3;i++) {
h.Fill(1.0);
}
for (int i=0;i<5;i++) {
h.Fill(2.0);
}
for (int i=0;i<4;i++) {
h.Fill(3.0);
}
TCanvas c1("c1", "c1", 400, 400);
h.Fit("poly");
h.Draw();
c1.Draw();
//
RooRealVar x ("x", "x", 0.0, 4.0);
RooRealVar a ("a", "a", 0.0, 3.0);
RooPolynomial poly ("poly", "poly", x, RooArgList(a));
double integral_lower = 0.0;
double integral_upper = 10000000.0;
RooRealVar poly_integral ("poly_integral", "poly_integral", integral_lower, integral_upper);
RooAddPdf pdf("pdf", "pdf", RooArgList(poly), RooArgList(poly_integral));
RooDataHist data ("data", "data", x, &h);
//
pdf.fitTo(data, Extended(kTRUE));
//
RooChi2Var chi2 ("chi2", "chi2", pdf, data, Extended(kTRUE));
Double_t chi2_val = chi2.getVal();
//
TCanvas c2("c2", "c2", 400, 400);
RooPlot* plot = x.frame();
data.plotOn(plot, Name("data"));
pdf.plotOn(plot, Name("pdf"));
pdf.plotOn(plot, RooFit::Components(poly), RooFit::LineColor(2));
plot->Draw();
//
cout << "reduced chi2 from RooChi2Var = " << chi2_val/(h.GetNbinsX()-2.0) << endl;
cout << "reduced chi2 from plot = " << plot->chiSquare("pdf", "data", (h.GetNbinsX()-2.0)) << endl;
}
