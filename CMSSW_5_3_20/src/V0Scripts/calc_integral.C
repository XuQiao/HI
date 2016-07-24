#include <TFile.h>
#include <TH1F.h>
double calc_integral( const TH1F* hist, float defBinWidth ) 
{
  double histoNbins = hist->GetXaxis()->GetNbins();

  double hist_integral = 0.;
  for( int ndx = 1; ndx <= histoNbins; ndx++ ) {
    //cout << hist->GetBinContent(ndx) << endl;
    hist_integral += 
      ( hist->GetBinContent(ndx)*hist->GetBinWidth(ndx) );
  }
  return hist_integral/defBinWidth;

}

