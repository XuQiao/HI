#ifndef __CINT__
# include <correlations/recursive/FromQVector.hh>
# include <iostream>
#endif

Bool_t   NextEvent();        // Load next event and return true
Bool_t   NextObservation();  // Load next observation and return true
Double_t GetPhi();           // Get current phi observation
Double_t GetWeight();        // Get current weight on observation

void Test() {
  correlations::Result  r;
  correlations::HarmonicVector h(n);
  correlations::QVector q(h);
  correlations::recursive::FromQVector c(q);

  while (NextEvent()) {
    q.reset();

    while (NextObservation()) {
      correlations::Real phi    = GetPhi();
      correlations::Real Weight = GetWeight();

      q.fill(phi, weight);
    }
    r += c.calculate(h);
  }
  std::cout << r.eval() << std::endl;
}
