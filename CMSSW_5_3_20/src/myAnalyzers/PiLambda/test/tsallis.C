//Macro myfunc.C
Double_t tsallis(Double_t *x, Double_t *par)
{
  const Double_t mass = 1.3217;
  const Double_t pi = 3.1415926536;
   Double_t pt =x[0];
   Double_t Et = sqrt(mass*mass+pt*pt)-mass;
   Double_t C = par[0];
   Double_t T = par[1];
   Double_t n = par[2];
   Double_t func = C*pt/pow((1+Et/(n*T)),n);
   return func;
}

Double_t tsallismean(Double_t *x, Double_t *par)
{
  const Double_t mass = 1.3217;
  const Double_t pi = 3.1415926536;
   Double_t pt =x[0];
   Double_t Et = sqrt(mass*mass+pt*pt)-mass;
   Double_t C = par[0];
   Double_t A = par[1];
   Double_t n = par[2];
   Double_t func = C*pt/pow((1+2.0*Et/(A*(n-3.0))),n);
   return func;
}

//2nT(m(n-2)+2nT) / (n^2-3n+2)
float tsallisInt(Double_t m, Double_t n, Double_t T, Double_t c=1) {
  return c*(n*T*((m*(n-2))+n*T) / ((n*n)-(3*n)+2));
}

