  // Macro myfunc.C
   Double_t tsallis(Double_t *x, Double_t *par)
   {
     const Double_t mass = 0.497614;
     const Double_t pi = 3.1415926536;
      Double_t pt =x[0];
      Double_t Et = sqrt(mass*mass+pt*pt)-mass;
      Double_t C = par[0];
      Double_t T = par[1];
      Double_t n = par[2];
      Double_t func = C*pt/pow((1+Et/(n*T)),n);
      return func;
   }
