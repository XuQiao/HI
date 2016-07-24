// http://www.hepforge.org/downloads/tglaubermc

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TNamed.h>
#include <TF1.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TString.h>
#include <TEllipse.h>
#endif

#ifndef _runglauber_
#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_
#endif

//---------------------------------------------------------------------------------
void runAndSaveNtuple(Int_t n,
                      Text_t *sysA="Au",
                      Text_t *sysB="Au",
                      Double_t signn=42,
                      Double_t mind=0.4,
                      const char *fname="glau_auau_ntuple.root");

//---------------------------------------------------------------------------------
void runAndSaveNucleons(Int_t n,                    
                        Text_t *sysA="Au",           
                        Text_t *sysB="Au",           
                        Double_t signn=42,           
                        Double_t mind=0.4,
                        Bool_t verbose=0,
                        const char *fname="glau_auau_nucleons.root");

//---------------------------------------------------------------------------------
class TGlauNucleon : public TNamed
{
   private:
      Double32_t fX;            //Position of nucleon
      Double32_t fY;            //Position of nucleon
      Double32_t fZ;            //Position of nucleon
      Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
      Int_t      fNColl;        //Number of binary collisions

   public:
      TGlauNucleon() : fX(0), fY(0), fZ(0), fInNucleusA(0), fNColl(0) {}
      virtual   ~TGlauNucleon() {}

      void       Collide()            {fNColl++;}
      Int_t      GetNColl()     const {return fNColl;}
      Double_t   GetX()         const {return fX;}
      Double_t   GetY()         const {return fY;}
      Double_t   GetZ()         const {return fZ;}
      Bool_t     IsInNucleusA() const {return fInNucleusA;}
      Bool_t     IsInNucleusB() const {return !fInNucleusA;}
      Bool_t     IsSpectator()  const {return !fNColl;}
      Bool_t     IsWounded()    const {return fNColl;}
      void       Reset()              {fNColl=0;}
      void       SetInNucleusA()      {fInNucleusA=1;}
      void       SetInNucleusB()      {fInNucleusA=0;}
      void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}

      ClassDef(TGlauNucleon,1)
};

//---------------------------------------------------------------------------------
class TGlauNucleus : public TNamed
{
   private:
      Int_t      fN;          //Number of nucleons
      Double_t   fR;          //Parameters of function
      Double_t   fA;          //Parameters of function
      Double_t   fW;          //Parameters of function
      Double_t   fMinDist;    //Minimum separation distance
      Int_t      fF;          //Type of radial distribution
      Int_t      fTrials;     //Store trials needed to complete nucleus
      TF1*       fFunction;   //Probability density function rho(r)
      TObjArray* fNucleons;   //Array of nucleons

      void       Lookup(Text_t* name);

   public:
      TGlauNucleus(Text_t* iname="Au", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
      virtual ~TGlauNucleus();

      using      TObject::Draw;
      void       Draw(Double_t xs, Int_t col);
      Int_t      GetN()             const {return fN;}
      Double_t   GetR()             const {return fR;}
      Double_t   GetA()             const {return fA;}
      Double_t   GetW()             const {return fW;}
      TObjArray *GetNucleons()      const {return fNucleons;}
      Int_t      GetTrials()        const {return fTrials;}
      void       SetN(Int_t in)           {fN=in;}
      void       SetR(Double_t ir);
      void       SetA(Double_t ia);
      void       SetW(Double_t iw);
      void       SetMinDist(Double_t min) {fMinDist=min;}
      void       ThrowNucleons(Double_t xshift=0.);

      ClassDef(TGlauNucleus,1)
};

//---------------------------------------------------------------------------------
class TGlauberMC : public TNamed
{
   private:
      TGlauNucleus fANucleus;       //Nucleus A
      TGlauNucleus fBNucleus;       //Nucleus B
      Double_t     fXSect;          //Nucleon-nucleon cross section
      TObjArray*   fNucleonsA;      //Array of nucleons in nucleus A
      TObjArray*   fNucleonsB;      //Array of nucleons in nucleus B
      Int_t        fAN;             //Number of nucleons in nucleus A
      Int_t        fBN;             //Number of nucleons in nucleus B
      TNtuple*     fnt;             //Ntuple for results (created, but not deleted)
      Double_t     fMeanX2;         //<x^2> of wounded nucleons
      Double_t     fMeanY2;         //<y^2> of wounded nucleons
      Double_t     fMeanXY;         //<xy> of wounded nucleons
      Double_t     fMeanXParts;     //<x> of wounded nucleons
      Double_t     fMeanYParts;     //<x> of wounded nucleons
      Double_t     fMeanXSystem;    //<x> of all nucleons
      Double_t     fMeanYSystem;    //<x> of all nucleons  
      Double_t     fMeanX_A;        //<x> of nucleons in nucleus A
      Double_t     fMeanY_A;        //<x> of nucleons in nucleus A
      Double_t     fMeanX_B;        //<x> of nucleons in nucleus B
      Double_t     fMeanY_B;        //<x> of nucleons in nucleus B
      Double_t     fB_MC;           //Impact parameter (b)
      Int_t        fEvents;         //Number of events with at least one collision
      Int_t        fTotalEvents;    //All events within selected impact parameter range
      Double_t     fBMin;           //Minimum impact parameter to be generated
      Double_t     fBMax;           //Maximum impact parameter to be generated
      Int_t        fMaxNpartFound;  //Largest value of Npart obtained
      Int_t        fNpart;          //Number of wounded (participating) nucleons in current event
      Int_t        fNcoll;          //Number of binary collisions in current event
      Double_t     fSx2;            //Variance of x of wounded nucleons
      Double_t     fSy2;            //Variance of y of wounded nucleons
      Double_t     fSxy;            //Covariance of x and y of wounded nucleons

      Bool_t       CalcResults(Double_t bgen);
      Bool_t       CalcEvent(Double_t bgen);

   public:
      TGlauberMC(Text_t* NA = "Au", Text_t* NB = "Au", Double_t xsect = 42);
      virtual     ~TGlauberMC() {Reset();}

      void         Draw(Option_t* option);
      Double_t     GetB()               const {return fB_MC;}
      Double_t     GetBMin()            const {return fBMin;}
      Double_t     GetBMax()            const {return fBMax;}
      Int_t        GetNcoll()           const {return fNcoll;}
      Int_t        GetNpart()           const {return fNpart;}
      Int_t        GetNpartFound()      const {return fMaxNpartFound;}
      TNtuple*     GetNtuple()          const {return fnt;}
      TObjArray   *GetNucleons();
      Double_t     GetTotXSect()        const;
      Double_t     GetTotXSectErr()     const;
      Bool_t       NextEvent(Double_t bgen=-1);
      void         Reset()                    {delete fnt; fnt=0; }
      void         Run(Int_t nevents);
      void         SetBmin(Double_t bmin)      {fBMin = bmin;}
      void         SetBmax(Double_t bmax)      {fBMax = bmax;}
      void         SetMinDistance(Double_t d)  {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}

      static 
        void       PrintVersion()             {cout << "TGlauberMC " << Version() << endl;}
      static 
        const char *Version()                 {return "v1.1";}

      ClassDef(TGlauberMC,1)
};

//---------------------------------------------------------------------------------
void runAndSaveNtuple(Int_t n,
                      Text_t *sysA,
                      Text_t *sysB,
                      Double_t signn,
                      Double_t mind,
                      const char *fname)
{
   TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn);
   mcg->SetMinDistance(mind);
   mcg->Run(n);
   TNtuple  *nt=mcg->GetNtuple();
   TFile out(fname,"recreate",fname,9);
   if(nt) nt->Write();
   out.Close();
}

//---------------------------------------------------------------------------------
void runAndSaveNucleons(Int_t n,                    
                        Text_t *sysA,           
                        Text_t *sysB,           
                        Double_t signn,
                        Double_t mind,
                        Bool_t verbose,
                        const char *fname)
{
   TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn);
   mcg->SetMinDistance(mind);
   TFile *out=0;
   if(fname) 
      out=new TFile(fname,"recreate",fname,9);

   for(Int_t ievent=0;ievent<n;ievent++){

      //get an event with at least one collision
      while(!mcg->NextEvent()) {}

      //access, save and (if wanted) print out nucleons
      TObjArray* nucleons=mcg->GetNucleons();
      if(!nucleons) continue;
      if(out)
         nucleons->Write(Form("nucleonarray%d",ievent),TObject::kSingleKey);

      if(verbose) {
         cout<<endl<<endl<<"EVENT NO: "<<ievent<<endl;
         cout<<"B = "<<mcg->GetB()<<"  Npart = "<<mcg->GetNpart()<<endl<<endl;
         printf("Nucleus\t X\t Y\t Z\tNcoll\n");
         Int_t nNucls=nucleons->GetEntries();
         for(Int_t iNucl=0;iNucl<nNucls;iNucl++) {
            TGlauNucleon *nucl=(TGlauNucleon *)nucleons->At(iNucl);
            Char_t nucleus='A';
            if(nucl->IsInNucleusB()) nucleus='B';
            Double_t x=nucl->GetX();
            Double_t y=nucl->GetY();
            Double_t z=nucl->GetZ();
            Int_t ncoll=nucl->GetNColl();
            printf("   %c\t%2.2f\t%2.2f\t%2.2f\t%3d\n",nucleus,x,y,z,ncoll);
         }
      }
   }
   if(out) delete out;
}

//---------------------------------------------------------------------------------
ClassImp(TGlauNucleus)
//---------------------------------------------------------------------------------

TGlauNucleus::TGlauNucleus(Text_t* iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
   fN(iN),fR(iR),fA(ia),fW(iw),fMinDist(-1),
   fF(0),fTrials(0),fFunction(ifunc),
   fNucleons(0)
{
   if (fN==0) {
      cout << "Setting up nucleus " << iname << endl;
      Lookup(iname);
   }
}

TGlauNucleus::~TGlauNucleus()
{
   if (fNucleons) {
      delete fNucleons;
   }
   delete fFunction;
}

void TGlauNucleus::Draw(Double_t xs, Int_t col)
{
   Double_t r = 0.5*sqrt(xs/TMath::Pi()/10.);
   TEllipse e;
   e.SetLineColor(col);
   e.SetFillColor(0);
   e.SetLineWidth(1);

   for (Int_t i = 0;i<fNucleons->GetEntries();++i) {
      TGlauNucleon* gn = (TGlauNucleon*) fNucleons->At(i);
      e.SetLineStyle(1);
      if (gn->IsSpectator()) e.SetLineStyle(3);
      e.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
   }
}

void TGlauNucleus::Lookup(Text_t* name)
{
   SetName(name);

   if      (TString(name) == "p")    {fN = 1;   fR = 0.6;   fA = 0;      fW =  0;      fF = 0;}
   else if (TString(name) == "d")    {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 1;}
   else if (TString(name) == "dh")   {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 3;}
   else if (TString(name) == "dhh")  {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 4;}
   else if (TString(name) == "O")    {fN = 16;  fR = 2.608; fA = 0.513;  fW = -0.051;  fF = 1;}
   else if (TString(name) == "Si")   {fN = 28;  fR = 3.34;  fA = 0.580;  fW = -0.233;  fF = 1;}
   else if (TString(name) == "S")    {fN = 32;  fR = 2.54;  fA = 2.191;  fW =  0.16;   fF = 2;}
   else if (TString(name) == "Ca")   {fN = 40;  fR = 3.766; fA = 0.586;  fW = -0.161;  fF = 1;}
   else if (TString(name) == "Ni")   {fN = 58;  fR = 4.309; fA = 0.517;  fW = -0.1308; fF = 1;}
   else if (TString(name) == "Cu")   {fN = 63;  fR = 4.2;   fA = 0.596;  fW =  0;      fF = 1;}
   else if (TString(name) == "W")    {fN = 186; fR = 6.58;  fA = 0.480;  fW =  0;      fF = 1;}
   else if (TString(name) == "Au")   {fN = 197; fR = 6.38;  fA = 0.535;  fW =  0;      fF = 1;}
   else if (TString(name) == "Pb")   {fN = 208; fR = 6.62;  fA = 0.546;  fW =  0;      fF = 1;}
   else if (TString(name) == "U")    {fN = 238; fR = 6.81;  fA = 0.6;    fW =  0;      fF = 1;}
   else {
      cout << "Could not find nucleus " << name << endl;
      return;
   }

   switch (fF)
   {
      case 0: // Proton
         fFunction = new TF1("prot","x*x*exp(-x/[0])",0,10);
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
         fFunction = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,15);
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 2: // 3pG
         fFunction = new TF1("3pg","x*x*(1+[2]*(x/[0])**2)/(1+exp((x**2-[0]**2)/[1]**2))",0,15);
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 3: // Hulthen
         fFunction = new TF1("f3","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,10);
         fFunction->SetParameters(1/4.38,1/.85);
         break;
      case 4: // Hulthen HIJING
         fFunction = new TF1("f4","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,20);
         fFunction->SetParameters(2/4.38,2/.85);
         break;
      default:
         cerr << "Could not find function type " << fF << endl;
         return;
   }
}

void TGlauNucleus::SetR(Double_t ir)
{
   fR = ir;
   switch (fF)
   {
      case 0: // Proton
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
         fFunction->SetParameter(0,fR);
         break;
      case 2: // 3pG
         fFunction->SetParameter(0,fR);
         break;
   }
}

void TGlauNucleus::SetA(Double_t ia)
{
   fA = ia;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(1,fA);
         break;
      case 2: // 3pG
         fFunction->SetParameter(1,fA);
         break;
   }
}

void TGlauNucleus::SetW(Double_t iw)
{
   fW = iw;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(2,fW);
         break;
      case 2: // 3pG
         fFunction->SetParameter(2,fW);
         break;
   }
}

void TGlauNucleus::ThrowNucleons(Double_t xshift)
{
   if (fNucleons==0) {
      fNucleons=new TObjArray(fN);
      fNucleons->SetOwner();
      for(Int_t i=0;i<fN;i++) {
	 TGlauNucleon *nucleon=new TGlauNucleon(); 
	 fNucleons->Add(nucleon); 
      }
   } 
   
   fTrials = 0;

   Double_t sumx=0;       
   Double_t sumy=0;       
   Double_t sumz=0;       

   Bool_t hulthen = (TString(GetName())=="dh");
   if (fN==2 && hulthen) { //special treatmeant for Hulten

      Double_t r = fFunction->GetRandom()/2;
      Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
      Double_t ctheta = 2*gRandom->Rndm() - 1 ;
      Double_t stheta = sqrt(1-ctheta*ctheta);
     
      TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
      TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
      nucleon1->Reset();
      nucleon1->SetXYZ(r * stheta * cos(phi) + xshift,
		       r * stheta * sin(phi),
		       r * ctheta);
      nucleon2->Reset();
      nucleon2->SetXYZ(-nucleon1->GetX() + 2*xshift,
		       -nucleon1->GetY(),
		       -nucleon1->GetZ());
      fTrials = 1;
      return;
   }

   for (Int_t i = 0; i<fN; i++) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while(1) {
         fTrials++;
         Double_t r = fFunction->GetRandom();
         Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
         Double_t ctheta = 2*gRandom->Rndm() - 1 ;
         Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
         Double_t x = r * stheta * cos(phi) + xshift;
         Double_t y = r * stheta * sin(phi);      
         Double_t z = r * ctheta;      
         nucleon->SetXYZ(x,y,z);
         if(fMinDist<0) break;
         Bool_t test=1;
         for (Int_t j = 0; j<i; j++) {
            TGlauNucleon *other=(TGlauNucleon*)fNucleons->At(j);
            Double_t xo=other->GetX();
            Double_t yo=other->GetY();
            Double_t zo=other->GetZ();
            Double_t dist = TMath::Sqrt((x-xo)*(x-xo)+
                                       (y-yo)*(y-yo)+
                                       (z-zo)*(z-zo));
	       
            if(dist<fMinDist) {
               test=0;
               break;
            }
         }
         if (test) break; //found nucleuon outside of mindist
      }
           
      sumx += nucleon->GetX();
      sumy += nucleon->GetY();
      sumz += nucleon->GetZ();
   }
      
   if(1) { // set the centre-of-mass to be at zero (+xshift)
      sumx = sumx/fN;  
      sumy = sumy/fN;  
      sumz = sumz/fN;  
      for (Int_t i = 0; i<fN; i++) {
         TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
         nucleon->SetXYZ(nucleon->GetX()-sumx-xshift,
                         nucleon->GetY()-sumy,
                         nucleon->GetZ()-sumz);
      }
   }
}

//---------------------------------------------------------------------------------
ClassImp(TGlauberMC)
//---------------------------------------------------------------------------------

TGlauberMC::TGlauberMC(Text_t* NA, Text_t* NB, Double_t xsect) :
   fANucleus(NA),fBNucleus(NB),
   fXSect(0),fNucleonsA(0),fNucleonsB(0),fnt(0),
   fMeanX2(0),fMeanY2(0),fMeanXY(0),fMeanXParts(0),
   fMeanYParts(0),fMeanXSystem(0),fMeanYSystem(0),  
   fMeanX_A(0),fMeanY_A(0),fMeanX_B(0),fMeanY_B(0),fB_MC(0),
   fEvents(0),fTotalEvents(0),fBMin(0),fBMax(0),fMaxNpartFound(0),
   fNpart(0),fNcoll(0),fSx2(0),fSy2(0),fSxy(0)
{
   fBMin = 0;
   fBMax = 20;
   fXSect = xsect;
   
   TString name(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
   TString title(Form("Glauber %s+%s Version",fANucleus.GetName(),fBNucleus.GetName()));
   SetName(name);
   SetTitle(title);
}

Bool_t TGlauberMC::CalcEvent(Double_t bgen)
{
   // prepare event
   fANucleus.ThrowNucleons(-bgen/2.);
   fNucleonsA = fANucleus.GetNucleons();
   fAN = fANucleus.GetN();
   for (Int_t i = 0; i<fAN; i++) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      nucleonA->SetInNucleusA();
   }
   fBNucleus.ThrowNucleons(bgen/2.);
   fNucleonsB = fBNucleus.GetNucleons();
   fBN = fBNucleus.GetN();
   for (Int_t i = 0; i<fBN; i++) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      nucleonB->SetInNucleusB();
   }

   // "ball" diameter = distance at which two balls interact
   Double_t d2 = (Double_t)fXSect/(TMath::Pi()*10); // in fm^2

   // for each of the A nucleons in nucleus B
   for (Int_t i = 0; i<fBN; i++) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      for (Int_t j = 0 ; j < fAN ;j++) {
	 TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(j));
         Double_t dx = nucleonB->GetX()-nucleonA->GetX();
         Double_t dy = nucleonB->GetY()-nucleonA->GetY();
         Double_t dij = dx*dx+dy*dy;
         if (dij < d2) {
            nucleonB->Collide();
            nucleonA->Collide();
         }
      }
  }
   
  return CalcResults(bgen);
}

Bool_t TGlauberMC::CalcResults(Double_t bgen)
{
   // calc results for the given event
   fNpart=0;
   fNcoll=0;
   fMeanX2=0;
   fMeanY2=0;
   fMeanXY=0;
   fMeanXParts=0;
   fMeanYParts=0;
   fMeanXSystem=0;
   fMeanYSystem=0;
   fMeanX_A=0;
   fMeanY_A=0;
   fMeanX_B=0;
   fMeanY_B=0;
  
   for (Int_t i = 0; i<fAN; i++) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      Double_t xA=nucleonA->GetX();
      Double_t yA=nucleonA->GetY();
      fMeanXSystem  += xA;
      fMeanYSystem  += yA;
      fMeanX_A  += xA;
      fMeanY_A  += yA;

      if(nucleonA->IsWounded()) {
         fNpart++;
         fMeanXParts  += xA;
         fMeanYParts  += yA;
         fMeanX2 += xA * xA;
         fMeanY2 += yA * yA;
         fMeanXY += xA * yA;
      }
   }

   for (Int_t i = 0; i<fBN; i++) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      Double_t xB=nucleonB->GetX();
      Double_t yB=nucleonB->GetY();
      fMeanXSystem  += xB;
      fMeanYSystem  += yB;
      fMeanX_B  += xB;
      fMeanY_B  += yB;

      if(nucleonB->IsWounded()) {
         fNpart++;
         fMeanXParts  += xB;
         fMeanYParts  += yB;
         fMeanX2 += xB * xB;
         fMeanY2 += yB * yB;
         fMeanXY += xB * yB;
	 fNcoll += nucleonB->GetNColl();
      }
   }

   if (fNpart>0) {
      fMeanXParts /= fNpart;
      fMeanYParts /= fNpart;
      fMeanX2 /= fNpart;
      fMeanY2 /= fNpart;
      fMeanXY /= fNpart;
   } else {
      fMeanXParts = 0;
      fMeanYParts = 0;
      fMeanX2 = 0;
      fMeanY2 = 0;
      fMeanXY = 0;
   }
   
   if(fAN+fBN>0) {
      fMeanXSystem /= (fAN + fBN);
      fMeanYSystem /= (fAN + fBN);
   } else {
      fMeanXSystem = 0;
      fMeanYSystem = 0;
   }
   if(fAN>0) {
      fMeanX_A /= fAN;
      fMeanY_A /= fAN;
   } else {
      fMeanX_A = 0;
      fMeanY_A = 0;
   }

   if(fBN>0) {
      fMeanX_B /= fBN;
      fMeanY_B /= fBN;
   } else {
      fMeanX_B = 0;
      fMeanY_B = 0;
   }

   fSx2=fMeanX2-(fMeanXParts*fMeanXParts);
   fSy2=fMeanY2-(fMeanYParts*fMeanYParts);
   fSxy=fMeanXY-fMeanXParts*fMeanYParts;

   fB_MC = bgen;

   fTotalEvents++;
   if (fNpart>0) fEvents++;

   if (fNpart==0) return kFALSE;
   if (fNpart > fMaxNpartFound) fMaxNpartFound = fNpart;

   return kTRUE;
}

void TGlauberMC::Draw(Option_t* /*option*/)
{
   fANucleus.Draw(fXSect, 2);
   fBNucleus.Draw(fXSect, 4);

   TEllipse e;
   e.SetFillColor(0);
   e.SetLineColor(1);
   e.SetLineStyle(2);
   e.SetLineWidth(1);
   e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
   e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
}

Double_t TGlauberMC::GetTotXSect() const
{
   return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBMax*fBMax/100;
}

Double_t TGlauberMC::GetTotXSectErr() const
{
   return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * 
      TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

TObjArray *TGlauberMC::GetNucleons() 
{
   if(!fNucleonsA || !fNucleonsB) return 0;
   fNucleonsA->SetOwner(0);
   fNucleonsB->SetOwner(0);
   TObjArray *allnucleons=new TObjArray(fAN+fBN);
   allnucleons->SetOwner();
   for (Int_t i = 0; i<fAN; i++) {
      allnucleons->Add(fNucleonsA->At(i));
   }
   for (Int_t i = 0; i<fBN; i++) {
      allnucleons->Add(fNucleonsB->At(i));
   }
   return allnucleons;
}

Bool_t TGlauberMC::NextEvent(Double_t bgen)
{
   if(bgen<0) 
      bgen = TMath::Sqrt((fBMax*fBMax-fBMin*fBMin)*gRandom->Rndm()+fBMin*fBMin);

   return CalcEvent(bgen);
}

void TGlauberMC::Run(Int_t nevents)
{
   TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
   TString title(Form("%s + %s (x-sect = %d mb)",fANucleus.GetName(),fBNucleus.GetName(),(Int_t) fXSect));
   if (fnt == 0) {
      fnt = new TNtuple(name,title,
                        "Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB");
      fnt->SetDirectory(0);
   }

   for (int i = 0;i<nevents;i++) {

      while(!NextEvent()) {}

      Float_t v[17];
      v[0]  = GetNpart();
      v[1]  = GetNcoll();
      v[2]  = fB_MC;
      v[3]  = fMeanXParts;
      v[4]  = fMeanYParts;
      v[5]  = fMeanX2;
      v[6]  = fMeanY2;
      v[7]  = fMeanXY;
      v[8]  = fSx2;
      v[9]  = fSy2;
      v[10] = fSxy;
      v[11] = fMeanXSystem;
      v[12] = fMeanYSystem;
      v[13] = fMeanX_A;
      v[14] = fMeanY_A;
      v[16] = fMeanX_B;
      v[17] = fMeanY_B;

      fnt->Fill(v);

      if (!(i%50)) cout << "Event # " << i << " x-sect = " << GetTotXSect() << " +- " << GetTotXSectErr() << " b        \r" << flush;  
   }
   cout << endl << "Done!" << endl;
}
#endif
