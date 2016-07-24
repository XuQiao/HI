//#include <libgen.h>

void anaphpythia(const int nevents = 0,
                 const int runnumber = 0,
                 const char *pyfname = "../../macros/phpythia.root")
{
  gSystem->Load("libfun4allfuncs.so");	// framework + reco modules
  gSystem->Load("libAnaPHPythia.so");

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",runnumber);

  /////////////////////////////////////////////////////////////////
  //  Server...
  Fun4AllServer *se = Fun4AllServer::instance();

  /////////////////////////////////////////////////////////////////
  //  Reconstruction Modules...
 
  SubsysReco *anaphpythia = new AnaPHPythia("anaphpythia.root");
  se->registerSubsystem(anaphpythia);

  Fun4AllDstInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  in1->AddFile( pyfname );
  se->registerInputManager(in1);

  se->run(nevents);  // run over all events
  se->End();
}

