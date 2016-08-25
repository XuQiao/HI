//#include <libgen.h>

void anaphhijing(const int nevents = 0,
                 const int runnumber = 0,
                 const char *pyfname = "phhijing.root")
{
  gSystem->Load("libfun4allfuncs.so");	// framework + reco modules
  gSystem->Load("/direct/phenix+plhf/xuq/phenix/MCgen/event_gen/lib/lib/libAnaPHPythia.so");

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",runnumber);

  /////////////////////////////////////////////////////////////////
  //  Server...
  Fun4AllServer *se = Fun4AllServer::instance();

  /////////////////////////////////////////////////////////////////
  //  Reconstruction Modules...
 
  SubsysReco *anaphhijing = new AnaPHHijing("anaphhijing.root");
  se->registerSubsystem(anaphhijing);

  Fun4AllDstInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  in1->AddFile( pyfname );
  se->registerInputManager(in1);

  se->run(nevents);  // run over all events
  se->End();
}

