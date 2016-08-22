
void
loadLib(const char* lib)
{
  std::cout << "Load lib " << lib << std::endl;
  int stat = gSystem->Load(lib);
  if ( stat < 0 ) 
    {
      std::string newlib = ".libs/";
      newlib += lib;
      gSystem->Load(newlib.c_str());
    }
}

void
phhijing(const int nevents = 10, const int outputnameindex = 0)
{ 
  //TString temp = Form("/phenix/scratch/xuq/Hijing/phhijing_dAu_%d.root",outputnameindex);
  TString temp = Form("/phenix/scratch/xuq/Hijing/phhijing_He3Au_%d.root",outputnameindex);
  //TString temp = Form("./phhijing_He3Au_%d.root",outputnameindex);
  char *outputname = temp.Data();
  gSystem->Load("libfun4allfuncs.so");	// framework + reco modules
  gSystem->Load("libsimreco.so");	// framework + reco modules
  loadLib("/direct/phenix+plhf/xuq/phenix/MCgen/event_gen/lib/lib/libPHHijing");

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",0);

  /////////////////////////////////////////////////////////////////
  //  Server...
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  /////////////////////////////////////////////////////////////////
  //  Reconstruction Modules...
  
  SubsysReco *sync = new SyncSimreco();
  se->registerSubsystem(sync);

  PHHijing* ss = new PHHijing();
  //  ss->SetSeed(99);
  //ss->SetOscarOut("pythia.osc");	// Set this to write to an oscar file
  //ss->SetSeed(1999);			// Set your own seed, otherwise, defaults to seed = time

  se->registerSubsystem(ss);

  // Example: Turn off all B and D decays
/*
  ss->SetDecay(411,0);
  ss->SetDecay(-411,0);
  ss->SetDecay(421,0);
  ss->SetDecay(-421,0);
  ss->SetDecay(431,0);
  ss->SetDecay(-431,0);
  
  ss->SetDecay(511,0);
  ss->SetDecay(-511,0);
  ss->SetDecay(521,0);
  ss->SetDecay(-521,0);
  ss->SetDecay(531,0);
  ss->SetDecay(-531,0);
*/

  //** Here you could put a module to trigger on selected events
  //** Select only one of the trigger types below
  //PHPyTrigger *phpytrigger = new PHPyTrigger();
  //phpytrigger->SetTrigger( PHPyTrigger::CHIC_MUONARM );
  //phpytrigger->SetTrigger( PHPyTrigger::PI0_CENTARM );
  //phpytrigger->SetEThreshold( 2.0 );
  //se->registerSubsystem(phpytrigger);
  //PHPyTrigger* trigger = new PHFvtxTrigger();
  //se->registerSubsystem(trigger);

  Fun4AllDummyInputManager *in1 = new Fun4AllDummyInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  TString OUTPUT = outputname;
  Fun4AllDstOutputManager *dst_output_mgr  = new Fun4AllDstOutputManager("PHHIJING",OUTPUT.Data());
  dst_output_mgr->AddNode("Sync");
  dst_output_mgr->AddNode("PHHijingHeader");
  dst_output_mgr->AddNode("PHHijing");
  se->registerOutputManager(dst_output_mgr);

  se->run(nevents);  // run over all events
  se->End();
}

