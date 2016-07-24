static  const char* clnames[] = {
  "LCGReflex/pat::HeavyIon",
  "LCGReflex/edm::Wrapper<pat::HeavyIon>",
  "LCGReflex/reco::Centrality",
  "LCGReflex/edm::Wrapper<reco::Centrality>",
  "LCGReflex/reco::EvtPlane",
  "LCGReflex/edm::Wrapper<reco::EvtPlane>",
  "LCGReflex/std::vector<reco::Centrality>",
  "LCGReflex/std::vector<reco::EvtPlane>",
  "LCGReflex/edm::Wrapper<std::vector<reco::EvtPlane> >",
  "LCGReflex/edm::Wrapper<std::vector<reco::Centrality> >",
};
extern "C" void SEAL_CAPABILITIES (const char**& names, int& n )
{names = clnames;n = sizeof(clnames)/sizeof(char*);}
