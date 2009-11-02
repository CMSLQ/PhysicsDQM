#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Leptoquarks/PhysicsDQM/interface/LQeDQM.h"

DEFINE_SEAL_MODULE();
// DEFINE_ANOTHER_FWK_MODULE(BPhysicsOniaDQM);
// DEFINE_ANOTHER_FWK_MODULE(EwkDQM);
// DEFINE_ANOTHER_FWK_MODULE(QcdPhotonsDQM);
// DEFINE_ANOTHER_FWK_MODULE(QcdHighPtDQM);
// DEFINE_ANOTHER_FWK_MODULE(TopDiLeptonDQM);
DEFINE_ANOTHER_FWK_MODULE(LQeDQM);
