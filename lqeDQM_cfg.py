import FWCore.ParameterSet.Config as cms

process = cms.Process("LQeDQM")
process.load("Leptoquarks.PhysicsDQM.lqeDQM_cfi")

process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.DQM.collectorHost = ''

process.dqmSaver.workflow = cms.untracked.string('/My/Test/DataSet')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#    '/store/mc/Fall08/QCD100to250-madgraph/GEN-SIM-RECO/IDEAL_V9_v4/0017/D6268CB8-F8E0-DD11-9D0F-00188B7AC621.root'
#    '/store/user/santanas/LQ_ue_400_10TeV_eejj/LQ_ue_400_10TeV_eejj/3ccb14828a2cff1a3f1e1cd9ebf9b6b4/PYTHIA6_Exotica_LQ_ue_400_10TeV_eejj_cff_py_GEN_FASTSIM_9.root'
    '/store/user/santanas/LQ_ue_400_10TeV_enuejj/LQ_ue_400_10TeV_enuejj/375c51e7a4b97bca5e29dd325fdfa9c1/PYTHIA6_Exotica_LQ_ue_400_10TeV_enuejjFilter_cff_py_GEN_FASTSIM_9.root'
    )
)
#process.MessageLogger = cms.Service("MessageLogger",
#    destinations = cms.untracked.vstring('detailedInfo',
#        'cout')
#)
#process.ana = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(process.lqeDQM+process.dqmSaver)

