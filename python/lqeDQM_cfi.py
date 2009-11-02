import FWCore.ParameterSet.Config as cms

# DQM monitor module for first generation leptoquarks
lqeDQM = cms.EDAnalyzer("LQeDQM",
            elecTriggerPathToPass    = cms.string("HLT_Ele10_LW_L1R"),
            muonTriggerPathToPass    = cms.string("HLT_Mu9"),
            jetPtCut                 = cms.double("50.0"),
            elePtCut                 = cms.double("50.0"),
            metCut                 = cms.double("50.0"),
            triggerResultsCollection = cms.InputTag("TriggerResults", "", "HLT"),
            muonCollection           = cms.InputTag("muons"),
            electronCollection       = cms.InputTag("gsfElectrons"),
            caloJetCollection        = cms.InputTag("sisCone5CaloJets"),
            caloMETCollection        = cms.InputTag("corMetGlobalMuons"),
            genParticleCollection    = cms.InputTag("genParticles")
)
