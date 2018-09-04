import FWCore.ParameterSet.Config as cms

dileptonPreseletor = cms.EDFilter("DileptonPreselector",
	muons = cms.InputTag("slimmedMuons"),
	nMuons = cms.double(2),
	ptCut = cms.double(40),	

)

dileptonPreseletor_AOD = cms.EDFilter("DileptonPreselector_AOD",
  muons = cms.InputTag("muons"),
  nMuons = cms.double(2),
  ptCut = cms.double(40),
)
