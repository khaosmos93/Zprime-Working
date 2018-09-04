import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction_MiniAOD

hardInteraction_MiniAOD.shutUp = cms.bool(True)

cut_trk_quality   = 'isGlobalMuon && ' \
                    'isTrackerMuon && ' \
                    'pt > 53 && ' \
                    'abs(eta) < 2.4 && ' \
                    'abs(dB) < 0.2 && ' \
                    'isolationR03.sumPt / innerTrack.pt < 0.10 && ' \
                    'globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && ' \
                    'globalTrack.hitPattern.numberOfValidPixelHits >= 1'

cut_for_valMuHits = 'isGlobalMuon && ' \
                    'isTrackerMuon && ' \
                    'pt > 53 && ' \
                    'abs(eta) < 2.4 && ' \
                    'abs(dB) < 0.2 && ' \
                    'isolationR03.sumPt / innerTrack.pt < 0.10 && ' \
                    'globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && ' \
                    'globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
                    '( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))'

cut_for_matchedSt = 'isGlobalMuon && ' \
                    'isTrackerMuon && ' \
                    'pt > 53 && ' \
                    'abs(eta) < 2.4 && ' \
                    'abs(dB) < 0.2 && ' \
                    'isolationR03.sumPt / innerTrack.pt < 0.10 && ' \
                    'globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && ' \
                    'globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ' \
                    '( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) )'


#             '( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && ' \
#             '( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))'

HistosGenMatching = cms.EDAnalyzer('HistosGenMatching',
                               muon_src = cms.InputTag('slimmedMuons'),
                               beamspot_src = cms.InputTag('offlineBeamSpot'),
                               vertex_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                               use_bs_and_pv = cms.bool(True),

                               hardInteraction = hardInteraction_MiniAOD,

                               muon_cut_DEN = cms.string(cut_trk_quality),
                               muon_dpt_over_pt_max_DEN = cms.double(0.3),
                               muon_dz_max_DEN = cms.double(0.5),

                               muon_cut_NUM = cms.string(cut_trk_quality),
                               muon_dpt_over_pt_max_NUM = cms.double(0.3),
                               muon_dz_max_NUM = cms.double(0.5),

                               useGenWeight = cms.bool(True),
                               usePileupWeight = cms.bool(False),
                               pileup_src = cms.InputTag('slimmedAddPileupInfo'),
                               vec_PileUpWeight = cms.vdouble(-1),

                               ShutUp = cms.bool(True)
)
