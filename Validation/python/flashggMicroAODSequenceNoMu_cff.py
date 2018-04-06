import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggTkVtxMapNoMu_cfi import flashggVertexMapUnique,flashggVertexMapNonUnique,flashggVertexMapForCHS,flashggVertexMapForPUPPI,flashggVertexMapUniqueNoMu,flashggVertexMapNonUniqueNoMu,flashggVertexMapForCHSNoMu,flashggVertexMapForPUPPINoMu
from flashgg.MicroAOD.flashggPhotons_cfi import flashggPhotons
from flashgg.MicroAOD.flashggMuons_cfi import flashggMuons

from flashgg.MicroAOD.flashggLeptonSelectors_cff import flashggSelectedMuons,flashggSelectedElectrons
from flashgg.MicroAOD.flashggMicroAODGenSequence_cff import *

from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import *


eventCount = cms.EDProducer("EventCountProducer")
weightsCount = cms.EDProducer("WeightsCountProducer",
                              generator=cms.InputTag("generator"),
                              pileupInfo=cms.InputTag("slimmedAddPileupInfo"),
                              doObsPileup=cms.untracked.bool(True),
                              minObsPileup=cms.double(-0.5),
                              maxObsPileup=cms.double(100.5),
                              nbinsObsPileup=cms.int32(101),
                              )

flashggMicroAODSequence = cms.Sequence(eventCount
                                       +flashggVertexMapUnique+flashggVertexMapNonUnique
                                       +flashggVertexMapUniqueNoMu+flashggVertexMapNonUniqueNoMu
                                       +flashggMuons*flashggSelectedMuons
                                       +flashggMicroAODGenSequence
                                       )
