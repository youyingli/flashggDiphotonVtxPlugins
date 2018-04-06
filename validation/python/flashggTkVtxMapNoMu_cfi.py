import FWCore.ParameterSet.Config as cms

flashggVertexMapUnique = cms.EDProducer('FlashggDzVertexMapProducer',
                                        PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                        VertexTag=cms.InputTag('offlinePrimaryVerticesWithMu'),
                                        MaxAllowedDz=cms.double(0.2),
                                        UseEachTrackOnce=cms.bool(True)
                                        )

flashggVertexMapNonUnique = cms.EDProducer('FlashggDzVertexMapProducer',
                                           PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                           VertexTag=cms.InputTag('offlinePrimaryVerticesWithMu'),
                                           MaxAllowedDz=cms.double(0.2), 
                                           UseEachTrackOnce=cms.bool(False)
                                           )

flashggVertexMapForCHS = cms.EDProducer('FlashggVertexMapFromCandidateProducer', # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#PV_Assignment
                                        PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                        VertexTag=cms.InputTag('offlinePrimaryVerticesWithMu'),
                                        FromPVCut=cms.uint32(0), # "The definition normally used for isolation calculations is fromPV() > 1; 
                                                                 # the definition used for CHS subtraction in jets is fromPV() > 0."
                                        FromPVCutIfPassDz=cms.uint32(0),
                                        DzCut=cms.double(999999.))


                                                   
flashggVertexMapForPUPPI = cms.EDProducer('FlashggVertexMapFromCandidateProducer',
                                          PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                          VertexTag=cms.InputTag('offlinePrimaryVerticesWithMu'),
                                          FromPVCut=cms.uint32(2), # if used in fit, keep it
                                          FromPVCutIfPassDz=cms.uint32(0), # if loose or tight, check dz cut
                                          DzCut=cms.double(0.3))


flashggVertexMapUniqueNoMu = cms.EDProducer('FlashggDzVertexMapProducer',
                                        PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                        VertexTag=cms.InputTag('offlinePrimaryVerticesNoMu'),
                                        MaxAllowedDz=cms.double(0.2),
                                        UseEachTrackOnce=cms.bool(True)
                                        )

flashggVertexMapNonUniqueNoMu = cms.EDProducer('FlashggDzVertexMapProducer',
                                           PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                           VertexTag=cms.InputTag('offlinePrimaryVerticesNoMu'),
                                           MaxAllowedDz=cms.double(0.2), 
                                           UseEachTrackOnce=cms.bool(False)
                                           )

flashggVertexMapForCHSNoMu = cms.EDProducer('FlashggVertexMapFromCandidateProducer', # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#PV_Assignment
                                        PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                        VertexTag=cms.InputTag('offlinePrimaryVerticesNoMu'),
                                        FromPVCut=cms.uint32(0), # "The definition normally used for isolation calculations is fromPV() > 1; 
                                                                 # the definition used for CHS subtraction in jets is fromPV() > 0."
                                        FromPVCutIfPassDz=cms.uint32(0),
                                        DzCut=cms.double(999999.))


                                                   
flashggVertexMapForPUPPINoMu = cms.EDProducer('FlashggVertexMapFromCandidateProducer',
                                          PFCandidatesTag=cms.InputTag('packedPFCandidates'),
                                          VertexTag=cms.InputTag('offlinePrimaryVerticesNoMu'),
                                          FromPVCut=cms.uint32(2), # if used in fit, keep it
                                          FromPVCutIfPassDz=cms.uint32(0), # if loose or tight, check dz cut
                                          DzCut=cms.double(0.3))

                                          
