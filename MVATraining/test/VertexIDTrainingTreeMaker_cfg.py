import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("VertexIDTrainingTreeMaker")

# geometry and global tag:
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v12')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring(
'/store/mc/RunIIFall17MiniAOD/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/A87F7D05-0403-E811-AE12-C0BFC0E567FE.root'
        )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("VertexIDTrainingTree.root")
)

#Sequence builder
#**************************************************************
process.load("flashggDiphotonVtxPlugins/MVATraining/flashggDiphotonVertexSequence_cff")
process.flashggDiPhotons.nVtxSaveInfo = cms.untracked.uint32(999)
process.flashggDiPhotons.useSingleLeg = cms.bool(True)
process.flashggDiPhotons.sigma1Pix               = cms.double( 0.00800379 )
process.flashggDiPhotons.sigma1Tib               = cms.double( 0.502127   )
process.flashggDiPhotons.sigma1Tob               = cms.double( 5.5138     )
process.flashggDiPhotons.sigma1PixFwd            = cms.double( 0.0318172  )
process.flashggDiPhotons.sigma1Tid               = cms.double( 0.325117   )
process.flashggDiPhotons.sigma1Tec               = cms.double( 1.19907    )
process.flashggDiPhotons.sigma2Pix               = cms.double( 0.0171381  )
process.flashggDiPhotons.sigma2Tib               = cms.double( 0.282616   )
process.flashggDiPhotons.sigma2Tob               = cms.double( 3.5737     )
process.flashggDiPhotons.sigma2PixFwd            = cms.double( 0.0923745  )
process.flashggDiPhotons.sigma2Tid               = cms.double( 0.355705   )
process.flashggDiPhotons.sigma2Tec               = cms.double( 0.863342   )
process.flashggDiPhotons.singlelegsigma1Pix      = cms.double( 0.00879849 )
process.flashggDiPhotons.singlelegsigma1Tib      = cms.double( 1.37155    )
process.flashggDiPhotons.singlelegsigma1Tob      = cms.double( 2.7242     )
process.flashggDiPhotons.singlelegsigma1PixFwd   = cms.double( 0.0596455  )
process.flashggDiPhotons.singlelegsigma1Tid      = cms.double( 0.479279   )
process.flashggDiPhotons.singlelegsigma1Tec      = cms.double( 2.02211    )
process.flashggDiPhotons.singlelegsigma2Pix      = cms.double( 0.0224474  )
process.flashggDiPhotons.singlelegsigma2Tib      = cms.double( 0.594662   )
process.flashggDiPhotons.singlelegsigma2Tob      = cms.double( 0.433137   )
process.flashggDiPhotons.singlelegsigma2PixFwd   = cms.double( 0.137922   )
process.flashggDiPhotons.singlelegsigma2Tid      = cms.double( 0.421378   )
process.flashggDiPhotons.singlelegsigma2Tec      = cms.double( 0.977421   )

process.commissioning = cms.EDAnalyzer('VertexIDTrainingTreeMaker',
                                       DiPhotonTag             = cms.InputTag('flashggDiPhotons'),
                                       VertexTag               = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       VertexCandidateMapTagDz = cms.InputTag('flashggVertexMapUnique'),
                                       BeamSpotTag             = cms.InputTag('offlineBeamSpot'),
                                       rhoTag                  = cms.InputTag('fixedGridRhoFastjetAll'),
                                       GenParticleTag          = cms.InputTag('prunedGenParticles'),
                                       GenEventInfo            = cms.InputTag('generator')
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,setupVIDPhotonSelection
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
#my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
process.flashggPhotons.effAreasConfigFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt")
process.flashggPhotons.egmMvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1Values")

process.p = cms.Path(process.flashggDiphotonVertexSequence
                    *process.egmPhotonIDSequence
                    *process.commissioning
                    )
