import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("ConvertedPhotonResoTreeMaker")

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
                                   fileName = cms.string("ConvertedPhotonResoTree.root")
)

#Sequence builder
#**************************************************************
process.load("flashggDiphotonVtxPlugins/MVATraining/flashggDiphotonVertexSequence_cff")
process.flashggDiPhotons.nVtxSaveInfo = cms.untracked.uint32(999)

process.commissioning = cms.EDAnalyzer('ConvertedPhotonResoTreeMaker',
                                       PhotonTag              = cms.InputTag('flashggPhotons'),
                                       VertexTag              = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       ConversionTag          = cms.InputTag('reducedEgamma', 'reducedConversions'),
                                       SingleLegConversionTag = cms.InputTag('reducedEgamma', 'reducedSingleLegConversions'),
                                       BeamSpotTag            = cms.InputTag('offlineBeamSpot'),
                                       rhoTag                 = cms.InputTag('fixedGridRhoFastjetAll'),
                                       GenParticleTag         = cms.InputTag('prunedGenParticles')
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
