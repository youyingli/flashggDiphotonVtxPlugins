import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("vertexProbTrainingTreeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")

# geometry and global tag:

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = '94X_mc2017_realistic_v10'

#**************************************************************

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

#chnage with input produced with CMSSW 7_4 release
readFiles.extend( [
       '/store/mc/RunIIFall17MiniAOD/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/A87F7D05-0403-E811-AE12-C0BFC0E567FE.root' 
       ] );

secFiles.extend( [
               ] )


#**************************************************************


process.load("flashgg/MicroAOD/flashggMicroAODSequence_cff")
process.flashggDiPhotons.nVtxSaveInfo = cms.untracked.uint32(999) 
process.flashggDiPhotons.convProbCut = cms.untracked.double(1.0E-6) 
#process.flashggDiPhotons.vertexIdMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_noSL_2015.xml")
process.flashggDiPhotons.vertexIdMVAweightfile = cms.FileInPath("flashgg/Validation/data/TMVAClassification_BDTVtxId.weights.xml")
process.flashggDiPhotons.vertexProbMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_noSL_2015.xml")
process.flashggDiPhotons.useSingleLeg=cms.bool(True)

process.commissioning = cms.EDAnalyzer('vertexProbTrainingTreeMaker',
                                       PhotonTag=cms.untracked.InputTag('flashggPhotons'),
                                       DiPhotonTag = cms.untracked.InputTag('flashggDiPhotons'),
                                       VertexTag=cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                       VertexCandidateMapTagDz=cms.InputTag('flashggVertexMapUnique'),
                                       ConversionTag=cms.untracked.InputTag("reducedEgamma","reducedConversions"), 
                                       BeamSpotTag=cms.untracked.InputTag('offlineBeamSpot'),
                                       evWeight = cms.untracked.double(1.0)
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,setupVIDPhotonSelection
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService")
process.RandomNumberGeneratorService.flashggRandomizedPhotons = cms.PSet(
                  initialSeed = cms.untracked.uint32(16253245)
                          )

#**************************************************************

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("vtxProbTest.root")
)

process.p = cms.Path(process.flashggMicroAODSequence * process.egmPhotonIDSequence * process.commissioning)
