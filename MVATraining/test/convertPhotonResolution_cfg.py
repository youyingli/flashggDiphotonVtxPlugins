import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("vertexTrainingTreeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")

# geometry and global tag:

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = '94X_mc2017_realistic_v10'

#**************************************************************

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [
       '/store/mc/RunIIFall17MiniAOD/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/A87F7D05-0403-E811-AE12-C0BFC0E567FE.root' ] );

secFiles.extend( [
               ] )

#**************************************************************

process.load("flashgg/MicroAOD/flashggMicroAODSequence_cff")

process.commissioning = cms.EDAnalyzer('convertPhotonResolution',
                                       PhotonTag=cms.untracked.InputTag('flashggPhotons'),
                                       ConversionTag=cms.untracked.InputTag("reducedEgamma","reducedConversions"), 
                                       SingleLegConversionTag=cms.untracked.InputTag("reducedEgamma","reducedSingleLegConversions"), 
                                       BeamSpotTag=cms.untracked.InputTag('offlineBeamSpot'),
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
                                   fileName = cms.string("tests.root")
)

process.p = cms.Path(process.flashggMicroAODSequence*process.egmPhotonIDSequence*process.commissioning)
#process.p = cms.Path(process.flashggMicroAODSequence*process.commissioning)
