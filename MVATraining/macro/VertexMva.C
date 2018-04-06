#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void VertexMva()
{
    TString outfileName( "outputTMVA_BDTVtxId.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    factory->AddVariable("ptasym", "ptasym", "", 'F'); 
    factory->AddVariable("ptbal", "ptbal", "", 'F');
    factory->AddVariable("logsumpt2", "logsumpt2", "", 'F'); 
    factory->AddVariable("limPullToConv", "limPullToConv", "", 'F');
    factory->AddVariable("nConv", "nConv", "", 'F');

	TFile* inputFile_ggf = TFile::Open("vtxid_training_samples/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8_94X_TrainV1.root");
	TFile* inputFile_vbf = TFile::Open("vtxid_training_samples/VBFHToGG_M125_13TeV_amcatnlo_pythia8_94X_TrainV1.root");
	TFile* inputFile_vh  = TFile::Open("vtxid_training_samples/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_94X_TrainV1.root");
	TFile* inputFile_tth = TFile::Open("vtxid_training_samples/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_94X_TrainV1.root");

  	TTree* signal_ggf  = (TTree*)inputFile_ggf->Get("commissioning/signalTree");
  	TTree* signal_vbf  = (TTree*)inputFile_vbf->Get("commissioning/signalTree");
  	TTree* signal_vh   = (TTree*)inputFile_vh ->Get("commissioning/signalTree");
  	TTree* signal_tth  = (TTree*)inputFile_tth->Get("commissioning/signalTree");
    TTree* bkg_ggf     = (TTree*)inputFile_ggf->Get("commissioning/backgroundTree");
    TTree* bkg_vbf     = (TTree*)inputFile_vbf->Get("commissioning/backgroundTree");
    TTree* bkg_vh      = (TTree*)inputFile_vh ->Get("commissioning/backgroundTree");
    TTree* bkg_tth     = (TTree*)inputFile_tth->Get("commissioning/backgroundTree");

	Double_t signalWeight_ggf = 1.0, signalWeight_vbf = 0.0779, signalWeight_vh = 0.0465, signalWeight_tth = 0.0104;
	Double_t backgroundWeight_ggf = 1.0, backgroundWeight_vbf = 0.0779, backgroundWeight_vh = 0.0465, backgroundWeight_tth = 0.0104;

	factory->AddSignalTree( signal_ggf, signalWeight_ggf );
	factory->AddSignalTree( signal_vbf, signalWeight_vbf );
	factory->AddSignalTree( signal_vh , signalWeight_vh  );
	factory->AddSignalTree( signal_tth, signalWeight_tth );

	factory->AddBackgroundTree( bkg_ggf, backgroundWeight_ggf );
	factory->AddBackgroundTree( bkg_vbf, backgroundWeight_vbf );
	factory->AddBackgroundTree( bkg_vh , backgroundWeight_vh  );
	factory->AddBackgroundTree( bkg_tth, backgroundWeight_tth );

    factory->SetWeightExpression( "genweight" );

	outputFile->cd();

    TCut mycuts = "ptbal < 500 && logsumpt2 >= -10 && dipho_index==0";
    TCut mycutb = "ptbal < 500 && logsumpt2 >= -10 && dipho_index==0"; 

    factory->PrepareTrainingAndTestTree( mycuts, mycutb,
//                                    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Alternate:NormMode=NumEvents:!V" );
//                                    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
//                                    "nTrain_Signal=200000:nTrain_Background=200000:SplitMode=Random:NormMode=NumEvents:!V" );
                                    "nTrain_Signal=200000:nTrain_Background=200000:nTest_Signal=200000:nTest_Background=200000:SplitMode=Random:NormMode=NumEvents:!V" );

    TMVA::MethodCategory* mcat = 0;

    TString theCat1Vars = "ptasym:ptbal:logsumpt2";
    TString theCat2Vars = "ptasym:ptbal:logsumpt2:limPullToConv";
 
    //TMVA::MethodBase* BDT_Cat = factory->BookMethod( TMVA::Types::kCategory, "BDT","" );
    //mcat = dynamic_cast<TMVA::MethodCategory*>(BDT_Cat);

    //mcat->AddMethod( "NConv<1", theCat1Vars, TMVA::Types::kBDT, "0_1_BDTGNewTrue","!H:!V:!CreateMVAPdfs:NTrees=1000:NNodesMax=5:BoostType=Grad:UseBaggedGrad:Shrinkage=0.30:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5");
    //mcat->AddMethod( "NConv>=1",  theCat2Vars, TMVA::Types::kBDT, "1_BDTGNewTrue","!H:!V:!CreateMVAPdfs:NTrees=1000:NNodesMax=5:BoostType=Grad:UseBaggedGrad:Shrinkage=0.30:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5" );

    TMVA::MethodBase* BDT_Cat = factory->BookMethod( TMVA::Types::kCategory, "BDTVtxId","" );
    mcat = dynamic_cast<TMVA::MethodCategory*>(BDT_Cat);

    mcat->AddMethod( "nConv<1", theCat1Vars, TMVA::Types::kBDT, "BDTVtxId_noconv","!H:!V:!CreateMVAPdfs:NTrees=1000:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=3:MinNodeSize=2:NegWeightTreatment=ignorenegweightsintraining");
    mcat->AddMethod( "nConv>=1",  theCat2Vars, TMVA::Types::kBDT, "BDTVtxId_conv","!H:!V:!CreateMVAPdfs:NTrees=1000:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=3:MinNodeSize=2:NegWeightTreatment=ignorenegweightsintraining");

    factory->TrainAllMethods();

    factory->TestAllMethods();

    factory->EvaluateAllMethods();

   // --------------------------------------------------------------

    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
}
