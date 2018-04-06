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

void VertexMvaProb()
{
    TString outfileName( "outputTMVA_BDTVtxProb.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    factory->AddVariable( "pt"    , "pt"   , 'F');
    factory->AddVariable( "NVert" , "NVert", 'F');
    factory->AddVariable( "MVA0"  , "MVA0" , 'F');
    factory->AddVariable( "MVA1"  , "MVA1" , 'F');
    factory->AddVariable( "DZ1"   , "DZ1"  , 'F');
    factory->AddVariable( "MVA2"  , "MVA2" , 'F');
    factory->AddVariable( "DZ2"   , "DZ2"  , 'F');
    factory->AddVariable( "NConv" , "NConv", 'F');


    TFile* inputFile1 = TFile::Open("input/GluGluHToGG_M120_13TeV_amcatnloFXFX_pythia8_Training2_s.root");
    TFile* inputFile2 = TFile::Open("input/GluGluHToGG_M120_13TeV_amcatnloFXFX_pythia8_Training2_b.root");
    //TTree* inputTree = (TTree*) inputFile->Get("commissioning/diphoTree");

    //TCut signalCut = "DZtrue < 1";
    //TCut backgrCut = "DZtrue >=1";
    //factory->SetInputTrees( inputTree, signalCut, backgrCut );

    TTree *signal1     = inputFile1->Get("diphoTree");
    TTree *bkg1        = inputFile2->Get("diphoTree");
                                                          
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
                                                          
    factory->AddSignalTree( signal1, signalWeight);
                                                          
    factory->AddBackgroundTree( bkg1, backgroundWeight);

	outputFile->cd();

    TCut mycuts = "MVA2 > -1 && dipho_index==0";
    TCut mycutb = "MVA2 > -1 && dipho_index==0"; 

    factory->PrepareTrainingAndTestTree( mycuts, mycutb,
//                                    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Alternate:NormMode=NumEvents:!V" );
                                    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );


    TMVA::MethodCategory* mcat = 0;

    TString theCat1Vars = "pt:NVert:MVA0:MVA1:DZ1:MVA2:DZ2";
    TString theCat2Vars = "pt:NVert:MVA0:MVA1:DZ1:MVA2:DZ2";
 
    //TMVA::MethodBase* BDT_Cat = factory->BookMethod( TMVA::Types::kCategory, "BDT","" );
    //mcat = dynamic_cast<TMVA::MethodCategory*>(BDT_Cat);

    //mcat->AddMethod( "NConv<1", theCat1Vars, TMVA::Types::kBDT, "0_1_BDTGNewTrue","!H:!V:!CreateMVAPdfs:NTrees=1000:NNodesMax=5:BoostType=Grad:UseBaggedGrad:Shrinkage=0.30:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5");
    //mcat->AddMethod( "NConv>=1",  theCat2Vars, TMVA::Types::kBDT, "1_BDTGNewTrue","!H:!V:!CreateMVAPdfs:NTrees=1000:NNodesMax=5:BoostType=Grad:UseBaggedGrad:Shrinkage=0.30:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5" );

    TMVA::MethodBase* BDT_Cat = factory->BookMethod( TMVA::Types::kCategory, "BDTVtxProb","" );
    mcat = dynamic_cast<TMVA::MethodCategory*>(BDT_Cat);

    mcat->AddMethod( "NConv<1", theCat1Vars, TMVA::Types::kBDT, "BDTVtxProb_noconv","!H:!V:!CreateMVAPdfs:NTrees=1000:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=3:MinNodeSize=2");
    mcat->AddMethod( "NConv>=1",  theCat2Vars, TMVA::Types::kBDT, "BDTVtxProb_conv","!H:!V:!CreateMVAPdfs:NTrees=1000:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=3:MinNodeSize=2");

    factory->TrainAllMethods();

    factory->TestAllMethods();

    factory->EvaluateAllMethods();

   // --------------------------------------------------------------

    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
}
