
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//for miniAOD:
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TTree.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees


// **********************************************************************

using namespace flashgg;

using namespace std;
using namespace edm;
using namespace reco;
using namespace math;

struct DimuonInfo {

    float BSz             ;
    float BSsigmaz        ;
    float mcWeight        ;
    float zRecoTrue       ;    
    float zRecoTrueNoMu   ;
    float zTrue           ;
    int   nPU             ;
    float dimuon_mass     ; 
    float dimuon_pt       ; 
    float muon1_pt        ; 
    float muon1_eta       ; 
    float muon2_pt        ; 
    float muon2_eta       ; 
    float zChosenWithMu   ; 
    float zChosenNoMu     ; 
    float zRecoFromMu     ;
    float zRecoFromMuNoMu ;

    float mvaProbWithMu   ; 
    float logSumPt2WithMu ; 
    float ptBalWithMu     ; 
    float ptAsymWithMu    ; 
    int   nVtxWithMu      ; 
    float mva0WithMu      ; 
    float mva1WithMu      ; 
    float mva2WithMu      ; 
    float dZ1WithMu       ; 
    float dZ2WithMu       ;

    float mvaProbNoMu     ; 
    float logSumPt2NoMu   ; 
    float ptBalNoMu       ; 
    float ptAsymNoMu      ; 
    int   nVtxNoMu        ; 
    float mva0NoMu        ; 
    float mva1NoMu        ; 
    float mva2NoMu        ; 
    float dZ1NoMu         ; 
    float dZ2NoMu         ; 

};

struct SignalInfo {

    float vtxIdmva;
    float LogSumPt2;
    float PtBal;
    float PtAsym;
    float DZtrue;

};

struct BackgroundInfo {

    float vtxIdmva;
    float LogSumPt2;
    float PtBal;
    float PtAsym;
    float DZtrue;

};

// **********************************************************************

class outputAnalyzerDiMuHighMass : public edm::EDAnalyzer
{
public:
    explicit outputAnalyzerDiMuHighMass( const edm::ParameterSet & );
    ~outputAnalyzerDiMuHighMass();
    
    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );
    
    
private:
    
    edm::Service<TFileService> fs_;
    
    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    void initEventStructure();
   
    float getMCTruthVertexZ( const std::vector<edm::Ptr<reco::GenParticle>>& gens );
    int getRecoClosestToTrueVertexIndex(  const std::vector<edm::Ptr<reco::GenParticle>>& gens, const vector<edm::Ptr<reco::Vertex> > & vertices, double dzMatch = 50. );
    int getRecoWithMuonsVertexIndex( const vector<edm::Ptr<reco::Vertex> > & vertices, Ptr<pat::Muon> mu1, Ptr<pat::Muon> mu2, double dzMatch = 0.2 );
    int sortedIndex( const unsigned int trueVtxIndex, const std::vector<int> VtxSortedIndexVector);
    std::vector<int> isTightMuonWithoutVtxCut(const vector<edm::Ptr<pat::Muon>> & patMuons); 

    edm::EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    edm::EDGetToken genEventInfoToken_;
    edm::EDGetTokenT<View<pat::Muon> > patMuonToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexTokenNoMu_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexTokenWithMu_;
    EDGetTokenT< VertexCandidateMap > vertexCandidateMapToken_;
    EDGetTokenT< VertexCandidateMap > vertexCandidateMapNoMuToken_;
    unique_ptr<VertexSelectorBase> vertexSelector_;
    edm::EDGetTokenT<reco::BeamSpot > beamSpotToken_;
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> >  PileUpToken_;

    TTree *vtxTree;
    TTree *signalTree;
    TTree *backgroundTree;

    DimuonInfo dimuonInfo;
    SignalInfo sigInfo;
    BackgroundInfo bkgInfo;
};

// ******************************************************************************************


//
// constructors and destructor
//
outputAnalyzerDiMuHighMass::outputAnalyzerDiMuHighMass( const edm::ParameterSet &iConfig ):
    genParticleToken_(consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag> ("GenParticleTag", InputTag("prunedGenParticles")))),
    genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter <edm::InputTag> ("genEventInfoProduct"))),
    patMuonToken_( consumes<View<pat::Muon> >( iConfig.getParameter<InputTag>( "patMuonTag" ) ) ),
    vertexTokenNoMu_( consumes<View<reco::Vertex> >( iConfig.getUntrackedParameter<InputTag> ( "VertexTagNoMu", InputTag( "offlinePrimaryVerticesNoMu" ) ) ) ),
    vertexTokenWithMu_( consumes<View<reco::Vertex> >( iConfig.getUntrackedParameter<InputTag> ( "VertexTagWithMu", InputTag( "offlinePrimaryVerticesWithMu" ) ) ) ),
    vertexCandidateMapToken_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTag" ) ) ),
    vertexCandidateMapNoMuToken_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTagNoMu" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getUntrackedParameter<InputTag>( "BeamSpotTag", InputTag( "offlineBeamSpot" ) ) ) ),
    PileUpToken_( consumes<View<PileupSummaryInfo> >( iConfig.getUntrackedParameter<InputTag> ( "PileUpTag", InputTag( "slimmedAddPileupInfo" ) ) ) )
{
        const std::string &VertexSelectorName = iConfig.getParameter<std::string>( "VertexSelectorName" );
        vertexSelector_.reset( FlashggVertexSelectorFactory::get()->create( VertexSelectorName, iConfig ) );
}


outputAnalyzerDiMuHighMass::~outputAnalyzerDiMuHighMass()
{


}



void
outputAnalyzerDiMuHighMass::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{
    Handle<View<reco::GenParticle> > genParticles;
    iEvent.getByToken(genParticleToken_,genParticles);

    Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(genEventInfoToken_,genEventInfo);

    Handle<View<pat::Muon> >  patMuons;
    iEvent.getByToken( patMuonToken_, patMuons );

    Handle<View<reco::Vertex> > primaryVerticesNoMu;
    iEvent.getByToken( vertexTokenNoMu_, primaryVerticesNoMu );

    Handle<View<reco::Vertex> > primaryVerticesWithMu;
    iEvent.getByToken( vertexTokenWithMu_, primaryVerticesWithMu );

    Handle<VertexCandidateMap> vertexCandidateMap;
    iEvent.getByToken( vertexCandidateMapToken_, vertexCandidateMap );

    Handle<VertexCandidateMap> vertexCandidateMapNoMu;
    iEvent.getByToken( vertexCandidateMapNoMuToken_, vertexCandidateMapNoMu );

    Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken( beamSpotToken_, recoBeamSpotHandle );
    math::XYZPoint beamSpot;
    double BSsigmaz = 0;
    if( recoBeamSpotHandle.isValid() ) {
        beamSpot = recoBeamSpotHandle->position();
        BSsigmaz = recoBeamSpotHandle->sigmaZ();

    } else {
        cout << " WARNING NO VALID BEAM SPOT: this should not happen!" << endl;
    }

    Handle<View< PileupSummaryInfo> > PileupInfos;
    iEvent.getByToken( PileUpToken_, PileupInfos ); 

    initEventStructure();

    if( !iEvent.isRealData() ) {

        dimuonInfo.mcWeight = genEventInfo->weight();

        int trueRecoVtxIdx = getRecoClosestToTrueVertexIndex( genParticles->ptrs(), primaryVerticesWithMu->ptrs());
        int trueRecoVtxIdxNoMu = getRecoClosestToTrueVertexIndex( genParticles->ptrs(), primaryVerticesNoMu->ptrs());

        if(trueRecoVtxIdx != -1)     dimuonInfo.zRecoTrue     = primaryVerticesWithMu->ptrAt( trueRecoVtxIdx )->z();    
        if(trueRecoVtxIdxNoMu != -1) dimuonInfo.zRecoTrueNoMu = primaryVerticesNoMu->ptrAt( trueRecoVtxIdxNoMu )->z();
        dimuonInfo.zTrue = getMCTruthVertexZ( genParticles->ptrs() );
        

        float nPu = -999.;
        
        // pileup info
        for( unsigned int PVI = 0; PVI < PileupInfos->size(); ++PVI ) {
            Int_t pu_bunchcrossing = PileupInfos->ptrAt( PVI )->getBunchCrossing();
            if( pu_bunchcrossing == 0 ) {
                nPu = PileupInfos->ptrAt( PVI )->getTrueNumInteractions();
                break;
            }
        }
        
        dimuonInfo.nPU = nPu;
    }
    
    
    vector<int> muon_index = isTightMuonWithoutVtxCut(patMuons->ptrs());
  
    if( muon_index.size() >=2 ){
        
        float Zmass = 91.19;
        Ptr<pat::Muon> pat_muon1 = patMuons->ptrAt(muon_index[0]);
        Ptr<pat::Muon> pat_muon2 = patMuons->ptrAt(muon_index[1]);
   
        XYZTLorentzVector dimuonP4 = pat_muon1->p4() + pat_muon2->p4();
        double dimuon_mass = dimuonP4.M();
        double dimuon_pt = dimuonP4.pt();

        if ( fabs(dimuon_mass - Zmass) < 30. ) {

            dimuonInfo.dimuon_mass = dimuon_mass;
            dimuonInfo.dimuon_pt = dimuon_pt;

            dimuonInfo.BSz = beamSpot.z();
            dimuonInfo.BSsigmaz = BSsigmaz;

            dimuonInfo.muon1_pt  = pat_muon1->pt();
            dimuonInfo.muon1_eta = pat_muon1->eta();
            dimuonInfo.muon2_pt  = pat_muon1->pt();
            dimuonInfo.muon2_eta = pat_muon1->eta();

            //Selected Vertex
            vector<float> infoWithMu;
            vector<float> infoNoMu;

            Ptr<reco::Vertex> selectedVtxWithMu = vertexSelector_->select( pat_muon1, pat_muon2, primaryVerticesWithMu->ptrs(), *vertexCandidateMap, beamSpot ); 
            vertexSelector_->getInfoFromLastSelection(infoWithMu);

            Ptr<reco::Vertex> selectedVtxNoMu = vertexSelector_->select( pat_muon1, pat_muon2, primaryVerticesNoMu->ptrs(), *vertexCandidateMapNoMu, beamSpot ); 
            vertexSelector_->getInfoFromLastSelection(infoNoMu);
            
            dimuonInfo.zChosenWithMu = selectedVtxWithMu->z();
            dimuonInfo.zChosenNoMu = selectedVtxNoMu->z();

            dimuonInfo.mvaProbWithMu   = infoWithMu[0]; 
            dimuonInfo.logSumPt2WithMu = infoWithMu[1];
            dimuonInfo.ptBalWithMu     = infoWithMu[2];
            dimuonInfo.ptAsymWithMu    = infoWithMu[3];
            dimuonInfo.nVtxWithMu      = infoWithMu[4];
            dimuonInfo.mva0WithMu      = infoWithMu[5];
            dimuonInfo.mva1WithMu      = infoWithMu[6];
            dimuonInfo.mva2WithMu      = infoWithMu[7];
            dimuonInfo.dZ1WithMu       = infoWithMu[8];
            dimuonInfo.dZ2WithMu       = infoWithMu[9];

            dimuonInfo.mvaProbNoMu     = infoNoMu[0]; 
            dimuonInfo.logSumPt2NoMu   = infoNoMu[1];
            dimuonInfo.ptBalNoMu       = infoNoMu[2];
            dimuonInfo.ptAsymNoMu      = infoNoMu[3];
            dimuonInfo.nVtxNoMu        = infoNoMu[4];
            dimuonInfo.mva0NoMu        = infoNoMu[5];
            dimuonInfo.mva1NoMu        = infoNoMu[6];
            dimuonInfo.mva2NoMu        = infoNoMu[7];
            dimuonInfo.dZ1NoMu         = infoNoMu[8];
            dimuonInfo.dZ2NoMu         = infoNoMu[9];


            //Reconstruction real vertex
            int iFromMuonsWithMu = getRecoWithMuonsVertexIndex(primaryVerticesWithMu->ptrs(), pat_muon1, pat_muon2);
            int iFromMuonsNoMu = getRecoWithMuonsVertexIndex(primaryVerticesNoMu->ptrs(), pat_muon1, pat_muon2);       
      
            if( iFromMuonsNoMu != -1 && iFromMuonsWithMu != -1 ) {
                dimuonInfo.zRecoFromMu = primaryVerticesWithMu->ptrAt(iFromMuonsWithMu)->z();
                dimuonInfo.zRecoFromMuNoMu = primaryVerticesNoMu->ptrAt(iFromMuonsNoMu)->z();
                vtxTree->Fill();

                //Right Vertex
                vector<int> VtxSortedIndexVector;
                for( int i = 0; i < (int)primaryVerticesNoMu->size(); i++) VtxSortedIndexVector.push_back( vertexSelector_->getSortedIndexFromLastSelection( i ) );
                unsigned int trueVtxIndex = iFromMuonsNoMu;
                int trueVtxSortedIndex = sortedIndex(trueVtxIndex, VtxSortedIndexVector);

                if (trueVtxSortedIndex < 0) {    
                    vector<float> info;
                    vertexSelector_->getInfoFromLastSelectionForVtxIdx( info , (unsigned int) trueVtxSortedIndex );
                    sigInfo.vtxIdmva  = info[0];
                    sigInfo.LogSumPt2 = info[1];
                    sigInfo.PtBal     = info[2];
                    sigInfo.PtAsym    = info[3];
                    sigInfo.DZtrue    = primaryVerticesNoMu->ptrAt(trueVtxIndex)->position().z() - primaryVerticesWithMu->ptrAt(iFromMuonsWithMu)->position().z();

                    signalTree->Fill();
                }

                //Wrong Vertex
                vector<int>	notrueVtxIndexVector;
                for( unsigned int i = 0 ; i < primaryVerticesNoMu->size() ; i++ ) {
                    if( i != trueVtxIndex ) { notrueVtxIndexVector.push_back( i ); }
                }

                int irand = -999;
                if( notrueVtxIndexVector.size() > 1 ) { irand = rand() % notrueVtxIndexVector.size(); }

                int randVtxIndex = -999;
                if( irand != -999 ) { randVtxIndex = notrueVtxIndexVector[irand]; }

                int randVtxSortedIndexI = sortedIndex(randVtxIndex, VtxSortedIndexVector);

                if( randVtxSortedIndexI != -1 ) {
                    unsigned int randVtxSortedIndex = randVtxSortedIndexI;

                    vector<float> info;
                    vertexSelector_->getInfoFromLastSelectionForVtxIdx( info , randVtxSortedIndex );
                    bkgInfo.vtxIdmva  = info[0];
                    bkgInfo.LogSumPt2 = info[1];
                    bkgInfo.PtBal     = info[2];
                    bkgInfo.PtAsym    = info[3];
                    bkgInfo.DZtrue    = primaryVerticesNoMu->ptrAt(randVtxIndex)->position().z() - primaryVerticesWithMu->ptrAt(iFromMuonsWithMu)->position().z();

                    backgroundTree->Fill();
                }
            }

        }//mass cut
    }//tight muon

}


void
outputAnalyzerDiMuHighMass::beginJob()
{

    vtxTree = fs_->make<TTree>( "vtxTree", "vtxTree" );   
    vtxTree->Branch( "BSz"              , &dimuonInfo.BSz              , "BSz/F"             );  
    vtxTree->Branch( "BSsigmaz"         , &dimuonInfo.BSsigmaz         , "BSsigmaz/F"        );  
    vtxTree->Branch( "mcWeight"         , &dimuonInfo.mcWeight         , "mcWeight/F"        );  
    vtxTree->Branch( "zRecoTrue"        , &dimuonInfo.zRecoTrue        , "zRecoTrue/F"       );  
    vtxTree->Branch( "zRecoTrueNoMu"    , &dimuonInfo.zRecoTrueNoMu    , "zRecoTrueNoMu/F"   );  
    vtxTree->Branch( "zTrue"            , &dimuonInfo.zTrue            , "zTrue/F"           );  
    vtxTree->Branch( "nPU"              , &dimuonInfo.nPU              , "nPU/I"             );  
    vtxTree->Branch( "dimuon_mass"      , &dimuonInfo.dimuon_mass      , "dimuon_mass/F"     );  
    vtxTree->Branch( "dimuon_pt"        , &dimuonInfo.dimuon_pt        , "dimuon_pt/F"       );  
    vtxTree->Branch( "muon1_pt"         , &dimuonInfo.muon1_pt         , "muon1_pt/F"        );  
    vtxTree->Branch( "muon1_eta"        , &dimuonInfo.muon1_eta        , "muon1_eta/F"       );  
    vtxTree->Branch( "muon2_pt"         , &dimuonInfo.muon2_pt         , "muon2_pt/F"        );  
    vtxTree->Branch( "muon2_eta"        , &dimuonInfo.muon2_eta        , "muon2_eta/F"       );  
    vtxTree->Branch( "zChosenWithMu"    , &dimuonInfo.zChosenWithMu    , "zChosenWithMu/F"   );  
    vtxTree->Branch( "zChosenNoMu"      , &dimuonInfo.zChosenNoMu      , "zChosenNoMu/F"     );  
    vtxTree->Branch( "zRecoFromMu"      , &dimuonInfo.zRecoFromMu      , "zRecoFromMu/F"     );  
    vtxTree->Branch( "zRecoFromMuNoMu"  , &dimuonInfo.zRecoFromMuNoMu  , "zRecoFromMuNoMu/F" );  
    vtxTree->Branch( "mvaProbWithMu"    , &dimuonInfo.mvaProbWithMu    , "mvaProbWithMu/F"   );  
    vtxTree->Branch( "logSumPt2WithMu"  , &dimuonInfo.logSumPt2WithMu  , "logSumPt2WithMu/F" );  
    vtxTree->Branch( "ptBalWithMu"      , &dimuonInfo.ptBalWithMu      , "ptBalWithMu/F"     );  
    vtxTree->Branch( "ptAsymWithMu"     , &dimuonInfo.ptAsymWithMu     , "ptAsymWithMu/F"    );  
    vtxTree->Branch( "nVtxWithMu"       , &dimuonInfo.nVtxWithMu       , "nVtxWithMu/I"      );  
    vtxTree->Branch( "mva0WithMu"       , &dimuonInfo.mva0WithMu       , "mva0WithMu/F"      );  
    vtxTree->Branch( "mva1WithMu"       , &dimuonInfo.mva1WithMu       , "mva1WithMu/F"      );  
    vtxTree->Branch( "mva2WithMu"       , &dimuonInfo.mva2WithMu       , "mva2WithMu/F"      );  
    vtxTree->Branch( "dZ1WithMu"        , &dimuonInfo.dZ1WithMu        , "dZ1WithMu/F"       );  
    vtxTree->Branch( "dZ2WithMu"        , &dimuonInfo.dZ2WithMu        , "dZ2WithMu/F"       );  
    vtxTree->Branch( "mvaProbNoMu"      , &dimuonInfo.mvaProbNoMu      , "mvaProbNoMu/F"     );  
    vtxTree->Branch( "logSumPt2NoMu"    , &dimuonInfo.logSumPt2NoMu    , "logSumPt2NoMu/F"   );  
    vtxTree->Branch( "ptBalNoMu"        , &dimuonInfo.ptBalNoMu        , "ptBalNoMu/F"       );  
    vtxTree->Branch( "ptAsymNoMu"       , &dimuonInfo.ptAsymNoMu       , "ptAsymNoMu/F"      );  
    vtxTree->Branch( "nVtxNoMu"         , &dimuonInfo.nVtxNoMu         , "nVtxNoMu/I"        );  
    vtxTree->Branch( "mva0NoMu"         , &dimuonInfo.mva0NoMu         , "mva0NoMu/F"        );  
    vtxTree->Branch( "mva1NoMu"         , &dimuonInfo.mva1NoMu         , "mva1NoMu/F"        );  
    vtxTree->Branch( "mva2NoMu"         , &dimuonInfo.mva2NoMu         , "mva2NoMu/F"        );  
    vtxTree->Branch( "dZ1NoMu"          , &dimuonInfo.dZ1NoMu          , "dZ1NoMu/F"         );  
    vtxTree->Branch( "dZ2NoMu"          , &dimuonInfo.dZ2NoMu          , "dZ2NoMu/F"         );  

    signalTree = fs_->make<TTree>( "signalTree", "per-diphoton tree" );
    signalTree->Branch( "BSz"        , &dimuonInfo.BSz       , "BSz/F"      );
    signalTree->Branch( "BSsigmaz"   , &dimuonInfo.BSsigmaz  , "BSsigmaz/F" );
    signalTree->Branch( "mcWeight"   , &dimuonInfo.mcWeight  , "mcWeight/F" );
    signalTree->Branch( "nPU"        , &dimuonInfo.nPU       , "nPU/I"      );
    signalTree->Branch( "nvertex"    , &dimuonInfo.nVtxNoMu  , "nvertex/I"  );
    signalTree->Branch( "vtxIdmva"   , &sigInfo.vtxIdmva     , "vtxIdmva/F" );
    signalTree->Branch( "LogSumPt2"  , &sigInfo.LogSumPt2    , "LogSumPt2/F");
    signalTree->Branch( "PtBal"      , &sigInfo.PtBal        , "PtBal/F"    );
    signalTree->Branch( "PtAsym"     , &sigInfo.PtAsym       , "PtAsym/F"   );
    signalTree->Branch( "DZtrue"     , &sigInfo.DZtrue       , "DZtrue/F"   );

    backgroundTree = fs_->make<TTree>( "backgroundTree", "per-diphoton tree" );
    backgroundTree->Branch( "BSz"        , &dimuonInfo.BSz       , "BSz/F"      );
    backgroundTree->Branch( "BSsigmaz"   , &dimuonInfo.BSsigmaz  , "BSsigmaz/F" );
    backgroundTree->Branch( "mcWeight"   , &dimuonInfo.mcWeight  , "mcWeight/F" );
    backgroundTree->Branch( "nPU"        , &dimuonInfo.nPU       , "nPU/I"      );
    backgroundTree->Branch( "nvertex"    , &dimuonInfo.nVtxNoMu  , "nvertex/I"  );
    backgroundTree->Branch( "vtxIdmva"   , &bkgInfo.vtxIdmva     , "vtxIdmva/F" );
    backgroundTree->Branch( "LogSumPt2"  , &bkgInfo.LogSumPt2    , "LogSumPt2/F");
    backgroundTree->Branch( "PtBal"      , &bkgInfo.PtBal        , "PtBal/F"    );
    backgroundTree->Branch( "PtAsym"     , &bkgInfo.PtAsym       , "PtAsym/F"   );
    backgroundTree->Branch( "DZtrue"     , &bkgInfo.DZtrue       , "DZtrue/F"   );

}

void
outputAnalyzerDiMuHighMass::endJob()
{
}

void
outputAnalyzerDiMuHighMass::initEventStructure()
{
    dimuonInfo.BSz             = -999.; 
    dimuonInfo.BSsigmaz        = -999.; 
    dimuonInfo.mcWeight        = -999.; 
    dimuonInfo.zRecoTrue       = -999.; 
    dimuonInfo.zRecoTrueNoMu   = -999.; 
    dimuonInfo.zTrue           = -999.; 
    dimuonInfo.nPU             = -999 ; 
    dimuonInfo.dimuon_mass     = -999.; 
    dimuonInfo.dimuon_pt       = -999.; 
    dimuonInfo.muon1_pt        = -999.; 
    dimuonInfo.muon1_eta       = -999.; 
    dimuonInfo.muon2_pt        = -999.; 
    dimuonInfo.muon2_eta       = -999.; 
    dimuonInfo.zChosenWithMu   = -999.; 
    dimuonInfo.zChosenNoMu     = -999.; 
    dimuonInfo.zRecoFromMu     = -999.; 
    dimuonInfo.zRecoFromMuNoMu = -999.; 
    dimuonInfo.mvaProbWithMu   = -999.; 
    dimuonInfo.logSumPt2WithMu = -999.; 
    dimuonInfo.ptBalWithMu     = -999.; 
    dimuonInfo.ptAsymWithMu    = -999.; 
    dimuonInfo.nVtxWithMu      = -999 ; 
    dimuonInfo.mva0WithMu      = -999.; 
    dimuonInfo.mva1WithMu      = -999.; 
    dimuonInfo.mva2WithMu      = -999.; 
    dimuonInfo.dZ1WithMu       = -999.; 
    dimuonInfo.dZ2WithMu       = -999.; 
    dimuonInfo.mvaProbNoMu     = -999.; 
    dimuonInfo.logSumPt2NoMu   = -999.; 
    dimuonInfo.ptBalNoMu       = -999.; 
    dimuonInfo.ptAsymNoMu      = -999.; 
    dimuonInfo.nVtxNoMu        = -999 ; 
    dimuonInfo.mva0NoMu        = -999.; 
    dimuonInfo.mva1NoMu        = -999.; 
    dimuonInfo.mva2NoMu        = -999.; 
    dimuonInfo.dZ1NoMu         = -999.; 
    dimuonInfo.dZ2NoMu         = -999.; 

    sigInfo.vtxIdmva   = -999.; 
    sigInfo.LogSumPt2  = -999.; 
    sigInfo.PtBal      = -999.; 
    sigInfo.PtAsym     = -999.; 
    sigInfo.DZtrue     = -999.; 

    bkgInfo.vtxIdmva   = -999.; 
    bkgInfo.LogSumPt2  = -999.; 
    bkgInfo.PtBal      = -999.; 
    bkgInfo.PtAsym     = -999.; 
    bkgInfo.DZtrue     = -999.; 

}

int outputAnalyzerDiMuHighMass::getRecoWithMuonsVertexIndex( const vector<edm::Ptr<reco::Vertex> > & vertices, Ptr<pat::Muon> mu1, Ptr<pat::Muon> mu2, double dzMatch)
{

    double IP1 = mu1->vz();
    double IP2 = mu2->vz(); 

    double average = 0.5 * (IP1 + IP2);
        
    int  ivMatch = 0;
    int  ivMatch1 = 0;
    int  ivMatch2 = 0;
    double dzMin = 999;
    double dzMin1 = 999;
    double dzMin2 = 999;
    
    for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {
        double dz = fabs( vertices[iv]->z() - average );
        double dz1 = fabs( vertices[iv]->z() - IP1 );
        double dz2 = fabs( vertices[iv]->z() - IP2 );
        if( dz < dzMin ) {
            ivMatch = iv;
            dzMin   = dz;
        } 
        if( dz1 < dzMin1 ) {
            ivMatch1 = iv;
            dzMin1   = dz1;
        }
        if( dz2 < dzMin2 ) {
            ivMatch2 = iv;
            dzMin2   = dz2;
        }
    }
     
    if(ivMatch == ivMatch1 && ivMatch == ivMatch2 && dzMin < dzMatch) return ivMatch;
    return -1;
}


float 
outputAnalyzerDiMuHighMass::getMCTruthVertexZ( const std::vector<edm::Ptr<reco::GenParticle>>& gens )
{
    float zTrue = 999.;
    for( unsigned int genLoop = 0 ; genLoop < gens.size(); genLoop++ ) {
        if( fabs( gens[genLoop]->pdgId() ) == 23 ) {
            zTrue = gens[genLoop]->vz();
            break;
        }
    }
    return zTrue;
}

int outputAnalyzerDiMuHighMass::getRecoClosestToTrueVertexIndex(  const std::vector<edm::Ptr<reco::GenParticle>>& gens , const vector<edm::Ptr<reco::Vertex> > & vertices, double dzMatch )
{

    reco::Vertex::Point hardVertex( 0, 0, 0 );

    for( unsigned int genLoop = 0 ; genLoop < gens.size(); genLoop++ ) {

        if( fabs( gens[genLoop]->pdgId() ) < 10 || fabs( gens[genLoop]->pdgId() ) == 23 ) {
            hardVertex.SetCoordinates( gens[genLoop]->vx(), gens[genLoop]->vy(), gens[genLoop]->vz() );
            break;
        }
    }

    int  ivMatch = 0;
    double dzMin = 999;
    for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {
        double dz = fabs( vertices[iv]->z() - hardVertex.z() );
        if( dz < dzMin ) {
            ivMatch = iv;
            dzMin   = dz;
        }
    }


    if( dzMin < dzMatch ) { return ivMatch; }

    return -1;
}

int
outputAnalyzerDiMuHighMass::sortedIndex( const unsigned int trueVtxIndex, const vector<int> VtxSortedIndexVector )
{
    for (unsigned int i = 0; i < VtxSortedIndexVector.size(); i++) {
        int index = VtxSortedIndexVector[i];
        if (index < 0) continue;
        if ( (unsigned int) index == trueVtxIndex ) return i;
    }
    return -1;
}

std::vector<int>
outputAnalyzerDiMuHighMass::isTightMuonWithoutVtxCut(const vector<edm::Ptr<pat::Muon>> & patMuons) 
{
    vector<int> muon_index;
    for ( unsigned int i = 0; i < patMuons.size(); i++) {
        Ptr<pat::Muon> pat_muon = patMuons[i];
        if(!pat_muon->isLooseMuon() || pat_muon->pt()<10. ) continue;
        if(!pat_muon->innerTrack().isNonnull()) continue;
        if(!pat_muon->globalTrack().isNonnull()) continue;
        if(pat_muon->globalTrack()->normalizedChi2() > 10. ) continue;
        if(pat_muon->globalTrack()->hitPattern().numberOfValidMuonHits() <= 0) continue;
        if(pat_muon->numberOfMatchedStations() <= 1) continue;
        if(pat_muon->innerTrack()->hitPattern().numberOfValidPixelHits() <= 0) continue;
        if(pat_muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= 5) continue;
        if(pat_muon->isolationR03().sumPt/pat_muon->pt()>0.05) continue;
        muon_index.push_back(i);
    }

    return muon_index;
}

void
outputAnalyzerDiMuHighMass::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}


DEFINE_FWK_MODULE( outputAnalyzerDiMuHighMass );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

