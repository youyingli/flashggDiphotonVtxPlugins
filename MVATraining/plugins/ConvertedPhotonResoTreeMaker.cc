
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
//#include "DataFormats/Common/interface/PtrVector.h"

//#include "DataFormats/VertexReco/interface/Vertex.h"

#include "flashgg/DataFormats/interface/Photon.h"

#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees

struct resolutionInfo {

    int ntracks;
    int nphotons;
    int photon_index;

    double dipho_mass;
    double pho_pt;
    double ztrue;
    double zvtx;

    double z1Pix;
    double z1Tib;
    double z1Tob;
    double z1PixFwd;
    double z1Tid;
    double z1Tec;
    double z2Pix;
    double z2Tib;
    double z2Tob;
    double z2PixFwd;
    double z2Tid;
    double z2Tec;

};

// **********************************************************************

using namespace std;
using namespace edm;
using namespace flashgg;


// **********************************************************************

class ConvertedPhotonResoTreeMaker : public edm::EDAnalyzer
{
public:
    explicit ConvertedPhotonResoTreeMaker( const edm::ParameterSet & );
    ~ConvertedPhotonResoTreeMaker();

    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


private:

    edm::Service<TFileService> fs_;



    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    void initEventStructure();
    std::vector<double> ZmcTruthVertexAndRecoVectex( const std::vector<edm::Ptr<reco::GenParticle> > &genParticles ,
            const std::vector<edm::Ptr<reco::Vertex> > &vertices, double dzMatch = 0.1);
    double vtxZFromConvOnly( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco:: Conversion> &conversion, const math::XYZPoint &beamSpot ) const;
    double vtxZFromConvSuperCluster( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco:: Conversion> &conversion, const math::XYZPoint &beamSpot ) const;
    vector<int> IndexMatchedConversion( const edm::Ptr<flashgg::Photon> &g, 
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,  const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
            bool useSingleLeg ) const;



    edm::EDGetTokenT<edm::View<flashgg::Photon> >       photonToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >          vertexToken_;
    edm::EDGetTokenT<View<reco::Conversion> >           conversionToken_;
    edm::EDGetTokenT<View<reco::Conversion> >           singlelegconversionToken_;
    edm::EDGetTokenT<reco::BeamSpot >                   beamSpotToken_;
    edm::EDGetTokenT<double>                            rhoTaken_;
    edm::EDGetTokenT<View<reco::GenParticle> >          genParticleToken_;


    TTree *resolutionTree;
    resolutionInfo resInfo;
};

// ******************************************************************************************


//
// constructors and destructor
//
ConvertedPhotonResoTreeMaker::ConvertedPhotonResoTreeMaker( const edm::ParameterSet &iConfig ):
    photonToken_(consumes<View<flashgg::Photon> >(iConfig.getUntrackedParameter<InputTag> ("PhotonTag", InputTag("flashggPhotons")))),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getUntrackedParameter<InputTag> ( "VertexTag", InputTag( "offlineSlimmedPrimaryVertices" ) ) ) ),
    conversionToken_(consumes<View<reco::Conversion> >(iConfig.getUntrackedParameter<InputTag>("ConversionTag",InputTag("reducedConversions")))),
    singlelegconversionToken_(consumes<View<reco::Conversion> >(iConfig.getUntrackedParameter<InputTag>("SingleLegConversionTag",InputTag("reducedSingleLegConversions")))),
    beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getUntrackedParameter<InputTag>( "BeamSpotTag", InputTag( "offlineBeamSpot" ) ) ) ),
    rhoTaken_( consumes<double>( iConfig.getUntrackedParameter<InputTag>( "rhoTag", InputTag( "fixedGridRhoFastjetAll" ) ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getUntrackedParameter<InputTag> ( "GenParticleTag", InputTag( "prunedGenParticles" ) ) ) )
{
}


ConvertedPhotonResoTreeMaker::~ConvertedPhotonResoTreeMaker()
{
}

void
ConvertedPhotonResoTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    // ********************************************************************************

    // access edm objects
    
    Handle<View<flashgg::Photon> > photons;
    iEvent.getByToken(photonToken_,photons);

    Handle<View<reco::Vertex> > primaryVertices;
    iEvent.getByToken( vertexToken_, primaryVertices );

    Handle<View<reco::Conversion> > conversions;
    iEvent.getByToken(conversionToken_,conversions);

    Handle<View<reco::Conversion> > singlelegconversions;
    iEvent.getByToken(singlelegconversionToken_,singlelegconversions);
    
    Handle<double>  rho;
    iEvent.getByToken(rhoTaken_,rho);
    double rho_    = *rho;

    Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken( beamSpotToken_, recoBeamSpotHandle );
    math::XYZPoint beamSpot;
    if( recoBeamSpotHandle.isValid() ) {
        beamSpot = recoBeamSpotHandle->position();
    } else {
        cout << " WARNING NO VALID BEAM SPOT: this should not happen!" << endl;
    }

    Handle<View<reco::GenParticle> > genParticles;
    iEvent.getByToken( genParticleToken_, genParticles );
    const std::vector<edm::Ptr<reco::GenParticle> > genParticlesPtrs = genParticles->ptrs();


    // ********************************************************************************

    initEventStructure();

    unsigned int iphotIsolnAreaValN_ = 7; 
    double photIsolnEAreaVal_[iphotIsolnAreaValN_] = {1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 999999.0};
    double photIsolnEAreaChgHad_[iphotIsolnAreaValN_] = {0.0385, 0.0468, 0.0435, 0.0378, 0.0338, 0.0314, 0.0269};
    double photIsolnEAreaNeuHad_[iphotIsolnAreaValN_] = {0.0636, 0.1103, 0.0759, 0.0236, 0.0151, 0.00007, 0.0132};
    double photIsolnEAreaPhot_[iphotIsolnAreaValN_] = {0.1240, 0.1093, 0.0631, 0.0779, 0.0999, 0.1155, 0.1373};

    // diphoton loop
    vector<int> photon_index;
    for( size_t ipho = 0; ipho < photons->size(); ipho++ ) {

        Ptr<flashgg::Photon> photon = photons->ptrAt(ipho);
        if (photon->pt() < 20) continue;
        if ( !photon->passElectronVeto() ) continue;
        double photPt_ = photon->pt();
        int iphotEA_ = iphotIsolnAreaValN_ - 1;
        for (int iEA = 0; iEA < (int)iphotIsolnAreaValN_; iEA++){
            if (fabs(photon->eta()) < photIsolnEAreaVal_[iEA]) {
                iphotEA_ = iEA;
                break;
            }
        }

        //Cut base photon ID
        if ( photon->isEB() ){
            if ( photon->hadronicOverEm() > 0.105 )  continue;
            if ( photon->full5x5_sigmaIetaIeta() > 0.0103 ) continue;
            double photChargedHadronIso_ = max(photon->chargedHadronIso() - rho_*photIsolnEAreaChgHad_[iphotEA_], 0.);
            if (photChargedHadronIso_ > 2.839 ) continue;
            double photNeutralHadronIso_ = max(photon->neutralHadronIso() - rho_*photIsolnEAreaNeuHad_[iphotEA_], 0.);
            if (photNeutralHadronIso_ - 0.0126*photPt_ - 0.000026*photPt_*photPt_ > 9.188 ) continue;
            double photPhotonIso_ = max(photon->photonIso() - rho_*photIsolnEAreaPhot_[iphotEA_], 0.);
            if (photPhotonIso_ - 0.0035*photPt_ > 2.956 ) continue;
        }
        if ( photon->isEE() ){
            if ( photon->hadronicOverEm() > 0.029 )  continue;
            if ( photon->full5x5_sigmaIetaIeta() > 0.0276 ) continue;
            double photChargedHadronIso_ = max(photon->chargedHadronIso() - rho_*photIsolnEAreaChgHad_[iphotEA_], 0.);
            if (photChargedHadronIso_ > 2.150 ) continue;
            double photNeutralHadronIso_ = max(photon->neutralHadronIso() - rho_*photIsolnEAreaNeuHad_[iphotEA_], 0.);
            if (photNeutralHadronIso_ - 0.0119*photPt_ - 0.000025*photPt_*photPt_ > 10.471 ) continue;
            double photPhotonIso_ = max(photon->photonIso() - rho_*photIsolnEAreaPhot_[iphotEA_], 0.);
            if (photPhotonIso_ - 0.0040*photPt_ > 4.895 ) continue;
        }
        photon_index.push_back(ipho);
    }

    vector<double> zposition = ZmcTruthVertexAndRecoVectex( genParticles->ptrs(), primaryVertices->ptrs() );

    if (photon_index.size() >= 2 && zposition.size() == 2) {
        Ptr<flashgg::Photon> photon1 = photons->ptrAt(photon_index[0]);
        Ptr<flashgg::Photon> photon2 = photons->ptrAt(photon_index[1]);
        resInfo.dipho_mass = (photon1->p4() + photon2->p4()).M();
        resInfo.nphotons = (int)photon_index.size();
        resInfo.ztrue = zposition[0];
        resInfo.zvtx = zposition[1];

        for (int ipho = 0; ipho < (int)photon_index.size(); ipho++) {
            Ptr<flashgg::Photon> photon = photons->ptrAt(photon_index[ipho]);
            resInfo.pho_pt = photon->pt();
            resInfo.photon_index = ipho;

            vector<int> indexConversion = IndexMatchedConversion(photon, conversions->ptrs(), singlelegconversions->ptrs(), true);
            if (indexConversion[1] == -1) continue;
            Ptr<reco:: Conversion> conversion = indexConversion[1] == 2? conversions->ptrAt(indexConversion[0]) : singlelegconversions->ptrAt(indexConversion[0]);
            double vtxZFromConv   = vtxZFromConvOnly(photon, conversion, beamSpot);
            double vtxZFromConvSC = vtxZFromConvSuperCluster(photon, conversion, beamSpot);

            if (fabs(photon->superCluster()->eta()) < 1.4442) {
                double perp = sqrt( conversion->conversionVertex().x() * conversion->conversionVertex().x() + conversion->conversionVertex().y() *
                                                    conversion->conversionVertex().y() );
                if (perp <= 15.0) {
                    resInfo.z1Pix = vtxZFromConv;
                    resInfo.z2Pix = vtxZFromConvSC;
                } else if (perp > 15 && perp <= 60.0 ) {
                    resInfo.z1Tib = vtxZFromConv;
                    resInfo.z2Tib = vtxZFromConvSC;
                } else {
                    resInfo.z1Tob = vtxZFromConv;
                    resInfo.z2Tob = vtxZFromConvSC;
                }
            } else if (fabs(photon->superCluster()->eta()) > 1.566) {
                if( fabs( conversion->conversionVertex().z() ) <= 50.0 ) {
                    resInfo.z1PixFwd = vtxZFromConv;
                    resInfo.z2PixFwd = vtxZFromConvSC;
                } else if (fabs( conversion->conversionVertex().z() ) > 50.0 && fabs( conversion->conversionVertex().z() ) <= 100.0) {
                    resInfo.z1Tid = vtxZFromConv;
                    resInfo.z2Tid = vtxZFromConvSC;
                } else {
                    resInfo.z1Tec = vtxZFromConv;
                    resInfo.z2Tec = vtxZFromConvSC;
                }
            }
            resInfo.ntracks = indexConversion[1];
            resolutionTree->Fill();
        }  // end photon candidate loop
    }
}


void
ConvertedPhotonResoTreeMaker::beginJob()
{
    resolutionTree = fs_->make<TTree>( "resolutionTree", "per-diphoton tree" );
    resolutionTree->Branch( "ntracks"     , &resInfo.ntracks,      "ntracks/I"      );
    resolutionTree->Branch( "nphotons"    , &resInfo.nphotons,     "nphotons/I"     );
    resolutionTree->Branch( "photon_index", &resInfo.photon_index, "photon_index/I" );
    resolutionTree->Branch( "dipho_mass"  , &resInfo.dipho_mass,   "dipho_mass/D"   );
    resolutionTree->Branch( "pho_pt"      , &resInfo.pho_pt,       "pho_pt/D"       );
    resolutionTree->Branch( "ztrue"       , &resInfo.ztrue,        "ztrue/D"        );
    resolutionTree->Branch( "zvtx"        , &resInfo.zvtx,         "zvtx/D"         );

    resolutionTree->Branch( "z1Pix"   , &resInfo.z1Pix,    "z1Pix/D"    );
    resolutionTree->Branch( "z1Tib"   , &resInfo.z1Tib,    "z1Tib/D"    );
    resolutionTree->Branch( "z1Tob"   , &resInfo.z1Tob,    "z1Tob/D"    );
    resolutionTree->Branch( "z1PixFwd", &resInfo.z1PixFwd, "z1PixFwd/D" );
    resolutionTree->Branch( "z1Tid"   , &resInfo.z1Tid,    "z1Tid/D"    );
    resolutionTree->Branch( "z1Tec"   , &resInfo.z1Tec,    "z1Tec/D"    );
    resolutionTree->Branch( "z2Pix"   , &resInfo.z2Pix,    "z2Pix/D"    );
    resolutionTree->Branch( "z2Tib"   , &resInfo.z2Tib,    "z2Tib/D"    );
    resolutionTree->Branch( "z2Tob"   , &resInfo.z2Tob,    "z2Tob/D"    );
    resolutionTree->Branch( "z2PixFwd", &resInfo.z2PixFwd, "z2PixFwd/D" );
    resolutionTree->Branch( "z2Tid"   , &resInfo.z2Tid,    "z2Tid/D"    );
    resolutionTree->Branch( "z2Tec"   , &resInfo.z2Tec,    "z2Tec/D"    );
}

void
ConvertedPhotonResoTreeMaker::endJob()
{
}

void
ConvertedPhotonResoTreeMaker::initEventStructure()
{
    resInfo.ntracks       = -999.;
    resInfo.nphotons      = -999.;
    resInfo.photon_index  = -999.;
    resInfo.dipho_mass    = -999.;
    resInfo.pho_pt        = -999.;
    resInfo.ztrue         = -999.;
    resInfo.zvtx          = -999.;

    resInfo.z1Pix     = -999.;
    resInfo.z1Tib     = -999.;
    resInfo.z1Tob     = -999.;
    resInfo.z1PixFwd  = -999.;
    resInfo.z1Tid     = -999.;
    resInfo.z1Tec     = -999.;
    resInfo.z2Pix     = -999.;
    resInfo.z2Tib     = -999.;
    resInfo.z2Tob     = -999.;
    resInfo.z2PixFwd  = -999.;
    resInfo.z2Tid     = -999.;
    resInfo.z2Tec     = -999.;
}

void
ConvertedPhotonResoTreeMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}

std::vector<double>
ConvertedPhotonResoTreeMaker::ZmcTruthVertexAndRecoVectex( const std::vector<edm::Ptr<reco::GenParticle> > &genParticles ,
        const std::vector<edm::Ptr<reco::Vertex> > &vertices, double dzMatch )
{

    reco::Vertex::Point hardVertex( 0, 0, 0 );

    for( unsigned int genLoop = 0 ; genLoop < genParticles.size(); genLoop++ ) {

        if( fabs( genParticles[genLoop]->pdgId() ) == 25 ) {
            hardVertex.SetCoordinates( genParticles[genLoop]->vx(), genParticles[genLoop]->vy(), genParticles[genLoop]->vz() );
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

    std::vector<double> zposition;
    if( dzMin < dzMatch ) { 
        zposition.push_back(hardVertex.z());
        zposition.push_back(vertices[ivMatch]->z()); 
    }

    return zposition;
}

double 
ConvertedPhotonResoTreeMaker::vtxZFromConvOnly( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco:: Conversion> &conversion,
        const math::XYZPoint &beamSpot ) const
{
    double dz = 0;
    if( conversion->nTracks() == 2 ) {
        double r = sqrt( conversion->refittedPairMomentum().perp2() );
        dz = ( conversion->conversionVertex().z() - beamSpot.z() )
             -
             ( ( conversion->conversionVertex().x() - beamSpot.x() ) * conversion->refittedPair4Momentum().x() + ( conversion->conversionVertex().y() - beamSpot.y() ) *
               conversion->refittedPair4Momentum().y() ) / r * conversion->refittedPair4Momentum().z() / r;
    }
    if( conversion->nTracks() == 1 ) {
        double r = sqrt( conversion->tracksPin()[0].x() * conversion->tracksPin()[0].x() + conversion->tracksPin()[0].y() * conversion->tracksPin()[0].y() );
        dz = ( conversion->conversionVertex().z() - beamSpot.z() )
             -
             ( ( conversion->conversionVertex().x() - beamSpot.x() ) * conversion->tracksPin()[0].x() + ( conversion->conversionVertex().y() - beamSpot.y() ) *
               conversion->tracksPin()[0].y() ) / r * conversion->tracksPin()[0].z() / r;
    }
    return dz + beamSpot.z();
}

double 
ConvertedPhotonResoTreeMaker::vtxZFromConvSuperCluster( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco:: Conversion> &conversion,
        const math::XYZPoint &beamSpot ) const
{
    // get the z from conversion plus SuperCluster
    double deltaX1 =  pho->caloPosition().x() - conversion->conversionVertex().x();
    double deltaY1 =  pho->caloPosition().y() - conversion->conversionVertex().y();
    double deltaZ1 =  pho->caloPosition().z() - conversion->conversionVertex().z();
    double R1 = sqrt( deltaX1 * deltaX1 + deltaY1 * deltaY1 );
    double tantheta = R1 / deltaZ1;
                                                                                                                                            
    double deltaX2 = conversion->conversionVertex().x() - beamSpot.x();
    double deltaY2 = conversion->conversionVertex().y() - beamSpot.y();
    double R2 = sqrt( deltaX2 * deltaX2 + deltaY2 * deltaY2 );
    double deltaZ2 = R2 / tantheta;
    double higgsZ =  pho->caloPosition().z() - deltaZ1 - deltaZ2;
    return higgsZ;
}

std::vector<int> 
ConvertedPhotonResoTreeMaker::IndexMatchedConversion( const edm::Ptr<flashgg::Photon> &g,
        const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,  const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg, 
        bool useSingleLeg ) const
{
    double mindR = 999;
    int nConvLegs = 0;
    bool doOneLeg = true;
    bool pureGeomConvMatching = true;
                                                                                                                                                        
    std::vector<int> result;
                                                                                                                                                        
    if(!pureGeomConvMatching) assert( g->hasConversionTracks() );
    int selected_conversion_index = -1;
                                                                                                                                                        
    if( (g->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching){
        
        for( unsigned int i = 0; i < conversionsVector.size(); i++ ) {
            edm::Ptr<reco::Conversion> conv = conversionsVector[i];
            if( conv->nTracks() == 2 ) {
                if( !conv->isConverted() ) { continue; }
                if( conv->refittedPair4Momentum().pt() < 10. ) { continue; }
                if( TMath::Prob( conv->conversionVertex().chi2(), conv->conversionVertex().ndof() ) < 1e-6 ) { continue; }
                
                TVector3 VtxtoSC;
                VtxtoSC.SetXYZ( g->superCluster()->position().x() - conv->conversionVertex().x(),
                                g->superCluster()->position().y() - conv->conversionVertex().y(),
                                g->superCluster()->position().z() - conv->conversionVertex().z() );
                TVector3 RefPairMo;
                RefPairMo.SetXYZ( conv->refittedPairMomentum().x(), conv->refittedPairMomentum().y(), conv->refittedPairMomentum().z() );
                double dR = 0;
                dR = VtxtoSC.DeltaR( RefPairMo );
                if( dR < mindR ) {
                    mindR = dR;
                    selected_conversion_index = i;
                }
            }
        }
        if( mindR < 0.1 ) {
            result.push_back( selected_conversion_index );
            nConvLegs = 2;
            result.push_back( nConvLegs );
            doOneLeg = false;
        }
        if( doOneLeg && useSingleLeg ) {
            mindR = 999;
            for( unsigned int j = 0; j < conversionsVectorSingleLeg.size(); j++ ) {
                edm::Ptr<reco::Conversion> conv = conversionsVectorSingleLeg[j];
                if( conv->nTracks() == 1 ) {
                    TVector3 VtxtoSC;
                    VtxtoSC.SetXYZ( g->superCluster()->position().x() - conv->conversionVertex().x(),
                                    g->superCluster()->position().y() - conv->conversionVertex().y(),
                                    g->superCluster()->position().z() - conv->conversionVertex().z() );
                    TVector3 RefPairMo;
                    float oneLegTrack_X = conv->tracksPin()[0].x();
                    float oneLegTrack_Y = conv->tracksPin()[0].y();
                    float oneLegTrack_Z = conv->tracksPin()[0].z();
                    
                    RefPairMo.SetXYZ( oneLegTrack_X, oneLegTrack_Y, oneLegTrack_Z );
                    double dR = 0;
                    dR = VtxtoSC.DeltaR( RefPairMo );
                    if( dR < mindR ) {
                        mindR = dR;
                        selected_conversion_index = j;
                    }                        
                }
            }
            if( mindR < 0.1 ) {
                result.push_back( selected_conversion_index );
                nConvLegs = 1;
                result.push_back( nConvLegs );
            }
        }
    }
    
    if( mindR < 0.1 )
        {return result;}
    else {
        result.push_back( -1 );
        result.push_back( -1 );
        return result;
    }
}

DEFINE_FWK_MODULE( ConvertedPhotonResoTreeMaker );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

