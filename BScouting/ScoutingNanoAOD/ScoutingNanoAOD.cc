// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>

// ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include <DataFormats/TrackReco/interface/TrackBase.h>

#include "DataFormats/Math/interface/libminifloat.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/TrackReco/interface/fillCovariance.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"

// Root include files
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// User include files

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

using namespace std;


class ScoutingNanoAOD : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
    explicit ScoutingNanoAOD(const edm::ParameterSet&);
    ~ScoutingNanoAOD();
		
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    virtual void clearVars();

    const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> >      muonsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> >  	electronsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> >  	photonsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  	pfcandsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet> >  	pfjetsToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingTrack> >  	tracksToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >  	pvToken;
    const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >  	dvToken;
    const edm::EDGetTokenT<double>  	rhoToken;
    const edm::EDGetTokenT<double>  	pfMetToken;
    const edm::EDGetTokenT<double>  	pfMetPhiToken;

    //const edm::EDGetTokenT<GenEventInfoProduct>             genEvtInfoToken;	
    bool doL1;       
    const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken;
	
    edm::InputTag                algInputTag_;       
    edm::InputTag                extInputTag_;       
    edm::EDGetToken              algToken_;
    //l1t::L1TGlobalUtil          *l1GtUtils_;
    std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
    std::vector<std::string>     l1Seeds_;
    std::vector<bool>            l1Result_;

    //Photon
    const static int 	max_pho = 1000;
    UInt_t n_pho;
    vector<Float_t> Photon_pt_;
    vector<Float_t> Photon_eta_;
    vector<Float_t> Photon_phi_;
    vector<Float_t> Photon_m_;
    vector<Float_t> Photon_sigmaIetaIeta_;
    vector<Float_t> Photon_hOverE_;
    vector<Float_t> Photon_ecalIso_;
    vector<Float_t> Photon_hcalIso_;
    vector<Float_t> Photon_trkIso_;
    vector<Float_t> Photon_r9_;
    vector<Float_t> Photon_sMin_;
    vector<Float_t> Photon_sMaj_;
    vector<UInt_t> Photon_seedId_;
    vector<vector<Float_t> > Photon_energyMatrix_;
    //vector<vector<UInt_t> > Photon_detIds_;
    vector<vector<Float_t> > Photon_timingMatrix_;
    //vector<bool> Photon_rechitZeroSuppression_;

    //Electron
    const static int 	max_ele = 1000;
    UInt_t n_ele;
    vector<Float_t> Electron_pt_;
    vector<Float_t> Electron_eta_;
    vector<Float_t> Electron_phi_;
    vector<Float_t> Electron_m_;
    vector<Float_t> Electron_d0_;
    vector<Float_t> Electron_dz_;
    vector<Float_t> Electron_dEtaIn_;
    vector<Float_t> Electron_dPhiIn_;
    vector<Float_t> Electron_sigmaIetaIeta_;
    vector<Float_t> Electron_hOverE_;
    vector<Float_t> Electron_ooEMOop_;
    vector<Int_t> Electron_missingHits_;
    vector<Int_t> Electron_charge_;
    vector<Float_t> Electron_ecalIso_;
    vector<Float_t> Electron_hcalIso_;
    vector<Float_t> Electron_trackIso_;
    vector<Float_t> Electron_r9_;
    vector<Float_t> Electron_sMin_;
    vector<Float_t> Electron_sMaj_;
    vector<UInt_t> Electron_seedId_;
    vector<vector<Float_t> > Electron_energyMatrix_;
    //vector<vector<UInt_t> > Electron_detIds_;
    vector<vector<Float_t> > Electron_timingMatrix_;
    //vector<bool> Electron_rechitZeroSuppression_;

    //Muon
    const static int 	max_mu = 1000;
    UInt_t n_mu;
    vector<Float_t> Muon_pt_;
    vector<Float_t> Muon_eta_;
    vector<Float_t> Muon_phi_;
    vector<Float_t> Muon_m_;
    vector<UInt_t> Muon_type_;
    vector<UInt_t> Muon_charge_;
    vector<Float_t> Muon_normalizedChi2_;
    vector<Float_t> Muon_ecalIso_;
    vector<Float_t> Muon_hcalIso_;
    vector<Float_t> Muon_trackIso_;
    vector<UInt_t> Muon_nValidStandAloneMuonHits_;
    vector<UInt_t> Muon_nStandAloneMuonMatchedStations_;
    vector<UInt_t> Muon_nValidRecoMuonHits_;
    vector<UInt_t> Muon_nRecoMuonChambers_;
    vector<UInt_t> Muon_nRecoMuonChambersCSCorDT_;
    vector<UInt_t> Muon_nRecoMuonMatches_;
    vector<UInt_t> Muon_nRecoMuonMatchedStations_;
    vector<UInt_t> Muon_nRecoMuonExpectedMatchedStations_;
    vector<UInt_t> Muon_recoMuonStationMask_;
    vector<UInt_t> Muon_nRecoMuonMatchedRPCLayers_;
    vector<UInt_t> Muon_recoMuonRPClayerMask_;
    vector<UInt_t> Muon_nValidPixelHits_;
    vector<UInt_t> Muon_nValidStripHits_;
    vector<UInt_t> Muon_nPixelLayersWithMeasurement_;
    vector<UInt_t> Muon_nTrackerLayersWithMeasurement_;
    vector<Float_t> Muon_trk_chi2_;
    vector<Float_t> Muon_trk_ndof_;
    vector<Float_t> Muon_trk_dxy_;
    vector<Float_t> Muon_trk_dz_;
    vector<Float_t> Muon_trk_qoverp_;
    vector<Float_t> Muon_trk_lambda_;
    vector<Float_t> Muon_trk_pt_;
    vector<Float_t> Muon_trk_phi_;
    vector<Float_t> Muon_trk_eta_;
    vector<Float_t> Muon_trk_dxyError_;
    vector<Float_t> Muon_trk_dzError_;
    vector<Float_t> Muon_trk_qoverpError_;
    vector<Float_t> Muon_trk_lambdaError_;
    vector<Float_t> Muon_trk_phiError_;
    vector<Float_t> Muon_trk_dsz_;
    vector<Float_t> Muon_trk_dszError_;
    vector<Float_t> Muon_trk_qoverp_lambda_cov_;
    vector<Float_t> Muon_trk_qoverp_phi_cov_;
    vector<Float_t> Muon_trk_qoverp_dxy_cov_;
    vector<Float_t> Muon_trk_qoverp_dsz_cov_;
    vector<Float_t> Muon_trk_lambda_phi_cov_;
    vector<Float_t> Muon_trk_lambda_dxy_cov_;
    vector<Float_t> Muon_trk_lambda_dsz_cov_;
    vector<Float_t> Muon_trk_phi_dxy_cov_;
    vector<Float_t> Muon_trk_phi_dsz_cov_;
    vector<Float_t> Muon_trk_dxy_dsz_cov_;
    vector<Float_t> Muon_trk_vx_;
    vector<Float_t> Muon_trk_vy_;
    vector<Float_t> Muon_trk_vz_;
    //vector<reco::HitPattern> trk_hitPattern_;
    vector<vector<Int_t> > Muon_vtxIndx_;


    //PFJets
    const static int 	max_jet = 1000;
    UInt_t n_jet;
    vector<Float_t> Jet_pt_;
    vector<Float_t> Jet_eta_;
    vector<Float_t> Jet_phi_;
    vector<Float_t> Jet_m_;
    vector<Float_t> Jet_area_;
    vector<Float_t> Jet_chargedHadronEnergy_;
    vector<Float_t> Jet_neutralHadronEnergy_;
    vector<Float_t> Jet_photonEnergy_;
    vector<Float_t> Jet_electronEnergy_;
    vector<Float_t> Jet_muonEnergy_;
    vector<Float_t> Jet_HFHadronEnergy_;
    vector<Float_t> Jet_HFEMEnergy_;
    vector<Int_t> Jet_chargedHadronMultiplicity_;
    vector<Int_t> Jet_neutralHadronMultiplicity_;
    vector<Int_t> Jet_photonMultiplicity_;
    vector<Int_t> Jet_electronMultiplicity_;
    vector<Int_t> Jet_muonMultiplicity_;
    vector<Int_t> Jet_HFHadronMultiplicity_;
    vector<Int_t> Jet_HFEMMultiplicity_;
    vector<Float_t> Jet_HOEnergy_;
    vector<Float_t> Jet_csv_;
    vector<Float_t> Jet_mvaDiscriminator_;
    vector<vector<Int_t> > Jet_constituents_;

    //PFCand
    const static int 	max_pfcand = 10000;
    UInt_t n_pfcand;
    vector<Float_t> PFCand_pt_;
    vector<Float_t> PFCand_eta_;
    vector<Float_t> PFCand_phi_;
    vector<Int_t> PFCand_pdgId_;
    vector<Int_t> PFCand_vertex_;

    //Track
    const static int max_track = 10000;
    UInt_t n_track;
    vector<Float_t> Track_pt_;
    vector<Float_t> Track_eta_;
    vector<Float_t> Track_phi_;
    vector<Float_t> Track_chi2_;
    vector<Float_t> Track_ndof_;
    vector<Int_t> Track_charge_;
    vector<Float_t> Track_dxy_;
    vector<Float_t> Track_dz_;
    vector<Int_t> Track_nValidPixelHits_;
    vector<Int_t> Track_nTrackerLayersWithMeasurement_;
    vector<Int_t> Track_nValidStripHits_;
    vector<Float_t> Track_qoverp_;
    vector<Float_t> Track_lambda_;
    vector<Float_t> Track_dxy_Error_;
    vector<Float_t> Track_dz_Error_;
    vector<Float_t> Track_qoverp_Error_;
    vector<Float_t> Track_lambda_Error_;
    vector<Float_t> Track_phi_Error_;
    vector<Float_t> Track_dsz_;
    vector<Float_t> Track_dsz_Error_;
    vector<Float_t> Track_qoverp_lambda_cov_;
    vector<Float_t> Track_qoverp_phi_cov_;
    vector<Float_t> Track_qoverp_dxy_cov_;
    vector<Float_t> Track_qoverp_dsz_cov_;
    vector<Float_t> Track_lambda_phi_cov_;
    vector<Float_t> Track_lambda_dxy_cov_;
    vector<Float_t> Track_lambda_dsz_cov_;
    vector<Float_t> Track_phi_dxy_cov_;
    vector<Float_t> Track_phi_dsz_cov_;
    vector<Float_t> Track_dxy_dsz_cov_;
    vector<Int_t> Track_vtxInd_;
    vector<Float_t> Track_vx_;
    vector<Float_t> Track_vy_;
    vector<Float_t> Track_vz_;

    //PrimaryVertex
    UInt_t n_primaryvtx;
    vector<Float_t> PrimaryVtx_x_;
    vector<Float_t> PrimaryVtx_y_;
    vector<Float_t> PrimaryVtx_z_;
    vector<Float_t> PrimaryVtx_zError_;
    vector<Float_t> PrimaryVtx_xError_;
    vector<Float_t> PrimaryVtx_yError_;
    vector<Int_t> PrimaryVtx_tracksSize_;
    vector<Float_t> PrimaryVtx_chi2_;
    vector<Int_t> PrimaryVtx_ndof_;
    vector<Bool_t> PrimaryVtx_isValidVtx_;

    //DisplacedVertex
    UInt_t n_displacedvtx;
    vector<Float_t> DisplacedVtx_x_;
    vector<Float_t> DisplacedVtx_y_;
    vector<Float_t> DisplacedVtx_z_;
    vector<Float_t> DisplacedVtx_zError_;
    vector<Float_t> DisplacedVtx_xError_;
    vector<Float_t> DisplacedVtx_yError_;
    vector<Int_t> DisplacedVtx_tracksSize_;
    vector<Float_t> DisplacedVtx_chi2_;
    vector<Int_t> DisplacedVtx_ndof_;
    vector<Bool_t> DisplacedVtx_isValidVtx_;

    //rho
    double rho;
    
    //PF Met
    double pfMet;
    double pfMetPhi;

    // TTree carrying the event weight information
    TTree* tree;

    //Run and lumisection
    int run;
    int lumSec;

};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig): 
    muonsToken               (consumes<std::vector<Run3ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))), 
    electronsToken           (consumes<std::vector<Run3ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))), 
    photonsToken           (consumes<std::vector<Run3ScoutingPhoton> >         (iConfig.getParameter<edm::InputTag>("photons"))), 
    pfcandsToken             (consumes<std::vector<Run3ScoutingParticle> >         (iConfig.getParameter<edm::InputTag>("pfcands"))), 
    pfjetsToken              (consumes<std::vector<Run3ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets"))),
    tracksToken              (consumes<std::vector<Run3ScoutingTrack> >            (iConfig.getParameter<edm::InputTag>("tracks"))), 
    pvToken              (consumes<std::vector<Run3ScoutingVertex> >            (iConfig.getParameter<edm::InputTag>("primaryVertices"))), 
    dvToken              (consumes<std::vector<Run3ScoutingVertex> >            (iConfig.getParameter<edm::InputTag>("displacedVertices"))), 
    rhoToken             (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
    pfMetToken             (consumes<double>(iConfig.getParameter<edm::InputTag>("pfMet"))),
    pfMetPhiToken             (consumes<double>(iConfig.getParameter<edm::InputTag>("pfMetPhi"))),
    doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false),
    particleTableToken       (esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>())
{
    usesResource("TFileService");
    if (doL1) {
        algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
        extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
        algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
        l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
        l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
    }
    else {
        l1Seeds_ = std::vector<std::string>();
        l1GtUtils_ = 0;
    }

    // Access the TFileService
    edm::Service<TFileService> fs;

    // Create the TTree
    tree = fs->make<TTree>("Events", "Events");

    // Event weights
    tree->Branch("lumSec"		, &lumSec			 , "lumSec/i" );
    tree->Branch("run"			, &run				 , "run/i" );
    //tree->Branch("nvtx"			, &nvtx				 , "nvtx/i" );
    
    // Triggers
    tree->Branch("l1Result_"		, "std::vector<bool>"             ,&l1Result_	, 32000, 0);		
   
    // Pileup info
    //tree->Branch("nvtx"                 , &nvtx                          , "nvtx/i"       );

    //Electrons
    tree->Branch("n_ele", &n_ele, "n_ele/i");
    tree->Branch("Electron_pt", &Electron_pt_ );
    tree->Branch("Electron_eta", &Electron_eta_ );
    tree->Branch("Electron_phi", &Electron_phi_ );
    tree->Branch("Electron_m", &Electron_m_ );
    tree->Branch("Electron_d0", &Electron_d0_ );
    tree->Branch("Electron_dz", &Electron_dz_ );
    tree->Branch("Electron_dEtaIn", &Electron_dEtaIn_ );
    tree->Branch("Electron_dPhiIn", &Electron_dPhiIn_ );
    tree->Branch("Electron_sigmaIetaIeta", &Electron_sigmaIetaIeta_ );
    tree->Branch("Electron_hOverE", &Electron_hOverE_ );
    tree->Branch("Electron_ooEMOop", &Electron_ooEMOop_ );
    tree->Branch("Electron_missingHits", &Electron_missingHits_ );
    tree->Branch("Electron_charge", &Electron_charge_ );
    tree->Branch("Electron_ecalIso", &Electron_ecalIso_ );
    tree->Branch("Electron_hcalIso", &Electron_hcalIso_ );
    tree->Branch("Electron_trackIso", &Electron_trackIso_ );
    tree->Branch("Electron_r9", &Electron_r9_ );
    tree->Branch("Electron_sMin", &Electron_sMin_ );
    tree->Branch("Electron_sMaj", &Electron_sMaj_ );
    tree->Branch("Electron_seedId", &Electron_seedId_ );
    //tree->Branch("Electron_detIds", &Electron_detIds_ );
    tree->Branch("Electron_energyMatrix", &Electron_energyMatrix_ );
    tree->Branch("Electron_timingMatrix", &Electron_timingMatrix_ );
    //tree->Branch("Electron_rechitZeroSuppression", &Electron_rechitZeroSuppression_ );

    //Muons
    tree->Branch("n_mu"            	   ,&n_mu 			, "n_mu/i"		);
    tree->Branch("Muon_pt", &Muon_pt_ );
    tree->Branch("Muon_eta", &Muon_eta_ );
    tree->Branch("Muon_phi", &Muon_phi_ );
    tree->Branch("Muon_m", &Muon_m_ );
    tree->Branch("Muon_type", &Muon_type_ );
    tree->Branch("Muon_charge", &Muon_charge_ );
    tree->Branch("Muon_normalizedChi2", &Muon_normalizedChi2_ );
    tree->Branch("Muon_ecalIso", &Muon_ecalIso_ );
    tree->Branch("Muon_hcalIso", &Muon_hcalIso_ );
    tree->Branch("Muon_trackIso", &Muon_trackIso_ );
    tree->Branch("Muon_nValidStandAloneMuonHits", &Muon_nValidStandAloneMuonHits_ );
    tree->Branch("Muon_nStandAloneMuonMatchedStations", &Muon_nStandAloneMuonMatchedStations_ );
    tree->Branch("Muon_nValidRecoMuonHits", &Muon_nValidRecoMuonHits_ );
    tree->Branch("Muon_nRecoMuonChambers", &Muon_nRecoMuonChambers_ );
    tree->Branch("Muon_nRecoMuonChambersCSCorDT", &Muon_nRecoMuonChambersCSCorDT_ );
    tree->Branch("Muon_nRecoMuonMatches", &Muon_nRecoMuonMatches_ );
    tree->Branch("Muon_nRecoMuonMatchedStations", &Muon_nRecoMuonMatchedStations_ );
    tree->Branch("Muon_nRecoMuonExpectedMatchedStations", &Muon_nRecoMuonExpectedMatchedStations_ );
    tree->Branch("Muon_recoMuonStationMask", &Muon_recoMuonStationMask_ );
    tree->Branch("Muon_nRecoMuonMatchedRPCLayers", &Muon_nRecoMuonMatchedRPCLayers_ );
    tree->Branch("Muon_recoMuonRPClayerMask", &Muon_recoMuonRPClayerMask_ );
    tree->Branch("Muon_nValidPixelHits", &Muon_nValidPixelHits_ );
    tree->Branch("Muon_nValidStripHits", &Muon_nValidStripHits_ );
    tree->Branch("Muon_nPixelLayersWithMeasurement", &Muon_nPixelLayersWithMeasurement_ );
    tree->Branch("Muon_nTrackerLayersWithMeasurement", &Muon_nTrackerLayersWithMeasurement_ );
    tree->Branch("Muon_trk_chi2", &Muon_trk_chi2_ );
    tree->Branch("Muon_trk_ndof", &Muon_trk_ndof_ );
    tree->Branch("Muon_trk_dxy", &Muon_trk_dxy_ );
    tree->Branch("Muon_trk_dz", &Muon_trk_dz_ );
    tree->Branch("Muon_trk_qoverp", &Muon_trk_qoverp_ );
    tree->Branch("Muon_trk_lambda", &Muon_trk_lambda_ );
    tree->Branch("Muon_trk_pt", &Muon_trk_pt_ );
    tree->Branch("Muon_trk_phi", &Muon_trk_phi_ );
    tree->Branch("Muon_trk_eta", &Muon_trk_eta_ );
    tree->Branch("Muon_trk_dxyError", &Muon_trk_dxyError_ );
    tree->Branch("Muon_trk_dzError", &Muon_trk_dzError_ );
    tree->Branch("Muon_trk_qoverpError", &Muon_trk_qoverpError_ );
    tree->Branch("Muon_trk_lambdaError", &Muon_trk_lambdaError_ );
    tree->Branch("Muon_trk_phiError", &Muon_trk_phiError_ );
    tree->Branch("Muon_trk_dsz", &Muon_trk_dsz_ );
    tree->Branch("Muon_trk_dszError", &Muon_trk_dszError_ );
    tree->Branch("Muon_trk_qoverp_lambda_cov", &Muon_trk_qoverp_lambda_cov_ );
    tree->Branch("Muon_trk_qoverp_phi_cov", &Muon_trk_qoverp_phi_cov_ );
    tree->Branch("Muon_trk_qoverp_dxy_cov", &Muon_trk_qoverp_dxy_cov_ );
    tree->Branch("Muon_trk_qoverp_dsz_cov", &Muon_trk_qoverp_dsz_cov_ );
    tree->Branch("Muon_trk_lambda_phi_cov", &Muon_trk_lambda_phi_cov_ );
    tree->Branch("Muon_trk_lambda_dxy_cov", &Muon_trk_lambda_dxy_cov_ );
    tree->Branch("Muon_trk_lambda_dsz_cov", &Muon_trk_lambda_dsz_cov_ );
    tree->Branch("Muon_trk_phi_dxy_cov", &Muon_trk_phi_dxy_cov_ );
    tree->Branch("Muon_trk_phi_dsz_cov", &Muon_trk_phi_dsz_cov_ );
    tree->Branch("Muon_trk_dxy_dsz_cov", &Muon_trk_dxy_dsz_cov_ );
    tree->Branch("Muon_trk_vx", &Muon_trk_vx_ );
    tree->Branch("Muon_trk_vy", &Muon_trk_vy_ );
    tree->Branch("Muon_trk_vz", &Muon_trk_vz_ );
    tree->Branch("Muon_vtxIndx", &Muon_vtxIndx_ );

    //Photons
    tree->Branch("n_pho"            	   ,&n_pho 			, "n_pho/i"		);
    tree->Branch("Photon_pt", &Photon_pt_ );
    tree->Branch("Photon_eta", &Photon_eta_ );
    tree->Branch("Photon_phi", &Photon_phi_ );
    tree->Branch("Photon_m", &Photon_m_ );
    tree->Branch("Photon_sigmaIetaIeta", &Photon_sigmaIetaIeta_ );
    tree->Branch("Photon_hOverE", &Photon_hOverE_ );
    tree->Branch("Photon_ecalIso", &Photon_ecalIso_ );
    tree->Branch("Photon_hcalIso", &Photon_hcalIso_ );
    tree->Branch("Photon_trkIso", &Photon_trkIso_ );
    tree->Branch("Photon_r9", &Photon_r9_ );
    tree->Branch("Photon_sMin", &Photon_sMin_ );
    tree->Branch("Photon_sMaj", &Photon_sMaj_ );
    tree->Branch("Photon_seedId", &Photon_seedId_ );
    //tree->Branch("Photon_detIds", &Photon_detIds_ );
    tree->Branch("Photon_energyMatrix", &Photon_energyMatrix_ );
    tree->Branch("Photon_timingMatrix", &Photon_timingMatrix_ );
    //tree->Branch("Photon_rechitZeroSuppression", &Photon_rechitZeroSuppression_ );

    //PFJets
    tree->Branch("n_jet"            	   	,&n_jet 			, "n_jet/i"		);
    tree->Branch("Jet_pt", &Jet_pt_ );
    tree->Branch("Jet_eta", &Jet_eta_ );
    tree->Branch("Jet_phi", &Jet_phi_ );
    tree->Branch("Jet_m", &Jet_m_ );
    tree->Branch("Jet_area", &Jet_area_ );
    tree->Branch("Jet_chargedHadronEnergy", &Jet_chargedHadronEnergy_ );
    tree->Branch("Jet_neutralHadronEnergy", &Jet_neutralHadronEnergy_ );
    tree->Branch("Jet_photonEnergy", &Jet_photonEnergy_ );
    tree->Branch("Jet_electronEnergy", &Jet_electronEnergy_ );
    tree->Branch("Jet_muonEnergy", &Jet_muonEnergy_ );
    tree->Branch("Jet_HFHadronEnergy", &Jet_HFHadronEnergy_ );
    tree->Branch("Jet_HFEMEnergy", &Jet_HFEMEnergy_ );
    tree->Branch("Jet_chargedHadronMultiplicity", &Jet_chargedHadronMultiplicity_ );
    tree->Branch("Jet_neutralHadronMultiplicity", &Jet_neutralHadronMultiplicity_ );
    tree->Branch("Jet_photonMultiplicity", &Jet_photonMultiplicity_ );
    tree->Branch("Jet_electronMultiplicity", &Jet_electronMultiplicity_ );
    tree->Branch("Jet_muonMultiplicity", &Jet_muonMultiplicity_ );
    tree->Branch("Jet_HFHadronMultiplicity", &Jet_HFHadronMultiplicity_ );
    tree->Branch("Jet_HFEMMultiplicity", &Jet_HFEMMultiplicity_ );
    tree->Branch("Jet_HOEnergy", &Jet_HOEnergy_ );
    tree->Branch("Jet_csv", &Jet_csv_ );
    tree->Branch("Jet_mvaDiscriminator", &Jet_mvaDiscriminator_ );
    tree->Branch("Jet_constituents", &Jet_constituents_ );

    //PFCands
    tree->Branch("n_pfcand"            	   ,&n_pfcand 		,"n_pfcand/i"		);	
    tree->Branch("PFCand_pt", &PFCand_pt_ );
    tree->Branch("PFCand_eta", &PFCand_eta_ );
    tree->Branch("PFCand_phi", &PFCand_phi_ );
    tree->Branch("PFCand_pdgId", &PFCand_pdgId_ );
    tree->Branch("PFCand_vertex", &PFCand_vertex_ );

    //Tracks
    //tree->Branch("n_track"            	   	,&n_track 			, "n_track/i"		);
    //tree->Branch("Track_pt", &Track_pt_ );
    //tree->Branch("Track_eta", &Track_eta_ );
    //tree->Branch("Track_phi", &Track_phi_ );
    //tree->Branch("Track_chi2", &Track_chi2_ );
    //tree->Branch("Track_ndof", &Track_ndof_ );
    //tree->Branch("Track_charge", &Track_charge_ );
    //tree->Branch("Track_dxy", &Track_dxy_ );
    //tree->Branch("Track_dz", &Track_dz_ );
    //tree->Branch("Track_nValidPixelHits", &Track_nValidPixelHits_ );
    //tree->Branch("Track_nTrackerLayersWithMeasurement", &Track_nTrackerLayersWithMeasurement_ );
    //tree->Branch("Track_nValidStripHits", &Track_nValidStripHits_ );
    //tree->Branch("Track_qoverp", &Track_qoverp_ );
    //tree->Branch("Track_lambda", &Track_lambda_ );
    //tree->Branch("Track_dxy_Error", &Track_dxy_Error_ );
    //tree->Branch("Track_dz_Error", &Track_dz_Error_ );
    //tree->Branch("Track_qoverp_Error", &Track_qoverp_Error_ );
    //tree->Branch("Track_lambda_Error", &Track_lambda_Error_ );
    //tree->Branch("Track_phi_Error", &Track_phi_Error_ );
    //tree->Branch("Track_dsz", &Track_dsz_ );
    //tree->Branch("Track_dsz_Error", &Track_dsz_Error_ );
    //tree->Branch("Track_qoverp_lambda_cov", &Track_qoverp_lambda_cov_ );
    //tree->Branch("Track_qoverp_phi_cov", &Track_qoverp_phi_cov_ );
    //tree->Branch("Track_qoverp_dxy_cov", &Track_qoverp_dxy_cov_ );
    //tree->Branch("Track_qoverp_dsz_cov", &Track_qoverp_dsz_cov_ );
    //tree->Branch("Track_lambda_phi_cov", &Track_lambda_phi_cov_ );
    //tree->Branch("Track_lambda_dxy_cov", &Track_lambda_dxy_cov_ );
    //tree->Branch("Track_lambda_dsz_cov", &Track_lambda_dsz_cov_ );
    //tree->Branch("Track_phi_dxy_cov", &Track_phi_dxy_cov_ );
    //tree->Branch("Track_phi_dsz_cov", &Track_phi_dsz_cov_ );
    //tree->Branch("Track_dxy_dsz_cov", &Track_dxy_dsz_cov_ );
    //tree->Branch("Track_vtxInd", &Track_vtxInd_ );
    //tree->Branch("Track_vx", &Track_vx_ );
    //tree->Branch("Track_vy", &Track_vy_ );
    //tree->Branch("Track_vz", &Track_vz_ );

    //Primary Vertex
    tree->Branch("PrimaryVtx_x", &PrimaryVtx_x_ );
    tree->Branch("PrimaryVtx_y", &PrimaryVtx_y_ );
    tree->Branch("PrimaryVtx_z", &PrimaryVtx_z_ );
    tree->Branch("PrimaryVtx_zError", &PrimaryVtx_zError_ );
    tree->Branch("PrimaryVtx_xError", &PrimaryVtx_xError_ );
    tree->Branch("PrimaryVtx_yError", &PrimaryVtx_yError_ );
    tree->Branch("PrimaryVtx_tracksSize", &PrimaryVtx_tracksSize_ );
    tree->Branch("PrimaryVtx_chi2", &PrimaryVtx_chi2_ );
    tree->Branch("PrimaryVtx_ndof", &PrimaryVtx_ndof_ );
    tree->Branch("PrimaryVtx_isValidVtx", &PrimaryVtx_isValidVtx_ );

    //Displaced Vertex
    tree->Branch("DisplacedVtx_x", &DisplacedVtx_x_ );
    tree->Branch("DisplacedVtx_y", &DisplacedVtx_y_ );
    tree->Branch("DisplacedVtx_z", &DisplacedVtx_z_ );
    tree->Branch("DisplacedVtx_zError", &DisplacedVtx_zError_ );
    tree->Branch("DisplacedVtx_xError", &DisplacedVtx_xError_ );
    tree->Branch("DisplacedVtx_yError", &DisplacedVtx_yError_ );
    tree->Branch("DisplacedVtx_tracksSize", &DisplacedVtx_tracksSize_ );
    tree->Branch("DisplacedVtx_chi2", &DisplacedVtx_chi2_ );
    tree->Branch("DisplacedVtx_ndof", &DisplacedVtx_ndof_ );
    tree->Branch("DisplacedVtx_isValidVtx", &DisplacedVtx_isValidVtx_ );

    //PF Met
    tree->Branch("PFMet_Pt", &pfMet );
    tree->Branch("PFMet_Phi", &pfMetPhi );

    //
    tree->Branch("rho", &rho );

}


ScoutingNanoAOD::~ScoutingNanoAOD() {
}

void ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;
    using namespace reco;
    using namespace fastjet;
    using namespace fastjet::contrib;
    
    Handle<vector<Run3ScoutingElectron> > electronsH;
    iEvent.getByToken(electronsToken, electronsH);

    Handle<vector<Run3ScoutingMuon> > muonsH;
    iEvent.getByToken(muonsToken, muonsH);

    Handle<vector<Run3ScoutingPhoton> > photonsH;
    iEvent.getByToken(photonsToken, photonsH);

    Handle<vector<Run3ScoutingPFJet> > pfjetsH;
    iEvent.getByToken(pfjetsToken, pfjetsH);
    
    Handle<vector<Run3ScoutingParticle> > pfcandsH;
    iEvent.getByToken(pfcandsToken, pfcandsH);

    Handle<vector<Run3ScoutingTrack> > tracksH;
    iEvent.getByToken(tracksToken, tracksH);

    Handle<vector<Run3ScoutingVertex> > pvH;
    iEvent.getByToken(pvToken, pvH);

    Handle<vector<Run3ScoutingVertex> > dvH;
    iEvent.getByToken(dvToken, dvH);

    Handle<double> rhoH;
    iEvent.getByToken(rhoToken, rhoH);
    rho = *rhoH;

    Handle<double> pfMetH;
    iEvent.getByToken(pfMetToken, pfMetH);
    pfMet = *pfMetH;

    Handle<double> pfMetPhiH;
    iEvent.getByToken(pfMetPhiToken, pfMetPhiH);
    pfMetPhi = *pfMetPhiH;
    
    run = iEvent.eventAuxiliary().run();
    lumSec = iEvent.eventAuxiliary().luminosityBlock();

    auto pdt = iSetup.getHandle(particleTableToken);
    const HepPDT::ParticleDataTable* pdTable = pdt.product();

    if (doL1) {
        l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
        //for( int r = 0; r<280; r++){
        //    string name ("empty");
        //    bool algoName_ = false;
        //    algoName_ = l1GtUtils_->getAlgNameFromBit(r,name);
        //    cout << "getAlgNameFromBit = " << algoName_  << endl;
        //    cout << "L1 bit number = " << r << " ; L1 bit name = " << name << endl;
        //}
        for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
            bool l1htbit = 0;	
			
            l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
            //cout<<string(l1Seeds_[iseed])<<"  "<<l1htbit<<endl;
            l1Result_.push_back( l1htbit );
        }
    }  

    n_ele = 0;
    for (auto iter = electronsH->begin(); iter != electronsH->end(); ++iter) 
    {

        Electron_pt_.push_back(iter->pt());
        Electron_eta_.push_back(iter->eta());
        Electron_phi_.push_back(iter->phi());
        Electron_m_.push_back(iter->m());
        Electron_d0_.push_back(iter->d0());
        Electron_dz_.push_back(iter->dz());
        Electron_dEtaIn_.push_back(iter->dEtaIn());
        Electron_dPhiIn_.push_back(iter->dPhiIn());
        Electron_sigmaIetaIeta_.push_back(iter->sigmaIetaIeta());
        Electron_hOverE_.push_back(iter->hOverE());
        Electron_ooEMOop_.push_back(iter->ooEMOop());
        Electron_missingHits_.push_back(iter->missingHits());
        Electron_charge_.push_back(iter->charge());
        Electron_ecalIso_.push_back(iter->ecalIso());
        Electron_hcalIso_.push_back(iter->hcalIso());
        Electron_trackIso_.push_back(iter->trackIso());
        Electron_r9_.push_back(iter->r9());
        Electron_sMin_.push_back(iter->sMin());
        Electron_sMaj_.push_back(iter->sMaj());
        Electron_seedId_.push_back(iter->seedId());
        //Electron_detIds_.push_back(vector<UInt_t>(iter->detIds()));
        Electron_energyMatrix_.push_back(vector<Float_t>(iter->energyMatrix()));
        Electron_timingMatrix_.push_back(vector<Float_t>(iter->timingMatrix()));
        //Electron_rechitZeroSuppression_.push_back(iter->rechitZeroSuppression());
        n_ele++;
    }

    n_mu=0;
    for (auto iter = muonsH->begin(); iter != muonsH->end(); ++iter) {
        Muon_pt_.push_back(iter->pt());
        Muon_eta_.push_back(iter->eta());
        Muon_phi_.push_back(iter->phi());
        Muon_m_.push_back(iter->m());
        Muon_type_.push_back(iter->type());
        Muon_charge_.push_back(iter->charge());
        Muon_normalizedChi2_.push_back(iter->normalizedChi2());
        Muon_ecalIso_.push_back(iter->ecalIso());
        Muon_hcalIso_.push_back(iter->hcalIso());
        Muon_trackIso_.push_back(iter->trackIso());
        Muon_nValidStandAloneMuonHits_.push_back(iter->nValidStandAloneMuonHits());
        Muon_nStandAloneMuonMatchedStations_.push_back(iter->nStandAloneMuonMatchedStations());
        Muon_nValidRecoMuonHits_.push_back(iter->nValidRecoMuonHits());
        Muon_nRecoMuonChambers_.push_back(iter->nRecoMuonChambers());
        Muon_nRecoMuonChambersCSCorDT_.push_back(iter->nRecoMuonChambersCSCorDT());
        Muon_nRecoMuonMatches_.push_back(iter->nRecoMuonMatches());
        Muon_nRecoMuonMatchedStations_.push_back(iter->nRecoMuonMatchedStations());
        Muon_nRecoMuonExpectedMatchedStations_.push_back(iter->nRecoMuonExpectedMatchedStations());
        Muon_recoMuonStationMask_.push_back(iter->recoMuonStationMask());
        Muon_nRecoMuonMatchedRPCLayers_.push_back(iter->nRecoMuonMatchedRPCLayers());
        Muon_recoMuonRPClayerMask_.push_back(iter->recoMuonRPClayerMask());
        Muon_nValidPixelHits_.push_back(iter->nValidPixelHits());
        Muon_nValidStripHits_.push_back(iter->nValidStripHits());
        Muon_nPixelLayersWithMeasurement_.push_back(iter->nPixelLayersWithMeasurement());
        Muon_nTrackerLayersWithMeasurement_.push_back(iter->nTrackerLayersWithMeasurement());
        Muon_trk_chi2_.push_back(iter->trk_chi2());
        Muon_trk_ndof_.push_back(iter->trk_ndof());
        Muon_trk_dxy_.push_back(iter->trk_dxy());
        Muon_trk_dz_.push_back(iter->trk_dz());
        Muon_trk_qoverp_.push_back(iter->trk_qoverp());
        Muon_trk_lambda_.push_back(iter->trk_lambda());
        Muon_trk_pt_.push_back(iter->trk_pt());
        Muon_trk_phi_.push_back(iter->trk_phi());
        Muon_trk_eta_.push_back(iter->trk_eta());
        Muon_trk_dxyError_.push_back(iter->trk_dxyError());
        Muon_trk_dzError_.push_back(iter->trk_dzError());
        Muon_trk_qoverpError_.push_back(iter->trk_qoverpError());
        Muon_trk_lambdaError_.push_back(iter->trk_lambdaError());
        Muon_trk_phiError_.push_back(iter->trk_phiError());
        Muon_trk_dsz_.push_back(iter->trk_dsz());
        Muon_trk_dszError_.push_back(iter->trk_dszError());
        Muon_trk_qoverp_lambda_cov_.push_back(iter->trk_qoverp_lambda_cov());
        Muon_trk_qoverp_phi_cov_.push_back(iter->trk_qoverp_phi_cov());
        Muon_trk_qoverp_dxy_cov_.push_back(iter->trk_qoverp_dxy_cov());
        Muon_trk_qoverp_dsz_cov_.push_back(iter->trk_qoverp_dsz_cov());
        Muon_trk_lambda_phi_cov_.push_back(iter->trk_lambda_phi_cov());
        Muon_trk_lambda_dxy_cov_.push_back(iter->trk_lambda_dxy_cov());
        Muon_trk_lambda_dsz_cov_.push_back(iter->trk_lambda_dsz_cov());
        Muon_trk_phi_dxy_cov_.push_back(iter->trk_phi_dxy_cov());
        Muon_trk_phi_dsz_cov_.push_back(iter->trk_phi_dsz_cov());
        Muon_trk_dxy_dsz_cov_.push_back(iter->trk_dxy_dsz_cov());
        Muon_trk_vx_.push_back(iter->trk_vx());
        Muon_trk_vy_.push_back(iter->trk_vy());
        Muon_trk_vz_.push_back(iter->trk_vz());
        Muon_vtxIndx_.push_back(vector<Int_t>(iter->vtxIndx()));
        n_mu++;
    }


    n_pho = 0;
    for (auto iter = photonsH->begin(); iter != photonsH->end(); ++iter) {
        Photon_pt_.push_back(iter->pt());
        Photon_eta_.push_back(iter->eta());
        Photon_phi_.push_back(iter->phi());
        Photon_m_.push_back(iter->m());
        Photon_sigmaIetaIeta_.push_back(iter->sigmaIetaIeta());
        Photon_hOverE_.push_back(iter->hOverE());
        Photon_ecalIso_.push_back(iter->ecalIso());
        Photon_hcalIso_.push_back(iter->hcalIso());
        Photon_trkIso_.push_back(iter->trkIso());
        Photon_r9_.push_back(iter->r9());
        Photon_sMin_.push_back(iter->sMin());
        Photon_sMaj_.push_back(iter->sMaj());
        Photon_seedId_.push_back(iter->seedId());
        //Photon_detIds_.push_back(vector<UInt_t>(iter->detIds()));
        Photon_energyMatrix_.push_back(vector<Float_t>(iter->energyMatrix()));
        Photon_timingMatrix_.push_back(vector<Float_t>(iter->timingMatrix()));    
        //Photon_rechitZeroSuppression_.push_back(iter->rechitZeroSuppression());
        n_pho++;
    }

    n_jet = 0;
    for (auto iter = pfjetsH->begin(); iter != pfjetsH->end(); ++iter) {
        Jet_pt_.push_back(iter->pt());
        Jet_eta_.push_back(iter->eta());
        Jet_phi_.push_back(iter->phi());
        Jet_m_.push_back(iter->m());
        Jet_area_.push_back(iter->jetArea());
        Jet_chargedHadronEnergy_.push_back(iter->chargedHadronEnergy());
        Jet_neutralHadronEnergy_.push_back(iter->neutralHadronEnergy());
        Jet_photonEnergy_.push_back(iter->photonEnergy());
        Jet_electronEnergy_.push_back(iter->electronEnergy());
        Jet_muonEnergy_.push_back(iter->muonEnergy());
        Jet_HFHadronEnergy_.push_back(iter->HFHadronEnergy());
        Jet_HFEMEnergy_.push_back(iter->HFEMEnergy());
        Jet_chargedHadronMultiplicity_.push_back(iter->chargedHadronMultiplicity());
        Jet_neutralHadronMultiplicity_.push_back(iter->neutralHadronMultiplicity());
        Jet_photonMultiplicity_.push_back(iter->photonMultiplicity());
        Jet_electronMultiplicity_.push_back(iter->electronMultiplicity());
        Jet_muonMultiplicity_.push_back(iter->muonMultiplicity());
        Jet_HFHadronMultiplicity_.push_back(iter->HFHadronMultiplicity());
        Jet_HFEMMultiplicity_.push_back(iter->HFEMMultiplicity());
        Jet_HOEnergy_.push_back(iter->HOEnergy());
        Jet_csv_.push_back(iter->csv());
        Jet_mvaDiscriminator_.push_back(iter->mvaDiscriminator());
        Jet_constituents_.push_back(vector<Int_t>(iter->constituents()));
        n_jet++;
    }

    n_pfcand = 0;
    for (auto iter = pfcandsH->begin(); iter != pfcandsH->end(); ++iter) {
        PFCand_pt_.push_back(iter->pt());
        PFCand_eta_.push_back(iter->eta());
        PFCand_phi_.push_back(iter->phi());
        PFCand_pdgId_.push_back(iter->pdgId());
        PFCand_vertex_.push_back(iter->vertex());
        n_pfcand++;
    }

    //n_track = 0;
    //for (auto iter = tracksH->begin(); iter != tracksH->end(); ++iter) {
    //    Track_pt_.push_back(iter->tk_pt());
    //    Track_eta_.push_back(iter->tk_eta());
    //    Track_phi_.push_back(iter->tk_phi());
    //    Track_chi2_.push_back(iter->tk_chi2());
    //    Track_ndof_.push_back(iter->tk_ndof());
    //    Track_charge_.push_back(iter->tk_charge());
    //    Track_dxy_.push_back(iter->tk_dxy());
    //    Track_dz_.push_back(iter->tk_dz());
    //    Track_nValidPixelHits_.push_back(iter->tk_nValidPixelHits());
    //    Track_nTrackerLayersWithMeasurement_.push_back(iter->tk_nTrackerLayersWithMeasurement());
    //    Track_nValidStripHits_.push_back(iter->tk_nValidStripHits());
    //    Track_qoverp_.push_back(iter->tk_qoverp());
    //    Track_lambda_.push_back(iter->tk_lambda());
    //    Track_dxy_Error_.push_back(iter->tk_dxy_Error());
    //    Track_dz_Error_.push_back(iter->tk_dz_Error());
    //    Track_qoverp_Error_.push_back(iter->tk_qoverp_Error());
    //    Track_lambda_Error_.push_back(iter->tk_lambda_Error());
    //    Track_phi_Error_.push_back(iter->tk_phi_Error());
    //    Track_dsz_.push_back(iter->tk_dsz());
    //    Track_dsz_Error_.push_back(iter->tk_dsz_Error());
    //    Track_qoverp_lambda_cov_.push_back(iter->tk_qoverp_lambda_cov());
    //    Track_qoverp_phi_cov_.push_back(iter->tk_qoverp_phi_cov());
    //    Track_qoverp_dxy_cov_.push_back(iter->tk_qoverp_dxy_cov());
    //    Track_qoverp_dsz_cov_.push_back(iter->tk_qoverp_dsz_cov());
    //    Track_lambda_phi_cov_.push_back(iter->tk_lambda_phi_cov());
    //    Track_lambda_dxy_cov_.push_back(iter->tk_lambda_dxy_cov());
    //    Track_lambda_dsz_cov_.push_back(iter->tk_lambda_dsz_cov());
    //    Track_phi_dxy_cov_.push_back(iter->tk_phi_dxy_cov());
    //    Track_phi_dsz_cov_.push_back(iter->tk_phi_dsz_cov());
    //    Track_dxy_dsz_cov_.push_back(iter->tk_dxy_dsz_cov());
    //    Track_vtxInd_.push_back(iter->tk_vtxInd());
    //    Track_vx_.push_back(iter->tk_vx());
    //    Track_vy_.push_back(iter->tk_vy());
    //    Track_vz_.push_back(iter->tk_vz());
    //    n_track++;
    //} 

    n_primaryvtx = 0;
    for (auto iter = pvH->begin(); iter != pvH->end(); ++iter) {
        PrimaryVtx_x_.push_back(iter->x());
        PrimaryVtx_y_.push_back(iter->y());
        PrimaryVtx_z_.push_back(iter->z());
        PrimaryVtx_zError_.push_back(iter->zError());
        PrimaryVtx_xError_.push_back(iter->xError());
        PrimaryVtx_yError_.push_back(iter->yError());
        PrimaryVtx_tracksSize_.push_back(iter->tracksSize());
        PrimaryVtx_chi2_.push_back(iter->chi2());
        PrimaryVtx_ndof_.push_back(iter->ndof());
        PrimaryVtx_isValidVtx_.push_back(iter->isValidVtx());
        n_primaryvtx++;
    }

    n_displacedvtx = 0;
    for (auto iter = pvH->begin(); iter != pvH->end(); ++iter) {
        DisplacedVtx_x_.push_back(iter->x());
        DisplacedVtx_y_.push_back(iter->y());
        DisplacedVtx_z_.push_back(iter->z());
        DisplacedVtx_zError_.push_back(iter->zError());
        DisplacedVtx_xError_.push_back(iter->xError());
        DisplacedVtx_yError_.push_back(iter->yError());
        DisplacedVtx_tracksSize_.push_back(iter->tracksSize());
        DisplacedVtx_chi2_.push_back(iter->chi2());
        DisplacedVtx_ndof_.push_back(iter->ndof());
        DisplacedVtx_isValidVtx_.push_back(iter->isValidVtx());
        n_displacedvtx++;
    }


    tree->Fill();	
    clearVars();
	
}

void ScoutingNanoAOD::clearVars(){
    l1Result_.clear();
    Photon_pt_.clear();
    Photon_eta_.clear();
    Photon_phi_.clear();
    Photon_m_.clear();
    Photon_sigmaIetaIeta_.clear();
    Photon_hOverE_.clear();
    Photon_ecalIso_.clear();
    Photon_hcalIso_.clear();
    Photon_trkIso_.clear();
    Photon_r9_.clear();
    Photon_sMin_.clear();
    Photon_sMaj_.clear();
    Photon_seedId_.clear();
    //Photon_detIds_.clear();
    Photon_energyMatrix_.clear();
    Photon_timingMatrix_.clear();
    //Photon_rechitZeroSuppression_.clear();
    Electron_pt_.clear();
    Electron_eta_.clear();
    Electron_phi_.clear();
    Electron_m_.clear();
    Electron_d0_.clear();
    Electron_dz_.clear();
    Electron_dEtaIn_.clear();
    Electron_dPhiIn_.clear();
    Electron_sigmaIetaIeta_.clear();
    Electron_hOverE_.clear();
    Electron_ooEMOop_.clear();
    Electron_missingHits_.clear();
    Electron_charge_.clear();
    Electron_ecalIso_.clear();
    Electron_hcalIso_.clear();
    Electron_trackIso_.clear();
    Electron_r9_.clear();
    Electron_sMin_.clear();
    Electron_sMaj_.clear();
    Electron_seedId_.clear();
    //Electron_detIds_.clear();
    Electron_energyMatrix_.clear();
    Electron_timingMatrix_.clear();
    //Electron_rechitZeroSuppression_.clear();
    Muon_pt_.clear();
    Muon_eta_.clear();
    Muon_phi_.clear();
    Muon_m_.clear();
    Muon_type_.clear();
    Muon_charge_.clear();
    Muon_normalizedChi2_.clear();
    Muon_ecalIso_.clear();
    Muon_hcalIso_.clear();
    Muon_trackIso_.clear();
    Muon_nValidStandAloneMuonHits_.clear();
    Muon_nStandAloneMuonMatchedStations_.clear();
    Muon_nValidRecoMuonHits_.clear();
    Muon_nRecoMuonChambers_.clear();
    Muon_nRecoMuonChambersCSCorDT_.clear();
    Muon_nRecoMuonMatches_.clear();
    Muon_nRecoMuonMatchedStations_.clear();
    Muon_nRecoMuonExpectedMatchedStations_.clear();
    Muon_recoMuonStationMask_.clear();
    Muon_nRecoMuonMatchedRPCLayers_.clear();
    Muon_recoMuonRPClayerMask_.clear();
    Muon_nValidPixelHits_.clear();
    Muon_nValidStripHits_.clear();
    Muon_nPixelLayersWithMeasurement_.clear();
    Muon_nTrackerLayersWithMeasurement_.clear();
    Muon_trk_chi2_.clear();
    Muon_trk_ndof_.clear();
    Muon_trk_dxy_.clear();
    Muon_trk_dz_.clear();
    Muon_trk_qoverp_.clear();
    Muon_trk_lambda_.clear();
    Muon_trk_pt_.clear();
    Muon_trk_phi_.clear();
    Muon_trk_eta_.clear();
    Muon_trk_dxyError_.clear();
    Muon_trk_dzError_.clear();
    Muon_trk_qoverpError_.clear();
    Muon_trk_lambdaError_.clear();
    Muon_trk_phiError_.clear();
    Muon_trk_dsz_.clear();
    Muon_trk_dszError_.clear();
    Muon_trk_qoverp_lambda_cov_.clear();
    Muon_trk_qoverp_phi_cov_.clear();
    Muon_trk_qoverp_dxy_cov_.clear();
    Muon_trk_qoverp_dsz_cov_.clear();
    Muon_trk_lambda_phi_cov_.clear();
    Muon_trk_lambda_dxy_cov_.clear();
    Muon_trk_lambda_dsz_cov_.clear();
    Muon_trk_phi_dxy_cov_.clear();
    Muon_trk_phi_dsz_cov_.clear();
    Muon_trk_dxy_dsz_cov_.clear();
    Muon_trk_vx_.clear();
    Muon_trk_vy_.clear();
    Muon_trk_vz_.clear();
    Muon_vtxIndx_.clear();
    Jet_pt_.clear();
    Jet_eta_.clear();
    Jet_phi_.clear();
    Jet_m_.clear();
    Jet_area_.clear();
    Jet_chargedHadronEnergy_.clear();
    Jet_neutralHadronEnergy_.clear();
    Jet_photonEnergy_.clear();
    Jet_electronEnergy_.clear();
    Jet_muonEnergy_.clear();
    Jet_HFHadronEnergy_.clear();
    Jet_HFEMEnergy_.clear();
    Jet_chargedHadronMultiplicity_.clear();
    Jet_neutralHadronMultiplicity_.clear();
    Jet_photonMultiplicity_.clear();
    Jet_electronMultiplicity_.clear();
    Jet_muonMultiplicity_.clear();
    Jet_HFHadronMultiplicity_.clear();
    Jet_HFEMMultiplicity_.clear();
    Jet_HOEnergy_.clear();
    Jet_csv_.clear();
    Jet_mvaDiscriminator_.clear();
    Jet_constituents_.clear();
    PFCand_pt_.clear();
    PFCand_eta_.clear();
    PFCand_phi_.clear();
    PFCand_pdgId_.clear();
    PFCand_vertex_.clear();
    //Track_pt_.clear();
    //Track_eta_.clear();
    //Track_phi_.clear();
    //Track_chi2_.clear();
    //Track_ndof_.clear();
    //Track_charge_.clear();
    //Track_dxy_.clear();
    //Track_dz_.clear();
    //Track_nValidPixelHits_.clear();
    //Track_nTrackerLayersWithMeasurement_.clear();
    //Track_nValidStripHits_.clear();
    //Track_qoverp_.clear();
    //Track_lambda_.clear();
    //Track_dxy_Error_.clear();
    //Track_dz_Error_.clear();
    //Track_qoverp_Error_.clear();
    //Track_lambda_Error_.clear();
    //Track_phi_Error_.clear();
    //Track_dsz_.clear();
    //Track_dsz_Error_.clear();
    //Track_qoverp_lambda_cov_.clear();
    //Track_qoverp_phi_cov_.clear();
    //Track_qoverp_dxy_cov_.clear();
    //Track_qoverp_dsz_cov_.clear();
    //Track_lambda_phi_cov_.clear();
    //Track_lambda_dxy_cov_.clear();
    //Track_lambda_dsz_cov_.clear();
    //Track_phi_dxy_cov_.clear();
    //Track_phi_dsz_cov_.clear();
    //Track_dxy_dsz_cov_.clear();
    //Track_vtxInd_.clear();
    //Track_vx_.clear();
    //Track_vy_.clear();
    //Track_vz_.clear();
    PrimaryVtx_x_.clear();
    PrimaryVtx_y_.clear();
    PrimaryVtx_z_.clear();
    PrimaryVtx_zError_.clear();
    PrimaryVtx_xError_.clear();
    PrimaryVtx_yError_.clear();
    PrimaryVtx_tracksSize_.clear();
    PrimaryVtx_chi2_.clear();
    PrimaryVtx_ndof_.clear();
    PrimaryVtx_isValidVtx_.clear();
    DisplacedVtx_x_.clear();
    DisplacedVtx_y_.clear();
    DisplacedVtx_z_.clear();
    DisplacedVtx_zError_.clear();
    DisplacedVtx_xError_.clear();
    DisplacedVtx_yError_.clear();
    DisplacedVtx_tracksSize_.clear();
    DisplacedVtx_chi2_.clear();
    DisplacedVtx_ndof_.clear();
    DisplacedVtx_isValidVtx_.clear(); 
}

void ScoutingNanoAOD::beginJob() {
  
}

void ScoutingNanoAOD::endJob() {
}

void ScoutingNanoAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}

void ScoutingNanoAOD::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingNanoAOD::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void ScoutingNanoAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void ScoutingNanoAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScoutingNanoAOD);
