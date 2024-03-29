#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

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

#include "Run3ScoutingJetTagging/Analysis/interface/FatJetMatching.h"

using namespace deepntuples;

FatJetMatching<fastjet::PseudoJet> ak8_pdseudojet_match;



class AK4JetNtupleProducer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit AK4JetNtupleProducer(const edm::ParameterSet&);
  ~AK4JetNtupleProducer();
		
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
  bool isNeutralPdg(int);

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> pfcand_token_;
  const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genjet_token_;
  const edm::EDGetTokenT<reco::GenJetCollection> fgenjet_token_;
  const edm::EDGetTokenT<reco::GenParticleCollection> gencand_token_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muon_token_;
  const edm::EDGetTokenT<double>  	pfMetToken;
  const edm::EDGetTokenT<double>  	pfMetPhiToken;
  const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particletable_token_;
  const edm::EDGetTokenT<std::vector<reco::GenJet>>  genjet_data_token_;
  
  edm::Handle<std::vector<Run3ScoutingParticle>> pfcands_;
  edm::Handle<reco::GenParticleCollection> gencands_;
  edm::Handle<reco::JetFlavourInfoMatchingCollection> genjets_;
  edm::Handle<reco::GenJetCollection> fgenjets_;
  edm::Handle<std::vector<Run3ScoutingMuon> > muon_;
  edm::Handle<double> pfMet_;
  edm::Handle<double> pfMetPhi_;
  edm::Handle<std::vector<reco::GenJet>> genjet_data_;

  TTree* tree;

  std::vector<Float16_t> pfcand_pt_log_nopuppi;
  std::vector<Float16_t> pfcand_e_log_nopuppi;
  std::vector<Float16_t> pfcand_absphi;
  std::vector<Float16_t> pfcand_abseta;
  std::vector<Int_t> pfcand_charge;
  std::vector<UInt_t> pfcand_isEl;
  std::vector<UInt_t> pfcand_isMu;
  std::vector<UInt_t> pfcand_isGamma;
  std::vector<UInt_t> pfcand_isChargedHad;
  std::vector<UInt_t> pfcand_isNeutralHad;
  std::vector<Float16_t> pfcand_lostInnerHits;
  std::vector<Float16_t> pfcand_normchi2;
  std::vector<Float16_t> pfcand_quality;
  std::vector<Float16_t> pfcand_dz;
  std::vector<Float16_t> pfcand_dzsig;
  std::vector<Float16_t> pfcand_dxy;
  std::vector<Float16_t> pfcand_dxysig;
  std::vector<Float16_t> pfcand_trk_pt;
  std::vector<Float16_t> pfcand_trk_eta;
  std::vector<Float16_t> pfcand_trk_phi;
  std::vector<std::vector<Float16_t>> pfcand_index;

  std::vector<Float16_t> j_pt;
  std::vector<Float16_t> j_eta;
  std::vector<Float16_t> j_phi;
  std::vector<Float16_t> j_mass;
  int j_no;
  int useless;
  std::vector<UInt_t> j_npfcands;

  std::vector<UInt_t> j_nCHadrons;
  std::vector<UInt_t> j_nBHadrons;
  std::vector<UInt_t> j_partonFlavour;
  std::vector<UInt_t> j_hadronFlavour;
  
  int sample_isQCD;
  int event_no;

  bool isQCD_ = false;

  float dR_= 0.4;
  float fdR_ = 0.8;



  std::vector<Float16_t> fgen_pt;
  std::vector<Float16_t> fgen_eta;
  std::vector<Float16_t> fgen_phi;
  std::vector<Float16_t> fgen_mass;
  int fgen_no;

  std::vector<Float16_t> gen_pt;
  std::vector<Float16_t> gen_eta;
  std::vector<Float16_t> gen_phi;
  std::vector<Float16_t> gen_mass;
  int gen_no;

  const static int 	max_mu = 1000;
  UInt_t n_mu;
  std::vector<Float_t> Muon_pt_;
  std::vector<Float_t> Muon_eta_;
  std::vector<Float_t> Muon_phi_;
  std::vector<Float_t> Muon_m_;
  std::vector<UInt_t> Muon_type_;
  std::vector<UInt_t> Muon_charge_;
  std::vector<Float_t> Muon_normalizedChi2_;
  std::vector<Float_t> Muon_ecalIso_;
  std::vector<Float_t> Muon_hcalIso_;
  std::vector<Float_t> Muon_trackIso_;
  std::vector<UInt_t> Muon_nValidStandAloneMuonHits_;
  std::vector<UInt_t> Muon_nStandAloneMuonMatchedStations_;
  std::vector<UInt_t> Muon_nValidRecoMuonHits_;
  std::vector<UInt_t> Muon_nRecoMuonChambers_;
  std::vector<UInt_t> Muon_nRecoMuonChambersCSCorDT_;
  std::vector<UInt_t> Muon_nRecoMuonMatches_;
  std::vector<UInt_t> Muon_nRecoMuonMatchedStations_;
  std::vector<UInt_t> Muon_nRecoMuonExpectedMatchedStations_;
  std::vector<UInt_t> Muon_recoMuonStationMask_;
  std::vector<UInt_t> Muon_nRecoMuonMatchedRPCLayers_;
  std::vector<UInt_t> Muon_recoMuonRPClayerMask_;
  std::vector<UInt_t> Muon_nValidPixelHits_;
  std::vector<UInt_t> Muon_nValidStripHits_;
  std::vector<UInt_t> Muon_nPixelLayersWithMeasurement_;
  std::vector<UInt_t> Muon_nTrackerLayersWithMeasurement_;
  std::vector<Float_t> Muon_trk_chi2_;
  std::vector<Float_t> Muon_trk_ndof_;
  std::vector<Float_t> Muon_trk_dxy_;
  std::vector<Float_t> Muon_trk_dz_;
  std::vector<Float_t> Muon_trk_qoverp_;
  std::vector<Float_t> Muon_trk_lambda_;
  std::vector<Float_t> Muon_trk_pt_;
  std::vector<Float_t> Muon_trk_phi_;
  std::vector<Float_t> Muon_trk_eta_;
  std::vector<Float_t> Muon_trk_dxyError_;
  std::vector<Float_t> Muon_trk_dzError_;
  std::vector<Float_t> Muon_trk_qoverpError_;
  std::vector<Float_t> Muon_trk_lambdaError_;
  std::vector<Float_t> Muon_trk_phiError_;
  std::vector<Float_t> Muon_trk_dsz_;
  std::vector<Float_t> Muon_trk_dszError_;
  std::vector<Float_t> Muon_trk_qoverp_lambda_cov_;
  std::vector<Float_t> Muon_trk_qoverp_phi_cov_;
  std::vector<Float_t> Muon_trk_qoverp_dxy_cov_;
  std::vector<Float_t> Muon_trk_qoverp_dsz_cov_;
  std::vector<Float_t> Muon_trk_lambda_phi_cov_;
  std::vector<Float_t> Muon_trk_lambda_dxy_cov_;
  std::vector<Float_t> Muon_trk_lambda_dsz_cov_;
  std::vector<Float_t> Muon_trk_phi_dxy_cov_;
  std::vector<Float_t> Muon_trk_phi_dsz_cov_;
  std::vector<Float_t> Muon_trk_dxy_dsz_cov_;
  std::vector<Float_t> Muon_trk_vx_;
  std::vector<Float_t> Muon_trk_vy_;
  std::vector<Float_t> Muon_trk_vz_;
  //vector<reco::HitPattern> trk_hitPattern_;
  std::vector<std::vector<Int_t> > Muon_vtxIndx_;
  //PF Met
  double pfMet;
  double pfMetPhi;

  std::vector<Float16_t> fj_pt;
  std::vector<Float16_t> fj_eta;
  std::vector<Float16_t> fj_phi;
  std::vector<Float16_t> fj_mass;
  std::vector<Float16_t> fj_msd;
  std::vector<Float16_t> fj_n2b1;
  int fj_no;
  
  std::vector<UInt_t> fj_npfcands;
  std::vector<std::vector<Float16_t> > f_pfcand_index;


  std::vector<Float16_t> fj_gen_mass;
  std::vector<Float16_t> fj_genjet_sdmass;

  std::vector<Float16_t> label_Top_bcq;
  std::vector<Float16_t> label_Top_bqq;
  std::vector<Float16_t> label_Top_bc;
  std::vector<Float16_t> label_Top_bq;
  std::vector<Float16_t> label_W_cq;
  std::vector<Float16_t> label_W_qq;
  std::vector<Float16_t> label_Z_bb;
  std::vector<Float16_t> label_Z_cc;
  std::vector<Float16_t> label_Z_qq;
  std::vector<Float16_t> label_H_bb;
  std::vector<Float16_t> label_H_cc;
  std::vector<Float16_t> label_H_qqqq;
  std::vector<Float16_t> label_H_tautau;
  std::vector<Float16_t> label_H_qq;
  std::vector<Float16_t> label_QCD_all;

};

AK4JetNtupleProducer::AK4JetNtupleProducer(const edm::ParameterSet& iConfig):
  pfcand_token_(consumes<std::vector<Run3ScoutingParticle>>(iConfig.getParameter<edm::InputTag>("pf_candidates"))),
  genjet_token_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("gen_jets"))),
  fgenjet_token_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("fgen_jets"))),
  gencand_token_(consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("gen_candidates"))),
  muon_token_(consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  pfMetToken             (consumes<double>(iConfig.getParameter<edm::InputTag>("pfMet"))),
  pfMetPhiToken             (consumes<double>(iConfig.getParameter<edm::InputTag>("pfMetPhi"))),
  particletable_token_(esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>()),
  genjet_data_token_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("gen_jet_data")))

{

  isQCD_ = iConfig.getUntrackedParameter<bool>("isQCD", false);

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("Events", "Events");

  tree->Branch("pfcand_pt_log_nopuppi", &pfcand_pt_log_nopuppi);
  tree->Branch("pfcand_e_log_nopuppi", &pfcand_e_log_nopuppi);
  tree->Branch("pfcand_absphi", &pfcand_absphi);
  tree->Branch("pfcand_abseta", &pfcand_abseta);
  tree->Branch("pfcand_charge", &pfcand_charge);
  tree->Branch("pfcand_isEl", &pfcand_isEl);
  tree->Branch("pfcand_isMu", &pfcand_isMu);
  tree->Branch("pfcand_isGamma", &pfcand_isGamma);
  tree->Branch("pfcand_isChargedHad", &pfcand_isChargedHad);
  tree->Branch("pfcand_isNeutralHad", &pfcand_isNeutralHad);
  tree->Branch("pfcand_lostInnerHits", &pfcand_lostInnerHits);
  tree->Branch("pfcand_normchi2", &pfcand_normchi2);
  tree->Branch("pfcand_quality", &pfcand_quality);
  tree->Branch("pfcand_dz", &pfcand_dz);
  tree->Branch("pfcand_dzsig", &pfcand_dzsig);
  tree->Branch("pfcand_dxy", &pfcand_dxy);
  tree->Branch("pfcand_dxysig", &pfcand_dxysig);
  tree->Branch("pfcand_trk_pt", &pfcand_trk_pt);
  tree->Branch("pfcand_trk_eta", &pfcand_trk_eta);
  tree->Branch("pfcand_trk_phi", &pfcand_trk_phi);
  tree->Branch("pfcand_index", &pfcand_index);

  tree->Branch("j_pt", &j_pt);
  tree->Branch("j_eta", &j_eta);
  tree->Branch("j_phi", &j_phi);
  tree->Branch("j_mass", &j_mass);
  tree->Branch("j_no", &j_no);
  tree->Branch("j_npfcands", &j_npfcands);

  tree->Branch("j_nCHadrons", &j_nCHadrons);
  tree->Branch("j_nBHadrons", &j_nBHadrons);
  tree->Branch("j_partonFlavour", &j_partonFlavour);
  tree->Branch("j_hadronFlavour", &j_hadronFlavour);
  tree->Branch("sample_isQCD", &sample_isQCD);

  tree->Branch("event_no", &event_no);

  tree->Branch("gen_pt", &gen_pt);
  tree->Branch("gen_eta", &gen_eta);
  tree->Branch("gen_phi", &gen_phi);
  tree->Branch("gen_mass", &gen_mass);
  tree->Branch("gen_no", &gen_no);

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

  tree->Branch("PFMet_Pt", &pfMet );
  tree->Branch("PFMet_Phi", &pfMetPhi );
  

  tree->Branch("fj_pt", &fj_pt);
  tree->Branch("fj_eta", &fj_eta);
  tree->Branch("fj_phi", &fj_phi);
  tree->Branch("fj_mass", &fj_mass);
  tree->Branch("fj_msd", &fj_msd);
  tree->Branch("fj_n2b1", &fj_n2b1);
  tree->Branch("fj_no", &fj_no);
  tree->Branch("fj_npfcands", &fj_npfcands);
  tree->Branch("f_pfcand_index", &f_pfcand_index);

  tree->Branch("fj_gen_mass", &fj_gen_mass);
  tree->Branch("fj_genjet_sdmass", &fj_genjet_sdmass);

  tree->Branch("fgen_pt", &fgen_pt);
  tree->Branch("fgen_eta", &fgen_eta);
  tree->Branch("fgen_phi", &fgen_phi);
  tree->Branch("fgen_mass", &fgen_mass);
  tree->Branch("fgen_no", &fgen_no);

  tree->Branch("label_Top_bcq", &label_Top_bcq);
  tree->Branch("label_Top_bqq", &label_Top_bqq);
  tree->Branch("label_Top_bc", &label_Top_bc);
  tree->Branch("label_Top_bq", &label_Top_bq);
  tree->Branch("label_W_cq", &label_W_cq);
  tree->Branch("label_W_qq", &label_W_qq);
  tree->Branch("label_Z_bb", &label_Z_bb);
  tree->Branch("label_Z_cc", &label_Z_cc);
  tree->Branch("label_Z_qq", &label_Z_qq);
  tree->Branch("label_H_bb", &label_H_bb);
  tree->Branch("label_H_cc", &label_H_cc);
  tree->Branch("label_H_qqqq", &label_H_qqqq);
  tree->Branch("label_H_tautau", &label_H_tautau);
  tree->Branch("label_H_qq", &label_H_qq);
  tree->Branch("label_QCD_all", &label_QCD_all);

}

AK4JetNtupleProducer::~AK4JetNtupleProducer() {
}

bool AK4JetNtupleProducer::isNeutralPdg(int pdgId) {
   const int neutralPdgs_array[] = {9, 21, 22, 23, 25, 12, 14, 16, 111, 130, 310, 311, 421, 511, 2112}; // gluon, gluon, gamma, Z0, higgs, electron neutrino, muon neutrino, tau neutrino, pi0, K0_L, K0_S; K0, neutron
   const std::vector<int> neutralPdgs(neutralPdgs_array, neutralPdgs_array + sizeof(neutralPdgs_array) / sizeof(int));
   if (std::find(neutralPdgs.begin(), neutralPdgs.end(), std::abs(pdgId)) == neutralPdgs.end())
     return false;
 
   return true;
}

void AK4JetNtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // PF candidates
  iEvent.getByToken(pfcand_token_, pfcands_);
  // GEN candidates
  iEvent.getByToken(gencand_token_, gencands_);
  // GEN jets
  iEvent.getByToken(genjet_token_, genjets_);
  //Fat Gen Jets
  iEvent.getByToken(fgenjet_token_, fgenjets_);
  // Particle table
  auto pdt = iSetup.getHandle(particletable_token_);
  const HepPDT::ParticleDataTable* particletable_ = pdt.product();

  //Muon
  iEvent.getByToken(muon_token_, muon_);

  //MeT
  iEvent.getByToken(pfMetToken, pfMet_);
  pfMet = *pfMet_;
  iEvent.getByToken(pfMetPhiToken, pfMetPhi_);
  pfMetPhi = *pfMetPhi_;

  iEvent.getByToken(genjet_data_token_, genjet_data_);

  //int pfcand_i = 0;
  for (auto pfcands_iter = pfcands_->begin(); pfcands_iter != pfcands_->end(); ++pfcands_iter) {
    //auto *reco_cand = dynamic_cast<const Run3ScoutingParticle*> (&pfcands_->at(cand.user_index()));

    //const reco::PFCandidate* pfcand = &(*pfcands_iter);

    //auto *reco_cand = dynamic_cast<const Run3ScoutingParticle*> (pfcands_iter->originalObject());
    
    //std::cout<<reco_cand->trk_pt()<<"\n";
    std::cout<<pfcands_iter->trk_eta()<<"\n";
    std::cout<<"=============="<<"\n";

    float pfcands_iter_p = pfcands_iter->pt() * cosh(pfcands_iter->eta());
    auto rcm = particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId())) != nullptr
                        ? particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId()))->mass()
                        : -99.f;
      if (rcm < -90) continue;

    pfcand_e_log_nopuppi.push_back(log(sqrt(pfcands_iter_p*pfcands_iter_p + rcm*rcm)));
    pfcand_pt_log_nopuppi.push_back(log(pfcands_iter->pt()));
    pfcand_absphi.push_back(pfcands_iter->phi());
    pfcand_abseta.push_back(pfcands_iter->eta());
    if (isNeutralPdg(pfcands_iter->pdgId())) {
      pfcand_charge.push_back(0);
    } else {
       pfcand_charge.push_back(abs(pfcands_iter->pdgId())/pfcands_iter->pdgId());
    }
    pfcand_isEl.push_back(abs(pfcands_iter->pdgId()) == 11);
    pfcand_isMu.push_back(abs(pfcands_iter->pdgId()) == 13);
    pfcand_isGamma.push_back(abs(pfcands_iter->pdgId()) == 22);
    pfcand_isChargedHad.push_back(abs(pfcands_iter->pdgId()) == 211);
    pfcand_isNeutralHad.push_back(abs(pfcands_iter->pdgId()) == 130);
    pfcand_lostInnerHits.push_back(pfcands_iter->lostInnerHits());
    pfcand_normchi2.push_back(pfcands_iter->normchi2());
    pfcand_quality.push_back(pfcands_iter->quality());
    pfcand_dz.push_back(pfcands_iter->dz());
    pfcand_dzsig.push_back(pfcands_iter->dzsig());
    pfcand_dxy.push_back(pfcands_iter->dxy());
    pfcand_dxysig.push_back(pfcands_iter->dxysig());
    pfcand_trk_pt.push_back(pfcands_iter->trk_pt() );
    pfcand_trk_eta.push_back(pfcands_iter->trk_eta());
    pfcand_trk_phi.push_back(pfcands_iter->trk_phi());
  }

  fgen_no = 0;
  for (auto fgenjets_iter = fgenjets_->begin(); fgenjets_iter != fgenjets_->end(); ++fgenjets_iter) {
    fgen_pt.push_back(fgenjets_iter->pt());
    fgen_eta.push_back(fgenjets_iter->eta());
    fgen_phi.push_back(fgenjets_iter->phi());
    fgen_mass.push_back(fgenjets_iter->mass());
    fgen_no ++;
  }


  gen_no = 0;
  for (auto genjets_iter = genjet_data_->begin(); genjets_iter != genjet_data_->end(); ++genjets_iter) {
    gen_pt.push_back(genjets_iter->pt());
    gen_eta.push_back(genjets_iter->eta());
    gen_phi.push_back(genjets_iter->phi());
    gen_mass.push_back(genjets_iter->mass());
    gen_no ++;
  }

  // Create jets
  std::vector<fastjet::PseudoJet> j_part;
  j_part.reserve(pfcands_->size());
  int pfcand_i = 0;
  for (auto pfcands_iter = pfcands_->begin(); pfcands_iter != pfcands_->end(); ++pfcands_iter) {

    auto pfm = particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId())) != nullptr
                        ? particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId()))->mass()
                        : -99.f;
    if (pfm < -90) continue;

    math::PtEtaPhiMLorentzVector p4(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfm);
    j_part.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
    j_part.back().set_user_index(pfcand_i);
    pfcand_i++;
  }

  fastjet::JetDefinition ak4_def = fastjet::JetDefinition(fastjet::antikt_algorithm, dR_);
  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  fastjet::ClusterSequenceArea ak4_cs(j_part, ak4_def, area_def);
  std::vector<fastjet::PseudoJet> ak4_jets = fastjet::sorted_by_pt(ak4_cs.inclusive_jets(15.0));

  // Match jet-gen jet
  std::map<int, reco::JetFlavourInfo> genmatch_resultmap;
  std::vector<int> genmatch_unmatched;
  std::vector<std::tuple<int, int, float> > pairlist;

  int ak4_jet_idx = 0;

  for(unsigned int i=0; i<ak4_jets.size(); i++) {
    bool found_match = false;
    for(unsigned int j=0; j<genjets_->size(); j++) {
      float dR = reco::deltaR(ak4_jets[i].eta(), ak4_jets[i].phi(), (*genjets_)[j].first.get()->eta(), (*genjets_)[j].first.get()->phi());
      if(dR < dR_) {
        pairlist.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    ak4_jet_idx++;
    if(!found_match) {
       genmatch_unmatched.push_back(i);
    }
  }

  std::sort(pairlist.begin(), pairlist.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(pairlist.size() > 0) {
    reco::JetFlavourInfo genjet_assn = (*genjets_)[std::get<1>(pairlist[0])].second;
    genmatch_resultmap[std::get<0>(pairlist[0])] = genjet_assn;
    for(unsigned int k=1; k<pairlist.size(); k++) {
      if(std::get<0>(pairlist[k]) == std::get<0>(pairlist[0]) ||
         std::get<1>(pairlist[k]) == std::get<1>(pairlist[0])) {
        pairlist.erase(pairlist.begin() + k);
      }
    }
    pairlist.erase(pairlist.begin());
  }

  // Create AK8 Jet
  std::vector<fastjet::PseudoJet> fj_part;
  fj_part.reserve(pfcands_->size());
  int fpfcand_i = 0;
  for (auto pfcands_iter = pfcands_->begin(); pfcands_iter != pfcands_->end(); ++pfcands_iter) {

    auto pfm = particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId())) != nullptr
                        ? particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId()))->mass()
                        : -99.f;
    if (pfm < -90) continue;

    math::PtEtaPhiMLorentzVector p4(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfm);
    fj_part.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
    fj_part.back().set_user_index(fpfcand_i);
    fpfcand_i++;
  }

  fastjet::JetDefinition ak8_def = fastjet::JetDefinition(fastjet::antikt_algorithm, fdR_);
  fastjet::GhostedAreaSpec farea_spec(5.0,1,0.01);
  fastjet::AreaDefinition farea_def(fastjet::active_area, farea_spec);

  // Substructure
  double sd_z_cut = 0.10;
  double sd_beta = 0;
  fastjet::contrib::SoftDrop sd_groomer = fastjet::contrib::SoftDrop(sd_beta, sd_z_cut, fdR_);
  fastjet::contrib::EnergyCorrelatorN2 N2 = fastjet::contrib::EnergyCorrelatorN2(1.0);

  fastjet::ClusterSequenceArea ak8_cs(fj_part, ak8_def, farea_def);
  std::vector<fastjet::PseudoJet> ak8_jets = fastjet::sorted_by_pt(ak8_cs.inclusive_jets(170.0));

  // Match jet-gen jet
  std::map<int, reco::GenJet> fgenmatch_resultmap;
  std::vector<int> fgenmatch_unmatched;
  std::vector<std::tuple<int, int, float> > parlist;

  int ak8_jet_idx = 0;
  //good from here

  for(unsigned int i=0; i<ak8_jets.size(); i++) {
    bool ffound_match = false;
    for(unsigned int j=0; j<fgenjets_->size(); j++) {
      float fdR = reco::deltaR(ak8_jets[i].eta(), ak8_jets[i].phi(), (*fgenjets_)[j].eta(), (*fgenjets_)[j].phi());
      if(fdR < fdR_) {
        parlist.push_back(std::make_tuple(i, j, fdR));
        ffound_match = true;
      }
    }
    ak8_jet_idx++;
    if(!ffound_match) {
       fgenmatch_unmatched.push_back(i);
    }
  }

  std::sort(parlist.begin(), parlist.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(parlist.size() > 0) {
    reco::GenJet fgenjet_assn = (*fgenjets_)[std::get<1>(parlist[0])];
    fgenmatch_resultmap[std::get<0>(parlist[0])] = fgenjet_assn;
    for(unsigned int k=1; k<parlist.size(); k++) {
      if(std::get<0>(parlist[k]) == std::get<0>(parlist[0]) ||
         std::get<1>(parlist[k]) == std::get<1>(parlist[0])) {
        parlist.erase(parlist.begin() + k);
      }
    }
    parlist.erase(parlist.begin());
  }


  n_mu=0;
  for (auto iter = muon_->begin(); iter != muon_->end(); ++iter) {
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
    Muon_vtxIndx_.push_back(std::vector<Int_t>(iter->vtxIndx()));
    n_mu++;
  }
  std::cout<<"***************"<<"\n";


  
  ak4_jet_idx = 0;
  j_no = 0;
  for(auto &j: ak4_jets) {

    if (abs(j.eta())< 2.5){
      j_no++;
      const std::vector<fastjet::PseudoJet> constituents = j.constituents();
      std::vector<Float16_t> user_index_vector;
      for (auto &cand : constituents) {
        user_index_vector.push_back(cand.user_index());
        auto *reco_cand = dynamic_cast<const Run3ScoutingParticle*> (&pfcands_->at(cand.user_index()));
        std::cout<<reco_cand->trk_eta()<<"\n";
      }
      pfcand_index.push_back(user_index_vector);
      user_index_vector.clear();

      j_pt.push_back(j.pt());
      j_eta.push_back(j.eta());
      j_phi.push_back(j.phi());
      j_mass.push_back(j.m());

      //should be pushed outside the loop to save ram
      //j_no = ak4_jets.size();
      j_npfcands.push_back(constituents.size());
      sample_isQCD = isQCD_;

      if (std::find(genmatch_unmatched.begin(), genmatch_unmatched.end(), ak4_jet_idx) != genmatch_unmatched.end()) {
        j_nCHadrons.push_back(50);
        j_nBHadrons.push_back(50);
        j_partonFlavour.push_back(50);
        j_hadronFlavour.push_back(50);
      } else {
        reco::JetFlavourInfo flavJet = genmatch_resultmap[ak4_jet_idx];
        j_nCHadrons.push_back(flavJet.getcHadrons().size());
        j_nBHadrons.push_back(flavJet.getbHadrons().size());
        j_partonFlavour.push_back(flavJet.getPartonFlavour());
        j_hadronFlavour.push_back(flavJet.getHadronFlavour());
      }

      event_no = iEvent.id().event();
    }
    else{
      useless = 0;
    }
    ak4_jet_idx++;
  }

  ak8_jet_idx = 0;
  fj_no = 0;
  for(auto &j: ak8_jets) {

    if (abs(j.eta())< 2.5){
      // Flavour matching
      fj_no++;
      auto ak8_label = ak8_pdseudojet_match.flavorLabel(j, *gencands_, fdR_);

      const std::vector<fastjet::PseudoJet> constituents = j.constituents();
      std::vector<Float16_t> f_user_index_vector;
      //std::vector<Float16_t> f_pfcand_trk_pt;
      for (auto &cand : constituents) {
        f_user_index_vector.push_back(cand.user_index());
        //auto *reco_cand = dynamic_cast<const Run3ScoutingParticle*> (&pfcands_->at(cand.user_index()));
        //f_pfcand_trk_pt.push_back(reco_cand->trk_pt() );

      }

      f_pfcand_index.push_back(f_user_index_vector);
      //pfcand_trk_pt.push_back(f_pfcand_trk_pt);
      f_user_index_vector.clear();

      
      fj_pt.push_back(j.pt());
      fj_eta.push_back(j.eta());
      fj_phi.push_back(j.phi());
      fj_mass.push_back(j.m());

      fastjet::PseudoJet sd_ak8 = sd_groomer(j);
      fj_msd.push_back(sd_ak8.m());
      fj_n2b1.push_back(N2(sd_ak8));

      //could be pushed out the loop
      //fj_no = ak8_jets.size();
      fj_npfcands.push_back(constituents.size());

      label_Top_bcq.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Top_bcq));
      label_Top_bqq.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Top_bqq));
      label_Top_bc.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Top_bc));
      label_Top_bq.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Top_bq));
      label_W_cq.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::W_cq));
      label_W_qq.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::W_qq));
      label_Z_bb.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Z_bb));
      label_Z_cc.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Z_cc));
      label_Z_qq.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Z_qq));
      label_H_bb.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_bb));
      label_H_cc.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_cc));
      label_H_qqqq.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_qqqq));
      label_H_tautau.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_tautau));
      label_H_qq.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_qq));
      label_QCD_all.push_back((ak8_label.first == FatJetMatching<fastjet::PseudoJet>::QCD_all));

      fj_gen_mass.push_back((ak8_label.first < FatJetMatching<fastjet::PseudoJet>::QCD_all && ak8_label.second) ? ak8_label.second->mass() : 0);

      if (std::find(fgenmatch_unmatched.begin(), fgenmatch_unmatched.end(), ak8_jet_idx) != fgenmatch_unmatched.end()) {
        fj_genjet_sdmass.push_back(-99);
      } else {
        fj_genjet_sdmass.push_back(fgenmatch_resultmap[ak8_jet_idx].mass());
      }
    }
    else {
      useless =0;
    }

    ak8_jet_idx++;
  }
  tree->Fill();
  clearVars();

}

void AK4JetNtupleProducer::clearVars(){
  pfcand_pt_log_nopuppi.clear();
  pfcand_e_log_nopuppi.clear();
  pfcand_absphi.clear();
  pfcand_abseta.clear();
  pfcand_charge.clear();
  pfcand_isEl.clear();
  pfcand_isMu.clear();
  pfcand_isGamma.clear();
  pfcand_isChargedHad.clear();
  pfcand_isNeutralHad.clear();
  pfcand_lostInnerHits.clear();
  pfcand_normchi2.clear();
  pfcand_quality.clear();
  pfcand_dz.clear();
  pfcand_dzsig.clear();
  pfcand_dxy.clear();
  pfcand_dxysig.clear();
  pfcand_index.clear();
  pfcand_trk_pt.clear();
  pfcand_trk_eta.clear();
  pfcand_trk_phi.clear();

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
  j_pt.clear();
  j_eta.clear();
  j_phi.clear();
  j_mass.clear();

  fgen_pt.clear();
  fgen_eta.clear();
  fgen_phi.clear();
  fgen_mass.clear();

  gen_pt.clear();
  gen_eta.clear();
  gen_phi.clear();
  gen_mass.clear();

  j_npfcands.clear();

  j_nCHadrons.clear();
  j_nBHadrons.clear();
  j_partonFlavour.clear();
  j_hadronFlavour.clear();

  fj_pt.clear();
  fj_eta.clear();
  fj_phi.clear();
  fj_mass.clear();
  fj_msd.clear();
  fj_n2b1.clear();

  fj_gen_mass.clear();
  fj_genjet_sdmass.clear();

  label_Top_bcq.clear();
  label_Top_bqq.clear();
  label_Top_bc.clear();
  label_Top_bq.clear();
  label_W_cq.clear();
  label_W_qq.clear();
  label_Z_bb.clear();
  label_Z_cc.clear();
  label_Z_qq.clear();
  label_H_bb.clear();
  label_H_cc.clear();
  label_H_qqqq.clear();
  label_H_tautau.clear();
  label_H_qq.clear();
  label_QCD_all.clear();
  f_pfcand_index.clear();
  fj_npfcands.clear();

}

void AK4JetNtupleProducer::beginJob() {
}

void AK4JetNtupleProducer::endJob() {
}

void AK4JetNtupleProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}

void AK4JetNtupleProducer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void AK4JetNtupleProducer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void AK4JetNtupleProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void AK4JetNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(AK4JetNtupleProducer);
