// -*- C++ -*-
//
// Package:    PickDilepton
// Class:      PickDilepton
// 
/**\class PickDilepton PickDilepton.cc SUSYBSMAnalysis/Zprime2muAnalysis/PickDilepton.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/



// system include files
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include <boost/foreach.hpp>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/Dumpers.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "TLorentzVector.h"
//
// class declaration
//

class PickDilepton : public edm::EDFilter {
public:
  explicit PickDilepton(const edm::ParameterSet&);
  ~PickDilepton();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  //edm::EDGetTokenT< std::vector< pat::Muon > >  muonToken_;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;

  bool dump_muons;
  edm::InputTag muon_src;
  bool dump_dileptons;
  bool dump_dileptons_2;
  edm::InputTag dilepton_src;
  edm::InputTag dilepton_src_2;
  bool dump_met;
  std::vector<edm::InputTag> met_srcs;
  bool dump_jets;
  edm::InputTag jet_src;

  double massCut_;
  double ptMin_;
  double ptMax_;
  double isoMin_;
  double isoMax_;

  const reco::Vertex*   vertex;
  int                   nVtx;

  double mu_mass;

  bool debug;
};

// constructors and destructor
PickDilepton::PickDilepton(const edm::ParameterSet& cfg)
  : beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),

    dump_muons(cfg.existsAs<edm::InputTag>("muon_src")),

    dump_dileptons(cfg.existsAs<edm::InputTag>("dilepton_src")),

    dump_dileptons_2(cfg.existsAs<edm::InputTag>("dilepton_src_2")),

    dump_met(cfg.existsAs<edm::InputTag>("met_srcs")),

    dump_jets(cfg.existsAs<edm::InputTag>("jet_src")),

    vertex(0),
    nVtx(0)
{
  mu_mass = 0.105658;

  massCut_    = cfg.getParameter < double > ("massCut");
  ptMin_      = cfg.getParameter < double > ("ptMin");
  ptMax_      = cfg.getParameter < double > ("ptMax");
  isoMin_     = cfg.getParameter < double > ("isoMin");
  isoMax_     = cfg.getParameter < double > ("isoMax");

  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertex_src);

  std::ostringstream out;
  out << "configuration:\n";

  out << "dump_muons: " << dump_muons << " ";
  if (dump_muons) {
    muon_src = cfg.getParameter<edm::InputTag>("muon_src");
    consumes<pat::MuonCollection>(muon_src);
    out << muon_src;
  }
  out << "\n";

  out << "dump_dileptons: " << dump_dileptons << " ";
  if (dump_dileptons) {
    dilepton_src = cfg.getParameter<edm::InputTag>("dilepton_src");
    consumes<pat::CompositeCandidateCollection>(dilepton_src);
    out << dilepton_src;
  }
  out << "\n";

  out << "dump_dileptons_2: " << dump_dileptons_2 << " ";
  if (dump_dileptons_2) {
    dilepton_src_2 = cfg.getParameter<edm::InputTag>("dilepton_src_2");
    consumes<pat::CompositeCandidateCollection>(dilepton_src_2);
    out << dilepton_src_2;
  }
  out << "\n";

  out << "dump_met: " << dump_met << " ";
  if (dump_met) {
    met_srcs = cfg.getParameter<std::vector<edm::InputTag> >("met_srcs");
    BOOST_FOREACH(const edm::InputTag& met_src, met_srcs) {
      consumes<pat::METCollection>(met_src);
      out << met_src.encode() << " ";
    }
  }
  out << "\n";

  out << "dump_jets: " << dump_jets << " ";
  if (dump_jets) {
    jet_src = cfg.getParameter<edm::InputTag>("jet_src");
    consumes<pat::JetCollection>(jet_src);
    out << jet_src;
  }
  out << "\n";

  std::cout << out.str();

  debug = cfg.getParameter < bool > ("debug");
}

PickDilepton::~PickDilepton(){}


// member functions
// ------------ method called on each new Event  ------------
bool
PickDilepton::filter(edm::Event& event, const edm::EventSetup& iSetup)
{
  bool filter = false;

  edm::Handle<reco::BeamSpot> beamSpot;
  event.getByLabel(beamspot_src, beamSpot);

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel(vertex_src, vertices);

  vertex = 0;
  int vertex_count = 0;
  for (reco::VertexCollection::const_iterator it = vertices->begin(), ite = vertices->end(); it != ite; ++it) {
    if (it->ndof() > 4 && fabs(it->z()) <= 24 && fabs(it->position().rho()) <= 2) {
      if (vertex == 0)
        vertex = &*it;
      ++vertex_count;
    }
  }
  nVtx = vertex_count;

  if( dump_dileptons && dump_dileptons_2 ) {
    edm::Handle<pat::CompositeCandidateCollection> dils;
    event.getByLabel(dilepton_src, dils);

    edm::Handle<pat::CompositeCandidateCollection> dils_2;
    event.getByLabel(dilepton_src_2, dils_2);

    if ( !dils.isValid())
      std::cout << "WARNING! tried to get dileptons using " << dilepton_src << " and failed!\n";

    if ( !dils_2.isValid())
      std::cout << "WARNING! tried to get dileptons using " << dilepton_src_2 << " and failed!\n";


    if( (dils->size() == 0) && (dils_2->size() > 0) ) {

      BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils_2) {

        if(dil.mass() < massCut_ )
          continue;

        const reco::CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
        const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
        const reco::TrackRef TunTrack0 = mu0->tunePMuonBestTrack();
        double relTkIso0 = mu0->isolationR03().sumPt / mu0->innerTrack()->pt();

        const reco::CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);
        const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
        const reco::TrackRef TunTrack1 = mu1->tunePMuonBestTrack();
        double relTkIso1 = mu1->isolationR03().sumPt / mu1->innerTrack()->pt();

        if( ( TunTrack0->pt() > ptMin_ && TunTrack0->pt() < ptMax_ &&
              relTkIso0 > isoMin_ && relTkIso0 < isoMax_ &&
              mu0->numberOfMatchedStations() <= 1 ) ||
            ( TunTrack1->pt() > ptMin_ && TunTrack1->pt() < ptMax_ &&
              relTkIso1 > isoMin_ && relTkIso1 < isoMax_ &&
              mu1->numberOfMatchedStations() <= 1 )
          )
          filter = true;

        if(debug) {
          if( TunTrack0->pt() > ptMin_ && TunTrack0->pt() < ptMax_ && mu0->numberOfMatchedStations() <= 1 ) {
            std::cout << "mu0: pT=" << TunTrack0->pt() << ", eta=" << TunTrack0->eta() << ", phi=" << TunTrack0->phi() << std::endl;
            std::cout << "\t MS=" << mu0->numberOfMatchedStations() << std::endl;
            std::cout << "\tiso=" << relTkIso0 << std::endl;
          }
          if( TunTrack1->pt() > ptMin_ && TunTrack1->pt() < ptMax_ && mu1->numberOfMatchedStations() <= 1 ) {
            std::cout << "mu1: pT=" << TunTrack1->pt() << ", eta=" << TunTrack1->eta() << ", phi=" << TunTrack1->phi() << std::endl;
            std::cout << "\t MS=" << mu1->numberOfMatchedStations() << std::endl;
            std::cout << "\tiso=" << relTkIso1 << std::endl;
          }
        }

      }

    }
  }

  return filter;
}

// ------------ method called once each job just before starting event loop  ------------
void 
PickDilepton::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PickDilepton::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PickDilepton);
