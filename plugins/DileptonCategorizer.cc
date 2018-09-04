// -*- C++ -*-
//
// Package:    DileptonCategorizer
// Class:      DileptonCategorizer
// 
/**\class DileptonCategorizer DileptonCategorizer.cc SUSYBSMAnalysis/Zprime2muAnalysis/DileptonCategorizer.cc

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

class DileptonCategorizer : public edm::EDFilter {
public:
  explicit DileptonCategorizer(const edm::ParameterSet&);
  ~DileptonCategorizer();
  
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
  edm::InputTag dilepton_src;
  bool dump_met;
  std::vector<edm::InputTag> met_srcs;
  bool dump_jets;
  edm::InputTag jet_src;

  double massCut_;
  double ProbePtCut_;

  const reco::Vertex*   vertex;
  int                   nVtx;

  double mu_mass;

  int nProbeT;
  int nProbeB;
  int nProbeO;
  int nProbeE;
  int nProbeF;

  int nIsTunepPass;
  int nIsDYTPass;
  int nIsPickyPass;
  int nIsTPFMSPass;

  int nIsNStationMoreThan1;

  bool debug;
};

// constructors and destructor
DileptonCategorizer::DileptonCategorizer(const edm::ParameterSet& cfg)
  : beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),

    dump_muons(cfg.existsAs<edm::InputTag>("muon_src")),
    //muon_src(dump_muons ? cfg.getParameter<edm::InputTag>("muon_src") : ),

    dump_dileptons(cfg.existsAs<edm::InputTag>("dilepton_src")),
    //dilepton_src(dump_dileptons ? cfg.getParameter<edm::InputTag>("dilepton_src") : ),

    dump_met(cfg.existsAs<edm::InputTag>("met_srcs")),
    //met_srcs(dump_met ? cfg.getParameter<std::vector<edm::InputTag>>("met_srcs") : ),

    dump_jets(cfg.existsAs<edm::InputTag>("jet_src")),
    //jet_src(dump_jets ? cfg.getParameter<edm::InputTag>("jet_src") : ),

    vertex(0),
    nVtx(0)
{
  mu_mass = 0.105658;

  massCut_    = cfg.getParameter < double > ("massCut");
  ProbePtCut_ = cfg.getParameter < double > ("ProbePtCut");

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

DileptonCategorizer::~DileptonCategorizer(){}


// member functions
// ------------ method called on each new Event  ------------
bool
DileptonCategorizer::filter(edm::Event& event, const edm::EventSetup& iSetup)
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

  if( dump_dileptons ) {
    std::ostringstream out;
    edm::Handle<pat::CompositeCandidateCollection> dils;
    event.getByLabel(dilepton_src, dils);

    if (!dils.isValid())
      std::cout << "WARNING! tried to get dileptons using " << dilepton_src << " and failed!\n";
    else if(dils->size()) {
      BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils) {
        
        if(dil.mass() < massCut_ )
          continue;

        const reco::CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
        const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
        const reco::CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);
        const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);

        const pat::Muon *TagMu;
        const pat::Muon *ProbeMu;
        if(mu0->globalTrack()->hitPattern().numberOfValidMuonHits() > 0) {
          TagMu = mu0;
          ProbeMu = mu1;
        }
        else {
          TagMu = mu1;
          ProbeMu = mu0;
        }

        if(ProbeMu->pt() < ProbePtCut_ )
          continue;

        const reco::TrackRef GloTrackTag = TagMu->globalTrack();
        const reco::TrackRef InnTrackTag = TagMu->innerTrack();
        const reco::TrackRef StaTrackTag = TagMu->standAloneMuon();
        const reco::TrackRef TunTrackTag = TagMu->tunePMuonBestTrack();
        const reco::TrackRef DytTrackTag = TagMu->dytTrack();
        const reco::TrackRef PicTrackTag = TagMu->pickyTrack();
        const reco::TrackRef TpfTrackTag = TagMu->tpfmsTrack();

        const reco::TrackRef GloTrackProbe = ProbeMu->globalTrack();
        const reco::TrackRef InnTrackProbe = ProbeMu->innerTrack();
        const reco::TrackRef StaTrackProbe = ProbeMu->standAloneMuon();
        const reco::TrackRef TunTrackProbe = ProbeMu->tunePMuonBestTrack();
        const reco::TrackRef DytTrackProbe = ProbeMu->dytTrack();
        const reco::TrackRef PicTrackProbe = ProbeMu->pickyTrack();
        const reco::TrackRef TpfTrackProbe = ProbeMu->tpfmsTrack();

        out << "Mass : " << dil.mass() << "\n\n";
        //out << "lep0(pT, eta, phi) : " << lep0->pt() << " " << lep0->eta() << " " << lep0->phi() << "\n";
        out << "TagMu(pT, eta, phi) : " << TagMu->pt() << " " << TagMu->eta() << " " << TagMu->phi() << "\n";
        out << "# Mat St : " << TagMu->numberOfMatchedStations() << "\n";
        if( GloTrackTag.isAvailable() ) {
          out << "\t"   << "Global(q*pT, eta, phi, dpT/pT) : " << GloTrackTag->charge() * GloTrackTag->pt() << " " << GloTrackTag->eta() << " " << GloTrackTag->phi() << " " << GloTrackTag->ptError() / GloTrackTag->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << GloTrackTag->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << GloTrackTag->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << GloTrackTag->hitPattern().numberOfValidMuonHits() << "\n";
        }
        if( InnTrackTag.isAvailable() ) {
          out << "\t"   << " Inner(q*pT, eta, phi, dpT/pT) : " << InnTrackTag->charge() * InnTrackTag->pt() << " " << InnTrackTag->eta() << " " << InnTrackTag->phi() << " " << InnTrackTag->ptError() / InnTrackTag->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << InnTrackTag->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << InnTrackTag->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << InnTrackTag->hitPattern().numberOfValidMuonHits() << "\n";
        }
        if( StaTrackTag.isAvailable() ) {
          out << "\t"   << "   STA(q*pT, eta, phi, dpT/pT) : " << StaTrackTag->charge() * StaTrackTag->pt() << " " << StaTrackTag->eta() << " " << StaTrackTag->phi() << " " << StaTrackTag->ptError() / StaTrackTag->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << StaTrackTag->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << StaTrackTag->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << StaTrackTag->hitPattern().numberOfValidMuonHits() << "\n";
        }
        if( TunTrackTag.isAvailable() ) {
          out << "\t"   << " TuneP(q*pT, eta, phi, dpT/pT) : " << TunTrackTag->charge() * TunTrackTag->pt() << " " << TunTrackTag->eta() << " " << TunTrackTag->phi() << " " << TunTrackTag->ptError() / TunTrackTag->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << TunTrackTag->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << TunTrackTag->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << TunTrackTag->hitPattern().numberOfValidMuonHits() << "\n";
        }
        if( DytTrackTag.isAvailable() ) {
          out << "\t"   << "   DYT(q*pT, eta, phi, dpT/pT) : " << DytTrackTag->charge() * DytTrackTag->pt() << " " << DytTrackTag->eta() << " " << DytTrackTag->phi() << " " << DytTrackTag->ptError() / DytTrackTag->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << DytTrackTag->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << DytTrackTag->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << DytTrackTag->hitPattern().numberOfValidMuonHits() << "\n";
        }
        if( PicTrackTag.isAvailable() ) {
          out << "\t"   << " Picky(q*pT, eta, phi, dpT/pT) : " << PicTrackTag->charge() * PicTrackTag->pt() << " " << PicTrackTag->eta() << " " << PicTrackTag->phi() << " " << PicTrackTag->ptError() / PicTrackTag->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << PicTrackTag->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << PicTrackTag->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << PicTrackTag->hitPattern().numberOfValidMuonHits() << "\n";
        }
        if( TpfTrackTag.isAvailable() ) {
          out << "\t"   << " TPFMS(q*pT, eta, phi, dpT/pT) : " << TpfTrackTag->charge() * TpfTrackTag->pt() << " " << TpfTrackTag->eta() << " " << TpfTrackTag->phi() << " " << TpfTrackTag->ptError() / TpfTrackTag->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << TpfTrackTag->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << TpfTrackTag->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << TpfTrackTag->hitPattern().numberOfValidMuonHits() << "\n";
        }

        //out << "lep1(pT, eta, phi) : " << lep1->pt() << " " << lep1->eta() << " " << lep1->phi() << "\n";
        out << "ProbeMu(pT, eta, phi) : " << ProbeMu->pt() << " " << ProbeMu->eta() << " " << ProbeMu->phi() << "\n";
        out << "# Mat St : " << ProbeMu->numberOfMatchedStations() << "\n";
        if( GloTrackProbe.isAvailable() ) {
          out << "\t"   << "Global(q*pT, eta, phi, dpT/pT) : " << GloTrackProbe->charge() * GloTrackProbe->pt() << " " << GloTrackProbe->eta() << " " << GloTrackProbe->phi() << " " << GloTrackProbe->ptError() / GloTrackProbe->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << GloTrackProbe->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << GloTrackProbe->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << GloTrackProbe->hitPattern().numberOfValidMuonHits() << "\n";
        }
        if( InnTrackProbe.isAvailable() ) {
          out << "\t"   << " Inner(q*pT, eta, phi, dpT/pT) : " << InnTrackProbe->charge() * InnTrackProbe->pt() << " " << InnTrackProbe->eta() << " " << InnTrackProbe->phi() << " " << InnTrackProbe->ptError() / InnTrackProbe->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << InnTrackProbe->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << InnTrackProbe->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << InnTrackProbe->hitPattern().numberOfValidMuonHits() << "\n";
        }
        if( StaTrackProbe.isAvailable() ) {
          out << "\t"   << "   STA(q*pT, eta, phi, dpT/pT) : " << StaTrackProbe->charge() * StaTrackProbe->pt() << " " << StaTrackProbe->eta() << " " << StaTrackProbe->phi() << " " << StaTrackProbe->ptError() / StaTrackProbe->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << StaTrackProbe->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << StaTrackProbe->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << StaTrackProbe->hitPattern().numberOfValidMuonHits() << "\n";
        }
        if( TunTrackProbe.isAvailable() ) {
          out << "\t"   << " TuneP(q*pT, eta, phi, dpT/pT) : " << TunTrackProbe->charge() * TunTrackProbe->pt() << " " << TunTrackProbe->eta() << " " << TunTrackProbe->phi() << " " << TunTrackProbe->ptError() / TunTrackProbe->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << TunTrackProbe->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << TunTrackProbe->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << TunTrackProbe->hitPattern().numberOfValidMuonHits() << "\n";

          if( TunTrackProbe->hitPattern().numberOfValidMuonHits() > 0 )
            nIsTunepPass += 1;
        }
        if( DytTrackProbe.isAvailable() ) {
          out << "\t"   << "   DYT(q*pT, eta, phi, dpT/pT) : " << DytTrackProbe->charge() * DytTrackProbe->pt() << " " << DytTrackProbe->eta() << " " << DytTrackProbe->phi() << " " << DytTrackProbe->ptError() / DytTrackProbe->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << DytTrackProbe->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << DytTrackProbe->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << DytTrackProbe->hitPattern().numberOfValidMuonHits() << "\n";

          if( DytTrackProbe->hitPattern().numberOfValidMuonHits() > 0 )
            nIsDYTPass += 1;
        }
        if( PicTrackProbe.isAvailable() ) {
          out << "\t"   << " Picky(q*pT, eta, phi, dpT/pT) : " << PicTrackProbe->charge() * PicTrackProbe->pt() << " " << PicTrackProbe->eta() << " " << PicTrackProbe->phi() << " " << PicTrackProbe->ptError() / PicTrackProbe->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << PicTrackProbe->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << PicTrackProbe->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << PicTrackProbe->hitPattern().numberOfValidMuonHits() << "\n";

          if( PicTrackProbe->hitPattern().numberOfValidMuonHits() > 0 )
            nIsPickyPass += 1;
        }
        if( TpfTrackProbe.isAvailable() ) {
          out << "\t"   << " TPFMS(q*pT, eta, phi, dpT/pT) : " << TpfTrackProbe->charge() * TpfTrackProbe->pt() << " " << TpfTrackProbe->eta() << " " << TpfTrackProbe->phi() << " " << TpfTrackProbe->ptError() / TpfTrackProbe->pt() << "\n";
          out << "\t\t" << "# Pix hits : " << TpfTrackProbe->hitPattern().numberOfValidPixelHits() << "\n";
          out << "\t\t" << "# Trk Lays : " << TpfTrackProbe->hitPattern().trackerLayersWithMeasurement() << "\n";
          out << "\t\t" << " # Mu hits : " << TpfTrackProbe->hitPattern().numberOfValidMuonHits() << "\n";

          if( TpfTrackProbe->hitPattern().numberOfValidMuonHits() > 0 )
            nIsTPFMSPass += 1;
        }
        out << "\n";

        reco::TrackBase::Point bs(beamSpot->x0(), beamSpot->y0(), beamSpot->z0());
        osprintf(out, "beamspot(x,y,z) : %f %f %f\n", bs.x(), bs.y(), bs.z());
        out << "offline PV(x,y,z) : " << vertex->x() << " " << vertex->y() << " " << vertex->z() << "  nVTX = " << nVtx << "\n";
        if (dil.hasUserFloat("vertexX"))
          out << "Common dimuon vertex(x, y, z) : "
              << dil.userFloat("vertexX") << ", " << dil.userFloat("vertexY") << ", " << dil.userFloat("vertexZ")
              << " chi2/dof: " << dil.userFloat("vertex_chi2") << "\n";
        if (dil.hasUserFloat("vertexM"))
          out << "Mass computed with the common-vertex constraint : "
              << dil.userFloat("vertexM") << " /- " << dil.userFloat("vertexMError") << "\n";

        nProbeT += 1;
        if( fabs(ProbeMu->eta()) <= 0.9 )
          nProbeB += 1;
        else if( fabs(ProbeMu->eta()) > 0.9 && fabs(ProbeMu->eta()) <= 1.2 )
          nProbeO += 1;
        else if( fabs(ProbeMu->eta()) > 1.2 && fabs(ProbeMu->eta()) <= 2.1 )
          nProbeE += 1;
        else if( fabs(ProbeMu->eta()) > 2.1 && fabs(ProbeMu->eta()) <= 2.4 )
          nProbeF += 1;

        if( ProbeMu->numberOfMatchedStations() > 1 )
          nIsNStationMoreThan1 += 1;

        filter = true;
      }
      out << "\n";
      std::cout << out.str();
    }
  }

  return filter;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DileptonCategorizer::beginJob()
{
  nProbeT = 0;
  nProbeB = 0;
  nProbeO = 0;
  nProbeE = 0;
  nProbeF = 0;

  nIsTunepPass = 0;
  nIsDYTPass = 0;
  nIsPickyPass = 0;
  nIsTPFMSPass = 0;

  nIsNStationMoreThan1 = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DileptonCategorizer::endJob()
{
  std::ostringstream out;

  out << "Event loop ended\n";

  out << "nProbeT : " << nProbeT << "\n";
  out << "nProbeB : " << nProbeB << "\n";
  out << "nProbeO : " << nProbeO << "\n";
  out << "nProbeE : " << nProbeE << "\n";
  out << "nProbeF : " << nProbeF << "\n";

  out << "nIsNStationMoreThan1 : " << nIsNStationMoreThan1 << "\n";

  out << "nIsTunepPass : " << nIsTunepPass << "\n";
  out << "  nIsDYTPass : " << nIsDYTPass << "\n";
  out << "nIsPickyPass : " << nIsPickyPass << "\n";
  out << "nIsTPFMSPass : " << nIsTPFMSPass << "\n";

  out << "\n\n";
  std::cout << out.str();
}

//define this as a plug-in
DEFINE_FWK_MODULE(DileptonCategorizer);
