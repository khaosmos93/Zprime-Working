#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/MuonReco/interface/MuonShower.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "DataFormats/DTRecHit/interface/DTSLRecSegment2D.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentUpdator.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentCleaner.h"
#include "RecoLocalMuon/DTSegment/src/DTHitPairForFit.h"

#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCRangeMapAccessor.h"

#include <vector>

bool approxEqual(LocalPoint a, LocalPoint b, const double tol=1.E-3)
{
  return ( (fabs(a.x()-b.x()) < tol) && (fabs(a.y()-b.y()) < tol) && (fabs(a.z()-b.z()) < tol) );
}

class ZprimeTnPPairSelector_FromAOD : public edm::EDProducer {
public:
  explicit ZprimeTnPPairSelector_FromAOD(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // Helper stuff.
  struct reverse_mass_sort {
    bool operator()(const pat::CompositeCandidate& lhs, const pat::CompositeCandidate& rhs) {
      return lhs.mass() > rhs.mass();
    }
  };

  struct lepton_pt_sort {
    bool operator()(const pat::CompositeCandidate& lhs, const pat::CompositeCandidate& rhs) {
      // Sort by pT: a pair with two highest-pT muons wins.
      const reco::CandidateBaseRef& lhs_mu0 = dileptonDaughter(lhs, 0);
      const reco::CandidateBaseRef& lhs_mu1 = dileptonDaughter(lhs, 1);
      const reco::CandidateBaseRef& rhs_mu0 = dileptonDaughter(rhs, 0);
      const reco::CandidateBaseRef& rhs_mu1 = dileptonDaughter(rhs, 1);
      // Get the unique ids of the dilepton daughters, as well as their pT's.
      int lhs_id0 = lhs_mu0.key();
      int lhs_id1 = lhs_mu1.key();
      int rhs_id0 = rhs_mu0.key();
      int rhs_id1 = rhs_mu1.key();
      double lhs_pt0 = lhs_mu0->pt();
      double lhs_pt1 = lhs_mu1->pt();
      double rhs_pt0 = rhs_mu0->pt();
      double rhs_pt1 = rhs_mu1->pt();
      //std::ostringstream out;
      //out << " lhs0: " << lhs_id0 << " " << lhs_pt0 << "\n";
      //out << " lhs1: " << lhs_id1 << " " << lhs_pt1 << "\n";
      //out << " rhs0: " << rhs_id0 << " " << rhs_pt0 << "\n";
      //out << " rhs1: " << rhs_id1 << " " << rhs_pt1 << "\n";
      //edm::LogDebug("Sort dileptons") << out.str();
      // Sort muons in each dimuon in decreasing order of pT
      if (lhs_pt0 < lhs_pt1) {
        double tmp_pt = lhs_pt0;
        lhs_pt0 = lhs_pt1;
        lhs_pt1 = tmp_pt;
        int tmp_id = lhs_id0;
        lhs_id0 = lhs_id1;
        lhs_id1 = tmp_id;
      }
      if (rhs_pt0 < rhs_pt1) {
        double tmp_pt = rhs_pt0;
        rhs_pt0 = rhs_pt1;
        rhs_pt1 = tmp_pt;
        int tmp_id = rhs_id0;
        rhs_id0 = rhs_id1;
        rhs_id1 = tmp_id;
      }
      // In sorting by pT, if there are more than four muons, only the
      // highest-ranking and lowest-ranking dimuons are well defined.
      // That's OK because we really care only about the
      // highest-ranking dimuon.
      if (lhs_id0 == rhs_id0 && lhs_id1 == rhs_id1) {
        edm::LogWarning("Sort dileptons") << "+++ Two identical dimuons found? +++";
        return false;
      }
      else if (lhs_id0 == rhs_id0 && lhs_pt1 > rhs_pt1)
        return true;
      else if (lhs_id1 == rhs_id1 && lhs_pt0 > rhs_pt0)
        return true;
      else if (lhs_id0 != rhs_id0 && lhs_id1 != rhs_id1 && lhs_pt0 > rhs_pt0 && lhs_pt1 > rhs_pt1)
        return true;
      else
        return false;
    }
  };

  void remove_overlap(pat::CompositeCandidateCollection&) const;
  std::vector<reco::TransientTrack> get_transient_tracks(const pat::CompositeCandidate&) const;

  // Evaluate cuts. Return values are pair<cut decision, variable or
  // variables to embed>. The cut decision should be made whether
  // we're actually going to drop the candidate because of this
  // decision or not (controlled by the cut_on_* variables); that will
  // be handled in the loop in produce.

  // FIXME : add muon around probe veto!!!
  bool                               TagAndProbeSelector(const pat::CompositeCandidate&, float, float, float, float) const;

  std::pair<bool, float>             back_to_back_cos_angle(const pat::CompositeCandidate&) const;
  std::pair<bool, CachingVertex<5> > vertex_constrained_fit(const pat::CompositeCandidate&) const;
  float                              dpt_over_pt(const reco::CandidateBaseRef&) const;

  std::pair<bool, float>             pt_ratio(const pat::CompositeCandidate&) const;
  std::pair<bool, float>             dil_deltaR(const pat::CompositeCandidate&) const;

  float                              veto_others_dphi(edm::Event&, const reco::CandidateBaseRef&);

  std::vector<int> countDTsegs(edm::Event&, reco::MuonRef);
  std::vector<int> countCSCsegs(edm::Event&, reco::MuonRef);
  void testCSChits(edm::Event&, reco::TrackRef);

  std::vector<int> nSegments(edm::Event&, reco::MuonRef);

  std::pair<bool, reco::MuonRef> getMuonRef(edm::Event&, const reco::CandidateBaseRef&);
  std::vector<int> nHits(edm::Event&, reco::MuonRef);

  // If the variable to embed in the methods above is a simple int or
  // float or is going to be embedded wholesale with the generic
  // userData mechanism, we'll do those explicitly in the loop in
  // produce. Otherwise if it's more complicated, the methods here
  // take care of it.
  void embed_vertex_constrained_fit(pat::CompositeCandidate&, const CachingVertex<5>& vtx) const;

  const edm::InputTag src;
  const edm::InputTag reco_muon_src;
  const edm::InputTag muonshower_src;
  const edm::InputTag dtseg_src;
  const edm::InputTag cscseg_src;

  edm::InputTag vertex_src;
  const reco::Vertex*   PV;
  StringCutObjectSelector<pat::CompositeCandidate> selector;
  StringCutObjectSelector<pat::Muon> tag_selector;
  StringCutObjectSelector<pat::Muon> probe_selector;

  const unsigned max_candidates;
  const bool sort_by_pt;
  const bool do_remove_overlap;

  const bool cut_on_back_to_back_cos_angle;
  const double back_to_back_cos_angle_min;

  const bool cut_on_vertex_chi2;
  const double vertex_chi2_max;

  const bool cut_on_tag_dpt_over_pt;
  const double tag_dpt_over_pt_max;
  const bool cut_on_tag_dz;
  const double tag_dz_max;

  const bool cut_on_probe_dpt_over_pt;
  const double probe_dpt_over_pt_max;

  const bool cut_on_pt_ratio;
  const double pt_ratio_max;

  const bool cut_on_dil_deltaR;
  const double dil_deltaR_min;

  const bool cut_on_veto_others_dphi;
  const double veto_others_dphi_min;

  const bool samePV;

  edm::ESHandle<TransientTrackBuilder> ttkb;

  const bool ShutUp;
};

ZprimeTnPPairSelector_FromAOD::ZprimeTnPPairSelector_FromAOD(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),
    reco_muon_src(cfg.getParameter<edm::InputTag>("reco_muon_src")),
    muonshower_src(cfg.getParameter<edm::InputTag>("muonshower_src")),
    dtseg_src(cfg.getParameter<edm::InputTag>("dtseg_src")),
    cscseg_src(cfg.getParameter<edm::InputTag>("cscseg_src")),

    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    PV(0),
    selector(cfg.getParameter<std::string>("cut")),

    ///
    tag_selector(cfg.getParameter<std::string>("tag_cut")),
    probe_selector(cfg.getParameter<std::string>("probe_cut")),

    max_candidates(cfg.getParameter<unsigned>("max_candidates")),
    sort_by_pt(cfg.getParameter<bool>("sort_by_pt")),
    do_remove_overlap(cfg.getParameter<bool>("do_remove_overlap")),
    cut_on_back_to_back_cos_angle(cfg.existsAs<double>("back_to_back_cos_angle_min")),
    back_to_back_cos_angle_min(cut_on_back_to_back_cos_angle ? cfg.getParameter<double>("back_to_back_cos_angle_min") : -2),
    cut_on_vertex_chi2(cfg.existsAs<double>("vertex_chi2_max")),
    vertex_chi2_max(cut_on_vertex_chi2 ? cfg.getParameter<double>("vertex_chi2_max") : 1e99),
    cut_on_tag_dpt_over_pt(cfg.existsAs<double>("tag_dpt_over_pt_max")),
    tag_dpt_over_pt_max(cut_on_tag_dpt_over_pt ? cfg.getParameter<double>("tag_dpt_over_pt_max") : 1e99),
    cut_on_tag_dz(cfg.existsAs<double>("tag_dz_max")),
    tag_dz_max(cut_on_tag_dz ? cfg.getParameter<double>("tag_dz_max") : 1e99),
    cut_on_probe_dpt_over_pt(cfg.existsAs<double>("probe_dpt_over_pt_max")),
    probe_dpt_over_pt_max(cut_on_probe_dpt_over_pt ? cfg.getParameter<double>("probe_dpt_over_pt_max") : 1e99),


    cut_on_pt_ratio(cfg.existsAs<double>("pt_ratio_max")),
    pt_ratio_max(cut_on_pt_ratio ? cfg.getParameter<double>("pt_ratio_max") : 1e99),
    cut_on_dil_deltaR(cfg.existsAs<double>("dil_deltaR_min")),
    dil_deltaR_min(cut_on_dil_deltaR ? cfg.getParameter<double>("dil_deltaR_min") : -2),

    cut_on_veto_others_dphi(cfg.existsAs<double>("veto_others_dphi_min")),
    veto_others_dphi_min(cut_on_veto_others_dphi ? cfg.getParameter<double>("veto_others_dphi_min") : -2.),

    samePV(cfg.getParameter<bool>("samePV")),

    ShutUp(cfg.getParameter<bool>("ShutUp"))
{
 consumes<pat::CompositeCandidateCollection>(src);
 consumes<std::vector< reco::Muon >>(reco_muon_src);
 consumes<edm::ValueMap<reco::MuonShower>>(muonshower_src);
 consumes<DTRecSegment4DCollection>(dtseg_src);
 consumes<CSCSegmentCollection>(cscseg_src);

 consumes<reco::VertexCollection>(vertex_src);
 produces<pat::CompositeCandidateCollection>();

}

void ZprimeTnPPairSelector_FromAOD::remove_overlap(pat::CompositeCandidateCollection& cands) const {
  // For the list of CompositeCandidates, find any that share leptons
  // and remove one of them. The sort order of the input is used to
  // determine which of the pair is to be removed: we keep the first
  // one.

  // Don't bother doing anything if there's just one candidate.
  if (cands.size() < 2) return;

  pat::CompositeCandidateCollection::iterator p, q;
  for (p = cands.begin(); p != cands.end() - 1; ) {
    for (q = p + 1; q != cands.end(); ++q) {
      // Check to see if any of the leptons in p is in q also. If so,
      // remove q (e.g. the one with lower invariant mass since we
      // have sorted the vector already), reset pointers and restart.

      // To do this we need the unique ids of the daughters, i.e. the
      // refs into the original lepton collections.
      typedef std::vector<reco::CandidateBaseRef> refs;
      refs prefs, qrefs;
      for (size_t i = 0; i < p->numberOfDaughters(); ++i)
        prefs.push_back(p->daughter(i)->masterClone());
      for (size_t i = 0; i < q->numberOfDaughters(); ++i)
        qrefs.push_back(q->daughter(i)->masterClone());

      // Compare every pair of (pref, qref) to check for any lepton
      // being shared.
      bool any_shared = false;
      refs::const_iterator pr = prefs.begin(), pre = prefs.end(),
        qr = qrefs.begin(), qre = qrefs.end();
      for ( ; pr != pre && !any_shared; ++pr)
        for ( ; qr != qre && !any_shared; ++qr)
          if (pr == qr)
            any_shared = true;

      if (any_shared) {
        cands.erase(q);
        p = cands.begin();
      }
      else
        ++p;
    }
  }
}

std::vector<reco::TransientTrack> ZprimeTnPPairSelector_FromAOD::get_transient_tracks(const pat::CompositeCandidate& dil) const {
  // Get TransientTracks (for use in e.g. the vertex fit) for each of
  // the muon tracks, using e.g. the cocktail momentum.

  std::vector<reco::TransientTrack> ttv;
  const size_t n = dil.numberOfDaughters();
  for (size_t i = 0; i < n; ++i) {
    const pat::Muon* mu = toConcretePtr<pat::Muon>(dileptonDaughter(dil, i));
    assert(mu);
    const reco::TrackRef& tk = patmuon::getPickedTrack(*mu);
    ttv.push_back(ttkb->build(tk));
  }

  return ttv;
}

std::pair<bool, float> ZprimeTnPPairSelector_FromAOD::back_to_back_cos_angle(const pat::CompositeCandidate& dil) const {
  // Back-to-back cut to kill cosmics.
  assert(dil.numberOfDaughters() == 2);
  const float cos_angle = dil.daughter(0)->momentum().Dot(dil.daughter(1)->momentum()) / dil.daughter(0)->p() / dil.daughter(1)->p();
  return std::make_pair(cos_angle >= back_to_back_cos_angle_min, cos_angle);
}

std::pair<bool, CachingVertex<5> > ZprimeTnPPairSelector_FromAOD::vertex_constrained_fit(const pat::CompositeCandidate& dil) const {
  // Loose common vertex chi2 cut.
  assert(dil.numberOfDaughters() == 2);
  if (abs(dil.daughter(0)->pdgId()) != 13 || abs(dil.daughter(1)->pdgId()) != 13)
    return std::make_pair(true, CachingVertex<5>()); // pass objects we don't know how to cut on, i.e. e-mu dileptons

  KalmanVertexFitter kvf(true);
  CachingVertex<5> v = kvf.vertex(get_transient_tracks(dil));

  return std::make_pair(v.isValid() && v.totalChiSquared()/v.degreesOfFreedom() <= vertex_chi2_max, v);
}

void ZprimeTnPPairSelector_FromAOD::embed_vertex_constrained_fit(pat::CompositeCandidate& dil, const CachingVertex<5>& vtx) const {
  if (!vtx.isValid()) {
    dil.addUserFloat("vertex_chi2", 1e8);
    return;
  }

  dil.addUserFloat("vertex_chi2", vtx.totalChiSquared()/vtx.degreesOfFreedom());
  dil.addUserFloat("vertex_ndof", vtx.degreesOfFreedom());
  dil.addUserFloat("vertexX", vtx.position().x());
  dil.addUserFloat("vertexY", vtx.position().y());
  dil.addUserFloat("vertexZ", vtx.position().z());
  dil.addUserFloat("vertexXError", sqrt(vtx.error().cxx()));
  dil.addUserFloat("vertexYError", sqrt(vtx.error().cyy()));
  dil.addUserFloat("vertexZError", sqrt(vtx.error().czz()));

  InvariantMassFromVertex imfv;
  static const double muon_mass = 0.1056583;
  InvariantMassFromVertex::LorentzVector p4 = imfv.p4(vtx, muon_mass);
  Measurement1D mass = imfv.invariantMass(vtx, muon_mass);

  dil.addUserFloat("vertexPX", p4.X());
  dil.addUserFloat("vertexPY", p4.Y());
  dil.addUserFloat("vertexPZ", p4.Z());

  dil.addUserFloat("vertexM",      mass.value());
  dil.addUserFloat("vertexMError", mass.error());
}

float ZprimeTnPPairSelector_FromAOD::dpt_over_pt(const reco::CandidateBaseRef& lep) const {
  // Cut on sigma(pT)/pT to reject grossly mismeasured tracks.
  double dpt_over_pt = 999;
  if (lep.isNonnull()) {
    const pat::Muon* mu = toConcretePtr<pat::Muon>(lep);
    if (mu) {
      const reco::Track* tk = patmuon::getPickedTrack(*mu).get();
      if (tk) {
        dpt_over_pt = ptError(tk)/tk->pt();
      }
    }
  }
  return dpt_over_pt;
}

///
std::pair<bool, float> ZprimeTnPPairSelector_FromAOD::pt_ratio(const pat::CompositeCandidate& dil) const {
  double pt0 = dil.daughter(0)->pt();
  double pt1 = dil.daughter(1)->pt();

  double the_pt_ratio = pt0 > pt1 ? (pt0/pt1) : (pt1/pt0);

  //if(!ShutUp)  std::cout << "ZprimeTnPPairSelector_FromAOD::pt_ratio[" << pt_ratio_max << "] : pt_ratio = " << the_pt_ratio << ", " << (the_pt_ratio < pt_ratio_max) << std::endl;

  return std::make_pair( the_pt_ratio < pt_ratio_max , the_pt_ratio );
}

///
std::pair<bool, float> ZprimeTnPPairSelector_FromAOD::dil_deltaR(const pat::CompositeCandidate& dil) const {

  double the_deltaR = reco::deltaR( *dil.daughter(0), *dil.daughter(1) );

  //if(!ShutUp)  std::cout << "ZprimeTnPPairSelector_FromAOD::dil_deltaR[" << dil_deltaR_min << "] : deltaR = " << the_deltaR << ", " << (the_deltaR > dil_deltaR_min) << std::endl;

  return std::make_pair( the_deltaR > dil_deltaR_min , the_deltaR );
}

///
float ZprimeTnPPairSelector_FromAOD::veto_others_dphi(edm::Event& event, const reco::CandidateBaseRef& lep) {

  edm::Handle< std::vector< reco::Muon > > recoMuons;
  event.getByLabel(reco_muon_src, recoMuons);

  float the_dphi = 999.;

  if (lep.isNonnull()) {
    const pat::Muon* mu = toConcretePtr<pat::Muon>(lep);

    for (std::vector<reco::Muon>::const_iterator imu = recoMuons->begin(); imu != recoMuons->end(); imu++) {
      if( imu->pt() < 10 || (imu->standAloneMuon()).isNull() )  continue;
      if( imu->standAloneMuon()->pt()  == mu->standAloneMuon()->pt() &&
          imu->standAloneMuon()->eta() == mu->standAloneMuon()->eta() &&
          imu->standAloneMuon()->phi() == mu->standAloneMuon()->phi() )
        continue;

      float temp_dphi = fabs( reco::deltaPhi(imu->phi(), mu->phi()) );
      if( temp_dphi < the_dphi )
        the_dphi = temp_dphi;
    }
  }

  return the_dphi;
}

///
bool ZprimeTnPPairSelector_FromAOD::TagAndProbeSelector(const pat::CompositeCandidate& dil,
                                                        float lep0_dpt_over_pt,      float lep1_dpt_over_pt,
                                                        float lep0_veto_others_dphi, float lep1_veto_others_dphi) const {
  bool isTag0Probe1 = false;
  bool isTag1Probe0 = false;

  const reco::CandidateBaseRef& lep0 = dileptonDaughter(dil, 0);
  const reco::CandidateBaseRef& lep1 = dileptonDaughter(dil, 1);
  if (lep0.isNonnull() && lep1.isNonnull()) {
    const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
    const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);

    if (mu0 && mu1) {

      isTag0Probe1 = ( tag_selector(*mu0)   && (lep0_dpt_over_pt<tag_dpt_over_pt_max) && (fabs(mu0->muonBestTrack()->dz( PV->position() )) < tag_dz_max) &&
                       probe_selector(*mu1) && (lep1_dpt_over_pt<probe_dpt_over_pt_max) );

      isTag1Probe0 = ( tag_selector(*mu1)   && (lep1_dpt_over_pt<tag_dpt_over_pt_max) && (fabs(mu1->muonBestTrack()->dz( PV->position() )) < tag_dz_max) &&
                       probe_selector(*mu0) && (lep0_dpt_over_pt<probe_dpt_over_pt_max) );

      if(cut_on_veto_others_dphi) {
        isTag0Probe1 = ( isTag0Probe1 && (lep1_veto_others_dphi>veto_others_dphi_min) );
        isTag1Probe0 = ( isTag1Probe0 && (lep0_veto_others_dphi>veto_others_dphi_min) );
      }

      if(!ShutUp && !cut_on_veto_others_dphi)  std::cout << "ZprimeTnPPairSelector_FromAOD::TagAndProbeSelector[" << (isTag0Probe1 || isTag1Probe0) << "] : isTag0=" 
        << ( tag_selector(*mu0)   && (lep0_dpt_over_pt<tag_dpt_over_pt_max) && (fabs(mu0->muonBestTrack()->dz( PV->position() )) < tag_dz_max) ) << " " << "isTag1=" 
        << ( tag_selector(*mu1)   && (lep1_dpt_over_pt<tag_dpt_over_pt_max) && (fabs(mu1->muonBestTrack()->dz( PV->position() )) < tag_dz_max) ) << " " << "isProbe0=" 
        << ( probe_selector(*mu0) && (lep0_dpt_over_pt<probe_dpt_over_pt_max) ) << " " << "isProbe1=" 
        << ( probe_selector(*mu1) && (lep1_dpt_over_pt<probe_dpt_over_pt_max) ) << std::endl;

      if(!ShutUp && cut_on_veto_others_dphi)  std::cout << "ZprimeTnPPairSelector_FromAOD::TagAndProbeSelector[" << (isTag0Probe1 || isTag1Probe0) << "] : isTag0=" 
        << ( tag_selector(*mu0)   && (lep0_dpt_over_pt<tag_dpt_over_pt_max) && (fabs(mu0->muonBestTrack()->dz( PV->position() )) < tag_dz_max) ) << " " << "isTag1=" 
        << ( tag_selector(*mu1)   && (lep1_dpt_over_pt<tag_dpt_over_pt_max) && (fabs(mu1->muonBestTrack()->dz( PV->position() )) < tag_dz_max) ) << " " << "isProbe0=" 
        << ( probe_selector(*mu0) && (lep0_dpt_over_pt<probe_dpt_over_pt_max) && (lep0_veto_others_dphi>veto_others_dphi_min) ) << " " << "isProbe1=" 
        << ( probe_selector(*mu1) && (lep1_dpt_over_pt<probe_dpt_over_pt_max) && (lep1_veto_others_dphi>veto_others_dphi_min) ) << std::endl;

      double DeltaX = fabs(mu0->innerTrack()->referencePoint().x() - mu1->innerTrack()->referencePoint().x());
      double DeltaY = fabs(mu0->innerTrack()->referencePoint().y() - mu1->innerTrack()->referencePoint().y());
      double DeltaZ = fabs(mu0->innerTrack()->referencePoint().z() - mu1->innerTrack()->referencePoint().z());
      if( samePV && !( (sqrt( DeltaX*DeltaX + DeltaY*DeltaY ) < 0.02) && (DeltaZ < 0.05) ) ) {
        if(!ShutUp) std::cout << "Not from same PV(|dxy|<0.02 && |dz|<0.05)" << " : " 
                              << sqrt( DeltaX*DeltaX + DeltaY*DeltaY ) << ", " << DeltaZ << std::endl;
        return false;
      }

    }
  }

  return ( isTag0Probe1 || isTag1Probe0 );
}


std::vector<int> ZprimeTnPPairSelector_FromAOD::countDTsegs(edm::Event& event, reco::MuonRef muon) {
  double DTCut = 30.;

  std::vector<int> stations={0,0,0,0};
  std::vector<int> removed={0,0,0,0};

  edm::Handle<DTRecSegment4DCollection> dtRecHits;
  event.getByLabel(dtseg_src, dtRecHits);

  if (!ShutUp) std::cout << std::endl << " *** DT Segment search" << std::endl;
  for (const auto &ch : muon->matches()) {
    if( ch.detector() != MuonSubdetId::DT )  continue;
    DTChamberId DTid( ch.id.rawId() );
    if (!ShutUp) std::cout << "   DT chamber in station " << ch.station() << "  DTChamberId:" << DTid << "  local position: (" << ch.x << ", " << ch.y << ", 0)" << std::endl;

    int nsegs_temp = 0;
    std::vector<float> nsegs_x_temp, nsegs_y_temp;

    for (auto seg = dtRecHits->begin(); seg!=dtRecHits->end(); ++seg) {
      DTChamberId myChamber((*seg).geographicalId().rawId());
      if (!(DTid==myChamber))  continue;
      LocalPoint posLocalSeg = seg->localPosition();
      if (!ShutUp) std::cout << "     Found segment in station " << ch.station() << "  DTChamberId:" << myChamber << "  local position: " << posLocalSeg << std::endl;

      /*if( ( (posLocalSeg.x()==0 && posLocalSeg.y()!=0) && (fabs(posLocalSeg.y()-ch.y)<DTCut) ) ||
          ( (posLocalSeg.x()!=0 && posLocalSeg.y()==0) && (fabs(posLocalSeg.x()-ch.x)<DTCut) ) ||
          ( (posLocalSeg.x()!=0 && posLocalSeg.y()!=0) && (sqrt( (posLocalSeg.x()-ch.x)*(posLocalSeg.x()-ch.x) + (posLocalSeg.y()-ch.y)*(posLocalSeg.y()-ch.y) )<DTCut) )
        ) {
        nsegs_temp++;
      }*/

      if( ( posLocalSeg.x()!=0 && ch.x!=0 ) && (fabs(posLocalSeg.x()-ch.x)<DTCut) ) {
        bool found = false;
        for( auto prev_x : nsegs_x_temp) {
          if( fabs(prev_x-posLocalSeg.x()) < 0.1 ) {
            found = true;
            break;
          }
        }
        if( !found )  nsegs_x_temp.push_back(posLocalSeg.x());
      }

      if( ( posLocalSeg.y()!=0 && ch.y!=0 ) && (fabs(posLocalSeg.y()-ch.y)<DTCut) ) {
        bool found = false;
        for( auto prev_y : nsegs_y_temp) {
          if( fabs(prev_y-posLocalSeg.y()) < 0.1 ) {
            found = true;
            break;
          }
        }
        if( !found )  nsegs_y_temp.push_back(posLocalSeg.y());
      }

    }

    nsegs_temp = (int)std::max(nsegs_x_temp.size(), nsegs_y_temp.size());

    //--- subtract best matched segment from given muon
    bool isBestMatched = false;
    for(std::vector<reco::MuonSegmentMatch>::const_iterator matseg = ch.segmentMatches.begin(); matseg != ch.segmentMatches.end(); matseg++) {
      if( matseg->isMask(reco::MuonSegmentMatch::BestInChamberByDR) ) {
        if (!ShutUp) std::cout << "     Found BestInChamberByDR in station " << ch.station() << "  DTChamberId:" << DTid << "  local position: (" << matseg->x << ", " << matseg->y << ", 0)" << std::endl;
        removed[ch.station()-1]++;
        isBestMatched = true;
        break;
      }
    }
    if (!ShutUp) std::cout << std::endl;

    if(isBestMatched) nsegs_temp = nsegs_temp-1;
    if(nsegs_temp>0)  stations[ch.station()-1] += nsegs_temp;

  }

  if (!ShutUp) {
    std::cout << " DT Shower pattern: ";
    int DTSum = 0;
    for (int i=0;i<4;i++) {
      std::cout << stations[i] << " ";
      DTSum += stations[i];
    }
    std::cout << std::endl;
    std::cout << " DTSum = " << DTSum << std::endl;

    std::cout << " DT removed pattern: ";
    for (int i=0;i<4;i++) {
      std::cout << removed[i] << " ";
    }
    std::cout << std::endl;
  }

  return stations;
}

std::vector<int> ZprimeTnPPairSelector_FromAOD::countCSCsegs(edm::Event& event, reco::MuonRef muon) {
  double CSCCut = 30.;

  std::vector<int> stations={0,0,0,0};
  std::vector<int> removed={0,0,0,0};

  edm::Handle<CSCSegmentCollection> cscRecHits;
  event.getByLabel(cscseg_src, cscRecHits);

  if (!ShutUp) std::cout << std::endl << " *** CSC Segment search" << std::endl;
  for (const auto &ch : muon->matches()) {
    if( ch.detector() != MuonSubdetId::CSC )  continue;
    CSCDetId CSCid( ch.id.rawId() );
    if (!ShutUp) std::cout << "   CSC chamber in station " << ch.station() << "  CSCDetId:" << CSCid << "  local position: (" << ch.x << ", " << ch.y << ", 0)" << std::endl;

    int nsegs_temp = 0;

    for (auto seg = cscRecHits->begin(); seg!=cscRecHits->end(); ++seg) {
      CSCDetId myChamber((*seg).geographicalId().rawId());
      if (!(CSCid==myChamber))  continue;
      LocalPoint posLocalSeg = seg->localPosition();
      if (!ShutUp) std::cout << "     Found segment in station " << ch.station() << "  CSCDetId:" << myChamber << "  local position: " << posLocalSeg << std::endl;

      if( (posLocalSeg.x()!=0 && posLocalSeg.y()!=0) && (sqrt( (posLocalSeg.x()-ch.x)*(posLocalSeg.x()-ch.x) + (posLocalSeg.y()-ch.y)*(posLocalSeg.y()-ch.y) )<CSCCut) )  { 
        nsegs_temp++;
      }
    }

    //--- subtract best matched segment from given muon
    bool isBestMatched = false;
    for(std::vector<reco::MuonSegmentMatch>::const_iterator matseg = ch.segmentMatches.begin(); matseg != ch.segmentMatches.end(); matseg++) {
      if( matseg->isMask(reco::MuonSegmentMatch::BestInChamberByDR) ) {
        if (!ShutUp) std::cout << "     Found BestInChamberByDR in station " << ch.station() << "  CSCDetId:" << CSCid << "  local position: (" << matseg->x << ", " << matseg->y << ", 0)" << std::endl;
        removed[ch.station()-1]++;
        isBestMatched = true;
        break;
      }
    }
    if (!ShutUp) std::cout << std::endl;

    if(isBestMatched) nsegs_temp = nsegs_temp-1;
    if(nsegs_temp>0)  stations[ch.station()-1] += nsegs_temp;

    /*for(std::vector<reco::MuonSegmentMatch>::const_iterator seg = ch.segmentMatches.begin(); seg != ch.segmentMatches.end(); seg++) {
      if( seg->isMask(reco::MuonSegmentMatch::BestInChamberByDR) ) {
        if (!ShutUp) std::cout << "     Found BestInChamberByDR in station " << ch.station() << "  local position: (" << seg->x << ", " << seg->y << ", 0)" << std::endl;
        removed[ch.station()-1]++;
      }
      else {
        if (!ShutUp) std::cout << "     Found segment in station " << ch.station() << "  local position: (" << seg->x << ", " << seg->y << ", 0)" << std::endl;
        if( (seg->x!=0 && seg->y!=0) && (sqrt( (seg->x-ch.x)*(seg->x-ch.x) + (seg->y-ch.y)*(seg->y-ch.y) )<CSCCut) )  { stations[ch.station()-1]++; }
      }
    }*/
  }

  if (!ShutUp) {
    std::cout << " CSC Shower pattern: ";
    int CSCSum = 0;
    for (int i=0;i<4;i++) {
      std::cout << stations[i] << " ";
      CSCSum += stations[i];
    }
    std::cout << std::endl;
    std::cout << " CSCSum = " << CSCSum << std::endl;

    std::cout << " CSC removed pattern: ";
    for (int i=0;i<4;i++) {
      std::cout << removed[i] << " ";
    }
    std::cout << std::endl;
  }

  return stations;
}

void ZprimeTnPPairSelector_FromAOD::testCSChits(edm::Event& event, reco::TrackRef muon) {

  int endcap = -1;
  std::vector<int> stations={0,0,0,0};

  edm::Handle<CSCSegmentCollection> cscRecHits;
  event.getByLabel(cscseg_src, cscRecHits);  

  std::cout << std::endl << " *** CSC muon recHit search" << std::endl;

  // Loop over muon recHits
  for(trackingRecHit_iterator muonHit = muon->recHitsBegin(); muonHit != muon->recHitsEnd(); ++muonHit) {
    if ( (*muonHit)->geographicalId().det() != DetId::Muon ) continue; 
    if ( (*muonHit)->geographicalId().subdetId() != MuonSubdetId::CSC ) continue;

    CSCDetId cscDetIdHitT((*muonHit)->geographicalId());
    LocalPoint posLocalMuon = (*muonHit)->localPosition();
    std::cout << "Found in " << std::endl;
    std::cout << "\t  Endcap " << cscDetIdHitT.endcap() << std::endl;
    std::cout << "\t Station " << cscDetIdHitT.station() << std::endl;
    std::cout << "\t    Ring " << cscDetIdHitT.ring() << std::endl;
    std::cout << "\t Chamber " << cscDetIdHitT.chamber() << std::endl;
    std::cout << "\t LocalPo " << posLocalMuon << std::endl;
    std::cout << std::endl;

    endcap = cscDetIdHitT.endcap();
    stations[cscDetIdHitT.station()-1]++;
  }
  std::cout << " Muon recHit in stations : ";
  for (int i=0;i<4;i++) {
    std::cout << stations[i] << " ";
  }
  std::cout << std::endl;

  for(int st=0; st<4; ++st) {
    if(stations[st] == 0) {
      std::cout << " No muon rechit in station " << (st+1) << std::endl;

      for (auto rechit = cscRecHits->begin(); rechit!=cscRecHits->end();++rechit) {
        CSCDetId myChamber((*rechit).geographicalId().rawId());
        if( ( (int)(myChamber.station()-1) == st ) && ( endcap == (int)(myChamber.endcap()) ) ) {
          LocalPoint posLocalHit = rechit->localPosition();
          std::cout << "\t  Endcap " << myChamber.endcap() << std::endl;
          std::cout << "\t Station " << myChamber.station() << std::endl;
          std::cout << "\t    Ring " << myChamber.ring() << std::endl;
          std::cout << "\t Chamber " << myChamber.chamber() << std::endl;
          std::cout << "\t LocalPo " << posLocalHit << std::endl;
        }
      }

    }
  }
}

std::vector<int> ZprimeTnPPairSelector_FromAOD::nSegments(edm::Event& event, reco::MuonRef MuRef) {

  std::vector<int> nsegments = {0,0,0,0};

  if( MuRef.isNonnull() ) {
    std::vector<int> dthits = countDTsegs(event, MuRef);
    std::vector<int> cschits = countCSCsegs(event, MuRef);
    for(int i=0; i<4; ++i) {
      nsegments[i] = dthits[i] + cschits[i];
    }
  }

  return nsegments;
}

std::pair<bool, reco::MuonRef> ZprimeTnPPairSelector_FromAOD::getMuonRef(edm::Event& event, const reco::CandidateBaseRef& lep) {

  edm::Handle< std::vector< reco::Muon > > recoMuons;
  event.getByLabel(reco_muon_src, recoMuons);

  bool isMatched = false;
  reco::MuonRef matchedMuRef = reco::MuonRef(recoMuons, 0);

  if (lep.isNonnull()) {
    const pat::Muon* mu = toConcretePtr<pat::Muon>(lep);

    int imucount = 0;
    for (std::vector<reco::Muon>::const_iterator imu = recoMuons->begin(); imu != recoMuons->end(); imu++) {
      if(imu->globalTrack()->pt() == mu->globalTrack()->pt() &&
         imu->globalTrack()->eta() == mu->globalTrack()->eta() &&
         imu->globalTrack()->phi() == mu->globalTrack()->phi() ) {
        isMatched = true;
        matchedMuRef = reco::MuonRef(recoMuons, imucount);
        break;
      }
      imucount++;
    }
  }
  return make_pair(isMatched, matchedMuRef);
}

std::vector<int> ZprimeTnPPairSelector_FromAOD::nHits(edm::Event& event, reco::MuonRef MuRef) {
  edm::Handle<edm::ValueMap<reco::MuonShower> > muonShowerInformationValueMap;
  event.getByLabel(muonshower_src, muonShowerInformationValueMap);

  std::vector<int> nhits = {0,0,0,0};
  reco::MuonShower muonShowerInformation = (*muonShowerInformationValueMap)[MuRef];
  for(int i=0; i<4; ++i) {
    nhits[i]  = (muonShowerInformation.nStationHits).at(i);
  }

  return nhits;
}

void ZprimeTnPPairSelector_FromAOD::produce(edm::Event& event, const edm::EventSetup& setup) {
  if(!ShutUp)  std::cout << "ZprimeTnPPairSelector_FromAOD::produce : Start!" << std::endl;

  edm::Handle<pat::CompositeCandidateCollection> cands;
  event.getByLabel(src, cands);

  //PV
  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel(vertex_src, vertices);
  PV = 0;
  for (reco::VertexCollection::const_iterator it = vertices->begin(), ite = vertices->end(); it != ite; ++it) {
    if (it->ndof() > 4 && fabs(it->z()) <= 24 && fabs(it->position().rho()) <= 2) {
      if (PV == 0){
        PV = &*it;
        break;
      }
    }
  }

  // does this get cached correctly? do we care?
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);

  std::unique_ptr<pat::CompositeCandidateCollection> new_cands(new pat::CompositeCandidateCollection);

  // Copy all the candidates that pass the specified cuts into the new
  // output vector. Also embed into the output dimuons any other
  // things that are best to just calculate once and for all.
  for (pat::CompositeCandidateCollection::const_iterator c = cands->begin(), ce = cands->end(); c != ce; ++c) {
    // Some cuts can be simply specified through the
    // StringCutSelector.
    if (!selector(*c))
      continue;

    // -- Selection for each Tag and Probe -- //

    //--- dpT/pT cut
    const reco::CandidateBaseRef& lep0 = dileptonDaughter(*c, 0);
    const reco::CandidateBaseRef& lep1 = dileptonDaughter(*c, 1);
    float lep0_dpt_over_pt = dpt_over_pt(lep0);
    float lep1_dpt_over_pt = dpt_over_pt(lep1);
    float lep0_veto_others_dphi = veto_others_dphi(event, lep0);
    float lep1_veto_others_dphi = veto_others_dphi(event, lep1);

    //if(!ShutUp)  std::cout << "ZprimeTnPPairSelector_FromAOD::produce dpt_over_pt : lep0 dpt_over_pt=" << lep0_dpt_over_pt << ", " << "lep1 dpt_over_pt=" << lep1_dpt_over_pt << std::endl;

    bool isTagAndProbe = TagAndProbeSelector(*c,
                                             lep0_dpt_over_pt,      lep1_dpt_over_pt,
                                             lep0_veto_others_dphi, lep1_veto_others_dphi);
    if( !isTagAndProbe )
      continue;


    // -- Selection for Tang and Probe pair -- //

    //---- Back-to-back cut to kill cosmics.
    std::pair<bool, float> cos_angle = back_to_back_cos_angle(*c);
    if (cut_on_back_to_back_cos_angle && !cos_angle.first)
      continue;

    //---- Loose common vertex chi2 cut
    std::pair<bool, CachingVertex<5> > vertex = vertex_constrained_fit(*c);
    if (cut_on_vertex_chi2 && !vertex.first)
      continue;

    //---- Cut on pT1 / pT2
    std::pair<bool, float> the_pt_ratio = pt_ratio(*c);
    if( cut_on_pt_ratio && !the_pt_ratio.first)
      continue;

    //---- Cut on deltaR(Tag, Probe)
    std::pair<bool, float> the_deltaR = dil_deltaR(*c);
    if( cut_on_dil_deltaR && !the_deltaR.first)
      continue;

    //--- Adding nShower variables
    // HERE
    bool testCSC = false;

    std::vector<int> lep0_nHits = {-1,-1,-1,-1};
    std::vector<int> lep0_nSegs = {-1,-1,-1,-1};
    std::pair<bool, reco::MuonRef> pair_lep0Ref = getMuonRef(event, lep0);
    int lep0_nMuonHits = -1;
    int lep0_nValidMuonHits = -1;

    if(pair_lep0Ref.first) {
      lep0_nMuonHits = (int)((pair_lep0Ref.second)->standAloneMuon()->hitPattern().numberOfMuonHits());
      lep0_nValidMuonHits = (int)((pair_lep0Ref.second)->standAloneMuon()->hitPattern().numberOfValidMuonHits());
      lep0_nHits = nHits(event, pair_lep0Ref.second);
      lep0_nSegs = nSegments(event, pair_lep0Ref.second);
      if(testCSC) {
        std::cout << "lep0 testCSC" << std::endl;
        testCSChits(event, (pair_lep0Ref.second)->standAloneMuon());
      }
    }
    else
      std::cout <<  "ZprimeTnPPairSelector_FromAOD::produce : No RECO MuonRef found for lep0" << std::endl;


    std::vector<int> lep1_nHits = {-1,-1,-1,-1};
    std::vector<int> lep1_nSegs = {-1,-1,-1,-1};
    std::pair<bool, reco::MuonRef> pair_lep1Ref = getMuonRef(event, lep1);
    int lep1_nMuonHits = -1;
    int lep1_nValidMuonHits = -1;

    if(pair_lep1Ref.first) {
      lep1_nMuonHits = (int)((pair_lep1Ref.second)->standAloneMuon()->hitPattern().numberOfMuonHits());
      lep1_nValidMuonHits = (int)((pair_lep1Ref.second)->standAloneMuon()->hitPattern().numberOfValidMuonHits());
      lep1_nHits = nHits(event, pair_lep1Ref.second);
      lep1_nSegs = nSegments(event, pair_lep1Ref.second);
      if(testCSC) {
        std::cout << "lep1 testCSC" << std::endl;
        testCSChits(event, (pair_lep1Ref.second)->standAloneMuon());
      }
    }
    else
      std::cout <<  "ZprimeTnPPairSelector_FromAOD::produce : No RECO MuonRef found for lep1" << std::endl;

    if(!ShutUp) {
      const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
      const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
      std::cout << "ZprimeTnPPairSelector_FromAOD::produce : RECO MuonRef" << std::endl;

      std::cout << "\tmu0" << std::endl;
      std::cout << "\t\t" << mu0->globalTrack()->pt() << ", " << mu0->globalTrack()->eta() << ", " << mu0->globalTrack()->phi() << std::endl;
      std::cout << "\tRECO Ref matched to mu0" << std::endl;
      std::cout << "\t\t" << pair_lep0Ref.second->globalTrack()->pt() << ", " << pair_lep0Ref.second->globalTrack()->eta() << ", " << pair_lep0Ref.second->globalTrack()->phi() << std::endl;

      std::cout << "\tlep0_nMuonHits : " << lep0_nMuonHits << std::endl;
      std::cout << "\tlep0_nValidMuonHits : " << lep0_nValidMuonHits << std::endl;

      std::cout << "\tlep0_nHits" << std::endl;
      std::cout << "\t\t" << lep0_nHits[0] << ", " << lep0_nHits[1] << ", " << lep0_nHits[2] << ", " << lep0_nHits[3] << std::endl;
      std::cout << "\tlep0_nSegs" << std::endl;
      std::cout << "\t\t" << lep0_nSegs[0] << ", " << lep0_nSegs[1] << ", " << lep0_nSegs[2] << ", " << lep0_nSegs[3] << std::endl;

      std::cout << "\n\tmu1" << std::endl;
      std::cout << "\t\t" << mu1->globalTrack()->pt() << ", " << mu1->globalTrack()->eta() << ", " << mu1->globalTrack()->phi() << std::endl;
      std::cout << "\tRECO Ref matched to mu1" << std::endl;
      std::cout << "\t\t" << pair_lep1Ref.second->globalTrack()->pt() << ", " << pair_lep1Ref.second->globalTrack()->eta() << ", " << pair_lep1Ref.second->globalTrack()->phi() << std::endl;

      std::cout << "\tlep1_nMuonHits : " << lep1_nMuonHits << std::endl;
      std::cout << "\tlep1_nValidMuonHits : " << lep1_nValidMuonHits << std::endl;

      std::cout << "\tlep1_nHits" << std::endl;
      std::cout << "\t\t" << lep1_nHits[0] << ", " << lep1_nHits[1] << ", " << lep1_nHits[2] << ", " << lep1_nHits[3] << std::endl;
      std::cout << "\tlep1_nSegs" << std::endl;
      std::cout << "\t\t" << lep1_nSegs[0] << ", " << lep1_nSegs[1] << ", " << lep1_nSegs[2] << ", " << lep1_nSegs[3] << std::endl;

      std::cout << std::endl;
    }

    // Save the dilepton since it passed the cuts, and store the cut
    // variables and other stuff for use later.
    new_cands->push_back(*c);
    new_cands->back().addUserFloat("lep0_dpt_over_pt",   lep0_dpt_over_pt);
    new_cands->back().addUserFloat("lep1_dpt_over_pt",   lep1_dpt_over_pt);

    new_cands->back().addUserFloat("lep0_veto_others_dphi", lep0_veto_others_dphi);
    new_cands->back().addUserFloat("lep1_veto_others_dphi", lep1_veto_others_dphi);

    new_cands->back().addUserFloat("cos_angle",   cos_angle.second);
    embed_vertex_constrained_fit(new_cands->back(), vertex.second);

    new_cands->back().addUserFloat("pt_ratio", the_pt_ratio.second);
    new_cands->back().addUserFloat("dil_deltaR", the_deltaR.second);

    new_cands->back().addUserInt("lep0_nMuonHits", lep0_nMuonHits);
    new_cands->back().addUserInt("lep0_nValidMuonHits", lep0_nValidMuonHits);

    new_cands->back().addUserInt("lep0_nHits_1", lep0_nHits[0]);
    new_cands->back().addUserInt("lep0_nHits_2", lep0_nHits[1]);
    new_cands->back().addUserInt("lep0_nHits_3", lep0_nHits[2]);
    new_cands->back().addUserInt("lep0_nHits_4", lep0_nHits[3]);

    new_cands->back().addUserInt("lep0_nSegs_1", lep0_nSegs[0]);
    new_cands->back().addUserInt("lep0_nSegs_2", lep0_nSegs[1]);
    new_cands->back().addUserInt("lep0_nSegs_3", lep0_nSegs[2]);
    new_cands->back().addUserInt("lep0_nSegs_4", lep0_nSegs[3]);

    new_cands->back().addUserInt("lep1_nMuonHits", lep1_nMuonHits);
    new_cands->back().addUserInt("lep1_nValidMuonHits", lep1_nValidMuonHits);

    new_cands->back().addUserInt("lep1_nHits_1", lep1_nHits[0]);
    new_cands->back().addUserInt("lep1_nHits_2", lep1_nHits[1]);
    new_cands->back().addUserInt("lep1_nHits_3", lep1_nHits[2]);
    new_cands->back().addUserInt("lep1_nHits_4", lep1_nHits[3]);

    new_cands->back().addUserInt("lep1_nSegs_1", lep1_nSegs[0]);
    new_cands->back().addUserInt("lep1_nSegs_2", lep1_nSegs[1]);
    new_cands->back().addUserInt("lep1_nSegs_3", lep1_nSegs[2]);
    new_cands->back().addUserInt("lep1_nSegs_4", lep1_nSegs[3]);


    if(!ShutUp) {
      if( (lep0_nSegs[0]+lep0_nSegs[1]+lep0_nSegs[2]+lep0_nSegs[3]) > 0 ||
          (lep1_nSegs[0]+lep1_nSegs[1]+lep1_nSegs[2]+lep1_nSegs[3]) > 0
       )
        std::cout << "Event ID  " << event.id().run() << ":" << event.id().luminosityBlock() << ":" << event.id().event() << std::endl;
    }
  }

  // Sort candidates so we keep either the ones with higher-pT
  // muons or the ones with larger invariant mass.
  if(sort_by_pt)
    sort(new_cands->begin(), new_cands->end(), lepton_pt_sort());
  else
    sort(new_cands->begin(), new_cands->end(), reverse_mass_sort());

  // Remove cands of lower invariant mass that are comprised of a
  // lepton that has been used by a higher invariant mass one.
  if (do_remove_overlap)
    remove_overlap(*new_cands);

  // Only return the maximum number of candidates specified.
  if (new_cands->size() > max_candidates)
    new_cands->erase(new_cands->begin() + max_candidates, new_cands->end());

  if(!ShutUp && new_cands->size()>0) {
    std::cout << "ZprimeTnPPairSelector_FromAOD::produce : TnP Pair found!" << std::endl;
    std::cout << "\tmass : " << new_cands->at(0).mass() << "\n"
              << "\tlep0_dpt_over_pt : " << new_cands->at(0).userFloat("lep0_dpt_over_pt") << "\n"
              << "\tlep1_dpt_over_pt : " << new_cands->at(0).userFloat("lep1_dpt_over_pt") << "\n"

              << "\tlep0_veto_others_dphi : " << new_cands->at(0).userFloat("lep0_veto_others_dphi") << "\n"
              << "\tlep1_veto_others_dphi : " << new_cands->at(0).userFloat("lep1_veto_others_dphi") << "\n"

              << "\tcos_angle : " << new_cands->at(0).userFloat("cos_angle") << "\n"
              << "\tpt_ratio : " << new_cands->at(0).userFloat("pt_ratio") << "\n"
              << "\tdil_deltaR : " << new_cands->at(0).userFloat("dil_deltaR") << "\n";
  }

  event.put(move(new_cands));
}

DEFINE_FWK_MODULE(ZprimeTnPPairSelector_FromAOD);
