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

class ZprimeTnPPairSelector : public edm::EDProducer {
public:
  explicit ZprimeTnPPairSelector(const edm::ParameterSet&);

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

  bool                               TagAndProbeSelector(const pat::CompositeCandidate&, float, float) const;

  std::pair<bool, float>             back_to_back_cos_angle(const pat::CompositeCandidate&) const;
  std::pair<bool, CachingVertex<5> > vertex_constrained_fit(const pat::CompositeCandidate&) const;
  float                              dpt_over_pt(const reco::CandidateBaseRef&) const;

  std::pair<bool, float>             pt_ratio(const pat::CompositeCandidate&) const;
  std::pair<bool, float>             dil_deltaR(const pat::CompositeCandidate&) const;

  // If the variable to embed in the methods above is a simple int or
  // float or is going to be embedded wholesale with the generic
  // userData mechanism, we'll do those explicitly in the loop in
  // produce. Otherwise if it's more complicated, the methods here
  // take care of it.
  void embed_vertex_constrained_fit(pat::CompositeCandidate&, const CachingVertex<5>& vtx) const;

  const edm::InputTag src;
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

  const bool samePV;

  edm::ESHandle<TransientTrackBuilder> ttkb;

  const bool ShutUp;
};

ZprimeTnPPairSelector::ZprimeTnPPairSelector(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),
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

    samePV(cfg.getParameter<bool>("samePV")),

    ShutUp(cfg.getParameter<bool>("ShutUp"))
{
 consumes<pat::CompositeCandidateCollection>(src);
 consumes<reco::VertexCollection>(vertex_src);
 produces<pat::CompositeCandidateCollection>();

}

void ZprimeTnPPairSelector::remove_overlap(pat::CompositeCandidateCollection& cands) const {
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

std::vector<reco::TransientTrack> ZprimeTnPPairSelector::get_transient_tracks(const pat::CompositeCandidate& dil) const {
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

std::pair<bool, float> ZprimeTnPPairSelector::back_to_back_cos_angle(const pat::CompositeCandidate& dil) const {
  // Back-to-back cut to kill cosmics.
  assert(dil.numberOfDaughters() == 2);
  const float cos_angle = dil.daughter(0)->momentum().Dot(dil.daughter(1)->momentum()) / dil.daughter(0)->p() / dil.daughter(1)->p();
  return std::make_pair(cos_angle >= back_to_back_cos_angle_min, cos_angle);
}

std::pair<bool, CachingVertex<5> > ZprimeTnPPairSelector::vertex_constrained_fit(const pat::CompositeCandidate& dil) const {
  // Loose common vertex chi2 cut.
  assert(dil.numberOfDaughters() == 2);
  if (abs(dil.daughter(0)->pdgId()) != 13 || abs(dil.daughter(1)->pdgId()) != 13)
    return std::make_pair(true, CachingVertex<5>()); // pass objects we don't know how to cut on, i.e. e-mu dileptons

  KalmanVertexFitter kvf(true);
  CachingVertex<5> v = kvf.vertex(get_transient_tracks(dil));

  return std::make_pair(v.isValid() && v.totalChiSquared()/v.degreesOfFreedom() <= vertex_chi2_max, v);
}

void ZprimeTnPPairSelector::embed_vertex_constrained_fit(pat::CompositeCandidate& dil, const CachingVertex<5>& vtx) const {
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

float ZprimeTnPPairSelector::dpt_over_pt(const reco::CandidateBaseRef& lep) const {
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
std::pair<bool, float> ZprimeTnPPairSelector::pt_ratio(const pat::CompositeCandidate& dil) const {
  double pt0 = dil.daughter(0)->pt();
  double pt1 = dil.daughter(1)->pt();

  double the_pt_ratio = pt0 > pt1 ? (pt0/pt1) : (pt1/pt0);

  //if(!ShutUp)  std::cout << "ZprimeTnPPairSelector::pt_ratio[" << pt_ratio_max << "] : pt_ratio = " << the_pt_ratio << ", " << (the_pt_ratio < pt_ratio_max) << std::endl;

  return std::make_pair( the_pt_ratio < pt_ratio_max , the_pt_ratio );
}

///
std::pair<bool, float> ZprimeTnPPairSelector::dil_deltaR(const pat::CompositeCandidate& dil) const {

  double the_deltaR = reco::deltaR( *dil.daughter(0), *dil.daughter(1) );

  //if(!ShutUp)  std::cout << "ZprimeTnPPairSelector::dil_deltaR[" << dil_deltaR_min << "] : deltaR = " << the_deltaR << ", " << (the_deltaR > dil_deltaR_min) << std::endl;

  return std::make_pair( the_deltaR > dil_deltaR_min , the_deltaR );
}

///
bool ZprimeTnPPairSelector::TagAndProbeSelector(const pat::CompositeCandidate& dil, float lep0_dpt_over_pt, float lep1_dpt_over_pt) const {
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

      if(!ShutUp)  std::cout << "ZprimeTnPPairSelector::TagAndProbeSelector[" << (isTag0Probe1 || isTag1Probe0) << "] : isTag0=" 
        << ( tag_selector(*mu0)   && (lep0_dpt_over_pt<tag_dpt_over_pt_max) && (fabs(mu0->muonBestTrack()->dz( PV->position() )) < tag_dz_max) ) << " " << "isTag1=" 
        << ( tag_selector(*mu1)   && (lep1_dpt_over_pt<tag_dpt_over_pt_max) && (fabs(mu1->muonBestTrack()->dz( PV->position() )) < tag_dz_max) ) << " " << "isProbe0=" 
        << ( probe_selector(*mu0) && (lep0_dpt_over_pt<probe_dpt_over_pt_max) ) << " " << "isProbe1=" 
        << ( probe_selector(*mu1) && (lep1_dpt_over_pt<probe_dpt_over_pt_max) ) << std::endl;

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

void ZprimeTnPPairSelector::produce(edm::Event& event, const edm::EventSetup& setup) {
  if(!ShutUp)  std::cout << "ZprimeTnPPairSelector::produce : Start!" << std::endl;

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

    //if(!ShutUp)  std::cout << "ZprimeTnPPairSelector::produce dpt_over_pt : lep0 dpt_over_pt=" << lep0_dpt_over_pt << ", " << "lep1 dpt_over_pt=" << lep1_dpt_over_pt << std::endl;

    bool isTagAndProbe = TagAndProbeSelector(*c, lep0_dpt_over_pt, lep1_dpt_over_pt);
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

    // Save the dilepton since it passed the cuts, and store the cut
    // variables and other stuff for use later.
    new_cands->push_back(*c);
    new_cands->back().addUserFloat("lep0_dpt_over_pt",   lep0_dpt_over_pt);
    new_cands->back().addUserFloat("lep1_dpt_over_pt",   lep1_dpt_over_pt);
    new_cands->back().addUserFloat("cos_angle",   cos_angle.second);
    embed_vertex_constrained_fit(new_cands->back(), vertex.second);

    new_cands->back().addUserFloat("pt_ratio", the_pt_ratio.second);
    new_cands->back().addUserFloat("dil_deltaR", the_deltaR.second);
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
    std::cout << "ZprimeTnPPairSelector::produce : TnP Pair found!" << std::endl;
    std::cout << "\tmass : " << new_cands->at(0).mass() << "\n"
              << "\tlep0_dpt_over_pt : " << new_cands->at(0).userFloat("lep0_dpt_over_pt") << "\n"
              << "\tlep1_dpt_over_pt : " << new_cands->at(0).userFloat("lep1_dpt_over_pt") << "\n"
              << "\tcos_angle : " << new_cands->at(0).userFloat("cos_angle") << "\n"
              << "\tpt_ratio : " << new_cands->at(0).userFloat("pt_ratio") << "\n"
              << "\tdil_deltaR : " << new_cands->at(0).userFloat("dil_deltaR") << "\n";
  }

  event.put(move(new_cands));
}

DEFINE_FWK_MODULE(ZprimeTnPPairSelector);
