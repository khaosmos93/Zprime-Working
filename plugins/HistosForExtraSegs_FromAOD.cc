#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class Zprime2muHistosForExtraSegs_FromAOD : public edm::EDAnalyzer {
 public:
  explicit Zprime2muHistosForExtraSegs_FromAOD(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void getBSandPV(const edm::Event&);
  std::pair<int, int> getNShowers(const pat::CompositeCandidate&);  // first : lep0,  second : lep1
  void fillTnPControlHistos(const pat::CompositeCandidate&, const reco::CandidateBaseRef&, const reco::CandidateBaseRef&, int, int, int, bool );

  edm::InputTag dilepton_src;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;
  const bool use_bs_and_pv;
  const reco::BeamSpot* beamspot;
  const reco::Vertex*   vertex;
  int                   nVtx;

  StringCutObjectSelector<pat::Muon> tag_selector;
  StringCutObjectSelector<pat::Muon> probe_selector;
  StringCutObjectSelector<pat::Muon> passing_probe_selector;
  StringCutObjectSelector<pat::Muon> comparison_probe_selector;

  double minMass;
  double maxMass;
  double tag_dpt_over_pt_max;
  double tag_dz_max;
  double probe_dpt_over_pt_max;
  double probe_pt_min;
  double passing_probe_dpt_over_pt_max;
  double passing_probe_dz_max;

  int nshowers_threshold_min;

  double comparison_probe_dpt_over_pt_max;
  double comparison_probe_dz_max;

  bool   _usePrescaleWeight;
  int    _prescaleWeight;
  double _eventWeight;
  bool   _useMadgraphWeight;
  double _madgraphWeight;

  edm::InputTag pileup_src;
  double _pileupWeight;
  std::vector<double> vec_PileUpWeight;

  double _totalWeight;

  const bool ShutUp;


  // -- Histograms -- //
  TH1F* NBeamSpot;
  TH1F* NVertices;
  TH1F* NDileptons;
  TH1F* WeightMadGraph;

  TH1F *TagPt;
  TH1F *TagEta;
  TH1F *TagPhi;

  /*TH1F* TagAbsTkIso;
  TH1F* TagRelTkIso;
  TH1F* TagChi2dof;
  TH1F* TagTrackDXYBS;
  TH1F* TagTrackDZBS;
  TH1F* TagTrackDXYPV;
  TH1F* TagTrackDZPV;
  TH1F* TagNPxHits;
  TH1F* TagNStHits;
  TH1F* TagNTkHits;
  TH1F* TagNMuHits;
  TH1F* TagNHits;
  TH1F* TagNPxLayers;
  TH1F* TagNStLayers;
  TH1F* TagNTkLayers;*/

  TH1F *ProbePt;
  TH1F *ProbeEta;
  TH1F *ProbePhi;
  TH1F *ProbeNVertices;
  TH2F *ProbeEtaPhi;

  /*TH1F* ProbeAbsTkIso;
  TH1F* ProbeRelTkIso;
  TH1F* ProbeChi2dof;
  TH1F* ProbeTrackDXYBS;
  TH1F* ProbeTrackDZBS;
  TH1F* ProbeTrackDXYPV;
  TH1F* ProbeTrackDZPV;
  TH1F* ProbeNPxHits;
  TH1F* ProbeNStHits;
  TH1F* ProbeNTkHits;
  TH1F* ProbeNMuHits;
  TH1F* ProbeNHits;
  TH1F* ProbeNPxLayers;
  TH1F* ProbeNStLayers;
  TH1F* ProbeNTkLayers;*/

  TH2F *ProbeNHitsBP;
  TH2F *ProbeNHitsEP;

  TH2F *ProbeNSegsBP;
  TH2F *ProbeNSegsBP1St;
  TH2F *ProbeNSegsBP2St;
  TH2F *ProbeNSegsBP3St;
  TH2F *ProbeNSegsBP4St;

  TH2F *ProbeNSegsEP;
  TH2F *ProbeNSegsEP1St;
  TH2F *ProbeNSegsEP2St;
  TH2F *ProbeNSegsEP3St;
  TH2F *ProbeNSegsEP4St;

  TH2F *ProbeNShowerBP;
  TH2F *ProbeNShowerEP;

  TH1F *PassingProbePt;
  TH1F *PassingProbeEta;
  TH1F *PassingProbePhi;
  TH1F *PassingProbeNVertices;
  TH2F *PassingProbeEtaPhi;

  TH2F *PassingProbeNHitsBP;
  TH2F *PassingProbeNHitsEP;
  TH2F *PassingProbeNSegsBP;
  TH2F *PassingProbeNSegsEP;
  TH2F *PassingProbeNShowerBP;
  TH2F *PassingProbeNShowerEP;

  TH1F *FailingProbePt;
  TH1F *FailingProbeEta;
  TH1F *FailingProbePhi;
  TH1F *FailingProbeNVertices;
  TH2F *FailingProbeEtaPhi;

  TH2F *FailingProbeNHitsBP;
  TH2F *FailingProbeNHitsEP;
  TH2F *FailingProbeNSegsBP;
  TH2F *FailingProbeNSegsEP;
  TH2F *FailingProbeNShowerBP;
  TH2F *FailingProbeNShowerEP;

  TH1F *PairNoPtMass;
  TH1F *PairNoPtPt;
  TH1F *PairNoPtEta;
  TH1F *PairNoPtRap;

  TH1F *PairMass;
  TH1F *PairPt;
  TH1F *PairEta;
  TH1F *PairRap;

  TH1F *PassingPairMass;
  TH1F *PassingPairPt;
  TH1F *PassingPairEta;
  TH1F *PassingPairRap;

  TH1F *FailingPairMass;
  TH1F *FailingPairPt;
  TH1F *FailingPairEta;
  TH1F *FailingPairRap;

};

Zprime2muHistosForExtraSegs_FromAOD::Zprime2muHistosForExtraSegs_FromAOD(const edm::ParameterSet& cfg)
  : dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    use_bs_and_pv(cfg.getParameter<bool>("use_bs_and_pv")),
    beamspot(0),
    vertex(0),
    nVtx(0),

    tag_selector(cfg.getParameter<std::string>("tag_cut")),
    probe_selector(cfg.getParameter<std::string>("probe_cut")),
    passing_probe_selector(cfg.getParameter<std::string>("passing_probe_cut")),
    comparison_probe_selector(cfg.getParameter<std::string>("comparison_probe_cut")),

    minMass(cfg.getParameter<double>("minMass")),
    maxMass(cfg.getParameter<double>("maxMass")),

    tag_dpt_over_pt_max(cfg.getParameter<double>("tag_dpt_over_pt_max")),
    tag_dz_max(cfg.getParameter<double>("tag_dz_max")),
    probe_dpt_over_pt_max(cfg.getParameter<double>("probe_dpt_over_pt_max")),
    probe_pt_min(cfg.getParameter<double>("probe_pt_min")),
    passing_probe_dpt_over_pt_max(cfg.getParameter<double>("passing_probe_dpt_over_pt_max")),
    passing_probe_dz_max(cfg.getParameter<double>("passing_probe_dz_max")),

    nshowers_threshold_min(cfg.getParameter<int>("nshowers_threshold_min")),

    comparison_probe_dpt_over_pt_max(cfg.getParameter<double>("comparison_probe_dpt_over_pt_max")),
    comparison_probe_dz_max(cfg.getParameter<double>("comparison_probe_dz_max")),

    _usePrescaleWeight(cfg.getUntrackedParameter<bool>("usePrescaleWeight",false)),
    _prescaleWeight(1),
    _useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    _madgraphWeight(1.),

    pileup_src(cfg.getParameter<edm::InputTag>("pileup_src")),
    _pileupWeight(1.),
    vec_PileUpWeight(cfg.getParameter<std::vector<double>>("vec_PileUpWeight")),

    _totalWeight(1.),

    ShutUp(cfg.getParameter<bool>("ShutUp"))
{
  consumes<pat::CompositeCandidateCollection>(dilepton_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertex_src);
  consumes<std::vector<PileupSummaryInfo> >(pileup_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  NBeamSpot = fs->make<TH1F>("NBeamSpot", "# beamspots/event",  2, 0,  2);
  NVertices = fs->make<TH1F>("NVertices", "# vertices/event",  200, 0, 200);

  // Dilepton multiplicity.
  NDileptons = fs->make<TH1F>("NDileptons", "# dileptons/event", 10, 0, 10);

  // Generator weights
  WeightMadGraph = fs->make<TH1F>("WeightMadGraph", "weight per event", 4, -2,2);

  // Tag
  TagPt = fs->make<TH1F>("TagPt", "Tag pT", 10000, 0, 10000);
  TagEta = fs->make<TH1F>("TagEta", "Tag #eta",    100, -5, 5);
  TagPhi = fs->make<TH1F>("TagPhi", "Tag #phi", 100, -TMath::Pi(), TMath::Pi());

  /*TagAbsTkIso = fs->make<TH1F>("TagAbsTkIso", "Tag Iso. (#Delta R < 0.3) #Sigma pT", 1000, 0, 1000);
  TagRelTkIso = fs->make<TH1F>("TagRelTkIso", "Tag Iso. (#Delta R < 0.3) #Sigma pT / tk. pT", 500, 0, 5);
  TagChi2dof = fs->make<TH1F>("TagChi2dof", "Tag #chi^{2}/dof", 500, 0, 50);
  TagTrackDXYBS = fs->make<TH1F>("TagTrackDXYBS", "Tag |dxy wrt BS|", 10000, 0, 2);
  TagTrackDZBS = fs->make<TH1F>("TagTrackDZBS", "Tag |dz wrt BS|", 10000, 0, 20);
  TagTrackDXYPV = fs->make<TH1F>("TagTrackDXYPV", "Tag |dxy wrt PV|", 10000, 0, 2);
  TagTrackDZPV = fs->make<TH1F>("TagTrackDZPV", "Tag |dz wrt PV|", 10000, 0, 20);
  TagNPxHits = fs->make<TH1F>("TagNPxHits", "Tag # pixel hits", 10, 0,  10);
  TagNStHits = fs->make<TH1F>("TagNStHits", "Tag # strip hits", 40, 0, 40);
  TagNTkHits = fs->make<TH1F>("TagNTkHits", "Tag # tracker hits", 50, 0, 50);
  TagNMuHits = fs->make<TH1F>("TagNMuHits", "Tag # muon hits", 60, 0, 60);
  TagNHits = fs->make<TH1F>("TagNHits", "Tag # hits", 80, 0, 80);
  TagNPxLayers = fs->make<TH1F>("TagNPxLayers", "Tag # pixel layers", 10, 0, 10);
  TagNStLayers = fs->make<TH1F>("TagNStLayers", "Tag # strip layers", 20, 0, 20);
  TagNTkLayers = fs->make<TH1F>("TagNTkLayers", "Tag # tracker layers", 30, 0, 30);*/

  // Probe
  ProbePt = fs->make<TH1F>("ProbePt", "Probe pT", 10000, 0, 10000);
  ProbeEta = fs->make<TH1F>("ProbeEta", "Probe #eta",    100, -5, 5);
  ProbePhi = fs->make<TH1F>("ProbePhi", "Probe #phi", 100, -TMath::Pi(), TMath::Pi());
  ProbeNVertices = fs->make<TH1F>("ProbeNVertices", "Probe # vertices/event",  200, 0, 200);
  ProbeEtaPhi = fs->make<TH2F>("ProbeEtaPhi", "Probe #eta #phi",    100, -5, 5, 100, -TMath::Pi(), TMath::Pi());

  /*ProbeAbsTkIso = fs->make<TH1F>("ProbeAbsTkIso", "Probe Iso. (#Delta R < 0.3) #Sigma pT", 1000, 0, 1000);
  ProbeRelTkIso = fs->make<TH1F>("ProbeRelTkIso", "Probe Iso. (#Delta R < 0.3) #Sigma pT / tk. pT", 500, 0, 5);
  ProbeChi2dof = fs->make<TH1F>("ProbeChi2dof", "Probe #chi^{2}/dof", 500, 0, 50);
  ProbeTrackDXYBS = fs->make<TH1F>("ProbeTrackDXYBS", "Probe |dxy wrt BS|", 10000, 0, 2);
  ProbeTrackDZBS = fs->make<TH1F>("ProbeTrackDZBS", "Probe |dz wrt BS|", 10000, 0, 20);
  ProbeTrackDXYPV = fs->make<TH1F>("ProbeTrackDXYPV", "Probe |dxy wrt PV|", 10000, 0, 2);
  ProbeTrackDZPV = fs->make<TH1F>("ProbeTrackDZPV", "Probe |dz wrt PV|", 10000, 0, 20);
  ProbeNPxHits = fs->make<TH1F>("ProbeNPxHits", "Probe # pixel hits", 10, 0,  10);
  ProbeNStHits = fs->make<TH1F>("ProbeNStHits", "Probe # strip hits", 40, 0, 40);
  ProbeNTkHits = fs->make<TH1F>("ProbeNTkHits", "Probe # tracker hits", 50, 0, 50);
  ProbeNMuHits = fs->make<TH1F>("ProbeNMuHits", "Probe # muon hits", 60, 0, 60);
  ProbeNHits = fs->make<TH1F>("ProbeNHits", "Probe # hits", 80, 0, 80);
  ProbeNPxLayers = fs->make<TH1F>("ProbeNPxLayers", "Probe # pixel layers", 10, 0, 10);
  ProbeNStLayers = fs->make<TH1F>("ProbeNStLayers", "Probe # strip layers", 20, 0, 20);
  ProbeNTkLayers = fs->make<TH1F>("ProbeNTkLayers", "Probe # tracker layers", 30, 0, 30);*/

  double arr_PBins[5] = {0., 200., 500., 1000., 10000.};

  ProbeNHitsBP = fs->make<TH2F>("ProbeNHitsBP", "Probe Barrel # hits |P|",    2000, 0, 2000, 4, arr_PBins);
  ProbeNHitsEP = fs->make<TH2F>("ProbeNHitsEP", "Probe Endcap # hits |P|",    2000, 0, 2000, 4, arr_PBins);
  
  ProbeNSegsBP = fs->make<TH2F>("ProbeNSegsBP", "Probe Barrel # segments |P|",    500, 0, 500, 4, arr_PBins);
  ProbeNSegsBP1St = fs->make<TH2F>("ProbeNSegsBP1St", "Probe Barrel 1St # segments |P|",    500, 0, 500, 4, arr_PBins);
  ProbeNSegsBP2St = fs->make<TH2F>("ProbeNSegsBP2St", "Probe Barrel 2St # segments |P|",    500, 0, 500, 4, arr_PBins);
  ProbeNSegsBP3St = fs->make<TH2F>("ProbeNSegsBP3St", "Probe Barrel 3St # segments |P|",    500, 0, 500, 4, arr_PBins);
  ProbeNSegsBP4St = fs->make<TH2F>("ProbeNSegsBP4St", "Probe Barrel 4St # segments |P|",    500, 0, 500, 4, arr_PBins);

  ProbeNSegsEP = fs->make<TH2F>("ProbeNSegsEP", "Probe Endcap # segments |P|",    500, 0, 500, 4, arr_PBins);
  ProbeNSegsEP1St = fs->make<TH2F>("ProbeNSegsEP1St", "Probe Endcap 1St # segments |P|",    500, 0, 500, 4, arr_PBins);
  ProbeNSegsEP2St = fs->make<TH2F>("ProbeNSegsEP2St", "Probe Endcap 2St # segments |P|",    500, 0, 500, 4, arr_PBins);
  ProbeNSegsEP3St = fs->make<TH2F>("ProbeNSegsEP3St", "Probe Endcap 3St # segments |P|",    500, 0, 500, 4, arr_PBins);
  ProbeNSegsEP4St = fs->make<TH2F>("ProbeNSegsEP4St", "Probe Endcap 4St # segments |P|",    500, 0, 500, 4, arr_PBins);

  ProbeNShowerBP = fs->make<TH2F>("ProbeNShowerBP", "Probe # showers Barrel |P|",    6, -1, 5, 4, arr_PBins);
  ProbeNShowerEP = fs->make<TH2F>("ProbeNShowerEP", "Probe # showers Endcap |P|",    6, -1, 5, 4, arr_PBins);

  PassingProbePt = fs->make<TH1F>("PassingProbePt", "PassingProbe pT", 10000, 0, 10000);
  PassingProbeEta = fs->make<TH1F>("PassingProbeEta", "PassingProbe #eta",    100, -5, 5);
  PassingProbePhi = fs->make<TH1F>("PassingProbePhi", "PassingProbe #phi", 100, -TMath::Pi(), TMath::Pi());
  PassingProbeNVertices = fs->make<TH1F>("PassingProbeNVertices", "PassingProbe # vertices/event",  200, 0, 200);
  PassingProbeEtaPhi = fs->make<TH2F>("PassingProbeEtaPhi", "PassingProbe #eta #phi",    100, -5, 5, 100, -TMath::Pi(), TMath::Pi());

  PassingProbeNHitsBP = fs->make<TH2F>("PassingProbeNHitsBP", "PassingProbe shower Barrel # hits |P|",    2000, 0, 2000, 4, arr_PBins);
  PassingProbeNHitsEP = fs->make<TH2F>("PassingProbeNHitsEP", "PassingProbe shower Endcap # hits |P|",    2000, 0, 2000, 4, arr_PBins);
  PassingProbeNSegsBP = fs->make<TH2F>("PassingProbeNSegsBP", "PassingProbe shower Barrel # segments |P|",    500, 0, 500, 4, arr_PBins);
  PassingProbeNSegsEP = fs->make<TH2F>("PassingProbeNSegsEP", "PassingProbe shower Endcap # segments |P|",    500, 0, 500, 4, arr_PBins);
  PassingProbeNShowerBP = fs->make<TH2F>("PassingProbeNShowerBP", "PassingProbe # showers Barrel |P|",    6, -1, 5, 4, arr_PBins);
  PassingProbeNShowerEP = fs->make<TH2F>("PassingProbeNShowerEP", "PassingProbe # showers Endcap |P|",    6, -1, 5, 4, arr_PBins);

  FailingProbePt = fs->make<TH1F>("FailingProbePt", "FailingProbe pT", 10000, 0, 10000);
  FailingProbeEta = fs->make<TH1F>("FailingProbeEta", "FailingProbe #eta",    100, -5, 5);
  FailingProbePhi = fs->make<TH1F>("FailingProbePhi", "FailingProbe #phi", 100, -TMath::Pi(), TMath::Pi());
  FailingProbeNVertices = fs->make<TH1F>("FailingProbeNVertices", "FailingProbe # vertices/event",  200, 0, 200);
  FailingProbeEtaPhi = fs->make<TH2F>("FailingProbeEtaPhi", "FailingProbe #eta #phi",    100, -5, 5, 100, -TMath::Pi(), TMath::Pi());

  FailingProbeNHitsBP = fs->make<TH2F>("FailingProbeNHitsBP", "FailingProbe shower Barrel # hits |P|",    2000, 0, 2000, 4, arr_PBins);
  FailingProbeNHitsEP = fs->make<TH2F>("FailingProbeNHitsEP", "FailingProbe shower Endcap # hits |P|",    2000, 0, 2000, 4, arr_PBins);
  FailingProbeNSegsBP = fs->make<TH2F>("FailingProbeNSegsBP", "FailingProbe shower Barrel # segments |P|",    500, 0, 500, 4, arr_PBins);
  FailingProbeNSegsEP = fs->make<TH2F>("FailingProbeNSegsEP", "FailingProbe shower Endcap # segments |P|",    500, 0, 500, 4, arr_PBins);
  FailingProbeNShowerBP = fs->make<TH2F>("FailingProbeNShowerBP", "FailingProbe # showers Barrel |P|",    6, -1, 5, 4, arr_PBins);
  FailingProbeNShowerEP = fs->make<TH2F>("FailingProbeNShowerEP", "FailingProbe # showers Endcap |P|",    6, -1, 5, 4, arr_PBins);

  // TnP pair
  PairNoPtMass = fs->make<TH1F>("PairNoPtMass", "TnP PairNoPt mass", 20000, 0, 20000);
  PairNoPtPt = fs->make<TH1F>("PairNoPtPt", "TnP PairNoPt pT", 10000, 0, 10000);
  PairNoPtEta = fs->make<TH1F>("PairNoPtEta", "TnP PairNoPt #eta",    100, -5, 5);
  PairNoPtRap = fs->make<TH1F>("PairNoPtRap", "TnP PairNoPt y", 100, -5, 5);

  PairMass = fs->make<TH1F>("PairMass", "TnP Pair mass", 20000, 0, 20000);
  PairPt = fs->make<TH1F>("PairPt", "TnP Pair pT", 10000, 0, 10000);
  PairEta = fs->make<TH1F>("PairEta", "TnP Pair #eta",    100, -5, 5);
  PairRap = fs->make<TH1F>("PairRap", "TnP Pair y", 100, -5, 5);

  PassingPairMass = fs->make<TH1F>("PassingPairMass", "TnP PassingPair mass", 20000, 0, 20000);
  PassingPairPt = fs->make<TH1F>("PassingPairPt", "TnP PassingPair pT", 10000, 0, 10000);
  PassingPairEta = fs->make<TH1F>("PassingPairEta", "TnP PassingPair #eta",    100, -5, 5);
  PassingPairRap = fs->make<TH1F>("PassingPairRap", "TnP PassingPair y", 100, -5, 5);

  FailingPairMass = fs->make<TH1F>("FailingPairMass", "TnP FailingPair mass", 20000, 0, 20000);
  FailingPairPt = fs->make<TH1F>("FailingPairPt", "TnP FailingPair pT", 10000, 0, 10000);
  FailingPairEta = fs->make<TH1F>("FailingPairEta", "TnP FailingPair #eta",    100, -5, 5);
  FailingPairRap = fs->make<TH1F>("FailingPairRap", "TnP FailingPair y", 100, -5, 5);

}

void Zprime2muHistosForExtraSegs_FromAOD::getBSandPV(const edm::Event& event) {
  // We store these as bare pointers. Should find better way, but
  // don't want to pass them around everywhere...
  edm::Handle<reco::BeamSpot> hbs;
  event.getByLabel(beamspot_src, hbs);
  beamspot = hbs.isValid() ? &*hbs : 0; // nice and fragile
  NBeamSpot->Fill(beamspot != 0);

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
  NVertices->Fill(vertex_count, _totalWeight );
}

std::pair<int, int> Zprime2muHistosForExtraSegs_FromAOD::getNShowers(const pat::CompositeCandidate& dil) {
  int lep0_nShowers = 0;
  if(dil.userInt("lep0_nSegs_1") >= nshowers_threshold_min)
    lep0_nShowers+=1;
  if(dil.userInt("lep0_nSegs_2") >= nshowers_threshold_min)
    lep0_nShowers+=1;
  if(dil.userInt("lep0_nSegs_3") >= nshowers_threshold_min)
    lep0_nShowers+=1;
  if(dil.userInt("lep0_nSegs_4") >= nshowers_threshold_min)
    lep0_nShowers+=1;

  int lep1_nShowers = 0;
  if(dil.userInt("lep1_nSegs_1") >= nshowers_threshold_min)
    lep1_nShowers+=1;
  if(dil.userInt("lep1_nSegs_2") >= nshowers_threshold_min)
    lep1_nShowers+=1;
  if(dil.userInt("lep1_nSegs_3") >= nshowers_threshold_min)
    lep1_nShowers+=1;
  if(dil.userInt("lep1_nSegs_4") >= nshowers_threshold_min)
    lep1_nShowers+=1;

  return std::make_pair(lep0_nShowers, lep1_nShowers);
}

void Zprime2muHistosForExtraSegs_FromAOD::fillTnPControlHistos(const pat::CompositeCandidate& dil,
                                                 const reco::CandidateBaseRef& TagMu,
                                                 const reco::CandidateBaseRef& ProbeMu,
                                                 int   probe_nHits,
                                                 int   probe_nSegs,
                                                 int   probe_nShowers,
                                                 bool  isPassNoPt ) {

  //for offline variables
  //const pat::Muon* TagPat = toConcretePtr<pat::Muon>(TagMu);
  //const pat::Muon* ProbePat = toConcretePtr<pat::Muon>(ProbeMu);

  TagPt->Fill( TagMu->pt(), _totalWeight );
  TagEta->Fill( TagMu->eta(), _totalWeight );
  TagPhi->Fill( TagMu->phi(), _totalWeight );

  /*if(TagPat) {  // fill offline
    const reco::MuonIsolation& iso = TagPat->isolationR03();
    TagAbsTkIso->Fill( iso.sumPt, _totalWeight );
    TagRelTkIso->Fill( iso.sumPt / TagPat->innerTrack()->pt(), _totalWeight );

    //const reco::TrackRef track = patmuon::getPickedTrack(*TagPat);
    const reco::TrackRef GlTrack = TagPat->globalTrack();
    const reco::TrackRef InTrack = TagPat->innerTrack();
    const reco::TrackRef BeTrack = TagPat->muonBestTrack();
    if (GlTrack.isAvailable() && InTrack.isAvailable() && BeTrack.isAvailable()) {
      TagChi2dof->Fill( GlTrack->normalizedChi2(), _totalWeight );

      if (beamspot != 0) {
        TagTrackDXYBS->Fill( fabs(BeTrack->dxy(beamspot->position())), _totalWeight );
        TagTrackDZBS->Fill( fabs(BeTrack->dz (beamspot->position())), _totalWeight );
      }

      if (vertex != 0) {
        TagTrackDXYPV->Fill( fabs(BeTrack->dxy(vertex->position())), _totalWeight );
        TagTrackDZPV->Fill( fabs(BeTrack->dz (vertex->position())), _totalWeight );
      }

      const reco::HitPattern& hp = InTrack->hitPattern();
      TagNPxHits->Fill( hp.numberOfValidPixelHits(), _totalWeight );
      TagNStHits->Fill( hp.numberOfValidStripHits(), _totalWeight );
      TagNTkHits->Fill( hp.numberOfValidTrackerHits(), _totalWeight );
      TagNMuHits->Fill( GlTrack->hitPattern().numberOfValidMuonHits(), _totalWeight );

      TagNHits->Fill( GlTrack->hitPattern().numberOfValidHits(), _totalWeight );

      TagNPxLayers->Fill( hp.pixelLayersWithMeasurement(), _totalWeight );
      TagNStLayers->Fill( hp.stripLayersWithMeasurement(), _totalWeight );
      TagNTkLayers->Fill( hp.trackerLayersWithMeasurement(), _totalWeight );
    }
  }*/

  ProbePt->Fill( ProbeMu->pt(), _totalWeight );
  if( isPassNoPt )
    PassingProbePt->Fill( ProbeMu->pt(), _totalWeight );
  else
    FailingProbePt->Fill( ProbeMu->pt(), _totalWeight );

  PairNoPtMass->Fill( dil.mass(), _totalWeight );
  PairNoPtPt->Fill( dil.pt(), _totalWeight );
  PairNoPtEta->Fill( dil.eta(), _totalWeight );
  PairNoPtRap->Fill( dil.rapidity(), _totalWeight );

  if(ProbeMu->pt() > probe_pt_min) {
    ProbeEta->Fill( ProbeMu->eta(), _totalWeight );
    ProbePhi->Fill( ProbeMu->phi(), _totalWeight );
    ProbeNVertices->Fill( nVtx, _totalWeight );
    ProbeEtaPhi->Fill( ProbeMu->eta(), ProbeMu->phi(), _totalWeight );

    if( fabs(ProbeMu->eta())<0.9 ) {
      ProbeNHitsBP->Fill( probe_nHits, ProbeMu->p(), _totalWeight );
      ProbeNSegsBP->Fill( probe_nSegs, ProbeMu->p(), _totalWeight );
      ProbeNShowerBP->Fill( probe_nShowers, ProbeMu->p(), _totalWeight );
    }
    else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.4 ) {
      ProbeNHitsEP->Fill( probe_nHits, ProbeMu->p(), _totalWeight );
      ProbeNSegsEP->Fill( probe_nSegs, ProbeMu->p(), _totalWeight );
      ProbeNShowerEP->Fill( probe_nShowers, ProbeMu->p(), _totalWeight );
    }

    PairMass->Fill( dil.mass(), _totalWeight );
    PairPt->Fill( dil.pt(), _totalWeight );
    PairEta->Fill( dil.eta(), _totalWeight );
    PairRap->Fill( dil.rapidity(), _totalWeight );

    if( isPassNoPt ) {
      PassingProbeEta->Fill( ProbeMu->eta(), _totalWeight );
      PassingProbePhi->Fill( ProbeMu->phi(), _totalWeight );
      PassingProbeNVertices->Fill( nVtx, _totalWeight );
      PassingProbeEtaPhi->Fill( ProbeMu->eta(), ProbeMu->phi(), _totalWeight );

      if( fabs(ProbeMu->eta())<0.9 ) {
        PassingProbeNHitsBP->Fill( probe_nHits, ProbeMu->p(), _totalWeight );
        PassingProbeNSegsBP->Fill( probe_nSegs, ProbeMu->p(), _totalWeight );
        PassingProbeNShowerBP->Fill( probe_nShowers, ProbeMu->p(), _totalWeight );
      }
      else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.4 ) {
        PassingProbeNHitsEP->Fill( probe_nHits, ProbeMu->p(), _totalWeight );
        PassingProbeNSegsEP->Fill( probe_nSegs, ProbeMu->p(), _totalWeight );
        PassingProbeNShowerEP->Fill( probe_nShowers, ProbeMu->p(), _totalWeight );
      }

      PassingPairMass->Fill( dil.mass(), _totalWeight );
      PassingPairPt->Fill( dil.pt(), _totalWeight );
      PassingPairEta->Fill( dil.eta(), _totalWeight );
      PassingPairRap->Fill( dil.rapidity(), _totalWeight );
    }
    else {
      FailingProbeEta->Fill( ProbeMu->eta(), _totalWeight );
      FailingProbePhi->Fill( ProbeMu->phi(), _totalWeight );
      FailingProbeNVertices->Fill( nVtx, _totalWeight );
      FailingProbeEtaPhi->Fill( ProbeMu->eta(), ProbeMu->phi(), _totalWeight );

      if( fabs(ProbeMu->eta())<0.9 ) {
        FailingProbeNHitsBP->Fill( probe_nHits, ProbeMu->p(), _totalWeight );
        FailingProbeNSegsBP->Fill( probe_nSegs, ProbeMu->p(), _totalWeight );
        FailingProbeNShowerBP->Fill( probe_nShowers, ProbeMu->p(), _totalWeight );
      }
      else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.4 ) {
        FailingProbeNHitsEP->Fill( probe_nHits, ProbeMu->p(), _totalWeight );
        FailingProbeNSegsEP->Fill( probe_nSegs, ProbeMu->p(), _totalWeight );
        FailingProbeNShowerEP->Fill( probe_nShowers, ProbeMu->p(), _totalWeight );
      }

      FailingPairMass->Fill( dil.mass(), _totalWeight );
      FailingPairPt->Fill( dil.pt(), _totalWeight );
      FailingPairEta->Fill( dil.eta(), _totalWeight );
      FailingPairRap->Fill( dil.rapidity(), _totalWeight );
    }

    /*if(ProbePat) {  // fill offline
      const reco::MuonIsolation& iso = ProbePat->isolationR03();
      ProbeAbsTkIso->Fill( iso.sumPt, _totalWeight );
      ProbeRelTkIso->Fill( iso.sumPt / ProbePat->innerTrack()->pt(), _totalWeight );

      //const reco::TrackRef track = patmuon::getPickedTrack(*ProbePat);
      const reco::TrackRef GlTrack = ProbePat->globalTrack();
      const reco::TrackRef InTrack = ProbePat->innerTrack();
      const reco::TrackRef BeTrack = ProbePat->muonBestTrack();
      if (GlTrack.isAvailable() && InTrack.isAvailable() && BeTrack.isAvailable()) {
        ProbeChi2dof->Fill( GlTrack->normalizedChi2(), _totalWeight );

        if (beamspot != 0) {
          ProbeTrackDXYBS->Fill( fabs(BeTrack->dxy(beamspot->position())), _totalWeight );
          ProbeTrackDZBS->Fill( fabs(BeTrack->dz (beamspot->position())), _totalWeight );
        }

        if (vertex != 0) {
          ProbeTrackDXYPV->Fill( fabs(BeTrack->dxy(vertex->position())), _totalWeight );
          ProbeTrackDZPV->Fill( fabs(BeTrack->dz (vertex->position())), _totalWeight );
        }

        const reco::HitPattern& hp = InTrack->hitPattern();
        ProbeNPxHits->Fill( hp.numberOfValidPixelHits(), _totalWeight );
        ProbeNStHits->Fill( hp.numberOfValidStripHits(), _totalWeight );
        ProbeNTkHits->Fill( hp.numberOfValidTrackerHits(), _totalWeight );
        ProbeNMuHits->Fill( GlTrack->hitPattern().numberOfValidMuonHits(), _totalWeight );

        ProbeNHits->Fill( GlTrack->hitPattern().numberOfValidHits(), _totalWeight );

        ProbeNPxLayers->Fill( hp.pixelLayersWithMeasurement(), _totalWeight );
        ProbeNStLayers->Fill( hp.stripLayersWithMeasurement(), _totalWeight );
        ProbeNTkLayers->Fill( hp.trackerLayersWithMeasurement(), _totalWeight );
      }
    }*/
  }
}


void Zprime2muHistosForExtraSegs_FromAOD::analyze(const edm::Event& event, const edm::EventSetup& setup) {

  //---- Prescales : Not using for now...
    //  edm::Handle<int> hltPrescale;
    //  edm::Handle<int> l1Prescale;
    //  event.getByLabel(edm::InputTag("getPrescales","HLTPrescale","Zprime2muAnalysis"), hltPrescale);
    //  event.getByLabel(edm::InputTag("getPrescales","L1Prescale","Zprime2muAnalysis"), l1Prescale);
    if (_usePrescaleWeight) {
      edm::Handle<int> totalPrescale;
      event.getByLabel(edm::InputTag("getPrescales","TotalPrescale","Zprime2muAnalysis"), totalPrescale);
      _prescaleWeight = *totalPrescale;
    }
    //  std::cout<<*hltPrescale<<std::endl;
    //  std::cout<<l1Prescale<<std::endl;
    //  std::cout<<totalPrescale<<std::endl;

  //---- Generator weights
  if (_useMadgraphWeight) {
    _eventWeight = 1.;
    _madgraphWeight = 1.;

    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid()){
      _eventWeight = gen_ev_info->weight();
      _madgraphWeight = ( _eventWeight > 0 ) ? 1.0 : -1.0;
    }
    WeightMadGraph->Fill( _madgraphWeight );
  }

  //---- Get PileUp Weights
  int thePU = -1;
  _pileupWeight = 1.;
  if( !event.isRealData() ){

    edm::Handle<std::vector< PileupSummaryInfo > > pileup;
    if ( event.getByLabel(pileup_src, pileup)){
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = pileup->begin(); PVI != pileup->end(); ++PVI)
      {
        if(PVI->getBunchCrossing()==0){
          thePU = PVI->getTrueNumInteractions();
          continue;
        }
      }

      int nPileUpBin = (int)vec_PileUpWeight.size();
      for(int iPU=0; iPU<nPileUpBin; ++iPU) {
        if(thePU == iPU) {
          _pileupWeight = vec_PileUpWeight[iPU];
        }
      }
    }
    else
      edm::LogError("") << "PU collection not found !!!";

    //if(!ShutUp)  std::cout << "Zprime2muHistosForExtraSegs_FromAOD::analyze : PU = " << thePU << "  PU weight = " << _pileupWeight << std::endl;
  }

  _totalWeight = 1.;
  if( _madgraphWeight != 1. || _pileupWeight != 1. )
    _totalWeight = _madgraphWeight * _pileupWeight;  // prescale weight is not considered yet
  //if(!ShutUp)  std::cout << "Zprime2muHistosForExtraSegs_FromAOD::analyze : Total weight = " << _totalWeight << std::endl;

  if (use_bs_and_pv)
    getBSandPV(event);

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);

  if( !dileptons.isValid() ) {
    std::cout << "Zprime2muHistosForExtraSegs_FromAOD::analyze : !dileptons.isValid() ---> return" << std::endl;
    return;
  }

  NDileptons->Fill(dileptons->size(), _totalWeight );

  pat::CompositeCandidateCollection::const_iterator dil = dileptons->begin(), dile = dileptons->end();
  for ( ; dil != dile; ++dil) {

    float lep0_dpt_over_pt = dil->userFloat("lep0_dpt_over_pt");
    float lep1_dpt_over_pt = dil->userFloat("lep1_dpt_over_pt");
    std::pair<int, int> pair_nShowers = getNShowers(*dil);

    //HERE
    int lep0_nHits = ( dil->userInt("lep0_nHits_1") + dil->userInt("lep0_nHits_2") + dil->userInt("lep0_nHits_3") + dil->userInt("lep0_nHits_4") );
    int lep0_nSegs = ( dil->userInt("lep0_nSegs_1") + dil->userInt("lep0_nSegs_2") + dil->userInt("lep0_nSegs_3") + dil->userInt("lep0_nSegs_4") );
    int lep0_nShowers = pair_nShowers.first;
    int lep1_nHits = ( dil->userInt("lep1_nHits_1") + dil->userInt("lep1_nHits_2") + dil->userInt("lep1_nHits_3") + dil->userInt("lep1_nHits_4") );
    int lep1_nSegs = ( dil->userInt("lep1_nSegs_1") + dil->userInt("lep1_nSegs_2") + dil->userInt("lep1_nSegs_3") + dil->userInt("lep1_nSegs_4") );
    int lep1_nShowers = pair_nShowers.second;

    const reco::CandidateBaseRef& lep0 = dileptonDaughter(*dil, 0);
    const reco::CandidateBaseRef& lep1 = dileptonDaughter(*dil, 1);

    if(lep0.isNonnull() && lep1.isNonnull()) {
      const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
      const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
      if(mu0 && mu1) {

        //---- Tag0 and Probe1
        if( tag_selector(*mu0) && (lep0_dpt_over_pt < tag_dpt_over_pt_max) && (fabs(mu0->muonBestTrack()->dz( vertex->position() )) < tag_dz_max) &&
            probe_selector(*mu1) && (lep1_dpt_over_pt < probe_dpt_over_pt_max) ) {
          const reco::CandidateBaseRef& TagMu = lep0;
          const reco::CandidateBaseRef& ProbeMu = lep1;

          int probe_nHits    = lep1_nHits;
          int probe_nSegs    = lep1_nSegs;
          int probe_nShowers = lep1_nShowers;

          if(!ShutUp)  std::cout << "Zprime2muHistosForExtraSegs_FromAOD::analyze : Tag0 and Probe1" << std::endl;
          if(!ShutUp)  std::cout << "                                              pT=" << ProbeMu->pt() << std::endl;
          if(!ShutUp)  std::cout << "                                             eta=" << ProbeMu->eta() << std::endl;
          if(!ShutUp)  std::cout << "                                             phi=" << ProbeMu->phi() << std::endl;
          if(!ShutUp)  std::cout << "                                            nVtx=" << nVtx << std::endl;

          if(!ShutUp)  std::cout << "                                            lep1_nSegs_1 =" << dil->userInt("lep1_nSegs_1") << std::endl;
          if(!ShutUp)  std::cout << "                                            lep1_nSegs_2 =" << dil->userInt("lep1_nSegs_2") << std::endl;
          if(!ShutUp)  std::cout << "                                            lep1_nSegs_3 =" << dil->userInt("lep1_nSegs_3") << std::endl;
          if(!ShutUp)  std::cout << "                                            lep1_nSegs_4 =" << dil->userInt("lep1_nSegs_4") << std::endl;
          if(!ShutUp)  std::cout << "                                            lep1_nSegs   =" << lep1_nSegs << std::endl;
          if(!ShutUp)  std::cout << "                                            lep1_nShowers=" << lep1_nShowers << std::endl;

          bool isPassingProbe = ( passing_probe_selector(*mu1) && (lep1_dpt_over_pt < passing_probe_dpt_over_pt_max) && (fabs(mu1->muonBestTrack()->dz( vertex->position() )) < passing_probe_dz_max) );
          if(!ShutUp)  std::cout << "                                  isPassingProbe = " << isPassingProbe << std::endl;

          fillTnPControlHistos(*dil, TagMu, ProbeMu, probe_nHits, probe_nSegs, probe_nShowers, isPassingProbe);

          if( ProbeMu->pt() > probe_pt_min ) {
            if( fabs(ProbeMu->eta())<0.9 ) {
              ProbeNSegsBP1St->Fill( dil->userInt("lep1_nSegs_1"), ProbeMu->p(), _totalWeight );
              ProbeNSegsBP2St->Fill( dil->userInt("lep1_nSegs_2"), ProbeMu->p(), _totalWeight );
              ProbeNSegsBP3St->Fill( dil->userInt("lep1_nSegs_3"), ProbeMu->p(), _totalWeight );
              ProbeNSegsBP4St->Fill( dil->userInt("lep1_nSegs_4"), ProbeMu->p(), _totalWeight );
            }
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.4 ) {
              ProbeNSegsEP1St->Fill( dil->userInt("lep1_nSegs_1"), ProbeMu->p(), _totalWeight );
              ProbeNSegsEP2St->Fill( dil->userInt("lep1_nSegs_2"), ProbeMu->p(), _totalWeight );
              ProbeNSegsEP3St->Fill( dil->userInt("lep1_nSegs_3"), ProbeMu->p(), _totalWeight );
              ProbeNSegsEP4St->Fill( dil->userInt("lep1_nSegs_4"), ProbeMu->p(), _totalWeight );
            }
          }

        } // Tag0 and Probe1

        //---- Tag1 and Probe0
        if( tag_selector(*mu1) && (lep1_dpt_over_pt < tag_dpt_over_pt_max) && (fabs(mu1->muonBestTrack()->dz( vertex->position() )) < tag_dz_max) &&
            probe_selector(*mu0) && (lep0_dpt_over_pt < probe_dpt_over_pt_max) ) {
          const reco::CandidateBaseRef& TagMu = lep1;
          const reco::CandidateBaseRef& ProbeMu = lep0;

          int probe_nHits    = lep0_nHits;
          int probe_nSegs    = lep0_nSegs;
          int probe_nShowers = lep0_nShowers;

          if(!ShutUp)  std::cout << "Zprime2muHistosForExtraSegs_FromAOD::analyze : Tag1 and Probe0" << std::endl;
          if(!ShutUp)  std::cout << "                                              pT=" << ProbeMu->pt() << std::endl;
          if(!ShutUp)  std::cout << "                                             eta=" << ProbeMu->eta() << std::endl;
          if(!ShutUp)  std::cout << "                                             phi=" << ProbeMu->phi() << std::endl;
          if(!ShutUp)  std::cout << "                                            nVtx=" << nVtx << std::endl;

          if(!ShutUp)  std::cout << "                                            lep0_nSegs_1 =" << dil->userInt("lep0_nSegs_1") << std::endl;
          if(!ShutUp)  std::cout << "                                            lep0_nSegs_2 =" << dil->userInt("lep0_nSegs_2") << std::endl;
          if(!ShutUp)  std::cout << "                                            lep0_nSegs_3 =" << dil->userInt("lep0_nSegs_3") << std::endl;
          if(!ShutUp)  std::cout << "                                            lep0_nSegs_4 =" << dil->userInt("lep0_nSegs_4") << std::endl;
          if(!ShutUp)  std::cout << "                                            lep0_nSegs   =" << lep0_nSegs << std::endl;
          if(!ShutUp)  std::cout << "                                            lep0_nShowers=" << lep0_nShowers << std::endl;

          bool isPassingProbe = ( passing_probe_selector(*mu0) && (lep0_dpt_over_pt < passing_probe_dpt_over_pt_max) && (fabs(mu0->muonBestTrack()->dz( vertex->position() )) < passing_probe_dz_max) );
          if(!ShutUp)  std::cout << "                                  isPassingProbe = " << isPassingProbe << std::endl;

          fillTnPControlHistos(*dil, TagMu, ProbeMu, probe_nHits, probe_nSegs, probe_nShowers, isPassingProbe);

          if( ProbeMu->pt() > probe_pt_min ) {
            if( fabs(ProbeMu->eta())<0.9 ) {
              ProbeNSegsBP1St->Fill( dil->userInt("lep0_nSegs_1"), ProbeMu->p(), _totalWeight );
              ProbeNSegsBP2St->Fill( dil->userInt("lep0_nSegs_2"), ProbeMu->p(), _totalWeight );
              ProbeNSegsBP3St->Fill( dil->userInt("lep0_nSegs_3"), ProbeMu->p(), _totalWeight );
              ProbeNSegsBP4St->Fill( dil->userInt("lep0_nSegs_4"), ProbeMu->p(), _totalWeight );
            }
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.4 ) {
              ProbeNSegsEP1St->Fill( dil->userInt("lep0_nSegs_1"), ProbeMu->p(), _totalWeight );
              ProbeNSegsEP2St->Fill( dil->userInt("lep0_nSegs_2"), ProbeMu->p(), _totalWeight );
              ProbeNSegsEP3St->Fill( dil->userInt("lep0_nSegs_3"), ProbeMu->p(), _totalWeight );
              ProbeNSegsEP4St->Fill( dil->userInt("lep0_nSegs_4"), ProbeMu->p(), _totalWeight );
            }
          }

        } // Tag1 and Probe0

        if(!ShutUp)  std::cout << std::endl;
      } // if(mu0 && mu1)
    } // if(lep0.isNonnull() && lep1.isNonnull())


  } // for ( ; dil != dile; ++dil)

}

DEFINE_FWK_MODULE(Zprime2muHistosForExtraSegs_FromAOD);
