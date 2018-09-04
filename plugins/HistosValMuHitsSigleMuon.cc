#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

//#include <boost/foreach.hpp>

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

//using namespace std;

class HistosValMuHitsSigleMuon : public edm::EDAnalyzer {
 public:
  explicit HistosValMuHitsSigleMuon(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void getBSandPV(const edm::Event&);

  edm::InputTag muon_src;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;
  const bool use_bs_and_pv;
  const reco::BeamSpot* beamspot;
  const reco::Vertex*   vertex;
  int                   nVtx;

  StringCutObjectSelector<pat::Muon> muon_selector0;
  StringCutObjectSelector<pat::Muon> muon_selector1;
  StringCutObjectSelector<pat::Muon> muon_selector2;

  double muon_dpt_over_pt_max;
  double muon_dz_max;

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
  TH1F* WeightMadGraph;

  TH1F *Muon0_Pt;
  TH1F *Muon0_PtB;
  TH1F *Muon0_PtO;
  TH1F *Muon0_PtE;
  TH1F *Muon0_PtF;

  TH1F *Muon0_AbsP;
  TH1F *Muon0_AbsPB;
  TH1F *Muon0_AbsPO;
  TH1F *Muon0_AbsPE;
  TH1F *Muon0_AbsPF;

  TH1F *Muon0_Eta;
  TH1F *Muon0_Phi;

  TH2F *Muon0_PtEta;

  TH1F *Muon1_Pt;
  TH1F *Muon1_PtB;
  TH1F *Muon1_PtO;
  TH1F *Muon1_PtE;
  TH1F *Muon1_PtF;

  TH1F *Muon1_AbsP;
  TH1F *Muon1_AbsPB;
  TH1F *Muon1_AbsPO;
  TH1F *Muon1_AbsPE;
  TH1F *Muon1_AbsPF;

  TH1F *Muon1_Eta;
  TH1F *Muon1_Phi;

  TH2F *Muon1_PtEta;

  TH1F *Muon2_Pt;
  TH1F *Muon2_PtB;
  TH1F *Muon2_PtO;
  TH1F *Muon2_PtE;
  TH1F *Muon2_PtF;

  TH1F *Muon2_AbsP;
  TH1F *Muon2_AbsPB;
  TH1F *Muon2_AbsPO;
  TH1F *Muon2_AbsPE;
  TH1F *Muon2_AbsPF;

  TH1F *Muon2_Eta;
  TH1F *Muon2_Phi;

  TH2F *Muon2_PtEta;
};

HistosValMuHitsSigleMuon::HistosValMuHitsSigleMuon(const edm::ParameterSet& cfg)
  : muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    use_bs_and_pv(cfg.getParameter<bool>("use_bs_and_pv")),
    beamspot(0),
    vertex(0),
    nVtx(0),

    muon_selector0(cfg.getParameter<std::string>("muon_cut0")),
    muon_selector1(cfg.getParameter<std::string>("muon_cut1")),
    muon_selector2(cfg.getParameter<std::string>("muon_cut2")),

    muon_dpt_over_pt_max(cfg.getParameter<double>("muon_dpt_over_pt_max")),
    muon_dz_max(cfg.getParameter<double>("muon_dz_max")),

    _useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    _madgraphWeight(1.),

    pileup_src(cfg.getParameter<edm::InputTag>("pileup_src")),
    _pileupWeight(1.),
    vec_PileUpWeight(cfg.getParameter<std::vector<double>>("vec_PileUpWeight")),

    _totalWeight(1.),

    ShutUp(cfg.getParameter<bool>("ShutUp"))
{
  consumes<pat::MuonCollection>(muon_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertex_src);
  consumes<std::vector<PileupSummaryInfo> >(pileup_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  NBeamSpot = fs->make<TH1F>("NBeamSpot", "# beamspots/event",  2, 0,  2);
  NVertices = fs->make<TH1F>("NVertices", "# vertices/event",  200, 0, 200);

  // Generator weights
  WeightMadGraph = fs->make<TH1F>("WeightMadGraph", "weight per event", 4, -2,2);


  //--- Muon0
  Muon0_Pt = fs->make<TH1F>("Muon0_Pt", "", 10000, 0, 10000);
  Muon0_PtB = fs->make<TH1F>("Muon0_PtB", "", 10000, 0, 10000);
  Muon0_PtO = fs->make<TH1F>("Muon0_PtO", "", 10000, 0, 10000);
  Muon0_PtE = fs->make<TH1F>("Muon0_PtE", "", 10000, 0, 10000);
  Muon0_PtF = fs->make<TH1F>("Muon0_PtF", "", 10000, 0, 10000);

  Muon0_AbsP = fs->make<TH1F>("Muon0_AbsP", "", 10000, 0, 10000);
  Muon0_AbsPB = fs->make<TH1F>("Muon0_AbsPB", "", 10000, 0, 10000);
  Muon0_AbsPO = fs->make<TH1F>("Muon0_AbsPO", "", 10000, 0, 10000);
  Muon0_AbsPE = fs->make<TH1F>("Muon0_AbsPE", "", 10000, 0, 10000);
  Muon0_AbsPF = fs->make<TH1F>("Muon0_AbsPF", "", 10000, 0, 10000);

  Muon0_Eta = fs->make<TH1F>("Muon0_Eta", "", 100, -5, 5);
  Muon0_Phi = fs->make<TH1F>("Muon0_Phi", "", 100, -TMath::Pi(), TMath::Pi());

  Muon0_PtEta = fs->make<TH2F>("Muon0_PtEta", "", 10000, 0, 10000, 100, -5, 5);

  //--- Muon1
  Muon1_Pt = fs->make<TH1F>("Muon1_Pt", "", 10000, 0, 10000);
  Muon1_PtB = fs->make<TH1F>("Muon1_PtB", "", 10000, 0, 10000);
  Muon1_PtO = fs->make<TH1F>("Muon1_PtO", "", 10000, 0, 10000);
  Muon1_PtE = fs->make<TH1F>("Muon1_PtE", "", 10000, 0, 10000);
  Muon1_PtF = fs->make<TH1F>("Muon1_PtF", "", 10000, 0, 10000);

  Muon1_AbsP = fs->make<TH1F>("Muon1_AbsP", "", 10000, 0, 10000);
  Muon1_AbsPB = fs->make<TH1F>("Muon1_AbsPB", "", 10000, 0, 10000);
  Muon1_AbsPO = fs->make<TH1F>("Muon1_AbsPO", "", 10000, 0, 10000);
  Muon1_AbsPE = fs->make<TH1F>("Muon1_AbsPE", "", 10000, 0, 10000);
  Muon1_AbsPF = fs->make<TH1F>("Muon1_AbsPF", "", 10000, 0, 10000);

  Muon1_Eta = fs->make<TH1F>("Muon1_Eta", "", 100, -5, 5);
  Muon1_Phi = fs->make<TH1F>("Muon1_Phi", "", 100, -TMath::Pi(), TMath::Pi());

  Muon1_PtEta = fs->make<TH2F>("Muon1_PtEta", "", 10000, 0, 10000, 100, -5, 5);

  // Muon2
  Muon2_Pt = fs->make<TH1F>("Muon2_Pt", "", 10000, 0, 10000);
  Muon2_PtB = fs->make<TH1F>("Muon2_PtB", "", 10000, 0, 10000);
  Muon2_PtO = fs->make<TH1F>("Muon2_PtO", "", 10000, 0, 10000);
  Muon2_PtE = fs->make<TH1F>("Muon2_PtE", "", 10000, 0, 10000);
  Muon2_PtF = fs->make<TH1F>("Muon2_PtF", "", 10000, 0, 10000);

  Muon2_AbsP = fs->make<TH1F>("Muon2_AbsP", "", 10000, 0, 10000);
  Muon2_AbsPB = fs->make<TH1F>("Muon2_AbsPB", "", 10000, 0, 10000);
  Muon2_AbsPO = fs->make<TH1F>("Muon2_AbsPO", "", 10000, 0, 10000);
  Muon2_AbsPE = fs->make<TH1F>("Muon2_AbsPE", "", 10000, 0, 10000);
  Muon2_AbsPF = fs->make<TH1F>("Muon2_AbsPF", "", 10000, 0, 10000);

  Muon2_Eta = fs->make<TH1F>("Muon2_Eta", "", 100, -5, 5);
  Muon2_Phi = fs->make<TH1F>("Muon2_Phi", "", 100, -TMath::Pi(), TMath::Pi());

  Muon2_PtEta = fs->make<TH2F>("Muon2_PtEta", "", 10000, 0, 10000, 100, -5, 5);
}

void HistosValMuHitsSigleMuon::getBSandPV(const edm::Event& event) {
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
  if(vertex == 0)
    vertex = &*vertices->begin();

  nVtx = vertex_count;
  NVertices->Fill(vertex_count, _totalWeight );
}

void HistosValMuHitsSigleMuon::analyze(const edm::Event& event, const edm::EventSetup& setup) {

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

    if(!ShutUp)  std::cout << "HistosValMuHitsSigleMuon::analyze : PU = " << thePU << "  PU weight = " << _pileupWeight << std::endl;
  }

  _totalWeight = 1.;
  if( _madgraphWeight != 1. || _pileupWeight != 1. )
    _totalWeight = _madgraphWeight * _pileupWeight;
  if(!ShutUp)  std::cout << "HistosValMuHitsSigleMuon::analyze : Total weight = " << _totalWeight << std::endl;

  if (use_bs_and_pv)
    getBSandPV(event);

  edm::Handle<pat::MuonCollection> muons;
  event.getByLabel(muon_src, muons);

  if( !muons.isValid() ) {
    std::cout << "HistosValMuHitsSigleMuon::analyze : !muons.isValid() ---> return" << std::endl;
    return;
  }

  //BOOST_FOREACH(const pat::Muon& mu, *muons) {
  for(pat::MuonCollection::const_iterator mu = muons->begin(), end = muons->end(); mu != end; ++mu) {

    const reco::TrackRef InnTrack = mu->innerTrack();
    const reco::TrackRef GloTrack = mu->globalTrack();
    const reco::TrackRef TunTrack = mu->tunePMuonBestTrack();
    if( !( InnTrack.isAvailable() && GloTrack.isAvailable() && TunTrack.isAvailable() ) )
      continue;

    bool is_dpt_over_pt = ( (TunTrack->ptError() / TunTrack->pt()) < muon_dpt_over_pt_max ); // ? true : false;
    bool is_dz          = ( fabs(InnTrack->dz( vertex->position() )) < muon_dz_max ); // ? true : false;

    if( !(is_dpt_over_pt && is_dz) )
      continue;

    bool is_Muon0 = (bool)muon_selector0(*mu);
    bool is_Muon1 = (bool)muon_selector1(*mu);
    bool is_Muon2 = (bool)muon_selector2(*mu);

    bool isB = ( fabs(mu->eta())<0.9 ); // ? true : false;
    bool isO = ( fabs(mu->eta())>=0.9 && fabs(mu->eta())<1.2 ); // ? true : false;
    bool isE = ( fabs(mu->eta())>=1.2 && fabs(mu->eta())<2.1 ); // ? true : false;
    bool isF = ( fabs(mu->eta())>=2.1 && fabs(mu->eta())<2.4 ); // ? true : false;

    if( is_Muon0 ) {
      Muon0_Pt->Fill( TunTrack->pt(), _totalWeight );
      Muon0_AbsP->Fill( TunTrack->p(), _totalWeight );
      if( isB ) {
        Muon0_PtB->Fill( TunTrack->pt(), _totalWeight );
        Muon0_AbsPB->Fill( TunTrack->p(), _totalWeight );
      }
      if( isO ) {
        Muon0_PtO->Fill( TunTrack->pt(), _totalWeight );
        Muon0_AbsPO->Fill( TunTrack->p(), _totalWeight );
      }
      if( isE ) {
        Muon0_PtE->Fill( TunTrack->pt(), _totalWeight );
        Muon0_AbsPE->Fill( TunTrack->p(), _totalWeight );
      }
      if( isF ) {
        Muon0_PtF->Fill( TunTrack->pt(), _totalWeight );
        Muon0_AbsPF->Fill( TunTrack->p(), _totalWeight );
      }
      Muon0_Eta->Fill( TunTrack->eta(), _totalWeight );
      Muon0_Phi->Fill( TunTrack->phi(), _totalWeight );
      Muon0_PtEta->Fill( TunTrack->pt(), TunTrack->eta(), _totalWeight );
    }

    if( is_Muon1 ) {
      Muon1_Pt->Fill( TunTrack->pt(), _totalWeight );
      Muon1_AbsP->Fill( TunTrack->p(), _totalWeight );
      if( isB ) {
        Muon1_PtB->Fill( TunTrack->pt(), _totalWeight );
        Muon1_AbsPB->Fill( TunTrack->p(), _totalWeight );
      }
      if( isO ) {
        Muon1_PtO->Fill( TunTrack->pt(), _totalWeight );
        Muon1_AbsPO->Fill( TunTrack->p(), _totalWeight );
      }
      if( isE ) {
        Muon1_PtE->Fill( TunTrack->pt(), _totalWeight );
        Muon1_AbsPE->Fill( TunTrack->p(), _totalWeight );
      }
      if( isF ) {
        Muon1_PtF->Fill( TunTrack->pt(), _totalWeight );
        Muon1_AbsPF->Fill( TunTrack->p(), _totalWeight );
      }
      Muon1_Eta->Fill( TunTrack->eta(), _totalWeight );
      Muon1_Phi->Fill( TunTrack->phi(), _totalWeight );
      Muon1_PtEta->Fill( TunTrack->pt(), TunTrack->eta(), _totalWeight );
    }

    if( is_Muon2 ) {
      Muon2_Pt->Fill( TunTrack->pt(), _totalWeight );
      Muon2_AbsP->Fill( TunTrack->p(), _totalWeight );
      if( isB ) {
        Muon2_PtB->Fill( TunTrack->pt(), _totalWeight );
        Muon2_AbsPB->Fill( TunTrack->p(), _totalWeight );
      }
      if( isO ) {
        Muon2_PtO->Fill( TunTrack->pt(), _totalWeight );
        Muon2_AbsPO->Fill( TunTrack->p(), _totalWeight );
      }
      if( isE ) {
        Muon2_PtE->Fill( TunTrack->pt(), _totalWeight );
        Muon2_AbsPE->Fill( TunTrack->p(), _totalWeight );
      }
      if( isF ) {
        Muon2_PtF->Fill( TunTrack->pt(), _totalWeight );
        Muon2_AbsPF->Fill( TunTrack->p(), _totalWeight );
      }
      Muon2_Eta->Fill( TunTrack->eta(), _totalWeight );
      Muon2_Phi->Fill( TunTrack->phi(), _totalWeight );
      Muon2_PtEta->Fill( TunTrack->pt(), TunTrack->eta(), _totalWeight );
    }

  }

}

DEFINE_FWK_MODULE(HistosValMuHitsSigleMuon);
