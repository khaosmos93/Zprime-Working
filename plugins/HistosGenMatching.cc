#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"

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

class HistosGenMatching : public edm::EDAnalyzer {
 public:
  explicit HistosGenMatching(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void getBSandPV(const edm::Event&);
  std::pair<const pat::Muon*, const pat::Muon*> getGenMatchedPair(HardInteraction *, const edm::Handle<pat::MuonCollection> &);

  edm::InputTag muon_src;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;
  const bool use_bs_and_pv;
  const reco::BeamSpot* beamspot;
  const reco::Vertex*   vertex;
  int                   nVtx;

  HardInteraction* hardInteraction;

  StringCutObjectSelector<pat::Muon> muon_selector_DEN;
  StringCutObjectSelector<pat::Muon> muon_selector_NUM;

  double muon_dpt_over_pt_max_DEN;
  double muon_dpt_over_pt_max_NUM;

  double muon_dz_max_DEN;
  double muon_dz_max_NUM;

  double _eventWeight;
  bool   _useGenWeight;
  double _genWeight;

  bool _usePileupWeight;
  edm::InputTag pileup_src;
  double _pileupWeight;
  std::vector<double> vec_PileUpWeight;

  double _totalWeight;

  const bool ShutUp;

  // -- Histograms -- //
  TH1F* NBeamSpot;
  TH1F* NVertices;
  TH1F* GenWeight;

  typedef std::pair<TH1F*, TH1F*> eff_th1;
  typedef std::pair<TH2F*, TH2F*> eff_th2;

  std::pair<TH1F*, TH1F*> make_th1_pair(TString name, double nbins, double minX, double maxX) {
    edm::Service<TFileService> fs;
    TH1F* a = fs->make<TH1F>(name+"_DEN", "", nbins, minX, maxX);  a->Sumw2();
    TH1F* b = fs->make<TH1F>(name+"_NUM", "", nbins, minX, maxX);  b->Sumw2();
    return std::make_pair(a, b);
  }

  std::pair<TH2F*, TH2F*> make_th2_pair(TString name, double nbinsX, double minX, double maxX, double nbinsY, double minY, double maxY) {
    edm::Service<TFileService> fs;
    TH2F* a = fs->make<TH2F>(name+"_DEN", "", nbinsX, minX, maxX, nbinsY, minY, maxY);  a->Sumw2();
    TH2F* b = fs->make<TH2F>(name+"_NUM", "", nbinsX, minX, maxX, nbinsY, minY, maxY);  b->Sumw2();
    return std::make_pair(a, b);
  }

  eff_th1 Mass;
  eff_th1 Mass_BB;
  eff_th1 Mass_BE;

  eff_th2 EtaPt;
  eff_th2 EtaAbsP;
  eff_th2 EtaPhi;
};

HistosGenMatching::HistosGenMatching(const edm::ParameterSet& cfg)
  : muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    use_bs_and_pv(cfg.getParameter<bool>("use_bs_and_pv")),
    beamspot(0),
    vertex(0),
    nVtx(0),

    hardInteraction(new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction"))),

    muon_selector_DEN(cfg.getParameter<std::string>("muon_cut_DEN")),
    muon_selector_NUM(cfg.getParameter<std::string>("muon_cut_NUM")),

    muon_dpt_over_pt_max_DEN(cfg.getParameter<double>("muon_dpt_over_pt_max_DEN")),
    muon_dpt_over_pt_max_NUM(cfg.getParameter<double>("muon_dpt_over_pt_max_NUM")),

    muon_dz_max_DEN(cfg.getParameter<double>("muon_dz_max_DEN")),
    muon_dz_max_NUM(cfg.getParameter<double>("muon_dz_max_NUM")),

    _useGenWeight(cfg.getParameter<bool>("useGenWeight")),
    _genWeight(1.),

    _usePileupWeight(cfg.getParameter<bool>("usePileupWeight")),
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
  consumes<std::vector<reco::GenParticle>>(hardInteraction->src);


  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  NBeamSpot = fs->make<TH1F>("NBeamSpot", "# beamspots/event",  2, 0,  2);
  NVertices = fs->make<TH1F>("NVertices", "# vertices/event",  200, 0, 200);

  // Generator weights
  GenWeight = fs->make<TH1F>("GenWeight", "weight per event", 4, -2,2);


  Mass    = make_th1_pair("Mass", 10000, 0, 10000);
  Mass_BB = make_th1_pair("Mass_BB", 10000, 0, 10000);
  Mass_BE = make_th1_pair("Mass_BE", 10000, 0, 10000);

  EtaPt     = make_th2_pair("EtaPt",   96, -4.8, 4.8, 5000, 0, 5000);
  EtaAbsP   = make_th2_pair("EtaAbsP", 96, -4.8, 4.8, 5000, 0, 5000);
  EtaPhi    = make_th2_pair("EtaPhi",  96, -4.8, 4.8, 100, -TMath::Pi(), TMath::Pi());
}

void HistosGenMatching::getBSandPV(const edm::Event& event) {
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

std::pair<const pat::Muon*, const pat::Muon*> HistosGenMatching::getGenMatchedPair(HardInteraction *hardInt, const edm::Handle<pat::MuonCollection> &muons) {

  double minDR = 0.2;

  if(hardInt->lepPlus == 0 || hardInt->lepMinus == 0) {
    edm::LogWarning("HistosGenMatching") << "hardInt->lepPlus == 0 || hardInt->lepMinus == 0";
    return std::make_pair(nullptr, nullptr);
  }

  const pat::Muon *matchedPlus = nullptr;
  const pat::Muon *matchedMinus = nullptr;

  double plus_dr_max  = minDR*2.;
  double minus_dr_max = minDR*2.;

  for(pat::MuonCollection::const_iterator mu = muons->begin(), end = muons->end(); mu != end; ++mu) {
    double plus_dr_temp = reco::deltaR(*mu, *(hardInt->lepPlus) );
    double minus_dr_temp = reco::deltaR(*mu, *(hardInt->lepMinus) );

    if(plus_dr_temp < plus_dr_max) {
      plus_dr_max = plus_dr_temp;
      matchedPlus = &*mu;
    }

    if(minus_dr_temp < minus_dr_max) {
      minus_dr_max = minus_dr_temp;
      matchedMinus = &*mu;
    }

  }

  if( !( (plus_dr_max < minDR) && (minus_dr_max < minDR) ) )
    return std::make_pair(nullptr, nullptr);

  return std::make_pair( matchedPlus, matchedMinus );
}


void HistosGenMatching::analyze(const edm::Event& event, const edm::EventSetup& setup) {

  static const double muon_mass = 0.1056583;

  //--- Generator weights
  if (_useGenWeight) {
    _eventWeight = 1.;
    _genWeight = 1.;

    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid()){
      _eventWeight = gen_ev_info->weight();
      _genWeight = ( _eventWeight > 0 ) ? 1.0 : -1.0;
    }
    GenWeight->Fill( _genWeight );
  }

  //--- Get PileUp Weights
  if(_usePileupWeight) {
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

      if(!ShutUp)  std::cout << "HistosGenMatching::analyze : PU = " << thePU << "  PU weight = " << _pileupWeight << std::endl;
    }
  }

  _totalWeight = 1.;
  if( _genWeight != 1. || _pileupWeight != 1. )
    _totalWeight = _genWeight * _pileupWeight;
  if(!ShutUp)  std::cout << "HistosGenMatching::analyze : Total weight = " << _totalWeight << std::endl;

  if (use_bs_and_pv)
    getBSandPV(event);

  //--- HardInteraction
  hardInteraction->Fill(event);
  if (!hardInteraction->IsValid() ) {
    edm::LogWarning("HistosGenMatching") << "!hardInteraction->IsValid()";
    return;
  }

  edm::Handle<pat::MuonCollection> muons;
  event.getByLabel(muon_src, muons);

  if( !muons.isValid() ) {
    std::cout << "HistosGenMatching::analyze : !muons.isValid() ---> return" << std::endl;
    return;
  }

  std::pair<const pat::Muon*, const pat::Muon*> pair = getGenMatchedPair( hardInteraction, muons );

  if( pair.first == nullptr || pair.second == nullptr )
    return;


  const reco::TrackRef InnTrack_first = pair.first->innerTrack();
  const reco::TrackRef GloTrack_first = pair.first->globalTrack();
  const reco::TrackRef TunTrack_first = pair.first->tunePMuonBestTrack();
  if( !( InnTrack_first.isAvailable() && GloTrack_first.isAvailable() && TunTrack_first.isAvailable() ) )
    return;

  bool is_dpt_over_pt_DEN_first = ( (TunTrack_first->ptError() / TunTrack_first->pt()) < muon_dpt_over_pt_max_DEN );
  bool is_dz_DEN_first          = ( fabs(InnTrack_first->dz( vertex->position() )) < muon_dz_max_DEN );

  bool is_dpt_over_pt_NUM_first = ( (TunTrack_first->ptError() / TunTrack_first->pt()) < muon_dpt_over_pt_max_NUM );
  bool is_dz_NUM_first          = ( fabs(InnTrack_first->dz( vertex->position() )) < muon_dz_max_NUM );

  bool isDEN_first = ( (bool)muon_selector_DEN(*pair.first) &&
                       is_dpt_over_pt_DEN_first &&
                       is_dz_DEN_first
                      );

  bool isNUM_first = ( isDEN_first &&
                       (bool)muon_selector_NUM(*pair.first) &&
                       is_dpt_over_pt_NUM_first &&
                       is_dz_NUM_first
                      );


  const reco::TrackRef InnTrack_second = pair.second->innerTrack();
  const reco::TrackRef GloTrack_second = pair.second->globalTrack();
  const reco::TrackRef TunTrack_second = pair.second->tunePMuonBestTrack();
  if( !( InnTrack_second.isAvailable() && GloTrack_second.isAvailable() && TunTrack_second.isAvailable() ) )
    return;

  bool is_dpt_over_pt_DEN_second = ( (TunTrack_second->ptError() / TunTrack_second->pt()) < muon_dpt_over_pt_max_DEN );
  bool is_dz_DEN_second          = ( fabs(InnTrack_second->dz( vertex->position() )) < muon_dz_max_DEN );

  bool is_dpt_over_pt_NUM_second = ( (TunTrack_second->ptError() / TunTrack_second->pt()) < muon_dpt_over_pt_max_NUM );
  bool is_dz_NUM_second          = ( fabs(InnTrack_second->dz( vertex->position() )) < muon_dz_max_NUM );

  bool isDEN_second = ( (bool)muon_selector_DEN(*pair.second) &&
                        is_dpt_over_pt_DEN_second &&
                        is_dz_DEN_second
                       );

  bool isNUM_second = ( isDEN_second &&
                        (bool)muon_selector_NUM(*pair.second) &&
                        is_dpt_over_pt_NUM_second &&
                        is_dz_NUM_second
                       );

  TLorentzVector mu_first, mu_second;
  mu_first.SetPtEtaPhiM(TunTrack_first->pt(), TunTrack_first->eta(), TunTrack_first->phi(), muon_mass);
  mu_second.SetPtEtaPhiM(TunTrack_second->pt(), TunTrack_second->eta(), TunTrack_second->phi(), muon_mass);
  double mass = (mu_first+mu_second).M();

  //--- Fill
  if( isDEN_first && isDEN_second ) {

    EtaPt.first->Fill( TunTrack_first->eta(), TunTrack_first->pt(), _totalWeight );
    EtaAbsP.first->Fill( TunTrack_first->eta(), TunTrack_first->p(), _totalWeight );
    EtaPhi.first->Fill( TunTrack_first->eta(), TunTrack_first->phi(), _totalWeight );

    EtaPt.first->Fill( TunTrack_second->eta(), TunTrack_second->pt(), _totalWeight );
    EtaAbsP.first->Fill( TunTrack_second->eta(), TunTrack_second->p(), _totalWeight );
    EtaPhi.first->Fill( TunTrack_second->eta(), TunTrack_second->phi(), _totalWeight );

    if( isNUM_first ) {
      EtaPt.second->Fill( TunTrack_first->eta(), TunTrack_first->pt(), _totalWeight );
      EtaAbsP.second->Fill( TunTrack_first->eta(), TunTrack_first->p(), _totalWeight );
      EtaPhi.second->Fill( TunTrack_first->eta(), TunTrack_first->phi(), _totalWeight );
    }

    if( isNUM_second ) {
      EtaPt.second->Fill( TunTrack_second->eta(), TunTrack_second->pt(), _totalWeight );
      EtaAbsP.second->Fill( TunTrack_second->eta(), TunTrack_second->p(), _totalWeight );
      EtaPhi.second->Fill( TunTrack_second->eta(), TunTrack_second->phi(), _totalWeight );
    }

  }


  bool isBB =  ( fabs(TunTrack_first->eta())<0.9 && fabs(TunTrack_second->eta())<0.9 );

  if( isDEN_first && isDEN_second ) {

    Mass.first->Fill( mass, _totalWeight );
    if( isBB )
      Mass_BB.first->Fill( mass, _totalWeight );
    else
      Mass_BE.first->Fill( mass, _totalWeight );

    if( isNUM_first && isNUM_second ) {
      Mass.second->Fill( mass, _totalWeight );
      if( isBB )
        Mass_BB.second->Fill( mass, _totalWeight );
      else
        Mass_BE.second->Fill( mass, _totalWeight );
    }

  }

}

DEFINE_FWK_MODULE(HistosGenMatching);
