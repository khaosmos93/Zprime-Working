#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TH1F.h"
#include "TTree.h"
#include "TMath.h"


class GenPileup : public edm::EDFilter {
 public:
  explicit GenPileup(const edm::ParameterSet&);
  ~GenPileup();

 private:
 virtual void beginJob();
 virtual bool filter(edm::Event&, const edm::EventSetup&);
 virtual void endJob();

  edm::InputTag pileup_src;
  edm::InputTag vertex_src;

  double eventWeight;
  double eventWeightSign;

  double SumOfWeight;
  double SumOfWeightSign;

  TH1F *GenWeightSign;
  TH1F *GenTruePileup;
  TH1F *GenTruePileup100;
  TH1F *OffNVTX;
};


GenPileup::GenPileup(const edm::ParameterSet& cfg) 
 :  pileup_src(cfg.getParameter<edm::InputTag>("pileup_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src"))
{
  consumes<reco::VertexCollection>(vertex_src);
  consumes<std::vector<PileupSummaryInfo> >(pileup_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));
}

GenPileup::~GenPileup(){}

bool GenPileup::filter(edm::Event& event, const edm::EventSetup&) {

  eventWeight = -99999.;
  eventWeightSign = -99999.;

  if(! event.isRealData() ) {

    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid() ){
      eventWeight = gen_ev_info->weight();
      eventWeightSign = ( eventWeight > 0 ) ? 1.0 : -1.0;

      SumOfWeight += eventWeight;
      SumOfWeightSign += eventWeightSign;
      GenWeightSign->Fill(eventWeightSign);

      //--- True PU
      double truePU = -999.;
      edm::Handle<std::vector< PileupSummaryInfo > > pileup;
      if ( event.getByLabel(pileup_src, pileup) ) {
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for(PVI = pileup->begin(); PVI != pileup->end(); ++PVI) {
          if(PVI->getBunchCrossing()==0){
            truePU = PVI->getTrueNumInteractions();
            continue;
          }
        }
        GenTruePileup->Fill(truePU, eventWeightSign);
        GenTruePileup100->Fill(truePU, eventWeightSign);
      }
      else
        edm::LogError("") << "! event.getByLabel(pileup_src, pileup)";

      //--- Offline VTX
      double nVTX = 0;
      edm::Handle<reco::VertexCollection> vertices;
      if( event.getByLabel(vertex_src, vertices) ) {
        for (reco::VertexCollection::const_iterator it = vertices->begin(), ite = vertices->end(); it != ite; ++it) {
          if (it->ndof() > 4 && fabs(it->z()) <= 24 && fabs(it->position().rho()) <= 2)
            nVTX += 1.;
        }
        OffNVTX->Fill(nVTX, eventWeightSign);
      }
      else
        edm::LogError("") << "! event.getByLabel(vertex_src, vertices)";

    }
  }

  return
    (eventWeight!=0.0);
}

void GenPileup::beginJob() {
  eventWeight = 0.;
  eventWeightSign = 0.;
  SumOfWeight = 0.;
  SumOfWeightSign = 0.;

  edm::Service<TFileService> fs;
  GenWeightSign = fs->make<TH1F>("GenWeightSign", "", 5, -2, 3);
  GenTruePileup = fs->make<TH1F>("GenTruePileup", "", 301, -1, 300);
  GenTruePileup100 = fs->make<TH1F>("GenTruePileup100", "", 100, 0, 100);
  OffNVTX = fs->make<TH1F>("OffNVTX", "", 500, 0, 500);
}

void GenPileup::endJob() {
  /*std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "    SumOfWeight : " << SumOfWeight << std::endl;
  std::cout << "SumOfWeightSign : " << SumOfWeightSign << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;*/

  //GenTruePileup->Scale(1./SumOfWeightSign);
  //GenTruePileup100->Scale(1./SumOfWeightSign);
}

DEFINE_FWK_MODULE(GenPileup);
