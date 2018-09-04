#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TH1F.h"
#include "TTree.h"
#include "TMath.h"


class GenWeight : public edm::EDFilter {
 public:
  explicit GenWeight(const edm::ParameterSet&);
  ~GenWeight();

 private:
 virtual void beginJob() ;
 virtual bool filter(edm::Event&, const edm::EventSetup&);
 virtual void endJob() ;

  double eventWeight;
  double eventWeightSign;

  double SumOfWeight;
  double SumOfWeightSign;

  TTree* tree;
  TH1F *GenWeightSign;
};


GenWeight::GenWeight(const edm::ParameterSet& cfg) {
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));
}

GenWeight::~GenWeight(){}

bool GenWeight::filter(edm::Event& event, const edm::EventSetup&) {

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
      tree->Fill();
    }
  }

  return
    (eventWeight!=0.0);
}

void GenWeight::beginJob() {
  eventWeight = 0.;
  eventWeightSign = 0.;
  SumOfWeight = 0.;
  SumOfWeightSign = 0.;

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("eventWeight",&eventWeight,"eventWeight/D");
  GenWeightSign = fs->make<TH1F>("GenWeightSign", "", 5, -2, 3);
}

void GenWeight::endJob() {
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "    SumOfWeight : " << SumOfWeight << std::endl;
  std::cout << "SumOfWeightSign : " << SumOfWeightSign << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

DEFINE_FWK_MODULE(GenWeight);
