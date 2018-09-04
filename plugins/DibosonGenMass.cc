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
#include "TMath.h"


class DibosonGenMass : public edm::EDFilter {
 public:
  explicit DibosonGenMass(const edm::ParameterSet&);

 private:
 virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::InputTag src;
  const double min_gen;
  const double max_gen;

  double eventWeight;
  bool useMadgraphWeight;
  double madgraphWeight;

  TH1F* pre_Gen_Weight;
  TH1F* pre_res_mass;
  TH1F* pre_res_pt;
  //TH1F* pre_res_eta;
  TH1F* pre_res_rap;
  TH1F* pre_res_phi;
  TH1F* Gen_Weight;
  TH1F* res_mass;
  TH1F* res_pt;
  //TH1F* res_eta;
  TH1F* res_rap;
  TH1F* res_phi;
};


DibosonGenMass::DibosonGenMass(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),
    min_gen(cfg.getParameter<double>("min_gen")),
    max_gen(cfg.getParameter<double>("max_gen")),
    eventWeight(1.0),
    useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    madgraphWeight(1.0)
{
  consumes<reco::GenParticleCollection>(src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));

  edm::Service<TFileService> fs;
  pre_Gen_Weight = fs->make<TH1F>("pre_Gen_Weight", "", 4, -2, 2);
  pre_res_mass = fs->make<TH1F>("pre_res_mass", "", 20000, 0, 20000);
  pre_res_pt = fs->make<TH1F>("pre_res_pt", "", 20000, 0, 20000);
  //pre_res_eta = fs->make<TH1F>("pre_res_eta", "", 100, -5, 5);
  pre_res_rap = fs->make<TH1F>("pre_res_rap", "", 100, -5, 5);
  pre_res_phi = fs->make<TH1F>("pre_res_phi", "", 100, -TMath::Pi(), TMath::Pi());
  Gen_Weight = fs->make<TH1F>("Gen_Weight", "", 4, -2, 2);
  res_mass = fs->make<TH1F>("res_mass", "", 20000, 0, 20000);
  res_pt = fs->make<TH1F>("res_pt", "", 20000, 0, 20000);
  //res_eta = fs->make<TH1F>("res_eta", "", 100, -5, 5);
  res_rap = fs->make<TH1F>("res_rap", "", 100, -5, 5);
  res_phi = fs->make<TH1F>("res_phi", "", 100, -TMath::Pi(), TMath::Pi());
}

bool DibosonGenMass::filter(edm::Event& event, const edm::EventSetup&) {

  eventWeight = 1.0;
  madgraphWeight = 1.0;

  if (useMadgraphWeight) {
    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid() ){
      eventWeight = gen_ev_info->weight();
      madgraphWeight = ( eventWeight > 0 ) ? 1.0 : -1.0;
    }
  }
  else {
    eventWeight = 1.0;
    madgraphWeight = 1.0;
  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByLabel(src, genParticles);

  reco::GenParticleCollection::const_iterator genp = genParticles->begin();

  bool isFind1 = false;
  bool isFind2 = false;

  float pt1 = 0;
  float pt2 = 0;
  float eta1 =0;
  float eta2 =0;
  float phi1 =0;
  float phi2 =0;
  float mass1 = 0;
  float mass2 = 0;
  TLorentzVector mu1, mu2;
  TLorentzVector WW;
  float massWW = -999.0;

  for(; genp != genParticles->end(); genp++) {

    /*if( genp->isHardProcess() ||
       (genp->pdgId() == 13 || genp->pdgId() == 11  || genp->pdgId() == 15) ||
       (genp->pdgId() == -13 || genp->pdgId() ==-11  || genp->pdgId() == -15) ||
       (genp->pdgId() == 24 || genp->pdgId() == -24) ||
       (genp->pdgId() == 23 || genp->pdgId() == -23) ||
       (genp->pdgId() == 6 || genp->pdgId() == -6) ) {
      const reco::Candidate* mom = genp->mother();
      //std::cout<<"    genp Id ="<<genp->pdgId()<<"  genp status = "<<genp->status()<<"  mother ="<<mom->pdgId() << "  isHardProcess = " << genp->isHardProcess() << std::endl;
    }*/

    if( (genp->pdgId() == 13 || genp->pdgId() == 11  || genp->pdgId() == 15) && (genp->isHardProcess()) ){
      const reco::Candidate* m = genp->mother();
      //std::cout<<"lep1 Id ="<<genp->pdgId()<<"  status="<<genp->status()<<"  mother="<< m->pdgId() << std::endl;

      if(m->pdgId()==-24){
        isFind1 = true;
        eta1 = genp->eta();
        phi1 = genp->phi();
        pt1 = genp->pt();
        mass1 = genp->mass();
        mu1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);  //0.105
      }
    }

    if( (genp->pdgId() == -13 || genp->pdgId() ==-11  || genp->pdgId() == -15) && (genp->isHardProcess()) ){
      const reco::Candidate* m2 = genp->mother();
      //std::cout<<"lep2 Id ="<<genp->pdgId()<<"  status="<<genp->status()<<"  mother="<< m2->pdgId() << std::endl;

      if(m2->pdgId()==24){
        isFind2 = true;
        eta2 = genp->eta();
        phi2 = genp->phi();
        pt2 = genp->pt();
        mass2 = genp->mass();
        mu2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);  //0.105
      }
    }

  }

  WW = mu1+ mu2;

  massWW = WW.M();

  //std::cout << "pt1=" << pt1 << ", eta1=" << eta1 << ", phi1=" << phi1 << ", mass1=" << mass1 << std::endl;
  //std::cout << "pt2=" << pt2 << ", eta2=" << eta2 << ", phi2=" << phi2 << ", mass2=" << mass2 << std::endl;
  //std::cout << "WW mass = "<<massWW << "    " << (massWW >= min_gen && massWW<=max_gen) << std::endl;

  if(isFind1 && isFind2){
    //std::cout << "Fill gen histograms" << std::endl;
    pre_Gen_Weight->Fill(madgraphWeight);
    pre_res_mass->Fill(WW.M(), madgraphWeight);
    pre_res_pt->Fill(WW.Pt(), madgraphWeight);
    //pre_res_eta->Fill(WW.Eta(), madgraphWeight);
    pre_res_phi->Fill(WW.Phi(), madgraphWeight);
    pre_res_rap->Fill(WW.Rapidity(), madgraphWeight);

    if(massWW >= min_gen && massWW<=max_gen){
      Gen_Weight->Fill(madgraphWeight);
      res_mass->Fill(WW.M(), madgraphWeight);
      res_pt->Fill(WW.Pt(), madgraphWeight);
      //res_eta->Fill(WW.Eta(), madgraphWeight);
      res_phi->Fill(WW.Phi(), madgraphWeight);
      res_rap->Fill(WW.Rapidity(), madgraphWeight);
    }
  }
  //std::cout << std::endl;

  return
    isFind1 && isFind2 && massWW >= min_gen && massWW <= max_gen;
}


DEFINE_FWK_MODULE(DibosonGenMass);
