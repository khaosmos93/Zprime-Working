// -*- C++ -*-
//
// Package:    DileptonMassFilter
// Class:      DileptonMassFilter
// 
/**\class DileptonMassFilter DileptonMassFilter.cc SUSYBSMAnalysis/Zprime2muAnalysis/DileptonMassFilter.cc

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
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TLorentzVector.h"
//#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//
// class declaration
//

class DileptonMassFilter : public edm::EDFilter {
public:
  explicit DileptonMassFilter(const edm::ParameterSet&);
  ~DileptonMassFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::EDGetTokenT< std::vector< pat::Muon > >  muonToken_;
  double ptCut_;
  int multiplicityCut_;
  double massCut_;

  double mu_mass;

  bool debug;
};

// constructors and destructor
DileptonMassFilter::DileptonMassFilter(const edm::ParameterSet& iConfig):
 muonToken_ (consumes< std::vector< pat::Muon > > (iConfig.getParameter<edm::InputTag>("muons")))
{
  ptCut_    = iConfig.getParameter < double > ("ptCut");
  multiplicityCut_    = iConfig.getParameter < double > ("nMuons");
  massCut_    = iConfig.getParameter < double > ("massCut");

  mu_mass = 0.105658;
  debug = iConfig.getParameter < bool > ("debug");
  std::cout << "debug = " << debug << std::endl;
}

DileptonMassFilter::~DileptonMassFilter(){}


// member functions
// ------------ method called on each new Event  ------------
bool
DileptonMassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< std::vector< pat::Muon > > muons;
  iEvent.getByToken(muonToken_, muons);
  bool filter = false;

  if(muons->size() < 2)
    return false;

  int nMuons = 0;
  double MaxMass = -999.;

  for (std::vector<pat::Muon>::const_iterator mu0 = muons->begin(); mu0 != muons->end(); mu0++){
    if ((*mu0).tunePMuonBestTrack()->pt() > ptCut_ && (*mu0).isGlobalMuon()) nMuons++;

    TLorentzVector vec4_mu0;
    vec4_mu0.SetPtEtaPhiM(mu0->pt(), mu0->eta(), mu0->phi(), mu_mass);

    if(mu0+1 != muons->end() ) {
      for (std::vector<pat::Muon>::const_iterator mu1 = mu0+1; mu1 != muons->end(); mu1++){
        TLorentzVector vec4_mu1;
        vec4_mu1.SetPtEtaPhiM(mu1->pt(), mu1->eta(), mu1->phi(), mu_mass);
        TLorentzVector vec4_pair = vec4_mu0+vec4_mu1;
        if(vec4_pair.M() > MaxMass){
          MaxMass = vec4_pair.M();
        }
      }
    }

  }

  if (nMuons >= multiplicityCut_ &&
      MaxMass >= massCut_) {
    if(debug)  std::cout << "nMuons : " << nMuons << std::endl;
    if(debug)  std::cout << "MaxMass : " << MaxMass << std::endl;
    filter = true;
  }
  return filter;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DileptonMassFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DileptonMassFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DileptonMassFilter);
