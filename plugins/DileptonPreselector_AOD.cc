// -*- C++ -*-
//
// Package:    DileptonPreselector_AOD
// Class:      DileptonPreselector_AOD
// 
/**\class DileptonPreselector_AOD DileptonPreselector_AOD.cc SUSYBSMAnalysis/Zprime2muAnalysis/DileptonPreselector_AOD.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/



// system include files
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
#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//
// class declaration
//

class DileptonPreselector_AOD : public edm::EDFilter {
public:
  explicit DileptonPreselector_AOD(const edm::ParameterSet&);
  ~DileptonPreselector_AOD();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::EDGetTokenT< std::vector< reco::Muon > >  muonToken_;
  int ptCut_;
  int multiplicityCut_;


  bool debug;
};

// constructors and destructor
DileptonPreselector_AOD::DileptonPreselector_AOD(const edm::ParameterSet& iConfig):
 muonToken_ (consumes< std::vector< reco::Muon > > (iConfig.getParameter<edm::InputTag>("muons")))
{
  ptCut_    = iConfig.getParameter < double > ("ptCut");
  multiplicityCut_    = iConfig.getParameter < double > ("nMuons");

  debug = false;
}

DileptonPreselector_AOD::~DileptonPreselector_AOD(){}


// member functions
// ------------ method called on each new Event  ------------
bool
DileptonPreselector_AOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< std::vector< reco::Muon > > muons;  
  iEvent.getByToken(muonToken_, muons);
  bool filter = false;

  int nMuons = 0;

  for (std::vector<reco::Muon>::const_iterator mu = muons->begin(); mu != muons->end(); mu++){

	if ((*mu).tunePMuonBestTrack()->pt() > ptCut_ && (*mu).isGlobalMuon()) nMuons++;


  }

  if (nMuons >= multiplicityCut_) filter = true;
  return filter;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DileptonPreselector_AOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DileptonPreselector_AOD::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DileptonPreselector_AOD);
