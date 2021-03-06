#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/Event.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerUtilities.h"

Zprime2muTriggerPathsAndFilters::Zprime2muTriggerPathsAndFilters(const edm::Event& event) {
  // Here we explicitly specify the HLT path/final filter name as a
  // function of run number, rather than relying on HLTConfigProvider
  // to give it to us, or using globbing on different path version
  // numbers, etc.

  // During running, online HLT menus change often so maintaining this
  // list is burdensome. Dealing with subtle bugs/inconsistencies from
  // assumptions made is worse.

  // If none of the data/MC/run number conditions below are
  // satisfied, path and filter will remain empty strings, and the
  // valid bit will be set false.

  valid = true; // assume until proven otherwise
  const unsigned run = event.id().run();
  // JMTBAD if there are different HLT menus used for different
  // samples (e.g. 51X/52X), the next line may not be sufficient, and
  // we will have to think of a better way to handle this.


  if (!event.isRealData())                { path = "HLT_Mu50_v10", filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q", prescaled_path = "HLT_Mu27_v10", prescaled_filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"; }  //92X MC

  else if (run >= 272007 and run < 289999)                 { path = "HLT_Mu50_v4", filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q", prescaled_path = "HLT_Mu27_v2", prescaled_filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"; }

  else if (run >= 290000 and run < 297099)                 { path = "HLT_Mu50_v6", filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q", prescaled_path = "HLT_Mu27_v6", prescaled_filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"; }
  else if (run >= 297099 and run < 297557)                 { path = "HLT_Mu50_v7", filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q", prescaled_path = "HLT_Mu27_v7", prescaled_filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"; }
  else if (run >= 297557 and run < 299368)                 { path = "HLT_Mu50_v8", filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q", prescaled_path = "HLT_Mu27_v8", prescaled_filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"; }
  else if (run >= 299368 and run < 299650)                 { path = "HLT_Mu50_v9", filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q", prescaled_path = "HLT_Mu27_v9", prescaled_filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"; }
  else if (run >= 299650 and run < 302020)                 { path = "HLT_Mu50_v10", filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q", prescaled_path = "HLT_Mu27_v10", prescaled_filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"; }
  else if (run >= 302020 and run < 999999)                 { path = "HLT_Mu50_v11", filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q", prescaled_path = "HLT_Mu27_v11", prescaled_filter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"; }

    else
        valid = false;
}

trigger::TriggerObjectCollection get_L3_muons(const edm::Event& event, const std::string& filter, const edm::InputTag& trigger_summary_src, const std::string& collection) {
  trigger::TriggerObjectCollection L3_muons;

  edm::Handle<trigger::TriggerEvent> trigEvent;

  event.getByLabel(trigger_summary_src, trigEvent);
  if (!trigEvent.isValid())
    throw cms::Exception("get_L3_muons") << "couldn't get hltTriggerSummaryAOD " << trigger_summary_src.encode() << " from event\n";

  // The TriggerEvent object keeps a list of all trigger-firing
  // objects, associated to the original collection module name
  // (usually "hltL3MuonCandidates"). Determine which entries are
  // muons, with "keys" being indices into the TriggerObjectCollection
  // below. If the collection is not found, the final loop will not
  // run thanks to the sentinel values.
  std::pair<int, int> collection_keys(-1,-2);
  int key_prev = 0;
  for (trigger::size_type iC = 0; iC < trigEvent->sizeCollections(); ++iC) {
    int key = trigEvent->collectionKey(iC);
    if (trigEvent->collectionTag(iC).label() == collection) {
      collection_keys = std::make_pair(key_prev, key - 1);
      break;
    }
    key_prev = key;
  }

  // Same idea for the filter name, but a slightly different way to
  // get the keys. If the filter is not found, the keys_passing_filter
  // vector will be empty and no muons will be kept in the loop below.
  trigger::Keys keys_passing_filter;
  for (trigger::size_type iF = 0; iF < trigEvent->sizeFilters(); ++iF) {
    if (trigEvent->filterTag(iF).label() == filter) {
      keys_passing_filter = trigEvent->filterKeys(iF);
      break;
    }
  }

  // Finally, run over the pre-selected range of keys to grab the
  // entries that also passed the filter.
  const trigger::TriggerObjectCollection& TOC(trigEvent->getObjects());
  for (trigger::size_type iO = collection_keys.first; iO <= collection_keys.second; ++iO)
    if (std::find(keys_passing_filter.begin(), keys_passing_filter.end(), iO) != keys_passing_filter.end())
      L3_muons.push_back(TOC.at(iO));

    //std::cout<<"filter: "<<filter<<"collection size: "<<L3_muons.size()<<std::endl;
  return L3_muons;
}

