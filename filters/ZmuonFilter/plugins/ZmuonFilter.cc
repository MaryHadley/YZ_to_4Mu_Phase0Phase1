// -*- C++ -*-
//
// Package:    filters/ZmuonFilter
// Class:      ZmuonFilter
// 
/**\class ZmuonFilter ZmuonFilter.cc filters/ZmuonFilter/plugins/ZmuonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mary Hill Hadley
//         Created:  Tue, 08 Sep 2020 18:57:22 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

//new includes  //need to specify these in the BuildFile.xml in the plugins directory, thanks to Kevin Pedro for explaining this
//By these I mean things like PhysicsTools, DataFormats, then what you use specifically you specify here in the .cc
//so for example in my BuildFile.xml I have: <use name="PhysicsTools/UtilAlgos"/> and then here I have stuff that comes out of PhysicsTools
//Ask John to elaborate on this because I have for instance DataFormats/VertexReco in my BuildFile.xml but then here I actually use DataFormats/PatCandidates and it doesn't complain...
//I think the answer at this point (after talking to John) is: this is done on a wing and a prayer, maybe ask K.P. later this week

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "FWCore/Common/interface/TriggerNames.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
//#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h" 


#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "FWCore/ServiceRegistry/interface/Service.h" //to make histos
#include "CommonTools/UtilAlgos/interface/TFileService.h"  //to make histos

//
// class declaration
//

class ZmuonFilter : public edm::stream::EDFilter<> {
   public:
      explicit ZmuonFilter(const edm::ParameterSet&);
      ~ZmuonFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      
      edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

      
      std::vector<std::string> triggerlist;
      double pTCut, etaCut, invMass4MuCut_low; // invMass4MuCut_high; //these come from the cfg
      bool verboseFilter = true; //instead of having it come from the cfg, doing it out a hacky way here
      double muon_mass = 0.1056583715; 
      
      std::unordered_map<std::string,TH1*> phase0_histContainer_; //for phase0_CutFlow histo
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZmuonFilter::ZmuonFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
//put needed tokens here
// define eta cuts, 4 object invariant mass cut, pT cut (if needed), triggercuts
// it looks like, based on Frank's example here: https://github.com/fojensen/ExcitedTau/blob/master/excitingAnalyzer_cfg.py#L137-L143, that this is how to add it to the cfg

muonsToken_        = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"));
triggerBits_       = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));
triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"));

pTCut = iConfig.getParameter<double>("pTCut");
etaCut = iConfig.getParameter<double>("etaCut");
invMass4MuCut_low = iConfig.getParameter<double>("invMass4MuCut_low");
//invMass4MuCut_high = iConfig.getParameter<double>("invMass4MuCut_high");
//verboseFilter      = iConfig.getParameter<bool>("verboseFilter"); //Not sure this will work when applyZmuonFilter is set to false, so trying it a hacky way first

edm::Service<TFileService> fs; //creating a TFileService instance

phase0_histContainer_ ["phase0_CutFlow"]   = fs->make<TH1F>("phase0_CutFlow",  ";phase0_CutFlow;Z+Upsi Candidate",4,0,4); //creating phase0_CutFlow histo in the phase0_histContainer_
for(std::unordered_map<std::string,TH1*>::iterator it=phase0_histContainer_.begin();   it!=phase0_histContainer_.end();   it++) it->second->Sumw2(); //call Sumw2 on all the hists in phase0_histContainer_  

}


ZmuonFilter::~ZmuonFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZmuonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{//I think this is the one I want and I can kill of the EVENTSETUP_EXAMPLE
   using namespace edm;
//since using namespace edm probably don't need to write edm before Handle below, leaving it for now
 
 //Trigger
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

//  edm::Handle<pat::PackedCandidateCollection> pfs; //do I need this
//  iEvent.getByToken(pfToken_, pfs); //do I need this
 
  
  //muons
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  
  //Flag initializations for each event 
  bool flagAtLeast4Mu = false; //set to true, commented out 4Mu code so I know it's not doing anything  
//  std::cout << "flagAtLeast4Mu is initialized to: " << flagAtLeast4Mu << std::endl;
  
  bool flagPassTrigger = false;//change to true so this cut is not doing anything right now
//  std::cout << "flagPassTrigger is initialized to: " << flagPassTrigger << std::endl;
  
  //Check if there are at least four mu in the event
  if ((int)muons->size() <= 3){
      if (verboseFilter){
         phase0_histContainer_["phase0_CutFlow"]->Fill(1);
      }
      return false; // If there are not at least four muons, the filter function will return false 
    }
  else { 
       flagAtLeast4Mu = true;
//        std::cout << "flatAtLeasat4Mu is true" << std::endl;
    
    }
   
   //Check triggers
   triggerlist.clear();
   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
//       std::cout << triggerBits->size() << std::endl;
       if (triggerBits->accept(i)) {
//           std::cout << "Got here" << std::endl;
           triggerlist.push_back(names.triggerName(i));
           std::string str (names.triggerName(i));
//           std::cout << "str is: " << str << std::endl;
           std::string str2 ("Mu");//confirmed that Mu is the right thing to search even for the DoubleMuon Dataset we will be using, see: https://twiki.cern.ch/twiki/bin/view/CMS/HLTPathsRunIIList
         
           std::size_t foundMu = str.find(str2);
//           std::cout << "Defined foundMu" <<  foundMu << std::endl;
           
           if (foundMu != std::string::npos) {
               flagPassTrigger = true;
 //              std::cout << "flagPassTrigger is true!" << std::endl;
               if (flagPassTrigger) { //could get rid of this if and just do a break after flagPassTrigger = true, can't decide right now which is more readable 
                   break;
               }
           }
                
        }
    

  }

 if (!flagPassTrigger) { //If, after going through all the trigger bits, we have not switched the flagPassTrigger to true and it is still false, reject the event 
//       std::cout << "DID NOT PASS TRIGGER" << std::endl;
       if (verboseFilter){
          phase0_histContainer_["phase0_CutFlow"]->Fill(2);
       }
       return false;
  }
  
  
  
  bool flagTotCharge_pT_Eta_InvMassOf4Mu = false; //set to false normally, set to true for debug  
//  std::cout << "flagTotCharge_pT_Eta_InvMassOf4Mu is initialized to: " << flagTotCharge_pT_Eta_InvMassOf4Mu << std::endl; 
  
   for (auto iM1 = muons->begin(); iM1 != muons->end(); ++iM1) {
      if  (flagTotCharge_pT_Eta_InvMassOf4Mu) break;  
      for (auto iM2 = iM1+1; iM2 != muons->end(); ++iM2) {
          if (flagTotCharge_pT_Eta_InvMassOf4Mu) break; 
          for (auto iM3 = iM2+1; iM3 != muons->end(); ++iM3) {
              if (flagTotCharge_pT_Eta_InvMassOf4Mu) break; 
              for (auto iM4 = iM3+1; iM4 != muons->end(); ++iM4) {
                  math::PtEtaPhiMLorentzVector lepton1(iM1->pt(), iM1->eta(), iM1->phi(), muon_mass);
                  math::PtEtaPhiMLorentzVector lepton2(iM2->pt(), iM2->eta(), iM2->phi(), muon_mass);
                  math::PtEtaPhiMLorentzVector lepton3(iM3->pt(), iM3->eta(), iM3->phi(), muon_mass);
                  math::PtEtaPhiMLorentzVector lepton4(iM4->pt(), iM4->eta(), iM4->phi(), muon_mass);
                  
                  if ((iM1->charge() + iM2->charge() + iM3->charge() + iM4->charge() == 0)
                    && iM1->pt() >= pTCut && iM2->pt() >= pTCut && iM3->pt() >= pTCut && iM4->pt() >= pTCut
                     && fabs(iM1->eta()) <= etaCut && fabs(iM2->eta()) <= etaCut && fabs(iM3->eta()) <= etaCut && fabs(iM4->eta()) <= etaCut
                    && (lepton1 + lepton2 + lepton3 + lepton4).mass() >= invMass4MuCut_low){
   //                  &&  (lepton1 + lepton2 + lepton3 + lepton4).mass() <= invMass4MuCut_high){
//  //                   std::cout << "invMass4Mu is  " <<(lepton1 + lepton2 + lepton3 + lepton4).mass() << std::endl;  //Recall that for events with more than one Y+Z candidate, the inv4MuMass might be over 120, but don't freak out, you are just guaranteeing there is at least 1 candidate in the event in the right range 
//                  //   if ((lepton1 + lepton2 + lepton3 + lepton4).mass() > 120) {
//                  //      std::cout << "BARK" << std::endl;
//                   //  }
                     flagTotCharge_pT_Eta_InvMassOf4Mu = true;
 //                    std::cout << "flagTotCharge_pT_Eta_InvMassOf4Mu is True" << std::endl;
                     break;
                    }
//                     else {
//                       std::cout << iM1->pt() << std::endl;
//                       std::cout << iM2->pt() << std::endl;
//                       std::cout << iM3->pt() << std::endl;
//                       std::cout << iM4->pt() << std::endl;
//                       std::cout << iM1->eta() << std::endl;
//                       std::cout << iM2->eta() << std::endl;
//                       std::cout << iM3->eta() << std::endl;
//                       std::cout << iM4->eta() << std::endl;
//                       std::cout << iM1->charge() << std::endl;
//                       std::cout << iM2->charge() << std::endl;
//                       std::cout << iM3->charge() << std::endl;
//                       std::cout << iM4->charge() << std::endl;
//                       
//                       std::cout << (lepton1 + lepton2 + lepton3 + lepton4).mass() << std::endl;
//                       
//                       
//                       //, iM2->pt(), iM3->pt(), iM4->pt(), iM1->eta(), iM2->eta(), iM3->eta(), iM4->eta() << std::endl;
//                   
//                   }
//                   
//                   
                  }
             }  
         }          
            
     }
 
  
   if (!flagTotCharge_pT_Eta_InvMassOf4Mu) { //If, after going through all the sets of 4 mu, we have not flipped our flagTotCharge_pT_Eta_InvMassOf4Mu to true, reject the event
//  //     std::cout << iM1->pt(), iM2->pt(), iM3->pt(), iM4->pt() << std::endl;
 //       std::cout << "FAILED flagTotCharge_pT_Eta_InvMassOf4Mu" << std::endl; 
       if (verboseFilter) {
          phase0_histContainer_["phase0_CutFlow"]->Fill(3);
       
       }
       return false;
   }
  
  if (flagAtLeast4Mu && flagPassTrigger && flagTotCharge_pT_Eta_InvMassOf4Mu) {
       return true; //If we have found at least four muons in the event, the event has passed our trigger of interest, and we can find at least one quartet of muons that satisfies our pT, eta, total charge, and total invariant mass of the four requirements, keep the event to pass to the analyzer!
  }
          
      
//Some examples of how to do things that might be helpful for reference.
    
// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif
// 
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif



   return false; //this bool defaults to false //I do not think it should ever get down here, I think it should return true or false before then, but I guess it is a protection
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ZmuonFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ZmuonFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ZmuonFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ZmuonFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ZmuonFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ZmuonFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZmuonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  
  
  edm::ParameterSetDescription desc;
//  desc.setUnknown(); // comment back in if you want to go to the default of not knowing what parameters are allowed so do no validation
//descriptions.addDefault(desc); // // comment back in if you want to go to the default of not knowing what parameters are allowed so do no validation
//If you wanted to do no validation, comment out all the stuff that is specific to this EDFilter, so muonCollection all the way through the  descriptions.add("ZmuonFilter", desc); line


//Needed parameters for the ZmuonFilter 
 //all values provided are default paramters that can be overwritten in the cfg (probably will never overwrite what muonCollection, bits, and objects are, will change and tune the <blah>Cut variables
 
 desc.add<edm::InputTag>("muonCollection",edm::InputTag("slimmedMuons"));
 desc.add<edm::InputTag>("bits",edm::InputTag("TriggerResults","", "HLT"));
 desc.add<edm::InputTag>("objects",edm::InputTag("selectedPatTrigger"));
 desc.add<double>("pTCut", 0); //these values provided for ptCut, etaCut, invMass4MuCut_low are defaults and will be overwritten by what you put in the cfg
 desc.add<double>("etaCut", 4); //I have put in really loose cuts that basically do nothing here
 desc.add<double>("invMass4MuCut_low", 0); //Thank you to Gabriele Benelli and Jan-Frederik Schulte for their tips on getting this fillDescription part working!
// desc.add<double>("invMass4MuCut_high", 10000);
 descriptions.add("ZmuonFilter", desc);


}
//define this as a plug-in
DEFINE_FWK_MODULE(ZmuonFilter);
