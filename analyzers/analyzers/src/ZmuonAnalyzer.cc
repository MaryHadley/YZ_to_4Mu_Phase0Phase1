// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// new includes
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>

#include <typeinfo>

//test comment
//test comment 2 //test comment 3

// for vertexing                                                                                                                                                                                        
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"                                                                                             
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
     
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"                                                                                                       
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"
             
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"                                                                                                 
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" //for pileup 

#include <math.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <algorithm>  // std::sort, std::swap
#include <iostream>  // std::cout, std::endl
#include <string>

#include "FWCore/Framework/interface/LuminosityBlock.h" //Needed for EventCountProducer 

#include "DataFormats/Common/interface/MergeableCounter.h" //EventCountProducer produces an instance edm::MergeableCounter, so we need this

#include<TH2.h> //to get histContainer2_ to work, thanks to Markus Seidel for input on this

class ZmuonAnalyzer : public edm::EDAnalyzer {
public:
	explicit ZmuonAnalyzer(const edm::ParameterSet&);
private:
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
   edm::EDGetTokenT<std::vector<reco::Vertex>> token_vertices;
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
//   edm::EDGetTokenT<pat::TriggerEvent> triggerEventToken_;
//   edm::EDGetTokenT<std::vector<std::string>> muonMatch_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
   
   edm::EDGetTokenT<edm::MergeableCounter> nTotalToken_; // need for EventCountProducer
   edm::EDGetTokenT<edm::MergeableCounter> nFilteredToken_; //need for EventCountProducer
   
   std::unordered_map<std::string,TH1*> histContainer_; //for counters, PU, Cutflow, trigger DR, etc
   
   std::unordered_map<std::string,TH2*> histContainer2_;
   
   edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_; //for pileup
   
   edm::EDGetTokenT<double> rhoToken_; //rho, median energy density, see: https://arxiv.org/pdf/0707.1378.pdf
   
   virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override; //If applyZmuonFilter=False, comment out!
   
   bool isMC;
   
   //bool use2018Triggers, use2017Triggers, use2016Triggers;
   
   int triggerYear;
   
   int PU, PU_True;
   
   double iTrkWeightOfRecoVtxTrk;

   TTree * tree;
   TTree * treemc;
   
   std::vector<bool> flagZplusY;
   std::vector<bool> flagZOnly;
   
   std::vector<double> lepton1_pt, lepton1_eta, lepton1_phi;
   std::vector<double> lepton1_charge;
   std::vector<double> lepton1_d0, lepton1_dz, lepton1_dxy; //recall d0 is the 3D IP
   std::vector<double> lepton1_iso03particle, lepton1_iso03hadron, lepton1_iso04hadron, lepton1_iso04particle; //these are the R = 0.4 or 0.3 isolation variables 
   std::vector<double> lepton1_iso03neutralHadron, lepton1_iso03photon, lepton1_iso03PU; 
   std::vector<double> lepton1_impactParameterSignificance;
   std::vector<bool> lepton1_isLooseMuon, lepton1_isSoftMuon, lepton1_isTightMuon, lepton1_isMediumMuon;
   std::vector<bool> lepton1_isPFMuon, lepton1_isGlobalMuon, lepton1_isTrackerMuon;
   std::vector<double> lepton1_deltaEta, lepton1_deltaPhi, lepton1_sigmaiEtaEta, lepton1_HoverE; //these are quantities for electrons, I'm leaving the branch here for ease of recycling when we modify the code to deal with the electron channel, but in this case, they will NOT be filled since we just have mu  //see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/EgammaPublicData#H_E
   std::vector<double> lepton1_OoEmOoP, lepton1_vtxFitProb, lepton1_missingHits; //OoEmOop is 1/E - 1/P, this is an electron quantity, leaving for ease of recycling  when we make a version of the code to do the el channel. The fitProb for a single mu seems nonsensical to us right now (28 Sept. 2020) so is not filled. Waiting for confirmation about right way to implement missing hits from LPC mattermost. //anyway these quantities are probably not that important 
   std::vector<double> lepton1_validHits;
   
   std::vector<double> lepton2_pt, lepton2_eta, lepton2_phi;
   std::vector<double> lepton2_charge;
   std::vector<double> lepton2_d0, lepton2_dz, lepton2_dxy;
   std::vector<double> lepton2_iso03particle, lepton2_iso03hadron, lepton2_iso04hadron, lepton2_iso04particle;
   std::vector<double> lepton2_iso03neutralHadron, lepton2_iso03photon, lepton2_iso03PU;
   std::vector<double> lepton2_impactParameterSignificance;
   std::vector<bool> lepton2_isLooseMuon, lepton2_isSoftMuon, lepton2_isTightMuon, lepton2_isMediumMuon;
   std::vector<bool> lepton2_isPFMuon, lepton2_isGlobalMuon, lepton2_isTrackerMuon;
   std::vector<double> lepton2_deltaEta, lepton2_deltaPhi, lepton2_sigmaiEtaEta, lepton2_HoverE; //these are quantities for electrons, I'm leaving the branch here for ease of recycling when we modify the code to deal with the electron channel, but in this case, they will NOT be filled since we just have mu  //see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/EgammaPublicData#H_E
   std::vector<double> lepton2_OoEmOoP, lepton2_vtxFitProb, lepton2_missingHits; //OoEmOop is 1/E - 1/P, this is an electron quantity, leaving for ease of recycling  when we make a version of the code to do the el channel. The fitProb for a single mu seems nonsensical to us right now (28 Sept. 2020) so is not filled. Waiting for confirmation about right way to implement missing hits from LPC mattermost. //anyway these quantities are probably not that important 
   std::vector<double> lepton2_validHits;
   
   std::vector<double> lepton3_pt, lepton3_eta, lepton3_phi;
   std::vector<double> lepton3_charge;
   std::vector<double> lepton3_d0, lepton3_dz, lepton3_dxy;
   std::vector<double> lepton3_iso03particle, lepton3_iso03hadron, lepton3_iso04hadron, lepton3_iso04particle;
   std::vector<double> lepton3_iso03neutralHadron, lepton3_iso03photon, lepton3_iso03PU;
   std::vector<double> lepton3_impactParameterSignificance;
   std::vector<bool> lepton3_isLooseMuon, lepton3_isSoftMuon, lepton3_isTightMuon, lepton3_isMediumMuon;
   std::vector<bool> lepton3_isPFMuon, lepton3_isGlobalMuon, lepton3_isTrackerMuon;
   std::vector<double> lepton3_deltaEta, lepton3_deltaPhi, lepton3_sigmaiEtaEta, lepton3_HoverE;  //these are quantities for electrons, I'm leaving the branch here for ease of recycling when we modify the code to deal with the electron channel, but in this case, they will NOT be filled since we just have mu  //see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/EgammaPublicData#H_E
   std::vector<double> lepton3_OoEmOoP, lepton3_vtxFitProb, lepton3_missingHits;//OoEmOop is 1/E - 1/P, this is an electron quantity, leaving for ease of recycling  when we make a version of the code to do the el channel. The fitProb for a single mu seems nonsensical to us right now (28 Sept. 2020) so is not filled. Waiting for confirmation about right way to implement missing hits from LPC mattermost. //anyway these quantities are probably not that important 
   std::vector<double> lepton3_validHits;

   std::vector<double> lepton4_pt, lepton4_eta, lepton4_phi;
   std::vector<double> lepton4_charge;
   std::vector<double> lepton4_d0, lepton4_dz, lepton4_dxy;
   std::vector<double> lepton4_iso03particle, lepton4_iso03hadron, lepton4_iso04hadron, lepton4_iso04particle;
   std::vector<double> lepton4_iso03neutralHadron, lepton4_iso03photon, lepton4_iso03PU;
   std::vector<double> lepton4_impactParameterSignificance;
   std::vector<bool> lepton4_isLooseMuon, lepton4_isSoftMuon, lepton4_isTightMuon, lepton4_isMediumMuon;
   std::vector<bool> lepton4_isPFMuon, lepton4_isGlobalMuon, lepton4_isTrackerMuon;
   std::vector<double> lepton4_deltaEta, lepton4_deltaPhi, lepton4_sigmaiEtaEta, lepton4_HoverE; //these are quantities for electrons, I'm leaving the branch here for ease of recycling when we modify the code to deal with the electron channel, but in this case, they will NOT be filled since we just have mu  //see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/EgammaPublicData#H_E
   std::vector<double> lepton4_OoEmOoP, lepton4_vtxFitProb, lepton4_missingHits; //OoEmOop is 1/E - 1/P, this is an electron quantity, leaving for ease of recycling  when we make a version of the code to do the el channel. The fitProb for a single mu seems nonsensical to us right now (28 Sept. 2020) so is not filled. Waiting for confirmation about right way to implement missing hits from LPC mattermost. //anyway these quantities are probably not that important 
   std::vector<double> lepton4_validHits;
   
   //fromPV and PVAssociationQuality info
   std::vector<double> lepton1_fromPV, lepton2_fromPV, lepton3_fromPV, lepton4_fromPV;
   std::vector<double> lepton1_pvassq, lepton2_pvassq, lepton3_pvassq, lepton4_pvassq;
   
  // std::vector<double> dimuon1mass, dimuon2mass;
   std::vector<std::vector<double>> dimuon1mass, dimuon2mass; //Bjorn
  // std::vector<double> dimuon1pt,   dimuon2pt ;
   std::vector<std::vector<double>> dimuon1pt,   dimuon2pt; //Bjorn
  // std::vector<double> dimuon1eta,  dimuon2eta;
   std::vector<std::vector<double>> dimuon1eta,  dimuon2eta; //Bjorn
   //std::vector<double> dimuon1phi,  dimuon2phi;
   std::vector<std::vector<double>>  dimuon1phi,  dimuon2phi; //Bjorn
  // std::vector<double> dimuon1vtx,  dimuon2vtx;
   std::vector<std::vector<double>> dimuon1vtx, dimuon2vtx; //Bjorn
   std::vector<double> dimuon1lxy,  dimuon2lxy;
   std::vector<double> dimuon1lxysig,     dimuon2lxysig  ;
   std::vector<double> dimuon1lxyctauPV,  dimuon2lxyctauPV;
   std::vector<double> sixmuonvtx;
   //std::vector<double> dimuon1vtx_xpos, dimuon2vtx_xpos;
   std::vector<std::vector<double>> dimuon1vtx_xpos, dimuon2vtx_xpos; //Bjorn
   //std::vector<double> dimuon1vtx_ypos, dimuon2vtx_ypos;
   std::vector<std::vector<double>> dimuon1vtx_ypos, dimuon2vtx_ypos; //Bjorn
   //std::vector<double> dimuon1vtx_zpos, dimuon2vtx_zpos;
   std::vector<std::vector<double>> dimuon1vtx_zpos, dimuon2vtx_zpos; //Bjorn
   //std::vector<double> dimuon1vtx_xposError, dimuon2vtx_xposError;
    std::vector<std::vector<double>> dimuon1vtx_xposError, dimuon2vtx_xposError; //Bjorn
  // std::vector<double> dimuon1vtx_yposError, dimuon2vtx_yposError;
   std::vector<std::vector<double>> dimuon1vtx_yposError, dimuon2vtx_yposError; //Bjorn
  // std::vector<double> dimuon1vtx_zposError, dimuon2vtx_zposError;
  std::vector<std::vector<double>> dimuon1vtx_zposError, dimuon2vtx_zposError; //Bjorn
   std::vector<double> rho;
   std::vector<unsigned int>  event_number;
   std::vector<unsigned int>  run_number;
   std::vector<unsigned int> lumi_section;
   std::vector<unsigned int> orbit_number;
   std::vector<unsigned int> bunch_crossing;
//   std::vector<unsigned int> process_ID;
   
   std::vector<unsigned int> save_event_count;

   std::vector<double> PVx, PVy, PVz;
   
   std::vector<int> quadComesFromEvtWithThisManyMu;

   std::vector<double> truth_Zmuon_pt, truth_Zmuon_eta, truth_Zmuon_phi;
   std::vector<double> truth_Z_pt, truth_Z_eta, truth_Z_phi, truth_Z_mass, truth_Z_pdgid;
   std::vector<double> mc_event_number;
   std::vector<double> mc_run_number;
   std::vector<unsigned int> mc_lumi_section;

   std::vector<double> truth_Upsimuon_pt, truth_Upsimuon_eta, truth_Upsimuon_phi;
   std::vector<double> truth_Upsi_pt, truth_Upsi_eta, truth_Upsi_phi, truth_Upsi_mass;
   
   std::vector<double> truth_Upsi2muon_pt, truth_Upsi2muon_eta, truth_Upsi2muon_phi;
   std::vector<double> truth_Upsi2_pt, truth_Upsi2_eta, truth_Upsi2_phi, truth_Upsi2_mass;
   
   std::vector<double> truth_Upsi3muon_pt, truth_Upsi3muon_eta, truth_Upsi3muon_phi;
   std::vector<double> truth_Upsi3_pt, truth_Upsi3_eta, truth_Upsi3_phi, truth_Upsi3_mass;
   
   //Chib0_1P state
   std::vector<double> truth_Chib0_1P_UPSI_muon_pt, truth_Chib0_1P_UPSI_muon_eta, truth_Chib0_1P_UPSI_muon_phi;
   std::vector<double> truth_Chib0_1P_pt, truth_Chib0_1P_eta, truth_Chib0_1P_phi;
   std::vector<double> truth_Chib0_1P_pdgid;
   std::vector<double>truth_Chib0_1P_mass;
  
   
   std::vector<double>truth_Chib0_1P_UPSI_pt, truth_Chib0_1P_UPSI_eta, truth_Chib0_1P_UPSI_phi;
   std::vector<double>truth_Chib0_1P_UPSI_pdgid;
   std::vector<double>truth_Chib0_1P_UPSI_mass;
   
   std::vector<double>truth_Chib0_1P_photon_pt, truth_Chib0_1P_photon_eta, truth_Chib0_1P_photon_phi;
   std::vector<double>truth_Chib0_1P_photon_pdgid;
   
   //Chib1_1P state
   std::vector<double>truth_Chib1_1P_UPSI_muon_pt, truth_Chib1_1P_UPSI_muon_eta, truth_Chib1_1P_UPSI_muon_phi;
   
   std::vector<double>truth_Chib1_1P_pt, truth_Chib1_1P_eta, truth_Chib1_1P_phi;
   
   std::vector<double>truth_Chib1_1P_pdgid;
   std::vector<double>truth_Chib1_1P_mass;
   
   std::vector<double>truth_Chib1_1P_UPSI_pt, truth_Chib1_1P_UPSI_eta, truth_Chib1_1P_UPSI_phi;
   std::vector<double>truth_Chib1_1P_UPSI_pdgid;
   std::vector<double>truth_Chib1_1P_UPSI_mass;
   
  std::vector<double>truth_Chib1_1P_photon_pt, truth_Chib1_1P_photon_eta, truth_Chib1_1P_photon_phi;
  std::vector<double>truth_Chib1_1P_photon_pdgid;
  
  //Chib2_1P state
  std::vector<double>truth_Chib2_1P_UPSI_muon_pt, truth_Chib2_1P_UPSI_muon_eta, truth_Chib2_1P_UPSI_muon_phi;
  
  std::vector<double>truth_Chib2_1P_pt,truth_Chib2_1P_eta, truth_Chib2_1P_phi;
  std::vector<double>truth_Chib2_1P_pdgid;
  std::vector<double>truth_Chib2_1P_mass;
  
  std::vector<double>truth_Chib2_1P_UPSI_pt, truth_Chib2_1P_UPSI_eta, truth_Chib2_1P_UPSI_phi;
  std::vector<double>truth_Chib2_1P_UPSI_pdgid;
  std::vector<double>truth_Chib2_1P_UPSI_mass;
  
  std::vector<double>truth_Chib2_1P_photon_pt,truth_Chib2_1P_photon_eta, truth_Chib2_1P_photon_phi;
  std::vector<double>truth_Chib2_1P_photon_pdgid;
   
   

   std::vector<double> truth_Upsi_pdgid;
   std::vector<double> truth_Upsi2_pdgid;
   std::vector<double> truth_Upsi3_pdgid;
   
   std::vector<double> loop_enter_check;

   std::vector<std::string> triggerlist;

   std::vector<double> numberOfVertices;
   std::vector<double> zOfVerticesError;
   std::vector<double> zOfVertices;

   std::vector<bool> pair_12_34_56;
   std::vector<bool> pair_13_24_56;
   std::vector<bool> pair_14_23_56;
   
   std::vector<bool> pair_12_34_ZOnly;
   std::vector<bool> pair_13_24_ZOnly;
   std::vector<bool> pair_14_23_ZOnly; 
   
   
   std::vector<double> inv4MuMass, big4MuEta, big4MuPhi, big4MuPt;
   
   std::vector<double> big4MuVtx;
   
   std::vector<int> eventHasZUpsiNTo4Mu_Count;
   
   //trigger matching variable
   std::vector<int> quadHasHowManyTrigMatches;
   
   //sanity check
   std::vector<int> sum4LepCharges;
   
//   std::vector<double> TrkWeightsRecoVtxTrks;

};  

ZmuonAnalyzer::ZmuonAnalyzer(const edm::ParameterSet& iConfig)
{
   muonsToken_        = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"));
   token_vertices     = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   triggerBits_       = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));
//   triggerEventToken_ = consumes<pat::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerEvent"));
//   muonMatch_         = consumes<std::vector<std::string>>(iConfig.getParameter< std::string >( "muonMatch" ) );
   triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"));
   genParticlesToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
   pfToken_ = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"));
   
   
   puToken_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));  //hard coding this because I doubt I will change it, for pileup HARDCODED 
   
   rhoToken_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
  
   nTotalToken_ = consumes<edm::MergeableCounter, // for EventCountProducer
   edm::InLumi>(edm::InputTag("nEventsTotal")); //for EventCountProducer 
   
   nFilteredToken_ = consumes<edm::MergeableCounter, // for EventCountProducer
   edm::InLumi>(edm::InputTag("nEventsFiltered")); //for EventCountProducer 

   
   isMC = iConfig.getParameter<bool>("isMC"); // ***** MC *** gets read in from ZmuonAnalyzer_cfg!
   triggerYear = iConfig.getParameter<int>("triggerYear");
   
   
   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree", "tree");
   
   //counter histograms 
   histContainer_["nEventsTotalCounter"]    = fs->make<TH1F>("nEventsTotalCounter", ";nEventsTotalCounter; nEvents",2,0,2);
   histContainer_["nEventsFilteredCounter"]  = fs->make<TH1F>("nEventsFilteredCounter", ";nEventsFilteredCounter; nEvents",2,0,2);
   histContainer_["PU"]                      =  fs->make<TH1F>("PU",      ";Pileup observed;Events",100,0,100);
   histContainer_["PU_True"]                 =  fs->make<TH1F>("PU_True",  ";Pileup true;Events",100,0,100);
   histContainer_["CutFlow"]                 = fs->make<TH1F>("CutFlow",  ";CutFlow;Z+Upsi Candidate",20,0,20);
   histContainer_["trigDR"]                  = fs->make<TH1F>("trigDR", ";trigMu-RecoMu_DR; Num of trigMu-RecoMu Pairs", 200, 0, 0.2); //May adjust binning!! //this is a histo of the trig Mu-Reco Mu pair in a given quad that has the smallest deltaR separation
   histContainer_["numOfGoodMuInEvt"]        = fs->make<TH1F>("numOfGoodMuInEvt", ";Num of Good Mu in Event; Events", 20, -0.5, 19.5);
   histContainer_["numOfGoodMuPlusInEvt"]    = fs->make<TH1F>("numOfGoodMuPlusInEvt", ";Num of Good Mu+ in Event; Events", 20, -0.5, 19.5);
   histContainer_["numOfGoodMuMinusInEvt"]   = fs->make<TH1F>("numOfGoodMuMinusInEvt", ";Num of Good Mu- in Event; Events", 20, -0.5, 19.5);
  
   histContainer2_["numOfGoodMuPlusVsGoodMuMinus"] = fs->make<TH2F>("numOfGoodMuPlusVsGoodMuMinus", "Num of Good Mu+ vs. Good Mu- in Event; Num of Good Mu+; Num of Good Mu-", 20, -0.5, 19.5, 20, -0.5, 19.5);
   
   for(std::unordered_map<std::string,TH1*>::iterator it=histContainer_.begin();   it!=histContainer_.end();   it++) it->second->Sumw2(); //call Sumw2 on all the hists in histContainer_   
   for(std::unordered_map<std::string,TH2*>::iterator it=histContainer2_.begin();   it!=histContainer2_.end();   it++) it->second->Sumw2(); //call Sumw2 on all the hists in histContainer2_
   
   tree->Branch("event_number", &event_number);
   tree->Branch("run_number", &run_number);
   tree->Branch("lumi_section", &lumi_section);
   tree->Branch("bunch_crossing", &bunch_crossing);
   tree->Branch("orbit_number", &orbit_number);
//   tree->Branch("process_ID", &process_ID);
   
   tree->Branch("rho", &rho);
   
   tree->Branch("flagZplusY", &flagZplusY);
   tree->Branch("flagZOnly", &flagZOnly);
   

   tree->Branch("lepton1_pt",  &lepton1_pt);
   tree->Branch("lepton1_eta", &lepton1_eta);
   tree->Branch("lepton1_phi", &lepton1_phi);
   tree->Branch("lepton1_charge", &lepton1_charge);
   tree->Branch("lepton1_d0",  &lepton1_d0);
   tree->Branch("lepton1_dxy", &lepton1_dxy);
   tree->Branch("lepton1_dz",  &lepton1_dz);
   tree->Branch("lepton1_iso03particle", &lepton1_iso03particle);
   tree->Branch("lepton1_iso04particle", &lepton1_iso04particle);
   tree->Branch("lepton1_iso03hadron", &lepton1_iso03hadron); //charged hadron iso
   tree->Branch("lepton1_iso04hadron", &lepton1_iso04hadron);
   tree->Branch("lepton1_iso03neutralHadron", &lepton1_iso03neutralHadron); //neutral hadron iso
   tree->Branch("lepton1_iso03photon", &lepton1_iso03photon); //photon iso 
   tree->Branch("lepton1_iso03PU", &lepton1_iso03PU);//PU iso 
   tree->Branch("lepton1_impactParameterSignificance", &lepton1_impactParameterSignificance);
   tree->Branch("lepton1_isLooseMuon", &lepton1_isLooseMuon);    // these are only for muons 
   tree->Branch("lepton1_isMediumMuon", &lepton1_isMediumMuon);
   tree->Branch("lepton1_isSoftMuon",  &lepton1_isSoftMuon);     // these are only for muons  
   tree->Branch("lepton1_isTightMuon", &lepton1_isTightMuon);    // these are only for muons  
   tree->Branch("lepton1_isPFMuon", &lepton1_isPFMuon);          // these are only for muons  
   tree->Branch("lepton1_isGlobalMuon", &lepton1_isGlobalMuon);
   tree->Branch("lepton1_isTrackerMuon", &lepton1_isTrackerMuon);
   tree->Branch("lepton1_fromPV", &lepton1_fromPV);
   tree->Branch("lepton1_pvassq", &lepton1_pvassq);

   tree->Branch("lepton2_pt",  &lepton2_pt);
   tree->Branch("lepton2_eta", &lepton2_eta);
   tree->Branch("lepton2_phi", &lepton2_phi);
   tree->Branch("lepton2_charge", &lepton2_charge);
   tree->Branch("lepton2_d0",  &lepton2_d0);
   tree->Branch("lepton2_dxy", &lepton2_dxy);
   tree->Branch("lepton2_dz",  &lepton2_dz);
   tree->Branch("lepton2_iso03particle", &lepton2_iso03particle); 
   tree->Branch("lepton2_iso04particle", &lepton2_iso04particle);
   tree->Branch("lepton2_iso03hadron",   &lepton2_iso03hadron); //charged hadron iso 
   tree->Branch("lepton2_iso04hadron",   &lepton2_iso04hadron);
   tree->Branch("lepton2_iso03neutralHadron", &lepton2_iso03neutralHadron); //neutral hadron iso 
   tree->Branch("lepton2_iso03photon", &lepton2_iso03photon); //photon iso 
   tree->Branch("lepton2_iso03PU", &lepton2_iso03PU); //PU iso 
   tree->Branch("lepton2_impactParameterSignificance", &lepton2_impactParameterSignificance);
   tree->Branch("lepton2_isLooseMuon", &lepton2_isLooseMuon);    // these are only for muons 
   tree->Branch("lepton2_isSoftMuon",  &lepton2_isSoftMuon);     // these are only for muons
   tree->Branch("lepton2_isMediumMuon", &lepton2_isMediumMuon);  
   tree->Branch("lepton2_isTightMuon", &lepton2_isTightMuon);    // these are only for muons  
   tree->Branch("lepton2_isPFMuon", &lepton2_isPFMuon);          // these are only for muons  
   tree->Branch("lepton2_isGlobalMuon", &lepton2_isGlobalMuon);
   tree->Branch("lepton2_isTrackerMuon", &lepton2_isTrackerMuon);
   tree->Branch("lepton2_fromPV", &lepton2_fromPV);
   tree->Branch("lepton2_pvassq", &lepton2_pvassq);

   tree->Branch("lepton3_pt",  &lepton3_pt);
   tree->Branch("lepton3_eta", &lepton3_eta);
   tree->Branch("lepton3_phi", &lepton3_phi);
   tree->Branch("lepton3_charge", &lepton3_charge);
   tree->Branch("lepton3_d0",  &lepton3_d0);
   tree->Branch("lepton3_dxy", &lepton3_dxy);
   tree->Branch("lepton3_dz",  &lepton3_dz);
   tree->Branch("lepton3_iso03particle", &lepton3_iso03particle);
   tree->Branch("lepton3_iso04particle", &lepton3_iso04particle);
   tree->Branch("lepton3_iso03hadron", &lepton3_iso03hadron); //charged hadron iso 
   tree->Branch("lepton3_iso04hadron", &lepton3_iso04hadron);
   tree->Branch("lepton3_iso03neutralHadron", &lepton3_iso03neutralHadron); //neutral hadron iso 
   tree->Branch("lepton3_iso03photon", &lepton3_iso03photon); //photon iso 
   tree->Branch("lepton3_iso03PU", &lepton3_iso03PU); //PU iso 
   tree->Branch("lepton3_impactParameterSignificance", &lepton3_impactParameterSignificance);
   tree->Branch("lepton3_isLooseMuon", &lepton3_isLooseMuon);    // these are only for muons 
   tree->Branch("lepton3_isMediumMuon", &lepton3_isMediumMuon);
   tree->Branch("lepton3_isSoftMuon",  &lepton3_isSoftMuon);     // these are only for muons  
   tree->Branch("lepton3_isTightMuon", &lepton3_isTightMuon);    // these are only for muons  
   tree->Branch("lepton3_isPFMuon", &lepton3_isPFMuon);          // these are only for muons  
   tree->Branch("lepton3_isGlobalMuon", &lepton3_isGlobalMuon);
   tree->Branch("lepton3_isTrackerMuon", &lepton3_isTrackerMuon);
   tree->Branch("lepton3_fromPV", &lepton3_fromPV);
   tree->Branch("lepton3_pvassq", &lepton3_pvassq);

   tree->Branch("lepton4_pt",  &lepton4_pt);
   tree->Branch("lepton4_eta", &lepton4_eta);
   tree->Branch("lepton4_phi", &lepton4_phi);
   tree->Branch("lepton4_charge", &lepton4_charge);
   tree->Branch("lepton4_d0",  &lepton4_d0);
   tree->Branch("lepton4_dxy", &lepton4_dxy);
   tree->Branch("lepton4_dz",  &lepton4_dz);
   tree->Branch("lepton4_iso03particle", &lepton4_iso03particle); 
   tree->Branch("lepton4_iso04particle", &lepton4_iso04particle);
   tree->Branch("lepton4_iso03hadron",   &lepton4_iso03hadron); //charged hadron iso 
   tree->Branch("lepton4_iso04hadron",   &lepton4_iso04hadron);
   tree->Branch("lepton4_iso03neutralHadron", &lepton4_iso03neutralHadron); //neutral hadron iso 
   tree->Branch("lepton4_iso03photon", &lepton4_iso03photon); //photon iso 
   tree->Branch("lepton4_iso03PU", &lepton4_iso03PU); //PU iso 
   tree->Branch("lepton4_impactParameterSignificance", &lepton4_impactParameterSignificance);
   tree->Branch("lepton4_isLooseMuon", &lepton4_isLooseMuon);    // these are only for muons 
   tree->Branch("lepton4_isMediumMuon", &lepton4_isMediumMuon);
   tree->Branch("lepton4_isSoftMuon",  &lepton4_isSoftMuon);     // these are only for muons  
   tree->Branch("lepton4_isTightMuon", &lepton4_isTightMuon);    // these are only for muons  
   tree->Branch("lepton4_isPFMuon", &lepton4_isPFMuon);          // these are only for muons  
   tree->Branch("lepton4_isGlobalMuon", &lepton4_isGlobalMuon);
   tree->Branch("lepton4_isTrackerMuon", &lepton4_isTrackerMuon);
   tree->Branch("lepton4_fromPV", &lepton4_fromPV);
   tree->Branch("lepton4_pvassq", &lepton4_pvassq);

   tree->Branch("dimuon1mass", &dimuon1mass);
   tree->Branch("dimuon2mass", &dimuon2mass);
   tree->Branch("dimuon1pt", &dimuon1pt);
   tree->Branch("dimuon2pt", &dimuon2pt);
   tree->Branch("dimuon1eta", &dimuon1eta);
   tree->Branch("dimuon2eta", &dimuon2eta);
   tree->Branch("dimuon1phi", &dimuon1phi);
   tree->Branch("dimuon2phi", &dimuon2phi);

   tree->Branch("dimuon1vtx", &dimuon1vtx);
   tree->Branch("dimuon2vtx", &dimuon2vtx);
   tree->Branch("dimuon1vtx_xpos", &dimuon1vtx_xpos);
   tree->Branch("dimuon2vtx_xpos", &dimuon2vtx_xpos);
   tree->Branch("dimuon1vtx_ypos", &dimuon1vtx_ypos);
   tree->Branch("dimuon2vtx_ypos", &dimuon2vtx_ypos);
   tree->Branch("dimuon1vtx_zpos", &dimuon1vtx_zpos);
   tree->Branch("dimuon2vtx_zpos", &dimuon2vtx_zpos);
   tree->Branch("dimuon1vtx_xposError", &dimuon1vtx_xposError);
   tree->Branch("dimuon2vtx_xposError", &dimuon2vtx_xposError);
   tree->Branch("dimuon1vtx_yposError", &dimuon1vtx_yposError);
   tree->Branch("dimuon2vtx_yposError", &dimuon2vtx_yposError);
   tree->Branch("dimuon1vtx_zposError", &dimuon1vtx_zposError);
   tree->Branch("dimuon2vtx_zposError", &dimuon2vtx_zposError);

   tree->Branch("dimuon1lxy", &dimuon1lxy);
   tree->Branch("dimuon2lxy", &dimuon2lxy);
   tree->Branch("dimuon1lxysig", &dimuon1lxysig);
   tree->Branch("dimuon2lxysig", &dimuon2lxysig);
   tree->Branch("dimuon1lxyctauPV", &dimuon1lxyctauPV);
   tree->Branch("dimuon2lxyctauPV", &dimuon2lxyctauPV);

   tree->Branch("sixmuonvtx", &sixmuonvtx);

   tree->Branch("triggerlist", &triggerlist);
   
   tree->Branch("save_event_count", &save_event_count);
   tree->Branch("numberOfVertices", &numberOfVertices);
   tree->Branch("zOfVertices",    &zOfVertices);
   tree->Branch("zOfVerticesError",    &zOfVerticesError);

   tree->Branch("pair_12_34_56", &pair_12_34_56);
   tree->Branch("pair_13_24_56", &pair_13_24_56);
   tree->Branch("pair_14_23_56", &pair_14_23_56);
   
   tree->Branch("pair_12_34_ZOnly", &pair_12_34_ZOnly);
   tree->Branch("pair_13_24_ZOnly", &pair_13_24_ZOnly);
   tree->Branch("pair_14_23_ZOnly", &pair_14_23_ZOnly);

   tree->Branch("PVx", &PVx);
   tree->Branch("PVy", &PVy);
   tree->Branch("PVz", &PVz);
   
   tree->Branch("quadComesFromEvtWithThisManyMu", &quadComesFromEvtWithThisManyMu);
   
   tree->Branch("inv4MuMass", &inv4MuMass);
   tree->Branch("big4MuVtx", &big4MuVtx);
   tree->Branch("big4MuEta", &big4MuEta);
   tree->Branch("big4MuPhi", &big4MuPhi);
   tree->Branch("big4MuPt", &big4MuPt);
   
   //trigger matching variable
   tree->Branch("quadHasHowManyTrigMatches", &quadHasHowManyTrigMatches);
   
   tree->Branch("sum4LepCharges", &sum4LepCharges);
   
//   tree->Branch("TrkWeightsRecoVtxTrks", &TrkWeightsRecoVtxTrks);
   
   tree->Branch("lepton1_validHits", &lepton1_validHits);
   tree->Branch("lepton2_validHits", &lepton2_validHits);
   tree->Branch("lepton3_validHits", &lepton3_validHits);
   tree->Branch("lepton4_validHits", &lepton4_validHits);

   treemc = fs->make<TTree>("treemc", "treemc");
   treemc->Branch("truth_Zmuon_pt",  &truth_Zmuon_pt);
   treemc->Branch("truth_Zmuon_eta", &truth_Zmuon_eta);
   treemc->Branch("truth_Zmuon_phi", &truth_Zmuon_phi);
   treemc->Branch("truth_Z_pt",   &truth_Z_pt);
   treemc->Branch("truth_Z_eta",  &truth_Z_eta);
   treemc->Branch("truth_Z_phi",  &truth_Z_phi);
   treemc->Branch("truth_Z_mass", &truth_Z_mass);
   treemc->Branch("truth_Z_pdgid",&truth_Z_pdgid);
   treemc->Branch("truth_Upsimuon_pt",  &truth_Upsimuon_pt);
   treemc->Branch("truth_Upsimuon_eta", &truth_Upsimuon_eta);
   treemc->Branch("truth_Upsimuon_phi", &truth_Upsimuon_phi);
   treemc->Branch("truth_Upsi_pt",   &truth_Upsi_pt);
   treemc->Branch("truth_Upsi_eta",  &truth_Upsi_eta);
   treemc->Branch("truth_Upsi_phi",  &truth_Upsi_phi);
   treemc->Branch("truth_Upsi_mass", &truth_Upsi_mass);
   treemc->Branch("truth_Upsi_pdgid", &truth_Upsi_pdgid);
   
   treemc->Branch("truth_Upsi2muon_pt",  &truth_Upsi2muon_pt);
   treemc->Branch("truth_Upsi2muon_eta", &truth_Upsi2muon_eta);
   treemc->Branch("truth_Upsi2muon_phi", &truth_Upsi2muon_phi);
   treemc->Branch("truth_Upsi2_pt",   &truth_Upsi2_pt);
   treemc->Branch("truth_Upsi2_eta",  &truth_Upsi2_eta);
   treemc->Branch("truth_Upsi2_phi",  &truth_Upsi2_phi);
   treemc->Branch("truth_Upsi2_mass", &truth_Upsi2_mass);
   treemc->Branch("truth_Upsi2_pdgid", &truth_Upsi2_pdgid);
   
   treemc->Branch("truth_Upsi3muon_pt",  &truth_Upsi3muon_pt);
   treemc->Branch("truth_Upsi3muon_eta", &truth_Upsi3muon_eta);
   treemc->Branch("truth_Upsi3muon_phi", &truth_Upsi3muon_phi);
   treemc->Branch("truth_Upsi3_pt",   &truth_Upsi3_pt);
   treemc->Branch("truth_Upsi3_eta",  &truth_Upsi3_eta);
   treemc->Branch("truth_Upsi3_phi",  &truth_Upsi3_phi);
   treemc->Branch("truth_Upsi3_mass", &truth_Upsi3_mass);
   treemc->Branch("truth_Upsi3_pdgid", &truth_Upsi3_pdgid);
   
   treemc->Branch("truth_Chib0_1P_UPSI_muon_pt", &truth_Chib0_1P_UPSI_muon_pt);
   treemc->Branch("truth_Chib0_1P_UPSI_muon_eta", &truth_Chib0_1P_UPSI_muon_eta);
   treemc->Branch("truth_Chib0_1P_UPSI_muon_phi", &truth_Chib0_1P_UPSI_muon_phi);
   
   treemc->Branch("truth_Chib0_1P_pt", &truth_Chib0_1P_pt);
   treemc->Branch("truth_Chib0_1P_eta", &truth_Chib0_1P_eta);
   treemc->Branch("truth_Chib0_1P_phi", &truth_Chib0_1P_phi);
   treemc->Branch("truth_Chib0_1P_pdgid", &truth_Chib0_1P_pdgid);
   treemc->Branch("truth_Chib0_1P_mass", &truth_Chib0_1P_mass);
   
   treemc->Branch("truth_Chib0_1P_UPSI_pt", &truth_Chib0_1P_UPSI_pt);
   treemc->Branch("truth_Chib0_1P_UPSI_phi", &truth_Chib0_1P_UPSI_phi);
   treemc->Branch("truth_Chib0_1P_UPSI_eta", &truth_Chib0_1P_UPSI_eta);
   treemc->Branch("truth_Chib0_1P_UPSI_pdgid", &truth_Chib0_1P_UPSI_pdgid);
   treemc->Branch("truth_Chib0_1P_UPSI_mass", &truth_Chib0_1P_UPSI_mass);
   
   treemc->Branch("truth_Chib0_1P_photon_pt", &truth_Chib0_1P_photon_pt);
   treemc->Branch("truth_Chib0_1P_photon_phi", &truth_Chib0_1P_photon_phi);
   treemc->Branch("truth_Chib0_1P_photon_eta", &truth_Chib0_1P_photon_eta);
   treemc->Branch("truth_Chib0_1P_photon_pdgid", &truth_Chib0_1P_photon_pdgid);
   
   treemc->Branch("truth_Chib1_1P_UPSI_muon_pt", &truth_Chib1_1P_UPSI_muon_pt);
   treemc->Branch("truth_Chib1_1P_UPSI_muon_eta", &truth_Chib1_1P_UPSI_muon_eta);
   treemc->Branch("truth_Chib1_1P_UPSI_muon_phi", &truth_Chib1_1P_UPSI_muon_phi);
   
   treemc->Branch("truth_Chib1_1P_pt", &truth_Chib1_1P_pt);
   treemc->Branch("truth_Chib1_1P_eta", &truth_Chib1_1P_eta);
   treemc->Branch("truth_Chib1_1P_phi", &truth_Chib1_1P_phi);
   treemc->Branch("truth_Chib1_1P_pdgid", &truth_Chib1_1P_pdgid);
   treemc->Branch("truth_Chib1_1P_mass", &truth_Chib1_1P_mass);
   
   treemc->Branch("truth_Chib1_1P_UPSI_pt", &truth_Chib1_1P_UPSI_pt);
   treemc->Branch("truth_Chib1_1P_UPSI_eta", &truth_Chib1_1P_UPSI_eta);
   treemc->Branch("truth_Chib1_1P_UPSI_phi", &truth_Chib1_1P_UPSI_phi);
   treemc->Branch("truth_Chib1_1P_UPSI_pdgid", &truth_Chib1_1P_UPSI_pdgid);
   treemc->Branch("truth_Chib1_1P_UPSI_mass", &truth_Chib1_1P_UPSI_mass);
   
   treemc->Branch("truth_Chib1_1P_photon_pt", &truth_Chib1_1P_photon_pt);
   treemc->Branch("truth_Chib1_1P_photon_eta", &truth_Chib1_1P_photon_eta);
   treemc->Branch("truth_Chib1_1P_photon_phi", &truth_Chib1_1P_photon_phi);
   treemc->Branch("truth_Chib1_1P_photon_pdgid", &truth_Chib1_1P_photon_pdgid);
   
   treemc->Branch("truth_Chib2_1P_UPSI_muon_pt", &truth_Chib2_1P_UPSI_muon_pt);
   treemc->Branch("truth_Chib2_1P_UPSI_muon_eta", &truth_Chib2_1P_UPSI_muon_eta);
   treemc->Branch("truth_Chib2_1P_UPSI_muon_phi", &truth_Chib2_1P_UPSI_muon_phi);
   
   treemc->Branch("truth_Chib2_1P_pt", &truth_Chib2_1P_pt);
   treemc->Branch("truth_Chib2_1P_eta", &truth_Chib2_1P_eta);
   treemc->Branch("truth_Chib2_1P_phi", &truth_Chib2_1P_phi);
   treemc->Branch("truth_Chib2_1P_pdgid", &truth_Chib2_1P_pdgid);
   treemc->Branch("truth_Chib2_1P_mass", &truth_Chib2_1P_mass);
   
   treemc->Branch("truth_Chib2_1P_UPSI_pt", &truth_Chib2_1P_UPSI_pt);
   treemc->Branch("truth_Chib2_1P_UPSI_eta", &truth_Chib2_1P_UPSI_eta);
   treemc->Branch("truth_Chib2_1P_UPSI_phi", &truth_Chib2_1P_UPSI_phi);
   treemc->Branch("truth_Chib2_1P_UPSI_pdgid", &truth_Chib2_1P_UPSI_pdgid);
   treemc->Branch("truth_Chib2_1P_UPSI_mass", &truth_Chib2_1P_UPSI_mass);
   
   treemc->Branch("truth_Chib2_1P_photon_pt", &truth_Chib2_1P_photon_pt);
   treemc->Branch("truth_Chib2_1P_photon_eta", &truth_Chib2_1P_photon_eta);
   treemc->Branch("truth_Chib2_1P_photon_phi", &truth_Chib2_1P_photon_phi);
   treemc->Branch("truth_Chib2_1P_photon_pdgid", &truth_Chib2_1P_photon_pdgid);
   
   treemc->Branch("loop_enter_check", &loop_enter_check);
   treemc->Branch("mc_event_number", &mc_event_number);
   treemc->Branch("mc_run_number", &mc_run_number);
   treemc->Branch("mc_lumi_section", &mc_lumi_section);
   
   treemc->Branch("eventHasZUpsiNTo4Mu_Count", &eventHasZUpsiNTo4Mu_Count);

}

void ZmuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<std::vector<pat::Muon>> muons;
   iEvent.getByToken(muonsToken_, muons);

   edm::ESHandle<TransientTrackBuilder> builder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

//Rho, the median energy desnity //

 //  edm::Handle< double > rhoH;
//   iEvent.getByToken(rhoToken_,rhoH);
// //  float rho=*rhoH;
//   rho.clear();
//   rho.push_back(*rhoH);
// *******
// TRIGGER
// *******
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_, triggerBits);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);

//   edm::Handle< pat::TriggerEvent > triggerEvent;
//   iEvent.getByToken( triggerEventToken_, triggerEvent );
//   const pat::helper::TriggerMatchHelper matchHelper;

// **************************************************
// FLAG FOR ASSOCIATED PRODUCTION OR RESONANCE SEARCH
// **************************************************
  // bool mc = true;
  
//  bool isMC = false; //testing 2018A SingleMu data
//   isMC = iConfig.getParameter<bool>("isMC");
//Recall that we get MC from the cfg setting, see:  isMC = iConfig.getParameter<bool>("isMC"); line above

   triggerlist.clear();
   save_event_count.clear();
   dimuon1mass.clear(); dimuon2mass.clear();
   dimuon1pt.clear(); dimuon2pt.clear();    
   dimuon1eta.clear(); dimuon2eta.clear();  
   dimuon1phi.clear(); dimuon2phi.clear();  
   dimuon1vtx.clear(); dimuon2vtx.clear();  
   dimuon1vtx_xpos.clear(); dimuon2vtx_xpos.clear(); 
   dimuon1vtx_ypos.clear(); dimuon2vtx_ypos.clear(); 
   dimuon1vtx_zpos.clear(); dimuon2vtx_zpos.clear(); 
   dimuon1vtx_xposError.clear(); dimuon2vtx_xposError.clear();
   dimuon1vtx_yposError.clear(); dimuon2vtx_yposError.clear();
   dimuon1vtx_zposError.clear(); dimuon2vtx_zposError.clear();
   sixmuonvtx.clear();

   dimuon1lxy.clear();       dimuon2lxy.clear();      
   dimuon1lxysig.clear();    dimuon2lxysig.clear();   
   dimuon1lxyctauPV.clear(); dimuon2lxyctauPV.clear();

   numberOfVertices.clear();
   zOfVerticesError.clear();
   zOfVertices.clear();

   PVx.clear();
   PVy.clear();
   PVz.clear();

   run_number.clear(); event_number.clear(); lumi_section.clear(); bunch_crossing.clear(); orbit_number.clear();
   
   quadHasHowManyTrigMatches.clear(); // branch passes a basic sanity check when  filled  this way
   quadComesFromEvtWithThisManyMu.clear();
   sum4LepCharges.clear();
   flagZplusY.clear();
   flagZOnly.clear();
   
   lepton1_pt                         .clear();  lepton2_pt                         .clear();
   lepton1_eta                        .clear();  lepton2_eta                        .clear();
   lepton1_phi                        .clear();  lepton2_phi                        .clear();
   lepton1_charge                     .clear();  lepton2_charge                     .clear();
   lepton1_dxy                        .clear();  lepton2_dxy                        .clear();
   lepton1_dz                         .clear();  lepton2_dz                         .clear();
   lepton1_d0                         .clear();  lepton2_d0                         .clear();
   lepton1_iso03particle              .clear();  lepton2_iso03particle              .clear();
   lepton1_iso03hadron                .clear();  lepton2_iso03hadron                .clear();
   lepton1_iso03neutralHadron         .clear();  lepton2_iso03neutralHadron         .clear();
   lepton1_iso03photon                .clear();  lepton2_iso03photon                .clear();
   lepton1_iso03PU                    .clear();  lepton2_iso03PU                    .clear();
   lepton1_iso04hadron                .clear();  lepton2_iso04hadron                .clear();
   lepton1_iso04particle              .clear();  lepton2_iso04particle              .clear();
   lepton1_impactParameterSignificance.clear();  lepton2_impactParameterSignificance.clear();
   lepton1_isLooseMuon                .clear();  lepton2_isLooseMuon                .clear();   
   lepton1_isMediumMuon               .clear();  lepton2_isMediumMuon               .clear();
   lepton1_isSoftMuon                 .clear();  lepton2_isSoftMuon                 .clear();        
   lepton1_isTightMuon                .clear();  lepton2_isTightMuon                .clear();
   lepton1_isPFMuon                   .clear();  lepton2_isPFMuon                   .clear();
   lepton1_isGlobalMuon               .clear();  lepton2_isGlobalMuon               .clear();
   lepton1_isTrackerMuon              .clear();  lepton2_isTrackerMuon              .clear();
   lepton1_deltaEta                   .clear();  lepton2_deltaEta                   .clear(); 
   lepton1_deltaPhi                   .clear();  lepton2_deltaPhi                   .clear(); 
   lepton1_sigmaiEtaEta               .clear();  lepton2_sigmaiEtaEta               .clear(); 
   lepton1_HoverE                     .clear();  lepton2_HoverE                     .clear(); 
   lepton1_OoEmOoP                    .clear();  lepton2_OoEmOoP                    .clear(); 
   lepton1_missingHits                .clear();  lepton2_missingHits                .clear(); 
   lepton1_fromPV                     .clear();  lepton2_fromPV                     .clear();
   lepton1_pvassq                     .clear();  lepton2_pvassq                     .clear();
   
   lepton3_pt                         .clear();  lepton4_pt                         .clear();
   lepton3_eta                        .clear();  lepton4_eta                        .clear();
   lepton3_phi                        .clear();  lepton4_phi                        .clear();
   lepton3_charge                     .clear();  lepton4_charge                     .clear();
   lepton3_dxy                        .clear();  lepton4_dxy                        .clear();
   lepton3_dz                         .clear();  lepton4_dz                         .clear();
   lepton3_d0                         .clear();  lepton4_d0                         .clear();
   lepton3_iso03particle              .clear();  lepton4_iso03particle              .clear();
   lepton3_iso03hadron                .clear();  lepton4_iso03hadron                .clear();
   lepton3_iso03neutralHadron         .clear();  lepton4_iso03neutralHadron         .clear();
   lepton3_iso03photon                .clear();  lepton4_iso03photon                .clear();
   lepton3_iso03PU                    .clear();  lepton4_iso03PU                    .clear();
   lepton3_iso04hadron                .clear();  lepton4_iso04hadron                .clear();
   lepton3_iso04particle              .clear();  lepton4_iso04particle              .clear();
   lepton3_impactParameterSignificance.clear();  lepton4_impactParameterSignificance.clear();
   lepton3_isLooseMuon                .clear();  lepton4_isLooseMuon                .clear(); 
   lepton3_isLooseMuon                .clear();  lepton4_isMediumMuon               .clear();   
   lepton3_isSoftMuon                 .clear();  lepton4_isSoftMuon                 .clear();        
   lepton3_isTightMuon                .clear();  lepton4_isTightMuon                .clear();
   lepton3_isPFMuon                   .clear();  lepton4_isPFMuon                   .clear();
   lepton3_isGlobalMuon               .clear();  lepton4_isGlobalMuon               .clear();
   lepton3_isTrackerMuon              .clear();  lepton4_isTrackerMuon              .clear();
   lepton3_deltaEta                   .clear();  lepton4_deltaEta                   .clear(); 
   lepton3_deltaPhi                   .clear();  lepton4_deltaPhi                   .clear(); 
   lepton3_sigmaiEtaEta               .clear();  lepton4_sigmaiEtaEta               .clear(); 
   lepton3_HoverE                     .clear();  lepton4_HoverE                     .clear(); 
   lepton3_OoEmOoP                    .clear();  lepton4_OoEmOoP                    .clear(); 
   lepton3_missingHits                .clear();  lepton4_missingHits                .clear();
   lepton3_fromPV                     .clear();  lepton4_fromPV                     .clear();
   lepton3_pvassq                     .clear();  lepton4_pvassq                     .clear(); 

   pair_12_34_56.clear();    
   pair_13_24_56.clear();    
   pair_14_23_56.clear();
   
   pair_12_34_ZOnly.clear();
   pair_13_24_ZOnly.clear();
   pair_14_23_ZOnly.clear();
   
   inv4MuMass.clear();
   big4MuVtx.clear();
   big4MuEta.clear();
   big4MuPhi.clear();
   big4MuPt.clear();
   
   lepton1_validHits.clear(); lepton2_validHits.clear(); lepton3_validHits.clear(); lepton4_validHits.clear();

   truth_Zmuon_pt.clear(); 
   truth_Zmuon_eta.clear();
   truth_Zmuon_phi.clear();
   truth_Z_pt.clear(); truth_Z_eta.clear(); truth_Z_phi.clear(); truth_Z_mass.clear(); truth_Z_pdgid.clear();
   truth_Upsimuon_pt.clear(); 
   truth_Upsimuon_eta.clear();
   truth_Upsimuon_phi.clear();
   truth_Upsi_pt.clear(); truth_Upsi_eta.clear(); truth_Upsi_phi.clear(); truth_Upsi_mass.clear();
   
   truth_Upsi2muon_pt.clear(); 
   truth_Upsi2muon_eta.clear();
   truth_Upsi2muon_phi.clear();
   truth_Upsi2_pt.clear(); truth_Upsi2_eta.clear(); truth_Upsi2_phi.clear(); truth_Upsi2_mass.clear();
   truth_Upsi2_pdgid.clear();
   
   truth_Upsi3muon_pt.clear(); 
   truth_Upsi3muon_eta.clear();
   truth_Upsi3muon_phi.clear();
   truth_Upsi3_pt.clear(); truth_Upsi3_eta.clear(); truth_Upsi3_phi.clear(); truth_Upsi3_mass.clear();
   truth_Upsi3_pdgid.clear();
   
   truth_Chib0_1P_UPSI_muon_pt.clear();
   truth_Chib0_1P_UPSI_muon_eta.clear();
   truth_Chib0_1P_UPSI_muon_phi.clear();
   
   truth_Chib0_1P_pt.clear();
   truth_Chib0_1P_eta.clear();
   truth_Chib0_1P_phi.clear();
   truth_Chib0_1P_pdgid.clear();
   truth_Chib0_1P_mass.clear();
   
   truth_Chib0_1P_UPSI_pt.clear();
   truth_Chib0_1P_UPSI_phi.clear();
   truth_Chib0_1P_UPSI_eta.clear();
   truth_Chib0_1P_UPSI_pdgid.clear();
   truth_Chib0_1P_UPSI_mass.clear();
   
   truth_Chib0_1P_photon_pt.clear();
   truth_Chib0_1P_photon_phi.clear();
   truth_Chib0_1P_photon_eta.clear();
   truth_Chib0_1P_photon_pdgid.clear();
   
   truth_Chib1_1P_UPSI_muon_pt.clear();
   truth_Chib1_1P_UPSI_muon_phi.clear();
   truth_Chib1_1P_UPSI_muon_eta.clear();
   
   truth_Chib1_1P_pt.clear();
   truth_Chib1_1P_eta.clear();
   truth_Chib1_1P_phi.clear();
   truth_Chib1_1P_pdgid.clear();
   truth_Chib1_1P_mass.clear();
   
   truth_Chib1_1P_UPSI_pt.clear();
   truth_Chib1_1P_UPSI_phi.clear();
   truth_Chib1_1P_UPSI_eta.clear();
   truth_Chib1_1P_UPSI_pdgid.clear(); 
   truth_Chib1_1P_UPSI_mass.clear();  
   
   truth_Chib1_1P_photon_pt.clear();
   truth_Chib1_1P_photon_eta.clear();
   truth_Chib1_1P_photon_phi.clear();
   truth_Chib1_1P_photon_pdgid.clear();
   
   truth_Chib2_1P_UPSI_muon_pt.clear();
   truth_Chib2_1P_UPSI_muon_eta.clear();
   truth_Chib2_1P_UPSI_muon_phi.clear();
   
   truth_Chib2_1P_pt.clear();
   truth_Chib2_1P_eta.clear();
   truth_Chib2_1P_phi.clear();
   truth_Chib2_1P_pdgid.clear();
   truth_Chib2_1P_mass.clear();
   
   truth_Chib2_1P_UPSI_pt.clear();
   truth_Chib2_1P_UPSI_eta.clear();
   truth_Chib2_1P_UPSI_phi.clear();
   truth_Chib2_1P_UPSI_pdgid.clear();
   truth_Chib2_1P_UPSI_mass.clear();
   
   truth_Chib2_1P_photon_pt.clear();
   truth_Chib2_1P_photon_eta.clear();
   truth_Chib2_1P_photon_phi.clear();
   truth_Chib2_1P_photon_pdgid.clear();
    
   
   mc_event_number.clear();
   mc_run_number.clear();
   mc_lumi_section.clear();

   truth_Upsi_pdgid.clear();  
   loop_enter_check.clear();
   rho.clear();
   
   eventHasZUpsiNTo4Mu_Count.clear();
   
   
 //  TrkWeightsRecoVtxTrks.clear();
   
 //  std::cout << "Am here at check 0" << std::endl;

////   bool nontriggeredevent = true;
//   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
//    if (triggerBits->accept(i)) {
//      if (names.triggerName(i).find("HLT")<100 && names.triggerName(i).find("Tau")>100 && names.triggerName(i).find("MET")>100 && names.triggerName(i).find("Jet")>100 && names.triggerName(i).find("isplaced")>100 && names.triggerName(i).find("SameSign")>100 && names.triggerName(i).find("Photon")>100) 
//        if (names.triggerName(i).find("Mu")<100 || names.triggerName(i).find("mu")<100) {
//          triggerlist.push_back(names.triggerName(i));
//  //      if (names.triggerName(i).find("Ele")<100 || names.triggerName(i).find("ele")<100)
////          nontriggeredevent = false;
//        }
//    }

   bool event_fails_trigger = true; //defaults to true 
   
   //bools for 2016 triggers
   bool singleMu2016Trig1Fired = false;
   bool singleMu2016Trig2Fired = false;
   bool doubleMu2016Trig1Fired = false;
   bool doubleMu2016Trig2Fired = false;
   bool tripleMu2016Trig1Fired = false; 
   
   //bools for 2017 triggers
   bool singleMu2017Trig1Fired = false;
   bool doubleMu2017Trig1Fired = false; 
   bool tripleMu2017Trig1Fired = false;
   bool tripleMu2017Trig2Fired = false;
   
   //bools for 2018 triggers
   bool singleMu2018Trig1Fired = false;
   bool doubleMu2018Trig1Fired = false; 
   bool doubleMu2018Trig2Fired = false; 
   bool tripleMu2018Trig1Fired = false;
   bool tripleMu2018Trig2Fired = false; 
   
   
   
 //   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
//     if (triggerBits->accept(i)) {
//       triggerlist.push_back(names.triggerName(i));
//       std::string str (names.triggerName(i));
//       std::cout << "str is :" << str << std::endl;
//       std::string str2 ("Ele");
//       std::string str3 ("Mu"); //confirmed that Mu is the right thing to search even for the DoubleMuon Dataset we will be using, see: https://twiki.cern.ch/twiki/bin/view/CMS/HLTPathsRunIIList
//       std::size_t foundEle = str.find(str2);
//       std::size_t foundMu  = str.find(str3);
//       if (foundEle!=std::string::npos && foundMu==std::string::npos){ //Do this as Stefanos did, full
//         event_fails_trigger = true;
//         std::cout << "event_fails_trigger is True" << std::endl;
//       }
//       
//       if (foundEle!=std::string::npos && foundMu!=std::string::npos){
//      
//         event_fails_trigger = false;
//         std::cout << "event_fails_trigger is False" << std::endl; 
//       }
//     }
// 
// 
//    } //trigger requirements 
//   if (event_fails_trigger) {
//      for (int itrig = 0; itrig < (int)triggerlist.size(); itrig++)
//        std::cout << triggerlist.at(itrig) << std::endl;
//      std::cout << "--------------\n";
//   }


//The above stuff is not the best way to deal with it, let's try this instead

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
     if (triggerBits->accept(i)) {
       triggerlist.push_back(names.triggerName(i));
      }
   }
   
  //  for (unsigned int i= 0, n = triggerlist.size(); i < n; i++) {
//        std::cout << triggerlist.at(i) << std::endl;
//        std::string str (triggerlist.at(i));
//        std::cout << "str is:" << str << std::endl;
//        std::string str2 ("Ele");
//        std::string str3 ("Mu");
//        std::size_t foundEle = str.find(str2);
//        std::size_t foundMu  = str.find(str3);
//        if (foundEle==std::string::npos && foundMu!=std::string::npos){ //avoid mixed e/mu triggers
//           event_fails_trigger=false;
//           std::cout << "Event passed trigger!" << std::endl;
//           break; //as soon as we found that it passed, get out of the loop!
//        
//        } 
//    
//    }
   
    for (unsigned int i= 0, n = triggerlist.size(); i < n; i++) {
       std::cout << triggerlist.at(i) << std::endl;
       std::string str (triggerlist.at(i));
       std::cout << "str is:" << str << std::endl;
       if (triggerYear == 2016){
         std::string str2 ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"); //call this DoubleMu2016Trig1
         std::string str3 ("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"); //call this DoubleMu2016Trig2
         std::string str4 ("HLT_IsoMu24_v"); //call this SingleMu2016Trig1
         std::string str5 ("HLT_IsoTkMu24_v"); //call this SingleMu2016Trig2
         std::string str5a ("HLT_TripleMu_12_10_5_v"); //call this TripleMu2016Trig1
         
         std::size_t foundDoubleMu2016Trig1 = str.find(str2);
         std::size_t foundDoubleMu2016Trig2 = str.find(str3);
         std::size_t foundSingleMu2016Trig1 = str.find(str4);
         std::size_t foundSingleMu2016Trig2 = str.find(str5);
         std::size_t foundTripleMu2016Trig1 = str.find(str5a);
   
         
         if (foundSingleMu2016Trig1 != std::string::npos){
           singleMu2016Trig1Fired = true;
           std::cout << "singleMu2016Trig1Fired:  "  << singleMu2016Trig1Fired << std::endl; 
         }
         
         if (foundSingleMu2016Trig2 != std::string::npos){
           singleMu2016Trig2Fired = true;
           std::cout << "singleMu2016Trig2Fired:  " << singleMu2016Trig2Fired << std::endl; 
         }
         
         if (foundDoubleMu2016Trig1 != std::string::npos){
           doubleMu2016Trig1Fired = true;
           std::cout << "doubleMu2016Trig1Fired:  " << doubleMu2016Trig1Fired << std::endl; 
         }
         
         if (foundDoubleMu2016Trig2 != std::string::npos){
           doubleMu2016Trig2Fired = true;
           std::cout << "doubleMu2016Trig2Fired:  " << doubleMu2016Trig2Fired << std::endl; 
         }
         
         if (foundTripleMu2016Trig1 != std::string::npos){
           tripleMu2016Trig1Fired = true;
           std::cout << "tripleMu2016Trig1Fired:  " << tripleMu2016Trig1Fired << std::endl; 
         }
         
//         if (foundDoubleMu2016Trig1 !=std::string::npos || foundDoubleMu2016Trig2 !=std::string::npos){ //OR of both DoubleMu trigs that we are considering
         if (foundSingleMu2016Trig1 != std::string::npos ||  foundSingleMu2016Trig2 != std::string::npos || foundDoubleMu2016Trig1 !=std::string::npos || foundDoubleMu2016Trig2 !=std::string::npos || foundTripleMu2016Trig1 != std::string::npos){
           event_fails_trigger=false;
           std::cout << "Event passed 2016 trigger!" << std::endl;
           break; //as soon as we found that it passed, get out of the loop!
         }
       } 
       if (triggerYear == 2017){
         std::string str6 ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"); //call this DoubleMu2017Trig1
         std::string str7 ("HLT_IsoMu27_v"); //call this SingleMu2017Trig1
         std::string str7a ("HLT_TripleMu_12_10_5_v"); //call this TripleMu2017Trig1
         std::string str7b ("HLT_TripleMu_10_5_5_DZ_v"); //call this TripleMu2017Trig2 
         
         std::size_t foundDoubleMu2017Trig1 = str.find(str6);
         std::size_t foundSingleMu2017Trig1 = str.find(str7);
         std::size_t foundTripleMu2017Trig1 = str.find(str7a);
         std::size_t foundTripleMu2017Trig2 = str.find(str7b);
         
         if (foundSingleMu2017Trig1 != std::string::npos){
           singleMu2017Trig1Fired = true;
           std::cout << "singleMu2017Trig1Fired:  " << singleMu2017Trig1Fired << std::endl; 
         }
         
         if (foundDoubleMu2017Trig1 != std::string::npos){
           doubleMu2017Trig1Fired = true;
           std::cout << "doubleMu2017Trig1Fired:  " << doubleMu2017Trig1Fired << std::endl; 
         }
         
         if (foundTripleMu2017Trig1 != std::string::npos){
           tripleMu2017Trig1Fired = true;
           std::cout << "tripleMu2017Trig1Fired:  " << tripleMu2017Trig1Fired << std::endl;
         }
         
          if (foundTripleMu2017Trig2 != std::string::npos){
           tripleMu2017Trig2Fired = true;
           std::cout << "tripleMu2017Trig2Fired:  " << tripleMu2017Trig2Fired << std::endl;
         }
         
         if (foundSingleMu2017Trig1 != std::string::npos || foundDoubleMu2017Trig1 != std::string::npos || foundTripleMu2017Trig1 != std::string::npos || foundTripleMu2017Trig2 != std::string::npos){
           event_fails_trigger=false;
           std::cout << "Event passed 2017 trigger!" << std::endl;
           break; 
         }
       
       }
       if (triggerYear == 2018){
         std::string str8 ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"); //call this DoubleMu2018Trig1
         std::string str8a ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"); // call this DoubleMu2018Trig2
         std::string str9 ("HLT_IsoMu24_v"); // call this SingleMu2018Trig1
         std::string str9a ("HLT_TripleMu_12_10_5_v"); // call this TripleMu2018Trig1
         std::string str9b ("HLT_TripleMu_10_5_5_DZ_v"); // call this TripleMu2018Trig2
         
         std::size_t foundDoubleMu2018Trig1 = str.find(str8);
         std::size_t foundDoubleMu2018Trig2 = str.find(str8a);
         std::size_t foundSingleMu2018Trig1 = str.find(str9);
         std::size_t foundTripleMu2018Trig1 = str.find(str9a);
         std::size_t foundTripleMu2018Trig2 = str.find(str9b);
         
         if (foundSingleMu2018Trig1 != std::string::npos){
           singleMu2018Trig1Fired = true;
           std::cout << "singleMu2018Trig1Fired:  " << singleMu2018Trig1Fired << std::endl; 
         }
         
         if (foundDoubleMu2018Trig1 != std::string::npos){
           doubleMu2018Trig1Fired = true;
           std::cout << "doubleMu2018Trig1Fired:  " << doubleMu2018Trig1Fired << std::endl; 
         }
         
         if (foundDoubleMu2018Trig2 != std::string::npos){
           doubleMu2018Trig2Fired = true;
           std::cout << "doubleMu2018Trig2Fired:  " << doubleMu2018Trig2Fired << std::endl;
         }
         
         if (foundTripleMu2018Trig1 != std::string::npos){
           tripleMu2018Trig1Fired = true;
           std::cout << "tripleMu2018Trig1Fired:  " << tripleMu2018Trig1Fired << std::endl; 
         }
         
         if (foundTripleMu2018Trig2 != std::string::npos){
           tripleMu2018Trig2Fired = true;
           std::cout << "tripleMu2018Trig2Fired:  " << tripleMu2018Trig2Fired << std::endl; 
         }
         
 //        if (foundDoubleMu2018Trig1 != std::string::npos){
         if (foundSingleMu2018Trig1 != std::string::npos || foundDoubleMu2018Trig1 != std::string::npos || foundDoubleMu2018Trig2 != std::string::npos || foundTripleMu2018Trig1 != std::string::npos || foundTripleMu2018Trig2 != std::string::npos){
           event_fails_trigger = false;
           std::cout << "Event passed 2018 trigger!" << std::endl;
           break; 
         
         }
    
    }
   }

// ******
// VERTEX
// ******
   edm::Handle<std::vector<reco::Vertex>> reco_vertices;
   iEvent.getByToken(token_vertices, reco_vertices);

   bool first_vertex = true; //first_vertex defaults to true 
   reco::Vertex primary_vertex;

   numberOfVertices.push_back(reco_vertices->size());
   for (unsigned int vertex=0; vertex < reco_vertices->size(); ++vertex)  { //loop over vertices, do some selections
     if (    // Criteria copied from twiki (mythical)
         !((*reco_vertices)[vertex].isFake()) //vertex must be true 
         && ((*reco_vertices)[vertex].ndof() > 4) //vertex must have ndof >4
         && (fabs((*reco_vertices)[vertex].z()) <= 24.0) //dz of the vertex must be <= 24 cm
         && ((*reco_vertices)[vertex].position().Rho() <= 2.0)//rho of the vertex must be less than or equal to 2 cm 
        ) {
//         reco_vert.num++;
//         reco_vert.x.push_back( (*reco_vertices)[vertex].x() );
//         reco_vert.y.push_back( (*reco_vertices)[vertex].y() );
//         reco_vert.z.push_back( (*reco_vertices)[vertex].z() );
//         // Store first good vertex as "primary"
//         // vertex is ordered by highest pt, in descending order //This is basically true, see email from Michal Bluj dated 7 Oct. 2020
//         //std::cout << iEvent.id().event() << "vertex pT: " << sumPtSquared((*reco_vertices)[vertex]) << std::endl;
//         sumOfVertices.push_back(sumPtSquared((*reco_vertices)[vertex]));

          /// <--------- to be tested
          /// <--------- to be tested
          reco::Vertex it_vertex = (*reco_vertices)[vertex]; //if the vertex meets the requirements, go ahead and do this loop
 //         std::cout << "Got here" << std::endl;
 
 //Confirmed with tracking experts that track info is NOT stored in miniAOD  (see email from Wolfram Erdmann from 12 Oct. 2020 )
   //        double sum = 0.;
//           double pT;
//           //It appears trackSize is always 0 and so this loop never fires
//           std::cout <<  "it_vertex.tracksSize() " << it_vertex.tracksSize();
    //       for (reco::Vertex::trackRef_iterator it = it_vertex.tracks_begin(); it != it_vertex.tracks_end(); it++) {
// //            for (auto it = it_vertex.tracks_begin(); it != it_vertex.tracks_end(); it++){
//  //           std::cout << "Also got here" << std::endl;
//             pT = (**it).pt();
//             sum += pT; //sum up the pT of all the tracks associated with that vertex
//             std::cout << pT << " " << sum << std::endl;
//             iTrkWeightOfRecoVtxTrk = it_vertex.trackWeight(*it);
//             TrkWeightsRecoVtxTrks.push_back(iTrkWeightOfRecoVtxTrk);
//             std::cout << it_vertex.trackWeight(*it) << std::endl;
//           } 
//           
                  

           
          zOfVerticesError.push_back((*reco_vertices)[vertex].zError());
          zOfVertices.push_back((*reco_vertices)[vertex].z());
//          std::cout << "Test" << std::endl;
          /// <--------- to be tested
          /// <--------- to be tested


         if (first_vertex) { //this ensures that the first vertex that passes the criteria is taken as the primary_vertex, assuming that miniaod does actually order its vertices, this is all good 
           first_vertex = false;
           primary_vertex = (*reco_vertices)[vertex];
           PVx.push_back((*reco_vertices)[vertex].x());
           PVy.push_back((*reco_vertices)[vertex].y());
           PVz.push_back((*reco_vertices)[vertex].z());
         }
     }
   }

   bool doit = false;
   if (iEvent.run() == 319077 || iEvent.run() == 319659 || iEvent.run() == 324205 || iEvent.run() == 322348)
      doit = true;
      

    
  
   
   // bypass for non-global muons check
   doit = true; //QUESTION! so doit is always true, so what is the point of the code above where doit was contingent on the event belonging to run NNN

//**** Lines for trigger efficiency studies*******
//   event_fails_trigger = false; // Comment in this line for getting denominator number for trigger efficiency studies
 //  std::cout << "event_fails_trigger ultimately is: " << event_fails_trigger << std::endl;
   
   bool save_event = false; 
   
   //Cut values 
   double muon_mass = 0.1056583715; 
   double upsi_mass_low  = 8.; //selections as from AN 
   double upsi_mass_high = 12.; //was 11 in AN but I don't think this really matters so much, can tighten in next step 
   double Z_mass_low  = 66.;
   double Z_mass_high = 116.;
   
   //double avoid_JPsi_ZOnly = 4;
  // double avoid_upsi_ZOnly = 12;
   double lowerBound_ZOnly = 1; 
   double ZOnly_mass_low = 80;
   double ZOnly_mass_high = 100; 
//   int save_event_count = 0;

//quick loop on muons collection too determine number of "good" muons in the event, where "good" is defined as passing pT_cut and having isTrackerMuon == true
   int numOfGoodMuInEvtCount = 0;
   int numOfGoodMuPlusInEvtCount = 0;
   int numOfGoodMuMinusInEvtCount = 0;
   
   double my_pT_cut = 3;
   for (auto iM = muons->begin(); iM != muons->end(); iM++){
     if (iM->pt() > my_pT_cut && iM->isTrackerMuon() == true){
       numOfGoodMuInEvtCount++;
     }
     if (iM->pt() > my_pT_cut && iM->isTrackerMuon() == true && iM->charge() ==1){
       numOfGoodMuPlusInEvtCount++;
     }
     if (iM->pt() > my_pT_cut && iM->isTrackerMuon() == true && iM->charge() == -1){
       numOfGoodMuMinusInEvtCount++;
     }
   }
    histContainer_["numOfGoodMuInEvt"]->Fill(numOfGoodMuInEvtCount);
    histContainer_["numOfGoodMuPlusInEvt"]->Fill(numOfGoodMuPlusInEvtCount);
    histContainer_["numOfGoodMuMinusInEvt"]->Fill(numOfGoodMuMinusInEvtCount);
    histContainer2_["numOfGoodMuPlusVsGoodMuMinus"]->Fill(numOfGoodMuPlusInEvtCount, numOfGoodMuMinusInEvtCount);

// high pt muons - low pt muons
//confirmed that in miniaod, muons are ordered from high to low pT, see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPATDataFormats#PatMuon
//  if (!nontriggeredevent)
    if ((int)muons->size()<4)
       std::cout << "LESS THAN FOUR"<< std::endl;
    std::cout << ((int)muons->size()) << std::endl;
    if (event_fails_trigger == true)
       std::cout << "FAILS TRIGGER" << std::endl;
   if ((int)muons->size()>3 && !event_fails_trigger && doit) {  //leave doit as homage to S.L.
//     bool save_event = false; //coment out for test
    // std::cout << "HERE" << std::endl;
   //  int numOfGoodMu = 0; //finish this later 12 Jan. 2021
     for (auto iM1 = muons->begin(); iM1 != muons->end(); ++iM1) {
      for (auto iM2 = iM1+1; iM2 != muons->end(); ++iM2) {
        for (auto iM3 = iM2+1; iM3 != muons->end(); ++iM3) {
          for (auto iM4 = iM3+1; iM4 != muons->end(); ++iM4) {
             
                 
                double pT_cut = 3; // we cut looser here and then harder (pT > 4) in phase 2 code
                histContainer_["CutFlow"]->Fill(1); //how many candidates we have 
   //             std::cout << " <<<<<<<<<<<<< ------------------------- entered " << std::endl;
    
                if (iM1->charge() + iM2->charge() + iM3->charge() + iM4->charge() !=0){
                 
                  std::cout << "FAILED AT CUT 1" << std::endl;
                  histContainer_["CutFlow"]->Fill(2);
                 continue;
                }
                // if (fabs(iM1->muonBestTrack()->dz(primary_vertex.position()))>20. || fabs(iM2->muonBestTrack()->dz(primary_vertex.position()))>20.){
//                  histContainer_["CutFlow"]->Fill(3);
//                  std::cout << "FAILED AT CUT 2" << std::endl;
//                   continue;
//                 }
//                 if (fabs(iM1->muonBestTrack()->dxy(primary_vertex.position()))>1. || fabs(iM2->muonBestTrack()->dxy(primary_vertex.position()))>1.){
//                   histContainer_["CutFlow"]->Fill(4);
//                   std::cout << "FAILED AT CUT 3" <<std::endl;
//                   continue;
//                 }
                if (fabs(iM1->eta())>2.5 || fabs(iM2->eta())>2.5){
                  histContainer_["CutFlow"]->Fill(5);
                  std::cout << "FAILED AT CUT 4" << std::endl;
                  continue;
                }
                if (iM1->pt()<pT_cut || iM2->pt()<pT_cut){
                  histContainer_["CutFlow"]->Fill(6);
                  std::cout << "FAILED AT CUT 5" << std::endl;
                  continue;
                }
                if (iM1->isTrackerMuon()==false || iM2->isTrackerMuon()==false){ //See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis#Tracker_muons for description of tracker muons
                  histContainer_["CutFlow"]->Fill(7);
                  std::cout << "FAILED AT CUT 6" << std::endl;
                  continue;
                }
  //              if (iM1->isSoftMuon(primary_vertex)==false || iM2->isSoftMuon(primary_vertex)==false){ //QUESTION!: I wouldn't have known isSoftMuon needed to (primary_vertex) argument, is this something you just knew? //Apparently you get compilation errors if you forget it
  //               histContainer_["CutFlow"]->Fill(8);
  //               std::cout << "FAILED AT CUT 7"  << std::endl;
   //              continue;
    //            }
           //     if (iM1->isGlobalMuon()==false || iM2->isGlobalMuon()==false)
             //     continue;
    
                // if (fabs(iM3->muonBestTrack()->dz(primary_vertex.position()))>20. || fabs(iM4->muonBestTrack()->dz(primary_vertex.position()))>20.){
//                   histContainer_["CutFlow"]->Fill(9);
//                   std::cout << "FAILED AT CUT 8" << std::endl;
//                   continue;
//                 }
//                 if (fabs(iM3->muonBestTrack()->dxy(primary_vertex.position()))>1. || fabs(iM4->muonBestTrack()->dxy(primary_vertex.position()))>1.){
//                   histContainer_["CutFlow"]->Fill(10);
//                   std::cout << "FAILED AT CUT 9" << std::endl;
//                   continue;
//                 }
                  
                if (fabs(iM3->eta())>2.5 || fabs(iM4->eta())>2.5){
                  histContainer_["CutFlow"]->Fill(11);
                  std::cout << "FAILED AT CUT 10" << std::endl;
                continue;
                }
                 
                if (iM3->pt()<pT_cut || iM4->pt()<pT_cut){
                  histContainer_["CutFlow"]->Fill(12);
                  std::cout << "FAILED AT CUT 11" << std::endl;
                  continue;
                }
                if (iM3->isTrackerMuon()==false || iM4->isTrackerMuon()==false){
                  histContainer_["CutFlow"]->Fill(13);
                  std::cout << "FAILED AT CUT 12" << std::endl; 
                  continue;
                }
                
                
   //             if (iM3->isSoftMuon(primary_vertex)==false || iM4->isSoftMuon(primary_vertex)==false){
   //                histContainer_["CutFlow"]->Fill(14);
    //               std::cout << "FAILED AT CUT 13" << std::endl;
    //               continue;
     //           }
                 
//                if (iM3->isGlobalMuon()==false || iM4->isGlobalMuon()==false)
//                  continue;
               //These cuts (dz = 20, dxy =1 )  come from the definition of the soft muon. S.L. pretty sure that if WE remove this cut, it won't make a difference...Probably obsolete.
               
               
                
                
                
                math::PtEtaPhiMLorentzVector lepton1(iM1->pt(), iM1->eta(), iM1->phi(), muon_mass);
                math::PtEtaPhiMLorentzVector lepton2(iM2->pt(), iM2->eta(), iM2->phi(), muon_mass);
                math::PtEtaPhiMLorentzVector lepton3(iM3->pt(), iM3->eta(), iM3->phi(), muon_mass);
                math::PtEtaPhiMLorentzVector lepton4(iM4->pt(), iM4->eta(), iM4->phi(), muon_mass);
    
                double mass_1_2=0.; double mass_1_3=0.; double mass_1_4 = 0.; 
                double mass_2_3=0.; double mass_2_4=0.;  
                double mass_3_4=0.;  
                double pt_1_2=0.;  double pt_1_3=0.;  double pt_1_4 = 0.;  
                double pt_2_3=0.;  double pt_2_4=0.;  
                double pt_3_4=0.; 
                double eta_1_2=0.; double eta_1_3=0.; double eta_1_4 = 0.; 
                double eta_2_3=0.; double eta_2_4=0.; 
                double eta_3_4=0.;
                double phi_1_2=0.; double phi_1_3=0.; double phi_1_4 = 0.; 
                double phi_2_3=0.; double phi_2_4=0.; 
                double phi_3_4=0.; 
                
                double mass_1_2_3_4 = 0;
                double pt_1_2_3_4 = 0;
                double eta_1_2_3_4 = 0;
                double phi_1_2_3_4 = 0;
                
                //build all the possible pairs of oppositely charged mu
                if (iM1->charge() + iM2->charge() == 0) { mass_1_2 = (lepton1 + lepton2).mass(); pt_1_2 = (lepton1 + lepton2).pt(); eta_1_2 = (lepton1 + lepton2).eta(); eta_1_2 = (lepton1 + lepton2).eta(); phi_1_2 = (lepton1 + lepton2).phi();}
                if (iM1->charge() + iM3->charge() == 0) { mass_1_3 = (lepton1 + lepton3).mass(); pt_1_3 = (lepton1 + lepton3).pt(); eta_1_3 = (lepton1 + lepton3).eta(); eta_1_3 = (lepton1 + lepton3).eta(); phi_1_3 = (lepton1 + lepton3).phi();}
                if (iM1->charge() + iM4->charge() == 0) { mass_1_4 = (lepton1 + lepton4).mass(); pt_1_4 = (lepton1 + lepton4).pt(); eta_1_4 = (lepton1 + lepton4).eta(); eta_1_4 = (lepton1 + lepton4).eta(); phi_1_4 = (lepton1 + lepton4).phi();}
                if (iM2->charge() + iM3->charge() == 0) { mass_2_3 = (lepton2 + lepton3).mass(); pt_2_3 = (lepton2 + lepton3).pt(); eta_2_3 = (lepton2 + lepton3).eta(); eta_2_3 = (lepton2 + lepton3).eta(); phi_2_3 = (lepton2 + lepton3).phi();}
                if (iM2->charge() + iM4->charge() == 0) { mass_2_4 = (lepton2 + lepton4).mass(); pt_2_4 = (lepton2 + lepton4).pt(); eta_2_4 = (lepton2 + lepton4).eta(); eta_2_4 = (lepton2 + lepton4).eta(); phi_2_4 = (lepton2 + lepton4).phi();}
                if (iM3->charge() + iM4->charge() == 0) { mass_3_4 = (lepton3 + lepton4).mass(); pt_3_4 = (lepton3 + lepton4).pt(); eta_3_4 = (lepton3 + lepton4).eta(); eta_3_4 = (lepton3 + lepton4).eta(); phi_3_4 = (lepton3 + lepton4).phi();}
                
                if (iM1->charge() + iM2->charge() + iM3->charge() + iM4->charge() == 0){ //this check should be redundant, because it's the first thing we check when we enter the quad loop
                  mass_1_2_3_4 = (lepton1 + lepton2 + lepton3 + lepton4).mass();
                  pt_1_2_3_4   = (lepton1 + lepton2 + lepton3 + lepton4).pt();
                  eta_1_2_3_4  = (lepton1 + lepton2 + lepton3 + lepton4).eta();
                  phi_1_2_3_4  = (lepton1 + lepton2 + lepton3 + lepton4).phi();
                
                }
                 
                bool isZplusY = false; 
                // ********************

                bool match_12_34_56 = false;

                if (mass_1_2 > upsi_mass_low && mass_1_2 < upsi_mass_high) 
                  if (mass_3_4 > Z_mass_low && mass_3_4 < Z_mass_high) 
                      match_12_34_56 = true;

                if (mass_1_2 > Z_mass_low && mass_1_2 < Z_mass_high) 
                  if (mass_3_4 > upsi_mass_low && mass_3_4 < upsi_mass_high) 
                      match_12_34_56 = true;

                // ********************

                bool match_13_24_56 = false;

                if (mass_1_3 > upsi_mass_low && mass_1_3 < upsi_mass_high) 
                  if (mass_2_4 > Z_mass_low && mass_2_4 < Z_mass_high) 
                      match_13_24_56 = true;

                if (mass_1_3 > Z_mass_low && mass_1_3 < Z_mass_high) 
                  if (mass_2_4 > upsi_mass_low && mass_2_4 < upsi_mass_high) 
                      match_13_24_56 = true;

                // ********************

                bool match_14_23_56 = false;

                if (mass_1_4 > upsi_mass_low && mass_1_4 < upsi_mass_high) 
                  if (mass_2_3 > Z_mass_low && mass_2_3 < Z_mass_high) 
                      match_14_23_56 = true;

                if (mass_1_4 > Z_mass_low && mass_1_4 < Z_mass_high) 
                  if (mass_2_3 > upsi_mass_low && mass_2_3 < upsi_mass_high) 
                      match_14_23_56 = true;
                      
            
            

                      
             //Remove cut on matching so that we allow quads into the sample that might work for Z->4 mu.

            //     if (!(match_12_34_56 || match_13_24_56 || 
//                       match_14_23_56)) {
//                    histContainer_["CutFlow"]->Fill(15);  
//                    std::cout << "FAILED At CUT 14" <<std::endl;
//                   continue;
//                 }
                   
                   if (match_12_34_56 || match_13_24_56 || match_14_23_56){
                     isZplusY = true;
                   }
                
                
                //******************************   
                bool isZOnly = false;
                //*******************
                
                bool match_12_34_ZOnly = false;
                //****************************
                
                if (mass_1_2_3_4 > ZOnly_mass_low && mass_1_2_3_4 < ZOnly_mass_high){
                  if (mass_1_2 > lowerBound_ZOnly && mass_3_4 > lowerBound_ZOnly){
                    match_12_34_ZOnly = true;
                   // std::cout << "match_12_34_ZOnly:  " << match_12_34_ZOnly << std::endl; 
                   }
                }
                 
             //      if (mass_1_2_3_4 > ZOnly_mass_low && mass_1_2_3_4 < ZOnly_mass_high){
//                     if (mass_3_4 > avoid_JPsi_ZOnly && mass_1_2 > avoid_upsi_ZOnly){
//                       match_12_34_ZOnly = true;
//                     //  std::cout << "match_12_34_ZOnly:  " << match_12_34_ZOnly << std::endl;
//                     }
//                   }
                 
                 bool match_13_24_ZOnly = false;
                 //****************************
                 if (mass_1_2_3_4 > ZOnly_mass_low && mass_1_2_3_4 < ZOnly_mass_high){
                   if (mass_1_3 > lowerBound_ZOnly && mass_2_4 > lowerBound_ZOnly){
                     match_13_24_ZOnly = true;
                   //  std::cout << "match_13_24_ZOnly:  " << match_13_24_ZOnly << std::endl;  
                   }
                 
                 }
                 
       //           if (mass_1_2_3_4 > ZOnly_mass_low && mass_1_2_3_4 < ZOnly_mass_high){
//                    if (mass_2_4 > avoid_JPsi_ZOnly && mass_1_3 > avoid_upsi_ZOnly){
//                      match_13_24_ZOnly = true; 
//                    //  std::cout << "match_13_24_ZOnly:  " << match_13_24_ZOnly << std::endl;
//                    }
//                  
//                  }
                 
                 bool match_14_23_ZOnly = false;
                 //****************************
                 if (mass_1_2_3_4 > ZOnly_mass_low && mass_1_2_3_4 < ZOnly_mass_high){
                   if (mass_1_4 > lowerBound_ZOnly && mass_2_3 > lowerBound_ZOnly){
                     match_14_23_ZOnly = true;
                     //std::cout << "match_14_23_ZOnly:  " << match_14_23_ZOnly << std::endl; 
                   }
                 }
                 
               //   if (mass_1_2_3_4 > ZOnly_mass_low && mass_1_2_3_4 < ZOnly_mass_high){
//                    if (mass_2_3 > avoid_JPsi_ZOnly && mass_1_4 > avoid_upsi_ZOnly){
//                      match_14_23_ZOnly = true;
//                    //  std::cout << "match_14_23_ZOnly:  " << match_14_23_ZOnly << std::endl;
//                    }
//                  }
                 
                 if (match_12_34_ZOnly || match_13_24_ZOnly || match_14_23_ZOnly){
                   isZOnly = true; 
                 }
                 
               //   //The statement below is normally what you want commented in 
                 if (! (isZplusY || isZOnly)){
                   histContainer_["CutFlow"]->Fill(16);
                   continue;
                   
                 }
                //***** TRIGGER EFFICIENCY STUDIES *****************

//              Comment this in if doing ZplusY trig eff studies
            //     if (!(isZplusY)){
//                   continue;
//                 }
                //Comment this in if doing ZOnly trig eff studies 
               //  if (!(isZOnly)){
//                   continue; 
//                 }
                 
                 

                  //At this point, we might have more than 1 Y+Z candidate in each event (e.g. from matching ambiguity within one set of four muons or from making 2 or more quads out of >4 muons). We will sort this out in the phase 2 code, not here.
                save_event = true; //if we found a good set, turn save_event to true 
                
                std::vector<double> dimuon1mass_temp, dimuon2mass_temp;
                std::vector<double> dimuon1pt_temp, dimuon2pt_temp;
                std::vector<double> dimuon1eta_temp, dimuon2eta_temp;
                std::vector<double> dimuon1phi_temp, dimuon2phi_temp;
                
                dimuon1mass_temp.clear(); dimuon2mass_temp.clear();
                dimuon1pt_temp.clear(); dimuon2pt_temp.clear();
                dimuon1eta_temp.clear(); dimuon2eta_temp.clear();
                dimuon1phi_temp.clear(); dimuon2phi_temp.clear();
                
                // match_12_34_56 case, aka the ZplusY case
                if (match_12_34_56) {dimuon1mass_temp.push_back(mass_1_2); dimuon2mass_temp.push_back(mass_3_4);  dimuon1pt_temp.push_back(pt_1_2); dimuon2pt_temp.push_back(pt_3_4);  dimuon1eta_temp.push_back(eta_1_2); dimuon2eta_temp.push_back(eta_3_4);  dimuon1phi_temp.push_back(phi_1_2); dimuon2phi_temp.push_back(phi_3_4); }
                else {
                  dimuon1mass_temp.push_back(-99);
                  dimuon2mass_temp.push_back(-99);
                  dimuon1pt_temp.push_back(-99);
                  dimuon2pt_temp.push_back(-99);
                  dimuon1eta_temp.push_back(-99);
                  dimuon2eta_temp.push_back(-99);
                  dimuon1phi_temp.push_back(-99);
                  dimuon2phi_temp.push_back(-99);
                }                                                                                                                                                                                                                                                                                                                                                                                                                                
                
                //match_13_24_56 case, aka the ZplusY case
                if (match_13_24_56) {dimuon1mass_temp.push_back(mass_1_3); dimuon2mass_temp.push_back(mass_2_4); dimuon1pt_temp.push_back(pt_1_3); dimuon2pt_temp.push_back(pt_2_4);dimuon1eta_temp.push_back(eta_1_3); dimuon2eta_temp.push_back(eta_2_4); dimuon1phi_temp.push_back(phi_1_3); dimuon2phi_temp.push_back(phi_2_4); }
                else {
                  dimuon1mass_temp.push_back(-99);
                  dimuon2mass_temp.push_back(-99);
                  dimuon1pt_temp.push_back(-99);
                  dimuon2pt_temp.push_back(-99);
                  dimuon1eta_temp.push_back(-99);
                  dimuon2eta_temp.push_back(-99);
                  dimuon1phi_temp.push_back(-99);
                  dimuon2phi_temp.push_back(-99);
                  
                }   
                 
                 //match_14_23_56 case, aka the ZplusY case                                                                                                                                                                                                                                                                                                  
                if (match_14_23_56) {dimuon1mass_temp.push_back(mass_1_4); dimuon2mass_temp.push_back(mass_2_3); dimuon1pt_temp.push_back(pt_1_4); dimuon2pt_temp.push_back(pt_2_3);dimuon1eta_temp.push_back(eta_1_4); dimuon2eta_temp.push_back(eta_2_3); dimuon1phi_temp.push_back(phi_1_4); dimuon2phi_temp.push_back(phi_2_3); }
                else {
                  dimuon1mass_temp.push_back(-99);
                  dimuon2mass_temp.push_back(-99);
                  dimuon1pt_temp.push_back(-99);
                  dimuon2pt_temp.push_back(-99);
                  dimuon1eta_temp.push_back(-99);
                  dimuon2eta_temp.push_back(-99);
                  dimuon1phi_temp.push_back(-99);
                  dimuon2phi_temp.push_back(-99);
                }  
                
                // match_12_34_ZOnly case
                if (match_12_34_ZOnly) {dimuon1mass_temp.push_back(mass_1_2); dimuon2mass_temp.push_back(mass_3_4);  dimuon1pt_temp.push_back(pt_1_2); dimuon2pt_temp.push_back(pt_3_4);  dimuon1eta_temp.push_back(eta_1_2); dimuon2eta_temp.push_back(eta_3_4);  dimuon1phi_temp.push_back(phi_1_2); dimuon2phi_temp.push_back(phi_3_4); }
                else {
                  dimuon1mass_temp.push_back(-99);
                  dimuon2mass_temp.push_back(-99);
                  dimuon1pt_temp.push_back(-99);
                  dimuon2pt_temp.push_back(-99);
                  dimuon1eta_temp.push_back(-99);
                  dimuon2eta_temp.push_back(-99);
                  dimuon1phi_temp.push_back(-99);
                  dimuon2phi_temp.push_back(-99);
                }                                                                                                                                                                                                                                                                                                                                                                                                                                
                
                // match_13_24_ZOnly case
                 if (match_13_24_ZOnly) {dimuon1mass_temp.push_back(mass_1_3); dimuon2mass_temp.push_back(mass_2_4); dimuon1pt_temp.push_back(pt_1_3); dimuon2pt_temp.push_back(pt_2_4);dimuon1eta_temp.push_back(eta_1_3); dimuon2eta_temp.push_back(eta_2_4); dimuon1phi_temp.push_back(phi_1_3); dimuon2phi_temp.push_back(phi_2_4); }
                 else {
                  dimuon1mass_temp.push_back(-99);
                  dimuon2mass_temp.push_back(-99);
                  dimuon1pt_temp.push_back(-99);
                  dimuon2pt_temp.push_back(-99);
                  dimuon1eta_temp.push_back(-99);
                  dimuon2eta_temp.push_back(-99);
                  dimuon1phi_temp.push_back(-99);
                  dimuon2phi_temp.push_back(-99);
                  
                } 
                
                //match_14_23_ZOnly case
                 if (match_14_23_ZOnly) {dimuon1mass_temp.push_back(mass_1_4); dimuon2mass_temp.push_back(mass_2_3); dimuon1pt_temp.push_back(pt_1_4); dimuon2pt_temp.push_back(pt_2_3);dimuon1eta_temp.push_back(eta_1_4); dimuon2eta_temp.push_back(eta_2_3); dimuon1phi_temp.push_back(phi_1_4); dimuon2phi_temp.push_back(phi_2_3); }
                 else {
                  dimuon1mass_temp.push_back(-99);
                  dimuon2mass_temp.push_back(-99);
                  dimuon1pt_temp.push_back(-99);
                  dimuon2pt_temp.push_back(-99);
                  dimuon1eta_temp.push_back(-99);
                  dimuon2eta_temp.push_back(-99);
                  dimuon1phi_temp.push_back(-99);
                  dimuon2phi_temp.push_back(-99);
                }  
                
                dimuon1mass.push_back(dimuon1mass_temp);
                dimuon2mass.push_back(dimuon2mass_temp); 
                dimuon1pt.push_back(dimuon1pt_temp);
                dimuon2pt.push_back(dimuon2pt_temp);
                dimuon1eta.push_back(dimuon1eta_temp);
                dimuon2eta.push_back(dimuon2eta_temp);
                dimuon1phi.push_back(dimuon1phi_temp);
                dimuon2phi.push_back(dimuon2phi_temp); 
                  
                
//**** Fresh attempt at trigger matching, using the triggered method of miniaod, see: https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_11_3_1/doc/html/d6/d13/classpat_1_1Muon.html, thanks to S.L. for finding this
              //   bool iM1Triggered = false;
//                 bool iM2Triggered = false;
//                 bool iM3Triggered = false;
//                 bool iM4Triggered = false;
//                 //note to self: reusing the same variables this way so that I take advantage of the fact that they have local scope and can just be recycled is not very good practice and I should avoid it from now on
//                 //check iM1
//                 int iM1TrigMatched = 0;
//                 for (unsigned int i = 0, n = triggerlist.size(); i < n; i++){
//                  // std::cout << "triggerlist.at(i)  " << triggerlist.at(i) << std::endl;
//                   std::string tempStr (triggerlist.at(i)); //need to convert string to char as S.L. showed me in email dated 19 Oct. 2021, subject = error when trying to use triggered in a loop over the triggerlist
//                   const char* c = tempStr.c_str(); 
//                   iM1Triggered = iM1->triggered(c);
//                  // std::cout << "tempStr  " << tempStr << std::endl; 
//                  // std::cout << "iM1Triggered is  " << iM1Triggered << std::endl; 
//                   if (iM1Triggered){
//                     iM1TrigMatched += 1;
//                     break;
//                   }
//                 }
//                   
//                 //check iM2
//                 int iM2TrigMatched = 0;
//                 for (unsigned int i = 0, n = triggerlist.size(); i < n; i++){
//                   std::string tempStr (triggerlist.at(i));
//                   const char* c = tempStr.c_str(); 
//                   iM2Triggered = iM2->triggered(c);
//                 //  std::cout << "tempStr  " << tempStr << std::endl;
//                //   std::cout << "iM2Triggered is  " << iM2Triggered << std::endl; 
//                   if (iM2Triggered){
//                     iM2TrigMatched +=1;
//                     break;
//                   }
//                 }
//                
//                //check iM3
//                 int iM3TrigMatched = 0;
//                 for (unsigned int i=0, n = triggerlist.size(); i < n; i++){
//                   std::string tempStr (triggerlist.at(i));
//                   const char* c = tempStr.c_str();
//                   iM3Triggered = iM3->triggered(c);
//                  // std::cout << "tempStr  " << tempStr << std::endl;
//                 //  std::cout << "iM3Triggered is  " << iM3Triggered << std::endl; 
//                   if (iM3Triggered){
//                     iM3TrigMatched += 1;
//                     break;
//                   }
//                 }
//                 
//                 //check iM4
//                 
//                 int iM4TrigMatched = 0;
//                 for (unsigned int i=0, n = triggerlist.size(); i < n; i++){
//                   std::string tempStr (triggerlist.at(i));
//                   const char* c = tempStr.c_str();
//                   iM4Triggered = iM4->triggered(c);
//                  // std::cout << "tempStr  " << tempStr << std::endl;
//                  // std::cout << "iM4Triggered is  " << iM4Triggered << std::endl; 
//                   if (iM4Triggered){
//                     iM4TrigMatched += 1;
//                     break;
//                   }
//                 }
//                 
//                 int totalTrigMatched = iM1TrigMatched + iM2TrigMatched + iM3TrigMatched + iM4TrigMatched;
//                 std::cout << "totalTrigMatched  " << totalTrigMatched << std::endl;                
//                
               
                
            //
                    

                
                 // iM1Triggered = iM1->triggered("HLT_diMu12_213_v4");
                 // std::cout << "iM1Triggered " << iM1Triggered << std::endl;
               
               
                 
                
                
 






//// ****************
//// trigger matching
//// ****************
//            pat::TriggerObjectStandAlone TO;
////            for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
////              std::string name = "HLT_IsoMu18_v5"; //names.triggerName(i);
//            std::string name2 = "HLT_Dimuon0_Jpsi3p5_Muon2_v2";
//            std::string name4 = "HLT_Dimuon0_Jpsi3p5_Muon2_v4";
//            std::string name5 = "HLT_Dimuon0_Jpsi3p5_Muon2_v5";
//            std::string name6 = "HLT_Dimuon0_Jpsi3p5_Muon2_v6";
//            int matched = 0;
//            for (unsigned int iTO = 0, n = triggerObjects->size(); iTO < n; ++iTO) {
//              TO = triggerObjects->at(iTO);
//              TO.unpackPathNames( names );
//              std::cout << "found trigger\n";
//              if (TO.hasPathName(name4)) {
////                double tempDR1 = TMath::Sqrt(TMath::Power(TO.eta() - iM1->eta(),2) + TMath::Power(TO.phi() - iM1->phi(),2));
////                double tempDR2 = TMath::Sqrt(TMath::Power(TO.eta() - iM2->eta(),2) + TMath::Power(TO.phi() - iM2->phi(),2));
////                double tempDR3 = TMath::Sqrt(TMath::Power(TO.eta() - iM3->eta(),2) + TMath::Power(TO.phi() - iM3->phi(),2));
////                double tempDR4 = TMath::Sqrt(TMath::Power(TO.eta() - iM4->eta(),2) + TMath::Power(TO.phi() - iM4->phi(),2));
////                std::cout << iTO << " " << tempDR1 << " " << tempDR2 << " " << tempDR3 << " " << tempDR4 << std::endl;
//                std::cout << TO.pt() << " " << iM1->pt() << " " << iM2->pt() <<  " " << iM3->pt() <<  " " << iM4->pt() << std::endl;
////                if (tempDR1<0.1 || tempDR2<0.1 || tempDR3<0.1 || tempDR4<0.1)
////                  matched++;
//              }
//            }
//            if (matched<2)
//              continue;
////            std::cout << "end of muon match thing: " << matched << std::endl;
////            }

  
                //****************
                //trigger matching //NEW PART HERE
                //****************
                
                //Credits to:  https://riptutorial.com/cplusplus/example/11151/find-max-and-min-element-and-respective-index-in-a-vector
                // https://stackoverflow.com/questions/571394/how-to-find-out-if-an-item-is-present-in-a-stdvector
                //https://www.techiedelight.com/check-vector-contains-given-element-cpp/
                 
                 double trigMatchDRCut = 0.01; //I set this based on looking at the output of the trigDR histo 0.1 is definitely too big, might still tune this 0.01
//               //  pat::TriggerObjectStandAlone TO;
//                 quadHasHowManyTrigMatches.clear(); //I think is where the bug was, I think this line always needs to be commented out and it should be cleared along with all the other stuff (around line 496 in this case, search "sanity check" to find it)
                 int trigMatched = 0;
                 std::vector<int> minElementIndexVec;
                 minElementIndexVec.clear();
                 std::cout << "Looking at trigger to reco matching module" << std::endl;
                 for (unsigned int iTO = 0; iTO < triggerObjects->size(); iTO++){
//                 //  std::cout << "triggerObjects->size()  " << triggerObjects->size() << std::endl;
                   pat::TriggerObjectStandAlone TO;
                   TO = triggerObjects->at(iTO);
                   TO.unpackPathNames(names);
    // }     //extra bracket that I put in now to keep all the nested loops happy, I may kill it when I eventually find out where the original enclosing bracket was       
//                  //Closely adapted from https://github.com/DryRun/BParkingNANO/blob/master/BParkingNano/plugins/MuonTriggerSelector.cc#L115-L123
                  TO.unpackFilterLabels(iEvent, *triggerBits);
                 
                 //  if(TMath::Abs(TO.pdgId()) != 13) continue; //make sure you are looking at TO muons!
//                   bool isTriggerMuon = false;
//                   for (unsigned h = 0; h < TO.filterIds().size(); ++h){
// 	                if(TO.filterIds()[h] == 83){ //83 = muon //see: https://github.com/cms-sw/cmssw/blob/master/DataFormats/HLTReco/interface/TriggerTypeDefs.h
// 	                  isTriggerMuon = true; 
// 	                  break;
// 	                }
// 	              }
 //Thanks to Sam Harper for his help with this section. See email and mattermost for discussions with him that greatly influenced the structure of this portion. Used his code, the printPathFilters.py, to get the filterNames. Remember to run the printPathFilters.py as "python3 printPathFilters.py <nameOfFileYouAreExaming>" (You need to expressly say python3)
                  bool isTriggerMuon = false;
                  for (unsigned h = 0; h < TO.filterLabels().size(); h++){
                    std::string filterName = TO.filterLabels()[h];
                    
                    if (triggerYear == 2016) {
                      if (singleMu2016Trig1Fired){
                        if (filterName.find("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09") != std::string::npos ){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for singleMu2016Trig1Fired" << std::endl; 
                          break;
                        }
                      }
                      if (singleMu2016Trig2Fired){
                        if (filterName.find("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09") != std::string::npos ){
                          isTriggerMuon = true; 
                          std::cout << "Found a filterName match for singleMu2016Trig2Fired" << std::endl;
                          break;
                        }
                      }
                      if (doubleMu2016Trig1Fired){
                        if (filterName.find("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2") != std::string::npos ){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for doubleMu2016Trig1Fired" << std::endl; 
                          break; 
                        }
                      }
                      if (doubleMu2016Trig2Fired){
                        if (filterName.find("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2") != std::string::npos ){
                          isTriggerMuon = true; 
                          std::cout << "Found a filterName match for doubleMu2016Trig2Fired" << std::endl; 
                          break; 
                        }
                      }
                     
                     if (tripleMu2016Trig1Fired){
                       if (filterName.find("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105") != std::string::npos){
                         isTriggerMuon = true;
                         std::cout << "Found a filterName match for tripleMu2016Trig1Fired" << std::endl;
                         break;
                       }
                     } 
                    
                    }
                    if (triggerYear == 2017){
                      if (singleMu2017Trig1Fired){
                        if (filterName.find("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07") != std::string::npos){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for singleMu2017Trig1Fired" << std::endl;
                          break;
                        }
                      }
                      if (doubleMu2017Trig1Fired){
                        if (filterName.find("hltDiMuon178Mass8Filtered") != std::string::npos){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for doubleMu2017Trig1Fired" << std::endl;
                          break; 
                        }
                      }
                      if (tripleMu2017Trig1Fired){
                        if (filterName.find("hltL3fL1TripleMu553f0Filtered12105") != std::string::npos){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for tripleMu2017Trig1Fired" << std::endl;
                          break;
                        }
                      }
                      if (tripleMu2017Trig2Fired){
                        if (filterName.find("hltTripleMu555TripleDZ0p2") != std::string::npos){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for tripleMu2017Trig2Fired" << std::endl;
                          break; 
                        }
                      }
                      
                    }
                    if (triggerYear == 2018){
                      if (singleMu2018Trig1Fired){
                        if (filterName.find("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07") != std::string::npos){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for singleMu2018Trig1Fired" << std::endl;
                          break;
                        }
                      }
                      if (doubleMu2018Trig1Fired){
                        if (filterName.find("hltDiMuon178Mass3p8Filtered") != std::string::npos){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for doubleMu2018Trig1Fired" << std::endl; 
                          break;
                        }
                      }
                      if (doubleMu2018Trig2Fired){
                        if (filterName.find("hltDiMuon178Mass8Filtered") != std::string::npos){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for doubleMu2018Trig2Fired" << std::endl; 
                          break;
                        }
                      } 
                      
                      if (tripleMu2018Trig1Fired){
                        if (filterName.find("hltL3fL1TripleMu553f0Filtered12105") != std::string::npos){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for tripleMu2018Trig1Fired" << std::endl;
                          break;
                        }
                      }
                      if (tripleMu2018Trig2Fired){
                        if (filterName.find("hltTripleMu555TripleDZ0p2") != std::string::npos){
                          isTriggerMuon = true;
                          std::cout << "Found a filterName match for tripleMu2018Trig2Fired" << std::endl; 
                          break; 
                        }
                      }
                      
                    }
                  
                  }
// 	                
 	              if(!isTriggerMuon) continue;  
 	               
//                  
//                   
                  std::vector<double> tempDRVec;
                  tempDRVec.clear();
                  
                  double PI = 3.14159265358979323846; //https://root.cern.ch/root/html524/TMath.html#TMath:Pi
                  double TWO_PI = 6.28318530717958647692; //two times the pi as defined above
                  
                  double tempDEta1Sq = TMath::Power(TO.eta() - iM1->eta(),2);
                  double tempDEta2Sq = TMath::Power(TO.eta() - iM2->eta(),2);
                  double tempDEta3Sq = TMath::Power(TO.eta() - iM3->eta(),2);
                  double tempDEta4Sq = TMath::Power(TO.eta() - iM4->eta(),2);
                  
                  //Thanks to David Yu for reminding me to put this protection in place
                  double tempDPhi1 = fabs(TO.phi() - iM1->phi());
                  while( tempDPhi1 > PI ){
                    tempDPhi1 -= TWO_PI;
                  }
                  
                  double tempDPhi2 = fabs(TO.phi() - iM2->phi());
                  while( tempDPhi2 > PI ){
                    tempDPhi2 -= TWO_PI;
                  }
                  
                  double tempDPhi3 = fabs(TO.phi() - iM3->phi());
                  while( tempDPhi3 > PI ){
                    tempDPhi3 -= TWO_PI;
                  }
                  
                  double tempDPhi4 = fabs(TO.phi() - iM4->phi());
                  while( tempDPhi4 > PI ){
                    tempDPhi4 -= TWO_PI;
                  }
                  
                  
                  double tempDPhi1Sq = TMath::Power(tempDPhi1,2);
                  double tempDPhi2Sq = TMath::Power(tempDPhi2,2);
                  double tempDPhi3Sq = TMath::Power(tempDPhi3,2);
                  double tempDPhi4Sq = TMath::Power(tempDPhi4,2);
                  
                  double tempDR1 = TMath::Sqrt(tempDEta1Sq + tempDPhi1Sq);
                  double tempDR2 = TMath::Sqrt(tempDEta2Sq + tempDPhi2Sq);
                  double tempDR3 = TMath::Sqrt(tempDEta3Sq + tempDPhi3Sq);
                  double tempDR4 = TMath::Sqrt(tempDEta4Sq + tempDPhi4Sq);
//                   
                  
                   
//            //       double tempDR1 = TMath::Sqrt(TMath::Power(TO.eta() - iM1->eta(),2) + TMath::Power(TO.phi() - iM1->phi(),2));
//             //      double tempDR2 = TMath::Sqrt(TMath::Power(TO.eta() - iM2->eta(),2) + TMath::Power(TO.phi() - iM2->phi(),2));
//               //    double tempDR3 = TMath::Sqrt(TMath::Power(TO.eta() - iM3->eta(),2) + TMath::Power(TO.phi() - iM3->phi(),2));
//                 //  double tempDR4 = TMath::Sqrt(TMath::Power(TO.eta() - iM4->eta(),2) + TMath::Power(TO.phi() - iM4->phi(),2));
//                  // std::cout << "Looking at trigger to reco matching module" << std::endl;
//                   
//                   
                   std::cout << iTO << " " << tempDR1 << " " << tempDR2 << " " << tempDR3 << " " << tempDR4 << std::endl;
                   std::cout << TO.pt() << " " << iM1->pt() << " " << iM2->pt() <<  " " << iM3->pt() <<  " " << iM4->pt() << std::endl;
    
                   tempDRVec.push_back(tempDR1); tempDRVec.push_back(tempDR2); tempDRVec.push_back(tempDR3); tempDRVec.push_back(tempDR4);
//             //      std::cout << "tempDRVec.at(0)  " << tempDRVec.at(0) << std::endl; 
                   int minElementIndex = std::min_element(tempDRVec.begin(),tempDRVec.end()) - tempDRVec.begin();
                   std::cout << "minElementIndex  " << minElementIndex << std::endl;
//                   
                  std::cout << "Collection  " << TO.collection() << std::endl; //this shows us that the duplicates are likely the same muon that has passed
//               //different triggers. The ones that print out as being the same are listed in each case as being part of a different collection each time
//               //To handle the fact that there are duplicates, when we do the matching, we use the more intensive matching I've done below as opposed to the simpler
//               // matching of  if (tempDR1<0.1 || tempDR2<0.1 || tempDR3<0.1 || tempDR4<0.1)
//                  //               matched++;
//                  //that was initially suggested by S.L.
//               
        
                  double minElement = *std::min_element(tempDRVec.begin(), tempDRVec.end());
                  std::cout << "minElement  " << minElement << std::endl;
                  histContainer_["trigDR"]->Fill(minElement); //trigDR is a histo of the trig Mu-Reco Mu pair for a given quad that has the smallest separation in deltaR 
                  
//                   
                  if (trigMatched == 0){
                    if (minElement < trigMatchDRCut){
                      trigMatched++;
                      minElementIndexVec.push_back(minElementIndex);
                      std::cout << "minElementIndexVec.at(0)  " << minElementIndexVec.at(0) << std::endl;
                      std::cout << "trigMatched is   " << trigMatched << std::endl;
                    }
                  
                  }
        //        } //I just commented this out 
//                   
//                 //  std::vector<int> testVec;
// //                  testVec.push_back(0);
// //                  std::cout << "testVec.at(0)  " << testVec.at(0) << std::endl; 
// //                  if ( std::count(testVec.begin(), testVec.end(), 1) ) {
// //                     std::cout << "true, 1 is found" << std::endl;
// //                  }
// //                  
// //                  if ( std::count(testVec.begin(), testVec.end(), 1) == 0 ) {
// //                     std::cout << "false, 1 is NOT found" << std::endl;
// //                  }
//                  
//            //      std::cout << "minElementIndex is:  " <<  minElementIndex << std::endl;
//            //      for (unsigned int i =0; i < minElementIndexVec.size(); i++){
//              //      std::cout << "minElementIndexVec.at(i)  " << minElementIndexVec.at(i) << std::endl; 
//              //    }
//                  
                 if (trigMatched > 0){
               //    std::cout << "Entered test area" << std::endl;
                   if ( std::count(minElementIndexVec.begin(), minElementIndexVec.end(), minElementIndex) == 0 ){
                 //    std::cout << "Entered here line 1077" << std::endl; 
                     if (minElement < trigMatchDRCut){
                       trigMatched++;
                       std::cout << "trigMatched > 1  " << trigMatched << std::endl; 
                       minElementIndexVec.push_back(minElementIndex);
                     }
                   
                   }
                 }
             }     
//                    
//                   
//          //          if (trigMatched > 0){
// //                      std::cout << "Entered test area" << std::endl;
// //                    // if ( (std::find(minElementIndexVec.begin(), minElementIndexVec.end(), minElementIndex) != minElementIndexVec.end() ) == false
// //                     // && minElement < trigMatchDRCut) {
// //                      if (minElement < trigMatchDRCut){
// //                        if (std::find(minElementIndexVec.begin(), minElementIndexVec.end(), minElementIndex) != minElementIndexVec.end() == false){
// //                         std::cout << "Got here line 1052" << std::endl; 
// //                         trigMatched++;
// //                         minElementIndexVec.push_back(minElementIndex);
// //                      }
// //                     }
// //                   
// //                   }
//                   
//             
//                // quadHasHowManyTrigMatches.push_back(trigMatched);
//                 }
//                 
//                 quadHasHowManyTrigMatches.push_back(trigMatched); //test comment //test comment 2

// *********     
// VERTEXING 
// *********
                const reco::TrackRef track1 = iM1->muonBestTrack();
                const reco::TrackRef track2 = iM2->muonBestTrack();
                const reco::TrackRef track3 = iM3->muonBestTrack();
                const reco::TrackRef track4 = iM4->muonBestTrack();
                
                //build all the track pairs
                std::vector<reco::TransientTrack> transient_tracks_muon12; transient_tracks_muon12.clear(); 
                std::vector<reco::TransientTrack> transient_tracks_muon13; transient_tracks_muon13.clear(); 
                std::vector<reco::TransientTrack> transient_tracks_muon14; transient_tracks_muon14.clear(); 
                std::vector<reco::TransientTrack> transient_tracks_muon23; transient_tracks_muon23.clear(); 
                std::vector<reco::TransientTrack> transient_tracks_muon24; transient_tracks_muon24.clear(); 
                std::vector<reco::TransientTrack> transient_tracks_muon34; transient_tracks_muon34.clear(); 
                
        
                transient_tracks_muon12.push_back((*builder).build(track1.get())); transient_tracks_muon12.push_back((*builder).build(track2.get()));
                transient_tracks_muon13.push_back((*builder).build(track1.get())); transient_tracks_muon13.push_back((*builder).build(track3.get()));
                transient_tracks_muon14.push_back((*builder).build(track1.get())); transient_tracks_muon14.push_back((*builder).build(track4.get()));
                transient_tracks_muon23.push_back((*builder).build(track2.get())); transient_tracks_muon23.push_back((*builder).build(track3.get()));
                transient_tracks_muon24.push_back((*builder).build(track2.get())); transient_tracks_muon24.push_back((*builder).build(track4.get()));
                transient_tracks_muon34.push_back((*builder).build(track3.get())); transient_tracks_muon34.push_back((*builder).build(track4.get()));
                
                //build big 4 mu Vertex
                
                std::vector<reco::TransientTrack> transient_tracks_muon1234; transient_tracks_muon1234.clear();
                transient_tracks_muon1234.push_back((*builder).build(track1.get())); transient_tracks_muon1234.push_back((*builder).build(track2.get())); 
                transient_tracks_muon1234.push_back((*builder).build(track3.get())); transient_tracks_muon1234.push_back((*builder).build(track4.get()));
                
                 //build the tracks for each muon // the stuff for single muons is currently non-compiling garbage and S.L. and I aren't even sure what we would want to build the vertex wrt, so commenting out and if we ever figure out what we meant, then we can fix and uncomment in 
             //    std::vector<reco::TransientTrack> transient_tracks_muon1; transient_tracks_muon1.clear();
//                 std::vector<reco::TransientTrack> transient_tracks_muon2; transient_tracks_muon2.clear();
//                 std::vector<reco::TransientTrack> transient_tracks_muon3; transient_tracks_muon3.clear();
//                 std::vector<reco::TransientTrack> transient_tracks_muon4; transient_tracks_muon4.clear();
//                 
//                 transient_tracks_muon1.push_back((*builder).build(track1.get())); 
//                 transient_tracks_muon2.push_back((*builder).build(track2.get()));
//                 transient_tracks_muon3.push_back((*builder).build(track3.get()));
//                 transient_tracks_muon4.push_back((*builder).build(track4.get()));
//                 
                

                std::cout << "DEBUG 1" << std::endl;
                KalmanVertexFitter kalman_fitter;
                //diumon vertices
                TransientVertex vtx12; vtx12 = kalman_fitter.vertex(transient_tracks_muon12);
                TransientVertex vtx13; vtx13 = kalman_fitter.vertex(transient_tracks_muon13);
                TransientVertex vtx14; vtx14 = kalman_fitter.vertex(transient_tracks_muon14);
                TransientVertex vtx23; vtx23 = kalman_fitter.vertex(transient_tracks_muon23);
                TransientVertex vtx24; vtx24 = kalman_fitter.vertex(transient_tracks_muon24);
                TransientVertex vtx34; vtx34 = kalman_fitter.vertex(transient_tracks_muon34);
                std::cout << "DEBUG 2" << std::endl;
                
                //big 4 mu vtx
                TransientVertex vtx1234; vtx1234 = kalman_fitter.vertex(transient_tracks_muon1234);
                std::cout << "DEBUG 3" << std::endl;
                
                if (vtx1234.isValid()) {
                   big4MuVtx.push_back(TMath::Prob(vtx1234.totalChiSquared(), int(vtx1234.degreesOfFreedom())));
                
                }
                
                else
                   big4MuVtx.push_back(-1000.);
                   
                   //more of the single mu stuff that is currently meaningless, maybe someday it will find its raison d'etre 
        //        TransientVertex vtx1; vtx1 = kalman_fitter.vertex(transient_tracks_muon1);
            //    TransientVertex vtx2; vtx2 = kalman_fitter.vertex(transient_tracks_muon2);
             //   TransientVertex vtx3; vtx3 = kalman_fitter.vertex(transient_tracks_muon3);
              //  TransientVertex vtx4; vtx4 = kalman_fitter.vertex(transient_tracks_muon4);
                
                //get vertex info for the pair vertices (assuming they are valid, else push back -1000. If no match, push back -10000)
               
                std::vector<double> dimuon1vtx_temp, dimuon2vtx_temp;
                dimuon1vtx_temp.clear(); dimuon2vtx_temp.clear();
            
                std::vector<double> dimuon1vtx_xpos_temp, dimuon2vtx_xpos_temp;
                dimuon1vtx_xpos_temp.clear(); dimuon2vtx_xpos_temp.clear();
                
                std::vector<double> dimuon1vtx_ypos_temp, dimuon2vtx_ypos_temp;
                dimuon1vtx_ypos_temp.clear(); dimuon2vtx_ypos_temp.clear();
                
                std::vector<double> dimuon1vtx_zpos_temp, dimuon2vtx_zpos_temp;
                dimuon1vtx_zpos_temp.clear(); dimuon2vtx_zpos_temp.clear();
                
                std::vector<double> dimuon1vtx_xposError_temp, dimuon2vtx_xposError_temp;
                dimuon1vtx_xposError_temp.clear(); dimuon2vtx_xposError_temp.clear();
                
                std::vector<double> dimuon1vtx_yposError_temp, dimuon2vtx_yposError_temp;
                dimuon1vtx_yposError_temp.clear(); dimuon2vtx_yposError_temp.clear();
                 
                std::vector<double> dimuon1vtx_zposError_temp, dimuon2vtx_zposError_temp;
                dimuon1vtx_zposError_temp.clear(); dimuon2vtx_zposError_temp.clear();
                 
                //same idea for dimuon1vtx_xpos, ypos, zpos, and for their corresponding errors
               
                std::cout << "DEBUG 4" << std::endl;
                
                //match_12_34_56 case, aka ZplusY case
                if (match_12_34_56) {
                  if (vtx12.isValid()) {
                    std::cout << "DEBUG 5" << std::endl;
                    dimuon1vtx_temp.push_back(TMath::Prob(vtx12.totalChiSquared(), int(vtx12.degreesOfFreedom())));
                 //   std::cout << "DEBUG 5A" << std::endl;
                 //   std::cout << TMath::Prob(vtx12.totalChiSquared(), int(vtx12.degreesOfFreedom())) << std::endl;
                    dimuon1vtx_xpos_temp.push_back(vtx12.position().x()); dimuon1vtx_xposError_temp.push_back(vtx12.positionError().cxx());
                //    std::cout << "DEBUG 5B" << std::endl;
                //    std::cout << vtx12.position().x() << std::endl; 
                    dimuon1vtx_ypos_temp.push_back(vtx12.position().y()); dimuon1vtx_yposError_temp.push_back(vtx12.positionError().cyy());
                //    std::cout << "DEBUG 5C" << std::endl;
                 //   std::cout << vtx12.position().y() << std::endl; 
                    dimuon1vtx_zpos_temp.push_back(vtx12.position().z()); dimuon1vtx_zposError_temp.push_back(vtx12.positionError().czz());
                  } 
                  else{
                    dimuon1vtx_temp.push_back(-1000.);
                    dimuon1vtx_xpos_temp.push_back(-1000.);
                    dimuon1vtx_ypos_temp.push_back(-1000.);
                    dimuon1vtx_zpos_temp.push_back(-1000.);
                    dimuon1vtx_xposError_temp.push_back(-1000.);
                    dimuon1vtx_yposError_temp.push_back(-1000.);
                    dimuon1vtx_zposError_temp.push_back(-1000.);
                    
                    }

                  if (vtx34.isValid()) {
                    std::cout << "DEBUG 6" << std::endl;
                    dimuon2vtx_temp.push_back(TMath::Prob(vtx34.totalChiSquared(), int(vtx34.degreesOfFreedom())));
                    dimuon2vtx_xpos_temp.push_back(vtx34.position().x()); dimuon2vtx_xposError_temp.push_back(vtx34.positionError().cxx());
                    dimuon2vtx_ypos_temp.push_back(vtx34.position().y()); dimuon2vtx_yposError_temp.push_back(vtx34.positionError().cyy());
                    dimuon2vtx_zpos_temp.push_back(vtx34.position().z()); dimuon2vtx_zposError_temp.push_back(vtx34.positionError().czz());
                  } 
                  else{
                    dimuon2vtx_temp.push_back(-1000.);
                    dimuon2vtx_xpos_temp.push_back(-1000.);
                    dimuon2vtx_ypos_temp.push_back(-1000.);
                    dimuon2vtx_zpos_temp.push_back(-1000.);
                    dimuon2vtx_xposError_temp.push_back(-1000.);
                    dimuon2vtx_yposError_temp.push_back(-1000.);
                    dimuon2vtx_zposError_temp.push_back(-1000.);
                   }

                }
                else {
                dimuon1vtx_temp.push_back(-10000);
                dimuon2vtx_temp.push_back(-10000);
                dimuon1vtx_xpos_temp.push_back(-10000.);
                dimuon2vtx_xpos_temp.push_back(-10000.);
                dimuon1vtx_ypos_temp.push_back(-10000.);
                dimuon2vtx_ypos_temp.push_back(-10000.);
                dimuon1vtx_zpos_temp.push_back(-10000.);
                dimuon2vtx_zpos_temp.push_back(-10000.);
                dimuon1vtx_xposError_temp.push_back(-10000.);
                dimuon2vtx_xposError_temp.push_back(-10000.);
                dimuon1vtx_yposError_temp.push_back(-10000.);
                dimuon2vtx_yposError_temp.push_back(-10000.);
                dimuon1vtx_zposError_temp.push_back(-10000.);
                dimuon2vtx_zposError_temp.push_back(-10000.);
                

                }
                // <<<<<<<<<-------------------------->>>>>>>>>>>>>
                //match_13_24_56 case, aka ZplusY case
                if (match_13_24_56) {
                  if (vtx13.isValid()) {
                    std::cout << "DEBUG 7" << std::endl;
                    dimuon1vtx_temp.push_back(TMath::Prob(vtx13.totalChiSquared(), int(vtx13.degreesOfFreedom())));
                    dimuon1vtx_xpos_temp.push_back(vtx13.position().x()); dimuon1vtx_xposError_temp.push_back(vtx13.positionError().cxx());
                    dimuon1vtx_ypos_temp.push_back(vtx13.position().y()); dimuon1vtx_yposError_temp.push_back(vtx13.positionError().cyy());
                    dimuon1vtx_zpos_temp.push_back(vtx13.position().z()); dimuon1vtx_zposError_temp.push_back(vtx13.positionError().czz());
                  } 
                  else{
                    dimuon1vtx_temp.push_back(-1000.);
                    dimuon1vtx_xpos_temp.push_back(-1000.);
                    dimuon1vtx_ypos_temp.push_back(-1000.);
                    dimuon1vtx_zpos_temp.push_back(-1000.);
                    dimuon1vtx_xposError_temp.push_back(-1000.);
                    dimuon1vtx_yposError_temp.push_back(-1000.);
                    dimuon1vtx_zposError_temp.push_back(-1000.);
                    
                   }

                  if (vtx24.isValid()) {
                    std::cout << "DEBUG 8" << std::endl;
                    dimuon2vtx_temp.push_back(TMath::Prob(vtx24.totalChiSquared(), int(vtx24.degreesOfFreedom())));
                    dimuon2vtx_xpos_temp.push_back(vtx24.position().x()); dimuon2vtx_xposError_temp.push_back(vtx24.positionError().cxx());
                    dimuon2vtx_ypos_temp.push_back(vtx24.position().y()); dimuon2vtx_yposError_temp.push_back(vtx24.positionError().cyy());
                    dimuon2vtx_zpos_temp.push_back(vtx24.position().z()); dimuon2vtx_zposError_temp.push_back(vtx24.positionError().czz());
                  } 
                  else{
                    dimuon2vtx_temp.push_back(-1000.);
                    dimuon2vtx_xpos_temp.push_back(-1000);
                    dimuon2vtx_ypos_temp.push_back(-1000);
                    dimuon2vtx_zpos_temp.push_back(-1000);
                    dimuon2vtx_xposError_temp.push_back(-1000.);
                    dimuon2vtx_yposError_temp.push_back(-1000.);
                    dimuon2vtx_zposError_temp.push_back(-1000.);
                    
                   }

                }
                else {
                dimuon1vtx_temp.push_back(-10000);
                dimuon2vtx_temp.push_back(-10000);
                dimuon1vtx_xpos_temp.push_back(-10000.);
                dimuon2vtx_xpos_temp.push_back(-10000.);
                dimuon1vtx_ypos_temp.push_back(-10000.);
                dimuon2vtx_ypos_temp.push_back(-10000.);
                dimuon1vtx_zpos_temp.push_back(-10000.);
                dimuon2vtx_zpos_temp.push_back(-10000.);
                dimuon1vtx_xposError_temp.push_back(-10000.);
                dimuon2vtx_xposError_temp.push_back(-10000.);
                dimuon1vtx_yposError_temp.push_back(-10000.);
                dimuon2vtx_yposError_temp.push_back(-10000.);
                dimuon1vtx_zposError_temp.push_back(-10000.);
                dimuon2vtx_zposError_temp.push_back(-10000.);
                }
                // <<<<<<<<<-------------------------->>>>>>>>>>>>>
                //match_14_23_56 case, aka ZplusY case
                if (match_14_23_56) {
                  if (vtx14.isValid()) {
                    std::cout << "DEBUG 9" << std::endl;
                    dimuon1vtx_temp.push_back(TMath::Prob(vtx14.totalChiSquared(), int(vtx14.degreesOfFreedom())));
                    dimuon1vtx_xpos_temp.push_back(vtx14.position().x()); dimuon1vtx_xposError_temp.push_back(vtx14.positionError().cxx());
                    dimuon1vtx_ypos_temp.push_back(vtx14.position().y()); dimuon1vtx_yposError_temp.push_back(vtx14.positionError().cyy());
                    dimuon1vtx_zpos_temp.push_back(vtx14.position().z()); dimuon1vtx_zposError_temp.push_back(vtx14.positionError().czz());
                  } 
                  else{
                    dimuon1vtx_temp.push_back(-1000.);
                    dimuon1vtx_xpos_temp.push_back(-1000.);
                    dimuon1vtx_ypos_temp.push_back(-1000.);
                    dimuon1vtx_zpos_temp.push_back(-1000.);
                    dimuon1vtx_xposError_temp.push_back(-1000.);
                    dimuon1vtx_yposError_temp.push_back(-1000.);
                    dimuon1vtx_zposError_temp.push_back(-1000.);
                    
                    }

                  if (vtx23.isValid()) {
                    std::cout << "DEBUG 10" << std::endl;
                    dimuon2vtx_temp.push_back(TMath::Prob(vtx23.totalChiSquared(), int(vtx23.degreesOfFreedom())));
                    dimuon2vtx_xpos_temp.push_back(vtx23.position().x()); dimuon2vtx_xposError_temp.push_back(vtx23.positionError().cxx());
                    dimuon2vtx_ypos_temp.push_back(vtx23.position().y()); dimuon2vtx_yposError_temp.push_back(vtx23.positionError().cyy());
                    dimuon2vtx_zpos_temp.push_back(vtx23.position().z()); dimuon2vtx_zposError_temp.push_back(vtx23.positionError().czz());
                  } 
                  else{
                    dimuon2vtx_temp.push_back(-1000.);
                    dimuon2vtx_xpos_temp.push_back(-1000.);
                    dimuon2vtx_ypos_temp.push_back(-1000.);
                    dimuon2vtx_zpos_temp.push_back(-1000.);
                    dimuon2vtx_xposError_temp.push_back(-1000.);
                    dimuon2vtx_yposError_temp.push_back(-1000.);
                    dimuon2vtx_zposError_temp.push_back(-1000.);
                    }

                }
                else {
                dimuon1vtx_temp.push_back(-10000);
                dimuon2vtx_temp.push_back(-10000);
                dimuon1vtx_xpos_temp.push_back(-10000);
                dimuon2vtx_xpos_temp.push_back(-10000);
                dimuon1vtx_ypos_temp.push_back(-10000);
                dimuon2vtx_ypos_temp.push_back(-10000);
                dimuon1vtx_zpos_temp.push_back(-10000);
                dimuon2vtx_zpos_temp.push_back(-10000);
                dimuon1vtx_xposError_temp.push_back(-10000.);
                dimuon2vtx_xposError_temp.push_back(-10000.);
                dimuon1vtx_yposError_temp.push_back(-10000.);
                dimuon2vtx_yposError_temp.push_back(-10000.);
                dimuon1vtx_zposError_temp.push_back(-10000.);
                dimuon2vtx_zposError_temp.push_back(-10000.);
                }
                // <<<<<<<<<-------------------------->>>>>>>>>>>>>
                //match_12_34_ZOnly case
                if (match_12_34_ZOnly) {
                  if (vtx12.isValid()) {
                  //  std::cout << "DEBUG 5" << std::endl;
                    dimuon1vtx_temp.push_back(TMath::Prob(vtx12.totalChiSquared(), int(vtx12.degreesOfFreedom())));
                    dimuon1vtx_xpos_temp.push_back(vtx12.position().x()); dimuon1vtx_xposError_temp.push_back(vtx12.positionError().cxx());
                    dimuon1vtx_ypos_temp.push_back(vtx12.position().y()); dimuon1vtx_yposError_temp.push_back(vtx12.positionError().cyy());
                    dimuon1vtx_zpos_temp.push_back(vtx12.position().z()); dimuon1vtx_zposError_temp.push_back(vtx12.positionError().czz());
                  } 
                  else{
                    dimuon1vtx_temp.push_back(-1000.);
                    dimuon1vtx_xpos_temp.push_back(-1000.);
                    dimuon1vtx_ypos_temp.push_back(-1000.);
                    dimuon1vtx_zpos_temp.push_back(-1000.);
                    dimuon1vtx_xposError_temp.push_back(-1000.);
                    dimuon1vtx_yposError_temp.push_back(-1000.);
                    dimuon1vtx_zposError_temp.push_back(-1000.);
                    
                    }

                  if (vtx34.isValid()) {
                //    std::cout << "DEBUG 6" << std::endl;
                    dimuon2vtx_temp.push_back(TMath::Prob(vtx34.totalChiSquared(), int(vtx34.degreesOfFreedom())));
                    dimuon2vtx_xpos_temp.push_back(vtx34.position().x()); dimuon2vtx_xposError_temp.push_back(vtx34.positionError().cxx());
                    dimuon2vtx_ypos_temp.push_back(vtx34.position().y()); dimuon2vtx_yposError_temp.push_back(vtx34.positionError().cyy());
                    dimuon2vtx_zpos_temp.push_back(vtx34.position().z()); dimuon2vtx_zposError_temp.push_back(vtx34.positionError().czz());
                  } 
                  else{
                    dimuon2vtx_temp.push_back(-1000.);
                    dimuon2vtx_xpos_temp.push_back(-1000.);
                    dimuon2vtx_ypos_temp.push_back(-1000.);
                    dimuon2vtx_zpos_temp.push_back(-1000.);
                    dimuon2vtx_xposError_temp.push_back(-1000.);
                    dimuon2vtx_yposError_temp.push_back(-1000.);
                    dimuon2vtx_zposError_temp.push_back(-1000.);
                   }

                }
                else {
                dimuon1vtx_temp.push_back(-10000);
                dimuon2vtx_temp.push_back(-10000);
                dimuon1vtx_xpos_temp.push_back(-10000.);
                dimuon2vtx_xpos_temp.push_back(-10000.);
                dimuon1vtx_ypos_temp.push_back(-10000.);
                dimuon2vtx_ypos_temp.push_back(-10000.);
                dimuon1vtx_zpos_temp.push_back(-10000.);
                dimuon2vtx_zpos_temp.push_back(-10000.);
                dimuon1vtx_xposError_temp.push_back(-10000.);
                dimuon2vtx_xposError_temp.push_back(-10000.);
                dimuon1vtx_yposError_temp.push_back(-10000.);
                dimuon2vtx_yposError_temp.push_back(-10000.);
                dimuon1vtx_zposError_temp.push_back(-10000.);
                dimuon2vtx_zposError_temp.push_back(-10000.);
                

                }
                // <<<<<<<<<-------------------------->>>>>>>>>>>>>
                //match_13_24_ZOnly case
                if (match_13_24_ZOnly) {
                  if (vtx13.isValid()) {
                   // std::cout << "DEBUG 7" << std::endl;
                    dimuon1vtx_temp.push_back(TMath::Prob(vtx13.totalChiSquared(), int(vtx13.degreesOfFreedom())));
                    dimuon1vtx_xpos_temp.push_back(vtx13.position().x()); dimuon1vtx_xposError_temp.push_back(vtx13.positionError().cxx());
                    dimuon1vtx_ypos_temp.push_back(vtx13.position().y()); dimuon1vtx_yposError_temp.push_back(vtx13.positionError().cyy());
                    dimuon1vtx_zpos_temp.push_back(vtx13.position().z()); dimuon1vtx_zposError_temp.push_back(vtx13.positionError().czz());
                  } 
                  else{
                    dimuon1vtx_temp.push_back(-1000.);
                    dimuon1vtx_xpos_temp.push_back(-1000.);
                    dimuon1vtx_ypos_temp.push_back(-1000.);
                    dimuon1vtx_zpos_temp.push_back(-1000.);
                    dimuon1vtx_xposError_temp.push_back(-1000.);
                    dimuon1vtx_yposError_temp.push_back(-1000.);
                    dimuon1vtx_zposError_temp.push_back(-1000.);
                    
                   }

                  if (vtx24.isValid()) {
                  //  std::cout << "DEBUG 8" << std::endl;
                    dimuon2vtx_temp.push_back(TMath::Prob(vtx24.totalChiSquared(), int(vtx24.degreesOfFreedom())));
                    dimuon2vtx_xpos_temp.push_back(vtx24.position().x()); dimuon2vtx_xposError_temp.push_back(vtx24.positionError().cxx());
                    dimuon2vtx_ypos_temp.push_back(vtx24.position().y()); dimuon2vtx_yposError_temp.push_back(vtx24.positionError().cyy());
                    dimuon2vtx_zpos_temp.push_back(vtx24.position().z()); dimuon2vtx_zposError_temp.push_back(vtx24.positionError().czz());
                  } 
                  else{
                    dimuon2vtx_temp.push_back(-1000.);
                    dimuon2vtx_xpos_temp.push_back(-1000);
                    dimuon2vtx_ypos_temp.push_back(-1000);
                    dimuon2vtx_zpos_temp.push_back(-1000);
                    dimuon2vtx_xposError_temp.push_back(-1000.);
                    dimuon2vtx_yposError_temp.push_back(-1000.);
                    dimuon2vtx_zposError_temp.push_back(-1000.);
                    
                   }

                }
                else {
                dimuon1vtx_temp.push_back(-10000);
                dimuon2vtx_temp.push_back(-10000);
                dimuon1vtx_xpos_temp.push_back(-10000.);
                dimuon2vtx_xpos_temp.push_back(-10000.);
                dimuon1vtx_ypos_temp.push_back(-10000.);
                dimuon2vtx_ypos_temp.push_back(-10000.);
                dimuon1vtx_zpos_temp.push_back(-10000.);
                dimuon2vtx_zpos_temp.push_back(-10000.);
                dimuon1vtx_xposError_temp.push_back(-10000.);
                dimuon2vtx_xposError_temp.push_back(-10000.);
                dimuon1vtx_yposError_temp.push_back(-10000.);
                dimuon2vtx_yposError_temp.push_back(-10000.);
                dimuon1vtx_zposError_temp.push_back(-10000.);
                dimuon2vtx_zposError_temp.push_back(-10000.);
                }
                // <<<<<<<<<-------------------------->>>>>>>>>>>>>
                //match_14_23_ZOnly case
                if (match_14_23_ZOnly) {
                  if (vtx14.isValid()) {
     //               std::cout << "DEBUG 9" << std::endl;
                    dimuon1vtx_temp.push_back(TMath::Prob(vtx14.totalChiSquared(), int(vtx14.degreesOfFreedom())));
                    dimuon1vtx_xpos_temp.push_back(vtx14.position().x()); dimuon1vtx_xposError_temp.push_back(vtx14.positionError().cxx());
                    dimuon1vtx_ypos_temp.push_back(vtx14.position().y()); dimuon1vtx_yposError_temp.push_back(vtx14.positionError().cyy());
                    dimuon1vtx_zpos_temp.push_back(vtx14.position().z()); dimuon1vtx_zposError_temp.push_back(vtx14.positionError().czz());
                  } 
                  else{
                    dimuon1vtx_temp.push_back(-1000.);
                    dimuon1vtx_xpos_temp.push_back(-1000.);
                    dimuon1vtx_ypos_temp.push_back(-1000.);
                    dimuon1vtx_zpos_temp.push_back(-1000.);
                    dimuon1vtx_xposError_temp.push_back(-1000.);
                    dimuon1vtx_yposError_temp.push_back(-1000.);
                    dimuon1vtx_zposError_temp.push_back(-1000.);
                    
                    }

                  if (vtx23.isValid()) {
                //    std::cout << "DEBUG 10" << std::endl;
                    dimuon2vtx_temp.push_back(TMath::Prob(vtx23.totalChiSquared(), int(vtx23.degreesOfFreedom())));
                    dimuon2vtx_xpos_temp.push_back(vtx23.position().x()); dimuon2vtx_xposError_temp.push_back(vtx23.positionError().cxx());
                    dimuon2vtx_ypos_temp.push_back(vtx23.position().y()); dimuon2vtx_yposError_temp.push_back(vtx23.positionError().cyy());
                    dimuon2vtx_zpos_temp.push_back(vtx23.position().z()); dimuon2vtx_zposError_temp.push_back(vtx23.positionError().czz());
                  } 
                  else{
                    dimuon2vtx_temp.push_back(-1000.);
                    dimuon2vtx_xpos_temp.push_back(-1000.);
                    dimuon2vtx_ypos_temp.push_back(-1000.);
                    dimuon2vtx_zpos_temp.push_back(-1000.);
                    dimuon2vtx_xposError_temp.push_back(-1000.);
                    dimuon2vtx_yposError_temp.push_back(-1000.);
                    dimuon2vtx_zposError_temp.push_back(-1000.);
                    }

                }
                else {
                dimuon1vtx_temp.push_back(-10000);
                dimuon2vtx_temp.push_back(-10000);
                dimuon1vtx_xpos_temp.push_back(-10000);
                dimuon2vtx_xpos_temp.push_back(-10000);
                dimuon1vtx_ypos_temp.push_back(-10000);
                dimuon2vtx_ypos_temp.push_back(-10000);
                dimuon1vtx_zpos_temp.push_back(-10000);
                dimuon2vtx_zpos_temp.push_back(-10000);
                dimuon1vtx_xposError_temp.push_back(-10000.);
                dimuon2vtx_xposError_temp.push_back(-10000.);
                dimuon1vtx_yposError_temp.push_back(-10000.);
                dimuon2vtx_yposError_temp.push_back(-10000.);
                dimuon1vtx_zposError_temp.push_back(-10000.);
                dimuon2vtx_zposError_temp.push_back(-10000.);
                }
                
                
                
                
                
                //Finally, fill the dimuon1vtx, dimuon2vtx vectors with vectors
                dimuon1vtx.push_back(dimuon1vtx_temp);
                dimuon2vtx.push_back(dimuon2vtx_temp);
                dimuon1vtx_xpos.push_back(dimuon1vtx_xpos_temp);
                dimuon2vtx_xpos.push_back(dimuon2vtx_xpos_temp);
                dimuon1vtx_ypos.push_back(dimuon1vtx_ypos_temp);
                dimuon2vtx_ypos.push_back(dimuon2vtx_ypos_temp);
                dimuon1vtx_zpos.push_back(dimuon1vtx_zpos_temp);
                dimuon2vtx_zpos.push_back(dimuon2vtx_zpos_temp);
                dimuon1vtx_xposError.push_back(dimuon1vtx_xposError_temp);
                dimuon2vtx_xposError.push_back(dimuon2vtx_xposError_temp);
                dimuon1vtx_yposError.push_back(dimuon1vtx_yposError_temp);
                dimuon2vtx_yposError.push_back(dimuon2vtx_yposError_temp);
                dimuon1vtx_zposError.push_back(dimuon1vtx_zposError_temp);
                dimuon2vtx_zposError.push_back(dimuon2vtx_zposError_temp);
                
                
                //protection here in case dimuon1vtx_xpos does not exist. without this protection, code fails with an unhelpful complaint about vector sizes. S.L. taught me with these kinds of errors, check the vectors first and see what is unprotected (like this was)
                //I also put in an analogous protection for diumuon2vtx_xpos
            //     if (dimuon1vtx_xpos.size()==0 || dimuon2vtx_xpos.size()==0){
//                     std::cout << "diumuon1 or 2 vtx failure!" << std::endl;
//                     continue;
//                 } //WARNING: DO NOT REINSTATE THIS CONTINUE STATEMENT, IT CAUSES A BUG WHERE SOMETIMES WE GET EVENTS WTIH NO QUADS IN THEM. DO IT AS DONE BELOW, WRAPPING EVERYTHING IN THE IF STATEMENT, INSTEAD!!
//              

              // if (dimuon1vtx_xpos.size() != 0 && dimuon2vtx_xpos.size() !=0 && dimuon1vtx_ypos.size() !=0 && dimuon2vtx_ypos.size() !=0){
// 
//                 // Lxy part
//                 // dimuon1lxyctauPV, dimuon1lxy, dimuon1lxysig
//                 //calculate pseudo proper decay time
//                 TVector3 pvtx;
//                 pvtx.SetXYZ(PVx.at(0),PVy.at(0),0);
//                 TVector3 vtx1, vtx2;
//                 vtx1.SetXYZ(dimuon1vtx_xpos.at(0),dimuon1vtx_ypos.at(0),0);
//                 vtx2.SetXYZ(dimuon2vtx_xpos.at(0),dimuon2vtx_ypos.at(0),0);
//                 TVector3 vdiff1 = vtx1 - pvtx;
//                 TVector3 vdiff2 = vtx2 - pvtx;
//                 VertexDistanceXY vdistXY1, vdistXY2; //look up VertexDistanceXY documentation here: https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/plugins/VertexTableProducer.cc, that is for NanoAOD but likely the right idea
//                 Measurement1D distXY1, distXY2;
//                 double cosAlphaXY1, cosAlphaXY2;
//                 if (match_13_24_56 && vtx13.isValid() && vtx24.isValid()) {
//                   TVector3 pperp1((lepton1+lepton3).px(), (lepton1+lepton3).py(), 0);
//                   TVector3 pperp2((lepton2+lepton4).px(), (lepton2+lepton4).py(), 0);
//                   std::cout << "DEBUG 11" << std::endl;
//                   distXY1 = vdistXY1.distance(vtx13.vertexState(), primary_vertex);
//                   std::cout << "DEBUG 12" << std::endl;
//                   distXY2 = vdistXY2.distance(vtx24.vertexState(), primary_vertex);
//                   cosAlphaXY1 = vdiff1.Dot(pperp1)/(vdiff1.Perp()*pperp1.Perp());
//                   cosAlphaXY2 = vdiff2.Dot(pperp2)/(vdiff2.Perp()*pperp2.Perp());
//                   std::cout << "DEBUG 13" << std::endl;
//                   dimuon1lxyctauPV.push_back (distXY1.value()*cosAlphaXY1 * (lepton1+lepton3).mass()/pperp1.Perp()); 
//                   dimuon2lxyctauPV.push_back (distXY2.value()*cosAlphaXY2 * (lepton2+lepton4).mass()/pperp2.Perp());
//                   dimuon1lxy.push_back(distXY1.value()); dimuon1lxysig.push_back(distXY1.value()/distXY1.error());
//                   dimuon2lxy.push_back(distXY2.value()); dimuon2lxysig.push_back(distXY2.value()/distXY2.error());
// 
//                 }
//                 if (match_14_23_56 && vtx14.isValid() && vtx23.isValid()) {
//                   TVector3 pperp1((lepton1+lepton4).px(), (lepton1+lepton4).py(), 0);
//                   TVector3 pperp2((lepton2+lepton3).px(), (lepton2+lepton3).py(), 0);
//                   std::cout << "DEBUG 14" << std::endl;
//                   distXY1 = vdistXY1.distance(vtx14.vertexState(), primary_vertex);
//                   std::cout << "DEBUG 15" << std::endl;
//                   distXY2 = vdistXY2.distance(vtx23.vertexState(), primary_vertex);
//                   std::cout << "DEBUG 16" << std::endl;
//                   cosAlphaXY1 = vdiff1.Dot(pperp1)/(vdiff1.Perp()*pperp1.Perp());
//                   cosAlphaXY2 = vdiff2.Dot(pperp2)/(vdiff2.Perp()*pperp2.Perp());
//                   dimuon1lxyctauPV.push_back (distXY1.value()*cosAlphaXY1 * (lepton1+lepton4).mass()/pperp1.Perp()); 
//                   dimuon2lxyctauPV.push_back (distXY2.value()*cosAlphaXY2 * (lepton2+lepton3).mass()/pperp2.Perp());
//                   dimuon1lxy.push_back(distXY1.value()); dimuon1lxysig.push_back(distXY1.value()/distXY1.error());
//                   dimuon2lxy.push_back(distXY2.value()); dimuon2lxysig.push_back(distXY2.value()/distXY2.error());
//                 }
//                 
//                 if (match_12_34_56 && vtx12.isValid() && vtx34.isValid()) {
//                   TVector3 pperp1((lepton1+lepton2).px(), (lepton1+lepton2).py(), 0); 
//                   TVector3 pperp2((lepton3+lepton4).px(), (lepton3+lepton4).py(), 0);
//                   std::cout << "DEBUG 17" << std::endl;
//                   distXY1 = vdistXY1.distance(vtx12.vertexState(), primary_vertex);
//                   std::cout << "DEBUG 18" << std::endl;
//                   distXY2 = vdistXY2.distance(vtx34.vertexState(), primary_vertex);
//                   std::cout << "DEBUG 19" << std::endl;
//                   cosAlphaXY1 = vdiff1.Dot(pperp1)/(vdiff1.Perp()*pperp1.Perp());
//                   cosAlphaXY2 = vdiff2.Dot(pperp2)/(vdiff2.Perp()*pperp2.Perp());
//                   dimuon1lxyctauPV.push_back (distXY1.value()*cosAlphaXY1 * (lepton1+lepton2).mass()/pperp1.Perp()); 
//                   dimuon2lxyctauPV.push_back (distXY2.value()*cosAlphaXY2 * (lepton3+lepton4).mass()/pperp2.Perp());
//                   dimuon1lxy.push_back(distXY1.value()); dimuon1lxysig.push_back(distXY1.value()/distXY1.error());
//                   dimuon2lxy.push_back(distXY2.value()); dimuon2lxysig.push_back(distXY2.value()/distXY2.error());
//                 }
//             }
            
             //loop on tracks
//            bool matchedPV1 = false; bool matchedPV2 = false; bool matchedPV3 = false; bool matchedPV4 = false;
            double lepton1frompv = -10; double lepton2frompv = -10; double lepton3frompv = -10; double lepton4frompv = -10;
            double lepton1pvassq = -10; double lepton2pvassq = -10; double lepton3pvassq = -10; double lepton4pvassq = -10;
    //        int lepton1PVIndex = -1; int lepton2PVIndex = -1; int lepton3PVIndex = -1;
            
           for (auto itrk = pfs->begin(); itrk != pfs->end(); ++itrk) {
            //  std::cout << "itrk type:  " << typeid(itrk).name();
            //learned how to check type from: https://stackoverflow.com/questions/11310898/how-do-i-get-the-type-of-a-variable
              if (TMath::Abs(itrk->pdgId())==13 && itrk->fromPV()>1) {
                if (TMath::Abs(iM1->pt() - itrk->pt()) <.5 && TMath::Abs(iM1->eta() - itrk->eta()  ) < .2) {
//                  std::cout << "iM1: " << itrk->pt() << " " << iM1->pt() << " " << itrk->eta() << " " << iM1->eta() << " " << itrk->fromPV() << " " << itrk->pvAssociationQuality() << "\n";
                  lepton1frompv = itrk->fromPV();
                  lepton1pvassq = itrk->pvAssociationQuality();
                 }
                if (TMath::Abs(iM2->pt() - itrk->pt()) <.5 && TMath::Abs( iM2->eta() - itrk->eta() ) < .2) {
//                  std::cout << "iM2: " << itrk->pt() << " " << iM2->pt() << " " << itrk->eta() << " " << iM2->eta() << " " << itrk->fromPV() << " " << itrk->pvAssociationQuality() << "\n";
                  lepton2frompv = itrk->fromPV();
                  lepton2pvassq = itrk->pvAssociationQuality();
                  
                }
                if (TMath::Abs(iM3->pt() - itrk->pt()) <.5 && TMath::Abs( iM3->eta() - itrk->eta() ) < .2) {
//                  std::cout << "iM3: " << itrk->pt() << " " << iM3->pt() << " " << itrk->eta() << " " << iM3->eta() << " " << itrk->fromPV() << " " << itrk->pvAssociationQuality() << "\n";
                  lepton3frompv = itrk->fromPV();
                  lepton3pvassq = itrk->pvAssociationQuality();
                }
                if (TMath::Abs(iM4->pt() - itrk->pt()) <.5 && TMath::Abs( iM4->eta() - itrk->eta() ) < .2) {
//                  std::cout << "iM4: " << itrk->pt() << " " << iM4->pt() << " " << itrk->eta() << " " << iM4->eta() << " " << itrk->fromPV() << " " << itrk->pvAssociationQuality() << "\n";
                  lepton4frompv = itrk->fromPV();
                  lepton4pvassq = itrk->pvAssociationQuality();
                }
              }  
            }
                
                
                lepton1_fromPV.push_back(lepton1frompv); lepton1_pvassq.push_back(lepton1pvassq); 
                lepton2_fromPV.push_back(lepton2frompv); lepton2_pvassq.push_back(lepton2pvassq);
                lepton3_fromPV.push_back(lepton3frompv); lepton3_pvassq.push_back(lepton3pvassq);
                lepton4_fromPV.push_back(lepton4frompv); lepton4_pvassq.push_back(lepton4pvassq);
                
                pair_12_34_56.push_back(match_12_34_56);
                pair_13_24_56.push_back(match_13_24_56);
                pair_14_23_56.push_back(match_14_23_56);
                
                pair_12_34_ZOnly.push_back(match_12_34_ZOnly);
                pair_13_24_ZOnly.push_back(match_13_24_ZOnly);
                pair_14_23_ZOnly.push_back(match_14_23_ZOnly);

                event_number.push_back(iEvent.id().event());
                run_number.push_back(iEvent.run());
                lumi_section.push_back(iEvent.luminosityBlock());
                bunch_crossing.push_back(iEvent.bunchCrossing());
                orbit_number.push_back(iEvent.orbitNumber());
  //              process_ID.push_back(iEvent.id());
  
                //Rho, the median energy desnity //

                edm::Handle< double > rhoH;
                iEvent.getByToken(rhoToken_,rhoH);
              //  float rho=*rhoH;
             // rho.clear();
                rho.push_back(*rhoH);
                
                //making sure this is right...
                //quadHasHowManyTrigMatches.push_back(totalTrigMatched);
                quadHasHowManyTrigMatches.push_back(trigMatched);
                quadComesFromEvtWithThisManyMu.push_back((int)muons->size()); // this isn't exactly the helpful number we want, we want to know how many GOOD (passing our cuts) muons are in a given event
                sum4LepCharges.push_back(iM1->charge() + iM2->charge() + iM3->charge() + iM4->charge());
                
                flagZplusY                         .push_back(isZplusY);
                flagZOnly                          .push_back(isZOnly);
                
                lepton1_pt                         .push_back(iM1->pt());  
                lepton1_eta                        .push_back(iM1->eta());
                lepton1_phi                        .push_back(iM1->phi());
                lepton1_charge                     .push_back(iM1->charge());
                lepton1_dxy                        .push_back(iM1->muonBestTrack()->dxy(primary_vertex.position()));
                lepton1_dz                         .push_back(iM1->muonBestTrack()->dz(primary_vertex.position()));
                lepton1_d0                         .push_back(iM1->muonBestTrack()->d0());
                lepton1_iso03particle              .push_back(iM1->pfIsolationR03().sumChargedParticlePt);
                lepton1_iso03hadron                .push_back(iM1->pfIsolationR03().sumChargedHadronPt);
                lepton1_iso03neutralHadron         .push_back(iM1->pfIsolationR03().sumNeutralHadronEt);
                lepton1_iso03photon                .push_back(iM1->pfIsolationR03().sumPhotonEt);
                lepton1_iso03PU                    .push_back(iM1->pfIsolationR03().sumPUPt);
                lepton1_iso04hadron                .push_back(iM1->pfIsolationR04().sumChargedHadronPt);
                lepton1_iso04particle              .push_back(iM1->pfIsolationR04().sumChargedParticlePt);
                
                lepton1_impactParameterSignificance.push_back(fabs(iM1->dB(pat::Muon::PV3D)/iM1->edB(pat::Muon::PV3D)));
                lepton1_isLooseMuon                .push_back(iM1->isLooseMuon());
                lepton1_isMediumMuon               .push_back(iM1->isMediumMuon());
                lepton1_isSoftMuon                 .push_back(iM1->isSoftMuon(primary_vertex));
                lepton1_isTightMuon                .push_back(iM1->isTightMuon(primary_vertex));
                lepton1_isPFMuon                   .push_back(iM1->isPFMuon());
                lepton1_isGlobalMuon               .push_back(iM1->isGlobalMuon());
                lepton1_isTrackerMuon              .push_back(iM1->isTrackerMuon());
            
                lepton2_pt                         .push_back(iM2->pt());  
                lepton2_eta                        .push_back(iM2->eta());
                lepton2_phi                        .push_back(iM2->phi());
                lepton2_charge                     .push_back(iM2->charge());
                lepton2_dxy                        .push_back(iM2->muonBestTrack()->dxy(primary_vertex.position()));
                lepton2_dz                         .push_back(iM2->muonBestTrack()->dz(primary_vertex.position()));
                lepton2_d0                         .push_back(iM2->muonBestTrack()->d0());
                lepton2_iso03particle              .push_back(iM2->pfIsolationR03().sumChargedParticlePt);
                lepton2_iso03hadron                .push_back(iM2->pfIsolationR03().sumChargedHadronPt);
                lepton2_iso03neutralHadron         .push_back(iM2->pfIsolationR03().sumNeutralHadronEt);
                lepton2_iso03photon                .push_back(iM2->pfIsolationR03().sumPhotonEt);
                lepton2_iso03PU                    .push_back(iM2->pfIsolationR03().sumPUPt);
                lepton2_iso04hadron                .push_back(iM2->pfIsolationR04().sumChargedHadronPt);
                lepton2_iso04particle              .push_back(iM2->pfIsolationR04().sumChargedParticlePt);
                lepton2_impactParameterSignificance.push_back(fabs(iM2->dB(pat::Muon::PV3D)/iM2->edB(pat::Muon::PV3D)));
                lepton2_isLooseMuon                .push_back(iM2->isLooseMuon());
                lepton2_isMediumMuon               .push_back(iM2->isMediumMuon());
                lepton2_isSoftMuon                 .push_back(iM2->isSoftMuon(primary_vertex));
                lepton2_isTightMuon                .push_back(iM2->isTightMuon(primary_vertex));
                lepton2_isPFMuon                   .push_back(iM2->isPFMuon());
                lepton2_isGlobalMuon               .push_back(iM2->isGlobalMuon());
                lepton2_isTrackerMuon              .push_back(iM2->isTrackerMuon());
        
                lepton3_pt                         .push_back(iM3->pt());  
                lepton3_eta                        .push_back(iM3->eta());
                lepton3_phi                        .push_back(iM3->phi());
                lepton3_charge                     .push_back(iM3->charge());
                lepton3_dxy                        .push_back(iM3->muonBestTrack()->dxy(primary_vertex.position()));
                lepton3_dz                         .push_back(iM3->muonBestTrack()->dz(primary_vertex.position()));
                lepton3_d0                         .push_back(iM3->muonBestTrack()->d0());
                lepton3_iso03particle              .push_back(iM3->pfIsolationR03().sumChargedParticlePt);
                lepton3_iso03hadron                .push_back(iM3->pfIsolationR03().sumChargedHadronPt);
                lepton3_iso03neutralHadron         .push_back(iM3->pfIsolationR03().sumNeutralHadronEt);
                lepton3_iso03photon                .push_back(iM3->pfIsolationR03().sumPhotonEt);
                lepton3_iso03PU                    .push_back(iM3->pfIsolationR03().sumPUPt);
                lepton3_iso04hadron                .push_back(iM3->pfIsolationR04().sumChargedHadronPt);
                lepton3_iso04particle              .push_back(iM3->pfIsolationR04().sumChargedParticlePt);
                lepton3_impactParameterSignificance.push_back(fabs(iM3->dB(pat::Muon::PV3D)/iM3->edB(pat::Muon::PV3D)));
                lepton3_isLooseMuon                .push_back(iM3->isLooseMuon());
                lepton3_isMediumMuon               .push_back(iM3->isMediumMuon());
                lepton3_isSoftMuon                 .push_back(iM3->isSoftMuon(primary_vertex));
                lepton3_isTightMuon                .push_back(iM3->isTightMuon(primary_vertex));
                lepton3_isPFMuon                   .push_back(iM3->isPFMuon());
                lepton3_isGlobalMuon               .push_back(iM3->isGlobalMuon());
                lepton3_isTrackerMuon              .push_back(iM3->isTrackerMuon());
            
                lepton4_pt                         .push_back(iM4->pt());  
                lepton4_eta                        .push_back(iM4->eta());
                lepton4_phi                        .push_back(iM4->phi());
                lepton4_charge                     .push_back(iM4->charge());
                lepton4_dxy                        .push_back(iM4->muonBestTrack()->dxy(primary_vertex.position()));
                lepton4_dz                         .push_back(iM4->muonBestTrack()->dz(primary_vertex.position()));
                lepton4_d0                         .push_back(iM4->muonBestTrack()->d0());
                lepton4_iso03particle              .push_back(iM4->pfIsolationR03().sumChargedParticlePt);
                lepton4_iso03hadron                .push_back(iM4->pfIsolationR03().sumChargedHadronPt);
                lepton4_iso03neutralHadron         .push_back(iM4->pfIsolationR03().sumNeutralHadronEt);
                lepton4_iso03photon                .push_back(iM4->pfIsolationR03().sumPhotonEt);
                lepton4_iso03PU                    .push_back(iM4->pfIsolationR03().sumPUPt);
                lepton4_iso04hadron                .push_back(iM4->pfIsolationR04().sumChargedHadronPt);
                lepton4_iso04particle              .push_back(iM4->pfIsolationR04().sumChargedParticlePt);
                lepton4_impactParameterSignificance.push_back(fabs(iM4->dB(pat::Muon::PV3D)/iM4->edB(pat::Muon::PV3D)));
                lepton4_isLooseMuon                .push_back(iM4->isLooseMuon());
                lepton4_isMediumMuon               .push_back(iM4->isMediumMuon());
                lepton4_isSoftMuon                 .push_back(iM4->isSoftMuon(primary_vertex));
                lepton4_isTightMuon                .push_back(iM4->isTightMuon(primary_vertex));
                lepton4_isPFMuon                   .push_back(iM4->isPFMuon());
                lepton4_isGlobalMuon               .push_back(iM4->isGlobalMuon());
                lepton4_isTrackerMuon              .push_back(iM4->isTrackerMuon());
                
                
                inv4MuMass                         .push_back((lepton1 + lepton2 + lepton3 + lepton4).mass());
                big4MuEta                          .push_back(eta_1_2_3_4);
                big4MuPhi                          .push_back(phi_1_2_3_4);
                big4MuPt                           .push_back(pt_1_2_3_4);
  
  
//Missing hit info, maybe not needed according to S.L., K.P. let me know that this NOT work because pat::Muon doesn't inherit from pat::PackedCandidate so this does not work            
  //              lepton1_missingHits                .push_back(iM1->muonBestTrack().lostInnerHits());
 //              lepton2_missingHits                .push_back(iM2->lostInnerHits());
//                lepton3_missingHits                .push_back(iM3->lostInnerHits());
 //               lepton4_missingHits                .push_back(iM4->lostInnerHits());
 
               lepton1_validHits                   .push_back(iM1->numberOfValidHits());
               lepton2_validHits                   .push_back(iM2->numberOfValidHits());
               lepton3_validHits                   .push_back(iM3->numberOfValidHits());
               lepton4_validHits                   .push_back(iM4->numberOfValidHits());

          }
        }
       }
     }
     if (save_event){
//       std::cout << "save_event for the data stuff: " << save_event << std::endl;
        save_event_count.push_back(1);
       tree->Fill();
      }
   }
//std::cout << "save_event_count: " << save_event_count << std::endl;
//std::cout << "save_event is: " << save_event << std::endl;
// **************
// MC starts here
// **************
   if (isMC){ // && save_event) { //change to (isMC && save_event) for testing 
//    std::cout << "Am here at check 2 " << std::endl; 
 
//*****
//PILEUP, observed and true 
//*****
   edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
   iEvent.getByToken(puToken_,PupInfo);
   std::vector<PileupSummaryInfo>::const_iterator ipu;   
   for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu){
      if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
      PU = ipu->getPU_NumInteractions();
      PU_True = ipu->getTrueNumInteractions();
       }
    histContainer_["PU"]->Fill(PU);
    histContainer_["PU_True"]->Fill(PU_True);
   
    
    loop_enter_check.push_back(1);
    mc_event_number.push_back(iEvent.id().event());
    mc_run_number.push_back(iEvent.run());
    mc_lumi_section.push_back(iEvent.luminosityBlock());
        
    edm::Handle<reco::GenParticleCollection> mc_particles;
    iEvent.getByToken(genParticlesToken_, mc_particles);

    int UPSI = 553;
    int Z = 23;
    int MUON = 13;
    
    int UPSI_2S = 100553;
    int UPSI_3S = 200553;
    int Z_to_fill_count = 0;
    int Upsi_to_fill_count =  0;
    
    
    //ChibJ_1P states
    int Chib0_1P = 10551;
    int Chib1_1P = 20553;
    int Chib2_1P = 555;
    
    //Photon
    int photon = 22;

    std::vector<double> upsi_pt_found;
    upsi_pt_found.clear();
    std::vector<double> Z_pt_found;
    Z_pt_found.clear();
    
    bool eventHasZToMuMu = false;
    std::cout << "eventHasZToMuMu initially is: " << eventHasZToMuMu << std::endl; 
    
    bool eventHasUpsiNToMuMu = false;
    std::cout << "eventHasUpsiNToMuMu initially is: " << eventHasUpsiNToMuMu << std::endl; 
    
    bool eventHasZUpsiNTo4Mu = false;
    std::cout << "eventHasZUpsiNTo4Mu initially is: " << eventHasZUpsiNTo4Mu << std::endl; 
    
//    std::cout << "upsi_pt_found.size() initially:  " << upsi_pt_found.size() <<std::endl; //Mary added for testing 
    //QUESTION!: I probably just missed something, but is there a reason you don't fill the upsi_pt_found until line 953 (similar question for Z_pt_found). I guess I'm not sure what this loop 
    //is doing, you said it was a protection against the pythia weirdos (when pythia has a lot of things with similar pT/eta) but I'm not sure I'm understanding what's going on here
    //wouldn't it make more sense to test the truth pT of the mc particle versus the reco pT of the mc particle?
 //   bool registered_upsi = false;
 //   bool registered_Z = false;

    for(unsigned int i = 0; i < mc_particles->size(); ++i) {
      const reco::GenParticle* gen_particle = &mc_particles->at(i);
//      if (gen_particle->pdgId() != UPSI && gen_particle->pdgId() != MUON && gen_particle->pdgId() != Z)
  //       std::cout << "NOT AN UPSI, MU, or Z" << std::endl;
  
         // get at final state muons that come from upsi1, upsi2, or upsi3 //verified that the statuts = 1 requirement is too harsh
         //also verified that when we relax it, we get the same answer working backwards from the muons 
         // as we do working downwards from the upsi
         
         
     //     if (TMath::Abs(gen_particle->pdgId()) == MUON){// && gen_particle->status() == 1){
//           // std::cout << "SECOND PLACEHOLDER" << std::endl;
//            bool hasUpsi1Ancestor = false;
//            for (size_t k = 0; k < gen_particle->numberOfMothers(); k++){
//            //  std::cout << "THIRD PLACE HOLDER" << std::endl;
//              if (gen_particle->mother(k)->pdgId() == UPSI){
//                hasUpsi1Ancestor = true;
//                std::cout << "hasUpsi1Ancestor  " << hasUpsi1Ancestor << std::endl;
//                std::cout << "gen_particle->mother(k)->pdgId()  " << gen_particle->mother(k)->pdgId() <<std::endl;
//                std::cout << "gen_particle->mother(k)->pt()     " << gen_particle->mother(k)->pt() << std::endl;
//                std::cout << "gen_particle->mother(k)->eta()    " << gen_particle->mother(k)->eta() << std::endl;
//                std::cout << "gen_particle->mother(k)->phi()    " << gen_particle->mother(k)->phi() << std::endl;
//                std::cout << "gen_particle->mother(k)->mass()   " << gen_particle->mother(k)->mass() << std::endl;
//                
//                std::cout << "gen_particle->pdgId()            " << gen_particle->pdgId() << std::endl;
//                std::cout << "gen_particle->pt()               " << gen_particle->pt() << std::endl;
//                std::cout << "gen_particle->eta()              " << gen_particle->eta() << std::endl;
//                std::cout << "gen_particle->phi()              " << gen_particle->phi() << std::endl;
//                std::cout << "gen_particle->mass()             " << gen_particle->mass() << std::endl;
//              
//              
//              }
//            }
//          }

         
         
         
         if (gen_particle->pdgId() == UPSI_2S ||gen_particle->pdgId() == UPSI_3S){
            std::cout << "FOUND UPSI_2S or UPSI_3S" << std::endl;
            std::cout << gen_particle->pdgId() << std::endl;
            std::cout << gen_particle->mass() << std::endl;
           // std::cout << gen_particle->daughter(0)->pdgId() << std::endl;
           // std::cout << gen_particle->daughter(1)->pdgId() << std::endl;
          //  std::cout << gen_particle->daughter(2)->pdgId() << std::endl;
          //test to induce new compilation
        }
        
         //UPSI state
         if (gen_particle->pdgId() == UPSI){
            std::cout << "gen_particle->pdgId() " << gen_particle->pdgId() << std::endl;
            std::cout << "gen_particle->status() " << gen_particle->status() << std::endl;
            
            int upsi_mu_daughter_counter = 0;
            for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
              if (TMath::Abs(gen_particle->daughter(j)->pdgId()) == MUON){ //&& gen_particle->daughter(j)->status() ==1  ){
//                std::cout << "Muon info coming down from the top!"<< std::endl;
                upsi_mu_daughter_counter +=1;
                std::cout << gen_particle->daughter(j)->status() << std::endl;
                truth_Upsimuon_pt .push_back(gen_particle->daughter(j)->pt());
                truth_Upsimuon_eta.push_back(gen_particle->daughter(j)->eta());
                truth_Upsimuon_phi.push_back(gen_particle->daughter(j)->phi());
                
               //  std::cout << "gen_particle->daughter(j)->pt() " << gen_particle->daughter(j)->pt() << std::endl;
//                 std::cout << "gen_particle->daughter(j)->eta() "<< gen_particle->daughter(j)->eta() << std::endl;
//                 std::cout << "gen_particle->daughter(j)->phi() "<< gen_particle->daughter(j)->phi() << std::endl;
//                 std::cout << "gen_particle->daughter(j)->mass() "<< gen_particle->daughter(j)->mass() << std::endl;
//                 std::cout << "gen_particle->daughter(j)->pdgId() " << gen_particle->daughter(j)->pdgId() << std::endl;
                
                if (upsi_mu_daughter_counter == 1){
                  std::cout << "Filling upsi 1 truth info" << std::endl; 
                  truth_Upsi_pt  .push_back(gen_particle->pt());
                  truth_Upsi_eta .push_back(gen_particle->eta());
                  truth_Upsi_phi .push_back(gen_particle->phi());
                  truth_Upsi_mass.push_back(gen_particle->mass());
                  truth_Upsi_pdgid.push_back(gen_particle->pdgId());
                  
                //   std::cout << "gen_particle->pt()  "<< gen_particle->pt() << std::endl;
//                   std::cout << "gen_particle->eta()  "<< gen_particle->eta() << std::endl;
//                   std::cout << "gen_particle->phi()  "<< gen_particle->phi() << std::endl;
//                   std::cout << "gen_particle->mass()  "<< gen_particle->mass() << std::endl;
//                   std::cout << "gen_particle->pdgId() " << gen_particle->pdgId() << std::endl; 
                  
                  
                
                
                }
                
                if (upsi_mu_daughter_counter == 2) {
                  eventHasUpsiNToMuMu = true; 
                }
             }
         }
        }
        
        
        //ChibJ_1P states that decay with a soft photon to upsi 1's
        
         //Chib0_1P state
         if (gen_particle->pdgId() == Chib0_1P){
           std::cout << "Found Chib0_1P" << std::endl;
           for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
             if (gen_particle->daughter(j)->pdgId() == UPSI){
                std::cout << "Found upsi 1 from Chib0_1P" << std::endl;
                int Chib0_1P_granddaughter_counter = 0;
                for (size_t k =0; k < gen_particle->daughter(j)->numberOfDaughters(); k++){
                 // std::cout << "ANOTHER PLACEHOLDER" << std::endl;
                   
                  if (TMath::Abs(gen_particle->daughter(j)->daughter(k)->pdgId()) == MUON){
                    Chib0_1P_granddaughter_counter++;
                    std::cout << "FILL muons that are granddaughters of Chib0_1P here" << std::endl;
                    truth_Chib0_1P_UPSI_muon_pt.push_back(gen_particle->daughter(j)->daughter(k)->pt());
                    truth_Chib0_1P_UPSI_muon_eta.push_back(gen_particle->daughter(j)->daughter(k)->eta());
                    truth_Chib0_1P_UPSI_muon_phi.push_back(gen_particle->daughter(j)->daughter(k)->phi());
                    
                    if (Chib0_1P_granddaughter_counter == 1){
                       std::cout << "Fill Chib0_1P, Upsi info here" << std::endl; 
                       truth_Chib0_1P_pt.push_back(gen_particle->pt());
                       truth_Chib0_1P_eta.push_back(gen_particle->eta());
                       truth_Chib0_1P_phi.push_back(gen_particle->phi());
                       truth_Chib0_1P_pdgid.push_back(gen_particle->pdgId());
                       truth_Chib0_1P_mass.push_back(gen_particle->mass());
                       
                       truth_Chib0_1P_UPSI_pt.push_back(gen_particle->daughter(j)->pt());
                       truth_Chib0_1P_UPSI_eta.push_back(gen_particle->daughter(j)->eta());
                       truth_Chib0_1P_UPSI_phi.push_back(gen_particle->daughter(j)->phi());
                       truth_Chib0_1P_UPSI_pdgid.push_back(gen_particle->daughter(j)->pdgId());
                       truth_Chib0_1P_UPSI_mass.push_back(gen_particle->daughter(j)->mass());
                    }
                  
                  }
                }
             }
             if (gen_particle->daughter(j)->pdgId() == photon){
                std::cout << "Found photon from Chib0_1P" << std::endl; 
                truth_Chib0_1P_photon_pt.push_back(gen_particle->daughter(j)->pt());
                truth_Chib0_1P_photon_eta.push_back(gen_particle->daughter(j)->eta());
                truth_Chib0_1P_photon_phi.push_back(gen_particle->daughter(j)->phi());
                truth_Chib0_1P_photon_pdgid.push_back(gen_particle->daughter(j)->pdgId());
             }
           
           }
         
         }
         
         //Chib1_1P state
         if (gen_particle->pdgId() == Chib1_1P){
           std::cout << "Found Chib1_1P" << std::endl;
           for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
             if (gen_particle->daughter(j)->pdgId() == UPSI){
                std::cout << "Found upsi 1 from Chib1_1P" << std::endl;
                int Chib1_1P_granddaughter_counter = 0;
                for (size_t k =0; k < gen_particle->daughter(j)->numberOfDaughters(); k++){
                 // std::cout << "ANOTHER PLACEHOLDER" << std::endl;
                   
                  if (TMath::Abs(gen_particle->daughter(j)->daughter(k)->pdgId()) == MUON){
                    Chib1_1P_granddaughter_counter++;
                    std::cout << "FILL muons that are granddaughters of Chib1_1P here" << std::endl;
                    truth_Chib1_1P_UPSI_muon_pt.push_back(gen_particle->daughter(j)->daughter(k)->pt());
                    truth_Chib1_1P_UPSI_muon_eta.push_back(gen_particle->daughter(j)->daughter(k)->eta());
                    truth_Chib1_1P_UPSI_muon_phi.push_back(gen_particle->daughter(j)->daughter(k)->phi());
                    
                    if (Chib1_1P_granddaughter_counter == 1){
                       std::cout << "Fill Chib1_1P, Upsi info here" << std::endl; 
                       truth_Chib1_1P_pt.push_back(gen_particle->pt());
                       truth_Chib1_1P_eta.push_back(gen_particle->eta());
                       truth_Chib1_1P_phi.push_back(gen_particle->phi());
                       truth_Chib1_1P_pdgid.push_back(gen_particle->pdgId());
                       truth_Chib1_1P_mass.push_back(gen_particle->mass());
                       
                       truth_Chib1_1P_UPSI_pt.push_back(gen_particle->daughter(j)->pt());
                       truth_Chib1_1P_UPSI_eta.push_back(gen_particle->daughter(j)->eta());
                       truth_Chib1_1P_UPSI_phi.push_back(gen_particle->daughter(j)->phi());
                       truth_Chib1_1P_UPSI_pdgid.push_back(gen_particle->daughter(j)->pdgId());
                       truth_Chib1_1P_UPSI_mass.push_back(gen_particle->daughter(j)->mass());
                    }
                  
                  }
                }
             }
             if (gen_particle->daughter(j)->pdgId() == photon){
                std::cout << "Found photon from Chib1_1P" << std::endl; 
                truth_Chib1_1P_photon_pt.push_back(gen_particle->daughter(j)->pt());
                truth_Chib1_1P_photon_eta.push_back(gen_particle->daughter(j)->eta());
                truth_Chib1_1P_photon_phi.push_back(gen_particle->daughter(j)->phi());
                truth_Chib1_1P_photon_pdgid.push_back(gen_particle->daughter(j)->pdgId());
             }
           
           }
         
         }
         
         //Chib2_1P state
         if (gen_particle->pdgId() == Chib2_1P){
           std::cout << "Found Chib2_1P" << std::endl;
           for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
             if (gen_particle->daughter(j)->pdgId() == UPSI){
                std::cout << "Found upsi 1 from Chib2_1P" << std::endl;
                int Chib2_1P_granddaughter_counter = 0;
                for (size_t k =0; k < gen_particle->daughter(j)->numberOfDaughters(); k++){
                 // std::cout << "ANOTHER PLACEHOLDER" << std::endl;
                
                 //Comments to  better understand what is happening in this sample
                 if (TMath::Abs(gen_particle->daughter(j)->daughter(k)->pdgId()) != MUON){
                    std::cout << "Upsi 1 daughter is not a muon" << std::endl;
                    std::cout << TMath::Abs(gen_particle->daughter(j)->daughter(k)->pdgId()) << std::endl;
                 }
                 
                 
                  if (TMath::Abs(gen_particle->daughter(j)->daughter(k)->pdgId()) == MUON){
                    Chib2_1P_granddaughter_counter++;
                    std::cout << "FILL muons that are granddaughters of Chib2_1P here" << std::endl;
                    truth_Chib2_1P_UPSI_muon_pt.push_back(gen_particle->daughter(j)->daughter(k)->pt());
                    truth_Chib2_1P_UPSI_muon_eta.push_back(gen_particle->daughter(j)->daughter(k)->eta());
                    truth_Chib2_1P_UPSI_muon_phi.push_back(gen_particle->daughter(j)->daughter(k)->phi());
                    
                    if (Chib2_1P_granddaughter_counter == 1){
                       std::cout << "Fill Chib2_1P, Upsi info here" << std::endl; 
                       truth_Chib2_1P_pt.push_back(gen_particle->pt());
                       truth_Chib2_1P_eta.push_back(gen_particle->eta());
                       truth_Chib2_1P_phi.push_back(gen_particle->phi());
                       truth_Chib2_1P_pdgid.push_back(gen_particle->pdgId());
                       truth_Chib2_1P_mass.push_back(gen_particle->mass());
                       
                       truth_Chib2_1P_UPSI_pt.push_back(gen_particle->daughter(j)->pt());
                       truth_Chib2_1P_UPSI_eta.push_back(gen_particle->daughter(j)->eta());
                       truth_Chib2_1P_UPSI_phi.push_back(gen_particle->daughter(j)->phi());
                       truth_Chib2_1P_UPSI_pdgid.push_back(gen_particle->daughter(j)->pdgId());
                       truth_Chib2_1P_UPSI_mass.push_back(gen_particle->daughter(j)->mass());
                    }
                  
                  }
                }
             }
             if (gen_particle->daughter(j)->pdgId() == photon){
                std::cout << "Found photon from Chib2_1P" << std::endl; 
                truth_Chib2_1P_photon_pt.push_back(gen_particle->daughter(j)->pt());
                truth_Chib2_1P_photon_eta.push_back(gen_particle->daughter(j)->eta());
                truth_Chib2_1P_photon_phi.push_back(gen_particle->daughter(j)->phi());
                truth_Chib2_1P_photon_pdgid.push_back(gen_particle->daughter(j)->pdgId());
             }
           
           }
         
         }
         
         /////////////////////////////////////
         
         //UPSI_2S state 
         if (gen_particle->pdgId() == UPSI_2S){
            std::cout << "gen_particle->pdgId() " << gen_particle->pdgId() << std::endl;
            std::cout << "gen_particle->status() " << gen_particle->status() << std::endl;
            
            int upsi2_mu_daughter_counter = 0;
            for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
              if (TMath::Abs(gen_particle->daughter(j)->pdgId()) == MUON){ //&& gen_particle->daughter(j)->status() ==1  ){
              //  std::cout << "PLACEHOLDER UPSI 2!"<< std::endl;
                upsi2_mu_daughter_counter +=1;
                std::cout << gen_particle->daughter(j)->status() << std::endl;
                truth_Upsi2muon_pt .push_back(gen_particle->daughter(j)->pt());
                truth_Upsi2muon_eta.push_back(gen_particle->daughter(j)->eta());
                truth_Upsi2muon_phi.push_back(gen_particle->daughter(j)->phi());
                
                if (upsi2_mu_daughter_counter == 1){
                  std::cout << "Filling upsi 2 truth info" << std::endl; 
                  truth_Upsi2_pt  .push_back(gen_particle->pt());
                  truth_Upsi2_eta .push_back(gen_particle->eta());
                  truth_Upsi2_phi .push_back(gen_particle->phi());
                  truth_Upsi2_mass.push_back(gen_particle->mass());
                  truth_Upsi2_pdgid.push_back(gen_particle->pdgId());
                
                
                }
                
                 if (upsi2_mu_daughter_counter == 2) {
                   eventHasUpsiNToMuMu = true;
                 }
             }
         }
        }
        
        //UPSI_3S state
         if (gen_particle->pdgId() == UPSI_3S){
            std::cout << "gen_particle->pdgId() " << gen_particle->pdgId() << std::endl;
            std::cout << "gen_particle->status() " << gen_particle->status() << std::endl;
            
            int upsi3_mu_daughter_counter = 0;
            for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
              if (TMath::Abs(gen_particle->daughter(j)->pdgId()) == MUON){ //&& gen_particle->daughter(j)->status() ==1  ){
              //  std::cout << "PLACEHOLDER UPSI 3!"<< std::endl;
                upsi3_mu_daughter_counter +=1;
                std::cout << gen_particle->daughter(j)->status() << std::endl;
                truth_Upsi3muon_pt .push_back(gen_particle->daughter(j)->pt());
                truth_Upsi3muon_eta.push_back(gen_particle->daughter(j)->eta());
                truth_Upsi3muon_phi.push_back(gen_particle->daughter(j)->phi());
                
                if (upsi3_mu_daughter_counter == 1){
                  std::cout << "Filling upsi 3 truth info" << std::endl; 
                  truth_Upsi3_pt  .push_back(gen_particle->pt());
                  truth_Upsi3_eta .push_back(gen_particle->eta());
                  truth_Upsi3_phi .push_back(gen_particle->phi());
                  truth_Upsi3_mass.push_back(gen_particle->mass());
                  truth_Upsi3_pdgid.push_back(gen_particle->pdgId());
                
                
                }
                 if (upsi3_mu_daughter_counter == 2) {
                   eventHasUpsiNToMuMu = true; 
                 }
             }
         }
        }
        
  //    registered_upsi = false;
//      POODLE=false;
 //      if (gen_particle->pdgId() == UPSI) { 
//         for (int iupsi = 0; iupsi < (int)upsi_pt_found.size(); iupsi++) {
//           if (TMath::Abs( gen_particle->pt() - upsi_pt_found.at(iupsi) ) < 0.1)
//             registered_upsi = true;
// //          std::cout << "diff " << (TMath::Abs( gen_particle->pt() - upsi_pt_found.at(iupsi) )) << std::endl; //mary added for testing/understanding this loop
//         }
//         for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
//           if (TMath::Abs(gen_particle->daughter(j)->pdgId()) == MUON && !registered_upsi) {
//             truth_Upsimuon_pt .push_back(gen_particle->daughter(j)->pt());
//             truth_Upsimuon_eta.push_back(gen_particle->daughter(j)->eta());
//             truth_Upsimuon_phi.push_back(gen_particle->daughter(j)->phi());
//             
//             //WARNING THE UPSI INFO NOW BEING FILLED FOR ALL MUONS, NEED TO CLEAN THIS UP WITH A MUON COUNTER TO DO 
//             truth_Upsi_pt  .push_back(gen_particle->pt());
//             truth_Upsi_eta .push_back(gen_particle->eta());
//             truth_Upsi_phi .push_back(gen_particle->phi());
//             truth_Upsi_mass.push_back(gen_particle->mass());
//             truth_Upsi_pdgid.push_back(gen_particle->pdgId());
//             Upsi_to_fill_count+=1;
//           }
//         }
//         upsi_pt_found  .push_back(gen_particle->pt());
// //        std::cout << "upsi_pt_found.size() now: " << upsi_pt_found.size() << std::endl; //mary added for testing 
// //        for (unsigned int i = 0; i < upsi_pt_found.size(); i++) {
//  //       std::cout<< "upsi_pt_found.at(i): " << upsi_pt_found.at(i) << std::endl; 
//  //       }
//       }
     
       
   //   registered_Z = false;
      if (gen_particle->pdgId() == Z) {
       
       //  std::cout << gen_particle->pt() << std::endl;
//         std::cout<< gen_particle->numberOfDaughters() << std::endl;
//         for (int iZ = 0; iZ < (int)Z_pt_found.size(); iZ++) {
//           if (TMath::Abs( gen_particle->pt() - Z_pt_found.at(iZ) ) < .05){
//   //          POODLE = true;
//             std::cout << "POODLE" << std::endl;
//           }
//           if (TMath::Abs( gen_particle->pt() - Z_pt_found.at(iZ) ) < 0.0000001) //tune this number! 0.000001 gives all events back
//             
//             registered_Z = true;
//  //           registered_Z_count+=1;
//           
//         //   if (!registered_Z){
// //             truth_Z_pt  .push_back(gen_particle->pt());
// //             truth_Z_eta .push_back(gen_particle->eta());
// //             truth_Z_phi .push_back(gen_particle->phi());
// //             truth_Z_mass.push_back(gen_particle->mass());
// //             truth_Z_pdgid.push_back(gen_particle->pdgId());
// //           
// //           }
// 
// 
//          }
     //    if (!registered_Z) {
//           truth_Z_pt  .push_back(gen_particle->pt());
//           truth_Z_eta .push_back(gen_particle->eta());
//            truth_Z_phi .push_back(gen_particle->phi());
//            truth_Z_mass.push_back(gen_particle->mass());
//            truth_Z_pdgid.push_back(gen_particle->pdgId());
//         
//         }
        int Z_mu_daughter_counter = 0;
        for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
          if (TMath::Abs(gen_particle->daughter(j)->pdgId()) == MUON ) {
            Z_mu_daughter_counter++;
            std::cout << "PLACEHOLDER 2" << std::endl; 
            std::cout << gen_particle->daughter(j)->status() << std::endl;
            truth_Zmuon_pt .push_back(gen_particle->daughter(j)->pt());
            truth_Zmuon_eta.push_back(gen_particle->daughter(j)->eta());
            truth_Zmuon_phi.push_back(gen_particle->daughter(j)->phi());
            Z_to_fill_count+=1;
            
           if (Z_mu_daughter_counter == 1){
             truth_Z_pt  .push_back(gen_particle->pt());
             truth_Z_eta .push_back(gen_particle->eta());
             truth_Z_phi .push_back(gen_particle->phi());
             truth_Z_mass.push_back(gen_particle->mass());
             truth_Z_pdgid.push_back(gen_particle->pdgId());
           }
           if (Z_mu_daughter_counter == 2){
             eventHasZToMuMu = true;
           }
          }
        }
 //       Z_pt_found  .push_back(gen_particle->pt());
      }
      
      
      }
    
    if (eventHasZToMuMu && eventHasUpsiNToMuMu){
        std::cout << "eventHasZToMuMu && eventHasUpsiNToMuMu" << std::endl; 
        eventHasZUpsiNTo4Mu = true;  
        eventHasZUpsiNTo4Mu_Count.push_back(eventHasZUpsiNTo4Mu);

    }
    treemc->Fill();
 //   std::cout << "Z_to_fill_count is:" << Z_to_fill_count << std::endl;
 //   std::cout << "Upsi_to_fill_count is:" << Upsi_to_fill_count << std::endl;
   } // end of mc

}


//If applyZmuonFilter=False, comment out this whole block below from void to last bracket above the //define this as a plug-in comment 
void
ZmuonAnalyzer::endLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup)
{
//  events seen in each lumi section
     edm::Handle<edm::MergeableCounter> nTotalHandle;
      iLumi.getByToken(nTotalToken_, nTotalHandle);
     int nTotInThisLS = nTotalHandle->value;
      
      histContainer_["nEventsTotalCounter"]->Fill(0); //the number of entries at 0 give us the number of LS total that were looked at at the end of the day
      histContainer_["nEventsTotalCounter"]->Fill(1, nTotInThisLS); //the number of entries at 1 will give us the total number of events looked at (each lumi section is weighted by the number of events seen in that lumi section)
      
      //Events seen in each LS that passed the filter
      edm::Handle<edm::MergeableCounter> nFilteredHandle;
     iLumi.getByToken(nFilteredToken_, nFilteredHandle);
      int nFilteredInThisLS = nFilteredHandle->value;
      
      histContainer_["nEventsFilteredCounter"]->Fill(0); //the number of entries at 0 give us the number of LS in total that passed the filter
      histContainer_["nEventsFilteredCounter"]->Fill(1, nFilteredInThisLS); //the number of entries at 1 will give us the total number of events that passed the filter 
    
    // Total number of events is the sum of the events in each of these luminosity blocks 
//    edm::Handle nEventsTotalCounter;
   // iLumi.getByLabel("nEventsTotal", nEventsTotalCounter); nEventsTotal +=nEventsTotalCounter->value;
//
//    edm::Handle nEventsFilteredCounter; iLumi.getByLabel("nEventsFiltered",
  //  nEventsFilteredCounter); nEventsFiltered +=nEventsFilteredCounter->value; 

}



//define this as a plug-in
DEFINE_FWK_MODULE(ZmuonAnalyzer);

