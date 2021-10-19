import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing #https://twiki.cern.ch/CMSPublic/SWGuideAboutPythonConfigFile#Passing_Command_Line_Arguments_T

process = cms.Process("analysis")

options = VarParsing.VarParsing("analysis")



#List of Options
options.register( "applyZmuonFilter",
    True,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "Shall we prefilter out some events and never bother giving them to the analyzer? Defaults to True. If you change to false, you will need to, for the moment, comment out the indicated lines in the ZmuonAnalyzer.cc. Will try to make a nice boolean to avoid this later. Will flag with a comment 'If applyZmuonFilter=False, comment out!'"
)

options.register("isMC",
    True, #set to false for crab test run
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "Is this simulation (aka MC) or is it data we have collected from the CMS detector? Defaults to True."
            
)

# This needs to go here, after we have defined all the options!
options.parseArguments()
print "Options are:", options

process.ZmuonAnalyzer = cms.EDAnalyzer("ZmuonAnalyzer",
   muonCollection = cms.InputTag("slimmedMuons"),
   electronCollection = cms.InputTag("slimmedElectrons"),
   bits = cms.InputTag("TriggerResults","", "HLT"),
   #objects = cms.InputTag("selectedPatTrigger"),
   objects = cms.InputTag("slimmedPatTrigger"), #corrected 11 Oct. 2021, see email from K.P.
   genParticles = cms.InputTag("prunedGenParticles"),
   vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
   metTag = cms.InputTag("slimmedMETs"),
   pfCands = cms.InputTag("packedPFCandidates"),
   rho     = cms.InputTag("fixedGridRhoFastjetAll"), #Median energy density, see: https://arxiv.org/pdf/0707.1378.pdf
  #Might add Pileup here later if I think it is something we might change ever, but at the moment is hard-coded in the analyzer
   isMC   = cms.bool(options.isMC),
)

if options.applyZmuonFilter:
   process.ZmuonFilter = cms.EDFilter("ZmuonFilter",
     muonCollection = cms.InputTag("slimmedMuons"),
     bits = cms.InputTag("TriggerResults","", "HLT"),
 #    objects = cms.InputTag("selectedPatTrigger"),
     objects = cms.InputTag("slimmedPatTrigger"),
#     genParticles = cms.InputTag("prunedGenParticles"), #not needed here
#     pfCands = cms.InputTag("packedPFCandidates")  #not needed here
     pTCut = cms.double(2.),
     etaCut = cms.double(3.),
     invMass4MuCut_low = cms.double(50), #usually 50
    #invMass4MuCut_high = cms.double(500.), #this is essentially not doing anything
  )
   process.nEventsTotal = cms.EDProducer("EventCountProducer")
   process.nEventsFiltered = cms.EDProducer("EventCountProducer")

#process.maxEvents.input = 1000

#process.maxEvents = cms.untracked.PSet(
#   input = cms.untracked.int32(40000)
#   input = cms.untracked.int32(-1)
 #  input = cms.untracked.int32(1000) #for crab test, just look at 1000 events

#)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#                                          "file:../miniAOD_01.root",
#                                          "file:../miniAOD_02.root",
#                                          "file:../miniAOD_03.root",
#                                          "file:../miniAOD_04.root",
#                                          "file:../miniAOD_05.root",
#                                          "file:../miniAOD_06.root"
#                                          "file:../../../miniAOD_lowStats.root",
#                                           "file:../SingleMu_Run2018A-17Sep2018-v2_F8CDAAA9-11C7-A34F-A059-409CF95EB82A.root",
#                                           "file:../SingleMu_2017B-31Mar2018-v1_FC2B7874-F538-E811-9C29-0025905A60A8.root",
#                                            "file:0079D4A1-71DE-AF4B-90E4-115037F02923.root", #DiMu Run18 A
 #                                             "file:FD881AE8-B6AC-8E4F-9453-2D5E6135D476.root" #DiMu Run18 A
 #                                            "file:E36DDD08-CC3D-6D49-8C71-C6E8C201A2E5.root", #this is DiMU 2018RunB data
 #                                            "file:../YZ_4MuFS_DPS_APV_MiniAOD_2016UL_Yfirst_Zsecond_03_10_2020.root",
 #                                           "file:../YZ_4MuFS_DPS_APV_MiniAOD_2016UL_Zfirst_Ysecond_03_10_2020.root",
#                                           "file:../YZ_4MuFS_DPS_MiniAOD_2016UL_Yfirst_Zsecond_03_10_2020.root",
#                                          "file:../YZ_4MuFS_DPS_MiniAOD_2016UL_Zfirst_Ysecond_03_10_2020.root",
#                                           "file:../YZ_4MuFS_DPS_MiniAOD_2016UL_Yfirst_Zsecond_06_11_2020.root",
#                                            "file:020353D5-EB7E-3A42-928B-64ABB6449999.root",
#                                            "file:../miniAOD_Marys_LHE.root",
                                            #  "file:Y1S3S11_Z_4Mu_HO_SPS_MiniAOD_shortTest_11_June_2021.root",
   #                                         "file:Y2S3S11_Z_4Mu_HO_SPS_MiniAOD.root",
                                           "file:MC_DPS_2016_YZ_00623223-2B20-AB42-A456-670F9B3875D5.root",
                                           #'file:14907725-2B14-B245-B076-5B04C5C36D55.root',
      #                                     'file:Y3S3S11_Z_4Mu_HO_SPS_MiniAOD.root',

                                    ),
   duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.maxEvents = cms.untracked.PSet(
#   input = cms.untracked.int32(40000)
  #input = cms.untracked.int32(5)
# input = cms.untracked.int32(1)
 input=cms.untracked.int32(-1)
  #  input = cms.untracked.int32(3000) #for crab test, just look at 1000 events

)

process.TFileService = cms.Service("TFileService",
   fileName = cms.string("testNewMCSection.root")
)

#process.maxEvents.input = 1000
#process.maxEvents = cms.untracked.PSet(
#   input = cms.untracked.int32(40000)
#   input = cms.untracked.int32(-1)
#    input = cms.untracked.int32(1000) #for crab test, just look at 1000 events

#)

process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
#   SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1#10000

#process.MessageLogger = cms.Service("MessageLogger",
#                    destinations   = cms.untracked.vstring('messages.txt')
#)

if options.applyZmuonFilter:
    process.p = cms.Path(
         process.nEventsTotal *process.ZmuonFilter *process.nEventsFiltered * process.ZmuonAnalyzer 
        )

else:


    process.p = cms.Path(
        process.ZmuonAnalyzer
)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

