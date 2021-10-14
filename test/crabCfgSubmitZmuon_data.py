#Usage: 
#python crabCfgSubmitZmuon_data.py (to submit jobs)
# ./multicrab.py (to check jobs) -w <workAreaName> -c <command to launch, e.g. status>
#Credit to Maximilian Heindl and this github repo for idea of how to write this crab cfg: https://github.com/CMSHCALCalib/RadDam/blob/master/HFmonitoring/nTuplizer/ggAnalysis/ggNtuplizer/test/crabConfig_data.py


#also credit to this tool for the multicrab.py file: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRABClientLibraryAPI#Multicrab_using_the_crabCommand

if __name__ == '__main__':

    


    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    from CRABClient.UserUtilities import config
    config = config()
    
    from multiprocessing import Process
    
    
    #Common configuration
    
    config.General.workArea     = 'Zmuon_DataJobs_UL2018ABCD_DoubleMu_22Mar2021_removePhase1SoftMuCut_5LSPerJob'
    config.General.transferLogs = False
#    config.JobType.maxMemoryMB = 5000 #Let's try the default to start and see if it works 
#    config.JobType.maxJobRuntimeMin = 2750 #Let's try the default to start and see if it works 
    config.JobType.pluginName   = 'Analysis' 
    config.JobType.psetName     = 'ZmuonAnalyzer_cfg.py'
    config.JobType.pyCfgParams = ["isMC=False"]
    config.JobType.allowUndistributedCMSSW = True
#    config.JobType.sendExternalFolder = True #I don't have an CMSSW_BASE/external so I don't think I need this 
    config.Data.inputDBS        = 'global'    #Checked, this is what we need 
    config.Data.splitting       = 'LumiBased' 
    config.Data.lumiMask        = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt' #roughly 58 fb^-1, see: https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis#13_TeV_pp_runs_ReReco (info only given for PromptReco)
    config.Data.unitsPerJob     = 7 #Try now with just 5 per job 
#    config.Data.totalUnits      = 100 #for testing purposes 
    config.Data.ignoreLocality  = False
    config.Data.publication     = False
    config.Data.allowNonValidInputDataset     = True
    config.Site.storageSite     = 'T3_US_FNALLPC' 
    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers
            
    
    
    #Dataset dependent Configuration
    
    
    
    # Run2018A SingleMu  
#     config.General.requestName = "SingleMuUL_Run2018A"
#     config.Data.inputDataset = "/SingleMuon/Run2018D-12Nov2019_UL2018-v4/MINIAOD"
#     config.Data.outLFNDirBase = "/store/user/mhadley/Zmuon_DataJobs_SingleMu_UL2018A_proofOfConcept_12Oct2020"
#     p = Process(target=submit, args=(config,))
#     p.start()
#     p.join()

    # Run2018A DoubleMu
    config.General.requestName = 'DoubleMuUL_Run2018A_22Mar2021_removePhase1SoftMuCut_7LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2018A-12Nov2019_UL2018-v2/MINIAOD' 
    config.Data.outLFNDirBase  = '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2018A_22Mar2021_removePhase1SoftMuCut_7LSPerJob'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
#     Run2018B DoubleMu
#     config.General.requestName = 'DoubleMuUL_Run2018B_22Mar2021_removePhase1SoftMuCut_5LSPerJob' 
#     config.Data.inputDataset   = '/DoubleMuon/Run2018B-12Nov2019_UL2018-v2/MINIAOD'
#     config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2018B_22Mar2021_removePhase1SoftMuCut_5LSPerJob'
#     p = Process(target=submit, args=(config,))
#     p.start()
#     p.join()
# 
#     Run2018C DoubleMu
#     config.General.requestName = 'DoubleMuUL_Run2018C_22Mar2021_removePhase1SoftMuCut_5LSPerJob' 
#     config.Data.inputDataset   = '/DoubleMuon/Run2018C-12Nov2019_UL2018-v2/MINIAOD' 
#     config.Data.outLFNDirBase  = '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2018C_22Mar2021_removePhase1SoftMuCut_5LSPerJob' 
#     p = Process(target=submit, args=(config,))
#     p.start()
#     p.join()
#     
#     Run2018D DoubleMu    
#     config.General.requestName = 'DoubleMuUL_Run2018D_22Mar2021_removePhase1SoftMuCut_5LSPerJob' 
#     config.Data.inputDataset   = '/DoubleMuon/Run2018C-12Nov2019_UL2018-v2/MINIAOD'
#     config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2018D_22Mar2021_removePhase1SoftMuCut_5LSPerJob' 
#     p = Process(target=submit, args=(config,))
#     p.start()
#     p.join()

    

