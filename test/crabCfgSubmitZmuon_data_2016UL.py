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
    
   
  # 2017 Common Configuration
    config.General.workArea      = 'Zmuon_DataJobs_UL2016BCDEFGH_DoubleMu_7Oct2021_10LSPerJob'
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
    config.Data.lumiMask        = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt:" #35.93 fb^-1, see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/DCUserPage#Legacy_2016_aka_UL2016_datasets

    config.Data.unitsPerJob     = 10 #5
#    config.Data.unitsPerJob     = 6
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
    
    
    

    # Run2016B DoubleMu HIPM ver1
    config.General.requestName = 'DoubleMuUL_Run2016B_HIPM_ver1_7Oct2021_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD' 
    config.Data.outLFNDirBase  = '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2016B_HIPM_ver1_7Oct2021_10LSPerJob'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
  # Run2016B DoubleMu HIPM ver2
    config.General.requestName = 'DoubleMuUL_Run2016B_HIPM_ver2_7Oct2021_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD'
    config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2016B_HIPM_ver2_7Oct2021_10LSPerJob'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

  # Run2016C DoubleMu HIPM
    config.General.requestName = 'DoubleMuUL_Run2016C_HIPM_7Oct2021_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD' 
    config.Data.outLFNDirBase  = '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2016C_HIPM_7Oct2021_10LSPerJob' 
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
#   Run2016D DoubleMu  HIPM  
    config.General.requestName = 'DoubleMuUL_Run2016D_HIPM_7Oct2021_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD'
    config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2016D_HIPM_7Oct2021_10LSPerJob' 
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
  #  Run2016E DoubleMu HIPM    
    config.General.requestName = 'DoubleMuUL_Run2016E_HIPM_7Oct2021_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD'
    config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2016E_HIPM_7Oct2021_10LSPerJob' 
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

#Start here next time, you are about to do 2016F

  #  Run2016F HIPM DoubleMu    
    config.General.requestName = 'DoubleMuUL_Run2016F_HIPM_7Oct2021_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD'
    config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2016F_HIPM_7Oct2021_10LSPerJob' 
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    #  Run2016F DoubleMu    
    config.General.requestName = 'DoubleMuUL_Run2016F_7Oct2021_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2016F-21Feb2020_UL2016-v1/MINIAOD'
    config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2016F_7Oct2021_10LSPerJob' 
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
     #  Run2016G DoubleMu    
    config.General.requestName = 'DoubleMuUL_Run2016G_7Oct2021_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2016G-21Feb2020_UL2016-v1/MINIAOD'
    config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2016G_7Oct2021_10LSPerJob' 
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    
     #  Run2016H DoubleMu    
    config.General.requestName = 'DoubleMuUL_Run2016H_7Oct2021_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2016H-21Feb2020_UL2016-v1/MINIAOD'
    config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2016H_7Oct2021_10LSPerJob' 
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    



    

