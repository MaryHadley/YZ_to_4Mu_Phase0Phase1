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
    config.General.workArea      = 'Zmuon_DataJobs_UL2017BC_DoubleMu_4May2022_10LSPerJob'
    config.General.transferLogs = False
#    config.JobType.maxMemoryMB = 5000 #Let's try the default to start and see if it works 
#    config.JobType.maxJobRuntimeMin = 2750 #Let's try the default to start and see if it works 
    config.JobType.pluginName   = 'Analysis' 
    config.JobType.psetName     = 'ZmuonAnalyzer_cfg.py'
    config.JobType.pyCfgParams = ["isMC=False", "triggerYear=2017"]
    config.JobType.allowUndistributedCMSSW = True
#    config.JobType.sendExternalFolder = True #I don't have an CMSSW_BASE/external so I don't think I need this 
    config.Data.inputDBS        = 'global'    #Checked, this is what we need 
    config.Data.splitting       = 'LumiBased' 
    config.Data.lumiMask        =  'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt' #Use https://blah blah blah when submitting from LPC, use /afs/blah blah blah when submitting from lxplus #'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt' #41.48 fb^-1, see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/DCUserPage#Legacy_2017_aka_UL2017_datasets

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
    
    
    

    # Run2017B DoubleMu
    config.General.requestName = 'DoubleMuUL_Run2017B_4May2022_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2017B-09Aug2019_UL2017-v1/MINIAOD' 
    config.Data.outLFNDirBase  = '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2017B_4May2022_10LSPerJob'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
  # Run2017C DoubleMu
    config.General.requestName = 'DoubleMuUL_Run2017C_4May2022_10LSPerJob' 
    config.Data.inputDataset   = '/DoubleMuon/Run2017C-09Aug2019_UL2017-v1/MINIAOD'
    config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2017C_4May2022_10LSPerJob'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

 #  # Run2017D DoubleMu
#     config.General.requestName = 'DoubleMuUL_Run2017D_3May2022_10LSPerJob' 
#     config.Data.inputDataset   = '/DoubleMuon/Run2017D-09Aug2019_UL2017-v1/MINIAOD' 
#     config.Data.outLFNDirBase  = '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2017D_3May2022_10LSPerJob' 
#     p = Process(target=submit, args=(config,))
#     p.start()
#     p.join()
#     
# #   Run2017E DoubleMu    
#     config.General.requestName = 'DoubleMuUL_Run2017E_3May2022_10LSPerJob' 
#     config.Data.inputDataset   = '/DoubleMuon/Run2017E-09Aug2019_UL2017-v1/MINIAOD'
#     config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2017E_3May2022_10LSPerJob' 
#     p = Process(target=submit, args=(config,))
#     p.start()
#     p.join()
#     
#   #  Run2017F DoubleMu    
#     config.General.requestName = 'DoubleMuUL_Run2017F_3May2022_10LSPerJob' 
#     config.Data.inputDataset   = '/DoubleMuon/Run2017F-09Aug2019_UL2017-v1/MINIAOD'
#     config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2017F_3May2022_10LSPerJob' 
#     p = Process(target=submit, args=(config,))
#     p.start()
#     p.join()


#Was a 5 TeV run, Do NOT use (if you try to submit it to Crab, the lumi mask saves you, as crab finds no events from Era G in the Golden JSON, so you just get a SUBMITFAILED message)
  #  Run2017G DoubleMu    
#     config.General.requestName = 'DoubleMuUL_Run2017G_3May2022_10LSPerJob' 
#     config.Data.inputDataset   = '/DoubleMuon/Run2017G-09Aug2019_UL2017-v1/MINIAOD'
#     config.Data.outLFNDirBase  =  '/store/user/mhadley/Zmuon_DataJobs_DiMu_UL2017G_3May2022_10LSPerJob' 
#     p = Process(target=submit, args=(config,))
#     p.start()
#     p.join()


    

