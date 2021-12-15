# YZ_to_4Mu_Phase0Phase1


mkdir < workArea >  
cd < workArea >  
cmsrel CMSSW_10_6_27 #10_6_27 is the correct release to use for UL analyses    
cd CMSSW_10_6_27/src  
  
git clone git@github.com:MaryHadley/YZ_to_4Mu_Phase0Phase1.git  
cd YZ_to_4Mu_Phase0Phase1  
mv * ..  
cd ..  
rm -rf YZ_to_4Mu_Phase0Phase1  
cmsenv  
scram b -j8 USER_CXXFLAGS="-g"  

