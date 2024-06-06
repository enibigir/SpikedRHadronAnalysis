# HSCP
This is a fork of the SUSYBSMAnalysis repository, a large repo that contains code for many different purposes. The focus of this fork will be in the HSCP directory.

## Setup
Once you are in the cmslpc cluster, run the following.
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_30
cd CMSSW_10_6_30/src
cmsenv
```
For the following step you should have a ssh key associated to your GitHub account. For more information, see [connecting-to-github-with-ssh-key.](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
```
git clone -b master git@github.com:ctdax/SUSYBSMAnalysis.git SUSYBSMAnalysis
```

## SimCaloHitAnalyzer.cc & GenSimEDMAnalyzer.cc

GenSimEDMAnalyzer.cc analyzes AOD level root files containing HSCP candidates in the tracker and muon chambers, it was written by Todd. SimCaloHitAnalyzer.cc is my revision of GenSimEDMAnalyzer to observe R-Hadron hits and saturation in the calorimiter.

## Information I found helpful when learning CMSSW

### CMSSW Directory Format

If this is your first time using or seeing a CMSSW directory I highly recommend you checkout [this brief presentation](https://www.hep.ph.ic.ac.uk/~dbauer/cms/tutorial_2011.pdf) and also [this twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter4#SampleCode), particularly section 4.1.2 'Writing your own EDAnalyzer' of the twiki. Ultimately you should be able to run the analyzer smoothly inside of either directory via the command
```
cmsRun python/config_cfg.py
```

### AOD Contents

The SimCaloHitAnalyzer currently analyzes an AOD gluino root file which is present in config_cfg.py. The contents of the AOD file can either be observed in a TBrowser
```
root -l file:/uscms_data/d2/tadams/hscp/fall22a/CMSSW_10_6_30/src/EXO-RunIISummer20UL18GENSIM-00010-v3.root
TBrowser b("Events")
```
or by using edmDump
```
edmDumpEventContent file:/uscms_data/d2/tadams/hscp/fall22a/CMSSW_10_6_30/src/EXO-RunIISummer20UL18GENSIM-00010-v3.root
```
