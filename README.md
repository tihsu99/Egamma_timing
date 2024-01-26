# Egamma timing study
To install the repository
```
cmsrel CMSSW_13_1_0
cd CMSSW_13_1_0/src
git clone https://github.com/tihsu99/Egamma_timing.git
cmsenv
scram b -j 8
cd Egamma_timing/EgammaTimingTools
```
The main code is stored in `plugins`. For `python`, it stored configuration file. To run it locally you first need the `reRECO.root` produced from [Generate_reReco](https://github.com/tihsu99/Generate_reReco):
```
cmsRun testElectronMVA_cfg_mod1.py electronLabel=ecalDrivenGsfElectronsHGC inputFiles=file:reRECO.root
```
And to submit it to condor (please remember to do `voms` first)
```
`mkdir -p $HOME/tmp
export X509_USER_PROXY=$HOME/tmp/x509up
python submit_condor.py --outdir [YOUR/OUTPUT/DIRECTORY]
```
To do the analysis, since `coffea` is used, it may be confict with `cmsenv`. Please restart the terminal and do not run the code under `cmsenv`
```
cd analysis
source script/env.sh
python summarize.py --indir [The/OUTPUT/DIRECTORY in previous stage] --fetch_distribution [To summarize the histograms and stored in root file]
python summarize.py --draw [To plot the distribution and AUC/ROC diagram]
```
