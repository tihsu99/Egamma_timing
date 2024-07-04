import subprocess, os
import optparse, argparse
import ROOT

cwd = os.getcwd()
dir_list = cwd.split('/')
cmssw_list = []
for dir_ in dir_list:
  cmssw_list.append(dir_)
  if 'CMSSW' in dir_: break
cmsswBase = '/'.join(cmssw_list)

def prepare_shell(shell_file, command, condor, FarmDir):
  cwd = os.getcwd()
  with open(os.path.join(FarmDir, shell_file), 'w') as shell:
    shell.write('#!/bin/bash\n')
    shell.write('WORKDIR=%s\n'%cwd)
    shell.write('cd %s\n'%cmsswBase)
    shell.write('eval `scram r -sh`\n')
    shell.write('cd ${WORKDIR}\n')
    shell.write(command)
  condor.write('cfgFile=%s\n'%shell_file)
  condor.write('queue 1\n')

def wrapper_condor(dataset_name, dataset_path, condor, outdir, farm_dir, barrel=False, check=False, particle='electron'):
  result = subprocess.Popen("gfal-ls " + dataset_path, shell=True, stdout=subprocess.PIPE)
  subprocess_return = result.stdout.read()
  file_list = subprocess_return.decode('utf-8').split("\n")
  n_file = len(file_list)
  for idx, inf in enumerate(file_list):
    if not ('.root' in inf): continue
    fout = os.path.join(outdir, dataset_name, 'ntupleTree_{}.root'.format(idx))
    #print(idx, inf)
    if check:
      try:
        fout_root = ROOT.TFile.Open(fout, "READ")
        if fout_root.IsZombie():
          print("{} is Zombie.".format(fout))
        h = fout_root.Get("ntuplizer/tree")
        if particle == 'electron':
          a = h.ele_pt
        else:
          a = h.eta
    #    print(h)
        fout_root.Close()
        continue
      except Exception as error:
        print(error)
        print("reproduce {}".format(fout))
    print(idx, inf)
    inf_full_path = os.path.join(dataset_path, inf)
    print(inf_full_path)
    region = 'barrel' if barrel else 'endcap'
    shell_file = '{}_{}_{}_{}.sh'.format(dataset_name, idx, region, particle)
    if particle == 'electron':
      if not barrel: 
        command = 'cmsRun testElectronMVA_cfg_mod1.py electronLabel=ecalDrivenGsfElectronsHGC inputFiles={} outDir={} outFileNumber={}'.format(inf_full_path, os.path.join(outdir, dataset_name), idx)
      else:
        command = 'cmsRun testElectronMVA_cfg_mod1.py inputFiles={} outDir={} outFileNumber={}'.format(inf_full_path, os.path.join(outdir, dataset_name), idx)
    else:
      if not barrel:
        command = 'cmsRun testPhotonMVA_cfg_mod1.py photonLabel=photonsHGC inputFiles={} outDir={} outFileNumber={}'.format(inf_full_path, os.path.join(outdir, dataset_name), idx)
      else:
        command = 'cmsRun testPhotonMVA_cfg_mod1.py inputFiles={} outDir={} outFileNumber={}'.format(inf_full_path, os.path.join(outdir, dataset_name), idx)


    prepare_shell(shell_file, command, condor, farm_dir)

if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--JobFlavour', dest='JobFlavour', help='espresso/microcentury/longlunch/workday/tomorrow', type=str, default='longlunch')
  parser.add_argument('--universe',   dest='universe',   help='vanilla/local', type=str, default='vanilla')
  parser.add_argument('--outdir',     dest='outdir',     help='output directory', type=str, default='./')
  parser.add_argument('--check',      action='store_true')
  parser.add_argument('--test',       dest='test',       action='store_true')
  parser.add_argument('--particle',   type=str, default = 'electron')
  args = parser.parse_args()
  args_dict = vars(args)

  ######################
  ##  Farm Directory  ##
  ######################

  farm_dir = os.path.join('./', 'Farm')
  os.system('mkdir -p %s '%farm_dir)

  ##############
  ##  Condor  ##  
  ##############

  condor = open(os.path.join(farm_dir, 'condor.sub'), 'w')
  condor.write('output = %s/job_common_$(Process).out\n'%farm_dir)
  condor.write('error  = %s/job_common_$(Process).err\n'%farm_dir)
  condor.write('log    = %s/job_common_$(Process).log\n'%farm_dir)
  condor.write('executable = %s/$(cfgFile)\n'%farm_dir)
  condor.write('universe = %s\n'%args.universe)
  condor.write('requirements = (OpSysAndVer =?= "CentOS7")\n')
  condor.write('on_exit_remove   = (ExitBySignal == False) && (ExitCode == 0)\n')
  condor.write('max_retries = 3\n')
  condor.write('+JobFlavour = "%s"\n'%args.JobFlavour)
  condor.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
  condor.write('use_x509userproxy = True\n')

  args.outdir = os.path.join(args.outdir, args.particle)

  if args.particle == 'electron':
    wrapper_condor('DY', "root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DY/240120_114733/0000", condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
    wrapper_condor('TT', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/TT_TuneCP5_14TeV-powheg-pythia8/crab_TT/240121_162232/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
    wrapper_condor('DY', "root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DY/240120_114733/0000", condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)
    wrapper_condor('TT', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/TT_TuneCP5_14TeV-powheg-pythia8/crab_TT/240121_162232/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)
    wrapper_condor('DY_noPU', "root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DY_noPU/240426_115232/0000", condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
    wrapper_condor('TT_noPU', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/TT_TuneCP5_14TeV-powheg-pythia8/crab_TT_noPU/240426_115307/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
    wrapper_condor('DY_noPU', "root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DY_noPU/240426_115232/0000", condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)
    wrapper_condor('TT_noPU', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/TT_TuneCP5_14TeV-powheg-pythia8/crab_TT_noPU/240426_115307/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)
  else:
    wrapper_condor('SinglePhoton2To200', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/SinglePhoton_Pt-2To200-gun/crab_SinglePhoton2to200/240303_225556/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('SinglePhoton200To500', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/SinglePhoton_Pt200To500-gun/crab_SinglePhoton200to500/240303_222110/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel = False, check=args.check, particle=args.particle)
    wrapper_condor('QCD15to20', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-15To20_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD15to20/240303_221527/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('QCD20to30', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD20to30/240303_221604/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('QCD30to50', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-30To50_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD30to50/240303_231625/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('QCD50to80', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-50To80_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD50to80/240303_221721/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('QCD80to120', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-80To120_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD80to120/240303_221759/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('QCD120to170', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-120To170_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD120to170/240303_221837/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('QCD170to300', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-170To300_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD170to300/240303_221915/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('QCD300toInf', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-300ToInf_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD300toInf/240303_221952/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)

    wrapper_condor('SinglePhoton2To200', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/SinglePhoton_Pt-2To200-gun/crab_SinglePhoton2to200/240303_225556/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    wrapper_condor('SinglePhoton200To500', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/SinglePhoton_Pt200To500-gun/crab_SinglePhoton200to500/240303_222110/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel = True, check=args.check, particle=args.particle)
    wrapper_condor('QCD15to20', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-15To20_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD15to20/240303_221527/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    wrapper_condor('QCD20to30', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD20to30/240303_221604/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    wrapper_condor('QCD30to50', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-30To50_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD30to50/240303_231625/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    wrapper_condor('QCD50to80', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-50To80_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD50to80/240303_221721/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    wrapper_condor('QCD80to120', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-80To120_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD80to120/240303_221759/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    wrapper_condor('QCD120to170', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-120To170_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD120to170/240303_221837/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    wrapper_condor('QCD170to300', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-170To300_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD170to300/240303_221915/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    wrapper_condor('QCD300toInf', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-300ToInf_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD300toInf/240303_221952/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)


  condor.close()
  if not args.test:
    os.system('condor_submit %s/condor.sub'%farm_dir)
