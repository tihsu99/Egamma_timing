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

def wrapper_condor(dataset_name, dataset_path, condor, outdir, farm_dir, barrel=False, check=False):
  result = subprocess.Popen("gfal-ls " + dataset_path, shell=True, stdout=subprocess.PIPE)
  subprocess_return = str(result.stdout.read())
  file_list = str(subprocess_return).split("\n")
  n_file = len(file_list)
  for idx, inf in enumerate(file_list):
    if not ('.root' in inf): continue
    fout = os.path.join(outdir, dataset_name, 'ntupleTree_{}.root'.format(idx))
    print(idx, inf)
    if check:
      try:
        fout_root = ROOT.TFile.Open(fout, "READ")
        if fout_root.IsZombie():
          print("{} is Zombie.".format(fout))
        h = fout_root.Get("ntuplizer/tree")
        a = h.ele_pt
        print(h)
        fout_root.Close()
        continue
      except Exception as error:
        print(error)
        print("reproduce {}".format(fout))
    print(idx, inf)
    inf_full_path = os.path.join(dataset_path, inf)
    print(inf_full_path)
    region = 'barrel' if barrel else 'endcap'
    shell_file = '{}_{}_{}.sh'.format(dataset_name, idx, region)
    if barrel: 
      command = 'cmsRun testElectronMVA_cfg_mod1.py electronLabel=ecalDrivenGsfElectronsHGC inputFiles={} outDir={} outFileNumber={}'.format(inf_full_path, os.path.join(outdir, dataset_name), idx)
    else:
      command = 'cmsRun testElectronMVA_cfg_mod1.py inputFiles={} outDir={} outFileNumber={}'.format(inf_full_path, os.path.join(outdir, dataset_name), idx)

    prepare_shell(shell_file, command, condor, farm_dir)

if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--JobFlavour', dest='JobFlavour', help='espresso/microcentury/longlunch/workday/tomorrow', type=str, default='longlunch')
  parser.add_argument('--universe',   dest='universe',   help='vanilla/local', type=str, default='vanilla')
  parser.add_argument('--outdir',     dest='outdir',     help='output directory', type=str, default='./')
  parser.add_argument('--check',      action='store_true')
  parser.add_argument('--test',       dest='test',       action='store_true')
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
  condor.write('+JobFlavour = "%s"\n'%args.JobFlavour)
  condor.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
  condor.write('use_x509userproxy = True\n')

  wrapper_condor('DY', "root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DY/240120_114733/0000", condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
  wrapper_condor('TT', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/TT_TuneCP5_14TeV-powheg-pythia8/crab_TT/240121_162232/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
  wrapper_condor('DY', "root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DY/240120_114733/0000", condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)
  wrapper_condor('TT', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/TT_TuneCP5_14TeV-powheg-pythia8/crab_TT/240121_162232/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)

  condor.close()
  if not args.test:
    os.system('condor_submit %s/condor.sub'%farm_dir)
