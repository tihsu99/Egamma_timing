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
#    shell.write('eval `scram r -sh`\n')
    shell.write('cd ${WORKDIR}\n')
    shell.write(command)
  condor.write('cfgFile=%s\n'%shell_file)
  condor.write('queue 1\n')

def wrapper_condor(dataset_name, dataset_path, condor, outdir, farm_dir, barrel=False, check=False, particle='electron'):
  print('start to process {}'.format(dataset_name))
  os.system('mkdir -p {}'.format(os.path.join(outdir, dataset_name)))
  print("xrdfs root://se01.grid.nchc.org.tw// ls " + dataset_path)
  os.system("xrdfs root://se01.grid.nchc.org.tw// ls " + dataset_path)
  result = subprocess.Popen("xrdfs root://se01.grid.nchc.org.tw// ls " + dataset_path, shell=True, stdout=subprocess.PIPE)
  subprocess_return = result.stdout.read()
  file_list = subprocess_return.decode('utf-8').split("\n")
  n_file = len(file_list)
  print(file_list)
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
          a = h.nEvent
        else:
          a = h.Pho_eta
    #    print(h)
        fout_root.Close()
        continue
      except Exception as error:
        print(error)
        print("reproduce {}".format(fout))
    print(idx, inf)
    inf_full_path = ('root://se01.grid.nchc.org.tw//' +  inf)
    print(inf_full_path)
    tmp_out_path = "/tmp/tihsu"
    os.system(f"mkdir -p {tmp_out_path}/{dataset_name}")
    region = 'barrel' if barrel else 'endcap'
    shell_file = '{}_{}_{}_{}.sh'.format(dataset_name, idx, region, particle)
    command = "eval `scram r -sh`\n"
    if particle == 'electron':
      if not barrel: 
        command += 'cmsRun ElectronTrackerNtuplizer.py electronLabel=ecalDrivenGsfElectronsHGC inputFiles={} outDir={} outFileNumber={}\n'.format(inf_full_path, os.path.join(tmp_out_path, dataset_name), idx)
      else:
        command += 'cmsRun ElectronTrackerNtuplizer.py inputFiles={} outDir={} outFileNumber={}\n'.format(inf_full_path, os.path.join(tmp_out_path, dataset_name), idx)
    else:
      if not barrel:
        command += 'cmsRun PhotonTrackerNtuplizer.py photonLabel=photonsHGC inputFiles={} outDir={} outFileNumber={}\n'.format(inf_full_path, os.path.join(tmp_out_path, dataset_name), idx)
      else:
        command += 'cmsRun testPhotonMVA_cfg_mod1.py inputFiles={} outDir={} outFileNumber={}\n'.format(inf_full_path, os.path.join(tmp_out_path, dataset_name), idx)

    command += f"mv {tmp_out_path}/{dataset_name}/ntupleTree_{idx}.root {fout}\n"
    prepare_shell(shell_file, command, condor, farm_dir)

if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--JobFlavour', dest='JobFlavour', help='espresso/microcentury/longlunch/workday/tomorrow', type=str, default='workday')
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
  condor.write('on_exit_remove   = (ExitBySignal == False) && (ExitCode == 0)\n')
  condor.write('max_retries = 3\n')
  condor.write('+JobFlavour = "%s"\n'%args.JobFlavour)
  condor.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
  condor.write('use_x509userproxy = True\n')

  args.outdir = os.path.join(args.outdir, args.particle)

  
  DY_path = "/cms/store/user/tihsu/ele_reRECO/2024-12-02/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DY_wTICL/241202_225444/0000"
#  DY_noPU_path = "root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DY_wTICL_noPU/240502_130306/0000"
  TT_path = "/cms/store/user/tihsu/ele_reRECO/2024-12-02/TT_TuneCP5_14TeV-powheg-pythia8/crab_TT_wTICL/241202_225507/0000"
#  TT_noPU_path = "root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/TT_TuneCP5_14TeV-powheg-pythia8/crab_TT_wTICL_noPU/240502_130344/0000"

  if args.particle == 'electron':
    wrapper_condor('DY', DY_path, condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
    wrapper_condor('TT', TT_path, condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
 #   wrapper_condor('DY', DY_path, condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)
 #   wrapper_condor('TT', TT_path, condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)
 #   wrapper_condor('DY_noPU', DY_noPU_path, condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
 #   wrapper_condor('TT_noPU', TT_noPU_path, condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check)
 #   wrapper_condor('DY_noPU', DY_noPU_path, condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)
 #   wrapper_condor('TT_noPU', TT_noPU_path, condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check)
  else:
#    wrapper_condor('SinglePhoton2To200_noPU', '/cms/store/user/tihsu/ele_reRECO/2024-09-06/SinglePhoton_Pt-2To200-gun/crab_SinglePhoton2to200_wTICL_noPU/240906_171249/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('SinglePhoton2To200', '/cms/store/user/tihsu/ele_reRECO/2025-03-24/SinglePhoton_Pt-2To200-gun/crab_SinglePhoton2to200_wTICL/250324_222449/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
    wrapper_condor('QCDEM', '/cms/store/user/tihsu/ele_reRECO/2025-03-24/QCD_Pt-15To3000_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCDEMEnriched_Pt-15To3000_wTICL/250324_222511/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
#    wrapper_condor('QCDEM_noPU', '/cms/store/user/tihsu/ele_reRECO/2024-09-06/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/crab_QCD_Pt-15To3000_wTICL_noPU/240906_171359/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
#    wrapper_condor('SinglePhoton200To500', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/SinglePhoton_Pt200To500-gun/crab_SinglePhoton200to500/240303_222110/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel = False, check=args.check, particle=args.particle)
#    wrapper_condor('QCD15to20', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-15To20_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD15to20/240303_221527/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
#    wrapper_condor('QCD20to30', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD20to30/240303_221604/0000', condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
#    wrapper_condor('QCD30to50', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-30To50_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD30to50/240303_231625/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
#    wrapper_condor('QCD50to80', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-50To80_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD50to80/240303_221721/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
   # wrapper_condor('QCD80to120', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-80To120_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD80to120/240303_221759/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
   # wrapper_condor('QCD120to170', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-120To170_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD120to170/240303_221837/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
   # wrapper_condor('QCD170to300', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-170To300_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD170to300/240303_221915/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)
   # wrapper_condor('QCD300toInf', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-300ToInf_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD300toInf/240303_221952/0000',  condor, os.path.join(args.outdir, 'endcap'), farm_dir, barrel=False, check=args.check, particle=args.particle)

    #wrapper_condor('SinglePhoton2To200', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/SinglePhoton_Pt-2To200-gun/crab_SinglePhoton2to200/240303_225556/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
   # wrapper_condor('SinglePhoton200To500', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/SinglePhoton_Pt200To500-gun/crab_SinglePhoton200to500/240303_222110/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel = True, check=args.check, particle=args.particle)
   # wrapper_condor('QCD15to20', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-15To20_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD15to20/240303_221527/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
   # wrapper_condor('QCD20to30', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD20to30/240303_221604/0000', condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
   # wrapper_condor('QCD30to50', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-30To50_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD30to50/240303_231625/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
   # wrapper_condor('QCD50to80', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-50To80_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD50to80/240303_221721/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    #wrapper_condor('QCD80to120', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-80To120_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD80to120/240303_221759/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    #wrapper_condor('QCD120to170', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-120To170_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD120to170/240303_221837/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    #wrapper_condor('QCD170to300', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-170To300_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD170to300/240303_221915/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)
    #wrapper_condor('QCD300toInf', 'root://se01.grid.nchc.org.tw//cms/store/user/tihsu/ele_reRECO/2024-01-19/QCD_Pt-300ToInf_EMEnriched_TuneCP5_14TeV-pythia8/crab_QCD300toInf/240303_221952/0000',  condor, os.path.join(args.outdir, 'barrel'), farm_dir, barrel=True, check=args.check, particle=args.particle)


  condor.close()
  if not args.test:
    os.system('condor_submit %s/condor.sub'%farm_dir)
