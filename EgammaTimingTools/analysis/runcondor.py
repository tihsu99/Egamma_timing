import os, sys
import optparse, argparse
import ROOT

cwd = os.getcwd()

def prepare_shell(shell_file, command, condor, FarmDir):
  cwd = os.getcwd()
  with open(os.path.join(FarmDir, shell_file), 'w') as shell:
    shell.write('#!/bin/bash\n')
    shell.write('WORKDIR=%s\n'%cwd)
    shell.write('cd ${WORKDIR}\n')
    shell.write('source script/env.sh\n')
    shell.write(command)
  condor.write('cfgFile=%s\n'%shell_file)
  condor.write('queue 1\n')

if __name__ == '__main__':
    usage = 'usage: %prog [options]'
    parser = argparse.ArgumentParser(description = usage)
    parser.add_argument('--indir', type=str)
    parser.add_argument('--outdir', type=str)
    parser.add_argument('--JobFlavour', dest='JobFlavour', help='espresso/microcentury/longlunch/workday/tomorrow', type=str, default='longlunch')
    parser.add_argument('--universe',   dest='universe',   help='vanilla/local', type=str, default='vanilla')
    parser.add_argument('--check', action = 'store_true')
    parser.add_argument('--test', action = 'store_true')
    parser.add_argument('--particle', type=str)
    args = parser.parse_args()
  
    farm_dir = "Farm"
    os.system('mkdir -p {}'.format(farm_dir))
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

    region = ['endcap']
    GreenLight = True
    for region_ in region:
      INDIR = os.path.join(args.indir, args.particle, region_)
      samples = os.listdir(INDIR)
      for sample_ in samples:
        SAMPLE_DIR = os.path.join(INDIR, sample_)
        files = os.listdir(SAMPLE_DIR)
        for file_ in files:
          if not "ntuple" in file_: continue
          input_file = os.path.join(SAMPLE_DIR, file_)
          output_file = os.path.join(args.outdir, args.particle, region_, sample_, file_)
          threshold   = 0.3 if "noPU" in sample_ else 0.7
          command = "python3 skim_photon_track.py --fin {} --fout {} --threshold {}\n".format(input_file, output_file, threshold)
          shell_name = "{}_{}_{}_{}.sh".format(region_, args.particle, sample_, file_)
          if args.check:
            try:
              fout_root = ROOT.TFile.Open(output_file, "READ")
              if fout_root.IsZombie():
                print("{} is Zombie.".format(output_file))
              h = fout_root.Get("ntuplizer/tree")
              if args.particle == 'electron':
                a = h.ele_pt
              else:
                a = h.Pho_eta
              fout_root.Close()
              continue
            except Exception as e:
              print(e)
              print("reproduce {}".format(output_file))
              GreenLight = False
          prepare_shell(shell_name, command, condor, farm_dir)

    condor.close()
    if not(args.test or (args.check and GreenLight)):
      os.system("condor_submit Farm/condor.sub")
