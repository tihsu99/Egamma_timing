import argparse
import os
import subprocess
from urllib.parse import urlparse

import yaml


def check_dir(path):
  if not os.path.exists(path):
    os.system("mkdir -p {}".format(path))


def resolve_path(path, base_dir):
  if path is None:
    return None
  if os.path.isabs(path):
    return os.path.abspath(path)
  return os.path.abspath(os.path.join(base_dir, path))


def resolve_cli_path(path):
  if path is None:
    return None
  if os.path.isabs(path):
    return os.path.abspath(path)
  return os.path.abspath(path)


def detect_cmssw_src():
  cmssw_base = os.environ.get("CMSSW_BASE")
  if cmssw_base:
    return os.path.join(cmssw_base, "src")

  current = os.getcwd()
  while True:
    if os.path.basename(current) == "src" and "CMSSW_" in os.path.basename(os.path.dirname(current)):
      return current
    parent = os.path.dirname(current)
    if parent == current:
      break
    current = parent
  raise RuntimeError("Unable to determine CMSSW src directory. Run inside a CMSSW area or set cmssw-dir in the config.")


def run_command(command):
  return subprocess.check_output(["bash", "-lc", command], text=True)


def write_condor_header(condor, farm_dir, log_dir, universe, job_flavour, proxy):
  condor.write("output = {}/job_common_$(Process).out\n".format(log_dir))
  condor.write("error  = {}/job_common_$(Process).err\n".format(log_dir))
  condor.write("log    = {}/job_common_$(Process).log\n".format(log_dir))
  condor.write("executable = {}/$(cfgFile)\n".format(farm_dir))
  condor.write("universe = {}\n".format(universe))
  condor.write("on_exit_remove   = (ExitBySignal == False) && (ExitCode == 0)\n")
  condor.write("max_retries = 3\n")
  condor.write('+JobFlavour = "{}"\n'.format(job_flavour))
  condor.write("x509userproxy = $ENV(X509_USER_PROXY)\n")
  condor.write("use_x509userproxy = True\n")
  if proxy is not None:
    condor.write("Proxy_path = {}\n".format(proxy))
    condor.write("arguments = $(Proxy_path)\n")


def detect_storage(dataset_path, storage_host):
  if dataset_path.startswith("root://"):
    parsed = urlparse(dataset_path)
    host = "{}://{}".format(parsed.scheme, parsed.netloc)
    remote_path = parsed.path
    return host, remote_path
  if dataset_path.startswith("/cms/store"):
    return storage_host, dataset_path
  return None, dataset_path


def is_das_dataset(sample_spec):
  if not isinstance(sample_spec, str):
    return False
  if sample_spec.startswith("/cms/store") or sample_spec.startswith("root://"):
    return False
  parts = [part for part in sample_spec.split("/") if part]
  return len(parts) == 3


def normalize_sample_spec(sample_spec):
  if isinstance(sample_spec, dict):
    return sample_spec
  if is_das_dataset(sample_spec):
    return {"dataset": sample_spec}
  return {"path": sample_spec}


def resolve_das_files(dataset_name, dbs_instance, redirector_host):
  query = "file dataset={} instance={}".format(dataset_name, dbs_instance)
  output = run_command("dasgoclient --query '{}'".format(query))
  host_prefix = redirector_host if redirector_host.endswith("//") else redirector_host + "//"
  files = []
  for line in output.splitlines():
    line = line.strip()
    if not line.endswith(".root"):
      continue
    if line.startswith("root://"):
      files.append(line)
    else:
      files.append("{}{}".format(host_prefix, line))
  files.sort()
  return files


def resolve_input_files(dataset_path, storage_host):
  host, path = detect_storage(dataset_path, storage_host)
  if host is None:
    output = run_command("ls {}".format(path))
    files = [os.path.join(path, item) for item in output.splitlines() if item.endswith(".root")]
    files.sort()
    return files

  xrdfs_host = host if host.endswith("//") else host + "//"
  output = run_command("xrdfs {} ls {}".format(xrdfs_host, path))
  files = []
  for line in output.splitlines():
    if not line.endswith(".root"):
      continue
    if line.startswith("/"):
      files.append("{}{}".format(xrdfs_host, line))
    else:
      files.append(line)
  files.sort()
  return files


def resolve_sample_files(sample_spec, das_redirector, storage_host):
  sample_cfg = normalize_sample_spec(sample_spec)
  if "dataset" in sample_cfg:
    dbs_instance = sample_cfg.get("dbs_instance", "prod/global")
    return resolve_das_files(sample_cfg["dataset"], dbs_instance, das_redirector)
  if "path" in sample_cfg:
    return resolve_input_files(sample_cfg["path"], storage_host)
  raise RuntimeError("Sample specification must contain either 'dataset' or 'path'")


def replace_placeholders(template, values):
  result = template
  for _ in range(5):
    updated = result
    for key, value in values.items():
      updated = updated.replace("{" + key + "}", str(value))
    if updated == result:
      break
    result = updated
  return result


def normalize_comma_separated(value):
  if value is None:
    return ""
  parts = [part.strip() for part in str(value).split(",")]
  return ",".join(part for part in parts if part)


def prepare_shell(shell_path, commands, cmssw_dir, workdir_expr):
  with open(shell_path, "w") as shell:
    shell.write("#!/bin/bash\n")
    shell.write("set -e\n")
    shell.write("if [ -n \"$1\" ]; then\n")
    shell.write("  export X509_USER_PROXY=$1\n")
    shell.write("fi\n")
    shell.write("WORKDIR={}\n".format(workdir_expr))
    shell.write("mkdir -p $WORKDIR\n")
    shell.write("echo 'Using WORKDIR=' $WORKDIR\n")
    shell.write("cd {}\n".format(cmssw_dir))
    shell.write("eval `scramv1 runtime -sh`\n")
    shell.write("cd $WORKDIR\n")
    for command in commands:
      shell.write(command)
      shell.write("\n")
    shell.write("rm -rf $WORKDIR\n")
  os.chmod(shell_path, 0o755)


def build_production_commands(config, args, particle, region, sample, variant, input_file, file_index, output_file):
  particle_cfg = config["particles"][particle]
  variant_cfg = config["variants"][variant]
  workflow_cfg = config["workflow"]
  offline_label = particle_cfg["offline_labels"][region]
  timing_dir = config["timing-dir"]
  comparison_cfg = particle_cfg["ntuplizer_cfg"]
  online_label = particle_cfg["online_label"]
  online_candidate_label = variant_cfg.get("online_candidate_label", "hltEgammaCandidatesUnseeded::HLTX")
  online_trackster_label = variant_cfg.get("online_trackster_label", "hltTiclCandidate::HLTX")
  comparison_input = workflow_cfg["ntuplizer_input"]

  placeholders = {
      "particle": particle,
      "region": region,
      "sample": sample,
      "variant": variant,
      "particle_type": particle.capitalize(),
      "proc_modifier": variant_cfg.get("proc_modifier", ""),
      "max_events": args.max_events if args.max_events is not None else -1,
      "max_events_arg": "-n {}".format(args.max_events) if args.max_events is not None else "",
      "input_file": input_file,
      "file_index": file_index,
      "timing_dir": timing_dir,
      "comparison_cfg": comparison_cfg,
      "offline_label": offline_label,
      "online_label": online_label,
      "online_candidate_label": online_candidate_label,
      "online_trackster_label": online_trackster_label,
      "comparison_input": comparison_input,
      "output_file": output_file,
      "workdir": "$WORKDIR",
      "l1_output": workflow_cfg.get("l1_output", "step_l1.root"),
      "hlt_output": workflow_cfg.get("hlt_output", "step_hlt.root"),
      "reco_output": workflow_cfg.get("reco_output", "step_reco.root"),
      "reco_output_commands": workflow_cfg.get("reco_output_commands", ""),
      "hlt_customise_commands": workflow_cfg.get("hlt_customise_commands", ""),
      "step1_customise": normalize_comma_separated(
          workflow_cfg.get("step1_customise", "SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000")
      ),
      }

  commands = []
  for command in workflow_cfg["commands"]:
    commands.append(replace_placeholders(command, placeholders))

  ntuplizer_input = replace_placeholders(comparison_input, placeholders)
  if particle == "electron":
    label_arg = "electronLabel={}".format(offline_label)
  else:
    label_arg = "photonLabel={}".format(offline_label)

  commands.append(
      "cmsRun {timing_dir}/python/{comparison_cfg} inputFiles={nt_input} outDir=$WORKDIR outFileNumber={file_index} "
      "{label_arg} offlineProcess=reRECO onlineLabel={online_label} onlineCandidateLabel={online_candidate_label} "
      "onlineTracksterLabel={online_trackster_label}".format(
          timing_dir=timing_dir,
          comparison_cfg=comparison_cfg,
          nt_input=ntuplizer_input,
          file_index=file_index,
          label_arg=label_arg,
          online_label=online_label,
          online_candidate_label=online_candidate_label,
          online_trackster_label=online_trackster_label,
      )
  )
  if args.max_events is not None:
    commands[-1] += " maxEvents={}".format(args.max_events)
  commands.append("mv $WORKDIR/comparisonNtuple_{}.root {}".format(file_index, output_file))
  return commands


def accept_sample(args, particle, region, sample):
  if args.particle is not None and particle != args.particle:
    return False
  if args.region is not None and region != args.region:
    return False
  if args.sample is not None and sample != args.sample:
    return False
  return True


def submit_production_jobs(config, condor, farm_dir, args):
  das_redirector = config.get("das-xrootd-host", "root://cms-xrd-global.cern.ch//")
  storage_host = config.get("storage-xrootd-host", config.get("xrootd-host", das_redirector))
  created_shells = []
  created_jobs = 0
  for particle, particle_cfg in config["samples"].items():
    for region, sample_map in particle_cfg.items():
      for sample, sample_spec in sample_map.items():
        if not accept_sample(args, particle, region, sample):
          continue
        input_files = resolve_sample_files(sample_spec, das_redirector, storage_host)
        for variant in args.variant:
          variant_outdir = os.path.join(args.outdir, variant, particle, region, sample)
          check_dir(variant_outdir)
          for file_index, input_file in enumerate(input_files):
            if args.max_files is not None and created_jobs >= args.max_files:
              return created_shells
            output_file = os.path.join(variant_outdir, "comparisonNtuple_{}.root".format(file_index))
            if args.check and os.path.exists(output_file):
              continue

            shell_name = "{}_{}_{}_{}_{}.sh".format(variant, particle, region, sample, file_index)
            shell_path = os.path.join(farm_dir, shell_name)
            workdir_expr = "${{TMPDIR:-/tmp}}/{}_{}_{}_{}_{}".format(args.tmp_tag, variant, particle, region, file_index)
            commands = build_production_commands(config, args, particle, region, sample, variant, input_file, file_index, output_file)
            prepare_shell(shell_path, commands, config["cmssw-dir"], workdir_expr)
            condor.write("cfgFile={}\n".format(shell_name))
            condor.write("queue 1\n")
            created_shells.append(shell_path)
            created_jobs += 1
  return created_shells


def discover_comparison_files(root_dir):
  mapping = dict()
  for particle in os.listdir(root_dir):
    particle_dir = os.path.join(root_dir, particle)
    if not os.path.isdir(particle_dir):
      continue
    for region in os.listdir(particle_dir):
      region_dir = os.path.join(particle_dir, region)
      if not os.path.isdir(region_dir):
        continue
      for sample in os.listdir(region_dir):
        sample_dir = os.path.join(region_dir, sample)
        if not os.path.isdir(sample_dir):
          continue
        for file_name in os.listdir(sample_dir):
          if file_name.endswith(".root"):
            mapping[(particle, region, sample, file_name)] = os.path.join(sample_dir, file_name)
  return mapping


def submit_merge_jobs(config, condor, farm_dir, args):
  v4_files = discover_comparison_files(args.v4_dir)
  v5_files = discover_comparison_files(args.v5_dir)
  if set(v4_files.keys()) != set(v5_files.keys()):
    raise RuntimeError("Mismatch between v4 and v5 comparison ntuples")

  timing_dir = config["timing-dir"]
  created_shells = []
  created_jobs = 0
  for (particle, region, sample, file_name), v4_file in sorted(v4_files.items()):
    if not accept_sample(args, particle, region, sample):
      continue
    if args.max_files is not None and created_jobs >= args.max_files:
      return created_shells
    v5_file = v5_files[(particle, region, sample, file_name)]
    output_dir = os.path.join(args.outdir, particle, region, sample)
    check_dir(output_dir)
    output_file = os.path.join(output_dir, file_name)
    if args.check and os.path.exists(output_file):
      continue

    shell_name = "merge_{}_{}_{}_{}.sh".format(particle, region, sample, file_name.replace(".root", ""))
    shell_path = os.path.join(farm_dir, shell_name)
    workdir_expr = "${{TMPDIR:-/tmp}}/{}_merge_{}_{}_{}".format(args.tmp_tag, particle, region, sample)
    commands = [
        "python3 {}/analysis/merge_comparison_trees.py --v4 {} --v5 {} --out $WORKDIR/{}".format(
            timing_dir, v4_file, v5_file, file_name
        ),
        "mv $WORKDIR/{} {}".format(file_name, output_file),
    ]
    prepare_shell(shell_path, commands, config["cmssw-dir"], workdir_expr)
    condor.write("cfgFile={}\n".format(shell_name))
    condor.write("queue 1\n")
    created_shells.append(shell_path)
    created_jobs += 1
  return created_shells


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Prepare Condor jobs for comparison ntuple production and merge with /tmp staging")
  parser.add_argument("--config", required=True, type=str)
  parser.add_argument("--sample-config", type=str, default=None)
  parser.add_argument("--mode", choices=["produce", "merge"], default="produce")
  parser.add_argument("--variant", nargs="+", choices=["v4", "v5"], default=["v4", "v5"])
  parser.add_argument("--farm", default="FarmComparison", type=str)
  parser.add_argument("--logdir", default=None, type=str)
  parser.add_argument("--outdir", required=True, type=str)
  parser.add_argument("--v4_dir", type=str)
  parser.add_argument("--v5_dir", type=str)
  parser.add_argument("--jobFlavour", default="workday", type=str)
  parser.add_argument("--universe", default="vanilla", type=str)
  parser.add_argument("--proxy", default=None, type=str)
  parser.add_argument("--tmp-tag", default="EGammaTimingComparison", type=str)
  parser.add_argument("--particle", choices=["electron", "photon"], default=None)
  parser.add_argument("--region", default=None, type=str)
  parser.add_argument("--sample", default=None, type=str)
  parser.add_argument("--max-files", default=None, type=int)
  parser.add_argument("--max-events", default=None, type=int)
  parser.add_argument("--check", action="store_true")
  parser.add_argument("--test", action="store_true")
  args = parser.parse_args()

  with open(args.config) as handle:
    config = yaml.safe_load(handle)

  config_dir = os.path.dirname(os.path.abspath(args.config))
  if config.get("cmssw-dir"):
    config["cmssw-dir"] = resolve_path(config["cmssw-dir"], config_dir)
  else:
    config["cmssw-dir"] = detect_cmssw_src()
  config["timing-dir"] = resolve_path(config.get("timing-dir", ".."), config_dir)
  if args.sample_config is not None:
    sample_config_path = resolve_cli_path(args.sample_config)
  else:
    sample_config_path = config.get("sample-config")
    if sample_config_path is not None:
      sample_config_path = resolve_path(sample_config_path, config_dir)
  if sample_config_path is not None:
    with open(sample_config_path) as handle:
      config["samples"] = yaml.safe_load(handle)["samples"]
  elif "samples" not in config:
    raise RuntimeError("No samples found. Provide --sample-config or set sample-config/samples in the workflow config.")

  args.farm = os.path.abspath(args.farm)
  if args.logdir is None:
    args.logdir = args.farm
  else:
    args.logdir = os.path.abspath(args.logdir)
  args.outdir = os.path.abspath(args.outdir)
  if args.v4_dir is not None:
    args.v4_dir = os.path.abspath(args.v4_dir)
  if args.v5_dir is not None:
    args.v5_dir = os.path.abspath(args.v5_dir)

  check_dir(args.farm)
  check_dir(args.logdir)
  condor_path = os.path.join(args.farm, "condor.sub")
  created_shells = []
  with open(condor_path, "w") as condor:
    write_condor_header(condor, args.farm, args.logdir, args.universe, args.jobFlavour, args.proxy)
    if args.mode == "produce":
      created_shells = submit_production_jobs(config, condor, args.farm, args)
    else:
      if args.v4_dir is None or args.v5_dir is None:
        raise RuntimeError("--v4_dir and --v5_dir are required in merge mode")
      created_shells = submit_merge_jobs(config, condor, args.farm, args)

  if len(created_shells):
    print("Generated {} shell script(s):".format(len(created_shells)))
    for shell_path in created_shells:
      print(shell_path)
  else:
    print("No shell scripts were generated.")

  if not args.test:
    os.system("condor_submit {}".format(condor_path))
