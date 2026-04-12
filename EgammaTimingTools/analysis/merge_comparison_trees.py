import argparse
import os

import awkward as ak
import uproot


def CheckDir(path):
  if not os.path.exists(path):
    os.system("mkdir -p {}".format(path))


def is_rootcompat(array):
  array_type = ak.type(array)
  if isinstance(array_type, ak._ext.ArrayType):
    if isinstance(array_type.type, ak._ext.PrimitiveType):
      return True
    if isinstance(array_type.type, ak._ext.ListType) and isinstance(array_type.type.type, ak._ext.PrimitiveType):
      return True
  return False


def uproot_writeable(events):
  output = dict()
  for branch_name in events.fields:
    if events[branch_name].fields:
      nested = dict()
      for field in events[branch_name].fields:
        if is_rootcompat(events[branch_name][field]):
          nested[field] = events[branch_name][field]
      output[branch_name] = ak.packed(ak.without_parameters(ak.zip(nested)))
    else:
      if "offset" in branch_name:
        continue
      output[branch_name] = ak.packed(ak.without_parameters(events[branch_name]))
  return output


def read_tree(file_name, tree_name):
  return uproot.open(file_name)[tree_name].arrays(library="ak")


def key_tuple(events, index):
  return (
      int(events["nRun"][index]),
      int(events["nLumi"][index]),
      int(events["nEvent"][index]),
      int(events["genIdx"][index]),
  )


def merge_gen_trees(v4_events, v5_events):
  v5_index = dict()
  for index in range(len(v5_events["genIdx"])):
    v5_index[key_tuple(v5_events, index)] = index

  output = dict()
  key_branches = {"nRun", "nLumi", "nEvent", "genIdx", "truthClass", "gen_pt", "gen_eta", "gen_phi"}

  for branch in v4_events.fields:
    if branch in key_branches:
      output[branch] = v4_events[branch]
    else:
      output["v4_" + branch] = v4_events[branch]

  for branch in v5_events.fields:
    if branch in key_branches:
      continue
    output["v5_" + branch] = ak.Array([v5_events[branch][v5_index[key_tuple(v4_events, idx)]] for idx in range(len(v4_events["genIdx"]))])

  if len(v4_events["genIdx"]) != len(v5_events["genIdx"]):
    raise RuntimeError("v4 and v5 gen trees have different number of entries")

  for idx in range(len(v4_events["genIdx"])):
    key = key_tuple(v4_events, idx)
    if key not in v5_index:
      raise RuntimeError("Missing v5 gen match for key {}".format(key))

  return output


def merge_file(v4_file, v5_file, output_file):
  v4_gen = read_tree(v4_file, "genTree")
  v5_gen = read_tree(v5_file, "genTree")
  v4_offline = read_tree(v4_file, "offlineTree")
  v5_offline = read_tree(v5_file, "offlineTree")
  v4_online = read_tree(v4_file, "onlineTree")
  v5_online = read_tree(v5_file, "onlineTree")

  output_dir = os.path.dirname(output_file)
  if len(output_dir):
    CheckDir(output_dir)
  fout = uproot.recreate(output_file)
  fout["gen_oriented"] = uproot_writeable(ak.zip(merge_gen_trees(v4_gen, v5_gen), depth_limit=1))
  fout["offline_v4"] = uproot_writeable(v4_offline)
  fout["offline_v5"] = uproot_writeable(v5_offline)
  fout["online_v4"] = uproot_writeable(v4_online)
  fout["online_v5"] = uproot_writeable(v5_online)


def discover_inputs(root_dir):
  mapping = dict()
  for region in os.listdir(root_dir):
    region_dir = os.path.join(root_dir, region)
    if not os.path.isdir(region_dir):
      continue
    for sample in os.listdir(region_dir):
      sample_dir = os.path.join(region_dir, sample)
      if not os.path.isdir(sample_dir):
        continue
      for file_name in os.listdir(sample_dir):
        if file_name.endswith(".root"):
          mapping[(region, sample, file_name)] = os.path.join(sample_dir, file_name)
  return mapping


def merge_directory(indir_v4, indir_v5, outdir):
  v4_files = discover_inputs(indir_v4)
  v5_files = discover_inputs(indir_v5)

  if set(v4_files.keys()) != set(v5_files.keys()):
    missing_v4 = sorted(set(v5_files.keys()) - set(v4_files.keys()))
    missing_v5 = sorted(set(v4_files.keys()) - set(v5_files.keys()))
    raise RuntimeError("File mismatch between v4 and v5 directories: missing_v4={}, missing_v5={}".format(missing_v4, missing_v5))

  for (region, sample, file_name), v4_file in sorted(v4_files.items()):
    v5_file = v5_files[(region, sample, file_name)]
    output_dir = os.path.join(outdir, region, sample)
    CheckDir(output_dir)
    merge_file(v4_file, v5_file, os.path.join(output_dir, file_name))


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Merge v4/v5 comparison ntuples into the final five-tree comparison file")
  parser.add_argument("--v4", type=str)
  parser.add_argument("--v5", type=str)
  parser.add_argument("--out", type=str)
  parser.add_argument("--indir_v4", type=str)
  parser.add_argument("--indir_v5", type=str)
  parser.add_argument("--outdir", type=str)
  args = parser.parse_args()

  if args.v4 and args.v5 and args.out:
    merge_file(args.v4, args.v5, args.out)
  elif args.indir_v4 and args.indir_v5 and args.outdir:
    merge_directory(args.indir_v4, args.indir_v5, args.outdir)
  else:
    raise RuntimeError("Provide either --v4/--v5/--out or --indir_v4/--indir_v5/--outdir")
