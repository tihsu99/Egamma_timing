#!/usr/bin/env python3
import argparse
import glob
import os

import awkward as ak
import numpy as np
import pandas as pd
import uproot


SCENARIOS = ["online_v4", "online_v5", "offline_v4", "offline_v5"]
COMMON_SCENARIOS = {"online_v4", "online_v5", "offline_v4", "offline_v5"}
ONLINE_SCENARIOS = {"online_v4", "online_v5"}
OFFLINE_SCENARIOS = {"offline_v4", "offline_v5"}

DT_CUTS = [round(0.01 * (index + 1), 3) for index in range(20)] + [9999.0]
VETO_MODES = ["cone", "seed", "cluster"]
COLLECTIONS = [
    ("HGCRecHit", "energy", 1.0),
    ("trackster", "pt", 1.0),
    ("layerCluster", "energy", 0.0),
]


def check_dir(path):
  if not os.path.exists(path):
    os.system("mkdir -p {}".format(path))


def object_prefix(particle):
  if particle == "electron":
    return "Ele", "matchedToGenEle"
  return "Pho", "matchedToGenPho"


def scenario_files(input_dir, particle, region, process_name):
  pattern = os.path.join(input_dir, particle, region, process_name, "*.root")
  return sorted(glob.glob(pattern))


def resolve_input_files(input_dir, input_file, particle, region, process_name):
  if input_file:
    return [input_file]
  return scenario_files(input_dir, particle, region, process_name)


def read_tree(files, tree_name):
  if not files:
    raise RuntimeError("No ROOT files found for tree {}".format(tree_name))
  print("[read] tree={} files={}".format(tree_name, len(files)))
  return uproot.concatenate([f"{name}:{tree_name}" for name in files], library="ak")


def absolute_time_difference(values, reference):
  broadcast_ref = ak.broadcast_arrays(reference, values)[0]
  valid = (values > -90) & (broadcast_ref > -90)
  return ak.where(valid, abs(values - broadcast_ref), 999.0)


def seed_time_from_tracksters(trackster_pt, trackster_dr, trackster_time, trackster_is_seed):
  seed_tracksters = trackster_is_seed & (trackster_pt > 1.0) & (trackster_dr < 0.3) & (trackster_time > -90.0)
  weighted = ak.mean(trackster_time[seed_tracksters], axis=1)
  return ak.fill_none(weighted, -99.0)


def compute_isolation(dr, values, is_seed, in_cluster, time_wrt_ref, ext_dr=0.15, int_dr=0.0,
                      energy_cut=0.0, dt_cut=9999.0, veto_seed=False, veto_cluster=False):
  selection = (dr < ext_dr) & (dr > int_dr) & (values > energy_cut)
  if time_wrt_ref is not None:
    selection = selection & (time_wrt_ref < dt_cut)
  if veto_seed:
    selection = selection & (~is_seed)
  if veto_cluster:
    selection = selection & (~in_cluster)
  return ak.sum(values[selection], axis=-1)


def flatten_valid(values, mask=None):
  if mask is None:
    return ak.to_numpy(ak.flatten(values))
  return ak.to_numpy(ak.flatten(values[mask]))


def scenario_collections(events, scenario, particle):
  prefix, matched_branch = object_prefix(particle)
  payload = {
      "scenario": scenario,
      "eg_pt": ak.to_numpy(events[f"{prefix}_pt"]),
      "eg_eta": ak.to_numpy(events[f"{prefix}_eta"]),
      "eg_phi": ak.to_numpy(events[f"{prefix}_phi"]),
      "matched": ak.to_numpy(events[matched_branch]),
      "gen_idx": ak.to_numpy(events["genIdx"]),
  }

  trackster_time = events["Trackster_Time"]
  seed_time = seed_time_from_tracksters(
      events["Trackster_pt"],
      events["Trackster_dr"],
      trackster_time,
      events["Trackster_isSeed"],
  )
  payload["seed_time"] = ak.to_numpy(seed_time)

  collections = {
      "trackster": {
          "dr": events["Trackster_dr"],
          "value": events["Trackster_pt"],
          "time": events["Trackster_Time"],
          "isSeed": events["Trackster_isSeed"],
          "inCluster": events["Trackster_inCluster"],
      },
      "HGCRecHit": {
          "dr": events["HGCRecHit_dr"],
          "value": events["HGCRecHit_energy"],
          "time": events["HGCRecHit_Time"],
          "isSeed": events["HGCRecHit_isSeed"],
          "inCluster": events["HGCRecHit_inCluster"],
      },
  }

  if scenario in ONLINE_SCENARIOS and "LayerCluster_Time" in events.fields:
    collections["layerCluster"] = {
        "dr": events["LayerCluster_dr"],
        "value": events["LayerCluster_energy"],
        "time": events["LayerCluster_Time"],
        "isSeed": events["LayerCluster_isSeed"],
        "inCluster": events["LayerCluster_inCluster"],
    }

  references = {
      "seed": seed_time,
  }
  if scenario in OFFLINE_SCENARIOS:
    references["pv"] = events["PV_Time"]
    references["sigTrk"] = events["sigTrkTime"]

  for collection_name, collection in collections.items():
    for ref_name, ref_values in references.items():
      collection[f"timeWrt_{ref_name}"] = absolute_time_difference(collection["time"], ref_values)

  if scenario in OFFLINE_SCENARIOS:
    payload["pv_time"] = ak.to_numpy(events["PV_Time"])
    payload["sig_trk_time"] = ak.to_numpy(events["sigTrkTime"])
    payload["sig_trk_time_err"] = ak.to_numpy(events["sigTrkTimeErr"])
    payload["sig_trk_mtdmva"] = ak.to_numpy(events["sigTrkMtdMva"])

  if scenario in ONLINE_SCENARIOS:
    for branch in [
        "hgcalPFIsol_default",
        "ecalPFIsol_default",
        "hcalPFIsol_default",
        "r9_default",
        "sigmaIEtaIEta_default",
        "hForHoverE_default",
    ]:
      if branch in events.fields:
        payload[branch] = ak.to_numpy(events[branch])
    for branch in ["trkIsol_default", "trkIsolV0_default", "trkIsolV6_default", "trkIsolV72_default"]:
      if branch in events.fields:
        payload[branch] = ak.to_numpy(events[branch])

  return payload, collections


def append_timing_summaries(payload, collections):
  for collection_name, collection in collections.items():
    for ref_name in [key.replace("timeWrt_", "") for key in collection.keys() if key.startswith("timeWrt_")]:
      time_name = f"timeWrt_{ref_name}"
      values = collection[time_name]
      valid = (collection["time"] > -90) & (values < 900)
      payload[f"{collection_name}-timeabsmean-ref_{ref_name}"] = ak.to_numpy(
          ak.fill_none(ak.mean(values[valid], axis=1), -1.0)
      )
      payload[f"{collection_name}-nhit-ref_{ref_name}"] = ak.to_numpy(ak.sum(valid, axis=1))


def append_isolation_columns(payload, collections, scenario):
  for collection_name, energy_name, base_energy_cut in COLLECTIONS:
    if collection_name not in collections:
      continue
    collection = collections[collection_name]
    reference_names = [key.replace("timeWrt_", "") for key in collection.keys() if key.startswith("timeWrt_")]
    for veto in VETO_MODES:
      inner_dr = 0.02 if veto == "cone" else 0.0
      energy_cut = 0.02 if (collection_name == "layerCluster" and veto == "cone") else base_energy_cut
      veto_seed = veto == "seed"
      veto_cluster = veto == "cluster"
      for ref_name in reference_names:
        ref_key = f"timeWrt_{ref_name}"
        for dt_cut in DT_CUTS:
          branch_name = (
              f"eg-{collection_name}_iso-{veto}_veto-ref_{ref_name}-dt_{dt_cut:.3f}"
          )
          payload[branch_name] = ak.to_numpy(
              compute_isolation(
                  collection["dr"],
                  collection["value"],
                  collection["isSeed"],
                  collection["inCluster"],
                  collection[ref_key],
                  ext_dr=0.15,
                  int_dr=inner_dr,
                  energy_cut=energy_cut,
                  dt_cut=dt_cut,
                  veto_seed=veto_seed,
                  veto_cluster=veto_cluster,
              )
          )

  if scenario in ONLINE_SCENARIOS:
    for branch in [
        "hgcalPFIsol_default",
        "ecalPFIsol_default",
        "hcalPFIsol_default",
        "trkIsol_default",
        "trkIsolV0_default",
        "trkIsolV6_default",
        "trkIsolV72_default",
    ]:
      if branch in payload:
        payload[f"eg-{branch}"] = payload[branch]


def build_dataframe(events, scenario, particle, process_name, is_signal):
  payload, collections = scenario_collections(events, scenario, particle)
  append_timing_summaries(payload, collections)
  append_isolation_columns(payload, collections, scenario)
  df = pd.DataFrame(payload)
  df["process"] = process_name
  df["is_signal"] = is_signal
  df["weight"] = 1.0

  if is_signal:
    df = df[df["matched"] == 1].copy()
  else:
    df = df[df["matched"] != 1].copy()

  df.rename(columns={
      "eg_pt": "eg-pt",
      "eg_eta": "eg-eta",
      "eg_phi": "eg-phi",
      "seed_time": "eg-seedTime",
  }, inplace=True)

  return df


def write_scenario_parquet(input_dir, input_file, out_dir, particle, region, scenario, process_name, is_signal):
  files = resolve_input_files(input_dir, input_file, particle, region, process_name)
  print("[build] scenario={} process={} signal={}".format(scenario, process_name, is_signal))
  events = read_tree(files, scenario)
  df = build_dataframe(events, scenario, particle, process_name, is_signal)
  scenario_dir = os.path.join(out_dir, particle, region, scenario)
  check_dir(scenario_dir)
  out_path = os.path.join(scenario_dir, f"{process_name}_iso.parquet")
  df.to_parquet(out_path, index=False)
  print("[saved] {} entries={}".format(out_path, len(df)))
  return out_path


def main():
  parser = argparse.ArgumentParser(description="Build HGCal timing/isolation parquet tables from merged comparison ntuples")
  parser.add_argument("--inputDir", type=str, help="Merged ROOT base directory")
  parser.add_argument("--inputFile", type=str, help="Single merged ROOT file")
  parser.add_argument("--outDir", required=True, type=str, help="Output parquet directory")
  parser.add_argument("--particle", required=True, choices=["electron", "photon"])
  parser.add_argument("--region", required=True, type=str)
  parser.add_argument("--signal", required=True, type=str)
  parser.add_argument("--background", required=True, type=str)
  parser.add_argument("--scenario", action="append", choices=SCENARIOS, help="Optional scenario filter")
  args = parser.parse_args()

  if not args.inputDir and not args.inputFile:
    raise RuntimeError("Provide either --inputDir or --inputFile")

  scenarios = args.scenario if args.scenario else SCENARIOS
  for scenario in scenarios:
    write_scenario_parquet(args.inputDir, args.inputFile, args.outDir, args.particle, args.region, scenario, args.signal, True)
    write_scenario_parquet(args.inputDir, args.inputFile, args.outDir, args.particle, args.region, scenario, args.background, False)


if __name__ == "__main__":
  main()
