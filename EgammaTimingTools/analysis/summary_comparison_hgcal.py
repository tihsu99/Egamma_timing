#!/usr/bin/env python3
import argparse
import glob
import json
import os
import re

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from sklearn.metrics import auc, roc_curve

try:
  from tqdm import tqdm
except ImportError:
  tqdm = None


SCENARIOS = ["online_v4", "online_v5", "offline_v4", "offline_v5"]
ONLINE_SCENARIOS = {"online_v4", "online_v5"}
OFFLINE_SCENARIOS = {"offline_v4", "offline_v5"}
COMMON_REFERENCES = ["seed"]
OFFLINE_ONLY_REFERENCES = ["pv", "sigTrk"]
COLLECTION_LABELS = {
    "trackster": "Trackster",
    "HGCRecHit": "HGCRecHit",
    "layerCluster": "LayerCluster",
}


def check_dir(path):
  if not os.path.exists(path):
    os.system("mkdir -p {}".format(path))


def load_frame(iso_dir, particle, region, scenario, process_name):
  file_name = os.path.join(iso_dir, particle, region, scenario, f"{process_name}_iso.parquet")
  print("[load] {}".format(file_name))
  return pd.read_parquet(file_name)


def compute_weights(sig_df, bkg_df, pt_bins, eta_bins):
  sig_hist, _, _ = np.histogram2d(sig_df["eg-pt"], sig_df["eg-eta"], bins=[pt_bins, eta_bins])
  bkg_hist, _, _ = np.histogram2d(bkg_df["eg-pt"], bkg_df["eg-eta"], bins=[pt_bins, eta_bins])
  weights_2d = np.divide(bkg_hist, sig_hist, out=np.zeros_like(bkg_hist), where=sig_hist > 0)

  pt_idx = np.clip(np.digitize(sig_df["eg-pt"], pt_bins) - 1, 0, len(pt_bins) - 2)
  eta_idx = np.clip(np.digitize(sig_df["eg-eta"], eta_bins) - 1, 0, len(eta_bins) - 2)
  return weights_2d[pt_idx, eta_idx]


def apply_reweighting(sig_df, bkg_df):
  pt_bins = np.linspace(15, 200, 20)
  eta_bins = np.linspace(-3.2, 3.2, 33)
  sig_df = sig_df.copy()
  bkg_df = bkg_df.copy()
  sig_df["weight"] = compute_weights(sig_df, bkg_df, pt_bins, eta_bins)
  bkg_df["weight"] = 1.0
  return sig_df, bkg_df


def parse_dt_from_column(column_name):
  match = re.search(r"-dt_([0-9.]+)$", column_name)
  if match is None:
    return None
  return float(match.group(1))


def compute_roc(sig_df, bkg_df, variable):
  sig_vals = sig_df[variable].to_numpy()
  bkg_vals = bkg_df[variable].to_numpy()
  sig_weights = sig_df["weight"].to_numpy()
  bkg_weights = bkg_df["weight"].to_numpy()

  score = np.concatenate([sig_vals, bkg_vals])
  labels = np.concatenate([np.ones_like(sig_vals), np.zeros_like(bkg_vals)])
  weights = np.concatenate([sig_weights, bkg_weights])
  fpr, tpr, thresholds = roc_curve(labels, -score, sample_weight=weights)
  return fpr, tpr, thresholds, auc(fpr, tpr)


def default_hlt_branch(frame):
  for branch in ["eg-hgcalPFIsol_default", "eg-trkIsol_default", "eg-trkIsolV6_default", "eg-trkIsolV72_default"]:
    if branch in frame.columns:
      return branch
  return None


def iso_family_prefix(collection_name, veto_mode, reference_name):
  return f"eg-{collection_name}_iso-{veto_mode}_veto-ref_{reference_name}-dt_"


def scenario_records(sig_df, bkg_df, collection_name, veto_mode, reference_name):
  prefix = iso_family_prefix(collection_name, veto_mode, reference_name)
  candidates = [column for column in sig_df.columns if column.startswith(prefix)]
  candidates = [column for column in candidates if column in bkg_df.columns]
  records = []
  for column in candidates:
    dt_value = parse_dt_from_column(column)
    if dt_value is None:
      continue
    fpr, tpr, _, roc_auc = compute_roc(sig_df, bkg_df, column)
    records.append({
        "column": column,
        "dt": dt_value,
        "fpr": fpr,
        "tpr": tpr,
        "auc": roc_auc,
    })
  return sorted(records, key=lambda item: item["dt"])


def best_record(records):
  if not records:
    return None
  return max(records, key=lambda item: item["auc"])


def default_record(records):
  for record in records:
    if abs(record["dt"] - 9999.0) < 1e-6:
      return record
  return None


def plot_roc_by_scenario(frames, collection_name, veto_mode, reference_name, plotdir, include_default_online=False):
  available = {}
  summary = {}
  for scenario, (sig_df, bkg_df) in frames.items():
    records = scenario_records(sig_df, bkg_df, collection_name, veto_mode, reference_name)
    best = best_record(records)
    default = default_record(records)
    if best is None:
      continue
    available[scenario] = (best, default)
    summary[scenario] = {
        "best_auc": best["auc"],
        "best_dt": best["dt"],
        "default_auc": default["auc"] if default is not None else None,
    }

  if not available:
    return None

  plt.figure(figsize=(8, 6))
  colors = plt.get_cmap("tab10").colors
  for index, scenario in enumerate(SCENARIOS):
    if scenario not in available:
      continue
    best, default = available[scenario]
    color = colors[index % len(colors)]
    label = f"{scenario} best ({best['dt']:.3f})"
    plt.plot(best["tpr"], 1.0 - best["fpr"], color=color, linewidth=2, label=label)
    if default is not None:
      plt.plot(default["tpr"], 1.0 - default["fpr"], color=color, linewidth=1.5, linestyle="--", label=f"{scenario} no cut")

    if include_default_online and scenario in ONLINE_SCENARIOS:
      branch = default_hlt_branch(frames[scenario][0])
      if branch is not None and branch in frames[scenario][1].columns:
        fpr, tpr, _, _ = compute_roc(frames[scenario][0], frames[scenario][1], branch)
        plt.plot(tpr, 1.0 - fpr, color=color, linewidth=1.2, linestyle=":", label=f"{scenario} HLT default")

  plt.xlabel("Signal efficiency")
  plt.ylabel("Background rejection")
  plt.title(f"{COLLECTION_LABELS[collection_name]} {veto_mode} veto ({reference_name} reference)")
  plt.grid(alpha=0.3)
  plt.legend(fontsize=8, ncol=2)
  plt.tight_layout()
  out_name = f"roc_{collection_name}_{veto_mode}_ref_{reference_name}.png"
  plt.savefig(os.path.join(plotdir, out_name))
  plt.savefig(os.path.join(plotdir, out_name.replace(".png", ".pdf")))
  plt.close()
  return summary


def plot_auc_scan(frames, collection_name, veto_mode, reference_name, plotdir):
  plt.figure(figsize=(8, 6))
  drew = False
  for scenario in SCENARIOS:
    if scenario not in frames:
      continue
    sig_df, bkg_df = frames[scenario]
    records = scenario_records(sig_df, bkg_df, collection_name, veto_mode, reference_name)
    filtered = [record for record in records if record["dt"] < 9998.0]
    if not filtered:
      continue
    dts = [record["dt"] for record in filtered]
    aucs = [record["auc"] for record in filtered]
    plt.plot(dts, aucs, marker="o", linewidth=2, label=scenario)
    drew = True
  if not drew:
    plt.close()
    return
  plt.xlabel("Timing cut")
  plt.ylabel("AUC")
  plt.title(f"AUC scan: {COLLECTION_LABELS[collection_name]} {veto_mode} veto ({reference_name})")
  plt.grid(alpha=0.3)
  plt.legend()
  plt.tight_layout()
  out_name = f"auc_scan_{collection_name}_{veto_mode}_ref_{reference_name}.png"
  plt.savefig(os.path.join(plotdir, out_name))
  plt.savefig(os.path.join(plotdir, out_name.replace(".png", ".pdf")))
  plt.close()


def read_merged_tree(files, tree_name):
  if not files:
    raise RuntimeError("No merged ROOT files found for tree {}".format(tree_name))
  print("[read] tree={} files={}".format(tree_name, len(files)))
  try:
    if len(files) == 1:
      return uproot.open(files[0])[tree_name].arrays(library="ak")
    return uproot.concatenate([f"{name}:{tree_name}" for name in files], library="ak")
  except ValueError as exc:
    message = str(exc)
    if "entries in normal baskets" in message and "don't add up to expected number of entries" in message:
      raise RuntimeError(
          "The merged ROOT file looks inconsistent for uproot reading, likely from `hadd` over uproot-written jagged TTrees. "
          "Use the original merged file directory instead of the hadd output, or pass non-hadded signal/background files."
      ) from exc
    raise


def merged_files(merged_dir, particle, region, process_name):
  return sorted(glob.glob(os.path.join(merged_dir, particle, region, process_name, "*.root")))


def resolve_merged_files(merged_dir, merged_file, particle, region, process_name):
  if merged_file:
    return [merged_file]
  return merged_files(merged_dir, particle, region, process_name)


def timing_arrays(events, scenario):
  arrays = {}
  if "Trackster_Time" in events.fields:
    seed_time = ak.fill_none(
        ak.mean(
            events["Trackster_Time"][
                events["Trackster_isSeed"] &
                (events["Trackster_pt"] > 1.0) &
                (events["Trackster_dr"] < 0.3) &
                (events["Trackster_Time"] > -90.0)
            ],
            axis=1,
        ),
        -99.0,
    )
    arrays["Trackster_TimeWrtSeed"] = ak.flatten(
        abs(events["Trackster_Time"] - ak.broadcast_arrays(seed_time, events["Trackster_Time"])[0])
    )
  if "HGCRecHit_Time" in events.fields and "Trackster_Time" in events.fields:
    seed_time = ak.fill_none(
        ak.mean(
            events["Trackster_Time"][
                events["Trackster_isSeed"] &
                (events["Trackster_pt"] > 1.0) &
                (events["Trackster_dr"] < 0.3) &
                (events["Trackster_Time"] > -90.0)
            ],
            axis=1,
        ),
        -99.0,
    )
    arrays["HGCRecHit_TimeWrtSeed"] = ak.flatten(
        abs(events["HGCRecHit_Time"] - ak.broadcast_arrays(seed_time, events["HGCRecHit_Time"])[0])
    )
  if scenario in ONLINE_SCENARIOS and "LayerCluster_Time" in events.fields and "Trackster_Time" in events.fields:
    seed_time = ak.fill_none(
        ak.mean(
            events["Trackster_Time"][
                events["Trackster_isSeed"] &
                (events["Trackster_pt"] > 1.0) &
                (events["Trackster_dr"] < 0.3) &
                (events["Trackster_Time"] > -90.0)
            ],
            axis=1,
        ),
        -99.0,
    )
    arrays["LayerCluster_TimeWrtSeed"] = ak.flatten(
        abs(events["LayerCluster_Time"] - ak.broadcast_arrays(seed_time, events["LayerCluster_Time"])[0])
    )
  if scenario in OFFLINE_SCENARIOS:
    arrays["Track_TimeWrtPV"] = ak.flatten(abs(events["Track_Time"] - ak.broadcast_arrays(events["PV_Time"], events["Track_Time"])[0]))
    arrays["Trackster_TimeWrtPV"] = ak.flatten(abs(events["Trackster_Time"] - ak.broadcast_arrays(events["PV_Time"], events["Trackster_Time"])[0]))
    arrays["HGCRecHit_TimeWrtPV"] = ak.flatten(abs(events["HGCRecHit_Time"] - ak.broadcast_arrays(events["PV_Time"], events["HGCRecHit_Time"])[0]))
    arrays["Track_TimeWrtSigTrk"] = ak.flatten(abs(events["Track_Time"] - ak.broadcast_arrays(events["sigTrkTime"], events["Track_Time"])[0]))
    arrays["Trackster_TimeWrtSigTrk"] = ak.flatten(abs(events["Trackster_Time"] - ak.broadcast_arrays(events["sigTrkTime"], events["Trackster_Time"])[0]))
    arrays["HGCRecHit_TimeWrtSigTrk"] = ak.flatten(abs(events["HGCRecHit_Time"] - ak.broadcast_arrays(events["sigTrkTime"], events["HGCRecHit_Time"])[0]))
  return {name: ak.to_numpy(values[(values > -90) & (values < 900)]) for name, values in arrays.items()}


def plot_timing_distributions(merged_dir, signal_merged_file, background_merged_file, particle, region, signal_process, background_process, plotdir):
  common_variables = ["Trackster_TimeWrtSeed", "HGCRecHit_TimeWrtSeed"]
  online_only = ["LayerCluster_TimeWrtSeed"]
  offline_only = [
      "Track_TimeWrtPV",
      "Track_TimeWrtSigTrk",
      "Trackster_TimeWrtPV",
      "Trackster_TimeWrtSigTrk",
      "HGCRecHit_TimeWrtPV",
      "HGCRecHit_TimeWrtSigTrk",
  ]

  scenario_payload = {}
  prefix = "Ele" if particle == "electron" else "Pho"
  matched_branch = "matchedToGenEle" if particle == "electron" else "matchedToGenPho"

  for scenario in SCENARIOS:
    scenario_payload[scenario] = {}
    for process_name, is_signal in [(signal_process, True), (background_process, False)]:
      merged_file = signal_merged_file if is_signal else background_merged_file
      files = resolve_merged_files(merged_dir, merged_file, particle, region, process_name)
      events = read_merged_tree(files, scenario)
      if is_signal:
        events = events[events[matched_branch] == 1]
      else:
        events = events[events[matched_branch] != 1]
      scenario_payload[scenario][process_name] = timing_arrays(events, scenario)

  variable_iter = common_variables + online_only + offline_only
  if tqdm is not None:
    variable_iter = tqdm(variable_iter, desc="Timing plots", unit="plot")
  for variable in variable_iter:
    valid_scenarios = []
    for scenario in SCENARIOS:
      if variable in scenario_payload[scenario].get(signal_process, {}):
        valid_scenarios.append(scenario)
    if not valid_scenarios:
      continue

    for process_name, label in [(signal_process, "signal"), (background_process, "background")]:
      plt.figure(figsize=(8, 6))
      drew = False
      for scenario in valid_scenarios:
        values = scenario_payload[scenario][process_name].get(variable)
        if values is None or len(values) == 0:
          continue
        plt.hist(values, bins=np.linspace(0.0, 0.25, 60), density=True, histtype="step", linewidth=2, label=scenario)
        drew = True
      if not drew:
        plt.close()
        continue
      plt.xlabel(variable)
      plt.ylabel("Density")
      plt.title(f"{variable} ({label})")
      plt.yscale("log")
      plt.grid(alpha=0.3)
      plt.legend()
      plt.tight_layout()
      out_name = f"timing_{variable}_{label}.png"
      plt.savefig(os.path.join(plotdir, out_name))
      plt.savefig(os.path.join(plotdir, out_name.replace(".png", ".pdf")))
      plt.close()
      print("[plot] {}".format(out_name))


def main():
  parser = argparse.ArgumentParser(description="Summarize HGCal timing/isolation performance for online/offline v4/v5 comparison")
  parser.add_argument("--isoDir", required=True, type=str)
  parser.add_argument("--mergedDir", type=str)
  parser.add_argument("--mergedFile", type=str)
  parser.add_argument("--signal-merged-file", type=str)
  parser.add_argument("--background-merged-file", type=str)
  parser.add_argument("--plotDir", required=True, type=str)
  parser.add_argument("--particle", required=True, choices=["electron", "photon"])
  parser.add_argument("--region", required=True, type=str)
  parser.add_argument("--signal", required=True, type=str)
  parser.add_argument("--background", required=True, type=str)
  args = parser.parse_args()

  if not args.mergedDir and not args.mergedFile:
    if not args.signal_merged_file or not args.background_merged_file:
      raise RuntimeError("Provide --mergedDir, --mergedFile, or both --signal-merged-file and --background-merged-file")

  check_dir(args.plotDir)

  frames = {}
  scenario_iter = SCENARIOS
  if tqdm is not None:
    scenario_iter = tqdm(SCENARIOS, desc="Loading scenarios", unit="scenario")
  for scenario in scenario_iter:
    print("[summary] scenario={}".format(scenario))
    sig_df = load_frame(args.isoDir, args.particle, args.region, scenario, args.signal)
    bkg_df = load_frame(args.isoDir, args.particle, args.region, scenario, args.background)
    frames[scenario] = apply_reweighting(sig_df, bkg_df)

  summary = {}
  roc_tasks = [(collection_name, veto_mode) for collection_name in ["trackster", "HGCRecHit"] for veto_mode in ["cone", "seed", "cluster"]]
  roc_iter = roc_tasks
  if tqdm is not None:
    roc_iter = tqdm(roc_tasks, desc="Common ROC scans", unit="scan")
  for collection_name, veto_mode in roc_iter:
      key = f"{collection_name}_{veto_mode}_seed"
      summary[key] = plot_roc_by_scenario(frames, collection_name, veto_mode, "seed", args.plotDir, include_default_online=True)
      plot_auc_scan(frames, collection_name, veto_mode, "seed", args.plotDir)
      print("[roc] {} {} ref=seed".format(collection_name, veto_mode))

  online_iter = ["cone", "seed", "cluster"]
  if tqdm is not None:
    online_iter = tqdm(online_iter, desc="Online layerCluster ROC", unit="scan")
  for veto_mode in online_iter:
    key = f"layerCluster_{veto_mode}_seed"
    summary[key] = plot_roc_by_scenario(
        {scenario: frames[scenario] for scenario in SCENARIOS if scenario in ONLINE_SCENARIOS},
        "layerCluster",
        veto_mode,
        "seed",
        args.plotDir,
    )
    plot_auc_scan({scenario: frames[scenario] for scenario in SCENARIOS if scenario in ONLINE_SCENARIOS}, "layerCluster", veto_mode, "seed", args.plotDir)
    print("[roc] layerCluster {} ref=seed".format(veto_mode))

  offline_frames = {scenario: frames[scenario] for scenario in SCENARIOS if scenario in OFFLINE_SCENARIOS}
  offline_tasks = [
      (reference_name, collection_name, veto_mode)
      for reference_name in OFFLINE_ONLY_REFERENCES
      for collection_name in ["trackster", "HGCRecHit"]
      for veto_mode in ["cone", "seed", "cluster"]
  ]
  offline_iter = offline_tasks
  if tqdm is not None:
    offline_iter = tqdm(offline_tasks, desc="Offline-only ROC scans", unit="scan")
  for reference_name, collection_name, veto_mode in offline_iter:
        key = f"{collection_name}_{veto_mode}_{reference_name}"
        summary[key] = plot_roc_by_scenario(offline_frames, collection_name, veto_mode, reference_name, args.plotDir)
        plot_auc_scan(offline_frames, collection_name, veto_mode, reference_name, args.plotDir)
        print("[roc] {} {} ref={}".format(collection_name, veto_mode, reference_name))

  with open(os.path.join(args.plotDir, "roc_summary.json"), "w") as handle:
    json.dump(summary, handle, indent=2, sort_keys=True)
  print("[saved] {}".format(os.path.join(args.plotDir, "roc_summary.json")))

  signal_merged_file = args.signal_merged_file if args.signal_merged_file else args.mergedFile
  background_merged_file = args.background_merged_file if args.background_merged_file else args.mergedFile
  plot_timing_distributions(
      args.mergedDir,
      signal_merged_file,
      background_merged_file,
      args.particle,
      args.region,
      args.signal,
      args.background,
      args.plotDir,
  )


if __name__ == "__main__":
  main()
