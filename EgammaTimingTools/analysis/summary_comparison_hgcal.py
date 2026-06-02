#!/usr/bin/env python3
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
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
  import pyarrow.parquet as pq
except ImportError:
  pq = None

try:
  from tqdm import tqdm
except ImportError:
  tqdm = None


SCENARIOS = ["online_v4", "online_v5", "offline_v4", "offline_v5"]
ONLINE_SCENARIOS = {"online_v4", "online_v5"}
OFFLINE_SCENARIOS = {"offline_v4", "offline_v5"}
COMMON_REFERENCES = ["seed"]
OFFLINE_ONLY_REFERENCES = ["pv", "sigTrk"]
DT_CUTS = [round(0.01 * (index + 1), 3) for index in range(20)] + [9999.0]
COLLECTION_LABELS = {
    "trackster": "Trackster",
    "HGCRecHit": "HGCRecHit",
    "layerCluster": "LayerCluster",
}


def check_dir(path):
  if not os.path.exists(path):
    os.system("mkdir -p {}".format(path))


def parquet_sources(iso_dir, particle, region, scenario, process_name):
  single_file = os.path.join(iso_dir, particle, region, scenario, f"{process_name}_iso.parquet")
  if os.path.exists(single_file):
    return [single_file]
  part_pattern = os.path.join(iso_dir, particle, region, scenario, f"{process_name}_iso_parts", "*.parquet")
  part_files = sorted(glob.glob(part_pattern))
  if part_files:
    return part_files
  raise RuntimeError("No parquet sources found for scenario={} process={}".format(scenario, process_name))


def apply_min_pt(frame, min_pt):
  if min_pt is None or min_pt <= 0.0 or "eg-pt" not in frame.columns:
    return frame
  return frame[frame["eg-pt"] >= min_pt].copy()


def load_frame(iso_dir, particle, region, scenario, process_name, min_pt):
  sources = parquet_sources(iso_dir, particle, region, scenario, process_name)
  print("[load] {} parquet file(s) for scenario={} process={}".format(len(sources), scenario, process_name))
  if len(sources) == 1:
    return apply_min_pt(pd.read_parquet(sources[0]), min_pt)
  return apply_min_pt(pd.concat((pd.read_parquet(source) for source in sources), ignore_index=True, copy=False), min_pt)


def load_scenario_frames(iso_dir, particle, region, scenario, signal_process, background_process, min_pt):
  print("[summary] scenario={}".format(scenario))
  sig_df = load_frame(iso_dir, particle, region, scenario, signal_process, min_pt)
  bkg_df = load_frame(iso_dir, particle, region, scenario, background_process, min_pt)
  return scenario, apply_reweighting(sig_df, bkg_df, min_pt)


def compute_weights(sig_df, bkg_df, pt_bins, eta_bins):
  sig_hist, _, _ = np.histogram2d(sig_df["eg-pt"], sig_df["eg-eta"], bins=[pt_bins, eta_bins])
  bkg_hist, _, _ = np.histogram2d(bkg_df["eg-pt"], bkg_df["eg-eta"], bins=[pt_bins, eta_bins])
  weights_2d = np.divide(bkg_hist, sig_hist, out=np.zeros_like(bkg_hist), where=sig_hist > 0)

  pt_idx = np.clip(np.digitize(sig_df["eg-pt"], pt_bins) - 1, 0, len(pt_bins) - 2)
  eta_idx = np.clip(np.digitize(sig_df["eg-eta"], eta_bins) - 1, 0, len(eta_bins) - 2)
  return weights_2d[pt_idx, eta_idx]


def apply_reweighting(sig_df, bkg_df, min_pt):
  pt_low = max(15.0, min_pt if min_pt is not None else 15.0)
  pt_bins = np.linspace(pt_low, 200, 20)
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
  sig_mask = np.isfinite(sig_vals) & np.isfinite(sig_weights) & (sig_vals > 0.0)
  bkg_mask = np.isfinite(bkg_vals) & np.isfinite(bkg_weights) & (bkg_vals > 0.0)
  sig_vals = sig_vals[sig_mask]
  bkg_vals = bkg_vals[bkg_mask]
  sig_weights = sig_weights[sig_mask]
  bkg_weights = bkg_weights[bkg_mask]
  if sig_vals.size == 0 or bkg_vals.size == 0:
    return None
  if np.sum(sig_weights) <= 0.0 or np.sum(bkg_weights) <= 0.0:
    return None

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
    roc = compute_roc(sig_df, bkg_df, column)
    if roc is None:
      continue
    fpr, tpr, _, roc_auc = roc
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
        roc = compute_roc(frames[scenario][0], frames[scenario][1], branch)
        if roc is not None:
          fpr, tpr, _, _ = roc
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


def parquet_schema_names(source):
  if pq is not None:
    return set(pq.read_schema(source).names)
  return set(pd.read_parquet(source).columns)


def parquet_has_column(sources, column):
  if not sources:
    return False
  return column in parquet_schema_names(sources[0])


def read_parquet_column(source, column):
  if not parquet_has_column([source], column):
    return np.asarray([], dtype=float)
  values = pd.read_parquet(source, columns=[column])[column].to_numpy()
  values = values[np.isfinite(values)]
  return values


def read_parquet_columns(source, columns):
  available = parquet_schema_names(source)
  if any(column not in available for column in columns):
    return pd.DataFrame(columns=columns)
  return pd.read_parquet(source, columns=columns)


def column_values_with_min_pt(source, column, min_pt):
  columns = [column]
  if min_pt is not None and min_pt > 0.0:
    columns.append("eg-pt")
  df = read_parquet_columns(source, columns)
  if df.empty or column not in df.columns:
    return np.asarray([], dtype=float)
  values = df[column].to_numpy()
  mask = np.isfinite(values)
  if min_pt is not None and min_pt > 0.0:
    if "eg-pt" not in df.columns:
      return np.asarray([], dtype=float)
    pt_values = df["eg-pt"].to_numpy()
    mask = mask & np.isfinite(pt_values) & (pt_values >= min_pt)
  return values[mask]


def column_quantile_range(sources, column, min_pt, max_values=200000):
  samples = []
  total = 0
  for source in sources:
    values = column_values_with_min_pt(source, column, min_pt)
    if values.size == 0:
      continue
    if total + values.size > max_values:
      keep = max(max_values - total, 0)
      if keep == 0:
        break
      values = values[:keep]
    samples.append(values)
    total += values.size
    if total >= max_values:
      break
  if not samples:
    return None
  combined = np.concatenate(samples)
  combined = combined[np.isfinite(combined)]
  if combined.size == 0:
    return None
  high = float(np.quantile(combined, 0.995))
  if not np.isfinite(high) or high <= 0.0:
    high = float(np.max(combined))
  if not np.isfinite(high) or high <= 0.0:
    high = 1.0
  return 0.0, high


def column_max(sources, column, min_pt):
  high = None
  for source in sources:
    values = column_values_with_min_pt(source, column, min_pt)
    if values.size == 0:
      continue
    source_high = float(np.max(values))
    if high is None or source_high > high:
      high = source_high
  if high is None or not np.isfinite(high) or high <= 0.0:
    return None
  return high


def histogram_column_sources(sources, column, bins, min_pt):
  counts = np.zeros(len(bins) - 1, dtype=float)
  n_values = 0
  for source in sources:
    values = column_values_with_min_pt(source, column, min_pt)
    if values.size == 0:
      continue
    hist, _ = np.histogram(values, bins=bins)
    counts += hist
    n_values += values.size
  return counts, n_values


def density_from_counts(counts, bins):
  total = np.sum(counts)
  widths = np.diff(bins)
  if total <= 0:
    return counts
  return counts / (total * widths)


def plot_step_counts(bins, counts, label):
  density = density_from_counts(counts, bins)
  x_values = np.repeat(bins, 2)[1:-1]
  y_values = np.repeat(density, 2)
  plt.plot(x_values, y_values, linewidth=2, label=label)


def plot_isolation_distribution_from_parquet(
    iso_dir,
    particle,
    region,
    signal_process,
    background_process,
    variable,
    scenario,
    collection_name,
    veto_mode,
    reference_name,
    tag,
    plotdir,
    min_pt,
):
  sig_sources = parquet_sources(iso_dir, particle, region, scenario, signal_process)
  bkg_sources = parquet_sources(iso_dir, particle, region, scenario, background_process)
  if not parquet_has_column(sig_sources, variable) or not parquet_has_column(bkg_sources, variable):
    return

  sig_range = column_quantile_range(sig_sources, variable, min_pt)
  bkg_range = column_quantile_range(bkg_sources, variable, min_pt)
  if sig_range is None or bkg_range is None:
    return
  high = max(sig_range[1], bkg_range[1])
  bins = np.linspace(0.0, high, 80)
  sig_counts, sig_entries = histogram_column_sources(sig_sources, variable, bins, min_pt)
  bkg_counts, bkg_entries = histogram_column_sources(bkg_sources, variable, bins, min_pt)
  if sig_entries == 0 or bkg_entries == 0:
    return

  plt.figure(figsize=(8, 6))
  plot_step_counts(bins, sig_counts, "signal")
  plot_step_counts(bins, bkg_counts, "background")
  plt.xlabel(variable)
  plt.ylabel("Density")
  plt.yscale("log")
  plt.title(f"{scenario}: {COLLECTION_LABELS[collection_name]} {veto_mode} veto ({reference_name}, {tag})")
  plt.grid(alpha=0.3)
  plt.legend()
  plt.tight_layout()
  safe_tag = tag.replace(" ", "_")
  out_name = f"iso_dist_{scenario}_{collection_name}_{veto_mode}_ref_{reference_name}_{safe_tag}.png"
  plt.savefig(os.path.join(plotdir, out_name))
  plt.savefig(os.path.join(plotdir, out_name.replace(".png", ".pdf")))
  plt.close()
  print("[plot] {}".format(out_name))


def accumulate_pt_eta_hist(sources, pt_bins, eta_bins, min_pt):
  counts = np.zeros((len(pt_bins) - 1, len(eta_bins) - 1), dtype=float)
  for source in sources:
    df = read_parquet_columns(source, ["eg-pt", "eg-eta"])
    if df.empty:
      continue
    pt_values = df["eg-pt"].to_numpy()
    eta_values = df["eg-eta"].to_numpy()
    mask = np.isfinite(pt_values) & np.isfinite(eta_values)
    if min_pt is not None and min_pt > 0.0:
      mask = mask & (pt_values >= min_pt)
    if not np.any(mask):
      continue
    hist, _, _ = np.histogram2d(pt_values[mask], eta_values[mask], bins=[pt_bins, eta_bins])
    counts += hist
  return counts


def streaming_weight_map(sig_sources, bkg_sources, min_pt):
  pt_low = max(15.0, min_pt if min_pt is not None else 15.0)
  pt_bins = np.linspace(pt_low, 200, 20)
  eta_bins = np.linspace(-3.2, 3.2, 33)
  sig_hist = accumulate_pt_eta_hist(sig_sources, pt_bins, eta_bins, min_pt)
  bkg_hist = accumulate_pt_eta_hist(bkg_sources, pt_bins, eta_bins, min_pt)
  weights_2d = np.divide(bkg_hist, sig_hist, out=np.zeros_like(bkg_hist), where=sig_hist > 0)
  return weights_2d, pt_bins, eta_bins


def lookup_streaming_weights(pt_values, eta_values, weights_2d, pt_bins, eta_bins):
  pt_idx = np.clip(np.digitize(pt_values, pt_bins) - 1, 0, len(pt_bins) - 2)
  eta_idx = np.clip(np.digitize(eta_values, eta_bins) - 1, 0, len(eta_bins) - 2)
  return weights_2d[pt_idx, eta_idx]


def weighted_score_histogram(sources, column, bins, weight_context=None, min_pt=None):
  counts = np.zeros(len(bins) - 1, dtype=float)
  n_values = 0
  for source in sources:
    columns = [column]
    if weight_context is not None:
      columns += ["eg-pt", "eg-eta"]
    elif min_pt is not None and min_pt > 0.0:
      columns.append("eg-pt")
    df = read_parquet_columns(source, columns)
    if df.empty or column not in df.columns:
      continue
    values = df[column].to_numpy()
    mask = np.isfinite(values) & (values > 0.0)
    if min_pt is not None and min_pt > 0.0:
      if "eg-pt" not in df.columns:
        continue
      pt_values_for_cut = df["eg-pt"].to_numpy()
      mask = mask & np.isfinite(pt_values_for_cut) & (pt_values_for_cut >= min_pt)
    weights = None
    if weight_context is not None:
      weights_2d, pt_bins, eta_bins = weight_context
      pt_values = df["eg-pt"].to_numpy()
      eta_values = df["eg-eta"].to_numpy()
      mask = mask & np.isfinite(pt_values) & np.isfinite(eta_values)
      if np.any(mask):
        weights = lookup_streaming_weights(pt_values[mask], eta_values[mask], weights_2d, pt_bins, eta_bins)
    if not np.any(mask):
      continue
    hist, _ = np.histogram(values[mask], bins=bins, weights=weights)
    counts += hist
    n_values += int(np.sum(mask))
  return counts, n_values


def binned_roc(sig_counts, bkg_counts):
  sig_total = np.sum(sig_counts)
  bkg_total = np.sum(bkg_counts)
  if sig_total <= 0 or bkg_total <= 0:
    return None
  tpr = np.concatenate([[0.0], np.cumsum(sig_counts) / sig_total])
  fpr = np.concatenate([[0.0], np.cumsum(bkg_counts) / bkg_total])
  return fpr, tpr, auc(fpr, tpr)


def streaming_contexts(iso_dir, particle, region, scenarios, signal_process, background_process, min_pt):
  contexts = {}
  for scenario in scenarios:
    sig_sources = parquet_sources(iso_dir, particle, region, scenario, signal_process)
    bkg_sources = parquet_sources(iso_dir, particle, region, scenario, background_process)
    print(
        "[stream] scenario={} signal_parts={} background_parts={}".format(
            scenario,
            len(sig_sources),
            len(bkg_sources),
        )
    )
    contexts[scenario] = {
        "sig_sources": sig_sources,
        "bkg_sources": bkg_sources,
        "weight_context": streaming_weight_map(sig_sources, bkg_sources, min_pt),
        "records_cache": {},
        "min_pt": min_pt,
    }
  return contexts


def streaming_candidate_columns(context, collection_name, veto_mode, reference_name):
  prefix = iso_family_prefix(collection_name, veto_mode, reference_name)
  sig_columns = parquet_schema_names(context["sig_sources"][0])
  bkg_columns = parquet_schema_names(context["bkg_sources"][0])
  columns = sorted((sig_columns & bkg_columns), key=lambda name: parse_dt_from_column(name) or -1.0)
  return [column for column in columns if column.startswith(prefix) and parse_dt_from_column(column) is not None]


def streaming_roc_record(context, column, roc_bins):
  high = max(
      column_max(context["sig_sources"], column, context["min_pt"]) or 0.0,
      column_max(context["bkg_sources"], column, context["min_pt"]) or 0.0,
  )
  if high <= 0.0:
    return None
  bins = np.linspace(0.0, high, roc_bins + 1)
  sig_counts, sig_entries = weighted_score_histogram(
      context["sig_sources"],
      column,
      bins,
      context["weight_context"],
      context["min_pt"],
  )
  bkg_counts, bkg_entries = weighted_score_histogram(
      context["bkg_sources"],
      column,
      bins,
      min_pt=context["min_pt"],
  )
  if sig_entries == 0 or bkg_entries == 0:
    return None
  roc = binned_roc(sig_counts, bkg_counts)
  if roc is None:
    return None
  fpr, tpr, roc_auc = roc
  return fpr, tpr, roc_auc


def streaming_scenario_records(context, collection_name, veto_mode, reference_name, roc_bins):
  records = []
  for column in streaming_candidate_columns(context, collection_name, veto_mode, reference_name):
    dt_value = parse_dt_from_column(column)
    roc = streaming_roc_record(context, column, roc_bins)
    if roc is None:
      continue
    fpr, tpr, roc_auc = roc
    records.append({
        "column": column,
        "dt": dt_value,
        "fpr": fpr,
        "tpr": tpr,
        "auc": roc_auc,
    })
  return sorted(records, key=lambda item: item["dt"])


def cached_streaming_scenario_records(context, collection_name, veto_mode, reference_name, roc_bins):
  cache_key = (collection_name, veto_mode, reference_name, roc_bins)
  if cache_key not in context["records_cache"]:
    context["records_cache"][cache_key] = streaming_scenario_records(
        context,
        collection_name,
        veto_mode,
        reference_name,
        roc_bins,
    )
  return context["records_cache"][cache_key]


def streaming_default_hlt_branch(context):
  sig_columns = parquet_schema_names(context["sig_sources"][0])
  bkg_columns = parquet_schema_names(context["bkg_sources"][0])
  for branch in ["eg-hgcalPFIsol_default", "eg-trkIsol_default", "eg-trkIsolV6_default", "eg-trkIsolV72_default"]:
    if branch in sig_columns and branch in bkg_columns:
      return branch
  return None


def plot_roc_by_scenario_streaming(contexts, collection_name, veto_mode, reference_name, plotdir, roc_bins, include_default_online=False):
  available = {}
  summary = {}
  for scenario, context in contexts.items():
    records = cached_streaming_scenario_records(context, collection_name, veto_mode, reference_name, roc_bins)
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
    plt.plot(best["tpr"], 1.0 - best["fpr"], color=color, linewidth=2, label=f"{scenario} best ({best['dt']:.3f})")
    if default is not None:
      plt.plot(default["tpr"], 1.0 - default["fpr"], color=color, linewidth=1.5, linestyle="--", label=f"{scenario} no cut")

    if include_default_online and scenario in ONLINE_SCENARIOS:
      branch = streaming_default_hlt_branch(contexts[scenario])
      if branch is not None:
        roc = streaming_roc_record(contexts[scenario], branch, roc_bins)
        if roc is not None:
          fpr, tpr, _ = roc
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


def plot_auc_scan_streaming(contexts, collection_name, veto_mode, reference_name, plotdir, roc_bins):
  plt.figure(figsize=(8, 6))
  drew = False
  for scenario in SCENARIOS:
    if scenario not in contexts:
      continue
    records = cached_streaming_scenario_records(contexts[scenario], collection_name, veto_mode, reference_name, roc_bins)
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


def plot_isolation_distributions_for_family_streaming(
    contexts,
    collection_name,
    veto_mode,
    reference_name,
    plotdir,
    roc_bins,
    iso_dir,
    particle,
    region,
    signal_process,
    background_process,
):
  for scenario, context in contexts.items():
    records = cached_streaming_scenario_records(context, collection_name, veto_mode, reference_name, roc_bins)
    best = best_record(records)
    default = default_record(records)
    if best is not None:
      plot_isolation_distribution_from_parquet(
          iso_dir,
          particle,
          region,
          signal_process,
          background_process,
          best["column"],
          scenario,
          collection_name,
          veto_mode,
          reference_name,
          f"best_dt_{best['dt']:.3f}",
          plotdir,
          context["min_pt"],
      )
    if default is not None and (best is None or default["column"] != best["column"]):
      plot_isolation_distribution_from_parquet(
          iso_dir,
          particle,
          region,
          signal_process,
          background_process,
          default["column"],
          scenario,
          collection_name,
          veto_mode,
          reference_name,
          "no_cut",
          plotdir,
          context["min_pt"],
      )


def run_streaming_roc_summary(args, scenarios):
  contexts = streaming_contexts(args.isoDir, args.particle, args.region, scenarios, args.signal, args.background, args.min_pt)
  summary = {}
  iso_distribution_args = (args.isoDir, args.particle, args.region, args.signal, args.background)

  roc_tasks = [(collection_name, veto_mode) for collection_name in ["trackster", "HGCRecHit"] for veto_mode in ["cone", "seed", "cluster"]]
  roc_iter = roc_tasks
  if tqdm is not None:
    roc_iter = tqdm(roc_tasks, desc="Common ROC scans", unit="scan")
  for collection_name, veto_mode in roc_iter:
    key = f"{collection_name}_{veto_mode}_seed"
    summary[key] = plot_roc_by_scenario_streaming(
        contexts,
        collection_name,
        veto_mode,
        "seed",
        args.plotDir,
        args.roc_bins,
        include_default_online=True,
    )
    plot_auc_scan_streaming(contexts, collection_name, veto_mode, "seed", args.plotDir, args.roc_bins)
    if not args.skip_iso_distributions:
      plot_isolation_distributions_for_family_streaming(
          contexts, collection_name, veto_mode, "seed", args.plotDir, args.roc_bins, *iso_distribution_args
      )
    print("[roc] {} {} ref=seed".format(collection_name, veto_mode))

  online_contexts = {scenario: contexts[scenario] for scenario in scenarios if scenario in ONLINE_SCENARIOS and scenario in contexts}
  online_iter = ["cone", "seed", "cluster"]
  if tqdm is not None:
    online_iter = tqdm(online_iter, desc="Online layerCluster ROC", unit="scan")
  for veto_mode in online_iter:
    key = f"layerCluster_{veto_mode}_seed"
    summary[key] = plot_roc_by_scenario_streaming(
        online_contexts,
        "layerCluster",
        veto_mode,
        "seed",
        args.plotDir,
        args.roc_bins,
    )
    plot_auc_scan_streaming(online_contexts, "layerCluster", veto_mode, "seed", args.plotDir, args.roc_bins)
    if not args.skip_iso_distributions:
      plot_isolation_distributions_for_family_streaming(
          online_contexts, "layerCluster", veto_mode, "seed", args.plotDir, args.roc_bins, *iso_distribution_args
      )
    print("[roc] layerCluster {} ref=seed".format(veto_mode))

  online_layer_seed_tasks = [
      (collection_name, veto_mode)
      for collection_name in ["trackster", "HGCRecHit", "layerCluster"]
      for veto_mode in ["cone", "seed", "cluster"]
  ]
  online_layer_seed_iter = online_layer_seed_tasks
  if tqdm is not None:
    online_layer_seed_iter = tqdm(online_layer_seed_tasks, desc="Online layerClusterSeed ROC", unit="scan")
  for collection_name, veto_mode in online_layer_seed_iter:
    key = f"{collection_name}_{veto_mode}_layerClusterSeed"
    summary[key] = plot_roc_by_scenario_streaming(
        online_contexts,
        collection_name,
        veto_mode,
        "layerClusterSeed",
        args.plotDir,
        args.roc_bins,
        include_default_online=(collection_name != "layerCluster"),
    )
    plot_auc_scan_streaming(online_contexts, collection_name, veto_mode, "layerClusterSeed", args.plotDir, args.roc_bins)
    if not args.skip_iso_distributions:
      plot_isolation_distributions_for_family_streaming(
          online_contexts, collection_name, veto_mode, "layerClusterSeed", args.plotDir, args.roc_bins, *iso_distribution_args
      )
    print("[roc] {} {} ref=layerClusterSeed".format(collection_name, veto_mode))

  online_cluster_ref_tasks = [
      (reference_name, collection_name, veto_mode)
      for reference_name in ["tracksterCluster", "layerClusterCluster"]
      for collection_name in ["trackster", "HGCRecHit", "layerCluster"]
      for veto_mode in ["cone", "seed", "cluster"]
  ]
  online_cluster_ref_iter = online_cluster_ref_tasks
  if tqdm is not None:
    online_cluster_ref_iter = tqdm(online_cluster_ref_tasks, desc="Online inCluster ROC", unit="scan")
  for reference_name, collection_name, veto_mode in online_cluster_ref_iter:
    key = f"{collection_name}_{veto_mode}_{reference_name}"
    summary[key] = plot_roc_by_scenario_streaming(
        online_contexts,
        collection_name,
        veto_mode,
        reference_name,
        args.plotDir,
        args.roc_bins,
        include_default_online=(collection_name != "layerCluster"),
    )
    plot_auc_scan_streaming(online_contexts, collection_name, veto_mode, reference_name, args.plotDir, args.roc_bins)
    if not args.skip_iso_distributions:
      plot_isolation_distributions_for_family_streaming(
          online_contexts, collection_name, veto_mode, reference_name, args.plotDir, args.roc_bins, *iso_distribution_args
      )
    print("[roc] {} {} ref={}".format(collection_name, veto_mode, reference_name))

  offline_contexts = {scenario: contexts[scenario] for scenario in scenarios if scenario in OFFLINE_SCENARIOS and scenario in contexts}
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
    summary[key] = plot_roc_by_scenario_streaming(
        offline_contexts,
        collection_name,
        veto_mode,
        reference_name,
        args.plotDir,
        args.roc_bins,
    )
    plot_auc_scan_streaming(offline_contexts, collection_name, veto_mode, reference_name, args.plotDir, args.roc_bins)
    if not args.skip_iso_distributions:
      plot_isolation_distributions_for_family_streaming(
          offline_contexts, collection_name, veto_mode, reference_name, args.plotDir, args.roc_bins, *iso_distribution_args
      )
    print("[roc] {} {} ref={}".format(collection_name, veto_mode, reference_name))

  with open(os.path.join(args.plotDir, "roc_summary.json"), "w") as handle:
    json.dump(summary, handle, indent=2, sort_keys=True)
  print("[saved] {}".format(os.path.join(args.plotDir, "roc_summary.json")))


def plot_isolation_distributions_for_family(
    frames,
    collection_name,
    veto_mode,
    reference_name,
    plotdir,
    iso_dir,
    particle,
    region,
    signal_process,
    background_process,
    min_pt,
):
  for scenario, (sig_df, bkg_df) in frames.items():
    records = scenario_records(sig_df, bkg_df, collection_name, veto_mode, reference_name)
    best = best_record(records)
    default = default_record(records)
    if best is not None:
      plot_isolation_distribution_from_parquet(
          iso_dir,
          particle,
          region,
          signal_process,
          background_process,
          best["column"],
          scenario,
          collection_name,
          veto_mode,
          reference_name,
          f"best_dt_{best['dt']:.3f}",
          plotdir,
          min_pt,
      )
    if default is not None and (best is None or default["column"] != best["column"]):
      plot_isolation_distribution_from_parquet(
          iso_dir,
          particle,
          region,
          signal_process,
          background_process,
          default["column"],
          scenario,
          collection_name,
          veto_mode,
          reference_name,
          "no_cut",
          plotdir,
          min_pt,
      )


def iso_column_name(collection_name, veto_mode, reference_name, dt_value):
  return "{}{:.3f}".format(iso_family_prefix(collection_name, veto_mode, reference_name), dt_value)


def histogram_family_tasks(scenarios):
  tasks = []
  for scenario in scenarios:
    for collection_name in ["trackster", "HGCRecHit"]:
      for veto_mode in ["cone", "seed", "cluster"]:
        tasks.append((scenario, collection_name, veto_mode, "seed"))
    if scenario in ONLINE_SCENARIOS:
      for veto_mode in ["cone", "seed", "cluster"]:
        tasks.append((scenario, "layerCluster", veto_mode, "seed"))
      for reference_name in ["layerClusterSeed", "tracksterCluster", "layerClusterCluster"]:
        for collection_name in ["trackster", "HGCRecHit", "layerCluster"]:
          for veto_mode in ["cone", "seed", "cluster"]:
            tasks.append((scenario, collection_name, veto_mode, reference_name))
    if scenario in OFFLINE_SCENARIOS:
      for reference_name in OFFLINE_ONLY_REFERENCES:
        for collection_name in ["trackster", "HGCRecHit"]:
          for veto_mode in ["cone", "seed", "cluster"]:
            tasks.append((scenario, collection_name, veto_mode, reference_name))
  return tasks


def plot_histogram_only_iso_distributions(
    iso_dir,
    particle,
    region,
    signal_process,
    background_process,
    scenarios,
    plotdir,
    dt_values,
    min_pt,
):
  tasks = histogram_family_tasks(scenarios)
  task_iter = tasks
  if tqdm is not None:
    task_iter = tqdm(tasks, desc="Isolation histograms", unit="family")
  for scenario, collection_name, veto_mode, reference_name in task_iter:
    for dt_value in dt_values:
      variable = iso_column_name(collection_name, veto_mode, reference_name, dt_value)
      tag = "no_cut" if abs(dt_value - 9999.0) < 1e-6 else f"dt_{dt_value:.3f}"
      plot_isolation_distribution_from_parquet(
          iso_dir,
          particle,
          region,
          signal_process,
          background_process,
          variable,
          scenario,
          collection_name,
          veto_mode,
          reference_name,
          tag,
          plotdir,
          min_pt,
      )


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
  seed_mask = None
  seed_time = None
  if "Trackster_Time" in events.fields:
    seed_mask = (
        events["Trackster_isSeed"] &
        (events["Trackster_pt"] > 1.0) &
        (events["Trackster_dr"] < 0.3) &
        (events["Trackster_Time"] > -90.0)
    )
    seed_time = ak.fill_none(ak.mean(events["Trackster_Time"][seed_mask], axis=1), -99.0)
    arrays["SeedTime"] = ak.to_numpy(seed_time[seed_time > -90.0])
    arrays["SeedMultiplicity"] = ak.to_numpy(ak.sum(seed_mask, axis=1))
    trackster_cluster_mask = (
        events["Trackster_inCluster"] &
        (events["Trackster_pt"] > 1.0) &
        (events["Trackster_dr"] < 0.3) &
        (events["Trackster_Time"] > -90.0)
    )
    trackster_cluster_weighted_sum = ak.sum(events["Trackster_Time"][trackster_cluster_mask] * events["Trackster_pt"][trackster_cluster_mask], axis=1)
    trackster_cluster_pt_sum = ak.sum(events["Trackster_pt"][trackster_cluster_mask], axis=1)
    trackster_cluster_time = ak.fill_none(ak.where(trackster_cluster_pt_sum > 0, trackster_cluster_weighted_sum / trackster_cluster_pt_sum, -99.0), -99.0)
    arrays["TracksterClusterTime"] = ak.to_numpy(trackster_cluster_time[trackster_cluster_time > -90.0])
    arrays["TracksterClusterMultiplicity"] = ak.to_numpy(ak.sum(trackster_cluster_mask, axis=1))
    arrays["Trackster_Time"] = ak.to_numpy(ak.flatten(events["Trackster_Time"][events["Trackster_Time"] > -90.0]))
    arrays["Trackster_TimeWrtSeed"] = ak.flatten(
        abs(events["Trackster_Time"] - ak.broadcast_arrays(seed_time, events["Trackster_Time"])[0])
    )
    arrays["Trackster_TimeWrtTracksterCluster"] = ak.flatten(
        abs(events["Trackster_Time"] - ak.broadcast_arrays(trackster_cluster_time, events["Trackster_Time"])[0])
    )
  if "HGCRecHit_Time" in events.fields and "Trackster_Time" in events.fields:
    arrays["HGCRecHit_Time"] = ak.to_numpy(ak.flatten(events["HGCRecHit_Time"][events["HGCRecHit_Time"] > -90.0]))
    arrays["HGCRecHit_TimeWrtSeed"] = ak.flatten(
        abs(events["HGCRecHit_Time"] - ak.broadcast_arrays(seed_time, events["HGCRecHit_Time"])[0])
    )
    arrays["HGCRecHit_TimeWrtTracksterCluster"] = ak.flatten(
        abs(events["HGCRecHit_Time"] - ak.broadcast_arrays(trackster_cluster_time, events["HGCRecHit_Time"])[0])
    )
  if scenario in ONLINE_SCENARIOS and "LayerCluster_Time" in events.fields and "Trackster_Time" in events.fields:
    layercluster_seed_mask = (
        events["LayerCluster_isSeed"] &
        (events["LayerCluster_energy"] > 0.5) &
        (events["LayerCluster_dr"] < 0.3) &
        (events["LayerCluster_Time"] > -90.0)
    )
    weighted_sum = ak.sum(events["LayerCluster_Time"][layercluster_seed_mask] * events["LayerCluster_energy"][layercluster_seed_mask], axis=1)
    energy_sum = ak.sum(events["LayerCluster_energy"][layercluster_seed_mask], axis=1)
    layercluster_seed_time = ak.fill_none(ak.where(energy_sum > 0, weighted_sum / energy_sum, -99.0), -99.0)
    arrays["LayerClusterSeedTime"] = ak.to_numpy(layercluster_seed_time[layercluster_seed_time > -90.0])
    arrays["LayerClusterSeedMultiplicity"] = ak.to_numpy(ak.sum(layercluster_seed_mask, axis=1))
    layercluster_cluster_mask = (
        events["LayerCluster_inCluster"] &
        (events["LayerCluster_energy"] > 0.5) &
        (events["LayerCluster_dr"] < 0.3) &
        (events["LayerCluster_Time"] > -90.0)
    )
    layercluster_cluster_weighted_sum = ak.sum(events["LayerCluster_Time"][layercluster_cluster_mask] * events["LayerCluster_energy"][layercluster_cluster_mask], axis=1)
    layercluster_cluster_energy_sum = ak.sum(events["LayerCluster_energy"][layercluster_cluster_mask], axis=1)
    layercluster_cluster_time = ak.fill_none(
        ak.where(layercluster_cluster_energy_sum > 0, layercluster_cluster_weighted_sum / layercluster_cluster_energy_sum, -99.0),
        -99.0,
    )
    arrays["LayerClusterClusterTime"] = ak.to_numpy(layercluster_cluster_time[layercluster_cluster_time > -90.0])
    arrays["LayerClusterClusterMultiplicity"] = ak.to_numpy(ak.sum(layercluster_cluster_mask, axis=1))
    arrays["LayerCluster_Time"] = ak.to_numpy(ak.flatten(events["LayerCluster_Time"][events["LayerCluster_Time"] > -90.0]))
    arrays["LayerCluster_TimeWrtSeed"] = ak.flatten(
        abs(events["LayerCluster_Time"] - ak.broadcast_arrays(seed_time, events["LayerCluster_Time"])[0])
    )
    arrays["Trackster_TimeWrtLayerClusterSeed"] = ak.flatten(
        abs(events["Trackster_Time"] - ak.broadcast_arrays(layercluster_seed_time, events["Trackster_Time"])[0])
    )
    arrays["LayerCluster_TimeWrtLayerClusterSeed"] = ak.flatten(
        abs(events["LayerCluster_Time"] - ak.broadcast_arrays(layercluster_seed_time, events["LayerCluster_Time"])[0])
    )
    arrays["HGCRecHit_TimeWrtLayerClusterSeed"] = ak.flatten(
        abs(events["HGCRecHit_Time"] - ak.broadcast_arrays(layercluster_seed_time, events["HGCRecHit_Time"])[0])
    )
    arrays["Trackster_TimeWrtLayerClusterCluster"] = ak.flatten(
        abs(events["Trackster_Time"] - ak.broadcast_arrays(layercluster_cluster_time, events["Trackster_Time"])[0])
    )
    arrays["LayerCluster_TimeWrtLayerClusterCluster"] = ak.flatten(
        abs(events["LayerCluster_Time"] - ak.broadcast_arrays(layercluster_cluster_time, events["LayerCluster_Time"])[0])
    )
    arrays["HGCRecHit_TimeWrtLayerClusterCluster"] = ak.flatten(
        abs(events["HGCRecHit_Time"] - ak.broadcast_arrays(layercluster_cluster_time, events["HGCRecHit_Time"])[0])
    )
  if "SuperCluster_energy" in events.fields:
    supercluster_mask = (
        (events["SuperCluster_energy"] > 0.0) &
        (events["SuperCluster_dr"] < 0.3)
    )
    arrays["SuperCluster_energy"] = ak.flatten(events["SuperCluster_energy"][supercluster_mask])
    arrays["SuperCluster_rawEnergy"] = ak.flatten(events["SuperCluster_rawEnergy"][supercluster_mask])
    arrays["SuperCluster_eta"] = ak.flatten(events["SuperCluster_eta"][supercluster_mask])
    arrays["SuperCluster_phi"] = ak.flatten(events["SuperCluster_phi"][supercluster_mask])
    arrays["SuperCluster_dr"] = ak.flatten(events["SuperCluster_dr"][supercluster_mask])
    arrays["SuperCluster_nClusters"] = ak.flatten(events["SuperCluster_nClusters"][supercluster_mask])
    arrays["SuperClusterMultiplicity"] = ak.sum(supercluster_mask, axis=1)
    candidate_mask = supercluster_mask & events["SuperCluster_isCandidate"]
    arrays["SuperClusterCandidateEnergy"] = ak.flatten(events["SuperCluster_energy"][candidate_mask])
    arrays["SuperClusterCandidateRawEnergy"] = ak.flatten(events["SuperCluster_rawEnergy"][candidate_mask])
    arrays["SuperClusterCandidateNClusters"] = ak.flatten(events["SuperCluster_nClusters"][candidate_mask])
    arrays["SuperClusterCandidateMultiplicity"] = ak.sum(candidate_mask, axis=1)
  if scenario in OFFLINE_SCENARIOS:
    arrays["Track_TimeWrtPV"] = ak.flatten(abs(events["Track_Time"] - ak.broadcast_arrays(events["PV_Time"], events["Track_Time"])[0]))
    arrays["Trackster_TimeWrtPV"] = ak.flatten(abs(events["Trackster_Time"] - ak.broadcast_arrays(events["PV_Time"], events["Trackster_Time"])[0]))
    arrays["HGCRecHit_TimeWrtPV"] = ak.flatten(abs(events["HGCRecHit_Time"] - ak.broadcast_arrays(events["PV_Time"], events["HGCRecHit_Time"])[0]))
    arrays["Track_TimeWrtSigTrk"] = ak.flatten(abs(events["Track_Time"] - ak.broadcast_arrays(events["sigTrkTime"], events["Track_Time"])[0]))
    arrays["Trackster_TimeWrtSigTrk"] = ak.flatten(abs(events["Trackster_Time"] - ak.broadcast_arrays(events["sigTrkTime"], events["Trackster_Time"])[0]))
    arrays["HGCRecHit_TimeWrtSigTrk"] = ak.flatten(abs(events["HGCRecHit_Time"] - ak.broadcast_arrays(events["sigTrkTime"], events["HGCRecHit_Time"])[0]))
  clean_arrays = {}
  for name, values in arrays.items():
    values = ak.to_numpy(values)
    mask = np.isfinite(values)
    if "Time" in name or "Wrt" in name:
      mask = mask & (values > -90) & (values < 900)
    elif "Multiplicity" in name or name.endswith("nClusters") or name.endswith("NClusters"):
      mask = mask & (values >= 0) & (values < 900)
    else:
      mask = mask & (values > -90)
    clean_arrays[name] = values[mask]
  return clean_arrays


def load_timing_payload_for_scenario(merged_dir, signal_merged_file, background_merged_file, particle, region, signal_process, background_process, scenario):
  scenario_payload = {}
  matched_branch = "matchedToGenEle" if particle == "electron" else "matchedToGenPho"
  print("[timing] scenario={}".format(scenario))
  for process_name, is_signal in [(signal_process, True), (background_process, False)]:
    merged_file = signal_merged_file if is_signal else background_merged_file
    files = resolve_merged_files(merged_dir, merged_file, particle, region, process_name)
    events = read_merged_tree(files, scenario)
    if is_signal:
      events = events[events[matched_branch] == 1]
    else:
      events = events[events[matched_branch] != 1]
    scenario_payload[process_name] = timing_arrays(events, scenario)
  return scenario, scenario_payload


def timing_plot_bins(variable, scenario_payload, valid_scenarios, signal_process, background_process):
  combined = []
  for scenario in valid_scenarios:
    for process_name in [signal_process, background_process]:
      values = scenario_payload[scenario].get(process_name, {}).get(variable)
      if values is None or len(values) == 0:
        continue
      combined.append(values)

  if not combined:
    return None

  combined = np.concatenate(combined)
  combined = combined[np.isfinite(combined)]
  if combined.size == 0:
    return None

  if variable in {"SeedMultiplicity", "LayerClusterSeedMultiplicity", "TracksterClusterMultiplicity", "LayerClusterClusterMultiplicity", "SuperClusterMultiplicity", "SuperClusterCandidateMultiplicity"} or variable.endswith("nClusters") or variable.endswith("NClusters"):
    upper = int(np.ceil(np.quantile(combined, 0.99)))
    upper = max(upper, 1)
    return np.arange(-0.5, upper + 1.5, 1.0)

  low = float(np.quantile(combined, 0.001))
  high = float(np.quantile(combined, 0.995))

  if variable.endswith("Time") or "Wrt" in variable or "SeedTime" in variable:
    low = min(low, 0.0)
  if not np.isfinite(low) or not np.isfinite(high) or high <= low:
    low = float(np.min(combined))
    high = float(np.max(combined))
  if high <= low:
    high = low + 1e-3

  return np.linspace(low, high, 70)


def plot_timing_distributions(merged_dir, signal_merged_file, background_merged_file, particle, region, signal_process, background_process, plotdir, n_workers=1, scenarios=None):
  scenarios = scenarios if scenarios else SCENARIOS
  common_variables = [
      "SeedTime",
      "SeedMultiplicity",
      "TracksterClusterTime",
      "TracksterClusterMultiplicity",
      "Trackster_Time",
      "Trackster_TimeWrtSeed",
      "Trackster_TimeWrtTracksterCluster",
      "HGCRecHit_Time",
      "HGCRecHit_TimeWrtSeed",
      "HGCRecHit_TimeWrtTracksterCluster",
  ]
  online_only = [
      "LayerClusterSeedTime",
      "LayerClusterSeedMultiplicity",
      "LayerClusterClusterTime",
      "LayerClusterClusterMultiplicity",
      "LayerCluster_Time",
      "LayerCluster_TimeWrtSeed",
      "Trackster_TimeWrtLayerClusterSeed",
      "LayerCluster_TimeWrtLayerClusterSeed",
      "HGCRecHit_TimeWrtLayerClusterSeed",
      "Trackster_TimeWrtLayerClusterCluster",
      "LayerCluster_TimeWrtLayerClusterCluster",
      "HGCRecHit_TimeWrtLayerClusterCluster",
      "SuperCluster_energy",
      "SuperCluster_rawEnergy",
      "SuperCluster_eta",
      "SuperCluster_phi",
      "SuperCluster_dr",
      "SuperCluster_nClusters",
      "SuperClusterMultiplicity",
      "SuperClusterCandidateEnergy",
      "SuperClusterCandidateRawEnergy",
      "SuperClusterCandidateNClusters",
      "SuperClusterCandidateMultiplicity",
  ]
  offline_only = [
      "Track_TimeWrtPV",
      "Track_TimeWrtSigTrk",
      "Trackster_TimeWrtPV",
      "Trackster_TimeWrtSigTrk",
      "HGCRecHit_TimeWrtPV",
      "HGCRecHit_TimeWrtSigTrk",
  ]

  scenario_payload = {}
  if n_workers <= 1:
    for scenario in scenarios:
      key, payload = load_timing_payload_for_scenario(
          merged_dir, signal_merged_file, background_merged_file,
          particle, region, signal_process, background_process, scenario,
      )
      scenario_payload[key] = payload
  else:
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
      futures = {
          executor.submit(
              load_timing_payload_for_scenario,
              merged_dir,
              signal_merged_file,
              background_merged_file,
              particle,
              region,
              signal_process,
              background_process,
              scenario,
          ): scenario
          for scenario in scenarios
      }
      future_iter = as_completed(futures)
      if tqdm is not None:
        future_iter = tqdm(future_iter, total=len(futures), desc="Loading timing trees", unit="scenario")
      for future in future_iter:
        key, payload = future.result()
        scenario_payload[key] = payload

  variable_iter = common_variables + online_only + offline_only
  if tqdm is not None:
    variable_iter = tqdm(variable_iter, desc="Timing plots", unit="plot")
  for variable in variable_iter:
    valid_scenarios = []
    for scenario in scenarios:
      if (
          variable in scenario_payload[scenario].get(signal_process, {}) or
          variable in scenario_payload[scenario].get(background_process, {})
      ):
        valid_scenarios.append(scenario)
    if not valid_scenarios:
      continue

    bins = timing_plot_bins(variable, scenario_payload, valid_scenarios, signal_process, background_process)
    if bins is None:
      continue

    for process_name, label in [(signal_process, "signal"), (background_process, "background")]:
      plt.figure(figsize=(8, 6))
      drew = False
      for scenario in valid_scenarios:
        values = scenario_payload[scenario][process_name].get(variable)
        if values is None or len(values) == 0:
          continue
        plt.hist(values, bins=bins, density=True, histtype="step", linewidth=2, label=scenario)
        drew = True
      if not drew:
        plt.close()
        continue
      plt.xlabel(variable)
      plt.ylabel("Density")
      plt.title(f"{variable} ({label})")
      if variable != "SeedMultiplicity":
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
  parser.add_argument("--scenario", action="append", choices=SCENARIOS, help="Optional scenario filter")
  parser.add_argument("--n-workers", default=1, type=int)
  parser.add_argument("--skip-timing", action="store_true", help="Skip merged ROOT timing distributions")
  parser.add_argument("--skip-iso-distributions", action="store_true", help="Skip isolation signal/background distributions")
  parser.add_argument("--histogram-only", action="store_true", help="Only make streaming isolation histograms from parquet columns")
  parser.add_argument("--histogram-roc", action="store_true", help="Use streaming high-granularity histograms for ROC/AUC")
  parser.add_argument("--roc-bins", default=2000, type=int, help="Number of score bins for --histogram-roc")
  parser.add_argument("--min-pt", default=30.0, type=float, help="Minimum eg-pt in GeV for summary plots")
  parser.add_argument(
      "--iso-distribution-dt",
      action="append",
      type=float,
      help="Timing cut value to use for histogram-only isolation plots. Defaults to 9999.0.",
  )
  args = parser.parse_args()

  if args.min_pt < 0.0:
    raise RuntimeError("--min-pt must be non-negative")

  if not args.skip_timing and not args.histogram_only and not args.mergedDir and not args.mergedFile:
    if not args.signal_merged_file or not args.background_merged_file:
      raise RuntimeError("Provide --mergedDir, --mergedFile, or both --signal-merged-file and --background-merged-file")

  check_dir(args.plotDir)

  scenarios = args.scenario if args.scenario else SCENARIOS
  iso_distribution_dt = args.iso_distribution_dt if args.iso_distribution_dt else [9999.0]

  if args.histogram_only:
    plot_histogram_only_iso_distributions(
        args.isoDir,
        args.particle,
        args.region,
        args.signal,
        args.background,
        scenarios,
        args.plotDir,
        iso_distribution_dt,
        args.min_pt,
    )
    return

  if args.histogram_roc:
    if args.roc_bins <= 0:
      raise RuntimeError("--roc-bins must be positive")
    run_streaming_roc_summary(args, scenarios)
    signal_merged_file = args.signal_merged_file if args.signal_merged_file else args.mergedFile
    background_merged_file = args.background_merged_file if args.background_merged_file else args.mergedFile
    if not args.skip_timing:
      plot_timing_distributions(
          args.mergedDir,
          signal_merged_file,
          background_merged_file,
          args.particle,
          args.region,
          args.signal,
          args.background,
          args.plotDir,
          args.n_workers,
          scenarios,
      )
    return

  frames = {}
  if args.n_workers <= 1:
    scenario_iter = scenarios
    if tqdm is not None:
      scenario_iter = tqdm(scenarios, desc="Loading scenarios", unit="scenario")
    for scenario in scenario_iter:
      key, frame = load_scenario_frames(
          args.isoDir,
          args.particle,
          args.region,
          scenario,
          args.signal,
          args.background,
          args.min_pt,
      )
      frames[key] = frame
  else:
    with ThreadPoolExecutor(max_workers=args.n_workers) as executor:
      futures = {
          executor.submit(
              load_scenario_frames,
              args.isoDir,
              args.particle,
              args.region,
              scenario,
              args.signal,
              args.background,
              args.min_pt,
          ): scenario
          for scenario in scenarios
      }
      future_iter = as_completed(futures)
      if tqdm is not None:
        future_iter = tqdm(future_iter, total=len(futures), desc="Loading scenarios", unit="scenario")
      for future in future_iter:
        key, frame = future.result()
        frames[key] = frame

  summary = {}
  iso_distribution_args = (args.isoDir, args.particle, args.region, args.signal, args.background, args.min_pt)
  roc_tasks = [(collection_name, veto_mode) for collection_name in ["trackster", "HGCRecHit"] for veto_mode in ["cone", "seed", "cluster"]]
  roc_iter = roc_tasks
  if tqdm is not None:
    roc_iter = tqdm(roc_tasks, desc="Common ROC scans", unit="scan")
  for collection_name, veto_mode in roc_iter:
      key = f"{collection_name}_{veto_mode}_seed"
      summary[key] = plot_roc_by_scenario(frames, collection_name, veto_mode, "seed", args.plotDir, include_default_online=True)
      plot_auc_scan(frames, collection_name, veto_mode, "seed", args.plotDir)
      if not args.skip_iso_distributions:
        plot_isolation_distributions_for_family(frames, collection_name, veto_mode, "seed", args.plotDir, *iso_distribution_args)
      print("[roc] {} {} ref=seed".format(collection_name, veto_mode))

  online_iter = ["cone", "seed", "cluster"]
  if tqdm is not None:
    online_iter = tqdm(online_iter, desc="Online layerCluster ROC", unit="scan")
  online_frames = {scenario: frames[scenario] for scenario in scenarios if scenario in ONLINE_SCENARIOS and scenario in frames}
  for veto_mode in online_iter:
    key = f"layerCluster_{veto_mode}_seed"
    summary[key] = plot_roc_by_scenario(
        online_frames,
        "layerCluster",
        veto_mode,
        "seed",
        args.plotDir,
    )
    plot_auc_scan(online_frames, "layerCluster", veto_mode, "seed", args.plotDir)
    if not args.skip_iso_distributions:
      plot_isolation_distributions_for_family(online_frames, "layerCluster", veto_mode, "seed", args.plotDir, *iso_distribution_args)
    print("[roc] layerCluster {} ref=seed".format(veto_mode))

  online_layer_seed_tasks = [
      (collection_name, veto_mode)
      for collection_name in ["trackster", "HGCRecHit", "layerCluster"]
      for veto_mode in ["cone", "seed", "cluster"]
  ]
  online_layer_seed_iter = online_layer_seed_tasks
  if tqdm is not None:
    online_layer_seed_iter = tqdm(online_layer_seed_tasks, desc="Online layerClusterSeed ROC", unit="scan")
  for collection_name, veto_mode in online_layer_seed_iter:
    key = f"{collection_name}_{veto_mode}_layerClusterSeed"
    summary[key] = plot_roc_by_scenario(
        online_frames,
        collection_name,
        veto_mode,
        "layerClusterSeed",
        args.plotDir,
        include_default_online=(collection_name != "layerCluster"),
    )
    plot_auc_scan(online_frames, collection_name, veto_mode, "layerClusterSeed", args.plotDir)
    if not args.skip_iso_distributions:
      plot_isolation_distributions_for_family(online_frames, collection_name, veto_mode, "layerClusterSeed", args.plotDir, *iso_distribution_args)
    print("[roc] {} {} ref=layerClusterSeed".format(collection_name, veto_mode))

  online_cluster_ref_tasks = [
      (reference_name, collection_name, veto_mode)
      for reference_name in ["tracksterCluster", "layerClusterCluster"]
      for collection_name in ["trackster", "HGCRecHit", "layerCluster"]
      for veto_mode in ["cone", "seed", "cluster"]
  ]
  online_cluster_ref_iter = online_cluster_ref_tasks
  if tqdm is not None:
    online_cluster_ref_iter = tqdm(online_cluster_ref_tasks, desc="Online inCluster ROC", unit="scan")
  for reference_name, collection_name, veto_mode in online_cluster_ref_iter:
    key = f"{collection_name}_{veto_mode}_{reference_name}"
    summary[key] = plot_roc_by_scenario(
        online_frames,
        collection_name,
        veto_mode,
        reference_name,
        args.plotDir,
        include_default_online=(collection_name != "layerCluster"),
    )
    plot_auc_scan(online_frames, collection_name, veto_mode, reference_name, args.plotDir)
    if not args.skip_iso_distributions:
      plot_isolation_distributions_for_family(online_frames, collection_name, veto_mode, reference_name, args.plotDir, *iso_distribution_args)
    print("[roc] {} {} ref={}".format(collection_name, veto_mode, reference_name))

  offline_frames = {scenario: frames[scenario] for scenario in scenarios if scenario in OFFLINE_SCENARIOS and scenario in frames}
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
        if not args.skip_iso_distributions:
          plot_isolation_distributions_for_family(offline_frames, collection_name, veto_mode, reference_name, args.plotDir, *iso_distribution_args)
        print("[roc] {} {} ref={}".format(collection_name, veto_mode, reference_name))

  with open(os.path.join(args.plotDir, "roc_summary.json"), "w") as handle:
    json.dump(summary, handle, indent=2, sort_keys=True)
  print("[saved] {}".format(os.path.join(args.plotDir, "roc_summary.json")))

  signal_merged_file = args.signal_merged_file if args.signal_merged_file else args.mergedFile
  background_merged_file = args.background_merged_file if args.background_merged_file else args.mergedFile
  if not args.skip_timing:
    plot_timing_distributions(
        args.mergedDir,
        signal_merged_file,
        background_merged_file,
        args.particle,
        args.region,
        args.signal,
        args.background,
        args.plotDir,
        args.n_workers,
        scenarios,
    )


if __name__ == "__main__":
  main()
