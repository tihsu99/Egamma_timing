import argparse
import os

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot


def CheckDir(path):
  if not os.path.exists(path):
    os.system("mkdir -p {}".format(path))


def get_efficiency(events, exists_branch, pt_bins):
  efficiency = []
  centers = []
  for low, high in zip(pt_bins[:-1], pt_bins[1:]):
    mask = (events["gen_pt"] >= low) & (events["gen_pt"] < high)
    total = ak.sum(mask)
    passed = ak.sum(events[exists_branch][mask] > 0)
    centers.append(0.5 * (low + high))
    efficiency.append(float(passed) / float(total) if total > 0 else 0.0)
  return np.array(centers), np.array(efficiency)


def get_response(events, pt_branch):
  mask = events[pt_branch] > 0
  return ak.to_numpy((events[pt_branch][mask] - events["gen_pt"][mask]) / events["gen_pt"][mask])


def draw_efficiency(events, plotdir):
  pt_bins = np.array([10, 15, 20, 30, 40, 50, 70, 100, 150, 200])
  labels = [
      ("v4_offline_exists", "offline v4"),
      ("v5_offline_exists", "offline v5"),
      ("v4_online_exists", "online v4"),
      ("v5_online_exists", "online v5"),
  ]

  plt.figure(figsize=(8, 6))
  for branch, label in labels:
    xvals, yvals = get_efficiency(events, branch, pt_bins)
    plt.plot(xvals, yvals, marker="o", label=label)

  plt.xlabel("Gen pT [GeV]")
  plt.ylabel("Efficiency")
  plt.ylim(0.0, 1.05)
  plt.grid(True)
  plt.legend()
  plt.tight_layout()
  plt.savefig(os.path.join(plotdir, "gen_efficiency.png"))
  plt.savefig(os.path.join(plotdir, "gen_efficiency.pdf"))
  plt.close()


def draw_response(events, plotdir):
  response_branches = [
      ("v4_offline_pt", "offline v4"),
      ("v5_offline_pt", "offline v5"),
      ("v4_online_pt", "online v4"),
      ("v5_online_pt", "online v5"),
  ]

  plt.figure(figsize=(8, 6))
  for branch, label in response_branches:
    response = get_response(events, branch)
    plt.hist(response, bins=np.linspace(-1.0, 1.0, 81), density=True, histtype="step", linewidth=2, label=label)

  plt.xlabel("(pT - gen pT) / gen pT")
  plt.ylabel("Density")
  plt.yscale("log")
  plt.grid(True)
  plt.legend()
  plt.tight_layout()
  plt.savefig(os.path.join(plotdir, "gen_response.png"))
  plt.savefig(os.path.join(plotdir, "gen_response.pdf"))
  plt.close()


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Summarize the merged gen-oriented comparison tree")
  parser.add_argument("--fin", required=True, type=str)
  parser.add_argument("--plotdir", default="plot_gen_compare", type=str)
  args = parser.parse_args()

  CheckDir(args.plotdir)
  events = uproot.open(args.fin)["gen_oriented"].arrays(library="ak")
  draw_efficiency(events, args.plotdir)
  draw_response(events, args.plotdir)
