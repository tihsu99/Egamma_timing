import awkward as ak
import numpy as np
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
import ROOT
import optparse, argparse
import os, sys
from coffea import processor
import hist
import matplotlib.pyplot as plt
import uproot
from array import array
ROOT.gStyle.SetOptStat(00000000)
ROOT.gStyle.SetPalette(ROOT.kSolar);
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.nanoevents.methods import vector
import scipy.stats as stats
from collections import OrderedDict
from sklearn.metrics import auc

def calculate_roc_auc(signal_histo, bkg_histo, signal_eff_threshold):
    """
    Calculate the ROC curve, AUC, and background efficiency at a specified signal efficiency.

    Parameters:
    - signal_hist (array-like): Histogram counts of the signal.
    - bkg_hist (array-like): Histogram counts of the background.
    - signal_eff_threshold (float): Signal efficiency threshold (in percentage) at which to extract the background efficiency.

    Returns:
    - fpr (array): False positive rate (background efficiency) at different thresholds.
    - tpr (array): True positive rate (signal efficiency) at different thresholds.
    - auc_value (float): The Area Under the ROC Curve (AUC).
    - bkg_eff_at_threshold (float): Interpolated background efficiency when signal efficiency achieves the specified threshold.
    """

    signal_hist, edges = signal_histo.to_numpy()
    bkg_hist,    edges = bkg_histo.to_numpy()
    # Cumulative sums to calculate efficiencies
    sig_cumsum = np.cumsum(signal_hist)
    bkg_cumsum = np.cumsum(bkg_hist)

    # Normalize to get efficiencies
    tpr = (sig_cumsum / sig_cumsum[-1])[::-1]  # Signal efficiency
    fpr = (bkg_cumsum / bkg_cumsum[-1])[::-1]  # Background efficiency

    # Calculate AUC using the trapezoidal rule
    tpr_ = list(tpr)
    fpr_ = list(fpr)
    tpr_.append(0.0)
    fpr_.append(0.0)
    auc_value = auc(fpr_, tpr_)

    # Find the background efficiency at the specified signal efficiency threshold using interpolation
    target_signal_eff = signal_eff_threshold

    if target_signal_eff > tpr[0]:  # If threshold is below the lowest signal efficiency
        bkg_eff_at_threshold = fpr[0]
    elif target_signal_eff <= tpr[-1]:  # If threshold is above the highest signal efficiency
        bkg_eff_at_threshold = fpr[-1]
    else:
        # Perform linear interpolation
        idx = len(tpr) - np.searchsorted(tpr[::-1], target_signal_eff) - 1
        x0, x1 = tpr[idx], tpr[idx + 1]
        y0, y1 = fpr[idx], fpr[idx + 1]
        bkg_eff_at_threshold = y0 + (y1 - y0) * (target_signal_eff - x0) / (x1 - x0)

    return {"fpr": fpr,
            "tpr": tpr,
            "auc": auc_value,
            "Eff(bkg)@0.95sig": bkg_eff_at_threshold}

def CheckDir(path):
  if not os.path.exists(path):
    os.system("mkdir -p {}".format(path))

def Get_hist(bins, array, weight):
    hist_ = (
        hist.Hist.new
        .Reg(bins[0], bins[1], bins[2], name = 'var')
        .Weight()
        .fill(var = array, weight = weight)
      )
    return hist_

def Draw_comparison(dataset_dict, histo_name, outputdir, fname, density = False,  logy=False, keyword = "Histogram", second_keyword = None):
  CheckDir(outputdir)
  fig, ax = plt.subplots()
  y_label = 'Density' if density else 'nEntry'

  PU_data = dict()
  noPU_data = dict()

  for dataset_ in dataset_dict:
    if 'noPU' in dataset_:
      label = dataset_.replace('_noPU', ' <PU>=0')
      linewidth = 3
      data_dict = noPU_data
    else:
      label = dataset_ + " <PU>=200"
      linewidth = 1
      data_dict = PU_data
    label = label.replace("DY", "Prompt").replace("TT", "NonPrompt")
    type_ = "NonPrompt" if "NonPrompt" in label else "Prompt"
    if second_keyword is None:
      histogram = dataset_dict[dataset_][keyword][histo_name]
    else:
      histogram = dataset_dict[dataset_][keyword][histo_name][second_keyword]

    histogram.plot1d(ax=ax, label = label, density = density, linewidth = linewidth)

    values, edges = histogram.to_numpy()
    center = (edges[:-1] + edges[1:]) / 2
    data_dict[type_] = np.repeat(center, (values/np.sum(values) * 1000).astype(int))

  plt.title(histo_name)
  plt.xlabel(histo_name)
  plt.ylabel(y_label)
  if logy:
    plt.yscale("log")
  ax.legend(title="Sample")
  plt.savefig(os.path.join(outputdir, fname + '.png'))
  plt.savefig(os.path.join(outputdir, fname + '.pdf'))
  plt.close()

def Draw_Trend(Data_Dict, InputVariable, PlotVariable, TimeSeries, baselinecut = 999.0, fname = '', second_InputVariable = None):

  x = []
  y = []
  x_second = []
  y_second = []

  plt.figure(figsize=(8, 6))
  for time_cut in TimeSeries:
      time_str = ("%.2f"%time_cut).replace(".","p")
      Dict = Data_Dict['{}{}'.format(InputVariable, time_str)]
      if PlotVariable == "ROC":
        fpr, tpr = Dict['fpr'], Dict['tpr']
        plt.plot(fpr[::-1], tpr[::-1], label = '{}{}'.format(InputVariable, time_str))
      else:
        x.append(time_cut)
        y.append(Dict[PlotVariable])

      if second_InputVariable is not None:
        Dict = Data_Dict['{}{}'.format(second_InputVariable, time_str)]
        if PlotVariable == "ROC":
          fpr, tpr = Dict['fpr'], Dict['tpr']
          plt.plot(fpr[::-1], tpr[::-1], '--', label = '{}{}'.format(second_InputVariable, time_str))
        else:
          x_second.append(time_cut)
          y_second.append(Dict[PlotVariable])

  baseline_time_cut = ("%.2f"%baselinecut).replace(".","p")
  baseline_name = '{}{}'.format(InputVariable, baseline_time_cut)
  baseline = Data_Dict[baseline_name]

  if second_InputVariable is not None:
    baseline_name_second = '{}{}'.format(second_InputVariable, baseline_time_cut)
    baseline_second = Data_Dict[baseline_name_second]

  if PlotVariable == 'ROC':
    base_fpr, base_tpr = baseline['fpr'], baseline['tpr']
    plt.plot(base_fpr[::-1], base_tpr[::-1], label = 'no Time Cut (Seed Cleaning)', linewidth = 3)
    if second_InputVariable is not None: 
      base_fpr_second, base_tpr_second = baseline_second['fpr'], baseline_second['tpr']
      plt.plot(base_fpr_second[::-1], base_tpr_second[::-1], '--', label = 'no Time Cut (Cluster Cleaning)', linewidth = 3)

    plt.xlabel('background efficiency')
    plt.ylabel('signal efficiency')
    plt.title('ROC curve')
    plt.legend(loc='lower right')
  else:
    base_value = baseline[PlotVariable]
    plt.axhline(y = base_value, color = 'r', linestyle = '-', linewidth = 3, label = 'no Time Cut (Seed Cleaning)')
    plt.plot(x, y, label=InputVariable, color='blue', marker='o')
    if second_InputVariable is not None:
      base_value_second = baseline_second[PlotVariable]
      plt.axhline(y = base_value_second, color = 'r', linestyle = '--', linewidth = 3, label = 'no Time Cut (Cluster Cleaning)')
      plt.plot(x_second, y_second, label=second_InputVariable, color='blue', marker='o', linestyle = '--')
    plt.xlabel('dT')
    plt.ylabel(PlotVariable)
    plt.title(InputVariable)
    plt.legend()
  plt.savefig(fname)
  plt.close()

class Accumulator(processor.ProcessorABC):
  def __init__(self, corr = None):
    self.corr = corr

  def process(self, events):

    histograms  = dict()

    dataset     = events.metadata['dataset']
    events = events[events.Pho.pt > 10]
    if 'NonPrompt' in dataset:
      events = events[~(events['matchedToGenPho'] == 1)]
    else:
      events = events[(events['matchedToGenPho'] == 1)]

    corr = self.corr['_noPU'] if '_noPU' in dataset else self.corr['']
    events["weight"] = corr(events.Pho.pt, events.Pho.eta) if not 'NonPrompt' in dataset else ak.ones_like(events.Pho.pt)


    histograms['Count']    = Get_hist([1,-0.5,0.5], ak.from_numpy(np.zeros(len(events))), events["weight"])
    Isolation_bin = [400,-0.1,50]

    time_cut_candidate = [0.04 * (i+1) for i in range(5)]
    time_cut_candidate.append(999.0)
    for time_cut in time_cut_candidate:
      time_str = ("%.2f"%time_cut).replace(".","p")
      histograms['TracksterIsolation_TimeWrtSig{}'.format(time_str)] = Get_hist(Isolation_bin, events.TracksterIsolation["TimeWrtSig{}".format(time_str)], events.weight)
      histograms['TracksterIsolation_ClusterCleaned_TimeWrtSig{}'.format(time_str)] = Get_hist(Isolation_bin, events.TracksterIsolation["ClusterCleaned_TimeWrtSig{}".format(time_str)], events.weight)

      for ring_ in range(5):
        events['HGCIsolationRing{}_TimeWrtSig{}'.format(ring_, time_str)] = ak.zeros_like(events.weight)
        events['HGCIsolationRing{}_ClusterCleaned_TimeWrtSig{}'.format(ring_, time_str)] = ak.zeros_like(events.weight)
        for sublayer_ in range(6):
          layer_ = ring_ * 6 + sublayer_
          events['HGCIsolationRing{}_TimeWrtSig{}'.format(ring_, time_str)]  = events['HGCIsolationRing{}_TimeWrtSig{}'.format(ring_, time_str)] + events.RecHitIsolation['TimeWrtSig{}_layer{}'.format(time_str, layer_)]
          events['HGCIsolationRing{}_ClusterCleaned_TimeWrtSig{}'.format(ring_, time_str)]  = events['HGCIsolationRing{}_ClusterCleaned_TimeWrtSig{}'.format(ring_, time_str)] + events.RecHitIsolation['ClusterCleaned_TimeWrtSig{}_layer{}'.format(time_str, layer_)]

        histograms['HGCIsolationRing{}_TimeWrtSig{}'.format(ring_, time_str)] = Get_hist(Isolation_bin, events['HGCIsolationRing{}_TimeWrtSig{}'.format(ring_, time_str)], events.weight)
        histograms['HGCIsolationRing{}_ClusterCleaned_TimeWrtSig{}'.format(ring_, time_str)] = Get_hist(Isolation_bin, events['HGCIsolationRing{}_ClusterCleaned_TimeWrtSig{}'.format(ring_, time_str)], events.weight)

      events['HGCIsolation_TimeWrtSig{}'.format(time_str)] = ak.zeros_like(events.weight)
      events['HGCIsolation_ClusterCleaned_TimeWrtSig{}'.format(time_str)] = ak.zeros_like(events.weight)
      for ring_ in range(5):
          events['HGCIsolation_TimeWrtSig{}'.format(time_str)] = events['HGCIsolation_TimeWrtSig{}'.format(time_str)] + events['HGCIsolationRing{}_TimeWrtSig{}'.format(ring_, time_str)]
          events['HGCIsolation_ClusterCleaned_TimeWrtSig{}'.format(time_str)] = events['HGCIsolation_ClusterCleaned_TimeWrtSig{}'.format(time_str)] + events['HGCIsolationRing{}_ClusterCleaned_TimeWrtSig{}'.format(ring_, time_str)]
      histograms['HGCIsolation_TimeWrtSig{}'.format(time_str)] = Get_hist(Isolation_bin, events['HGCIsolation_TimeWrtSig{}'.format(time_str)], events.weight)
      histograms['HGCIsolation_ClusterCleaned_TimeWrtSig{}'.format(time_str)] = Get_hist(Isolation_bin, events['HGCIsolation_ClusterCleaned_TimeWrtSig{}'.format(time_str)], events.weight)

    return{
      dataset: {
        'Histogram': histograms,
      }
    }

  def postprocess(self, accumulator):
    pass
def analysis(region, config):
  
  input_path = os.path.join(config.indir, config.particle, region)
  store_path = os.path.join(config.plotdir, config.particle, region)

  CheckDir(input_path)
  CheckDir(store_path)

  dataset = {
    'PromptPhoton': [os.path.join(input_path, 'SinglePhoton2To200', 'SinglePhoton2To200.root')],
    'PromptPhoton_noPU': [os.path.join(input_path, 'SinglePhoton2To200_noPU', 'SinglePhoton2To200_noPU.root')],
    'NonPromptPhoton': [os.path.join(input_path, 'QCDEM', 'QCDEM.root')],
    'NonPromptPhoton_noPU': [os.path.join(input_path, 'QCDEM_noPU', 'QCDEM_noPU.root')]
  }

  process = ['SinglePhoton2To200', 'QCDEM']
  PU      = ['', '_noPU']

  Kinematics  = dict()
  Kinematics_Weight = dict()
  for pileup in PU:
    for process_ in process:
      process_name = process_ + pileup
      events = NanoEventsFactory.from_root(
        os.path.join(input_path, process_name, process_name + ".root"),
        schemaclass = BaseSchema,
        treepath = "ntuplizer/tree"
      ).events()

      events = events[events.Pho_pt > 10]
      dists = (
        hist.Hist.new
        .Reg(100, 0, 200, name = "pt")
        .Reg(16, -3, 3,   name = "eta")
        .Weight()
        .fill(
          pt = events["Pho_pt"],
          eta = events["Pho_eta"]
        )
      )
      Kinematics[process_name] = dists

    num = Kinematics['QCDEM'+pileup][:, :].values()
    den = Kinematics["SinglePhoton2To200"+pileup][:, :].values()
    sf  = np.where(
            (den > 0),
            num / np.maximum(den, 1) * den.sum() / num.sum(),
            1.0,
          )
    corr = dense_lookup(sf, [ax.edges for ax in Kinematics['QCDEM'+pileup].axes[:]])
    Kinematics_Weight[pileup] = corr


  run = processor.Runner(
    executor = processor.FuturesExecutor(compression=None, workers=4),
    schema=NanoAODSchema,
    chunksize = 20000,
    maxchunks = config.maxchunks,
  )

  output = run(
             dataset,
             "ntuplizer/tree",
             processor_instance = Accumulator(Kinematics_Weight)
           )

  Comparison = []
  for dataset_ in output:
    for Histogram_name in output[dataset_]["Histogram"]:
      if not Histogram_name in Comparison: Comparison.append(Histogram_name)

  PU_data = dict()
  
  for variable_ in Comparison:
    if "Count" in variable_:
      Draw_comparison(output, variable_, store_path, variable_, density=False, logy = ("Isolation" in variable_))
    else:
      Draw_comparison(output, variable_, store_path, variable_, density=True, logy = ("Isolation" in variable_))
    PU_data[variable_] = calculate_roc_auc(output['PromptPhoton']['Histogram'][variable_], output['NonPromptPhoton']['Histogram'][variable_], 0.95)

  time_cut_candidate = [0.04 * (i+1) for i in range(5)]
  target_variable = [['TracksterIsolation_TimeWrtSig', 'TracksterIsolation_ClusterCleaned_TimeWrtSig'], ['HGCIsolation_TimeWrtSig', 'HGCIsolation_ClusterCleaned_TimeWrtSig']]
  for ring_ in range(5):
    target_variable.append(['HGCIsolationRing{}_TimeWrtSig'.format(ring_), 'HGCIsolationRing{}_ClusterCleaned_TimeWrtSig'.format(ring_)])

  for target_ in target_variable:
    for plot_variable_ in ['auc', 'Eff(bkg)@0.95sig', 'ROC']:
      Draw_Trend(PU_data, target_[0], plot_variable_, time_cut_candidate, baselinecut = 999.0, fname = os.path.join(store_path, target_[0] + '_' + plot_variable_ + ".png"), second_InputVariable = target_[1])

if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--indir', default = './', type = str)
  parser.add_argument('--maxchunks', default = None, type = int)
  parser.add_argument('--plotdir', default = 'plot/', type = str)
  parser.add_argument('--region', default = ['all'], nargs = '+', type = str)
  parser.add_argument('--particle', default = 'photon', type = str)
  parser.add_argument('--method', default = 'analysis', type=str)
  config = parser.parse_args()

  if 'all' in config.region:
    config.region = ['barrel', 'endcap']

  for region_ in config.region:
      analysis(region_, config)
