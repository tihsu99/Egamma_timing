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

def Draw_hist(hist, title, outputdir, fname, density = False, logy=False):
  CheckDir(outputdir)
  fig, ax = plt.subplots()
  hist.plot1d(ax=ax)
  y_label = 'Density' if density else 'nEntry'
  plt.title(title)
  if logy:
     plt.yscale("log")
  plt.xlabel(title)
  plt.ylabel(y_label)
  plt.savefig(os.path.join(outputdir, fname + '.png'))
  plt.savefig(os.path.join(outputdir, fname + '.pdf'))
  plt.close()

def Draw_comparison(dataset_dict, histo_name, outputdir, fname, density = False,  logy=False):
  CheckDir(outputdir)
  fig, ax = plt.subplots()
  y_label = 'Density' if density else 'nEntry'
  for dataset_ in dataset_dict:
    if 'noPU' in dataset_:
      label = dataset_.replace('_noPU', ' <PU>=0')
      linewidth = 3
    else:
      label = dataset_ + " <PU>=200"
      linewidth = 1
    label = label.replace("DY", "Prompt").replace("TT", "NonPrompt")
    dataset_dict[dataset_]["Histogram"][histo_name].plot1d(ax=ax, label = label, density = density, linewidth = linewidth)
  plt.title(histo_name)
  plt.xlabel(histo_name)
  plt.ylabel(y_label)
  if logy:
    plt.yscale("log")
  ax.legend(title="Sample")
  plt.savefig(os.path.join(outputdir, fname + '.png'))
  plt.savefig(os.path.join(outputdir, fname + '.pdf'))
  plt.close()

def is_rootcompat(a):
    """Is it a flat or 1-d jagged array?"""
    t = ak.type(a)
    if isinstance(t, ak._ext.ArrayType):
        if isinstance(t.type, ak._ext.PrimitiveType):
            return True
        if isinstance(t.type, ak._ext.ListType) and isinstance(t.type.type, ak._ext.PrimitiveType):
            return True
    return False


def uproot_writeable(events, black_list = []):
    """Restrict to columns that uproot can write compactly"""
    out = {}
    for bname in events.fields:
        if bname in black_list: continue
        if events[bname].fields:
            sub_out = {}
            for n in events[bname].fields:
                if is_rootcompat(events[bname][n]):
                  n_store = n
                  sub_out[n_store] = events[bname][n] 
            out[bname] = ak.packed(ak.without_parameters(ak.zip(sub_out)))
        else: 
            if 'offset' in bname: continue # offset makes TTree entries inconsistent
            out[bname] = ak.packed(ak.without_parameters(events[bname]))
    return out 


class Accumulator(processor.ProcessorABC):
  def __init__(self):
    pass

  def process(self, events):

    histograms  = dict()
    dataset     = events.metadata['dataset']
    events = events[events['Ele_pt'] > 15]
    if 'DY' in events:
      events = events[events['matchedToGenEle'] == 1]
    else:
      events = events[~(events['matchedToGenEle'] == 1)]
    histograms['Count'] = Get_hist([1,-0.5,0.5], ak.from_numpy(np.zeros(len(events))), events["Weight"])
    #events = events[events['sigTrkMtdMva'] > -1]
    histograms['Count_afterEventCut'] = Get_hist([1,-0.5,0.5], ak.from_numpy(np.zeros(len(events))), events["Weight"])
    histograms['Group']       = Get_hist([2, -0.5, 1.5], events['Group'], events["Weight"])
    histograms['Train']       = Get_hist([2, -0.5, 1.5], events['Train'], events['Weight'])
    histograms['Ele_pt']      = Get_hist([50, 15, 100], events['Ele_pt'], events["Weight"])
    histograms['Ele_eta']      = Get_hist([10, -2.5, 2.5], events['Ele_eta'], events["Weight"])
    histograms['sigTrkTime']  = Get_hist([50, 0, 1], events['sigTrkTime'] - events['PV_Time'], events["Weight"])
    histograms['sigTrkTimeErr'] = Get_hist([50, 0, 0.1], events['sigTrkTimeErr'], events["Weight"])
    histograms['sigTrkMtdMva']  = Get_hist([20, 0, 1], events['sigTrkMtdMva'], events["Weight"])
    TrackWeight = ak.broadcast_arrays(events["Weight"], events['Track_pt'])
    TrackWeight = ak.flatten(TrackWeight[0])
    histograms['Track_pt'] = Get_hist([45, 1, 10], ak.flatten(events['Track_pt']), TrackWeight)
    histograms['Track_Time'] = Get_hist([50, 0.0, 0.1], abs(ak.flatten(events['Track_Time'] - ak.unflatten(events['PV_Time'],counts = 1))), TrackWeight)
    histograms['Track_Time_zoom'] = Get_hist([50, 0.02, 0.1], abs(ak.flatten(events['Track_Time'] - ak.unflatten(events['PV_Time'],counts = 1))), TrackWeight)
    histograms['Track_TimeErr'] = Get_hist([50, 0, 0.1], ak.flatten(events['Track_TimeErr']), TrackWeight)
    histograms['Track_MtdMva'] = Get_hist([20, 0, 1], ak.flatten(events['Track_MtdMva']), TrackWeight)
    histograms['Track_dz']   = Get_hist([25, 0, 0.5], ak.flatten(events['Track_dz']), TrackWeight)
    histograms['Track_dxy']  = Get_hist([20, 0, 0.2], ak.flatten(events['Track_dxy']), TrackWeight)
    histograms['Track_TimeWrtSig'] = Get_hist([50, 0.00, 0.1], abs(ak.flatten(events['Track_Time'] - ak.unflatten(events['sigTrkTime'],counts = 1))),TrackWeight)
    histograms['Track_TimeWrtSig_zoom'] = Get_hist([50, 0.02, 0.1], abs(ak.flatten(events['Track_Time'] - ak.unflatten(events['sigTrkTime'],counts = 1))), TrackWeight)
    histograms['Isolation_bfcut'] = Get_hist([20, 0.0, 20], ak.sum(events['Track_pt'], axis = -1), events["Weight"])

    if 'TrackMva' in events.fields:
      histograms['Track_Mva'] = Get_hist([20, 0, 1.0], ak.flatten(events['TrackMva']), TrackWeight)
      histograms['Isolation'] = Get_hist([20, 0, 20], ak.sum(events['Track_pt'][events['TrackMva'] > 0.4], axis = -1), events["Weight"])
      histograms['Isolation_Cut0p5'] = Get_hist([20, 0, 20], ak.sum(events['Track_pt'][events['TrackMva'] > 0.5], axis = -1), events["Weight"])
      histograms['Isolation_dzCut']  =  Get_hist([20, 0, 20], ak.sum(events['Track_pt'][events['Track_dz'] < 0.15], axis = -1), events["Weight"])


    events["Trackster_idx"] = ak.argsort(ak.argsort(events["Trackster_pt"]), ascending=False)

    # Finding signal Trackster
    InnerTracksterFilter = (events['Trackster_pt'] > 1) & (events['Trackster_dr'] < 0.15)
    histograms['InnerTracksterFilter'] = Get_hist([6, -0.5, 5.5], ak.num(events['Trackster_pt'][InnerTracksterFilter]), events["Weight"])
    signal_tracker_idx   = events['Trackster_idx'][InnerTracksterFilter][ak.argsort(events['Trackster_pt'][InnerTracksterFilter], ascending = False)]
    signal_tracker_idx   = ak.from_regular(ak.pad_none(signal_tracker_idx, 1, axis = -1, clip=True))

    print(signal_tracker_idx)
    print(events['Trackster_pt'][signal_tracker_idx])

    events["sigTrackster_Time"] = ak.flatten(ak.fill_none(events["Trackster_Time"][signal_tracker_idx], -1))

    signal_tracker_idx = ak.fill_none(signal_tracker_idx, -1)
    print(signal_tracker_idx)
    TracksterFilter = (events['Trackster_pt'] > 1) & ((events['Trackster_dr'] > 0.1) & ~(events["Trackster_idx"] == ak.flatten(signal_tracker_idx)))
    TracksterWeight = ak.broadcast_arrays(events["Weight"], events["Trackster_pt"][TracksterFilter])
    TracksterWeight = ak.flatten(TracksterWeight[0])

    histograms['Trackster_pt'] = Get_hist([20, 0, 20], ak.flatten(events['Trackster_pt'][TracksterFilter]), TracksterWeight)
    histograms['Trackster_dr_raw'] = Get_hist([20, 0, 0.3], ak.flatten(events['Trackster_dr'][events['Trackster_pt'] > 1]), ak.flatten(ak.broadcast_arrays(events["Weight"], events["Trackster_pt"][events['Trackster_pt'] > 1])[0]))
    histograms['Trackster_dr'] = Get_hist([20, 0, 0.3], ak.flatten(events['Trackster_dr'][TracksterFilter]), TracksterWeight)
    histograms['Trackster_Time'] = Get_hist([50, 0, 0.2], abs(ak.flatten(events['Trackster_Time'][TracksterFilter] - ak.unflatten(events['PV_Time'], counts = 1))), TracksterWeight)
    histograms['Trackster_Time_zoom'] = Get_hist([50, 0.02, 2.0], abs(ak.flatten(events['Trackster_Time'][TracksterFilter] - ak.unflatten(events['PV_Time'], counts = 1))), TracksterWeight)
    histograms['Trackster_TimeErr'] = Get_hist([50, 0, 0.1], ak.flatten(events['Trackster_TimeErr'][TracksterFilter]), TracksterWeight)
    histograms['Trackster_TimeWrtSig'] = Get_hist([50, 0.00, 0.2], abs(ak.flatten(events['Trackster_Time'][TracksterFilter] - ak.unflatten(events['sigTrkTime'], counts=1))), TracksterWeight)
    histograms['Trackster_TimeWrtSig_zoom'] = Get_hist([50, 0.02, 2.0], abs(ak.flatten(events['Trackster_Time'][TracksterFilter] - ak.unflatten(events['sigTrkTime'], counts=1))), TracksterWeight)
    histograms['Trackster_TimeWrtSigHGCal'] = Get_hist([50, 0.00, 0.1], abs(ak.flatten(events['Trackster_Time'][TracksterFilter] - ak.unflatten(events['sigTrackster_Time'], counts=1))), TracksterWeight)
    histograms['Trackster_TimeWrtSigHGCal_zoom'] = Get_hist([50, 0.02, 2.0], abs(ak.flatten(events['Trackster_Time'][TracksterFilter] - ak.unflatten(events['sigTrackster_Time'], counts=1))), TracksterWeight)
    if 'Trackster_MVA' in events.fields:
      histograms['Trackster_MVA']  = Get_hist([20, 0, 1], ak.flatten(events['Trackster_MVA'][TracksterFilter]), TracksterWeight)
      histograms['HGCalIsolation_Cut0p3'] = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter & (events['Trackster_MVA'] >0.3)], axis = -1), events["Weight"])
      histograms['HGCalIsolation_Cut0p4'] = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter & (events['Trackster_MVA'] >0.4)], axis = -1), events["Weight"])
      histograms['HGCalIsolation_Cut0p5'] = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter & (events['Trackster_MVA'] >0.5)], axis = -1), events["Weight"])
    histograms['HGCalIsolation_timeCut']  = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter & ((events['Trackster_Time'] - ak.unflatten(events['sigTrackster_Time'], counts=1)) < 0.04)], axis = -1), events["Weight"])
    histograms['HGCalIsolation'] = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter], axis = -1), events["Weight"])

    return{
      dataset: {
        'Histogram': histograms
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
    'DY_noPU': [os.path.join(input_path, 'DY_noPU', 'DY_noPU.root')],
    'TT_noPU': [os.path.join(input_path, 'TT_noPU', 'TT_noPU.root')],
    'DY':      [os.path.join(input_path, 'DY', 'DY.root')],
    'TT':      [os.path.join(input_path, 'TT', 'TT.root')]
  }

  run = processor.Runner(
    executor = processor.FuturesExecutor(compression=None, workers=4),
    schema=BaseSchema,
    chunksize = 5000000,
    maxchunks = config.maxchunks,
  )

  output = run(
             dataset,
             "ntuplizer/tree",
             processor_instance = Accumulator()
           )

  Comparison = []
  for dataset_ in output:
    for Histogram_name in output[dataset_]["Histogram"]:
      if not Histogram_name in Comparison: Comparison.append(Histogram_name)

  for variable_ in Comparison:
    if "Count" in variable_:
      Draw_comparison(output, variable_, store_path, variable_, density=False, logy = ("Isolation" in variable_))
    else:
      Draw_comparison(output, variable_, store_path, variable_, density=True, logy = ("Isolation" in variable_))

def produce_ntuple(region, config):
  input_path  = os.path.join(config.indir, config.particle, region)
  store_path = os.path.join(config.plotdir, config.particle, region)

  CheckDir(store_path)

  process = ['DY', 'TT']
  PU      = ['', '_noPU']
  for pileup in PU:
    Events  = dict()
    for process_ in process:
      process_name = process_ + pileup
      events = NanoEventsFactory.from_root(
        os.path.join(input_path, process_name, process_name + ".root"),
        schemaclass = BaseSchema,
        treepath = "ntuplizer/tree"
      ).events()

      events = events[events["Ele_pt"] > 15]
      Events[process_] = events
    dists = (
      hist.Hist.new
      .StrCat(["DY", "TT"], name="dataset")
      .Reg(40, 0, 200, name = "pt")
      .Reg(8, -3, 3,   name = "eta")
      .Weight()
      .fill(
        dataset = "DY",
        pt = Events["DY"]["Ele_pt"],
        eta = Events["DY"]["Ele_eta"]
      )
      .fill(
        dataset = "TT",
        pt = Events["TT"]["Ele_pt"],
        eta = Events["TT"]["Ele_eta"]
      )
    )
    num = dists["TT", :, :].values()
    den = dists["DY", :, :].values()
    sf  = np.where(
            (den > 0),
            num / np.maximum(den, 1) * den.sum() / num.sum(),
            1.0,
          )
    corr = dense_lookup(sf, [ax.edges for ax in dists.axes[1:]])
    for process_ in process:
      Events[process_]["Weight"] = corr(Events[process_]["Ele_pt"], Events[process_]["Ele_eta"]) if process_ == "DY" else ak.from_numpy(np.ones(len(Events[process_])))#Kinematic weight 
      threshold = 0.3 if pileup == '_noPU' else 0.7
      group = np.where(np.random.rand(len(Events[process_])) > threshold, 1, 0)
      train = np.where(np.random.rand(len(Events[process_])) > 0.5, 1, 0)
      Events[process_]["Group"]  = ak.from_numpy(group)
      Events[process_]["Train"]  = ak.from_numpy(train)

      Events[process_]["Trackster_idx"] = ak.argsort(ak.argsort(Events[process_]["Trackster_pt"]), ascending=False)
      InnerTracksterFilter = (Events[process_]['Trackster_pt'] > 1) & (Events[process_]['Trackster_dr'] < 0.15)
      signal_tracker_idx   = Events[process_]['Trackster_idx'][InnerTracksterFilter][ak.argsort(Events[process_]['Trackster_pt'][InnerTracksterFilter], ascending = False)]
      signal_tracker_idx   = ak.from_regular(ak.pad_none(signal_tracker_idx, 1, axis = -1, clip=True))

      Events[process_]["sigTrackster_Time"] = ak.flatten(ak.fill_none(Events[process_]["Trackster_Time"][signal_tracker_idx], -1))
      signal_tracker_idx = ak.fill_none(signal_tracker_idx, -1)

      process_name = process_ + pileup
      print(process_name)
      output_path = os.path.join(store_path, process_name)
      CheckDir(output_path)
      file = uproot.recreate(os.path.join(output_path, process_name + ".root"))
      file["ntuplizer/tree"] = uproot_writeable(Events[process_])


if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--indir', default = './', type = str)
  parser.add_argument('--maxchunks', default = None, type = int)
  parser.add_argument('--plotdir', default = 'plot/', type = str)
  parser.add_argument('--region', default = ['all'], nargs = '+', type = str)
  parser.add_argument('--particle', default = 'electron', type = str)
  parser.add_argument('--method', default = 'analysis', type=str)
  config = parser.parse_args()

  if 'all' in config.region:
    config.region = ['barrel', 'endcap']

  if config.method == 'analysis':
    for region_ in config.region:
      analysis(region_, config)
  elif config.method == 'ntuple':
    for region_ in config.region:
      produce_ntuple(region_, config)
