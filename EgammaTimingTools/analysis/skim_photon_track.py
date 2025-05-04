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


  KS_test = dict()
  y_space = 0.0
  for type_ in PU_data:
    ks_stat, ks_p_value = stats.ks_2samp(PU_data[type_], noPU_data[type_])
    KS_test[type_] = ks_stat
    ks_text = f'K-S test({type_}): {ks_stat:.3f}'

    plt.text(0.05, 0.95 - y_space, ks_text, transform=ax.transAxes, fontsize=12,
                 verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    y_space -= 0.08
  
  plt.title(histo_name)
  plt.xlabel(histo_name)
  plt.ylabel(y_label)
  if logy:
    plt.yscale("log")
  ax.legend(title="Sample")
  plt.savefig(os.path.join(outputdir, fname + '.png'))
  plt.savefig(os.path.join(outputdir, fname + '.pdf'))
  plt.close()

  return KS_test

def Draw_scan(dataset_dict, histo_name, outputdir, fname, density = False,  logy=False):

  scan_value = []
  y_value     = []
  for dataset_ in dataset_dict:
    for scan_ in dataset_dict[dataset_]['scanHistograms'][histo_name]:
      KS_dict = Draw_comparison(dataset_dict, histo_name, outputdir, fname + "%.2f"%scan_, density,  logy, keyword = 'scanHistograms', second_keyword = scan_)
      scan_value.append(scan_)
      y_value.append(KS_dict['Prompt'])
    break

  sorted_indices = np.argsort(scan_value)
  sorted_scan_value = np.array(scan_value)[sorted_indices]
  sorted_y_value = np.array(y_value)[sorted_indices]

  plt.figure(figsize=(8, 6))
  plt.plot(sorted_scan_value, sorted_y_value, marker='o', linestyle='-', color='b', label='Prompt')

  # Add labels and title
  plt.xlabel(histo_name)
  plt.ylabel('KS stat')
  plt.legend()


  plt.show()
  # Save the plot to a file
  plt.savefig(os.path.join(outputdir, fname + '.png'))
  plt.savefig(os.path.join(outputdir, fname + '.png'))

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


def Get_hist2D(binsX, binsY, arrayX, arrayY, scale = 1.0):
  hist_ = (
    hist.Hist.new
    .Reg(binsX[0], binsX[1], binsX[2], name = 'varX')
    .Reg(binsY[0], binsY[1], binsY[2], name = 'varY')
    .Weight()
    .fill(varX = arrayX, varY = arrayY, weight = ak.from_numpy(np.ones(len(arrayX))*scale))
  )
  return hist_

def Draw_hist2D(hist, title, xlabel, ylabel, outputdir, fname):

  CheckDir(outputdir)
  fig, ax = plt.subplots()
  hist.plot2d(ax = ax)
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.savefig(os.path.join(outputdir, fname + '.png'))
  plt.savefig(os.path.join(outputdir, fname + '.pdf'))
  plt.close()

class Accumulator(processor.ProcessorABC):
  def __init__(self):
    pass

  def process(self, events):

    histograms  = dict()
    histograms2D = dict()
    scanHistograms = dict()

    dataset     = events.metadata['dataset']
    events = events[events.Pho.pt > 10]
    if 'NonPrompt' in dataset:
      events = events[~(events['matchedToGenPho'] == 1)]
    else:
      events = events[(events['matchedToGenPho'] == 1)]

    events["Weight"] = ak.ones_like(events.Pho.pt)
    histograms['Count']    = Get_hist([1,-0.5,0.5], ak.from_numpy(np.zeros(len(events))), events["Weight"])
    histograms['Pho_pt']   = Get_hist([50, 15, 200], events.Pho.pt, events["Weight"])
    histograms['Pho_eta']  = Get_hist([10, -2.5, 2.5], events.Pho.eta, events["Weight"])
    histograms2D['Pho_pt_v_Pho_eta'] = Get_hist2D([20, 15, 50], [10, -2.5, 2.5], events.Pho.pt, events.Pho.eta)


    Track = events.Track
    PV    = events.PV
    Trackster = events.Trackster
    Photon = events.Pho

    Photon["energy"] = Photon.pt * np.cosh(Photon.eta)

    Track["weight"] = ak.broadcast_arrays(events["Weight"], Track.pt)[0]
    histograms['Track_pt'] = Get_hist([45, 1, 10], ak.flatten(Track.pt), ak.flatten(Track.weight))
    histograms['Track_Time'] = Get_hist([50, 0.0, 0.1], abs(ak.flatten(Track.Time - ak.unflatten(PV.Time, counts = 1))), ak.flatten(Track.weight))
    histograms['Track_Time_zoom'] = Get_hist([50, 0.02, 0.1], abs(ak.flatten(Track.Time - ak.unflatten(PV.Time ,counts = 1))), ak.flatten(Track.weight))
    histograms['Track_TimeErr'] = Get_hist([50, 0, 0.1], ak.flatten(Track.TimeErr), ak.flatten(Track.weight))
    histograms['Track_MtdMva'] = Get_hist([20, 0, 1], ak.flatten(Track.MtdMva), ak.flatten(Track.weight))
    histograms['Track_dz']   = Get_hist([25, 0, 0.5], ak.flatten(abs(Track.dz)), ak.flatten(Track.weight))
    histograms['Track_dxy']  = Get_hist([20, 0, 0.2], ak.flatten(abs(Track.dxy)), ak.flatten(Track.weight))
    histograms['TrackIsolation_bfcut'] = Get_hist([20, 0.0, 20], ak.sum(Track.pt, axis = -1), events.Weight)

    

    Trackster["idx"]  = ak.argsort(ak.argsort(Trackster.pt), ascending=False)
    Trackster["pt_ratio"] = Trackster.pt / ak.unflatten(events.Pho.pt, counts = 1)
    Trackster["weight"]   = ak.broadcast_arrays(events["Weight"], Trackster.pt)[0]

    seedTrackster = Trackster[(Trackster.pt > 1) & (Trackster.dr < 0.30) &  (Trackster.isSeed)]
    Trackster["seedTime"] = ak.broadcast_arrays(ak.fill_none(ak.mean(seedTrackster.Time, axis = 1), -99), Trackster.pt)[0]
    Trackster     = Trackster[(Trackster.pt > 1) & (Trackster.dr < 0.30) & ~(Trackster.isSeed)]

    histograms['Trackster_pt'] = Get_hist([20, 0, 20], ak.flatten(Trackster.pt), ak.flatten(Trackster.weight))
    histograms['Trackster_dr'] = Get_hist([20, 0, 0.3],ak.flatten(Trackster.dr), ak.flatten(Trackster.weight))
    histograms['Trackster_Time'] = Get_hist([50, 0, 0.2], abs(ak.flatten(Trackster.Time - ak.unflatten(PV.Time, counts = 1))), ak.flatten(Trackster.weight))
    histograms['Trackster_Time_zoom'] = Get_hist([49, 0.02, 0.2], abs(ak.flatten(Trackster.Time - ak.unflatten(PV.Time, counts = 1))), ak.flatten(Trackster.weight))
    histograms['Trackster_TimeErr']   = Get_hist([49, 0.02, 0.1], ak.flatten(Trackster.TimeErr), ak.flatten(Trackster.weight))
    histograms['Trackster_HGCsigTime'] = Get_hist([50, 0, 0.1], ak.flatten(abs(Trackster.Time - Trackster.seedTime)), ak.flatten(Trackster.weight))
    histograms['Trackster_pt_ratio']  = Get_hist([50, 0.00, 1.0], ak.flatten(Trackster.pt_ratio), ak.flatten(Trackster.weight))
    histograms['Trackster_Count']     = Get_hist([10, -0.5, 9.5], ak.num(Trackster.pt), events.Weight)
    histograms['Trackster_match_Count'] = Get_hist([10, -0.5, 9.5], ak.num(Trackster.pt[Trackster.Trackpt > 0]), events.Weight)
    histograms['Trackster_match_rate']  = Get_hist([10, 0, 0.5], ak.num(Trackster.pt[Trackster.Trackpt > 0])/ak.num(Trackster.pt), events.Weight)
    histograms['Trackster_Isolation']    = Get_hist([20, 0.0, 20.0], ak.sum(Trackster.pt, axis = -1), events.Weight)
    histograms['Trackster_relIsolation'] = Get_hist([20, 0.0, 1.0],  ak.sum(Trackster.pt, axis = -1)/events.Pho.pt, events.Weight)

    scanHistograms['Trackster_relIsolation_scan_dr'] = OrderedDict()
    for dr_ in np.linspace(0.05, 0.3, 6):
      skimmedTrackster = Trackster[Trackster.dr < dr_]
      scanHistograms['Trackster_relIsolation_scan_dr'][dr_] = Get_hist([20, 0.0, 1.0],  ak.sum(skimmedTrackster.pt, axis = -1)/events.Pho.pt, events.Weight)


    histograms2D['Trackster_Count_v_Pho_pt'] = Get_hist2D([10, -0.5, 9.5], [10, 15, 50], ak.num(Trackster.pt), events.Pho.pt)
    histograms2D['Trackster_dr_v_Trackster_pt'] = Get_hist2D([20, 0, 0.3], [20, 0, 5], ak.flatten(Trackster.dr), ak.flatten(Trackster.pt))
    histograms2D['Trackster_dr_v_Trackster_pt_ratio'] = Get_hist2D([20, 0, 0.3], [20, 0.0, 0.3], ak.flatten(Trackster.dr), ak.flatten(Trackster.pt_ratio))
    histograms2D['Trackster_dr_v_Trackster_Time'] = Get_hist2D([20, 0, 0.3], [20, 0, 0.2], ak.flatten(Trackster.dr), abs(ak.flatten(Trackster.Time - ak.unflatten(PV.Time, counts = 1))))
    histograms2D['Trackster_dr_v_Trackster_HGCsigTime'] = Get_hist2D([20, 0, 0.3], [20, 0, 0.2], ak.flatten(Trackster.dr), abs(ak.flatten(Trackster.Time - Trackster.seedTime)))

    histograms['seedTrackster_pt'] = Get_hist([20, 0, 20], ak.flatten(seedTrackster.pt), ak.flatten(seedTrackster.weight))
    histograms['seedTrackster_dr'] = Get_hist([20, 0, 0.3],ak.flatten(seedTrackster.dr), ak.flatten(seedTrackster.weight))
    histograms['seedTrackster_Time'] = Get_hist([20, 0, 0.3], ak.flatten(seedTrackster.Time - ak.unflatten(PV.Time, counts = 1)), ak.flatten(seedTrackster.weight))
    histograms['seedTrackster_pt_ratio'] = Get_hist([50, 0., 1.0], ak.flatten(seedTrackster.pt_ratio), ak.flatten(seedTrackster.weight))
    histograms['seedTrackster_Count'] = Get_hist([10, -0.5, 9.5], ak.num(seedTrackster.pt), events.Weight)
    histograms['seedTrackster_match_rate'] = Get_hist([10, 0, 0.5], ak.num(seedTrackster.pt[seedTrackster.Trackpt > 0])/ak.num(seedTrackster.pt), events.Weight)

    RecHit     = events.HGCRecHit
    RecHit["weight"] = ak.broadcast_arrays(events["Weight"], RecHit.energy)[0]
    seedRecHit = RecHit[RecHit.isSeed & (RecHit.dr < 0.15)]
    RecHit     = RecHit[~RecHit.isSeed & (RecHit.dr < 0.15)]
    RecHit["seedTime"] = ak.broadcast_arrays(ak.fill_none(ak.mean(seedTrackster.Time,axis = 1), -99), RecHit.energy)[0]
    #RecHit["RecHitSeedTime"] = ak.broadcast_arrays(ak.mean(seedRecHit.Time, weight = seedRecHit.energy, axis = 1), RecHit.energy)[0]

#    histograms['RecHit_energy'] = Get_hist([20, 0, 5],   ak.flatten(RecHit.energy), ak.flatten(RecHit.weight))
#    histograms['RecHit_dr']     = Get_hist([20, 0, 0.3], ak.flatten(RecHit.dr),     ak.flatten(RecHit.weight))
#    histograms['seedRecHit_dr'] = Get_hist([20, 0, 0.3], ak.flatten(seedRecHit.dr), ak.flatten(seedRecHit.weight))
#    histograms['RecHit_layer']  = Get_hist([32, -0.5, 31.5], ak.flatten(RecHit.layer), ak.flatten(RecHit.weight))
#    histograms['seedRecHit_layer'] = Get_hist([32, -0.5, 31.5], ak.flatten(seedRecHit.layer), ak.flatten(seedRecHit.weight))
#    histograms['HGCal_Isolation']  = Get_hist([20, 0.0, 20.0], ak.sum(RecHit.energy, axis = -1), events.Weight)
#    histograms['HGCal_relIsolation'] = Get_hist([20, 0.0, 1.0], ak.sum(RecHit.energy, axis = -1)/Photon.energy, events.Weight)
#    histograms['RecHit_Time']      = Get_hist([49, 0.02, 0.1], abs(ak.flatten(RecHit.Time - ak.unflatten(PV.Time, counts = 1))), ak.flatten(RecHit.weight))
#    histograms['RecHit_HGCsigTime'] = Get_hist([49, 0.02, 0.1], abs(ak.flatten(RecHit.Time - RecHit.seedTime)), ak.flatten(RecHit.weight))
#    histograms['RecHit_HGCRecHitsigTime'] = Get_hist([49, 0.02, 0.1], abs(ak.flatten(RecHit.Time - RecHit.RecHitSeedTime)), ak.flatten(RecHit.weight)) 
#    highest_PT_trackster_idx = ak.argmax(events['Trackster_pt'][TracksterFilter], axis = 1)
#    histograms['highest_PT_Trackster_dr'] = Get_hist([20, 0, 0.3], ak.flatten(events['Trackster_dr'][TracksterFilter][ak.local_index(events['Trackster_pt'][TracksterFilter], axis=1) == highest_PT_trackster_idx[:, None]]), 1.0)

#    nearest_trackster_idx = ak.argmin(events['Trackster_dr'][TracksterFilter], axis = 1)
#    histograms['nearest_Trackster_pt_ratio'] = Get_hist([20, 0, 0.3], ak.flatten(events['Trackster_pt_ratio'][TracksterFilter][ak.local_index(events['Trackster_pt'][TracksterFilter], axis=1) == nearest_trackster_idx[:, None]]), 1.0)

#    if 'Trackster_MVA' in events.fields:
#      histograms['Trackster_MVA']  = Get_hist([20, 0, 1], ak.flatten(events['Trackster_MVA'][TracksterFilter]), TracksterWeight)
#      histograms['HGCalIsolation_Cut0p3'] = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter & (events['Trackster_MVA'] >0.3)], axis = -1), events["Weight"])
#      histograms['HGCalIsolation_Cut0p4'] = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter & (events['Trackster_MVA'] >0.4)], axis = -1), events["Weight"])
#      histograms['HGCalIsolation_Cut0p5'] = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter & (events['Trackster_MVA'] >0.5)], axis = -1), events["Weight"])
#    histograms['HGCalIsolation_timeCut']  = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter & ((events['Trackster_Time'] - ak.unflatten(events['sigTrackster_Time'], counts=1)) < 0.04)], axis = -1), events["Weight"])
#    histograms['HGCalIsolation'] = Get_hist([20, 0, 20], ak.sum(events['Trackster_pt'][TracksterFilter], axis = -1), events["Weight"])




    return{
      dataset: {
        'Histogram': histograms,
        'Histogram2D': histograms2D,
        'scanHistograms': scanHistograms
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

  run = processor.Runner(
    executor = processor.FuturesExecutor(compression=None, workers=4),
    schema=NanoAODSchema,
    chunksize = 2000,
    maxchunks = config.maxchunks,
  )

  output = run(
             dataset,
             "ntuplizer/tree",
             processor_instance = Accumulator()
           )

  Comparison = []
  scanName   = []
  for dataset_ in output:
    for Histogram_name in output[dataset_]["Histogram"]:
      if not Histogram_name in Comparison: Comparison.append(Histogram_name)
    for Histogram_name in output[dataset_]['Histogram2D']:
      xlabel = Histogram_name.split('_v_')[0]
      ylabel = Histogram_name.split('_v_')[1]
      Draw_hist2D(output[dataset_]["Histogram2D"][Histogram_name], Histogram_name, xlabel, ylabel, store_path, Histogram_name + "_" + dataset_)
    for Histogram_name in output[dataset_]['scanHistograms']:
      if not Histogram_name in scanName: scanName.append(Histogram_name)


  for variable_ in Comparison:
    if "Count" in variable_:
      Draw_comparison(output, variable_, store_path, variable_, density=False, logy = ("Isolation" in variable_))
    else:
      Draw_comparison(output, variable_, store_path, variable_, density=True, logy = ("Isolation" in variable_))

  for variable_ in scanName:
    Draw_scan(output, variable_, os.path.join(store_path, variable_), variable_, density = True, logy = ("Isolation" in variable_))

def produce_ntuple(config):

  store_path = '/'.join(config.fout.split('/')[:-1])
  CheckDir(store_path)

  events = NanoEventsFactory.from_root(
        config.fin,
        schemaclass = NanoAODSchema,
        treepath = "ntuplizer/tree"
      ).events()


  events = events[events.Pho.pt > 10]
  threshold = config.threshold
  group = np.where(np.random.rand(len(events)) > threshold, 1, 0)
  train = np.where(np.random.rand(len(events)) > 0.5, 1, 0)
  events["Group"]  = ak.from_numpy(group)
  events["Train"]  = ak.from_numpy(train)

  Trackster = events.Trackster
  PV    = events.PV
  seedTrackster = Trackster[(Trackster.pt > 1) & (Trackster.dr < 0.30) &  (Trackster.isSeed)]
  RecHit     = events.HGCRecHit


  Trackster["seedTime"] = ak.broadcast_arrays(ak.fill_none(ak.mean(seedTrackster.Time, axis = 1), -99), Trackster.pt)[0]
  Trackster["TimeWrtSig"] = abs(Trackster.Time - Trackster.seedTime)
  Trackster["TimeWrtPV"]  = abs(Trackster.Time - ak.unflatten(PV.Time, counts = 1))

  RecHit["seedTime"] = ak.broadcast_arrays(ak.fill_none(ak.mean(seedTrackster.Time,axis = 1), -99), RecHit.energy)[0]
  RecHit["TimeWrtSig"]    = abs(RecHit.Time - RecHit.seedTime)
  RecHit["TimeWrtPV"]     = abs(RecHit.Time - ak.unflatten(PV.Time, counts = 1))

  time_cut_candidate = [0.01 * (i+1) for i in range(20)]
  time_cut_candidate.append(999.0)

  # TracksterIsolation
  for time_cut in time_cut_candidate:
    general_cut = ((~Trackster.isSeed) & (Trackster.dr < 0.15)) & (Trackster.pt > 1)
    Trackster_TimeWrtSig_selected = Trackster[(Trackster["TimeWrtSig"] < time_cut) & general_cut]
    Trackster_TimeWrtPV_selected  = Trackster[(Trackster["TimeWrtPV"]  < time_cut) & general_cut]
    time_str = ("%.2f"%time_cut).replace(".", "p")
    events['TracksterIsolation_TimeWrtSig{}'.format(time_str)]    = ak.sum(Trackster_TimeWrtSig_selected.pt, axis = -1)
    events['TracksterIsolation_TimeWrtPV{}'.format(time_str)]     = ak.sum(Trackster_TimeWrtPV_selected.pt,  axis = -1)

  # TracksterIsolation without all tracksters relevant to cluster.
  for time_cut in time_cut_candidate:
    general_cut = ((~Trackster.inCluster) & (Trackster.dr < 0.15)) & (Trackster.pt > 1)
    Trackster_TimeWrtSig_selected = Trackster[(Trackster["TimeWrtSig"] < time_cut) & general_cut]
    Trackster_TimeWrtPV_selected  = Trackster[(Trackster["TimeWrtPV"]  < time_cut) & general_cut]
    time_str = ("%.2f"%time_cut).replace(".", "p")
    events['TracksterIsolation_ClusterCleaned_TimeWrtSig{}'.format(time_str)]    = ak.sum(Trackster_TimeWrtSig_selected.pt, axis = -1)
    events['TracksterIsolation_ClusterCleaned_TimeWrtPV{}'.format(time_str)]     = ak.sum(Trackster_TimeWrtPV_selected.pt,  axis = -1)


  # HGCRecHitIsolation
  for time_cut in time_cut_candidate:
    general_cut = (~RecHit.isSeed) & (RecHit.dr < 0.15)
    for layer in range(31):
      RecHit_TimeWrtSig_selected = RecHit[general_cut & (RecHit["TimeWrtSig"] < time_cut) & (RecHit.layer == layer)]
      RecHit_TimeWrtPV_selected  = RecHit[general_cut & (RecHit["TimeWrtPV"] < time_cut)  & (RecHit.layer == layer)]
      time_str = ("%.2f"%time_cut).replace(".","p")
      events['RecHitIsolation_TimeWrtSig{}_layer{}'.format(time_str, layer)] = ak.sum(RecHit_TimeWrtSig_selected.energy, axis = -1)
      events['RecHitIsolation_TimeWrtPV{}_layer{}'.format(time_str, layer)] = ak.sum(RecHit_TimeWrtPV_selected.energy, axis = -1)

  # HGCRecHitIsolation without all hits relevant to SC.
  for time_cut in time_cut_candidate:
    general_cut = (~RecHit.inCluster) & (RecHit.dr < 0.15)
    for layer in range(31):
      RecHit_TimeWrtSig_selected = RecHit[general_cut & (RecHit["TimeWrtSig"] < time_cut) & (RecHit.layer == layer)]
      RecHit_TimeWrtPV_selected  = RecHit[general_cut & (RecHit["TimeWrtPV"] < time_cut)  & (RecHit.layer == layer)]
      time_str = ("%.2f"%time_cut).replace(".","p")
      events['RecHitIsolation_ClusterCleaned_TimeWrtSig{}_layer{}'.format(time_str, layer)] = ak.sum(RecHit_TimeWrtSig_selected.energy, axis = -1)
      events['RecHitIsolation_ClusterCleaned_TimeWrtPV{}_layer{}'.format(time_str, layer)] = ak.sum(RecHit_TimeWrtPV_selected.energy, axis = -1)

  events["Trackster"] = Trackster
  events["HGCRecHit"] = RecHit
  file = uproot.recreate(config.fout)
  file["ntuplizer/tree"] = uproot_writeable(events)


if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--fin', type = str)
  parser.add_argument('--fout', type = str)
  parser.add_argument('--threshold', type = float)
  config = parser.parse_args()

  produce_ntuple(config)
