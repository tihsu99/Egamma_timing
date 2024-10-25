import awkward as ak
import numpy as np
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import ROOT
import optparse, argparse
import os, sys
from coffea import processor
import hist
import matplotlib.pyplot as plt

def Get_hist(bins, array):
  hist_ = (
    hist.Hist.new
    .Reg(bins[0], bins[1], bins[2], name = 'var')
    .Weight()
    .fill(var = array)
  )
  return hist_

def CheckDir(path):
  if not os.path.exists(path):
    os.system('mkdir -p {}'.format(path))

def Draw_comparison(dataset_dict, histo_name, outputdir, fname, density = True):
  CheckDir(outputdir)
  fig, ax = plt.subplots()
  y_label = 'Density' if density else 'nEntry'
  for dataset_ in dataset_dict:
    linewidth = 3 if '<PU>=0' in dataset_ else 1
    dataset_dict[dataset_]["Histogram"][histo_name].plot1d(ax=ax, label=dataset_, density=density, linewidth=linewidth)
  plt.xlabel('Isolation')
  plt.ylabel(y_label)
  plt.title('PU comparison')
  ax.legend()
  plt.yscale("log")
  plt.savefig(os.path.join(outputdir, fname + '.png'))
  plt.savefig(os.path.join(outputdir, fname + '.pdf'))
  plt.close()

class Accumulator(processor.ProcessorABC):

  def __init__(self, plotVar):
    self.plotVar = plotVar

  def process(self, events):
    histograms = dict()
    dataset = events.metadata['dataset']
    nElectron = len(events)

    # Basic selection
    events = events[events['ele_pt'] > 15]
    if 'DY' in dataset:
      events = events[events['matchedToGenEle']==1]
    else:
      events = events[~(events['matchedToGenEle']==1)]

    # Fill histogram
    for var_ in self.plotVar:
      hist_ = (
        hist.Hist.new
        .Variable(np.arange(0, 20.125, 1.0), name='var')
        .Weight()
      )
      hist_.fill(
        var = events[var_],
        weight = ak.ones_like(events[var_])
      )
      histograms[var_] = hist_

    return {
      dataset:{
        'nEvents': len(events),
        'Histogram': histograms
      }
    }
  def postprocess(self, accumulator):
    pass

def Draw(config):
  
  region = config.region
  dataset = {
    'Prompt(<PU>=0)': ['/eos/user/t/tihsu/database/EGM_timing_w_noPU/electron/{}/DY_noPU/DY_noPU.root'.format(region)],
    'NonPrompt(<PU>=0)': ['/eos/user/t/tihsu/database/EGM_timing_w_noPU/electron/{}/TT_noPU/TT_noPU.root'.format(region)],
    'Prompt(<PU>=200)': ['/eos/user/t/tihsu/database/EGM_timing_w_noPU/electron/{}/DY/DY.root'.format(region)],
    'NonPrompt(<PU>=200)':['/eos/user/t/tihsu/database/EGM_timing_w_noPU/electron/{}/TT/TT.root'.format(region)]
  }

  Var_list = ["eleTrkIsoDtWrtPVMax9999p0dzMax0p15"]


  run = processor.Runner(
        executor = processor.FuturesExecutor(compression = None, workers = 4),
        schema = BaseSchema,
        chunksize = 500000,
        maxchunks = config.maxchunks,
  )

  output = run(
           dataset,
           "ntuplizer/tree",
           processor_instance=Accumulator(Var_list)
  )

  outputdir = 'plot_PU_comparison/{}'.format(region)

  for var_ in Var_list:
    Draw_comparison(output, var_, outputdir, var_)
  

if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--region', default = 'barrel', type=str)
  parser.add_argument('--maxchunks', default = None, type=int)
  config = parser.parse_args()

  Draw(config)
