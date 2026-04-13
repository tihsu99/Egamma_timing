#!/usr/bin/env python3
import argparse
import glob
import os

import awkward as ak
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pyarrow as pa

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.processor import ProcessorABC, Runner, FuturesExecutor
from coffea import processor
import uproot


def compute_isolation(
    objs, energy_name, energy_cut = 0.0,
    ext_dr = 0.15, int_dr = 0.0,
    cut_feature=None, cut_value = 99999,
    veto_seed=False, veto_cluster=False, Ring=False):

    dr_selection = (objs["dr"] < ext_dr) & (objs["dr"] > int_dr)
    energy_selection = (objs[energy_name] > energy_cut)
    objs_select = objs[dr_selection & energy_selection]

    if cut_feature is not None:
      cut_selection = objs_select[cut_feature] < cut_value
      objs_select = objs_select[cut_selection]

    if veto_seed:
        objs_select = objs_select[~objs_select.isSeed]
    if veto_cluster:
        objs_select = objs_select[~objs_select.inCluster]

    if not Ring:
        isolation = ak.sum(objs_select[energy_name], axis=-1)
        return isolation
    else:
        isolations = []
        for ring in range(5):
          for layer in range(6):
            layer_idx = ring * 6 + layer
            layer_hit = objs_select[objs_select.layer == layer_idx]
            isolation = ak.sum(layer_hit[energy_name], axis=-1)
            isolations.append(isolation)
        return isolations


def apply_cut(collection, selection):

  for key in collection:
    collection[key] = collection[key][selection]
  return collection

class EGProcessor(ProcessorABC):
    def __init__(self, process_name, is_signal=True):
        self.process_name = process_name
        self.is_signal = is_signal

        self.iso_collections = {
           f"{var}{ring}_iso-{veto}_veto-dt_{dtcut}": {
              "variable": var,
              "energy": var_energy,
              "energy-cut": (0.02 if veto == "cone" else 0.0) if var == "layerCluster" else 1.0,
              "dr-min": 0.02 if veto == "cone" else 0.0,
              "veto-seed": veto == "seed",
              "veto-cluster": veto == "cluster",
              "dt-cut": dtcut,
              "ring": -1 if ring == "" else int(ring.replace("ring", ""))
           }
           for veto in ["cone", "seed", "cluster"]
           for var, var_energy, rings  in zip(["HGCRecHit", "layerCluster", "trackster"], ["energy", "energy", "pt"], [[f"ring{i}" for i in range(5)] + list([""]), [""], [""]])
           for ring in rings
           for dtcut in [0.1 * (i+1) for i in range(50)] + list([9999.0])
        }

        self.store_values = {
            "eg": ["eta", "phi", "et", "energy", "seedTime"]  + list(self.iso_collections)
        }


        self._accumulator = processor.dict_accumulator({
            f"{obj}-{key}": processor.column_accumulator(np.array([]))
            for obj, keys in self.store_values.items()
            for key in keys
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        collections = dict()
        collections["eg"] = events.eg
        collections["layerCluster"] = events.layerCluster
        collections["trackster"] = events.trackster
        collections["HGCRecHit"] = events.HGCRecHit

        # 0: UNMATCHED, 1:TRUE_PROMPT_ELECTRON, 2:TRUE_ELECTRON_FROM_TAU, 3:TRUE_NON_PROMPT_ELECTRON
        if self.is_signal:
            if self.process_name == "SinglePhoton":
                selection = (collections["eg"].status == 3)
            else:
                selection = (collections["eg"].status == 1)
        else:
            selection = ~(collections["eg"].status == 1)

        collections = apply_cut(collections, selection)
        
        seedTrackster = collections["trackster"][collections["trackster"].isSeed]
        seedTrackster = seedTrackster[seedTrackster.time > -90] # select only nice trackster
        seedTrackster = seedTrackster[seedTrackster.pt > 1]
        seedTrackster = seedTrackster[seedTrackster.dr < 0.3]
        #weighted_time = ak.sum(seedTrackster.time * seedTrackster.pt, axis=1) / ak.sum(seedTrackster.pt, axis=1)
        #avg_seed_time = ak.fill_none(weighted_time, -99)
        avg_seed_time = ak.fill_none(ak.mean(seedTrackster.time, axis = 1), -99)
        collections["eg"]["seedTime"] = avg_seed_time


        trackster_time = collections["trackster"]["time"]
        seed_time = ak.broadcast_arrays(collections["eg"]["seedTime"], trackster_time)[0]
        #mask = (trackster_time > -90) & (seed_time >  -90)
        mask = (seed_time >  -90)
        collections["trackster"]["TimeWrtSig"] = ak.where(mask,
                                                 abs(trackster_time - seed_time),
                                                 0)

        layercluster_time = collections["layerCluster"]["time"]
        seed_time = ak.broadcast_arrays(collections["eg"]["seedTime"], layercluster_time)[0]
 #       mask = (layercluster_time > -90) & (seed_time > -90)
        mask = (seed_time > -90)
        collections["layerCluster"]["TimeWrtSig"] = ak.where(mask,
                                                    abs(layercluster_time - seed_time),
                                                    0)


        rechit_time = collections["HGCRecHit"]["time"]
        seed_time = ak.broadcast_arrays(collections["eg"]["seedTime"], rechit_time)[0]
#        mask = (rechit_time > -90) & (seed_time > -90)
        mask = (seed_time > -90)
        collections["HGCRecHit"]["TimeWrtSig"] = ak.where(mask,
                                                 abs(rechit_time - seed_time),
                                                 0)


        for iso, config in self.iso_collections.items():
          if (config["ring"] > -0.5): 
             continue

          ring_based = "HGCRecHit" in iso

          result = compute_isolation(
            collections[config["variable"]], config["energy"], energy_cut = config["energy-cut"],
            ext_dr = 0.15, int_dr = config["dr-min"],
            veto_seed = config["veto-seed"], veto_cluster=config["veto-cluster"],
            cut_feature="TimeWrtSig", cut_value = config["dt-cut"], Ring = ring_based
          )

          if not ring_based:
            collections["eg"][iso] = result
          else:
            overall = ak.zeros_like(result[0])
            for ring in range(5):
              isoname = iso.replace("HGCRecHit", f"HGCRecHitring{ring}")
              collections["eg"][isoname] = result[ring]
              overall = overall + result[ring]
            collections["eg"][iso] = overall


        for group, keys in self.store_values.items():
            for key in keys:
              output[f"{group}-{key}"] += processor.column_accumulator(ak.to_numpy(collections[group][key]))
        return output

    def postprocess(self, accumulator):
        pass

def save_to_parquet(results, outname):
    """
    Save awkward array to parquet for future analysis.
    """

    df = pd.DataFrame({k: v.value for k, v in results.items()})
    df.to_parquet(outname, index=False)
    print(f"✅ Saved {outname}")
    print(results.keys())


def main():
    parser = argparse.ArgumentParser(description="Coffea script with eg_status selection and isolation output")
    parser.add_argument("--inputDir", required=True, help="Input directory containing ROOT files")
    parser.add_argument("--signal", required=True, help="Signal process name")
    parser.add_argument("--background", required=True, help="Background process name")
    parser.add_argument("--outDir", required=True, help="Output directory for parquet files")
    parser.add_argument("--nCPU", default = 1, type = int)
    parser.add_argument("--max_file", default = -1, type = int)
    args = parser.parse_args()

    os.makedirs(args.outDir, exist_ok=True)


    files = {}
    for proc in [args.signal, args.background]:
      files[proc] = []
      for f in glob.glob(os.path.join(args.inputDir, proc, "output_*.root")):
          if (args.max_file > 0) and (len(files[proc]) > args.max_file):
            continue
          files[proc].append(f)
#          try:
#            with uproot.open(f) as uf:  # ensures file is closed after this block
#              keys = uf.keys()
#            if any(name.split(";")[0] == "eg_endcap" for name in keys):
#              files[proc].append(f)
#          except:
#            print("skip", f)
      if args.max_file > 0:
          num_file = min(len(files[proc]), args.max_file)
          files[proc] = files[proc][:num_file]


    executor = FuturesExecutor(workers=args.nCPU)
    runner = Runner(
        executor=executor,
        schema=NanoAODSchema,
        chunksize=50000,
    )

    signal_output = runner(
        fileset={args.signal: files[args.signal]},
        treename="eg_endcap",
        processor_instance=EGProcessor(args.signal, is_signal=True),
    )

    background_output = runner(
        fileset={args.background: files[args.background]},
        treename="eg_endcap",
        processor_instance=EGProcessor(args.background, is_signal=False),
    )

    # Save outputs
    save_to_parquet(signal_output, os.path.join(args.outDir, f"{args.signal}_iso.parquet"))
    save_to_parquet(background_output, os.path.join(args.outDir, f"{args.background}_iso.parquet"))


if __name__ == "__main__":
    main()

