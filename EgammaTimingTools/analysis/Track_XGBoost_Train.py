import awkward as ak
import numpy as np
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import optparse, argparse
import os, sys
import uproot
from array import array
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
import pandas as pd
import matplotlib.pyplot as plt

def CheckDir(path):
  if not os.path.exists(path):
    os.system("mkdir -p {}".format(path))

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


def get_xgboost_array(events, type_ = 'Track', Trackster_Filter = False):
  X = []
  #X.append(ak.to_numpy(abs(ak.flatten(events['Track_Time'] - ak.unflatten(events['PV_Time'],counts = 1)))))
  #X.append(ak.to_numpy(ak.flatten(events['Track_TimeErr'])))
  #X.append(ak.to_numpy(ak.flatten(events['Track_MtdMva'])))
  if type_ == 'Track':
    X.append(ak.to_numpy(abs(ak.flatten(events['Track_Time'] - ak.unflatten(events['sigTrkTime'],counts = 1)))))
    sigTrkMtdMva = ak.flatten(ak.broadcast_arrays(events["sigTrkMtdMva"], events["Track_pt"])[0])
    X.append(ak.to_numpy(sigTrkMtdMva))
    X.append(ak.to_numpy(ak.flatten(events['Track_dz'])))
    X.append(ak.to_numpy(ak.flatten(events['Track_dxy'])))
    X    = (np.stack(X, axis = -1))
    Y    = ak.to_numpy(ak.flatten(ak.broadcast_arrays(events['Group'], events['Track_Time'])[0]))
    Feature_map = {"f{}".format(idx): name for idx,name in enumerate(['Track_TimeWrtSig', 'Track_MTDMVA', 'Track_dz', 'Track_dxy'])}
  if type_ == 'Trackster':
    if Trackster_Filter:
         TracksterFilter = (events['Trackster_pt'] > 1) & ((events['Trackster_dr'] > 0.1) & ~(events["Trackster_idx"] == ak.unflatten(events['sigTrackster_Time'], counts=1)))
    else:
         TracksterFilter = (ak.broadcast_arrays(ak.from_numpy(np.ones(len(events))), events['Trackster_pt'])[0])==1
         
    X.append(ak.to_numpy(ak.flatten(abs(events['Trackster_Time'][TracksterFilter] - ak.unflatten(events['sigTrackster_Time'], counts=1)))))
    X.append(ak.to_numpy(ak.flatten(events['Trackster_TimeErr'][TracksterFilter])))
    X.append(ak.to_numpy(abs(ak.flatten(events['Trackster_Time'][TracksterFilter] - ak.unflatten(events['PV_Time'], counts = 1)))))
    X.append(ak.to_numpy(abs(ak.flatten(events['Trackster_Time'][TracksterFilter] - ak.unflatten(events['sigTrkTime'], counts = 1)))))
    X.append(ak.to_numpy(ak.flatten(events['Trackster_dr'][TracksterFilter])))
    X    = (np.stack(X, axis = -1))
    Y    = ak.to_numpy(ak.flatten(ak.broadcast_arrays(events['Group'], events['Trackster_Time'][TracksterFilter])[0]))
    Feature_map = {"f{}".format(idx): name for idx,name in enumerate(['Trackster_TimeWrtSigHGCal', 'Trackster_TimeErr', 'Trackster_TimeWrtPV', 'Trackster_TimeWrtSig', 'Trackster_dr'])}
  return X, Y, Feature_map
def train_xgboost(region, config):

  # Track MVA
  input_path = os.path.join(config.indir, config.particle, region)

  process = ['DY', 'TT', 'DY_noPU', 'TT_noPU']

  input_files = [os.path.join(input_path, process_, process_ + '.root') for process_ in process]

  X = None
  Y = None

  for file_ in input_files:
    print(file_)
    events = NanoEventsFactory.from_root(
          file_,
          schemaclass = BaseSchema,
          treepath = "ntuplizer/tree"
          ).events()

    events = events[(events['Train'] == 1) & (events['Ele_pt'] > 15)]
    if X is None:
      X, Y, FeatureMap = get_xgboost_array(events)
    else:
      X_, Y_,_ = get_xgboost_array(events)
      X = np.concatenate([X,X_], axis = 0)
      Y = np.concatenate([Y,Y_], axis = 0)
  print(np.shape(X), np.shape(Y))
  X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = .1, shuffle=True)
 # bst = XGBClassifier(n_estimators = 5, max_depth = 6, learning_rate = 0.3, objective = 'binary:logistic', verbose = 1)
  bst = XGBClassifier(n_estimators = 500, max_depth = 6, learning_rate = 0.3, objective = 'binary:logistic', verbose = 1)
  bst.fit(X_train, Y_train)
  preds = bst.predict_proba(X_test)
  CheckDir(config.outdir)
  bst.save_model(os.path.join(config.outdir, 'track_xgb.json'))

  feature_important = bst.get_booster().get_score(importance_type='weight')
  keys = list(feature_important.keys())
  values = list(feature_important.values())
  keys_name = [FeatureMap[key] for key in keys]
  data = pd.DataFrame(data=values, index=keys_name, columns=["score"]).sort_values(by = "score", ascending=False)
  print(data)
  data.nlargest(10, columns="score").plot(kind='barh', figsize = (10,10)) ## plot top 40 features
  plt.xlabel("score")
  plt.ylabel("feature")
  plt.title("Variable Importance")
  plt.savefig(os.path.join(config.outdir, "Track_Variable_Importance.png"))


  for file_ in input_files:
     events = NanoEventsFactory.from_root(
          file_,
          schemaclass = BaseSchema,
          treepath = "ntuplizer/tree"
          ).events()

     events = events[events['Train'] == 0]
     X, Y,_ = get_xgboost_array(events)
     Y_pred = ak.from_numpy(bst.predict_proba(X))[:,1]
     length = ak.num(events['Track_pt'])
     print(Y_pred)
     Y_pred = ak.unflatten(Y_pred, counts = length)
     print(Y_pred)
     events['TrackMva'] = Y_pred
     CheckDir('/'.join(file_.replace(config.indir, config.outdir).split('/')[:-1]))
     file = uproot.recreate(file_.replace(config.indir, config.outdir))
     file["ntuplizer/tree"] = uproot_writeable(events)


  #################
  ##  Trackster  ##
  #################

  X = None
  Y = None

  for file_ in input_files:
    print(file_)
    events = NanoEventsFactory.from_root(
          file_,
          schemaclass = BaseSchema,
          treepath = "ntuplizer/tree"
          ).events()

    events = events[(events['Train'] == 1) & (events['Ele_pt'] > 15)]
    if X is None:
      X, Y, FeatureMap = get_xgboost_array(events,  type_ = 'Trackster', Trackster_Filter = True)
    else:
      X_, Y_,_ = get_xgboost_array(events, type_ = 'Trackster', Trackster_Filter = True)
      X = np.concatenate([X,X_], axis = 0)
      Y = np.concatenate([Y,Y_], axis = 0)
  print(np.shape(X), np.shape(Y))
  X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = .1, shuffle=True)
#  bst = XGBClassifier(n_estimators = 5, max_depth = 6, learning_rate = 0.3, objective = 'binary:logistic', verbose = 1)
  bst = XGBClassifier(n_estimators = 2000, max_depth = 6, learning_rate = 0.3, objective = 'binary:logistic', verbose = 1)
  bst.fit(X_train, Y_train)
  preds = bst.predict_proba(X_test)
  CheckDir(config.outdir)
  bst.save_model(os.path.join(config.outdir, 'trackster_xgb.json'))

  feature_important = bst.get_booster().get_score(importance_type='weight')
  keys = list(feature_important.keys())
  values = list(feature_important.values())
  keys_name = [FeatureMap[key] for key in keys]
  data = pd.DataFrame(data=values, index=keys_name, columns=["score"]).sort_values(by = "score", ascending=False)
  print(data)
  data.nlargest(10, columns="score").plot(kind='barh', figsize = (10,10)) ## plot top 40 features
  plt.xlabel("score")
  plt.ylabel("feature")
  plt.title("Variable Importance")
  plt.savefig(os.path.join(config.outdir, "TracksterVariable_Importance.png"))


  for file_ in input_files:
     events = NanoEventsFactory.from_root(
          file_.replace(config.indir, config.outdir), #Update result
          schemaclass = BaseSchema,
          treepath = "ntuplizer/tree"
          ).events()

     events = events[events['Train'] == 0]
     X, Y,_ = get_xgboost_array(events, type_ = 'Trackster')
     Y_pred = ak.from_numpy(bst.predict_proba(X))[:,1]
     length = ak.num(events['Trackster_pt'])
     Y_pred = ak.unflatten(Y_pred, counts = length)
     events['Trackster_MVA'] = Y_pred
     CheckDir('/'.join(file_.replace(config.indir, config.outdir).split('/')[:-1]))
     file = uproot.update(file_.replace(config.indir, config.outdir))
     file["ntuplizer/tree"] = uproot_writeable(events)


if __name__ == '__main__':
  usage = 'usage: %prog [otpions]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--indir', default = './', type=str)
  parser.add_argument('--outdir', default = 'xgboost_result/', type=str)
  parser.add_argument('--region', default = ['all'], nargs = '+', type = str)
  parser.add_argument('--particle', default = 'electron', type = str)
  config = parser.parse_args()

  if 'all' in config.region:
    config.region = ['barrel', 'endcap']

  config.outdir = config.outdir + "/"
  for region_ in config.region:
    train_xgboost(region_, config)
