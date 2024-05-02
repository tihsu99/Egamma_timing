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
from coffea.lookup_tools import extractor


class Accumulator(processor.ProcessorABC):
  def __init__(self, plotVar, weight_file = 'data/kinematic_weight.root'):
    self.plotVar = plotVar
    ext = extractor()
    ext.add_weight_sets(["sigSF2d h_sig_weight {}".format(weight_file)])
    ext.add_weight_sets(["bkgSF2d h_bkg_weight {}".format(weight_file)])
    ext.finalize()
    self.evaluator = ext.make_evaluator()

  def process(self, events):
    dataset = events.metadata['dataset']

    events = events[events['pT'] > 15]
    if dataset in ['SinglePhoton200To500', 'SinglePhoton2To200']: # Signal
      events = events[events['matchedToGenPh']==1]
      weight_name = 'sigSF2d'
    else: # Background
      events = events[~(events['matchedToGenPh']==1)]
      weight_name = 'bkgSF2d'
    events_barrel = events[abs(events['eta']) < 1.44]
    events_endcap = events[(abs(events['eta']) > 1.57) & (abs(events['eta']) < 2.5)] 
    Histograms = dict()
    
    for var_ in self.plotVar:
        hist_ = (
            hist.Hist.new
            .StrCat(['barrel', 'endcap'], name='region')\
            .Variable(np.arange(0, 20.125, 0.125), name='var')
            .Weight()
        )
        hist_.fill(
            region = 'barrel',
            var = events_barrel[var_],
            weight = self.evaluator[weight_name](events_barrel.pT, abs(events_barrel.eta))
        )\
        .fill(
            region = 'endcap',
            var = events_endcap[var_],
            weight = self.evaluator[weight_name](events_endcap.pT, abs(events_endcap.eta))
        )
        Histograms[var_] = hist_

    for var_ in self.plotVar:
        hist_ = (
            hist.Hist.new
            .StrCat(['barrel', 'endcap'], name='region')\
            .Variable(np.arange(0, 5.0, 0.0125), name='var')
            .Weight()
        )
        hist_.fill(
            region = 'barrel',
            var = events_barrel[var_]/events_barrel['pT'],
            weight = self.evaluator[weight_name](events_barrel.pT, abs(events_barrel.eta))
        )\
        .fill(
            region = 'endcap',
            var = events_endcap[var_]/events_endcap['pT'],
            weight = self.evaluator[weight_name](events_endcap.pT, abs(events_endcap.eta))
        )
        Histograms['relative_' + var_] = hist_


    hist_ = (
      hist.Hist.new
      .StrCat(['barrel', 'endcap'], name = 'region')\
      .Variable(np.arange(15, 1000, 5), name='pt')\
      .Variable([0.0, 0.8, 1.44, 1.57, 2.0, 2.5], name='eta')
      .Weight()
    )

    hist_.fill(
      region = 'barrel',
      pt = events_barrel['pT'],
      eta = abs(events_barrel['eta']),
      weight = self.evaluator[weight_name](events_barrel['pT'], abs(events_barrel['eta']))
    )

    hist_.fill(
      region = 'endcap',
      pt = events_endcap['pT'],
      eta = abs(events_endcap['eta']),
      weight = self.evaluator[weight_name](events_endcap['pT'], abs(events_endcap['eta']))
    )
    Histograms['pT_v_eta'] = hist_

    return {
      dataset: {
          "entries": len(events),
          "Histograms": Histograms
      }
    }
  def postprocess(self, accumulator):
    pass


def Integral_AUC(signal_eff, bkg_eff):
  integral_ = 0.0
  for ibin in range(len(signal_eff)-1):
    integral_ += (bkg_eff[ibin] + bkg_eff[ibin+1])/2. * abs(signal_eff[ibin+1] - signal_eff[ibin])
  return 1.0 - integral_

def obtain_ROC(Dt_list, Other_list, OtherPrefix, region, fin, plotdir='plot', other_ref=0.0, signal_list=None, background_list=None, IsoPrefix = ''):

  fin_root = ROOT.TFile.Open(fin, "READ")
  c = ROOT.TCanvas()

  color_template_forROC = [ROOT.kAzure-5, ROOT.kTeal+5, ROOT.kOrange-3, ROOT.kViolet-4, ROOT.kGreen+2, ROOT.kCyan+2, ROOT.kAzure-7, ROOT.kBlue-7, ROOT.kViolet-3, ROOT.kRed-5, ROOT.kRed+1, ROOT.kGreen-10, ROOT.kAzure-9]

  h_AUC = dict()
  mg = ROOT.TMultiGraph()


  ROOT.gStyle.SetPaintTextFormat(".3f")

  for type_ in ["Dt", "DtSigni"]:
    for ref_ in ["PV", "SIG"]:
      h_AUC[type_+ref_] = ROOT.TH2D("",";{};{}".format("%sWrt%sMax"%(type_,ref_), OtherPrefix), len(Dt_list), 0, len(Dt_list), len(Other_list), 0, len(Other_list))
      h_AUC[type_+ref_].SetDirectory(0)


  for iOther, Other in enumerate(Other_list):
    AUC_best = 0.0
    graph = None
    for type_ in ["Dt", "DtSigni"]:
      for ref_ in ["PV", "SIG"]:

        for iDt, dtMax in enumerate(Dt_list):
          dtMaxStr = "%sWrt%sMax%s" %(type_, ref_, str(dtMax).replace(".", "p"))
          OtherStr = "{}{}".format(OtherPrefix, str(Other).replace(".", "p"))
          h_AUC[type_+ref_].GetXaxis().SetBinLabel(iDt+1, str(dtMax))
          h_AUC[type_+ref_].GetYaxis().SetBinLabel(iOther+1, str(Other))
          var_ = (IsoPrefix + 'photonTrkIso{}{}'.format(dtMaxStr, OtherStr))

          hist_sig = None
          hist_bkg = None

          for signal_ in signal_list:
            if hist_sig is None:
              hist_sig = fin_root.Get('{}_{}_{}'.format(signal_, var_, region)).Clone()
              hist_sig.SetDirectory(0)
            else:
              hist_sig.Add(fin_root.Get('{}_{}_{}'.format(signal_, var_, region)).Clone())
          for background_ in background_list:
            if hist_bkg is None:
              hist_bkg = fin_root.Get('{}_{}_{}'.format(background_, var_, region)).Clone()
              hist_bkg.SetDirectory(0)
            else:
              hist_bkg.Add(fin_root.Get('{}_{}_{}'.format(background_, var_, region)).Clone())

          signal_eff = array('d')
          bkg_eff = array('d')

          nbinx = hist_sig.GetNbinsX()
          for ibinx in range(nbinx + 2):
            signal_eff.append(hist_sig.Integral(-1, 1 + nbinx-ibinx)/float(hist_sig.Integral(-1,nbinx+1)+1e-6))
            bkg_eff.append(float(hist_bkg.Integral(-1, 1 + nbinx-ibinx))/float(hist_bkg.Integral(-1, nbinx+1)+1e-6))
          AUC = Integral_AUC(signal_eff, bkg_eff)
          h_AUC[type_+ref_].SetBinContent(iDt+1, iOther+1, AUC)


          bkg_rej = array('d', [1. - bkg_eff_ for bkg_eff_ in bkg_eff])
          if Other == other_ref and dtMax == 9999.0 and type_=="Dt" and ref_ == "PV":
            graph = ROOT.TGraph(len(signal_eff), signal_eff, bkg_rej)
            graph.SetTitle('default'.format(OtherStr, dtMaxStr))
            graph.SetLineColor(ROOT.kBlack)
            graph.SetLineWidth(5)
            graph.SetMarkerColor(ROOT.kBlack)
            graph.SetMarkerStyle(8)
            graph.SetMarkerSize(0.6)
            mg.Add(graph, "l")
          if AUC > AUC_best:
            AUC_best = AUC
            graph = ROOT.TGraph(len(signal_eff), signal_eff, bkg_rej)
            graph.SetTitle(('{}, {}'.format(OtherStr, dtMaxStr)).replace('pt','Pt').replace('p','.').replace('Max','<').replace('Min','>'))
            graph.SetLineColor(color_template_forROC[iOther+1])
            graph.SetLineWidth(2)
            graph.SetMarkerColor(color_template_forROC[iOther+1])
            graph.SetMarkerStyle(8)
            graph.SetMarkerSize(0.6)
    mg.Add(graph, "l")

  for type_ in ["Dt", "DtSigni"]:
    for ref_ in ["PV", "SIG"]:
      h_AUC[type_+ref_].Draw("COLZ TEXT")
      c.SaveAs(os.path.join(plotdir, "AUC_{}{}Wrt{}_{}_{}.png".format(IsoPrefix, type_, ref_, OtherPrefix, region))) 
      c.SaveAs(os.path.join(plotdir, "AUC_{}{}Wrt{}_{}_{}.pdf".format(IsoPrefix, type_, ref_, OtherPrefix, region))) 

  c = ROOT.TCanvas('','',650, 800)
  c.SetLeftMargin(0.15)
  c.cd()
  mg.SetTitle("ROC;Signal Efficiency; Bkg Rejection Rate")
  c.DrawFrame(0.8,0.3,1.0,0.75)
  c.Modified()
  mg.Draw('alp')
  mg.GetXaxis().SetLimits(0.85, 1.0);
  c.Modified()
  if region == 'endcap':
    mg.SetMinimum(0.05);
  else:
    mg.SetMinimum(0.15)
  mg.SetMaximum(1.0)
  c.Modified()
  Legend = c.BuildLegend(0.25, 0.1, 0.65, 0.5)
  Legend.SetTextSize(0.02)
  c.SaveAs(os.path.join(plotdir, "ROC_{}{}_{}.png".format(IsoPrefix, OtherPrefix, region))) 
  c.SaveAs(os.path.join(plotdir, "ROC_{}{}_{}.pdf".format(IsoPrefix, OtherPrefix, region))) 

def obtain_distribution(maxchunks=None, fout='record.root', indir='./', sig_list = None, bkg_list=None):


  process_list = sig_list + bkg_list

  fileset = dict()

  for process in process_list:
    endcap_dir = os.path.join(indir, 'endcap/{}'.format(process))
    barrel_dir = os.path.join(indir, 'barrel/{}'.format(process))
    fileset[process]  = [os.path.join(barrel_dir, "{}.root".format(process))]
    fileset[process] += [os.path.join(endcap_dir, "{}.root".format(process))]


  output = dict()

  l_dtMax = [9999.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.3, 0.5, 1.0, 3.0, 5.0]
  l_dzMax = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1.0, 1.5, 2.0]
  l_MtdMvaMin = [0.1, 0.3, 0.5, 0.7, 0.9]
  l_ptMin = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]
  var_list = []
  for ref_ in ["PV", "SIG"]:
    for type_ in ["Dt", "DtSigni"]:
      for iVal, dtMax in enumerate(l_dtMax) :
        dtMaxStr = "%sWrt%sMax%s" %(type_, ref_, str(dtMax).replace(".", "p"))

        for dzMax in l_dzMax:
            dzMaxStr = "dzMax%s" %(str(dzMax).replace(".","p"))
            var_list.append('photonTrkIso{}{}'.format(dtMaxStr, dzMaxStr))
        for iMtdMva, MtdMvaMin in enumerate(l_MtdMvaMin):
            MtdMvaMinStr = "MtdMvaMin%s" %(str(MtdMvaMin).replace(".","p"))
            var_list.append('photonTrkIso{}{}'.format(dtMaxStr, MtdMvaMinStr))
        
        for iptMin, ptMin in enumerate(l_ptMin):
            ptMinStr = "ptMin%s" %(str(ptMin).replace(".","p"))
            var_list.append("photonTrkIso%s%s" %(dtMaxStr, ptMinStr))
  run = processor.Runner(
          executor = processor.FuturesExecutor(compression=None, workers=8),
          schema = BaseSchema,
          chunksize = 20_000,
          maxchunks = maxchunks
  )

  output = run(fileset, 
               "ntuplizer/tree",
               processor_instance=Accumulator(plotVar = var_list)
               )

  f = uproot.recreate(fout)
  for dataset_ in fileset:
    for var_ in var_list:
      f['{}_rel{}_barrel'.format(dataset_,var_)] = output[dataset_]['Histograms']['relative_' + var_]['barrel',:]
      f['{}_rel{}_endcap'.format(dataset_,var_)] = output[dataset_]['Histograms']['relative_' + var_]['endcap',:]
      f['{}_{}_barrel'.format(dataset_,var_)] = output[dataset_]['Histograms'][var_]['barrel',:]
      f['{}_{}_endcap'.format(dataset_,var_)] = output[dataset_]['Histograms'][var_]['endcap',:]
    f['{}_pT_v_eta_barrel'.format(dataset_)] = output[dataset_]['Histograms']['pT_v_eta']['barrel',:,:]
    f['{}_pT_v_eta_endcap'.format(dataset_)] = output[dataset_]['Histograms']['pT_v_eta']['endcap',:,:]


def draw_distribution(fin_name = 'record.root', sig_list=None, bkg_list=None, relative = False):

  IsoPrefix = 'rel' if relative else ''

  l_dtMax = [0.03, 0.12, 0.3, 9999.0]
  for ref_ in ["PV", "SIG"]:
    for type_ in ["Dt", "DtSigni"]:
      dtstr_list = []
      var_list = []
      for iVal, dtMax in enumerate(l_dtMax) :
        dtMaxStr = "%sWrt%sMax%sdzMax0p15" %(type_, ref_, str(dtMax).replace(".", "p"))
        dtstr_list.append(dtMaxStr)
        var_list.append(IsoPrefix + 'photonTrkIso{}'.format(dtMaxStr))

  ################
  ##  Plotting  ##
  ################

      color_template = [ROOT.kAzure-5, ROOT.kTeal+5, ROOT.kOrange-3, ROOT.kViolet-4, ROOT.kGreen+4, ROOT.kCyan+3, ROOT.kAzure-7, ROOT.kBlue+1, ROOT.kViolet+9, ROOT.kViolet+3]

      for region in ['endcap', 'barrel']:
        fin = ROOT.TFile.Open(fin_name, 'READ')
        c = ROOT.TCanvas("","",620,600)
        pad1 = ROOT.TPad('pad1', '', 0.00, 0.35, 0.99, 0.99)
        pad2 = ROOT.TPad('pad2', '', 0.00, 0.00, 0.99, 0.35)
        pad1.SetBottomMargin(0.01)
        pad1.SetTicks(1,1)
        pad1.SetGrid(1,1)
        pad2.SetTopMargin(0.035)
        pad2.SetBottomMargin(0.45)
        pad2.SetTicks(1,1)
        pad2.SetGrid(1,1)
        pad1.Draw()
        pad2.Draw()
        pad1.cd()

        c.SetGrid(1,1)
        c.SetLeftMargin(0.12)
        c.SetRightMargin(0.08)

        legend = ROOT.TLegend(.65, .5, .9, .9)
        already_draw = False
        color_idx = 0
        pad1.cd()
        pad1.SetLogy()

        hist_dict = dict()


        process_list = sig_list + bkg_list
        hist_dict['sig'] = dict()
        hist_dict['bkg'] = dict()

        for process_ in process_list:
          hist_dict[process_] = dict()

        
        for data_type_ in ['bkg', 'sig']:
          for idx, var_ in enumerate(var_list):
            sample_list = None
            if data_type_ == 'sig':
              sample_list = sig_list
            else:
              sample_list = bkg_list

            hist = None
            for sample in sample_list:
              hist_tmp = fin.Get('{}_{}_{}'.format(sample, var_, region))
              hist_tmp.Rebin(8)
              print(region, sample, var_, hist_tmp.Integral(-1, hist_tmp.GetNbinsX()+1), color_idx)
              if hist is None:
                hist = hist_tmp.Clone()
              else:
                hist.Add(hist_tmp)

            hist.Scale(1./float(hist.Integral(-1, hist.GetNbinsX()+1)+1e-9))
            hist.SetLineWidth(2)
            hist.SetLineColor(color_template[color_idx])
            hist.SetMarkerStyle(8)
            hist.SetMarkerSize(0.8)
            hist.SetMarkerColor(color_template[color_idx])
            color_idx += 1
            if data_type_ == 'sig':
              hist_dict['sig'][var_] = hist.Clone()
              legend.AddEntry(hist, 'Sig(part gun), {}'.format(dtstr_list[idx]),'PL')
            else:
              hist_dict['bkg'][var_] = hist.Clone()
              legend.AddEntry(hist, 'Bkg(QCD), {}'.format(dtstr_list[idx]),'PL')

            if not already_draw:
              hist.GetYaxis().SetTitle("Events/Integral")
              hist.GetYaxis().SetRangeUser(1e-3, 1)
              hist.Draw('PE')
              already_draw = True
            else:
              hist.Draw('SAME PE')
        
        legend.Draw('SAME')
        pad2.cd()
        already_draw_ratio = False
        for idx, var_ in enumerate(var_list):
          sig = hist_dict['sig'][var_]
          bkg = hist_dict['bkg'][var_]
          nbinx = bkg.GetNbinsX()
          eff = sig.Clone(sig.GetName())
          eff.SetDirectory(sig.GetDirectory())
          print(sig.GetDirectory(), eff.GetDirectory())
          for ibinx in range(bkg.GetNbinsX()):
            eff.SetBinContent(ibinx+1, float(sig.Integral(0, ibinx+1))/float((bkg.Integral(0, ibinx+1))**0.5+ 0.00001))
          for ibinx in range(eff.GetNbinsX()):
            sig.SetBinContent(ibinx+1, eff.GetBinContent(ibinx+1))
          sig.GetYaxis().SetRangeUser(1.05, 3.0)
          if not already_draw_ratio:
            sig.GetYaxis().SetTitle("sig/#sqrt{bkg} for upper cut")
            sig.GetYaxis().SetNdivisions(4)
            sig.GetYaxis().SetTitleSize(0.1)
            sig.GetYaxis().SetLabelSize(0.08)
            sig.GetYaxis().SetTitleOffset(0.3)
            sig.GetXaxis().SetTitle(IsoPrefix+'Isolation')
            sig.GetXaxis().SetLabelSize(0.08)
            sig.GetXaxis().SetTitleSize(0.1)
            sig.Draw('PE')
            already_draw_ratio = True
          else:
            sig.Draw('SAME E')
        c.SaveAs('plot/{}Isolation_{}_{}_{}.png'.format(IsoPrefix, region, ref_, type_))  
        c.SaveAs('plot/{}Isolation_{}_{}_{}.pdf'.format(IsoPrefix, region, ref_, type_))  
        fin.Close()


  ###########
  ##  ROC  ##
  ###########
  l_dtMax_forROC = [9999.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.3, 0.5, 1.0, 3.0, 5.0]
  l_dzMax = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1.0, 1.5, 2.0]
  l_MtdMvaMin = [0.1, 0.3, 0.5, 0.7, 0.9]
  l_ptMin = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]

  for region in ['endcap', 'barrel']:
        obtain_ROC(Dt_list=l_dtMax_forROC, Other_list=l_dzMax, OtherPrefix="dzMax", region=region, fin=fin_name, other_ref=0.15, signal_list = sig_list, background_list = bkg_list, IsoPrefix = IsoPrefix)
        obtain_ROC(Dt_list=l_dtMax_forROC, Other_list=l_MtdMvaMin, OtherPrefix="MtdMvaMin", region=region, fin=fin_name, other_ref=0.5, signal_list = sig_list, background_list = bkg_list, IsoPrefix = IsoPrefix)
        obtain_ROC(Dt_list=l_dtMax_forROC, Other_list=l_ptMin, OtherPrefix="ptMin", region=region, fin=fin_name, other_ref=1.0, signal_list = sig_list, background_list = bkg_list, IsoPrefix = IsoPrefix)

def plot_kinematic(fin_name, sig_list, bkg_list, fout_name):
  
  fin = ROOT.TFile.Open(fin_name, 'READ')
  h_sig = None
  h_bkg = None
  for region_ in ['endcap', 'barrel']:
    for sig_ in sig_list:
      if h_sig is None:
        h_sig = fin.Get('{}_pT_v_eta_{}'.format(sig_, region_)).Clone()
      else:
        h_sig.Add(fin.Get('{}_pT_v_eta_{}'.format(sig_, region_)).Clone())
    for bkg_ in bkg_list:
      if h_bkg is None:
        h_bkg = fin.Get('{}_pT_v_eta_{}'.format(bkg_, region_)).Clone()
      else:
        h_bkg.Add(fin.Get('{}_pT_v_eta_{}'.format(bkg_, region_)).Clone())

  c = ROOT.TCanvas()
  h_sig.Draw('COLZ TEXT')
  c.SaveAs('plot/sig_kinematic.png')
  c.SaveAs('plot/sig_kinematic.pdf')
  h_bkg.Draw('COLZ TEXT')
  c.SaveAs('plot/bkg_kinematic.png')
  c.SaveAs('plot/bkg_kinematic.pdf')

  h_sig_weight = h_sig.Clone()
  h_sig_weight.SetDirectory(0)
  h_bkg_weight = h_bkg.Clone()
  h_bkg_weight.SetDirectory(0)

  h_sig_density = h_sig.Clone()
  h_bkg_density = h_bkg.Clone()

  for xbin in range(h_sig_weight.GetNbinsX()):
    for ybin in range(h_sig_weight.GetNbinsY()):
      n_sig = h_sig.GetBinContent(xbin + 1, ybin + 1)
      n_bkg = h_bkg.GetBinContent(xbin + 1, ybin + 1)
      density_sig = n_sig / (h_sig.GetXaxis().GetBinWidth(xbin+1) * h_sig.GetYaxis().GetBinWidth(ybin+1))
      density_bkg = n_bkg / (h_bkg.GetXaxis().GetBinWidth(xbin+1) * h_bkg.GetYaxis().GetBinWidth(ybin+1))
      h_sig_density.SetBinContent(xbin + 1, ybin + 1, density_sig)
      h_bkg_density.SetBinContent(xbin + 1, ybin + 1, density_bkg)
      h_sig_weight.SetBinContent(xbin + 1, ybin + 1, 100./density_sig if density_sig > 0 else 1.0)
      h_bkg_weight.SetBinContent(xbin + 1, ybin + 1, 100./density_bkg if density_bkg > 0 else 1.0)

  h_sig_weight.Draw('COLZ TEXT')
  c.SaveAs('plot/sig_kinematic_weight.png')
  c.SaveAs('plot/sig_kinematic_weight.pdf')
  h_bkg_weight.Draw('COLZ TEXT')
  c.SaveAs('plot/bkg_kinematic_weight.png')
  c.SaveAs('plot/bkg_kinematic_weight.pdf')

  h_sig_density.Draw('COLZ TEXT')
  c.SaveAs('plot/sig_density.png')
  c.SaveAs('plot/sig_density.pdf')
  h_bkg_density.Draw('COLZ TEXT')
  c.SaveAs('plot/bkg_density.png')
  c.SaveAs('plot/bkg_density.pdf')

  fin.Close()
  fout = ROOT.TFile.Open(str(fout_name), 'RECREATE')
  h_sig_weight.Write('h_sig_weight')
  h_bkg_weight.Write('h_bkg_weight')
  fout.Close()

  

if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--maxchunks', dest='maxchunks', default=-1, type=int)
  parser.add_argument('--fout', dest='fout', default='record.root', type=str)
  parser.add_argument('--fetch_distribution', action ='store_true')
  parser.add_argument('--indir', default = './', type=str)
  parser.add_argument('--draw', action = 'store_true')
  parser.add_argument('--plot_kinematic', action = 'store_true')
  args = parser.parse_args()

  os.system('mkdir -p plot')
  sig_list = ['SinglePhoton2To200', 'SinglePhoton200To500']
  bkg_list = ['QCD15to20', 'QCD20to30', 'QCD30to50', 'QCD50to80', 'QCD80to120', 'QCD120to170', 'QCD170to300', 'QCD300toInf']
  if args.fetch_distribution:
    if args.maxchunks == -1: args.maxchunks = None
    obtain_distribution(maxchunks = args.maxchunks, fout = args.fout, indir=args.indir, sig_list=sig_list, bkg_list=bkg_list)
  if args.draw:
    draw_distribution(fin_name = args.fout, sig_list = sig_list, bkg_list = bkg_list)
    draw_distribution(fin_name = args.fout, sig_list = sig_list, bkg_list = bkg_list, relative=True)
  if args.plot_kinematic:
    plot_kinematic(fin_name = args.fout, sig_list = sig_list, bkg_list = bkg_list, fout_name = 'kinematic_weight.root')
