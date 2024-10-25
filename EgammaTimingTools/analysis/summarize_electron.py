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
 

class Accumulator(processor.ProcessorABC):
  def __init__(self, plotVar):
    self.plotVar = plotVar
  def process(self, events):
    dataset = events.metadata['dataset']

    events = events[events['ele_pt'] > 15]
    if dataset == 'DY': # Signal
      events = events[events['matchedToGenEle']==1]
    else: # Background
      events = events[~(events['matchedToGenEle']==1)]
    events_barrel = events[events['ele_isEB']]
    events_endcap = events[events['ele_isEE']] 
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
            weight = ak.ones_like(events_barrel[var_])
        )\
        .fill(
            region = 'endcap',
            var = events_endcap[var_],
            weight = ak.ones_like(events_endcap[var_])
        )
        Histograms[var_] = hist_
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

def obtain_ROC(Dt_list, Other_list, OtherPrefix, region, fin, plotdir='plot', other_ref=0.0):

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
          var_ = ('eleTrkIso{}{}'.format(dtMaxStr, OtherStr))

          hist_DY = fin_root.Get('DY_{}_{}'.format(var_, region)).Clone()
          hist_TT = fin_root.Get('TT_{}_{}'.format(var_, region)).Clone()
          hist_DY.SetDirectory(0)
          hist_TT.SetDirectory(0)

          signal_eff = array('d')
          bkg_eff = array('d')

          nbinx = hist_DY.GetNbinsX()
          for ibinx in range(nbinx + 2):
            signal_eff.append(hist_DY.Integral(-1, 1 + nbinx-ibinx)/float(hist_DY.Integral(-1,nbinx+1)))
            bkg_eff.append(float(hist_TT.Integral(-1, 1 + nbinx-ibinx))/float(hist_TT.Integral(-1, nbinx+1)))
          AUC = Integral_AUC(signal_eff, bkg_eff)
          h_AUC[type_+ref_].SetBinContent(iDt+1, iOther+1, AUC)


          if Other == other_ref and dtMax == 9999.0 and type_=="Dt" and ref_ == "PV":
            graph = ROOT.TGraph(len(signal_eff), signal_eff, bkg_eff)
            graph.SetTitle('default'.format(OtherStr, dtMaxStr))
            graph.SetLineColor(ROOT.kBlack)
            graph.SetLineWidth(5)
            graph.SetMarkerColor(ROOT.kBlack)
            graph.SetMarkerStyle(8)
            graph.SetMarkerSize(0.6)
            mg.Add(graph, "l")
          if AUC > AUC_best:
            AUC_best = AUC
            graph = ROOT.TGraph(len(signal_eff), signal_eff, bkg_eff)
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
      c.SaveAs(os.path.join(plotdir, "AUC_{}Wrt{}_{}_{}.png".format(type_, ref_, OtherPrefix, region))) 
      c.SaveAs(os.path.join(plotdir, "AUC_{}Wrt{}_{}_{}.pdf".format(type_, ref_, OtherPrefix, region))) 

  c = ROOT.TCanvas('','',650, 800)
  c.SetLeftMargin(0.15)
  c.cd()
  mg.SetTitle("ROC;Signal Efficiency; Bkg Efficiency")
  c.DrawFrame(0.8,0.3,1.0,0.75)
  c.Modified()
  mg.Draw('alp')
  mg.GetXaxis().SetLimits(0.7, 1.0);
  c.Modified()
  if region == 'endcap':
    mg.SetMinimum(0.25);
  else:
    mg.SetMinimum(0.15)
  mg.SetMaximum(0.75)
  c.Modified()
  Legend = c.BuildLegend(0.15, 0.5, 0.5, 0.9)
  Legend.SetTextSize(0.02)
  c.SaveAs(os.path.join(plotdir, "ROC_{}_{}.png".format(OtherPrefix, region))) 
  c.SaveAs(os.path.join(plotdir, "ROC_{}_{}.pdf".format(OtherPrefix, region))) 

def obtain_distribution(maxchunks=None, fout='record.root', indir='./'):


  DY_dir = os.path.join(indir, 'endcap/DY')
  TT_dir = os.path.join(indir, 'endcap/TT')
  fileset = {
    'DY': [os.path.join(DY_dir, file_) for file_ in os.listdir(DY_dir)],
    'TT': [os.path.join(TT_dir, file_) for file_ in os.listdir(TT_dir)]
  }

  DY_dir_barrel = os.path.join(indir, 'barrel/DY')
  TT_dir_barrel = os.path.join(indir, 'barrel/TT')
  fileset['DY'] += [os.path.join(DY_dir_barrel, file_) for file_ in os.listdir(DY_dir_barrel)]
  fileset['TT'] += [os.path.join(TT_dir_barrel, file_) for file_ in os.listdir(TT_dir_barrel)]

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
            var_list.append('eleTrkIso{}{}'.format(dtMaxStr, dzMaxStr))
        for iMtdMva, MtdMvaMin in enumerate(l_MtdMvaMin):
            MtdMvaMinStr = "MtdMvaMin%s" %(str(MtdMvaMin).replace(".","p"))
            var_list.append('eleTrkIso{}{}'.format(dtMaxStr, MtdMvaMinStr))
        
        for iptMin, ptMin in enumerate(l_ptMin):
            ptMinStr = "ptMin%s" %(str(ptMin).replace(".","p"))
            var_list.append("eleTrkIso%s%s" %(dtMaxStr, ptMinStr))
  run = processor.Runner(
          executor = processor.FuturesExecutor(compression=None, workers=8),
          schema = BaseSchema,
          chunksize = 2_000,
          maxchunks = maxchunks
  )

  output = run(fileset, 
               "ntuplizer/tree",
               processor_instance=Accumulator(plotVar = var_list)
               )

  f = uproot.recreate(fout)
  for dataset_ in fileset:
    for var_ in var_list:
      f['{}_{}_barrel'.format(dataset_,var_)] = output[dataset_]['Histograms'][var_]['barrel',:]
      f['{}_{}_endcap'.format(dataset_,var_)] = output[dataset_]['Histograms'][var_]['endcap',:]

def draw_distribution(fin_name = 'record.root'):

  l_dtMax = [0.03, 0.12, 0.3, 9999.0]
  for ref_ in ["PV", "SIG"]:
    for type_ in ["Dt", "DtSigni"]:
      dtstr_list = []
      var_list = []
      for iVal, dtMax in enumerate(l_dtMax) :
        dtMaxStr = "%sWrt%sMax%sdzMax0p15" %(type_, ref_, str(dtMax).replace(".", "p"))
        dtstr_list.append(dtMaxStr)
        var_list.append('eleTrkIso{}'.format(dtMaxStr))

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
        hist_dict['DY'] = dict()
        hist_dict['TT'] = dict()
        
        for sample in ['DY', 'TT']:
          for idx, var_ in enumerate(var_list):
            hist = fin.Get('{}_{}_{}'.format(sample, var_, region))
            hist.Rebin(8)
            print(region, sample, var_, hist.Integral(-1, hist.GetNbinsX()+1), color_idx)
            hist.Scale(1./float(hist.Integral(-1, hist.GetNbinsX()+1)))
            hist.SetLineWidth(2)
            hist.SetLineColor(color_template[color_idx])
            hist.SetMarkerStyle(8)
            hist.SetMarkerSize(0.8)
            hist.SetMarkerColor(color_template[color_idx])
            color_idx += 1
            if sample == 'DY':
              hist_dict['DY'][var_] = hist.Clone()
              legend.AddEntry(hist, 'Sig(ZEE), {}'.format(dtstr_list[idx]),'PL')
            else:
              hist_dict['TT'][var_] = hist.Clone()
              legend.AddEntry(hist, 'Bkg(QCD), {}'.format(dtstr_list[idx]),'PL')

            if not already_draw:
              hist.GetYaxis().SetTitle("Events/Integral")
              hist.Draw('PE')
              already_draw = True
            else:
              hist.Draw('SAME PE')
        
        legend.Draw('SAME')
        pad2.cd()
        already_draw_ratio = False
        for idx, var_ in enumerate(var_list):
          sig = hist_dict['DY'][var_]
          bkg = hist_dict['TT'][var_]
          nbinx = bkg.GetNbinsX()
          eff = sig.Clone(sig.GetName())
          eff.SetDirectory(sig.GetDirectory())
          print(sig.GetDirectory(), eff.GetDirectory())
          for ibinx in range(bkg.GetNbinsX()):
            eff.SetBinContent(ibinx+1, float(sig.Integral(-1, ibinx+1))/float((bkg.Integral(-1, ibinx+1))**0.5+ 0.00001))
          for ibinx in range(eff.GetNbinsX()):
            sig.SetBinContent(ibinx+1, eff.GetBinContent(ibinx+1))
          sig.GetYaxis().SetRangeUser(1.05, 2.0)
          if not already_draw_ratio:
            sig.GetYaxis().SetTitle("sig/#sqrt{bkg} for upper cut")
            sig.GetYaxis().SetNdivisions(4)
            sig.GetYaxis().SetTitleSize(0.1)
            sig.GetYaxis().SetLabelSize(0.08)
            sig.GetYaxis().SetTitleOffset(0.3)
            sig.GetXaxis().SetTitle('Isolation')
            sig.GetXaxis().SetLabelSize(0.08)
            sig.GetXaxis().SetTitleSize(0.1)
            sig.Draw('PE')
            already_draw_ratio = True
          else:
            sig.Draw('SAME E')
        c.SaveAs('plot/Isolation_{}_{}_{}.png'.format(region, ref_, type_))  
        c.SaveAs('plot/Isolation_{}_{}_{}.pdf'.format(region, ref_, type_))  
        fin.Close()


  ###########
  ##  ROC  ##
  ###########
  l_dtMax_forROC = [9999.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.3, 0.5, 1.0, 3.0, 5.0]
  l_dzMax = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1.0, 1.5, 2.0]
  l_MtdMvaMin = [0.1, 0.3, 0.5, 0.7, 0.9]
  l_ptMin = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]

  for region in ['endcap', 'barrel']:
        obtain_ROC(Dt_list=l_dtMax_forROC, Other_list=l_dzMax, OtherPrefix="dzMax", region=region, fin=fin_name, other_ref=0.15)
        obtain_ROC(Dt_list=l_dtMax_forROC, Other_list=l_MtdMvaMin, OtherPrefix="MtdMvaMin", region=region, fin=fin_name, other_ref=0.5)
        obtain_ROC(Dt_list=l_dtMax_forROC, Other_list=l_ptMin, OtherPrefix="ptMin", region=region, fin=fin_name, other_ref=1.0)
       

if __name__ == '__main__':
  usage = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--maxchunks', dest='maxchunks', default=-1, type=int)
  parser.add_argument('--fout', dest='fout', default='record.root', type=str)
  parser.add_argument('--fetch_distribution', action ='store_true')
  parser.add_argument('--indir', default = './', type=str)
  parser.add_argument('--draw', action = 'store_true')
  args = parser.parse_args()

  os.system('mkdir -p plot')
  if args.fetch_distribution:
    if args.maxchunks == -1: args.maxchunks = None
    obtain_distribution(maxchunks = args.maxchunks, fout = args.fout, indir=args.indir)
  if args.draw:
    draw_distribution(fin_name = args.fout)
