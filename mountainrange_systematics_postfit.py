import argparse
import pandas as pd
import numpy as np
import re
import ROOT as pyr
import root_numpy
import os

pyr.gROOT.SetBatch(True)
debug = False

parser = argparse.ArgumentParser()
parser.add_argument("fitinput")
parser.add_argument("--cardname", default="datacard")

args = parser.parse_args()

os.system("mkdir fitDiagnostics_plots_{}".format(args.cardname))
fit_file = pyr.TFile(args.fitinput)

process_prefit = {}
total_signal_prefit = []
total_background_prefit = []
total_prefit = []

process_fit_b = {}
total_signal_fit_b = []
total_background_fit_b = []
total_fit_b = []

process_fit_s = {}
total_signal_fit_s = []
total_background_fit_s = []
total_fit_s = []

prefit_dir = fit_file.Get("shapes_prefit")
fit_b_dir = fit_file.Get("shapes_fit_b")
fit_s_dir = fit_file.Get("shapes_fit_s")

def extract_hist(directory, process_dict, total_signal_list, total_background_list, total_list):
    for key in directory.GetListOfKeys():
        print(key.GetName())
        subdir = directory.Get(key.GetName())
        for subkey in subdir.GetListOfKeys():
            print(subkey.GetName())
            print(subkey.GetName() == "total_signal")
            if subkey.GetName() == "data": continue
            hist = subdir.Get(subkey.GetName())
            if "total" not in subkey.GetName():
                if subkey.GetName() not in process_dict.keys(): process_dict[subkey.GetName()] = []
                process_dict[subkey.GetName()] = np.append(process_dict[subkey.GetName()], root_numpy.hist2array(hist))
            elif subkey.GetName() == "total_signal": total_signal_list.append(list(root_numpy.hist2array(hist)))
            elif subkey.GetName() == "total_background": total_background_list.append(list(root_numpy.hist2array(hist)))
            elif subkey.GetName() == "total": total_list.append(list(root_numpy.hist2array(hist)))

def plot_hist(prefit_list, fit_b_list, fit_s_list, hist_name):
    print("Processing {}".format(hist_name))
    print(prefit_list)
    print(fit_b_list)
    print(fit_s_list)
    canvas = pyr.TCanvas("canvas", "Histogram", 1800, 800)

    mountainrange_prefit = pyr.TH1F("mountainrange_prefit", hist_name, len(prefit_list), 0, len(prefit_list))
    mountainrange_fit_b  = pyr.TH1F("mountainrange_fit_b", hist_name, len(fit_b_list), 0, len(fit_b_list))
    mountainrange_fit_s  = pyr.TH1F("mountainrange_fit_s", hist_name, len(fit_s_list), 0, len(fit_s_list))

    root_numpy.array2hist(prefit_list, mountainrange_prefit)
    root_numpy.array2hist(fit_b_list, mountainrange_fit_b)
    root_numpy.array2hist(fit_s_list, mountainrange_fit_s)

    max_boundary = max([max(prefit_list), max(fit_b_list), max(fit_s_list)])*10.
    min_boundary = min([min(prefit_list), min(fit_b_list), min(fit_s_list)])/10.
    if min_boundary == 0: min_boundary = 1e-7

    mountainrange_prefit.Draw("HIST")
    mountainrange_fit_b.Draw("HIST SAME")
    mountainrange_fit_s.Draw("HIST SAME")

    mountainrange_prefit.SetLineColor(pyr.kBlack)
    mountainrange_fit_b.SetLineColor(pyr.kRed)
    mountainrange_fit_s.SetLineColor(pyr.kOrange)

    mountainrange_prefit.SetStats(0)
    mountainrange_fit_b.SetStats(0)
    mountainrange_fit_s.SetStats(0)

    mountainrange_prefit.GetXaxis().SetTitle("Histogram bins")
    pyr.gPad.SetLogy(True)
    mountainrange_prefit.GetYaxis().SetRangeUser(min_boundary, max_boundary)
    mountainrange_fit_b.GetYaxis().SetRangeUser(min_boundary, max_boundary)
    mountainrange_fit_s.GetYaxis().SetRangeUser(min_boundary, max_boundary)
    
    legend = pyr.TLegend(0.8, 0.8, 0.9, 0.9)
    legend.AddEntry(mountainrange_prefit, "prefit")
    legend.AddEntry(mountainrange_fit_b, "fit_b")
    legend.AddEntry(mountainrange_fit_s, "fit_s")
    legend.Draw()
    
    canvas.SaveAs("fitDiagnostics_plots_{}/{}.pdf".format(args.cardname, hist_name))
    canvas.SaveAs("fitDiagnostics_plots_{}/{}.png".format(args.cardname, hist_name))

extract_hist(prefit_dir, process_prefit, total_signal_prefit, total_background_prefit, total_prefit)
extract_hist(fit_b_dir, process_fit_b, total_signal_fit_b, total_background_fit_b, total_fit_b)
extract_hist(fit_s_dir, process_fit_s, total_signal_fit_s, total_background_fit_s, total_fit_s)

total_signal_prefit = [l for sublist in total_signal_prefit for l in sublist]
total_background_prefit = [l for sublist in total_background_prefit for l in sublist]
total_prefit = [l for sublist in total_prefit for l in sublist]

total_signal_fit_b = [l for sublist in total_signal_fit_b for l in sublist]
total_background_fit_b = [l for sublist in total_background_fit_b for l in sublist]
total_fit_b = [l for sublist in total_fit_b for l in sublist]

total_signal_fit_s = [l for sublist in total_signal_fit_s for l in sublist]
total_background_fit_s = [l for sublist in total_background_fit_s for l in sublist]
total_fit_s = [l for sublist in total_fit_s for l in sublist]

for process_name in process_prefit.keys(): plot_hist(process_prefit[process_name], process_fit_b[process_name], process_fit_s[process_name], process_name)
plot_hist(total_signal_prefit, total_signal_fit_b, total_signal_fit_s, "total_signal")
plot_hist(total_background_prefit, total_background_fit_b, total_background_fit_s, "total_background")
plot_hist(total_prefit, total_fit_b, total_fit_s, "total")