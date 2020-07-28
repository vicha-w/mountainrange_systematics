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
parser.add_argument("filename")
parser.add_argument("--mass", default=690)
parser.add_argument("--cardname", default="datacard")
parser.add_argument("--ratio", action="store_true")

args = parser.parse_args()

infile = open(args.filename, 'r')

if args.ratio: os.system("mkdir mountainrange_ratio_{}".format(args.cardname))
else: os.system("mkdir mountainrange_{}".format(args.cardname))

# Initial section. Skip for now.
for line_num, line in enumerate(infile):
    if line_num == 0: continue
    if debug: print(line.strip())
    if not (line.startswith("imax") or line.startswith("jmax") or line.startswith("kmax")): break

# Shapes section
shape_lines = []
for line_num, line in enumerate(infile):
    if not line.startswith("shape"): break
    shape_lines.append(line.strip().split()[1:])
if debug: print(shape_lines)
shape_df = pd.DataFrame(shape_lines, columns=["process", "channel", "file", "central_hist", "syst_hist"])

# Observation section. Skip for now.
for line_num, line in enumerate(infile):
    if debug: print(line.strip())
    if line.startswith("-"): 
        break

process_lines = []
for line_num, line in enumerate(infile):
    if line.startswith("-"): break
    process_lines.append(line.strip().split())
process_df = pd.DataFrame(np.asarray(process_lines).transpose()[1:], columns=["category", "process_name", "process_num", "rate"])

systematics_lines = []
systematics_lines_length = -1
for line_num, line in enumerate(infile):
    if line.startswith("-"): break
    if systematics_lines_length == -1: systematics_lines_length = len(line.strip().split())
    elif systematics_lines_length != len(line.strip().split()): break
    systematics_lines.append(line.strip().split())
    
systematics_df = pd.DataFrame(np.asarray(systematics_lines), columns=["systematic", "systematic_type"] + list(range(len(process_df))))
if debug: print(systematics_df)
systematics_df_melted = pd.melt(systematics_df, id_vars=["systematic", "systematic_type"], value_vars=list(range(len(process_df))), value_name="systematic_size", var_name="process")
systematics_df_melted = systematics_df_melted[systematics_df_melted["systematic_size"]!="-"]
systematics_df_melted = systematics_df_melted.sort_values(["systematic", "process"]).reset_index(drop=True)
if debug: print(systematics_df_melted)

systematics_with_process_df = pd.merge(systematics_df_melted, process_df, how="left", left_on="process", right_index=True)
if debug: print(systematics_with_process_df.columns)
systematics_with_process_df_unique = systematics_with_process_df[["systematic", "process_num", "process_name", "systematic_type", "systematic_size"]].drop_duplicates()

histogram_collection_df = pd.DataFrame(columns=["process_name", "process_num", "category", "histogram_type", "histogram_bins"])

for index, systematic_row in systematics_with_process_df_unique.iterrows():
    systematic_name = systematic_row["systematic"]
    systematic_type = systematic_row["systematic_type"]
    systematic_size = systematic_row["systematic_size"]
    #if systematic_size == '-': continue
    process_num = systematic_row["process_num"]
    process_name = systematic_row["process_name"]
    print("Starting mountainrange plot for {} in {}".format(systematic_name, process_name))

    search_results_df = systematics_with_process_df[(systematics_with_process_df["systematic"] == systematic_name) & (systematics_with_process_df["process_num"] == process_num)]

    hist_central_array = np.array([])
    hist_upper_array = np.array([])
    hist_lower_array = np.array([])

    separator_locations = []
    label_locations = [0.]
    category_names = []

    if debug: print("Harvesting histograms")
    for search_result_index, search_result_row in search_results_df.iterrows():
        shape_search_df = shape_df[(shape_df["channel"] == search_result_row["category"]) & (shape_df["process"] == search_result_row["process_name"])]
        if len(shape_search_df) == 0:
            #print("No explicit file list for category {} with process {}. Searching for process widcard.".format(search_result_row["category"], search_result_row["process_name"]))
            shape_search_df = shape_df[(shape_df["channel"] == search_result_row["category"]) & (shape_df["process"] == "*")]
            if len(shape_search_df) != 1:
                print("Cannot resolve search result. Skipping.")
                continue
        hist_filename = os.path.dirname(args.filename) + '/' + shape_search_df.iloc[0]["file"]
        hist_name_central = shape_search_df.iloc[0]["central_hist"]
        hist_name_central = re.sub(r"\$PROCESS", process_name, hist_name_central)
        hist_name_central = re.sub(r"\$MASS", str(args.mass), hist_name_central)

        if systematic_type == "shape":
            hist_name_upper = shape_search_df.iloc[0]["syst_hist"]
            hist_name_upper = re.sub(r"\$PROCESS", process_name, hist_name_upper)
            hist_name_upper = re.sub(r"\$SYSTEMATIC", systematic_name + "Up", hist_name_upper)
            hist_name_upper = re.sub(r"\$MASS", str(args.mass), hist_name_upper)

            hist_name_lower = shape_search_df.iloc[0]["syst_hist"]
            hist_name_lower = re.sub(r"\$PROCESS", process_name, hist_name_lower)
            hist_name_lower = re.sub(r"\$SYSTEMATIC", systematic_name + "Down", hist_name_lower)
            hist_name_lower = re.sub(r"\$MASS", str(args.mass), hist_name_lower)
        if debug: print(hist_filename)
        if debug: print(hist_name_central)
        if debug: print(systematic_name)
        if debug: print(systematic_size)
        input_rootfile = pyr.TFile(hist_filename)
        hist_central = input_rootfile.Get(hist_name_central)
        if debug: print(hist_central)
        if systematic_size != '-':
            if search_result_row["systematic_type"] == "shape":
                if debug: print("Reading upper histogram")
                hist_upper = input_rootfile.Get(hist_name_upper)
                if debug: print(hist_upper)
                if debug: print("Reading lower histogram")
                hist_lower = input_rootfile.Get(hist_name_lower)
                if debug: print(hist_lower)
            elif search_result_row["systematic_type"] == "lnN":
                if debug: print("Copying upper histogram (for lnN uncertainty)")
                hist_upper = hist_central.Clone()
                hist_upper.Scale(float(systematic_size))
                if debug: print("Copying lower histogram (for lnN uncertainty)")
                hist_lower = hist_central.Clone()
                hist_lower.Scale(1./float(systematic_size))
        
        if debug: print("Saving central hist")
        hist_central_array = np.concatenate((hist_central_array, root_numpy.hist2array(hist_central)))
        if debug: print("Saving to DataFrame")
        histogram_collection_df.loc[len(histogram_collection_df)] = [process_name, process_num, search_result_row["category"], "central", root_numpy.hist2array(hist_central)]

        if systematic_size != '-':
            separator_locations.append(len(hist_central_array))
            label_locations.append(len(hist_central_array))
            category_names.append(search_result_row["category"])
            if debug: print("Saving upper hist")
            hist_upper_array = np.concatenate((hist_upper_array, root_numpy.hist2array(hist_upper)))
            if debug: print("Saving to DataFrame")
            histogram_collection_df.loc[len(histogram_collection_df)] = [process_name, process_num, search_result_row["category"], systematic_name+"Up", root_numpy.hist2array(hist_upper)]
            if debug: print("Saving lower hist")
            hist_lower_array = np.concatenate((hist_lower_array, root_numpy.hist2array(hist_lower)))
            if debug: print("Saving to DataFrame")
            histogram_collection_df.loc[len(histogram_collection_df)] = [process_name, process_num, search_result_row["category"], systematic_name+"Down", root_numpy.hist2array(hist_upper)]
    
    print("Constructing mountainrange")
    if systematic_size == '-': continue

    mountainrange_central = pyr.TH1F("mountainrange_central", "{} in {}".format(systematic_name, process_name), len(hist_central_array), 0, len(hist_central_array))
    mountainrange_upper = pyr.TH1F("mountainrange_upper", "{} in {}".format(systematic_name, process_name), len(hist_upper_array), 0, len(hist_upper_array))
    mountainrange_lower = pyr.TH1F("mountainrange_lower", "{} in {}".format(systematic_name, process_name), len(hist_lower_array), 0, len(hist_lower_array))

    #for bin_num in range(len(hist_central_array)):
    #    mountainrange_central.SetBinContent(1+bin_num, hist_central_array[bin_num])
    #    mountainrange_upper.SetBinContent(1+bin_num, hist_upper_array[bin_num])
    #    mountainrange_lower.SetBinContent(1+bin_num, hist_lower_array[bin_num])

    #mountainrange_central_ratio = pyr.TH1F("mountainrange_central_ratio", "{} in {}".format(systematic_name, process_name), len(hist_central_array), 0, len(hist_central_array))
    #mountainrange_upper_ratio = pyr.TH1F("mountainrange_upper_ratio", "{} in {}".format(systematic_name, process_name), len(hist_upper_array), 0, len(hist_upper_array))
    #mountainrange_lower_ratio = pyr.TH1F("mountainrange_lower_ratio", "{} in {}".format(systematic_name, process_name), len(hist_lower_array), 0, len(hist_lower_array))

    if args.ratio:
        max_boundary = 10.
        min_boundary = 0.00001
    else:
        max_boundary = max([max(hist_central_array[np.nonzero(hist_central_array)]), max(hist_upper_array[np.nonzero(hist_upper_array)]), max(hist_lower_array[np.nonzero(hist_lower_array)])])*10.
        min_boundary = min([min(hist_central_array[np.nonzero(hist_central_array)]), min(hist_upper_array[np.nonzero(hist_upper_array)]), min(hist_lower_array[np.nonzero(hist_lower_array)])])/10.

    root_numpy.array2hist(hist_central_array, mountainrange_central)
    root_numpy.array2hist(hist_upper_array, mountainrange_upper)
    root_numpy.array2hist(hist_lower_array, mountainrange_lower)
    #mountainrange_central = root_numpy.array2hist(hist_central_array)
    #mountainrange_upper = root_numpy.array2hist(hist_upper_array)
    #mountainrange_lower = root_numpy.array2hist(hist_lower_array)
    #root_numpy.array2hist(hist_central_array, mountainrange_central_ratio)
    #root_numpy.array2hist(hist_upper_array, mountainrange_upper_ratio)
    #root_numpy.array2hist(hist_lower_array, mountainrange_lower_ratio)

    if args.ratio:
        mountainrange_upper.Divide(mountainrange_central)
        mountainrange_lower.Divide(mountainrange_central)
        mountainrange_central.Divide(mountainrange_central)

    if debug: print("Drawing mountainrange")
    canvas = pyr.TCanvas("canvas", "Histogram", 1800, 800)
    #pad1 = pyr.TPad("pad1", "", 0.0, 0.2, 1.0, 1.0)
    #pad2 = pyr.TPad("pad2", "", 0.0, 0.0, 1.0, 0.2)

    #pad1.cd()
    if debug: print("Drawing histograms")
    mountainrange_central.Draw("HIST")
    mountainrange_upper.Draw("HIST SAME")
    mountainrange_lower.Draw("HIST SAME")
    
    if debug: print("Setting line colour")
    mountainrange_central.SetLineColor(1)
    mountainrange_upper.SetLineColor(632)
    mountainrange_lower.SetLineColor(600)

    mountainrange_central.SetStats(0)
    mountainrange_upper.SetStats(0)
    mountainrange_lower.SetStats(0)

    mountainrange_central.GetXaxis().SetTitle("Histogram bins")
    pyr.gPad.SetLogy(True)
    mountainrange_central.GetYaxis().SetRangeUser(min_boundary, max_boundary)
    mountainrange_upper.GetYaxis().SetRangeUser(min_boundary, max_boundary)
    mountainrange_lower.GetYaxis().SetRangeUser(min_boundary, max_boundary)

    if args.ratio: mountainrange_central.GetYaxis().SetTitle("Ratio to central")
    else: mountainrange_central.GetYaxis().SetTitle("Event yields")

    if debug: print(separator_locations)
    separator_lines = []
    for seploc in separator_locations:
        separator_lines.append(pyr.TLine(seploc, min_boundary, seploc, max_boundary))
    for line in separator_lines:
        line.SetLineStyle(pyr.kDashed)
        line.SetLineWidth(2)
        line.SetLineColor(pyr.kBlack)
        line.Draw()

    #print("Drawing legend")
    legend = pyr.TLegend(0.8, 0.8, 0.9, 0.9)
    legend.AddEntry(mountainrange_central, "central")
    legend.AddEntry(mountainrange_upper, "upper")
    legend.AddEntry(mountainrange_lower, "lower")
    legend.Draw()

    #pad2.cd()
    #mountainrange_central_ratio.Draw("HIST")
    #mountainrange_upper_ratio.Draw("HIST SAME")
    #mountainrange_lower_ratio.Draw("HIST SAME")
    
    #print("Setting line colour")
    #mountainrange_central_ratio.SetLineColor(1)
    #mountainrange_upper_ratio.SetLineColor(632)
    #mountainrange_lower_ratio.SetLineColor(600)

    #mountainrange_central_ratio.SetStats(0)
    #mountainrange_upper_ratio.SetStats(0)
    #mountainrange_lower_ratio.SetStats(0)

    #mountainrange_central_ratio.GetXaxis().SetTitle("Histogram bins")
    #mountainrange_central_ratio.GetYaxis().SetRangeUser(0., 2.)
    #mountainrange_upper_ratio.GetYaxis().SetRangeUser(0., 2.)
    #mountainrange_lower_ratio.GetYaxis().SetRangeUser(0., 2.)

    #mountainrange_central_ratio.GetYaxis().SetTitle("Ratio to central")

    #print(separator_locations)
    #separator_lines = []
    #for seploc in separator_locations:
    #    if args.ratio: separator_lines.append(pyr.TLine(seploc, 0., seploc, 2.))
    #    else: separator_lines.append(pyr.TLine(seploc, min_boundary, seploc, max_boundary))
    #for line in separator_lines:
    #    line.SetLineStyle(pyr.kDashed)
    #    line.SetLineWidth(2)
    #    line.SetLineColor(pyr.kBlack)
    #    line.Draw()

    if debug: print("Saving")
    if args.ratio:
        canvas.SaveAs("mountainrange_ratio_{}/{}_{}.pdf".format(args.cardname, systematic_name, process_name))
        canvas.SaveAs("mountainrange_ratio_{}/{}_{}.png".format(args.cardname, systematic_name, process_name))
    else:
        canvas.SaveAs("mountainrange_{}/{}_{}.pdf".format(args.cardname, systematic_name, process_name))
        canvas.SaveAs("mountainrange_{}/{}_{}.png".format(args.cardname, systematic_name, process_name))

    #pad1.cd()
    category_nametags = []
    label_locations.pop(-1)
    for catname, seploc in zip(category_names, label_locations):
        category_nametags.append(pyr.TLatex(seploc, max_boundary/2, catname))
    for cattag in category_nametags:
        cattag.Draw()
    if args.ratio: canvas.SaveAs("mountainrange_ratio_{}/{}_{}.root".format(args.cardname, systematic_name, process_name))
    else: canvas.SaveAs("mountainrange_{}/{}_{}.root".format(args.cardname, systematic_name, process_name))

    del mountainrange_central
    del mountainrange_upper
    del mountainrange_lower

    #del mountainrange_central_ratio
    #del mountainrange_upper_ratio
    #del mountainrange_lower_ratio

histogram_collection_df = histogram_collection_df.iloc[histogram_collection_df[["process_name", "category", "histogram_type"]].drop_duplicates().index].reset_index(drop=True)
print("Finalising")

category_series = histogram_collection_df["category"].drop_duplicates().reset_index(drop=True)
process_series = histogram_collection_df["process_name"].drop_duplicates().reset_index(drop=True)

central_histogram_pd = pd.DataFrame(index = pd.MultiIndex.from_product([category_series, process_series], names=["category", "process"])).reset_index()

for central_histogram_index, central_histogram_row in central_histogram_pd.iterrows():
    shape_search_df = shape_df[(shape_df["channel"] == central_histogram_row["category"]) & (shape_df["process"] == central_histogram_row["process"])]
    if len(shape_search_df) == 0:
        shape_search_df = shape_df[(shape_df["channel"] == central_histogram_row["category"]) & (shape_df["process"] == "*")]
        if len(shape_search_df) != 1:
            print("Cannot resolve search result. Skipping.")
            continue
    hist_filename = os.path.dirname(args.filename) + '/' + shape_search_df.iloc[0]["file"]
    hist_name_central = shape_search_df.iloc[0]["central_hist"]
    hist_name_central = re.sub(r"\$PROCESS", process_name, hist_name_central)
    hist_name_central = re.sub(r"\$MASS", str(args.mass), hist_name_central)
    input_rootfile = pyr.TFile(hist_filename)
    hist_central = input_rootfile.Get(hist_name_central)
    hist_central_array = np.concatenate((hist_central_array, root_numpy.hist2array(hist_central)))
    histogram_collection_df.loc[len(histogram_collection_df)] = [process_name, process_num, search_result_row["category"], "central", root_numpy.hist2array(hist_central)]

histogram_collection_df = histogram_collection_df.iloc[histogram_collection_df[["process_name", "category", "histogram_type"]].drop_duplicates().index].reset_index()
histogram_collection_df.to_pickle("histogram_collection_{}.p".format(args.cardname))