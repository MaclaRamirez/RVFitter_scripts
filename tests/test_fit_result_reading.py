import unittest

import numpy as np

from RVFitter import RVFitter
import pkg_resources
import os
import copy
import sys
import matplotlib.pyplot as plt
import pandas as pd

def get_tmp_file(filename):
    with open(filename, "r") as f:
        data = f.read()

    tmp_specsfilelist = filename.replace(".txt", "_tmp.txt")
    with open(tmp_specsfilelist, "w") as f:
        for fileline in data.splitlines():
            tmp_dir = os.path.dirname(filename)
            tmp_data = os.path.join(tmp_dir, fileline)
            f.write(tmp_data + "\n")
    return tmp_specsfilelist


class TestRVFitter(unittest.TestCase):
    def setUp(self):
        line_list = pkg_resources.resource_filename(
            "RVFitter",
            "tests/test_data/debug_spectral_lines_RVmeasurement.txt")
        self.specsfilelist = pkg_resources.resource_filename(
            "RVFitter", "tests/test_data/debug_specfile_list.txt")

        tmp_specsfilelist = get_tmp_file(self.specsfilelist)

        myfitter = RVFitter.from_specsfilelist_name_flexi(
            specsfilelist_name=tmp_specsfilelist,
            line_list=line_list)
        self.myfitter = myfitter

    #  def test_fitting(self):
    #
    #      filename = os.path.join(os.path.dirname(self.specsfilelist),
    #                              "B275_speclist.pkl")
    #      self.myfitter = RVFitter.load_from_df_file(filename=filename)
    #
    #      collected_fitters = []
    #      for shape_profile in ["gaussian", "lorentzian"]:
    #          this_fitter = fit_without_constraints(self.myfitter, shape_profile=shape_profile)
    #          this_fitter.label = shape_profile + " without constraints"
    #          collected_fitters.append(this_fitter)
    #
    #          this_fitter = fit_with_constraints(self.myfitter, shape_profile=shape_profile)
    #          this_fitter.label = shape_profile + " with constraints"
    #          collected_fitters.append(this_fitter)
    #
    #      color_dict = {0: "red", 1: "blue", 2: "green", 3: "orange"}
    #      fig, axes = self.myfitter.get_fig_and_axes()
    #      for idx, this_fitter in enumerate(collected_fitters):
    #          if idx == 0:
    #              this_fitter.plot_data(fig=fig, axes=axes)
    #          this_fitter.plot_fit(fig=fig, axes=axes, plot_dict={"zorder": 2.5, "color": color_dict[idx], "label": this_fitter.label})
    #      handles, labels = axes[-1, -1].get_legend_handles_labels()
    #      fig.legend(handles, labels, ncol=2, loc='lower center')
    #      fig.savefig("test_constraints.pdf")
    #      # plt.show()
    #
    #      #Compare RV measurements
    #      dates = self.myfitter.df['date'].unique()
    #      date = dates[0]
    #      # line_hashes = this_epoch['line_hash'].unique()
    #      fig, ax = plt.subplots()
    #      colors = ["red", "blue", "green"]
    #      for j, date in enumerate(dates):
    #          this_epoch = self.myfitter.df.query("date == '{}'".format(date))
    #          my_dict = {}
    #          for fitter in collected_fitters:
    #              print(fitter.label)
    #              my_dict[fitter.label] = {}
    #              for idx, row in this_epoch.iterrows():
    #                  hash = row['line_hash']
    #                  variable_formater = 'cen_line_{}_epoch_{}'.format(hash, date)
    #                  print(variable_formater)
    #                  print(fitter.result.params[variable_formater].value)
    #                  my_dict[fitter.label][hash] = fitter.result.params[variable_formater].value
    #
    #
    #          ax.plot(list(my_dict['gaussian without constraints'].values()),
    #                  list(my_dict['lorentzian without constraints'].values()), 'o', color=colors[j])#, label='without constraints')
    #          ax.plot(np.linspace(-1, 60, 1000), np.linspace(-1, 60, 1000), '--',  color=colors[j], label=date)
    #          ax.plot(list(my_dict['gaussian with constraints'].values()),
    #                  list(my_dict['lorentzian with constraints'].values()), 'o',  color=colors[j], alpha=0.05)# label='with constraints')
    #
    #      ax.legend()
    #      #  plt.show()


    def test_table(self):

        filename = os.path.join(os.path.dirname(self.specsfilelist),
                                "B275_speclist.pkl")
        self.myfitter = RVFitter.load_from_df_file(filename=filename)

        collected_fitters = {}
        for shape_profile in ["gaussian", "lorentzian"]:
            this_fitter = fit_without_constraints(self.myfitter, shape_profile=shape_profile)
            this_fitter.label = shape_profile + " without constraints"
            collected_fitters["no_constraint_" + shape_profile] = this_fitter

            this_fitter = fit_with_constraints(self.myfitter, shape_profile=shape_profile)
            this_fitter.label = shape_profile + " with constraints"
            collected_fitters["constraint_" + shape_profile] = this_fitter




        # create dataframe from collected_fitters
        df = pd.DataFrame()
        for key, fitter in collected_fitters.items():
            df = add_parameters_to_df(df=df, fitter=fitter, suffix="_" + key)
        print(df)

def add_parameters_to_df(df, fitter, suffix):
    l_amp = []
    l_cen = []
    l_sig = []
    l_error_amp = []
    l_error_cen = []
    l_error_sig = []
    for parameter in fitter.df["parameters"]:
        amp = parameter["amp"]
        cen = parameter["cen"]
        sig = parameter["sig"]

        l_amp.append(fitter.result.params[amp].value)
        l_cen.append(fitter.result.params[cen].value)
        l_sig.append(fitter.result.params[sig].value)

        l_error_amp.append(fitter.result.params[amp].stderr)
        l_error_cen.append(fitter.result.params[cen].stderr)
        l_error_sig.append(fitter.result.params[sig].stderr)

    this_df = pd.DataFrame({"date" + suffix: fitter.df["date"],
        "line_profile" + suffix: fitter.df["line_profile"],
        "amp" + suffix: l_amp,
        "cen" + suffix: l_cen,
        "sig" + suffix: l_sig,
        "error_amp" + suffix: l_error_amp,
        "error_cen" + suffix: l_error_cen,
        "error_sig" + suffix: l_error_sig
        })

    for col in df.columns:
        if "date" in col:
            check = df[col] == this_df["date" + suffix]
            if np.sum(check) != len(check):
                raise Exception("Error: date mismatch")
        if "line_profile" in col:
            check = df[col] == this_df["line_profile" + suffix]
            if np.sum(check) != len(check):
                raise Exception("Error: line_profile mismatch")
    df = pd.concat([df, this_df], axis=1)
    return df



def fit_with_constraints(myfitter, shape_profile="gaussian"):
    # prepare fitting
    this_fitter = copy.deepcopy(myfitter)
    this_fitter.shape_profile = shape_profile
    this_fitter.constrain_parameters(group="cen", constraint_type="epoch")
    this_fitter.constrain_parameters(group="amp", constraint_type="line_profile")
    this_fitter.constrain_parameters(group="sig", constraint_type="line_profile")
    this_fitter.run_fit()
    return this_fitter

def fit_without_constraints(myfitter, shape_profile="gaussian"):
    # prepare fitting
    this_fitter = copy.deepcopy(myfitter)
    this_fitter.shape_profile = shape_profile
    this_fitter.run_fit()
    return this_fitter

if __name__ == "__main__":
    unittest.main()
