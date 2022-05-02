from RVFitter import RVFitter
#  import pkg_resources
import argparse
import os
import copy
import sys
import matplotlib.pyplot as plt

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

def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Execute the fitting function from RVFitter')
    parser.add_argument('--debug', action='store_true', help="Debug mode")
    parser.add_argument(
        '--processed_spectra',
        type=str,
        required=True,
        help=
        "Path to the pickle-file which holds the dataframe."
    )
    return dict(vars(parser.parse_args()))

def main(args):
    parsed_args = parse_args(args)

    myfitter = RVFitter.load_from_df_file(parsed_args["processed_spectra"])

    collected_fitters = []
    for shape_profile in ["gaussian", "lorentzian"]:
        this_fitter = fit_without_constraints(myfitter, shape_profile=shape_profile)
        this_fitter.label = shape_profile + " without constraints"
        collected_fitters.append(this_fitter)
        #  myfitter.plot_fit()
        #  myfitter.plot_residuals()
        #  myfitter.plot_fit_with_residuals()

        this_fitter = fit_with_constraints(myfitter, shape_profile=shape_profile)
        this_fitter.label = shape_profile + " with constraints"
        collected_fitters.append(this_fitter)

    color_dict = {0: "red", 1: "blue", 2: "green", 3: "orange"}
    fig, axes = myfitter.get_fig_and_axes()
    for idx, this_fitter in enumerate(collected_fitters):
        if idx == 0:
            this_fitter.plot_data(fig=fig, axes=axes)
        this_fitter.plot_fit(fig=fig, axes=axes, plot_dict={"zorder": 2.5, "color": color_dict[idx], "label": this_fitter.label})
    handles, labels = axes[-1, -1].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=2, loc='lower center')
    fig.savefig("test_constraints.pdf")
    plt.show()


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
    main(sys.argv[1:])

