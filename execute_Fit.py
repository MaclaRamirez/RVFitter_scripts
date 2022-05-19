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

    output_file = parsed_args["processed_spectra"].replace(".pkl", "_{suffix}.pkl")

    collected_fitters = []
    for shape_profile in ["gaussian", "lorentzian"]:
        this_fitter = myfitter.fit_without_constraints(shape_profile=shape_profile)
        this_output_file = output_file.format(suffix=shape_profile + "_without_constraints")
        this_fitter.save_fit_result(this_output_file)
        collected_fitters.append(this_fitter)

        this_output_file = output_file.format(suffix=shape_profile + "_with_constraints")
        this_fitter = myfitter.fit_with_constraints(shape_profile=shape_profile)
        collected_fitters.append(this_fitter)
        this_fitter.save_fit_result(this_output_file)

    #  color_dict = {0: "red", 1: "blue", 2: "green", 3: "orange"}
    #  fig, axes = myfitter.get_fig_and_axes()
    #  for idx, this_fitter in enumerate(collected_fitters):
    #      if idx == 0:
    #          this_fitter.plot_data(fig=fig, axes=axes)
    #      this_fitter.plot_fit(fig=fig, axes=axes, plot_dict={"zorder": 2.5, "color": color_dict[idx], "label": this_fitter.label})
    #  handles, labels = axes[-1, -1].get_legend_handles_labels()
    #  fig.legend(handles, labels, ncol=2, loc='lower center')
    #  fig.savefig("test_constraints.pdf")
    #  plt.show()

    fig, ax_dict = this_fitter.get_fig_and_ax_dict()

    color_dict = {0: "red", 1: "blue", 2: "green", 3: "orange"}
    for idx, this_fitter in enumerate(collected_fitters):
        if idx == 0:
            this_fitter.plot_data_and_residuals(fig=fig, ax_dict=ax_dict)
        this_fitter.plot_fit_and_residuals(fig=fig,
                                           ax_dict=ax_dict,
                                           add_legend_label=True,
                                           plot_dict={"zorder": 2.5,
                                                      "color": color_dict[idx],
                                                      "markersize": "1",
                                                      },
                                           plot_dict_res={"color": color_dict[idx],
                                                          "marker": ".",
                                                          "linestyle": "None",
                                                          "markersize": "2"})
    handles, labels = ax_dict[list(ax_dict.items())[0][0]].get_legend_handles_labels()
    labels = [this_fitter.label for this_fitter in collected_fitters]
    fig.legend(handles, labels, ncol=2, loc='lower center')
    plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])

