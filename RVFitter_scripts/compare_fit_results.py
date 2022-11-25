import argparse
import yaml
import sys
import copy
import os
import matplotlib.pyplot as plt
from RVFitter import RVFitter#,  RVFitter_comparison
from RVFitter.RVFitter_comparison import RVFitter_comparison

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
    parser.add_argument(
        '--fit_results',
        type=str,
        required=True,
        help=
        "YAML config to compare fit-results."
    )
    parsed_args = dict(vars(parser.parse_args()))
    with open(parsed_args["fit_results"], 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            raise Exception(exc)

    return parsed_args, config

def main(args):
    parsed_args, config = parse_args(args)

    myfitter = RVFitter.load_from_df_file(parsed_args["processed_spectra"])
    fits_to_compare = config["fits_to_compare"]
    output_folder = config["output_folder"]

    prefix = os.path.basename(parsed_args["processed_spectra"]).replace(".pkl", "")
    dirnames = [os.path.dirname(fit_file) for fit_file in fits_to_compare]
    dirnames = list(set(dirnames))
    if len(dirnames) != 1:
        raise Exception("All fits must be in the same directory.")
    output_dir = dirnames[0]

    object_list = []
    for fit_file in fits_to_compare:
        myfitter.load_fit_result(fit_file)
        myfitter.setup_parameters()
        object_list.append(copy.deepcopy(myfitter))
    comparer = RVFitter_comparison(object_list, output_folder)
    comparer.create_overview_df()

    for variable in ["cen"]:
         # no plots for "amp" and "sig" as they are not comparable between lines
         comparer.compare_fit_results_1D(variable=variable, suffix=suffix)
    
    figname = os.path.join(
                output_folder, prefix + "_fits_and_residuals.png")
    comparer.plot_fits_and_residuals(figname=figname)
    plt.show()
    
    for variable in ["cen"]:#, "amp", "sig"]:
        comparer.write_overview_table(variable=variable)

if __name__ == "__main__":
    main(sys.argv[1:])
