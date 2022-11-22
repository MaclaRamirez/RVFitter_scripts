import copy
import glob
import os
import pandas as pd
import argparse
import sys

from RVFitter import RVFitter, RVFitter_comparison

def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Script to aggregate results of different stars')
    parser.add_argument('--debug', action='store_true', help="Debug mode")
    parser.add_argument(
        '--pattern',
        type=str,
        required=True,
        help=
        "pattern handed to glob.glob to look for processed_spectra."
    )
    parser.add_argument('--add_spectra_search_suffix', action='store_true', help="Debug mode")
    parser.add_argument(
        '--tablename',
        type=str,
        required=True,
        help=
        "Name of the table file written to disk."
    )

    parsed_args = dict(vars(parser.parse_args()))
    if parsed_args["add_spectra_search_suffix"]:
        parsed_args["pattern"] = os.path.join(parsed_args["pattern"], "*/*_processed_spectra.pkl")

    return parsed_args

def get_available_fit_results(processed_spectrum):
    dirname = os.path.dirname(processed_spectrum)
    pattern = os.path.join(dirname, "*_processed_spectra_*.pkl")
    available_fits = glob.glob(pattern)
    return available_fits

def main(args):
    parsed_args = parse_args(args)
    processed_spectra = glob.glob(parsed_args["pattern"])

    print(processed_spectra)
    print("Found {n} processed spectra".format(n=len(processed_spectra)))
    # processed_spectra = glob.glob('/Users/ramirez/Dropbox/Projects/M17_binaries/RVFitter_results/*/')

    list_of_dfs = []
    for processed_spectrum in processed_spectra:
        try:
            l_objects = []
            print("loading:", processed_spectrum)
            myfitter = RVFitter.load_from_df_file(processed_spectrum)

            available_fits = get_available_fit_results(processed_spectrum=processed_spectrum)

            for available_fit in available_fits:
                myfitter.load_fit_result(available_fit)
                myfitter.setup_parameters()
                l_objects.append(copy.deepcopy(myfitter))
        except:
            print("\tfailed loading:", processed_spectrum)
            continue
        this_comparer = RVFitter_comparison(l_objects)
        this_comparer.create_overview_df()
        this_df = this_comparer.get_df_for_latex(variable="cen", add_starname=True)
        list_of_dfs.append(this_df)

    master_df = pd.concat(list_of_dfs)
    table_name = parsed_args["tablename"]
    RVFitter_comparison.write_df_to_table(input_df=master_df,
                                          filename=table_name)

if __name__ == "__main__":
    main(sys.argv[1:])