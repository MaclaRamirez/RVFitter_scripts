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

    #  myfitter = RVFitter.load_from_df_file(parsed_args["processed_spectra"])
    fits_to_compare = config["fits_to_compare"]
    output_folder = config["output_folder"]

    object_list = []
    for fit_config in fits_to_compare:
        processed_spectra = fit_config["processed_spectra"]
        fit_file = fit_config["fit_file"]

        myfitter = RVFitter.load_from_df_file(processed_spectra)
        myfitter.load_fit_result(fit_file)
        myfitter.setup_parameters()
        object_list.append(copy.deepcopy(myfitter))

    for fitter in object_list:
        print("Constraints applied:", fitter.constraints_applied)
        for _, row in fitter.df.iterrows():
            parameter = row["parameters"]
            cen = parameter["cen"]
            res = fitter.result.params[cen].value
            print(row["line_name"], row["line_profile"], row["line_hash"], row["date"], res)
        print("")

    __import__('ipdb').set_trace()

if __name__ == "__main__":
    main(sys.argv[1:])
