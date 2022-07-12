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
    for shape_profile in ["voigt", "gaussian", "lorentzian"]:
        # if shape_profile != "voigt":
        this_fitter = myfitter.fit_without_constraints(shape_profile=shape_profile)
        this_output_file = output_file.format(suffix=shape_profile + "_without_constraints")
        this_fitter.save_fit_result(this_output_file)
        collected_fitters.append(this_fitter)

        this_output_file = output_file.format(suffix=shape_profile + "_with_constraints")
        this_fitter = myfitter.fit_with_constraints(shape_profile=shape_profile)
        collected_fitters.append(this_fitter)
        this_fitter.save_fit_result(this_output_file)

if __name__ == "__main__":
    main(sys.argv[1:])

