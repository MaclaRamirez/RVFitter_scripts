from RVFitter import RVFitter
from RVFitter import utils
#  import pkg_resources
import argparse
import os
import sys
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


def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Execute the fitting function from RVFitter')
    parser.add_argument('--debug', action='store_true', help="Debug mode")
    parser.add_argument(
        '--processed_spectra',
        type=str,
        required=True,
        help="Path to the pickle-file which holds the dataframe.")
    parser.add_argument(
        '--line_list',
        type=str,
        required=False,
        help="Path to the line-list which should be used in the fit.",
        default=None)
    return dict(vars(parser.parse_args()))


def main(args):
    parsed_args = parse_args(args)
    if parsed_args["line_list"] is not None:
        print("Reading line list from file.")
        line_list = utils.read_line_list(parsed_args["line_list"])
    else:
        line_list = None

    df = pd.read_pickle(parsed_args["processed_spectra"])

    skimmed_df = utils.manipulate_df_by_line_list(df, line_list)

    myfitter = RVFitter.load_from_df(skimmed_df)

    if line_list is not None:
        suffix = os.path.basename(parsed_args["line_list"]).replace(".txt", "")

        directory, filename = os.path.split(parsed_args["processed_spectra"])
        output_processed_df = os.path.join(directory, filename.split(".")[0] + "_" + suffix + ".pkl")

        info = f"""
        Writing out dataframe with hash for reproduction.
        Hash: {suffix}
        Line list: {line_list}
        Output file: {output_processed_df}
        """
        print(info)
        logfile = os.path.join(directory, "log_" + suffix + ".txt")
        with open(logfile, "w") as f:
            f.write(info)
        myfitter.save_df(output_processed_df)
    else:
        suffix = ""

    output_file = parsed_args["processed_spectra"].replace(
        ".pkl", "_{suffix}.pkl")

    for shape_profile in ["gaussian", "lorentzian"]:#"voigt", 
        # if shape_profile != "voigt":
        this_fitter = myfitter.fit_without_constraints(
            shape_profile=shape_profile)
        if suffix != "":
            this_output_file = output_file.format(suffix=shape_profile +
                                                  "_without_constraints_" + suffix)
        else:
            this_output_file = output_file.format(suffix=shape_profile +
                                                  "_without_constraints")
        this_fitter = myfitter.fit_without_constraints(
            shape_profile=shape_profile)
        this_fitter.save_fit_result(this_output_file)

        if suffix != "":
            this_output_file = output_file.format(suffix=shape_profile +
                                                  "_with_constraints_" + suffix)
        else:
            this_output_file = output_file.format(suffix=shape_profile +
                                                  "_with_constraints")
        this_fitter = myfitter.fit_with_constraints(
            shape_profile=shape_profile)
        this_fitter.save_fit_result(this_output_file)


if __name__ == "__main__":
    main(sys.argv[1:])
