# RVFitter_scripts

## Installation

Clone the repository and install the package in development mode:

```
pip install -r requirements.txt
```

## Idea of the project
Scripts to run RVFitter (https://github.com/WolfXeHD/RVFitter) in order to normalize and clip spectra, fit spectral lines, plot fitting results and compare the results with different fitting functions.

## Usage

RVFitter_scripts has curently 3 steps:

### 1. runRVFitter.py:
 Look trough the spectral lines and choose which ones are to be fitted, then re-normalize and clip each line. This is done for all spectra corresponding to each epoch for a single star. The lines are chosen for the first spectrum, for the consequent epochs only the re-normalization and clipping is possible.

INPUT:

  * `--specfile_list`: txt file containing the path to the spectra corresponding to each epoch. The spectra should be text files containing two (or three) columns with the wavelength in nm and the normalized flux (and error).

  * `--output_file`: path and name of the pkl file to be written by runRVFitter containing the chosen, re-normalized and clipped lines

  * `--line_list`: txt file containing three columns: name; line_center (in Anstroms); plotting_window for each of the lines that one wants to look through

OUTPUT:

* `output_file`: pkl file containing the re-normalized and clipped spectra

EXAMPLE:

```python runRVfitter.py --specfile_list specfile_list.txt --line_list configs/metalic_lines.txt --output_file B111_processed_spectra.pkl```

where **specfile_list.txt** looks like:
```
/data/normalized/B111_UVB_20190411T09.nspec
/data/normalized/B111_UVB_20190413T08.nspec
/data/normalizedB111_UVB_20190605T05.nspec
```

and **metalic_lines.txt**:

```NII 3995 8
HeI 5015.6783 8
HeII 5411.53 12
HeI 5875.621 20
HeI 6678.151 20
```

### 2. execute_fit.py:

Fits Gaussian, Lorenzian and Voigt profiles for all the selected lines. It executes two kinds of fits for each function:

- without constrains: Each line is fit individually, resulting in an amplitude, sigma and centroid for each lines

- with constrains: All lines for a single epoch are forced to have the same velocity, while a given line for all epochs is forced to have the same amplitude and sigma.

INPUT:

- **--processed_spectra**: pkl file containing the re-normalized and clipped spectra

OUTPUT:

- **fit_results**: one pkl file per fit executed

EXAMPLE:

```
python execute_Fit.py --processed_spectra B111_processed_spectra.pkl
```

### 3. compare_fit_results.py:

Makes two plots, one with an overview of the fitted lines and one showing the results for the 6 different fits

INPUT:

- **--processed_spectra** from step 1.

- **--fit_results**: a yaml file containing the paths to the output folder where the plots and tables need to be saved and to the fit results from step 2.

OUTPUT:

- **results_table.tex**: table summarizing the fitting results. When the fits are done with constraints (v_{.wc}), the errors are those given by lmfit. When the fir is done without constrains (v_{.nc}) the velocity given is the average between all lines used and the error corresponds to the standard deviation.

- **fits_and_residuals.png**: Plot showing the fit results for all lines

- **compare_results_velocity.png**: Plot showing the velocity measurements in all epochs for a single star

EXAMPLE:

```
python compare_fit_results.py --processed_spectra B111_processed_spectra.pkl --fit_results fits_to_compare.yaml
```

where **fits_to_compare.yaml** looks like:

```
output_folder: RVFitter_results/B111/ <br>
fits_to_compare:
  - B111_processed_spectra_gaussian_with_constraints.pkl
  - B111_processed_spectra_gaussian_without_constraints.pkl
  - B111_processed_spectra_lorentzian_with_constraints.pkl
  - B111_processed_spectra_lorentzian_without_constraints.pkl
  - B111_processed_spectra_voigt_with_constraints.pkl
  - B111_processed_spectra_voigt_without_constraints.pkl
```
