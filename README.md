# NalaFit

These are scripts for processing NMR relaxation peak intensities into relaxation time constants (i.e. T1 or T2 for amides)
This is the successor to NalaFit running with python 2.7

Errors in peak intensities for a time point are estimated by measuring the median squared difference between 
corresponding peak heights from a duplicate. A duplicate is usually collected at the first, middle and last 
(or near last) points (3 duplicates). Thus an error in peak height is collected as a function of relaxation delay. Errors
for non duplicated time points are estimated by interpolation.

This is an improvement over the method presented in: Skelton NJ, Palmer AG 3rd, Akke M, Kordel J, Rance M, and Chazin WJ JMR B 102 1993
which can lead to distortion due to a few large outlying differences
### Future Features

I may include support for things like exponential fits with baseline or biexponential fits should the need arise


## Installing NalaFit


### Create a Conda Environment (or use an existing)
```bash
    conda create --name NalaFit python=3.9
```

### Activate Your Conda Environment

```bash
    conda activate NalaFit #Or whatever conda environment you want to install it in
```

### Install NalaFit
```bash
    pip install git+https://github.com/anthony-bishop-tamu/NalaFit.git@main
```
Be sure your target conda environment is active!
pip should automatically install the following dependencies
- numpy>=1.24.0,<2 #numpy < 2 is enforced out of caution due to compatibility issues of dependencies on certain platforms, 
If all dependencies support numpy 2.0 on your platform, you may upgrade
- pandas>=2.0,
- matplotlib>=3.2,
- scipy >= 1.10

### Test NalaFit
Do a quick test by running in your terminal.
```bash
   NalaFit
```
You should see an error message with usage instructions
## NalaFit usage
```bash
    NalaFit --input relaxation_data.txt --output output_folder
```
where relaxation_data.txt and output_folder are the data source and destination respectively. 

### NalaFit Input
In TestDataSets there are two example input files. NalaFit has no requirement as to the number of data points in your
relaxation curve except that at least two time points were collected in duplicate.

The input file is just whitespace delimited (any number of tabs and/or spaces will work)

1st column is just an arbitrary name for each relaxation site (usually something like S2N-H)
Each other column is the headed by the relaxation delay (In seconds) with each entry being a _raw_ peak height
Time points may appear in any order your choose. 

## NalaFit Outputs
### NalaFit outputs the following files

#### errorValues.csv - This is just a list of what errors were used for each time point
#### fit_Chi2_distribution.png - This the distribution of the chi2 values for each fit. It is overlayed with the theoretical distribution
#### fittedTaus.txt - This contains the mean and standard deviation for the time constant (obtained by monte carlo). These are the values you will use downstream
#### runParameters.txt - A quick summary of the run parameters, containing most importantly the version of the software

### TwoParameterFit Directory
This directory contains all of the two parameter fit plots you will find a folder for each relaxation
site. These folders contain the histograms of the fitted parameters from the Monte Carlo run overlayed with 
a Gaussian distribution of the same mean and variance as simulation results

You will also find plots of the fit for each relaxation site, with a thin strip plot indicating residuals
Data are normalized prior to fitting


    

