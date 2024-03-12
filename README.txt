Release Notes:
The file format has changed, Masfit now takes unnormalized data. 
Mas fit also outputs the distribution of chi2 values for each fit, against the expected distribution. This is a way of assessing the quality of the errors

call Masfit in the following way

python3 Masfit.py data_file.txt False

Use "False" if you do not want to try 3-parameter fits. There is generally little reason to try 3-parameter fits. I recommend this mode for virtually all cases
Use "True" if you want to perform BIC model selection between 2 and 3 parameter fits. Masfit will automatically choose which model to use the tau values from in the final output file


Masfit has been tested on python3.8 and requires the following packages, I recommend creating a conda virtual environment for Masfit

matplotlib
scipy
numpy
pandas
