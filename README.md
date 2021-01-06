# SVAR
MATLAB code for quick n dirty structural VAR modelling to check impulse response functions

The main code is called TEMPLATE. Other files are functions called in the main code. Code is written for Matlab, should easily be convertable to Python (using python packages for cointegration)

Section 1 sets up the data, using a txt file.

Section 2 plots the data, often useful when the type of (S)VAR/Cointegration model to be used is intuitively obvious.

Section 3 conducts ADF tests for unit roots (work in progress).

Section 4 conducts tests on a regular VAR for optimal leg length choice.

Section 5 Choose your lag length at the top using economic intuition and discretion, then it plots/prints output for you

Section 6-8 are tests for misspecification on your VAR model

Section 9-10 uses the matlab built-in Johansen tests for cointegration and VEC modelling.

Section 11-12 conducts tests on the chosen VECM

Structural VAR modelling (work in progress):

Section 13-15 sets up three typical types of structural VAR models, Cholesky, common trends and long run restrictions. Economic understanding is needed to impose proper restrictions in the restrictions function attached.

Section 16 conducts bootstrapping

Section 17 plots impulse response functions with confidence intervals

Section 18 prints FEVD results
