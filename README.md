# Gene-environment interactions in COVID-19 patients
code necessary to reproduce the main analyses published in:

Widespread gene-environment interactions shape the immune response to SARS-CoV-2 infection in hospitalized COVID-19 patients (2025)

HE Randolph, R Aguirre-Gamboa, E Brunet-Ratnasingham, T Nakanishi, V Locher, E Ketter, C Brandolino, C Larochelle, A Prat, N Arbour, A Dumaine, A Finzi, M Durand, JB Richards, DE Kaufmann, and LB Barreiro

# General dependencies
R (tested in versions 4.1.0 and 4.3.1)

gcc or LLMV (tested with Apple clang version 14.0.3 (clang-1403.0.22.14.1))

# Usage
1. place the folder `inputs`, available for download on Zenodo (10.5281/zenodo.10928039), in the desired working directory
2. place the folders `main_analyses` and `common_functions` in the same working directory
3. check individual dependencies at the header of each script stored in `main_analyses`, and install the necessary CRAN and Bioconductor packages
4. change the working directory in the script to your desired working directory
```
current = [desired working directory]
```
5. run the scripts stored in `main_analyses` in R -- if `main_analyses` and `inputs` are located in the same directory, results will populate in a folder named `outputs`