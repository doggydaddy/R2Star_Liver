# R2Star_Liver
liver r2star iron quantification for central radiology, karolinska university hospital, solna.
Date: 2017-02-28
Author: Yanlu Wang

Assumes input files to be 5 datasets (magnitude images only). Each dataset has 3 TEs and 5 slices. 15 TEs in total are interleaved between the 5 datasets. for "ultra-rapid" multi-gradient-echo acquisition. The script works best on ABNORMAL liver with high iron concentration and might not work well for normal tissue. 

r2star.py -- main file, run to see usage: <r2star.py> -i <input directory> -o <output directory>

r2star_viewData.py -- script only for visualizing input data. Note the input data directory is hard-coded into the script

r2star_modelInspector.py -- script for visualzing the model used for fitting

r2star_fittingInpsector.py -- script for visualizing fitting results, acutally reads the data according to r2star_viewData.py and does the fitting, so this might take some time. Like r2star_viewData.py, the input data directory is hard-coded into the script.
