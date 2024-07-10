# FCR-GLM-AED-surrogate-runs

Repository to hold model files and output for multiple model runs of the General Lake Model - Aquatic EcoDynamics (GLM-AED) coupled hydrodynamic-water quality model.

Model runs are to inform development of a surrogate model to assist with GLM-AED model calibration and forecasting.

Each of 1000 model runs is parameterized with a different combination of phytoplankton parameters.

## Guide to folders and files

1) **Data:** contains chlorophyll-a observations against which to compare model output as well as parameters and parameter ranges that are altered across model runs
2) **Scripts:** contains R scripts to run 1000 iterations of GLM-AED with different parameter combinations, as well as helper scripts to wrangle data, etc.
3) **GLM-AED:** contains model input, output, and namelist files for the GLM-AED model
