# Sensory-Gating-Analysis
MATLAB scripts for fitting cumulative Gaussian psychometric functions to behavioral data, extracting PSE and JND

## Psychometric Function Fitting (MATLAB)

This repository provides MATLAB scripts to fit cumulative Gaussian psychometric functions to behavioral data collected from active, passive, control, and sensory attenuation (SA) conditions.  
It extracts the **Point of Subjective Equality (PSE)** and **Just Noticeable Difference (JND)** for each condition, with optional computation of **Weber fraction** and **McFadden’s R²** as measures of goodness of fit.

### Overview
The script loads participant CSV files, merges data across blocks, filters trials based on latency and duration quality criteria, and fits psychometric curves using the `anamax` fitting function (cumulative Gaussian model).

By default, the output table contains PSE and JND values per participant and condition.  
Lines computing Weber fraction and McFadden’s R² are included but commented out for flexibility and clarity.
 
### Requirements
- MATLAB R2021a or later  
- The `anamax.m` function available in your MATLAB path  
