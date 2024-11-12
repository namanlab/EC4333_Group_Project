# EC4333 Group Project: Yield Curve Modeling in Negative Interest Rate Environments

## Abstract

The persistence of negative interest rates in several economies, especially in the Eurozone and Japan, has presented significant challenges for traditional financial modeling frameworks. Negative rates disrupt conventional models that assume positive rates, increasing uncertainty in projecting the term structure of interest rates. Understanding and accurately modeling these negative interest rates is essential, as it directly impacts bond pricing, investment decisions, and broader economic forecasts.

This project investigates the effectiveness of the Nelson-Siegel-Svensson (NSS) model, the Nelson-Siegel (NS) model, and cubic splines in fitting yield curves across environments where interest rates may be negative, positive, or close to zero. These models are widely used by central banks due to their adaptability and precision in capturing diverse yield curve shapes, particularly in negative interest rate regimes. Our primary objectives are to replicate key findings from academic studies on yield curve modeling under negative rates, evaluate model performance, and apply these models to real-world data from the Eurozone.

Through this project, we aim to understand how well these models perform in negative interest rate environments, identify the challenges of dynamic fitting, and assess whether the prevailing interest rate environment influences the quality of model fit.

---

## Project Structure

Below is an overview of the project repository's structure, along with detailed descriptions of each folder and file.

```plaintext
.
├── Cubic_Splines.R               # Code for applying cubic splines on Eurozone data, including analysis on knots
├── NS_Model.R                    # Code for NS and NSS models on Eurozone data
├── Process_Data.R                # Code for processing Eurozone data
├── README.md                     # Project documentation
├── app                           # Directory containing code for a web application to explore yield curve shapes
│   ├── df_spot_rates.csv         # Sample data for the web application
│   ├── rsconnect                 # Deployment files for the Shiny app
│   │   └── documents
│   │       └── simple_app.R      # Shiny app script
│   │           └── shinyapps.io
│   │               └── 3vxc23-naman-agrawal
│   │                   └── NS_Sample_App_EC4333.dcf
│   └── simple_app.R              # Main R script for the web app
├── categorisation.R              # Code for categorizing Eurozone yield curves by rate environment
├── data                          # Directory containing all datasets
│   ├── data_aaa.csv              # Eurozone data for AAA bonds
│   ├── data_all.csv              # Eurozone data for all bonds
│   ├── df_results_exns.csv       # NSS results across all time periods for Eurozone
│   ├── df_results_exns_params.csv # NSS model parameters across all time periods for Eurozone
│   ├── df_results_ns.csv         # NS results across all time periods for Eurozone
│   ├── df_results_ns_params.csv  # NS model parameters across all time periods for Eurozone
│   ├── df_results_spline.csv     # Spline model results across all time periods for Eurozone
│   ├── df_spot_rates_aaa.csv     # Cleaned Eurozone data for AAA bonds
│   ├── df_spot_rates_all.csv     # Cleaned Eurozone data for all bonds
│   ├── japan_df_results_exns.csv # NSS results across all time periods for Japan
│   ├── japan_df_results_exns_params.csv # NSS model parameters across all time periods for Japan
│   ├── japan_df_results_ns.csv   # NS results across all time periods for Japan
│   ├── japan_df_results_ns_params.csv # NS model parameters across all time periods for Japan
│   ├── japan_df_results_spline.csv # Spline model results across all time periods for Japan
│   └── toy_data.csv              # Toy data based on Eurozone data for model testing
├── imgs                          # Directory for Eurozone results and graphs
│   ├── NSSException_1.png        # Exception plot in NSS model for Eurozone
│   ├── NSSException_2.png
│   ├── abnormal_eg1.png
│   ├── abormal_eg_2.png
│   ├── adjr2_ncs.png             # Adjusted R2 plot for cubic spline model (Eurozone)
│   ├── adjr2_ncs_bp.png
│   ├── adjr2_ncs_heatmap.png
│   ├── adjr2_nss.png             # Adjusted R2 plot for NSS model (Eurozone)
│   ├── adjr2_nss_bp.png
│   ├── beta0nschange.png         # Change in beta0 parameter for NSS model
│   ├── beta10nschange.png        # Change in beta1 parameter for NSS model
│   ├── beta1nschange.png
│   ├── categorisation_plot.png   # Yield curve categorization plot
│   ├── model_comp_adjr2.png      # Model comparison plot (Eurozone)
│   ├── model_error_comparison.png # Error comparison between models
│   ├── test_mse_heatmap.png
│   ├── test_mse_ncs.png          # MSE heatmap for cubic spline model (Eurozone)
│   ├── test_mse_ncs_bp.png
│   ├── test_mse_nss.png          # MSE heatmap for NSS model (Eurozone)
│   └── test_mse_nss_bp.png
├── japan_results                 # Directory for analysis on Japanese data
│   ├── Cubic_Splines.R           # Code for cubic splines on Japanese data, including knot analysis
│   ├── NS_Model.R                # Code for NS and NSS models on Japanese data
│   ├── categorisation.R          # Code for categorizing Japanese yield curves by rate environment
│   ├── df_japan_spot_rates.csv   # Cleaned Japanese data for spot rates
│   ├── imgs                      # Directory for Japanese results and graphs
│   │   ├── adj22_ns_bp.png
│   │   ├── adj_ncs.png           # Adjusted R2 plot for cubic spline model (Japan)
│   │   ├── adj_ncs_knots.png
│   │   ├── adjr2_ns.png          # Adjusted R2 plot for NS model (Japan)
│   │   ├── beta0.png             # Beta parameter plot for NS/NSS model (Japan)
│   │   ├── beta1.png
│   │   ├── beta10.png
│   │   ├── categorisation.png     # Categorization plot (Japan)
│   │   ├── test_mse_ncs.png      # MSE heatmap for cubic spline model (Japan)
│   │   ├── test_mse_ncs_knots.png
│   │   ├── test_mse_ns.png       # MSE heatmap for NS model (Japan)
│   │   └── test_mse_ns_bp.png
│   └── model_compare.R           # Code to compare NS, NSS, and cubic spline models in Japan
├── model_compare.R               # Code to compare NS, NSS, and cubic spline models in Eurozone
├── old                           # Directory for older, alternative models
│   ├── Varsieck_Model.R          # Code for the Varsieck model (initially tested)
│   ├── model_cir.R               # Code for Cox-Ingersoll-Ross (CIR) model
│   └── other                     # Additional references and notes
│       ├── EC4333_Proposal_Group_6.pdf # Project proposal document
│       └── Model_NS.R            # Initial NS model testing script
└── papers                        # Directory for referenced research papers
    └── A cross-sectional application of the Nelson-Siegel-Svensson model to several negative yield cases.pdf

```


## File Descriptions

### Main Project Files

- `Cubic_Splines.R`: Code for applying cubic splines to Eurozone yield curve data, with an analysis of optimal knot placement.  
- `NS_Model.R`: Code for fitting the Nelson-Siegel (NS) and Nelson-Siegel-Svensson (NSS) models to Eurozone data.  
- `Process_Data.R`: Script to clean and preprocess Eurozone yield curve data.  
- `categorisation.R`: Script for categorizing Eurozone yield curves based on the rate environment (negative, positive, or near-zero).
- `model_compare.R`: Compares the NS, NSS, and cubic spline models on Eurozone data, evaluating metrics such as adjusted R² and MSE.


### Japan Codes Directory

Contains model analysis scripts and visualizations specific to Japan:
- `japan_results/Cubic_Splines.R`: Code for cubic spline modeling of Japanese data.
- `japan_results/NS_Model.R`: Code for NS and NSS models on Japanese data.
- `japan_results/categorisation.R`: Script for categorizing Japan yield curves based on the rate environment (negative, positive, or near-zero).




### Data Directory

The project makes use of various datasets for the Eurozone and Japan analyses. These datasets are critical for fitting yield curves and assessing model performance across different time periods and interest rate regimes.

### Eurozone Data

- `data_aaa.csv`: This dataset contains Eurozone bond data for AAA-rated bonds. It includes historical spot rates and maturities for high-quality bonds, used in modeling yield curves with the NS and NSS models.
- `data_all.csv`: This dataset includes Eurozone bond data for all bond types (not limited to AAA-rated bonds). It is used to analyze the overall yield curve behavior across various bond categories.
- `df_results_exns.csv`: This file contains the results from the NSS model applied to Eurozone data, covering all time periods considered in the analysis.
- `df_results_exns_params.csv`: This file holds the parameters from the NSS model applied to Eurozone data, also covering all time periods. It includes the fitted parameters for each period.
- `df_results_ns.csv`: This dataset contains the results from the NS model applied to Eurozone data across all time periods.
- `df_results_ns_params.csv`: This file contains the parameters from the NS model for Eurozone data across all time periods.
- `df_results_spline.csv`: This file holds the results from the cubic spline model applied to Eurozone data, covering all time periods.
- `df_spot_rates_aaa.csv`: This is the cleaned Eurozone data for AAA-rated bonds, which is the result of preprocessing applied to the raw spot rates data.
- `df_spot_rates_all.csv`: This is the cleaned Eurozone data for all bond types, also resulting from preprocessing.

### Japan Data

- `japan_df_results_exns.csv`: This dataset contains the results of applying the NSS model to Japanese bond data across all time periods.
- `japan_df_results_exns_params.csv`: This file contains the parameters from the NSS model applied to Japanese data across all time periods.
- `japan_df_results_ns.csv`: This dataset contains the results from the NS model applied to Japanese data across all time periods.
- `japan_df_results_ns_params.csv`: This file holds the parameters from the NS model applied to Japanese data across all time periods.
- `japan_df_results_spline.csv`: This file contains the cubic spline model results applied to Japanese data across all time periods.


### App Directory

- `app/simple_app.R`: Shiny app that provides an interactive interface to explore different yield curve shapes. The app allows users to visualize the fit of various models to the Eurozone yield data.



## Installation Guide

To run this project, you’ll need to install several R packages listed in `requirements.txt`. Follow the steps below to set up the required packages:

**Prerequisites**

Ensure you have R installed. You can download R from CRAN.

Installing Packages
1. Save the Package List: The packages required for this project are listed in requirements.txt.  
2. Install Packages: Use the following code snippet in R to install all the packages listed in requirements.txt.

```r
packages <- readLines("requirements.txt")
install.packages(packages)
```


This code will read the list of packages from `requirements.txt` and install any missing packages in your R environment.

**Verifying Installation**

After running the above commands, load the packages in your R session to confirm they are installed correctly:

```r
# Load packages
library(tidyverse)
library(minpack.lm)  
library(forecast)
library(YieldCurve)
library(gridExtra)
```

If there are no errors, you are ready to proceed with the project!



