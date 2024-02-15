# EEG-preprocessing
 # EEG Cleaning Pipeline

This repository contains code for an EEG (Electroencephalography) cleaning pipeline designed to preprocess resting state raw EEG data. The pipeline includes several steps to enhance the quality of EEG signals for further analysis.

## Required toolboxes
Parallel Computing Toolbox (https://it.mathworks.com/products/parallel-computing.html), Brainstorm, \\
c) NoiseTools\\
d) PrepPipeline\\
e) clean_rawdata\\


## Overview

The EEG cleaning pipeline consists of the following steps:

1. **Data Loading**: Raw EEG data is loaded into the pipeline.

2. **Flat segments' Detection**: a portion of the signal is marked as BAD if the difference between consecutive samples is lower than an empirical threshold (1e-08).

3. **Bad channels detection based on Power Spectral Densities**: channels with PSDs in the range of 1 - 45 Hz outside the confidence interval are marked as bad. The confidence interval is equal to the median ± 5 * standard_deviation 

4. **Noisy segments' detection**: portions of the signal are identified and marked as bad by the algorithm according to the following steps: i) band-pass filter between 60 - F_{Nyquist} - 1 Hz (Muscular activity typical range); ii) mean across channels of the channels' envelopes, computed using the absolute value of the Hilbert-transformed signals; iii) definition of a threshold as 10 times the standard deviation of the mean obtained in ii). A robust estimation of the mean's standard deviation is computed as 1.4826 multiplied by the median absolute deviation (MAD). iv) segments displaying values above this threshold are marked as bad.

5. **Notch Filtering via Zapline**: Power line noise is removed with the Zapline algorithm [1]

6. **Robust referencing**: based on the process shown in the Prep pipeline [2]

7. **Data Export**: Cleaned EEG data is exported for further analysis.



## Usage

Follow the instructions in the provided scripts to run each step of the pipeline. Make sure to configure input and output paths as needed.

## References
[1] Alain de Cheveigné, "ZapLine: A simple and effective method to remove power line artifacts", NeuroImage, Volume 207, 15 February 2020, 116356, https://doi.org/10.1016/j.neuroimage.2019.116356
[2] Nima Bigdely-Shamlo, Tim Mullen, Christian Kothe, Kyung-Min Su, and Kay A. Robbins, "The PREP pipeline: standardized preprocessing for large-scale EEG analysis", Frontiers in Neuroinformatics, 2015, 9, pp:16, DOI=10.3389/fninf.2015.00016 
## Contributors

...

Feel free to contribute by opening issues, suggesting improvements, or submitting pull requests!


