# EEG-preprocessing
 # EEG Cleaning Pipeline

This repository contains code for an EEG (Electroencephalography) cleaning pipeline designed to preprocess resting state raw EEG data. The pipeline includes several steps to enhance the quality of EEG signals for further analysis.

## Required toolboxes
Parallel Computing, EEGLAB, Brainstorm, NoiseTools, PrepPipeline, clean_rawdata.


## Overview

The EEG cleaning pipeline consists of the following steps:

1. **Data Loading**: Raw EEG data is loaded into the pipeline.

2. **Flat segments' Detection**: a portion of the signal is marked as BAD if the difference between consecutive samples is lower than an empirical threshold (1e-08).

3. **Bad channels detection based on Power Spectral Densities**: channels with PSDs in the 1 - 45 Hz range outside the confidence interval are marked as bad. The confidence interval is equal to the median ± 5 * standard_deviation 

4. **Noisy segments' detection**: portions of the signal are identified and marked as bad by the algorithm according to the following steps: i) band-pass filter between 60 - F_{Nyquist} - 1 Hz (Muscular activity typical range); ii) means across channels of the channels' envelopes, computed using the absolute value of the Hilbert-transformed signals; iii) definition of a threshold as the mean ± 10 * standard_deviation. A robust estimation of the mean's standard deviation is computed as 1.4826 multiplied by the median absolute deviation (MAD). iv) segments displaying values above this threshold are marked as bad.

5. **Notch filtering via Zapline**: Power line noise is removed with the Zapline algorithm [1]

6. **Robust referencing**: based on the process shown in the Prep pipeline [2]

7. ** EOG regression**: two methods are implanted: i) the frontal channel with the best blink fitting is selected and subjected to singular value decomposition. EOG artifacts are removed using singular value decomposition (SVD) and ICLabel classifier [3] ii) linear regression

8. **High pass filtering**: Butterworth, 5th order, cutoff frequency = 80 Hz, zero-phase

9. **Clean raw data**: channel rejection using clean_rawata EEGLAB Plugin with the following parameters: 'ChannelCriterion', 'off', 'LineNoiseCriterion', 'off', 'BurstCriterion', 25, 'WindowCriterion', 0.25, 'Highpass', 'off' [4]

10.  **Artifacts removal via ICA**: ICs are visually identified by the user in terms of time series, topoplot, and power spectrum (similarly to step 7). An automatic classification is provided to help the user in the process based on the ICLabel classifier[3].

11. **Bad channels interpolation via cubic spline**

12. **Low pass filtering**: Butterworth, 5th order, cutoff frequency = 80 Hz, zero-phase

13. **Data Export**: Cleaned EEG data is exported for further analysis. If specified in the DefaultParameters.m file, an additional table EEGsteps is saved with EEG data after each implemented cleaning step. Moreover, plots comparing raw and cleaned EEG are obtained.



## Usage

Follow the instructions in the provided scripts to run each pipeline step. Make sure to configure input and output paths as needed.

## References
[1] Alain de Cheveigné, "ZapLine: A simple and effective method to remove power line artifacts", NeuroImage, Volume 207, 15 February 2020, 116356, https://doi.org/10.1016/j.neuroimage.2019.116356
[2] Nima Bigdely-Shamlo, Tim Mullen, Christian Kothe, Kyung-Min Su, and Kay A. Robbins, "The PREP pipeline: standardized preprocessing for large-scale EEG analysis", Frontiers in Neuroinformatics, 2015, 9, pp:16, DOI=10.3389/fninf.2015.00016 
[3] Luca Pion-Tonachini, Ken Kreutz-Delgado, Scott Makeig, "ICLabel: An automated electroencephalographic independent component classifier, dataset, and website", NeuroImage
Volume 198, September 2019, Pages 181-197, https://doi.org/10.1016/j.neuroimage.2019.05.026
[4] Kothe, Miyakoshi, Delorme, C. M. A. (2019). clean_rawdata (Version 2.7) [Computer software]


## Contributors

Jasmin Del Vecchio Del Vecchio, Andrea Canessa, Chiara Palmisano, Franziska Pellegrini, Stefan Haufe

Feel free to contribute by opening issues, suggesting improvements, or submitting pull requests!


