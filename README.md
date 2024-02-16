# EEG Cleaning Pipeline

This repository contains code for an EEG (Electroencephalography) cleaning pipeline designed to preprocess resting state raw EEG data. The pipeline includes several steps to enhance the quality of EEG signals for further analysis.

## Required Toolboxes

- Parallel Computing
- EEGLAB
- Brainstorm
- NoiseTools
- PrepPipeline
- clean_rawdata

## Overview

This EEG cleaning pipeline is based on [Pedroni et al., 2019](https://doi.org/10.1016/j.neuroimage.2019.06.046). It consists of the following steps:

1. **Data Loading**: Raw EEG data is loaded into the pipeline.

2. **Flat Segments' Detection**: Portions of the signal are marked as BAD if the difference between consecutive samples is lower than an empirical threshold (1e-08).

3. **Bad Channels Detection Based on Power Spectral Densities**: Channels with PSDs in the 1 - 45 Hz range outside the confidence interval are marked as bad.

4. **Noisy Segments' Detection**: Portions of the signal are identified and marked as bad by the algorithm according to the following steps: i) band-pass filter between 60 - (F_{Nyquist} - 5 Hz) (muscular activity’s typical range); ii) computation of means across channels of the channels' envelopes; iii) definition of a threshold as the mean ± 5 * standard_deviation.

5. **Notch Filtering via Zapline**: Power line noise is removed with the Zapline algorithm [Alain de Cheveigné, 2020](https://doi.org/10.1016/j.neuroimage.2019.116356).

6. **Robust Referencing**: Based on the process shown in the Prep pipeline [Bigdely-Shamlo et al., 2015](https://doi.org/10.3389/fninf.2015.00016).

7. **EOG Regression**: Two methods are implemented: i) singular value decomposition (SVD) and ICLabel classifier; ii) linear regression [Pion-Tonachini et al., 2019](https://doi.org/10.1016/j.neuroimage.2019.05.026).

8. **High Pass Filtering**: Butterworth, 5th order, cutoff frequency = 1 Hz, zero-phase.

9. **Clean Raw Data**: Channel rejection using clean_rawdata EEGLAB Plugin with specific parameters [Kothe et al., 2019].

10. **Artifacts Removal via ICA**: ICs are visually identified by the user, and automatic classification is provided to help in the process.

11. **Bad Channels Interpolation via Cubic Spline**.

12. **Low Pass Filtering**: Butterworth, 5th order, cutoff frequency = 80 Hz, zero-phase.

13. **Data Export**: Cleaned EEG data is exported for further analysis. If specified, an additional table EEGsteps is saved with EEG data after each implemented cleaning step. Moreover, plots comparing raw and cleaned EEG are obtained.

## Usage

Follow the instructions in the provided scripts to run each pipeline step. Configure input and output paths as needed.

## References

1. Pedroni, A., Bahreini, A., & Langer, N. (2019). Automagic: Standardized preprocessing of big EEG data. Neuroimage. doi: 10.1016/j.neuroimage.2019.06.046
2. Alain de Cheveigné, "ZapLine: A simple and effective method to remove power line artifacts", NeuroImage, Volume 207, 15 February 2020, 116356, https://doi.org/10.1016/j.neuroimage.2019.116356
3. Nima Bigdely-Shamlo, Tim Mullen, Christian Kothe, Kyung-Min Su, and Kay A. Robbins, "The PREP pipeline: standardized preprocessing for large-scale EEG analysis", Frontiers in Neuroinformatics, 2015, 9, pp:16, DOI=10.3389/fninf.2015.00016 
4. Luca Pion-Tonachini, Ken Kreutz-Delgado, Scott Makeig, "ICLabel: An automated electroencephalographic independent component classifier, dataset, and website", NeuroImage, Volume 198, September 2019, Pages 181-197, https://doi.org/10.1016/j.neuroimage.2019.05.026
5. clean_rawdata (Version 2.7) [Computer software]

## Contributors

- Jasmin Del Vecchio
- Andrea Canessa
- Chiara Palmisano
- Franziska Pellegrini
- Stefan Haufe

Feel free to contribute by opening issues, suggesting improvements, or submitting pull requests!
