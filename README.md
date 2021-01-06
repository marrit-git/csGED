# Introduction
This MATLAB code accompanies the following publication:

Zuure, M.B., & Cohen, M. X (2021). Narrowband multivariate source separation for semi-blind discovery of experiment contrasts. _Journal of Neuroscience Methods, 350_, 109063. DOI: https://doi.org/10.1016/j.jneumeth.2020.109063

This code is suited only to analyze the dataset made available at the [Donders Institute data repository](https://doi.org/10.34973/pjq2-et81).

## Note on provided data
The data available [here](https://doi.org/10.34973/pjq2-et81) consists of two folders.
- The `Data` folder, containing:
	-  the MEG data for each subject
	-  a [Brainstorm](https://neuroimage.usc.edu/brainstorm)-generated forward model, for use in the simulations
- The `Results` folder, containing:
	- csGED results (eigenvectors and eigenvalues at each frequency)
	-  simulation results

The materials provided in the `Results` folder are the results of intermediate processing steps, provided to save time (as some of the analysis steps are computationally intensive). The code checks whether these intermediate results are present and loads them. If they are not present, they will be generated from the materials in `Data`.<sup>[1](#fn1)</sup> <sup>[2](#fn2)</sup>

## Dependencies
- MATLAB; R2017b was used but earlier versions may work.
- EEGlab toolbox, available [here](https://sccn.ucsd.edu/eeglab/index.php). We used v14.1.2.
 - The data set to be analyzed, available [here](https://doi.org/10.34973/pjq2-et81).

## How to use
 1. Clone this repository (or download the files) to a location of your liking. 
 2. Download the [data set](https://doi.org/10.34973/pjq2-et81) and put it in the same location, so that on the same level as the `Scripts` folder there is a  `Data` folder containing files files `Sxx_data.mat`, etc., and a `Results` folder containing `Sxx_GED.mat`, etc.
 3. Change the hardcoded paths in`setdirs.m` to point to the different directories and dependencies listed above.
 4. In MATLAB, navigate to the `Scripts` folder and execute the csGED scripts in order. Scripts with `sim` generate and analyze the simulated data, scripts with `emp` analyze the empirical data. Whether you run the `sim` or `emp`  scripts first doesn't matter.
 5. After everything has been run, figure panels from the manuscript have been output to a `Figures` folder, and many additional plots have been output to the `Plots` folder.

## Summary of scripts
- `csGED_sim_01` simulates MEG data from a number of dipoles, as described in the manuscript, and compares different methods of reconstructing the ground truth. There are three different simulations: one with fixed SNR and csGED scanning a single frequency (6 Hz), one with fixed SNR and csGED scanning the full range of frequencies, and one with variable SNR and csGED scanning the full range of frequencies.
- `csGED_emp_01` performs csGED on the empirical data for each subject, with two contrasts: Perception > Imagery and Imagery > Perception. This script can be skipped when using the intermediate csGED results files.
- `csGED_emp_02` constructs the eigenvalue spectrum per subject and identifies spectral peaks.
- `csGED_emp_03`pools data across subjects and generates manuscript figure panels.

## License
This code is released under the MIT license.

## Contact
For further questions regarding this code, please refer to the manuscript or contact me (Marrit Zuure) at my github-associated e-mail address.

## Footnotes
<a name="fn1">1</a>: The recomputed GEDs may be subtly different from the ones provided. This is because the data was converted from double to single precision prior to distribution.

<a name="fn2">2</a>: The recomputed simulations may also differ from the ones provided. The default random generator was changed in MATLAB 2019a, which may be the cause. Alternatively, if you've changed any of the simulation code, the code path may have diverged after setting the random seed, yielding different results. The overall conclusions should be robust to these changes.
