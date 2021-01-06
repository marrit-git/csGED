function [dirs] = setdirs()
%SETDIRS adds necessary directories to the MATLAB path. This is defined as 
% a function so only one file needs to be updated when running code from a
% new location.

dirs.eeglab = '/home/mzuure/Documents/eeglab/';

dirs.home = '/home/mzuure/Data/csGED_fresh/'; % needs to exist
dirs.data = [dirs.home '/Data/']; % needs to exist
dirs.results = [dirs.home '/Results/']; % does not need to exist, but is distributed with the data

dirs.emp_plots = [dirs.home 'Plots/Empirical data/']; % created automatically
dirs.sim_plots = [dirs.home 'Plots/Simulated data/']; % created automatically
dirs.figs = [dirs.home 'Figures/']; % created automatically

addpath(genpath(dirs.eeglab));

end