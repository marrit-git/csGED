%% From each subject's GED results, extract peaks in the eigenvalue spectrum
% Additionaly, compute topoplots, component time series, and TF
% decompositions for components found at the peaks. Output is saved to file
% (to dirs.results) for later group-level analysis and plotting.
%
% Requires EEGlab.
%
% Author: Marrit Zuure

%% Preamble

close all; clear;

dirs = setdirs();

for subno = 1:30
    sublist{subno} = ['S' num2str(subno,'%02.f')];
end

% Reject subjects - here rejecting directly from sublist, to prevent gaps
% in multi-subject plots
sub2reject = [2 5 10 12 13]; % same as in Dijkstra et al., 2018
subidx = ones(size(sublist));
subidx(sub2reject) = 0;
sublist = sublist(logical(subidx));

%% Loop around subjects
for subno = 1:length(sublist)
    disp(['Processing subject ' sublist{subno} '...']);
    
    % Clear previous data
    clearvars -except dirs sublist subno data;
    
    %% Load data
    temp = load([dirs.results sublist{subno} '_GED.mat']); % GED results
    
    data(subno).GED_P = temp.GED_P;
    data(subno).GED_IM = temp.GED_IM;
    
    %% Extract subject-specific MEG sensor coordinates into chanlocs using EEGlab
    trialdata = load([dirs.data sublist{subno} '_data.mat']);
    trialdata = trialdata.data;
    
    chanlocs = struct([]);
    % Create EEGlab-compatible chanlocs (struct array)
    for i = 1:270
        chanlocs(i).labels = trialdata.grad.label{i}; % grad = gradiometer, type of MEG
        chanlocs(i).X = trialdata.grad.chanpos(i,1);
        chanlocs(i).Y = trialdata.grad.chanpos(i,2);
        chanlocs(i).Z = trialdata.grad.chanpos(i,3);
    end
    
    % Use EEGlab to repopulate all non-XYZ coordinates
    data(subno).chanlocs = pop_chanedit(chanlocs, 'convert', 'cart2all');
    
    % Sanity check: plot electrode locations
    % topoplot([], chanlocs, 'emarker', {'.', 'k', 4, 1});
    
    %% Extract max eigenvalues at each frequency, for P > IM and IM > P
    
    % Find subject-specific peaks in eigenspectra
    [data(subno).peaks_P.height, data(subno).peaks_P.idx, ...
        data(subno).peaks_P.width, data(subno).peaks_P.prominence] = findpeaks(data(subno).GED_P.evals(:,1));
    [data(subno).peaks_IM.height, data(subno).peaks_IM.idx, ...
        data(subno).peaks_IM.width, data(subno).peaks_IM.prominence] = findpeaks(data(subno).GED_IM.evals(:,1));
    
    % Get frequencies at which peaks occur
    data(subno).peaks_P.frex = data(subno).GED_P.frex(data(subno).peaks_P.idx);
    data(subno).peaks_IM.frex = data(subno).GED_IM.frex(data(subno).peaks_IM.idx);
    
    %% Create channels x time x trials matrix, for computing component time series
    disp('Creating MEG data matrix...');
    MEG = reshape(cell2mat(trialdata.trial), length(trialdata.label), [], length(trialdata.trial));
    disp('Finished creating MEG matrix.');
    
    % Generate topoplots by multiplying eigenvectors by covP, covIM
    % Generate component time series by weighting sensors by first eigenvector
    for i = 1:length(data(subno).peaks_P.idx)
        data(subno).peaks_P.topo(i,:) = squeeze(data(subno).GED_P.covP(data(subno).peaks_P.idx(i),:,:)) * data(subno).GED_P.evecs(data(subno).peaks_P.idx(i), :, 1)';
        data(subno).peaks_P.compts(i,:,:) = mean(squeeze(data(subno).GED_P.evecs(data(subno).peaks_P.idx(i),:,1))' .* MEG,1);
    end

    for i = 1:length(data(subno).peaks_IM.idx)
        data(subno).peaks_IM.topo(i,:) = squeeze(data(subno).GED_IM.covIM(data(subno).peaks_IM.idx(i),:,:)) * data(subno).GED_IM.evecs(data(subno).peaks_IM.idx(i), :, 1)';
        data(subno).peaks_IM.compts(i,:,:) = mean(squeeze(data(subno).GED_IM.evecs(data(subno).peaks_IM.idx(i),:,1))' .* MEG,1);
    end
        
      %% Time-frequency decompose component time series
      
      metadata = struct('pnts', length(trialdata.time{1}), 'times', trialdata.time{1}, ...
          'trials', length(trialdata.trial), 'srate', trialdata.fsample);
      
      times2save = [-0.5:0.2:9];
      min_freq = 4;
      max_freq = 80;
      num_frex = 40;
      tf_frex = logspace(log10(min_freq), log10(max_freq), num_frex);
      baseline = [-0.8 -0.2]; % 0.6 s baseline - can't go down to -1 on account of edge artifacts
      
      disp('TF decomposing component time series...');
      data(subno).peaks_P.tfd.tf = tfdecomp(data(subno).peaks_P.compts, metadata, min_freq, max_freq, num_frex, 'log', 'means', times2save, baseline);
      data(subno).peaks_IM.tfd.tf = tfdecomp(data(subno).peaks_IM.compts, metadata, min_freq, max_freq, num_frex, 'log', 'means', times2save, baseline);
      disp('Finished TF decomposing.');
      
      data(subno).peaks_P.tfd.times2save = times2save;
      data(subno).peaks_P.tfd.frex = tf_frex;
      data(subno).peaks_P.tfd.baseline = baseline;
      data(subno).peaks_P.tfd.metadata = metadata;
      data(subno).peaks_IM.tfd.times2save = times2save;
      data(subno).peaks_IM.tfd.frex = tf_frex;
      data(subno).peaks_IM.tfd.baseline = baseline;
      data(subno).peaks_IM.tfd.metadata = metadata;
      
      %% Interpolate for subject average eigenspectra
      % Frex differ per subject and for averaging they need to be all the same. 
      data(subno).interp_P.frex = min(data(subno).GED_P.frex):0.2:max(data(subno).GED_P.frex);
      data(subno).interp_P.evals = interp1(data(subno).GED_P.frex, data(subno).GED_P.evals(:,1), data(subno).interp_P.frex);
      data(subno).interp_IM.frex = min(data(subno).GED_IM.frex):0.2:max(data(subno).GED_IM.frex);
      data(subno).interp_IM.evals = interp1(data(subno).GED_IM.frex, data(subno).GED_IM.evals(:,1), data(subno).interp_IM.frex);
end

%% Save data to file
disp('Saving data to file...');
save([dirs.data 'csGED_spectrum_peaks.mat'], 'data', '-v7.3');
disp('Done.');
quit;