%% Perform comparison scanning GED on empirical data from visual imagery task
% Contrast "perception" and "imagery" time windows at a range of
% frequencies, using generalized eigendecomposition (GED). GED results are
% saved to file (to dirs.results) for further analysis.
% Data supplied by Nadine Dijkstra.
%
% Author: Marrit Zuure

%% Set preliminaries
close all; clear;

dirs = setdirs;

%% Loop around subjects
for subno = 1:30
    sublist{subno} = ['S' num2str(subno,'%02.f')];
end

% Reject subjects
sub2reject = [2 5 10 12 13]; % same as in Dijkstra et al., 2018

for subno = 1:length(sublist)
    % Clear variables from previous subject, if applicable
    clearvars -except dirs sublist subno sub2reject;
    
    if ismember(subno, sub2reject)
        warning([ sublist{subno} 'has been rejected. Skipping...']);
        continue;
    end
    
    if exist([dirs.results sublist{subno} '_GED.mat'], 'file')
        msg = sprintf(['A GED file already exists for %s .\n' ...
            'These intermediate results files were distributed together with the data, ' ...
            'to save processing time. You don''t need to run csGED_emp_01.m ' ...
            '(but you can remove the GED file if you want to verify the GED results).'], sublist{subno});
        warning(msg);
        continue;
    end
    
    disp(['Processing subject ' num2str(subno) ' of ' num2str(length(sublist)) '...']);
    
    %% Initialize array structs holding GED results at the end
    GED_IM = struct('evals', [], 'evecs', [], 'covIM', []);
    GED_P = struct('evals', [], 'evecs', [], 'covP', []);
    
    %% Load subject data
    data = load([dirs.data sublist{subno} '_data.mat']);
    data = data.data;
    disp('Finished loading data.');
    
    %% Convert data back to double precision (was stored as single)
    data.trial = cellfun(@double, data.trial, 'un', 0);
    
    %% R2R: perform example analysis on subject 01, removing evoked activity
    if subno == 1
        disp('Computing ERP...');
        % Compute ERP (per-sensor grand average per trial)
        ERP_sum = zeros(size(data.trial{1}));
        for i = 1:length(data.trial)
            ERP_sum = ERP_sum + data.trial{i};
        end
        ERP = ERP_sum / length(data.trial);
        
        % Subtract ERP from each trial
        disp('Subtracting ERP...');
        for i = 1:length(data.trial)
            data.trial{i} = data.trial{i} - ERP;
       end
    end
    
    %% Create plot output directory
    subplotdir = [dirs.emp_plots sublist{subno} '/'];
    if ~exist(subplotdir, 'dir')
        mkdir(subplotdir);
    end
    
    %% Extract some useful variables
    srate = data.fsample;
    pnts = length(data.trial{1});
    nbchan = length(data.label);
    trials = length(data.trialnumbers);
    
    %% Extract data matrix
    datamtx = cell2mat(data.trial);
    % Data scaling is now done on data import
    
    %% Find FP, IM windows
    % Using 200 ms offset in order to ignore initial ERP
    FPtime = data.trialinfo(:,10);
    FPwindow = [ FPtime + 0.2, FPtime + 0.8 ];
    
    cuetime = data.trialinfo(:,12);
    IMwindow = [ cuetime + 0.7, cuetime + 4.0 ]; % cue lasts 500 ms
    
    SPtime = data.trialinfo(:,11);
    SPwindow = [ SPtime + 0.2, SPtime + 0.8 ];
    
    for i = 1:length(data.trial)
        FPwidx(i,:) = dsearchn(data.time{1}(:), FPwindow(i,:)');
        IMwidx(i,:) = dsearchn(data.time{1}(:), IMwindow(i,:)');
        SPwidx(i,:) = dsearchn(data.time{1}(:), SPwindow(i,:)');
    end
    
    % Force all windows to be the same length (rounding error can shift
    % them one sample)
    FPwidx(diff(FPwidx') < mean(diff(FPwidx')),2) = FPwidx(diff(FPwidx') < mean(diff(FPwidx')),2) + 1;
    FPwidx(diff(FPwidx') > mean(diff(FPwidx')),2) = FPwidx(diff(FPwidx') > mean(diff(FPwidx')),2) - 1;
    IMwidx(diff(IMwidx') < mean(diff(IMwidx')),2) = IMwidx(diff(IMwidx') < mean(diff(IMwidx')),2) + 1;
    IMwidx(diff(IMwidx') > mean(diff(IMwidx')),2) = IMwidx(diff(IMwidx') > mean(diff(IMwidx')),2) - 1;
    SPwidx(diff(SPwidx') < mean(diff(SPwidx')),2) = SPwidx(diff(SPwidx') < mean(diff(SPwidx')),2) + 1;
    SPwidx(diff(SPwidx') > mean(diff(SPwidx')),2) = SPwidx(diff(SPwidx') > mean(diff(SPwidx')),2) - 1;
    
    %% Construct Morlet wavelets for low-resolution scan
    % Frequencies
    min_freq = 2;
    max_freq = 80;
    num_frex = 40;
    frex_lowres = logspace(log10(min_freq), log10(max_freq), num_frex);
    
    % Gaussian width
    min_width = 4; % cycles
    max_width = 10; % cycles
    s = 2 * ( (logspace(log10(min_width), log10(max_width), num_frex))./(2*pi*frex_lowres) ).^2;
    
    % Create wavelets
    wt = -2:1/srate:2+1/srate;
    nWave = length(wt);
    halfw = floor(nWave/2);
    nConv = pnts * trials + nWave - 1;
    
    waveX = zeros(num_frex, nConv);
    
    for f = 1:num_frex
        temp = exp(1i * 2 * pi * frex_lowres(f) * wt) .* exp( -wt.^2 / s(f) );
        waveX(f,:) = fft(real(temp ./ max(temp)), nConv); % real wavelet, not complex
    end
    
    %% Perform low resolution scan
    % Loop over frequencies
    disp('Spectral scanning at low resolution...');
    for f = 1:num_frex
        freq = frex_lowres(f);
        
        %% Filter entire trial through wavelet convolution
        disp(['Filtering frequency ' num2str(f) ' of ' num2str(num_frex) '...']);
        
        % Filter data
        temp = 2 * real(ifft( bsxfun(@times, fft(datamtx, nConv, 2), waveX(f,:)), [], 2) );
        fdata = reshape(temp(:, halfw:end-halfw), nbchan, pnts, trials);
        
        %% Extract P and IM window
        for i = 1:trials
            P(:,:,i) = cat(2, fdata(:, FPwidx(i,1):FPwidx(i,2),i), fdata(:, SPwidx(i,1):SPwidx(i,2),i));
            IM(:,:,i) = fdata(:, IMwidx(i,1):IMwidx(i,2),i);
        end
        
        %% Perform GEDs:
        [GED_P_lowres.evecs(f,:,:), GED_P_lowres.evals(f,:), GED_P_lowres.covP(f,:,:), ~] = performGED(P, IM); % P > IM
        [GED_IM_lowres.evecs(f,:,:), GED_IM_lowres.evals(f,:), GED_IM_lowres.covIM(f,:,:), ~] = performGED(IM, P); % IM > P
        
    end % end frequency loop
    
    %% Identify peaks in spectrum
    
    spectrum_lowres_P = GED_P_lowres.evals(:,1);
    spectrum_lowres_IM = GED_IM_lowres.evals(:,1);
    
    peakthres = 0.005; % Very low threshold; would rather have false positives than false negatives
    [peak_P, peakidx_P] = findpeaks(spectrum_lowres_P, 'MinPeakProminence', peakthres);
    [peak_IM, peakidx_IM] = findpeaks(spectrum_lowres_IM, 'MinPeakProminence', peakthres);
    
    %% Plot to see whether peakfinding went well
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(211);
    plot(frex_lowres, spectrum_lowres_P);
    hold on;
    plot(frex_lowres(peakidx_P), peak_P, 'sr');
    xlabel('Frequency (Hz)')
    ylabel('\lambda');
    title([sublist{subno} ' P > IM, low resolution scan']);
    subplot(212);
    plot(frex_lowres, spectrum_lowres_IM);
    hold on;
    plot(frex_lowres(peakidx_IM), peak_IM, 'sr');
    xlabel('Frequency (Hz)')
    ylabel('\lambda');
    title([sublist{subno} ' IM > P, low resolution scan']);
    
    saveas(gca, [subplotdir '00_eigenspectra_lowres.png']);
    close;
    
    %% Construct Morlet wavelets for high-resolution scans
    % Frequencies 
    centerfrex_P = frex_lowres(peakidx_P);
    centerfrex_IM = frex_lowres(peakidx_IM);
    width = 4; % How far around each peak frequency to scan, in Hz
    step = 0.5; % The actual frequency resolution, in Hz
    
    frex_P = [];
    for i = 1:length(centerfrex_P)
        if centerfrex_P(i)-width >= min_freq
            frex_P = [frex_P, centerfrex_P(i)-width:step:centerfrex_P(i)+width];
        end
    end

    frex_IM = [];
    for i = 1:length(centerfrex_IM)
        if centerfrex_IM(i)-width >= min_freq
            frex_IM = [frex_IM, centerfrex_IM(i)-width:step:centerfrex_IM(i)+width];
        end
    end
   
    % Combine frex_IM and frex_P (data from both conditions need to be 
    % filtered at each frequency for the condition comparison; might as 
    % well perform the GED both ways)
    % Merge with frex_lowres; these need to be rerun with the new Gaussian
    % width, otherwise the eigenspectra don't match up (= look choppy)
    frex = unique([frex_P, frex_IM, frex_lowres]);
    num_frex = length(frex);
    
    % Gaussian width
    min_width = 4; % cycles
    max_width = 10; % cycles
    s = 2 * ( (logspace(log10(min_width), log10(max_width), num_frex))./(2*pi*frex) ).^2;
    
    % Create wavelets
    wt = -2:1/srate:2+1/srate;
    nWave = length(wt);
    halfw = floor(nWave/2);
    nConv = pnts * trials + nWave - 1;
    
    waveX = zeros(num_frex, nConv);
    for f = 1:num_frex
        temp = exp(1i * 2 * pi * frex(f) * wt) .* exp( -wt.^2 / s(f) );
        waveX(f,:) = fft(real(temp ./ max(temp)), nConv); % real wavelet, not complex
    end
    
    %% Perform high resolution scan
    % Loop over frequencies
    disp('Spectral scanning at high resolution...');
    for f = 1:num_frex 
        freq = frex(f);
        
        %% Filter entire trial through wavelet convolution
        disp(['Filtering frequency ' num2str(f) ' of ' num2str(num_frex) '...']);
        
        % Filter data
        temp = 2 * real(ifft( bsxfun(@times, fft(datamtx, nConv, 2), waveX(f,:)), [], 2) );
        fdata = reshape(temp(:, halfw:end-halfw), nbchan, pnts, trials);
        
        %% Extract P and IM window
        for i = 1:trials
            P(:,:,i) = cat(2, fdata(:, FPwidx(i,1):FPwidx(i,2),i), fdata(:, SPwidx(i,1):SPwidx(i,2),i));
            IM(:,:,i) = fdata(:, IMwidx(i,1):IMwidx(i,2),i);
        end
        
        %% Peform GEDs:
        [GED_P.evecs(f,:,:), GED_P.evals(f,:), GED_P.covP(f,:,:), ~] = performGED(P, IM); % P > IM
        [GED_IM.evecs(f,:,:), GED_IM.evals(f,:), GED_IM.covIM(f,:,:), ~] = performGED(IM, P); % IM > P
        
    end % end frequency loop
    
    %% Save frequencies to variable
    GED_IM.frex = frex;
    GED_P.frex = frex;
    
    %% Sanity check: plot spectra
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(211);
    plot(GED_P.frex, GED_P.evals(:,1));
    xlabel('Frequency (Hz)')
    ylabel('\lambda');
    title([sublist{subno} ' P > IM']);
    subplot(212);
    plot(GED_IM.frex, GED_IM.evals(:,1));
    xlabel('Frequency (Hz)')
    ylabel('\lambda');
    title([sublist{subno} ' IM > P']); 
    
    saveas(gca, [subplotdir '00_eigenspectra_highres.png']);
    close;
    
    %% Save results to file    
    save([dirs.results sublist{subno} '_GED.mat'], 'GED_P', 'GED_IM');
    
end % end subject loop