function [ data ] = MEEGsim( dip, MEEG, lf, SNR, stim, randdip, waveform, GED, tf, fwhm, plotting, plotdir, rngseed, doICA)
%MEEGsim( dip, MEEG, lf, SNR, stim, randdip, waveform, GED, tf, fwhm, ...
% plotting, plotdir, rngseed) applies condition-specific spectral scanning
% using GED to simulated dipole data.
%
% FUNCTION INPUTS
%   dip         Struct with fields "dip", "frex", "labels",
%                   describing dipole index (dip x 1 vector), 
%                   signal frequencies (dip x 1 vector),
%                   labels for plotting (dip x 1 cell).
%   MEEG        Empty EEGlab-style MEEG struct with sensor locations 
%                   matching leadfield model. Also needs at least fields 
%                   "trials" (int), "pnts" (int), "srate" (int), "times"
%                   (pnts x 1 vector of sample times), "nbchan" (int)
%   lf          Brainstorm-style leadfield model, matching MEEG sensor locs
%   SNR         Struct with fields "SNR.sigamp", "SNR.randamp", "SNR.noiseamp" (ints).
%   stim        1 x 2 vector with signal dipole onset and offset. Make sure
%                   this fits MEEG.times and leaves enough time for baselines.
%   randdip        Struct with fields "num" (number of random dipoles, int),
%                   "prop" (proportion of those dipoles to be active on
%                   each trial), and "max_freq" (max frequency to generate
%                   dipoles; 
%   waveform    'sine' or 'bio'. If 'bio', also pass fwhm.
%   GED         Struct with fields "frex" (1 x n, frequencies at which to
%                   scan).
%   tf          Struct with fields "baseline" (1 x 2, onset and offset), 
%                   "times2save" (1 x n, sample times), and "frex" (1 x n), 
%                   frequencies at which to decompose. Only used for TF
%                   decompositions for plots, can be left empty if plotting is 0.
%   fwhm        Integer value for full width half maximum for dipole
%                    waveforms. Smaller is more sine-like.
%   plotting     1 or 0.
%   plotdir     Directory to save plots to. Can be left empty  if
%                   plotting is 0.
%   rngseed     Seed for the RNG. If not set, is left random.
%   doICA       0, 1, or 2. If 0, do not perform ICAs. If 1, perform
%               broadband, as well as narrowband at all frequencies. If 2,
%               perform broadband, as well as narrowband at 6 and 10 Hz.
%    
%
% FUNCTION OUTPUTS
%   data        Struct containing most of the relevant variables for
%                   plotting, including the following fields:
%   R2_comps   The R^2s between the original dipole signal and the
%                   GED-reconstructed time course.
%   R2_sensors   The R^2s between the original dipole signal and the
%                   best single sensor.
%   SNR         The computed SNR value (signal power / noise power)

%% Setup

dip.loc = lf.GridLoc(dip.dip,:);

% Do a few consistency checks
if mod(MEEG.trials, length(dip.frex)) ~= 0
    error(['Number of trials (' num2str(MEG.trials) ') is not divisible by number of conditions (' num2str(length(dip.frex)) ').']);
end

if max(stim) > MEEG.times(end) || min(stim) < MEEG.times(1)
    error('Stimulus on- or offset is outside time vector.');
end

if diff(stim) < 2*(1/min(GED.frex))
    warning(['Stimulus duration (' num2str(diff(stim)) ' s) covers less than 2 cycles (' num2str(diff(stim)/(1/min(GED.frex))) ' cycles) of minimum scanning frequency (' num2str(min(GED.frex)) ' Hz)']);
end

if diff(tf.baseline) < 2*(1/min(tf.frex))
    warning(['TF baseline (' num2str(diff(tf.baseline)) ' s) covers less than 2 cycles (' num2str(diff(tf.baseline)/(1/min(tf.frex))) ' cycles) of minimum decomposition frequency (' num2str(min(tf.frex)) ' Hz)']);
end


% Construct trial condition membership index
trialtype = ones(MEEG.trials,1);
trialtype(MEEG.trials/2+1:end) = 2;

% Set RNG seed
if nargin >= 14
    rng(rngseed);
end

if strcmpi(waveform, 'sine')
    maxjitter = 300; % max temporal jitter (sine phase offset) in simulation (only if waveform is "sine")
end

% Compute SNR
SNR_value = (SNR.sigamp / (SNR.randamp * (randdip.num * randdip.prop)))^2;


%% Create linear weighting of orientations for each dipole (matching sensor orientation)
dipweight = squeeze(lf.Gain(:,1,:)) .* lf.GridOrient(:,1)' + squeeze(lf.Gain(:,2,:)) .* lf.GridOrient(:,2)' ...
    + squeeze(lf.Gain(:,3,:)) .* lf.GridOrient(:,3)';

% Normalize dipweight so that all dipoles project equally strongly to scalp
dipweight=bsxfun(@rdivide,dipweight,max(dipweight));

%% Generate condition-specific transients for selected dipoles

% Find index for "stimulus" onset and offset
tidx = dsearchn(MEEG.times(:), stim');


if strcmpi(waveform, 'sine')
    waves = zeros(length(dip.frex), length(MEEG.times));
    for f = 1:length(dip.frex)
        % Generate "inner part" of wavelet (allowing for temporal jitter later)
        waves(f,:) = 2*pi*dip.frex(f)*MEEG.times';
    end
elseif strcmpi(waveform, 'bio')
    waves = zeros(length(dip.frex), length(MEEG.times), MEEG.trials); % assumes equal numbers of trials for each condition, which is the case by default
    for f = 1:length(dip.frex)
        % Generate biological-looking waveforms, by multiplying
        % frequency-domain noise with a Gaussian
        waves(f,:,trialtype==f) = filterFGx(randn(sum(trialtype==1), length(MEEG.times)), MEEG.srate, dip.frex(f), fwhm)';
        waves(f,:,trialtype==f) = squeeze(waves(f,:,trialtype==f)) ./ sqrt((squeeze(max(abs(waves(f,:,trialtype==f))))').^2 / 2); % scale RMS amplitude to +- 1
    end
end

% Generate tapered window to have signal fade in and out
taperwin = zeros(1, MEEG.pnts);
taperwin(tidx(1):tidx(2)) = tukeywin(tidx(2)-tidx(1)+1, .6);

%% Generate random shared transients at other random dipoles

% Randomize on- and offset times
min_duration = 0.3; % 300 ms
for i = 1:randdip.num
    diff_r = 0;
    while diff_r <= min_duration
        r = rand(1,2) * (max(MEEG.times) - min(MEEG.times)) + min(MEEG.times);
        diff_r = diff(r);
    end
    ridx(i,:) = dsearchn(MEEG.times(:), [min(r) max(r)]');
end

% Randomize frequencies and create waves
randdip.frex = randi(randdip.max_freq, randdip.num);

if strcmpi(waveform, 'sine')
    rwaves = zeros(randdip.num, length(MEEG.times));
    for f = 1:randdip.num
        % Generate "inner part" of wavelet (allowing for temporal jitter later)
        rwaves(f,:) = 2*pi*randdip.frex(f)*MEEG.times';
    end
elseif strcmpi(waveform, 'bio')
    disp('Generating biological noise...');
    rwaves = zeros(randdip.num, length(MEEG.times), MEEG.trials);
    % Generate biological-looking waveforms, by multiplying
    % frequency-domain noise with a Gaussian (new noise on each trial)
    for f = 1:randdip.num
        rwaves(f,:,:) = filterFGx(randn(MEEG.trials, length(MEEG.times)), MEEG.srate, randdip.frex(f), fwhm)';
        rwaves(f,:,:) = squeeze(rwaves(f,:,:)) ./ sqrt((squeeze(max(abs(rwaves(f,:,:))))').^2 / 2); % scale RMS amplitude to +- 1
    end
end

% Randomize dipole locations
randdip.dip = randi(length(lf.GridLoc), randdip.num, 1);
while sum(ismember(randdip.dip, dip.dip)) > 0 % if random dipole is also signal dipole
    randdip.dip(ismember(randdip.dip, dip.dip)) = randi(length(lf.GridLoc), sum(ismember(randdip.dip, dip.dip)), 1); % re-randomize until all replaced
end

% Generate tapered window to have signals fade in and out
randdip.taperwin = zeros(randdip.num, MEEG.pnts);
for i = 1:randdip.num
    randdip.taperwin(i,ridx(i,1):ridx(i,2)) = tukeywin(ridx(i,2)-ridx(i,1)+1, .6);
end

%% Simulate trials
disp('Simulating trials...');

orig_signal = zeros(MEEG.pnts, MEEG.trials);

for c = 1:length(dip.frex)
    for ti = (MEEG.trials / length(dip.frex)) * (c-1)+1:(MEEG.trials / length(dip.frex)) * c
        dipole_data = SNR.noiseamp * randn(MEEG.pnts, size(lf.Gain,3)); % Generate noise at all dipoles
        
        % Generate transients at shared dipoles
        for f = 1:randdip.num
            if rand < randdip.prop % Randomize how many transients manifest on each trial
                if strcmpi(waveform, 'sine')
                    swave = SNR.randamp * sin(rwaves(f,:) + randi(maxjitter,1)); % adding temporal jitter to randomize phase offset
                elseif strcmpi(waveform, 'bio')
                    swave = SNR.randamp * rwaves(f,:,ti);
                end
                swave = swave .* randdip.taperwin(f,:); % apply tapered window to signal
                dipole_data(:,randdip.dip(f)) = dipole_data(:,randdip.dip(f)) + swave';
            end
        end
        
        % Generate signal at condition-specific dipoles
%         for f = 1:length(dip.frex)
            if strcmpi(waveform, 'sine')
                swave = SNR.sigamp * sin( waves(c,:) + randi(maxjitter,1)); % adding temporal jitter to randomize phase offset
            elseif strcmpi(waveform, 'bio')
                swave = SNR.sigamp * waves(c,:,ti);
            end
            swave = swave .* taperwin; % apply tapered window to signal
            
            % Save swave to variable, for TF decomposition + to see how
            % well GED reconstructs it
            orig_signal(:, ti) = dipole_data(:,dip.dip(c)) + swave';
%             orig_signal(:, ti) = swave'; % without white noise
            
            dipole_data(:,dip.dip(c)) = dipole_data(:,dip.dip(c)) + swave';
            
%         end
        
        % Project to scalp (should be chans X time)
        MEEG.data(:,:,ti) = (dipole_data * dipweight')';
    end
    condition(c,:,:) = dipole_data; % for later plotting
end
disp('Finished simulating trials.');

%% Apply spectral scanning GED; A vs. B and B vs. A

% Initialize structs holding GED and ICA results
GED_A = ([]); % A > B
GED_B = ([]); % B > A
ICA = struct('bb_A', [], 'bb_B', [], 'nb_A', [], 'nb_B', []); % broadband/narrowband, A>B/B>A

%% Construct Morlet wavelets
% GED.frex is taken from inputs

% Gaussian width
min_width = 4; % cycles
max_width = 10; % cycles
s = 2 * ( (logspace(log10(min_width), log10(max_width), length(GED.frex)))./(2*pi*GED.frex) ).^2;

% Create wavelets
wt = -2:1/MEEG.srate:2+1/MEEG.srate;
nWave = length(wt);
halfw = floor(nWave/2);
nConv = MEEG.pnts * MEEG.trials + nWave - 1;

waveX = zeros(length(GED.frex), nConv);

for f = 1:length(GED.frex)
    temp = exp(1i * 2 * pi * GED.frex(f) * wt) .* exp( -wt.^2 / s(f) );
    waveX(f,:) = fft(real(temp), nConv); % real wavelet, not complex
    waveX(f,:) = waveX(f,:) / max(waveX(f,:)); % normalize frequency-domain representation
end

%% Loop over frequencies
for f = 1:length(GED.frex)
    freq = GED.frex(f);
    
    %% Filter entire trial through wavelet convolution
    disp(['Filtering + GED on frequency ' num2str(f) ' of ' num2str(length(GED.frex)) '...']);
    
    % Reshape to channels x time
    data = reshape(MEEG.data, MEEG.nbchan, MEEG.pnts * MEEG.trials);
    
    % Filter data
    temp = 2 * real(ifft( bsxfun(@times, fft(data, nConv, 2), waveX(f,:)), [], 2) ); % multiply with wavelet in the frequency domain; convert to power (2*amp) in the time domain
    fdata = reshape(temp(:, halfw:end-halfw), MEEG.nbchan, MEEG.pnts, MEEG.trials); % convert back to original data dimensions
    
    %% Extract window of interest for each condition
    A = fdata(:, tidx(1):tidx(2), 1:MEEG.trials/2);
    B = fdata(:, tidx(1):tidx(2), MEEG.trials/2+1:MEEG.trials);
    
    %% Peform GEDs
    [GED_A.evecs(f,:,:), GED_A.evals(f,:), GED_A.covA(f,:,:), ~] = performGED(A, B, 1); % A > B
    [GED_B.evecs(f,:,:), GED_B.evals(f,:), GED_B.covB(f,:,:), ~] = performGED(B, A, 1); % B > A

    if doICA == 1 % at each frequency
        % Morlet-filtered is too narrow to successfully perform an ICA;
        % applying a wider filter
        disp(['Applying ' num2str(GED.frex(f)) ' Hz filter with a FWHM of 2 for narrowband ICA...']);

        % Filter data
        fdata = filterFGx(MEEG.data, MEEG.srate, GED.frex(f), 2, 0);

         %% Extract window of interest for each condition
        A = fdata(:, tidx(1):tidx(2), 1:MEEG.trials/2);
        B = fdata(:, tidx(1):tidx(2), MEEG.trials/2+1:MEEG.trials);

        %% Perform ICAs
        EEG = MEEG;
        EEG.data = A;
        % Look for ICs in top 30 PCs; compute R2 only for top source
        [EEG, eucnorm, ~] = pop_runica(EEG, 'icatype', 'jader', 'dataset', 1, 'options', {30});

        % Store ICA results in struct
        ICA.nb_A(f).icawinv = EEG.icawinv;
        ICA.nb_A(f).icaweights = EEG.icaweights;
        ICA.nb_A(f).icaact = EEG.icaact;
        ICA.nb_A(f).icachansind = EEG.icachansind;
        ICA.nb_A(f).eucnorm = eucnorm;

        clear EEG;
        EEG = MEEG;
        EEG.data = B;
        
        [EEG, eucnorm, ~] = pop_runica(EEG, 'icatype', 'jader', 'dataset', 1, 'options', {30});

        % Store ICA results in struct
        ICA.nb_B(f).icawinv = EEG.icawinv;
        ICA.nb_B(f).icaweights = EEG.icaweights;
        ICA.nb_B(f).icaact = EEG.icaact;
        ICA.nb_B(f).icachansind = EEG.icachansind;
        ICA.nb_B(f).eucnorm = eucnorm;
    elseif doICA == 2  % 6 and 10 Hz only
        if f == dsearchn(GED.frex(:), dip.frex(1)) % Only do ICAs at 6 and 10 Hz (providing ICA with prior knowledge about frequencies)
            % Morlet-filtered is too narrow to successfully perform an ICA;
            % applying a wider filter
            disp(['Applying ' num2str(GED.frex(f)) ' Hz filter with a FWHM of 2 for narrowband ICA...']);

            % Filter data
            fdata = filterFGx(MEEG.data, MEEG.srate, GED.frex(f), 2, 0);

            % Extract window of interest for condition A
            A = fdata(:, tidx(1):tidx(2), 1:MEEG.trials/2);
            
            clear EEG;
            EEG = MEEG;
            EEG.data = A;
            % Look for ICs in top 30 PCs; compute R2 only for top source
            [EEG, eucnorm, ~] = pop_runica(EEG, 'icatype', 'jader', 'dataset', 1, 'options', {30});
            
            % Store ICA results in struct
            ICA.nb_A.icawinv = EEG.icawinv;
            ICA.nb_A.icaweights = EEG.icaweights;
            ICA.nb_A.icaact = EEG.icaact;
            ICA.nb_A.icachansind = EEG.icachansind;
            ICA.nb_A.eucnorm = eucnorm;

        elseif f == dsearchn(GED.frex(:), dip.frex(2))
            % Morlet-filtered is too narrow to successfully perform an ICA;
            % applying a wider filter
            disp(['Applying ' num2str(GED.frex(f)) ' Hz filter with a FWHM of 2 for narrowband ICA...']);

            % Filter data
            fdata = filterFGx(MEEG.data, MEEG.srate, GED.frex(f), 2, 0);

            % Extract window of interest for condition B
            B = fdata(:, tidx(1):tidx(2), MEEG.trials/2+1:MEEG.trials);
            
            clear EEG;
            EEG = MEEG;
            EEG.data = B;
            [EEG, eucnorm, ~] = pop_runica(EEG, 'icatype', 'jader', 'dataset', 1, 'options', {30});

            % Store ICA results in struct
            ICA.nb_B.icawinv = EEG.icawinv;
            ICA.nb_B.icaweights = EEG.icaweights;
            ICA.nb_B.icaact = EEG.icaact;
            ICA.nb_B.icachansind = EEG.icachansind;
            ICA.nb_B.eucnorm = eucnorm;
        end
    end
end % end frequency loop

%% Compute R squared
%% Identify the "best component" for each signal dipole
% Find peaks in eigenspectra
if size(GED_A.evals,1) > 2
    [Apeakheight, Apeaks] = findpeaks(GED_A.evals(:,1));
    [Bpeakheight, Bpeaks] = findpeaks(GED_B.evals(:,1));
else
    [Apeakheight, Apeaks] = max(GED_A.evals(:,1));
    [Bpeakheight, Bpeaks] = max(GED_B.evals(:,1));
end

% If no peaks are found (can happen for low SNR), just grab highest value overall
if length(Apeaks) < 1
   [Apeakheight, Apeaks] = max(GED_A.evals(:,1));
end
if length(Bpeaks) < 1
   [Bpeakheight, Bpeaks] = max(GED_B.evals(:,1));
end

Afrex = GED.frex(Apeaks);
Bfrex = GED.frex(Bpeaks);

% Sort peak indices and frequencies by peak height
[~, Asidx] = sort(Apeakheight, 'descend');
[~, Bsidx] = sort(Bpeakheight, 'descend');

Apeaks = Apeaks(Asidx);
Afrex = Afrex(Asidx);
Bpeaks = Bpeaks(Bsidx);
Bfrex = Bfrex(Bsidx);

% First, select spectral scanning peaks as components.
% Grab number of highest peaks matching the number of dipoles per
% condition.
compA = GED_A.evecs(Apeaks(1),:,1); % first sorted peak, first eigenvector
compB = GED_B.evecs(Bpeaks(1),:,1); % first sorted peak, first eigenvector

% Construct component time courses from sensor data
compAts = zeros(size(compA,1), size(MEEG.data,2), size(MEEG.data,3));
compBts = zeros(size(compB,1), size(MEEG.data,2), size(MEEG.data,3));
for i = 1:size(compA,1)
    compAts(i,:,:) = sum(compA(i,:)' .* MEEG.data);   % components x time x trials
end
for i = 1:size(compB,1)
    compBts(i,:,:) = sum(compB(i,:)' .* MEEG.data);   % components x time x trials
end

%% Identify the "best sensor" for each signal dipole
% Look for ICs in top 30 PCs; compute R2 only for top source for comparison sensor-level and GED analysis
% Find the sensor that each signal dipole projects most strongly to. This
% represents the assumption that we have a priori (correct) ideas about
% where the phenomenon under study should manifest (i.e., the best-case
% scenario for sensor-level data analysis).

[~, bestsensor] = max(dipweight(:,dip.dip));

%% R2R: apply ICA to broadband sensor data

if doICA ~= 0
    % Convert simulation results to be compatible with EEGlab
    EEG = MEEG;
    % Broadband, A trials:
    EEG.data = MEEG.data(:, tidx(1):tidx(2), 1:MEEG.trials/2); % ICA on broadband A trials

    % Apply ICA algorithm
    EEG = pop_runica(EEG, 'icatype', 'jader', 'dataset', 1, 'options', {30});
    % Look for ICs in top 30 PCs; compute R2 only for top source

    % Store ICA results in struct
    ICA.bb_A.icawinv = EEG.icawinv;
    ICA.bb_A.icaweights = EEG.icaweights;
    ICA.bb_A.icaact = EEG.icaact;
    ICA.bb_A.icachansind = EEG.icachansind;

    % Broadband, B trials:
    EEG.data = MEEG.data(:, tidx(1):tidx(2), MEEG.trials/2+1:MEEG.trials); % ICA on broadband B trials

    % Apply ICA algorithm
    EEG = pop_runica(EEG, 'icatype', 'jader', 'dataset', 1, 'options', {30});
    % Look for ICs in top 30 PCs; compute R2 only for top source

    % Store ICA results in struct
    ICA.bb_B.icawinv = EEG.icawinv;
    ICA.bb_B.icaweights = EEG.icaweights;
    ICA.bb_B.icaact = EEG.icaact;
    ICA.bb_B.icachansind = EEG.icachansind;
end

%% Compute R2s
% Compute R2s for raw component time series, raw sensor-level data for that
% condition:
for i = 1:2
    temp = corrcoef(squeeze(MEEG.data(bestsensor(i),:,trialtype==i)), orig_signal(:,trialtype==i));
    R2_sensors(i) = temp(1,2)^2;
end

temp = corrcoef(compAts(1,:,trialtype==1), orig_signal(:,trialtype==1));
R2_comps(1) = temp(1,2)^2;

temp = corrcoef(compBts(1,:,trialtype==2), orig_signal(:,trialtype==2));
R2_comps(2) = temp(1,2)^2;

if doICA ~= 0
    % Compute ICA activations and compute R2s for ICA - broadband
    temp = ICA.bb_A.icaweights(1,:) * reshape(MEEG.data, size(MEEG.data,1), []);
    ICA.bb_A.icaact = reshape(temp, 1, size(MEEG.data,2), size(MEEG.data,3));
    temp = corrcoef(ICA.bb_A.icaact(1,:,trialtype==1), orig_signal(:,trialtype==1));
    ICA.bb_A.R2 = temp(1,2)^2;

    temp = ICA.bb_B.icaweights(1,:) * reshape(MEEG.data, size(MEEG.data,1), []);
    ICA.bb_B.icaact = reshape(temp, 1, size(MEEG.data,2), size(MEEG.data,3));
    temp = corrcoef(ICA.bb_B.icaact(1,:,trialtype==2), orig_signal(:,trialtype==2));
    ICA.bb_B.R2 = temp(1,2)^2;

    % Compute ICA activations and compute R2s for ICA - narrowband
    if doICA == 1 % ICA at every frequency
        for f = 1:length(GED.frex)
            temp = ICA.nb_A(f).icaweights(1,:) * reshape(MEEG.data, size(MEEG.data,1), []);
            ICA.nb_A(f).icaact = reshape(temp, 1, size(MEEG.data,2), size(MEEG.data,3));
            temp = corrcoef(ICA.nb_A(f).icaact(1,:,trialtype==1), orig_signal(:,trialtype==1));
            ICA.nb_A(f).R2 = temp(1,2)^2;

            temp = ICA.nb_B(f).icaweights(1,:) * reshape(MEEG.data, size(MEEG.data,1), []);
            ICA.nb_B(f).icaact = reshape(temp, 1, size(MEEG.data,2), size(MEEG.data,3));
            temp = corrcoef(ICA.nb_B(f).icaact(1,:,trialtype==2), orig_signal(:,trialtype==2));
            ICA.nb_B(f).R2 = temp(1,2)^2;
        end
    elseif doICA == 2 % ICA at 6 and 10 Hz
        % For ICA just at 6 and 10 Hz:
        temp = ICA.nb_A.icaweights(1,:) * reshape(MEEG.data, size(MEEG.data,1), []);
        ICA.nb_A.icaact = reshape(temp, 1, size(MEEG.data,2), size(MEEG.data,3));
        temp = corrcoef(ICA.nb_A.icaact(1,:,trialtype==1), orig_signal(:,trialtype==1));
        ICA.nb_A.R2 = temp(1,2)^2;

        temp = ICA.nb_B.icaweights(1,:) * reshape(MEEG.data, size(MEEG.data,1), []);
        ICA.nb_B.icaact = reshape(temp, 1, size(MEEG.data,2), size(MEEG.data,3));
        temp = corrcoef(ICA.nb_B.icaact(1,:,trialtype==2), orig_signal(:,trialtype==2));
        ICA.nb_B.R2 = temp(1,2)^2;
    end
end

%% Saving results starts here - useful to later load and create plots
% Gather relevant variables into a struct
data = struct('dip', dip, 'randdip', randdip, 'GED_A', GED_A, 'GED_B', GED_B, 'fwhm', fwhm, 'waveform', waveform, ...
    'bestsensor', bestsensor, 'dipweight', dipweight, 'MEEG', MEEG, 'GED', GED, 'orig_signal', orig_signal, ...
    'compAevec', compA, 'compBevec', compB, 'compAts', compAts, 'compBts', compBts, 'Xaf', A, 'Xbf', B, 'f', f, ...
    'trialtype', trialtype, 'R2_sensors', R2_sensors, 'R2_comps', R2_comps, 'SNR_value', SNR_value, ...
    'ICA', ICA);

%% Plotting code starts here
if plotting == 1
    
    % Set plotdirectory to save plots
    if ~exist(plotdir, 'dir')
        mkdir(plotdir);
    end
    
    %% Save config file to plots dir, storing simulation parameters
    fileID = fopen([plotdir 'config.txt'],'w');
    fprintf(fileID, 'Simulation parameters\n\n');
    fprintf(fileID, 'dipole labels: %s %s\n', dip.labels{:});
    fprintf(fileID,'waveform: %s\n',waveform);
    fprintf(fileID,'fwhm: %.2f\n',fwhm);
    fprintf(fileID,'stim: %.2f to %.2f\n',stim);
    fprintf(fileID,'dip.frex: %.2f Hz, %.2f Hz\n',dip.frex);
    fprintf(fileID,'GED.frex: %i\n',GED.frex);
    fprintf(fileID,'randdip.num: %i\n',randdip.num);
    fprintf(fileID,'randdip.prop: %.2f\n',randdip.prop);
    fprintf(fileID,'randdip.max_freq: %i\n',randdip.max_freq);
    fprintf(fileID,'SNR.sigamp: %.2f\n',SNR.sigamp);
    fprintf(fileID,'SNR.randamp: %.2f\n',SNR.randamp);
    fprintf(fileID,'SNR.noiseamp: %.2f\n',SNR.noiseamp);
    fprintf(fileID,'SNR: %.2f\n',SNR_value);
    
    fclose(fileID);
    
    %% Set plotting defaults
    setPlotDefaults();
    
    %% View the selected dipoles
    figure;
    hold on;
    plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), '.');
    colors = {'r', 'g', 'b', 'k'};
    for i = 1:length(dip.frex)
        plot3(lf.GridLoc(dip.dip(i),1), lf.GridLoc(dip.dip(i),2), lf.GridLoc(dip.dip(i),3), 'rs' ,'markerfacecolor',colors{i},'markersize',10);
        text(lf.GridLoc(dip.dip(i),1), lf.GridLoc(dip.dip(i),2), lf.GridLoc(dip.dip(i),3), dip.labels{i});
    end
    ax = gca;
    ax.SortMethod = 'childorder';
    rotate3d on; axis equal;
    title('Brain dipole locations');
    
    filename = [plotdir '01_selected_dipoles'];
    saveas(gca, [filename '.png']);
    savefig([filename '.fig'])
    close;
    
    %% View the dipole projections (all three orientations, plus sensor orientation)
    
    setPlotDefaults();
    
    % Project the dipoles to the scalp
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    clim = 1;
    p = 0;
    for c = 1:length(dip.frex)
        for orient = 1:3
            p = p+1;
            subplot(4,4,p);
            topoplot(lf.Gain(:,orient,dip.dip(c)), MEEG.chanlocs,'maplimits',[-1 1]*clim, 'numcontour',0,'electrodes','off','shading','interp');
            colorbar;
            title([dip.labels{c} ', orientation ' num2str(orient) '/3']);
        end
        p = p+1;
        subplot(4,4,p);
        topoplot(dipweight(:,dip.dip(c)), MEEG.chanlocs,'maplimits',[-1 1]*clim, 'numcontour',0,'electrodes','off','shading','interp');
        colorbar;
        title([dip.labels{c} ', sensor orientation']);
    end
    
    saveas(gca, [plotdir '02_dipole_projections.png']);
    close;
    
    %% Sanity check: plot randomized dipole locations
    %     figure;
    %     hold on;
    %     plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), '.');
    %     for i = 1:randdip.num
    %         plot3(lf.GridLoc(randdip.dip(i,:),1), lf.GridLoc(randdip.dip(i,:),2), lf.GridLoc(randdip.dip(i,:),3), 's', 'markerfacecolor', 'y', 'markersize', 10);
    %     end
    %     rotate3d on, axis equal;
    %     title(['Randomized dipole locations (n = ' num2str(randdip.num) ')']);
    %
    %     filename = [plotdir '03_shared_dipoles'];
    %     saveas(gca, [filename '.png']);
    %     savefig([filename '.fig']);
    %     close;
    
    %% Sanity check: inspect dipole activity on a single trial in each condition
    d = 2; % data distance
    s = 0.4; % data scale
    
    colors = {'red', 'blue', 'green', 'yellow'};
    
    for c = 1:length(dip.frex) % create a plot for each condition
        setPlotDefaults();
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        hold on;
        for i = 1:length(dip.dip) % for each dipole
            plot(MEEG.times, condition(c,:,dip.dip(i)) * s - i*d, colors{i}, 'linewidth', 1); % plot dipole activity in color belonging to condition
        end
        for j = 1:min(randdip.num, 50) % for a max of 50 dipoles
            plot(MEEG.times, condition(c,:,randdip.dip(j)) * s - (i+j)*d, 'k', 'linewidth', 1); % plot shared transients in black
        end
        ylim([-(i+50) 0]);
        title(['Condition ' num2str(c) '/' num2str(length(dip.frex)) '; single trial; 50 dipoles shown (' num2str(randdip.num) ' total, avg ' num2str(randdip.num * randdip.prop) ' per trial)']);
        saveas(gca, [plotdir '04_cond' num2str(c) '_singletrial.png']);
        close;
    end
    
    %% Sanity check: Inspect simulated data topography, single trials
    %     % Note: we're looking at a single time point here, and it's entirely
    %     % possible that the signals of interest are in a trough or near 0 rather
    %     % than a peak. This is not meant to investigate the power from each dipole,
    %     % just to check that the simulation went fine.
    %
    %     figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    %     idx = round(mean(tidx)) + 10; % grab a time idx at the center of stimulus time
    %     for c = 1:length(dip.frex)
    %         trial = MEEG.trials/length(dip.frex) * (c-1) + 1; % trial to show, depending on condition
    %         subplot(1,length(dip.frex),c);
    %         topoplot(squeeze(MEEG.data(:,idx,trial)), MEEG.chanlocs);
    %         title(['Condition ' num2str(c) '/' num2str(length(dip.frex)) '; single trial topographic projection']);
    %         colorbar;
    %     end
    %
    %     saveas(gca, [plotdir '05_cond_singletrial_topoplot.png']);
    %     close;
    
    %% Sanity check: inspect sensor-level data over time during a single trial
    
    s = 0.1 * (1 - SNR.noiseamp); % scale
    d = 4; % distance
    
    trial2plot = 11;
    
    setPlotDefaults();
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    hold on;
    for i = 1:50 % plot first 50 sensors
        plot(MEEG.times, MEEG.data(i,:,trial2plot) * s - i*d, 'k', 'linewidth', 1);
    end
    ylim([-51*d 0]);
    title(['Sensor-level data over time (first 50 sensors, trial ' num2str(trial2plot) ')']);
    
    saveas(gca, [plotdir '06_cond_singletrial_sensorlevel.png']);
    close;
    
    %% Sanity check: plot ERP (mean over trials, mean over sensors) per condition
    
    setPlotDefaults();
    
    for i = 1:2
        erp(i,:) = mean(mean(MEEG.data(:,:,trialtype==i),3));
    end
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(2,1,1);
    plot(MEEG.times, erp(1,:));
    title('Condition 1 ERP', 'linewidth', 1);
    xlabel('Time (s)');
    subplot(2,1,2);
    plot(MEEG.times, erp(2,:))
    title('Condition 2 ERP', 'linewidth', 1);
    xlabel('Time (s)');
    
    saveas(gca, [plotdir '07_cond_ERP.png']);
    close;
    
    %% Apply FFT to mean sensor data, to check whether frequency attenuation is present in the leadfield model
    
    meansens(1,:) = mean(MEEG.data(:,:,find(trialtype==1,1)),1); % mean over sensors, single trial
    meansens(2,:) = mean(MEEG.data(:,:,find(trialtype==2,1)),1); % mean over sensors, single trial
    nyquist = MEEG.srate/2;
    frex2plot = linspace(0,nyquist,1000);
    fftres(1,:) = abs(fft(meansens(1,:), length(frex2plot)*2)).^2;
    fftres(2,:) = abs(fft(meansens(2,:), length(frex2plot)*2)).^2;
    
    setPlotDefaults();
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(2,2,1);
    plot(MEEG.times, meansens(1,:), 'linewidth', 1);
    xlabel('Time (s)');
    title('Condition 1 trial 1, mean over sensors');
    subplot(2,2,2);
    plot(MEEG.times, meansens(2,:), 'linewidth', 1);
    xlabel('Time (s)');
    title('Condition 2 trial 1, mean over sensors');
    subplot(2,2,3);
    plot(frex2plot, fftres(1,1:length(frex2plot)));
    xlim([0 50]); xticks(2:2:50);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(['Condition 1 frequency spectrum (generated: ' num2str(dip.frex(1)) ' Hz)']);
    subplot(2,2,4);
    plot(frex2plot, fftres(2,1:length(frex2plot)));
    xlim([0 50]); xticks(2:2:50);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(['Condition 2 frequency spectrum (generated: ' num2str(dip.frex(2)) ' Hz)']);
    
    saveas(gca, [plotdir '08_sensor_level_FFT.png']);
    close;
    
    %% Plot max eigenvalue frequency spectra for A > B and B > A
    
    setPlotDefaults();
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(2,1,1);
    hold on;
    plot(GED.frex, GED_A.evals(:,1));
    plot(Afrex(1), GED_A.evals(Apeaks(1),1), 'rs');
    xlabel('Frequency (Hz)');
    ylabel('Power ratio (\lambda)');
    title(['A > B (simulated transients: ' num2str(dip.frex(1)) ' Hz)']);
    xlim([min(GED.frex) max(GED.frex)]);
    xticks(min(GED.frex):2:max(GED.frex));
    
    subplot(2,1,2);
    hold on;
    plot(GED.frex, GED_B.evals(:,1));
    plot(Bfrex(1), GED_B.evals(Bpeaks(1),1), 'rs');
    xlabel('Frequency (Hz)');
    ylabel('Power ratio (\lambda)');
    title(['B > A (simulated transients: ' num2str(dip.frex(2)) ' Hz)']);
    xlim([min(GED.frex) max(GED.frex)]);
    xticks(min(GED.frex):2:max(GED.frex));
    
    suptitle('Frequency-specific peak eigenvalues per GED contrast');
    
    saveas(gca, [plotdir '09_GED_spectrum_peaks.png']);
    close;
    
    %% Plot topoplots at peaks
    
    setPlotDefaults();
    
    % Generate topoplots by multiplying eigenvectors by covA, covB
    for i = 1:length(Apeaks)
        Atopo(i,:) = squeeze(GED_A.covA(Apeaks(i),:,:)) * GED_A.evecs(Apeaks(i), :, 1)';
    end
    for i = 1:length(Bpeaks)
        Btopo(i,:) = squeeze(GED_B.covB(Bpeaks(i),:,:)) * GED_B.evecs(Bpeaks(i), :, 1)';
    end
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    for i = 1:length(Apeaks)
        subplot(2,max(length(Apeaks), length(Bpeaks)), i);
        topoplot(Atopo(i,:), MEEG.chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
        title(['A, ' num2str(Afrex(i),3) ' Hz, \lambda ' num2str(GED_A.evals(Apeaks(i),1),3)]);
    end
    for i = 1:length(Bpeaks)
        subplot(2,max(length(Apeaks), length(Bpeaks)), i+max(length(Apeaks), length(Bpeaks)));
        topoplot(Btopo(i,:), MEEG.chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
        title(['B, ' num2str(Bfrex(i),3) ' Hz, \lambda ' num2str(GED_B.evals(Bpeaks(i),1),3)]);
    end
    
    suptitle('Topoplots of components at spectrum peaks');
    
    saveas(gca, [plotdir '10_GED_spectrum_peaks_topoplots.png']);
    close;
    
    %% TF-decompose single sensors and GED components
    
    % Decompose sensor-level data
    [sensortf, ~] = tfdecomp(MEEG.data(bestsensor,:,:), MEEG, min(tf.frex), max(tf.frex), length(tf.frex), 'lin', 'means', tf.times2save, tf.baseline, trialtype);
    
    % Decompose component time courses
    [compAtf, ~] = tfdecomp(compAts, MEEG, min(tf.frex), max(tf.frex), length(tf.frex), 'lin', 'means', tf.times2save, tf.baseline, trialtype);
    [compBtf, frex] = tfdecomp(compBts, MEEG, min(tf.frex), max(tf.frex), length(tf.frex), 'lin', 'means', tf.times2save, tf.baseline, trialtype);
    
    %% Plot and compare TF decompositions for single (best) sensors and (best) GED components
    
    % Concatenate plot contents, to prevent repeating plotting code
    things2plot = cell(length(bestsensor) + size(compA,1) + size(compB,1),1);
    
    for i = 1:length(bestsensor)
        clim = max(abs(sensortf(:)));
        clim_diff = max(max(max(abs(sensortf(:,:,:,2) - sensortf(:,:,:,3)))));
        things2plot{i} = {
            {sensortf(i, :, :, 2)} % mean over condition 1 trials
            {sensortf(i, :, :, 3)} % mean over condition 2 trials
            {sensortf(i, :, :, 2) - sensortf(i, :, :, 3)} % difference
            {sensortf(i, :, :, 1)} % mean over all trials
            {['dB-normalized power at best MEEG sensor (' MEEG.chanlocs(bestsensor(i)).labels ') for dipole ' num2str(i) ' (' num2str(dip.frex(i)) ' Hz)']} % plot title
            {['11_sensorlevel_TF_sensor' num2str(i)]} % save as
            {[clim clim clim_diff clim]} % color limit
            };
    end
    
    things2plot{2}{3} = {sensortf(i, :, :, 3) - sensortf(i, :, :, 2)}; % difference (other contrast)
    
    for j = 1:size(compA,1)
        clim = max([abs(compAtf(:)); abs(compBtf(:))]);
        clim_diff = max([ max(max(max(abs(compAtf(:,:,:,2) - compAtf(:,:,:,3))))), ...
            max(max(max(abs(compBtf(:,:,:,2) - compBtf(:,:,:,3))))) ]);
        things2plot{i+j} = {
            {compAtf(j, :, :, 2)} % mean over condition 1 trials
            {compAtf(j, :, :, 3)} % mean over condition 2 trials
            {compAtf(j, :, :, 2) - compAtf(j, :, :, 3)} % difference
            {compAtf(j, :, :, 1)} % mean over all trials
            {'dB-normalized power for A > B component ' num2str(j) ' (peak at ' num2str(Afrex(j),3) 'Hz, \lambda = ' num2str(GED_A.evals(Apeaks(j),1),3) ')'} % title
            {['12_complevel_TF_compA_' num2str(j)]} % save as
            %             {ceil(max([abs(compAtf(:)); abs(compBtf(:))]))} % color limit
            {[clim clim clim_diff clim]} % color limit
            };
    end
    
    for k = 1:size(compB,1)
        clim = floor(max([abs(compAtf(:)); abs(compBtf(:))]));
        clim_diff = floor(max([ max(max(max(abs(compAtf(:,:,:,2) - compAtf(:,:,:,3))))), ...
            max(max(max(abs(compBtf(:,:,:,2) - compBtf(:,:,:,3))))) ]));
        things2plot{i+j+k} = {
            {compBtf(k, :, :, 2)} % mean over condition 1 trials
            {compBtf(k, :, :, 3)} % mean over condition 2 trials
            {compBtf(k, :, :, 3) - compBtf(k, :, :, 2)} % difference
            {compBtf(k, :, :, 1)} % mean over all trials
            {'dB-normalized power for B > A component ' num2str(k) ' (peak at ' num2str(Bfrex(k),3) 'Hz, \lambda = ' num2str(GED_B.evals(Bpeaks(k),1),3) ')'} %title
            {['13_complevel_TF_compB_' num2str(k)]} % save as
            %             {ceil(max([abs(compAtf(:)); abs(compBtf(:))]))} % color limit
            {[clim clim clim_diff clim]} % color limit
            };
    end
    
    subtitles = {'Mean over condition 1 trials', 'Mean over condition 2 trials', 'Difference', 'Mean over all trials'}; % the same for all plots
    
    for i = 1:length(things2plot)
        setPlotDefaults();
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        for j = 1:4
            subplot(2,2,j)
            contourf(tf.times2save, tf.frex, squeeze(cell2mat(things2plot{i}{j})), length(tf.frex), 'linecolor', 'none');
            %             yticks(log10(frex(1:2:end))); yticklabels(sprintfc('%.1f', frex(1:2:end)));
            title(subtitles{j});
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            c = cell2mat(things2plot{i}{7});
            caxis([-1 1] * c(j)); colormap jet; colorbar;
        end
        suptitle(strjoin(things2plot{i}{5}, ''));
        
        saveas(gca, strjoin([plotdir things2plot{i}{6} '.png'], ''));
        close;
    end
    
    %% TF decompose original dipole signal, for comparison
    
    orig_signal = reshape(orig_signal, 1, MEEG.pnts, MEEG.trials); % reshape to sources x time x trials
    
    [diptf, frex] = tfdecomp(orig_signal, MEEG, min(tf.frex), max(tf.frex), length(tf.frex), 'lin', 'means', tf.times2save, tf.baseline, trialtype);
    
    setPlotDefaults();
    
    % Plot TF decomposition
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(2,2,1);
    contourf(tf.times2save, frex, squeeze(diptf(1,:,:,2)), length(frex), 'linecolor', 'none');
    title(['Condition A signal dipole (' num2str(dip.frex(1)) ' Hz)']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    caxis([-1 1]*max(abs(caxis))+3); colorbar;
    subplot(2,2,2);
    contourf(tf.times2save, frex, squeeze(diptf(1,:,:,3)), length(frex), 'linecolor', 'none');
    title(['Condition B signal dipole (' num2str(dip.frex(2)) ' Hz)']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    caxis([-1 1]*max(abs(caxis))+3); colorbar;
    colormap jet;
    
    saveas(gca, [plotdir '15_signal_tf_decomp.png']);
    close;
    
    if doICA == 1
    %% TF decompose and plot ICA signals
    % Decompose ICA activations
    if doICA == 1
        [ICA_bb_A_tf, frex] = tfdecomp(ICA.bb_A.icaact, MEEG, min(tf.frex), max(tf.frex), length(tf.frex), 'lin', 'means', tf.times2save, tf.baseline, trialtype);
        [ICA_bb_B_tf, ~] = tfdecomp(ICA.bb_B.icaact, MEEG, min(tf.frex), max(tf.frex), length(tf.frex), 'lin', 'means', tf.times2save, tf.baseline, trialtype);
        [ICA_nb_A_tf, ~] = tfdecomp(ICA.nb_A.icaact, MEEG, min(tf.frex), max(tf.frex), length(tf.frex), 'lin', 'means', tf.times2save, tf.baseline, trialtype);
        [ICA_nb_B_tf, ~] = tfdecomp(ICA.nb_B.icaact, MEEG, min(tf.frex), max(tf.frex), length(tf.frex), 'lin', 'means', tf.times2save, tf.baseline, trialtype);
    end
    
    setPlotDefaults();
    
    % Plot TF decomposition
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(2,2,1);
    contourf(tf.times2save, frex, squeeze(ICA_bb_A_tf(1,:,:,2)), length(frex), 'linecolor', 'none');
    title('Condition A broadband IC');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    caxis([-1 1]*max(abs(caxis))+3); colorbar;
    
    subplot(2,2,2);
    contourf(tf.times2save, frex, squeeze(ICA_bb_B_tf(1,:,:,3)), length(frex), 'linecolor', 'none');
    title('Condition B broadband IC');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    caxis([-1 1]*max(abs(caxis))+3); colorbar;
    colormap jet;
    
    subplot(2,2,3);
    contourf(tf.times2save, frex, squeeze(ICA_nb_A_tf(1,:,:,2)), length(frex), 'linecolor', 'none');
    title('Condition A narrowband IC (6 Hz)');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    caxis([-1 1]*max(abs(caxis))+3); colorbar;
    
    subplot(2,2,4);
    contourf(tf.times2save, frex, squeeze(ICA_nb_B_tf(1,:,:,3)), length(frex), 'linecolor', 'none');
    title('Condition B narrowband IC (10 Hz)');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    caxis([-1 1]*max(abs(caxis))+3); colorbar;
    colormap jet;
    
    saveas(gca, [plotdir '16_ICs_tf_decomp.png']);
    close;
    
    end
    
end % end plotting code

end

