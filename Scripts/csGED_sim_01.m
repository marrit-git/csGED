%% Simulate EEG data containing transient oscillations from different dipoles.
% Then, compare ground truth reconstruction techniques (GED, ICA, and best
% sensor).
% Code based on Mike X Cohen's EEG simulation code from
% the Linear Algebra for Neuroscientists course. We use a Brainstorm-
% generated MEG lead field model (overlapping spheres, CTF275 sensors 
% from empirical dataset, cortical surface only, default anatomy, 
% manually aligned).
%
% Script saves simulation results if not already present (to dirs.results), 
% produces plots (in dirs.sim_plots) and generates semi-polished panels for 
% manuscript figures (in dirs.figs).

%% Preamble

close all; clear;

dirs = setdirs();

if(~exist(dirs.sim_plots, 'dir'))
    mkdir(dirs.sim_plots);
end

if(~exist(dirs.results, 'dir'))
    mkdir(dirs.results);
end

% Load EEGlab
eeglab;

sim_exist_msg = sprintf(['\nThese files were distributed with the data to save processing time. ' ...
    'csGED_sim_01.m will load them instead of rerunning the simulation. ' ...
    'You can delete them if you want to verify the simulation results.']);

%% First do all the simulations, so that all data is ready for plotting

% General settings:
% Load MEG forward model and empty MEG struct containing associated sensor locations
% lf is forward model (leadfield model); MEG is struct
load([dirs.data 'emptyMEG_CTF275.mat']);

% Reshape gain matrix to separate dipole locations from dipole orientations
lf.Gain = reshape(lf.Gain, size(lf.Gain,1), 3, size(lf.GridLoc,1));

% Drop dipoles to speed up code: from 15000 to 3000
lf.Gain = lf.Gain(:,:,1:5:end);
lf.GridLoc = lf.GridLoc(1:5:end,:);
lf.GridOrient = lf.GridOrient(1:5:end,:);

% Normalize leadfield so that all dipoles project equally strongly to scalp
lf.Gain=bsxfun(@rdivide,lf.Gain,max(lf.Gain));

%% Dipole selection
% Create linear weighting of orientations for each dipole (matching sensor orientation)
dipweight = squeeze(lf.Gain(:,1,:)) .* lf.GridOrient(:,1)' + squeeze(lf.Gain(:,2,:)) .* lf.GridOrient(:,2)' ...
    + squeeze(lf.Gain(:,3,:)) .* lf.GridOrient(:,3)';
% Normalize dipweight so that all dipoles project equally strongly to scalp
dipweight=bsxfun(@rdivide,dipweight,max(dipweight));

% Compute spatial overlap between each pair of dipoles
R2s = corrcoef(dipweight).^2;

% % Plot histogram of spatial overlap, to get a feel for variability
% idx = logical(triu(ones(size(R2s)),1));
% R2_vector = R2s(idx);
% hist(R2_vector,100);

% Pick two dipoles with a specified amount of overlap
R2tochoose = 0.4;
[~,pos] = min(abs(R2s(:) - R2tochoose));
[dip.dip(1),dip.dip(2)]=ind2sub(size(R2s),pos);

dip.frex = [6 10];
dip.labels = {'dipole 1', 'dipole 2'};

% % Plot dipoles in 3D, to inspect selection
% setPlotDefaults();
% figure;
% hold on;
% plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), '.');
% colors = {'r', 'g', 'b', 'k'};
% for i = 1:length(dip.frex)
%     plot3(lf.GridLoc(dip.dip(i),1), lf.GridLoc(dip.dip(i),2), lf.GridLoc(dip.dip(i),3), 'rs' ,'markerfacecolor',colors{i},'markersize',10);
%     text(lf.GridLoc(dip.dip(i),1), lf.GridLoc(dip.dip(i),2), lf.GridLoc(dip.dip(i),3), dip.labels{i});
% end
% ax = gca;
% ax.SortMethod = 'childorder';
% rotate3d on; axis equal;
% title('Brain dipole locations');
% 
% % Plot dipole projections
% figure;
% clim = 1;
% for c = 1:length(dip.dip)
%     subplot(2,1,c);
%     topoplot(dipweight(:,dip.dip(c)), MEG.chanlocs,'maplimits',[-1 1]*clim, 'numcontour',0,'electrodes','off','shading','interp');
%     colorbar;
%     title(dip.labels{c});
% end
% suptitle(['Spatial R^{2} = ' num2str(R2s(dip.dip(1), dip.dip(2)))]);

%% Simulation 1 (scan at just 6 Hz):
fname = [dirs.results 'sim1.mat'];

% Initialize structure of data to simulate
MEG.trials = 200; % trials to simulate
MEG.srate = 1000; % sampling rate in Hz
MEG.times = 0:1/MEG.srate:3; % time vector
MEG.pnts = length(MEG.times); % time points per trial

% General settings
waveform = 'bio'; % options: 'sine' for pure sine waves, or 'bio' for more biological-looking signals
fwhm = 2; % full width half maximum for "bio" bandpass filter. Ignored if waveform is "sine".

stim = [1 2]; % "stimulus" onset and offset, in seconds

randdip.num = 1000; % the number of random dipole signals to generate (shared across conditions)
randdip.prop = 0.3; % on which proportion of trials each random dipole is active
randdip.max_freq = 25; % maximum frequency in Hz for random dipoles

SNR.sigamp = 3; % amplitude of dipole signals
SNR.randamp = 3; % amplitude of randomized dipole signals
SNR.noiseamp = 0.5;

% % Frequencies at which to perform spectral scanning
min_freq = 6;
max_freq = 6;
num_frex = 1; % note: on account of Morlet wavelet width spacing in MEEGsim, simulation will be run with a wavelet width of 10 cycles at 6 Hz
GED.frex = 6;

% Frequencies at which to perform TF decomposition (for plotting)
tf.frex = linspace(min_freq, max_freq, 1);
tf.baseline = [0.2 0.8];
tf.times2save = [0.5:0.01:2.5]; % set to right around stim

if exist(fname, 'file')
    warning(['Simulated data file ' fname ' already exists.' sim_exist_msg]);
    load(fname);
    disp('Loaded data file.')
else
%     Run simulation (with RNG seed set to 0; deterministic)
    sim1 = MEEGsim(dip, MEG, lf, SNR, stim, randdip, waveform, GED, tf, fwhm, 0, [], 0, 0);
    save(fname, 'sim1', '-v7.3');
end


%% Simulation 2 (full range of frequencies)

% Update frequencies: first only 6 Hz, now whole spectrum
% % Frequencies at which to perform spectral scanning
min_freq = 2;
max_freq = 20;
num_frex = 50;
GED.frex = linspace(min_freq, max_freq, num_frex);

tf.frex = linspace(min_freq, max_freq, 40);

doICA = 1; % broadband + all frequencies

SNR.sigamp = 3; % amplitude of dipole signals
SNR.randamp = 3; % amplitude of randomized dipole signals
SNR.noiseamp = 0.5;

fname = [dirs.results 'sim2.mat'];
if exist(fname, 'file')
    warning(['Simulated data file ' fname ' already exists.' sim_exist_msg]);
    load(fname);
    disp('Loaded data file.')
else
    % Run simulation (with RNG seed set to 1; deterministic)
    sim2 = MEEGsim(dip, MEG, lf, SNR, stim, randdip, waveform, GED, tf, fwhm, 0, [], 1, doICA);
    
    save(fname, 'sim2', '-v7.3');
end

%% Extra analysis added post-hoc: R2 values between original signal and
% sensors, components are very low. Is this because of all the distractors
% in other frequencies? Filter both at 6 and 10 Hz and see whether R2s turn
% out a lot higher.

% Filter original signal at 6 & 10 Hz

fwhm = sim2.fwhm;
srate = MEG.srate;

disp('Filtering time series...');
orig_6Hz = filterFGx(sim2.orig_signal', srate, 6, fwhm, 0);
orig_10Hz = filterFGx(sim2.orig_signal', srate, 10, fwhm, 0);

% Filter sensor-level signal at 6 & 10 Hz
sensorA_6Hz = filterFGx(squeeze(sim2.MEEG.data(sim2.bestsensor(1),:,:))', srate, 6, fwhm, 0);
sensorB_10Hz = filterFGx(squeeze(sim2.MEEG.data(sim2.bestsensor(2),:,:))', srate, 10, fwhm, 0);

% Filter component signal at 6 & 10 Hz
compA_6Hz = filterFGx(squeeze(sim2.compAts)', srate, 6, fwhm, 0);
compB_10Hz = filterFGx(squeeze(sim2.compBts)', srate, 10, fwhm, 0);

% Filter ICA signal at 6 & 10 Hz
ICAA_6Hz = filterFGx(squeeze(sim2.ICA.bb_A.icaact)', srate, 6, fwhm, 0);
ICAB_10Hz = filterFGx(squeeze(sim2.ICA.bb_B.icaact)', srate, 10, fwhm, 0);
disp('Done.');


% Compute new R2 values
temp = corrcoef(sim2.orig_signal', sensorA_6Hz);
R2_sensors(1) = temp(1,2)^2;
temp = corrcoef(orig_10Hz, sensorB_10Hz);
R2_sensors(2) = temp(1,2)^2;

temp = corrcoef(orig_6Hz, compA_6Hz);
R2_comps(1) = temp(1,2)^2;
temp = corrcoef(orig_10Hz, compB_10Hz);
R2_comps(2) = temp(1,2)^2;

temp = corrcoef(orig_6Hz, ICAA_6Hz);
R2_ICA(1) = temp(1,2)^2;
temp = corrcoef(orig_10Hz, ICAB_10Hz);
R2_ICA(2) = temp(1,2)^2;


disp('R2 between original signal and sensors, narrowband-filtered:');
disp(R2_sensors);
disp('R2 between original signal and GED components, narrowband-filtered:');
disp(R2_comps);
disp('R2 between original signal and ICA components, narrowband-filtered:');
disp(R2_ICA);

%% Simulation 3: different SNR values, full range of frequencies; ICA at 6 and 10 Hz; save modularly so run doesn't need to take too long

fname = [dirs.results 'sim3.mat'];

min_freq = 2;
max_freq = 20;
num_frex = 50;
GED.frex = linspace(min_freq, max_freq, num_frex);

tf.frex = linspace(min_freq, max_freq, 40);

iterations = 5; % how often to rerun, to average over runs afterwards
sigamps = logspace(log10(0.9), log10(900), 30);

doICA = 2; % broadband + 6 and 10 Hz

if ~exist(fname) % if file does not yet exist
    m = matfile(fname, 'Writable', true);
    m.results = zeros(length(sigamps), iterations, 10); % Initialize m.results shape.
end

m = matfile(fname);
if any(m.results(length(sigamps),iterations,:)) % if results of the last run are nonzero, this simulation is done and we can proceed. 
    warning(['Simulated data file ' fname ' already exists.' sim_exist_msg]);
    disp('Loaded data file.')
    sim3 = m.results;
    
else %Otherwise, start running simulations. Note: the current code will overwrite any existing results at i,j, change i and j to start indexing at nonzero when resuming.
    m = matfile(fname, 'Writable', true);
    counter = 0;
    for j = 1:iterations
        for i = 1:length(sigamps)
            counter = counter+1;
            disp(['Run ' num2str(counter) ' of ' num2str(length(sigamps) * iterations) '...']);
            SNR.sigamp = sigamps(i);
            sim3 = MEEGsim(dip, MEG, lf, SNR, stim, randdip, 'bio', GED, tf, fwhm, 0, [], j, doICA);
            
            if i == 1 && j == 1 % on the first run, output all settings to file
                m.settings = sim3;
                m.fields = {'sigamp', 'SNR', 'R2_sensors A', 'R2_sensors B', 'R2_comps A', 'R2_comps B', ...
                    'R2_ICA bb A', 'R2_ICA bb B', 'R2_ICA nb A', 'R2_ICA nb B'}; % headers for m.results
            end
            
            % Create variable to add to .mat file: [sigamp, SNR, R2s]
            m.results(i,j,:) = reshape([sigamps(i), sim3.SNR_value, sim3.R2_sensors, sim3.R2_comps, ...
                sim3.ICA.bb_A.R2, sim3.ICA.bb_B.R2, sim3.ICA.nb_A.R2, sim3.ICA.nb_B.R2], 1, 1, 10);
            disp(['Saved results from run i ' num2str(i) ', j ' num2str(j) ' to file.']);
        end
    end
    quit;
end

%% Figure generation starts here

%% Figure 1
setPlotDefaults();
set(0,'defaultAxesFontSize',10); 
close;

% Illustration of methods, using simulated data filtered at 6 Hz
fig = 1;
status = mkdir([dirs.sim_plots 'Figure ' num2str(fig)]);

panel = 'A1'; % Condition A and B frequency spectra, averaged over trials

trials_A = mean(sim2.MEEG.data(:,:,sim2.trialtype==1),1); % mean over sensors; all trials for that condition
trials_B = mean(sim2.MEEG.data(:,:,sim2.trialtype==2),1); % mean over sensors; all trials for that condition

nyquist = sim2.MEEG.srate/2;
frex2plot = linspace(0,nyquist,1000);

for ti = 1:MEG.trials/2
    freqspec_A(:, ti) = abs(fft(trials_A(:,:,ti), length(frex2plot)*2)).^2;
    freqspec_B(:, ti) = abs(fft(trials_B(:,:,ti), length(frex2plot)*2)).^2;
end

avgfreqspec_A = mean(freqspec_A,2);
avgfreqspec_B = mean(freqspec_B,2);

ax = axgrid(5, 4, 1, .5, .5, .5, 2, 1, 'centimeters'); % ref column widths: 8.5, 11.6, 17.6 cm (JNeurosci)

axes(ax(1));
plot(frex2plot, avgfreqspec_A(1:length(frex2plot)));
xlim([0 30]); xticks(0:5:30);
% xlabel('Frequency (Hz)');
ylabel('Power'); % not dB, not baseline-normalized. pT^2.
ylim([0 3*1E5]);

set(gca,'box','off');

axes(ax(2));
plot(frex2plot, avgfreqspec_B(1:length(frex2plot)));
xlim([0 30]); xticks(0:5:30);
xlabel('Frequency (Hz)');
ylabel('Power'); % not dB, not baseline-normalized
ylim([0 3*1E5]);

set(gca,'box','off');

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
setPlotDefaults();
set(0,'defaultAxesFontSize',10); 
close;

panel = 'A2'; % Condition A and B filtered Xaf, Xbf for a single frequency (here 6 Hz)

sensors2plot = [21:26]; % random selection
trial2plot = 21; % random choice

s = 0.05; % scale
d = 1.5; % distance

ax = axgrid(4, 4, .5, .5, .5, .5, 1, 1, 'centimeters'); % ref column widths: 8.5, 11.6, 17.6 cm (JNeurosci)

axes(ax(1));
hold on;
for si = 1:length(sensors2plot)
    plot(sim1.Xaf(sensors2plot(si),:,trial2plot) * s - d*si, 'r');
    plot(sim1.Xbf(sensors2plot(si),:,trial2plot) * s - d*si -(length(sensors2plot)*d), 'b');
end
% use sim1; sim3 Xaf, Xbf are filtered at 20 Hz (last frequency scanned)

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
setPlotDefaults();

thetaidx = dsearchn(GED.frex(:), 6);

panel = 'A3'; % Covariance matrices Af and Bf

clim = 10;

figure;
imagesc(squeeze(sim2.GED_A.covA(thetaidx,:,:)));
axis square;
caxis([-1 1]*clim);
% colormap bluewhitered;

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

panel = 'A4';

figure;
imagesc(squeeze(sim2.GED_B.covB(thetaidx,:,:)));
axis square;
caxis([-1 1]*clim);
% colormap bluewhitered;

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;
%%
setPlotDefaults();

panel = 'B1'; % Generalized eigenvalue equation; Wf and Lambda_f
% Covariance matrices Af and Bf are in panel A3

figure;
imagesc(squeeze(sim2.GED_A.evecs(thetaidx,:,:)));
axis square;
caxis([-1 1]*.3);
% colormap bluewhitered;

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

panel = 'B2';

figure;
imagesc(diag(sim2.GED_A.evals(thetaidx,:)));
axis square;
caxis([0 0.1]); % limited for visualization purposes
% colormap bluewhitered;

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
setPlotDefaults();
set(0,'defaultLineMarkerSize',5);

panel = 'C1'; % Eigenspectrum for a single GED

ax = axgrid(6.5, 8.5, 1, .5, 1, .5, 1, 1, 'centimeters');

axes(ax(1));
scatter(1:sim2.MEEG.nbchan, sim2.GED_A.evals(thetaidx,:), 'filled');
xlabel('Components');
ylabel('Power ratio (\lambda)');
xlim([0 30]);

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
setPlotDefaults();

panel = 'C2'; % Topoplot for max eigenvalue

figure;
topo = squeeze(sim2.GED_A.covA(thetaidx,:,:)) * squeeze(sim2.GED_A.evecs(thetaidx,:,1)');
topoplot(topo, sim2.MEEG.chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
axis([axis * 1.2]); % prevent topoplot nose being cut off

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
setPlotDefaults();

panel = 'C3'; % Component time series for max eigenvalue

ax = axgrid(4, 4, .5, .5, .5, .5, 1, 1, 'centimeters');

axes(ax(1));
plot(sim2.MEEG.times, mean(sim2.compAts,3));
xlim([0.8, 2.2]);
ylim([-1 1]*5);

set(gca,'box','off');

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
setPlotDefaults();

set(0,'defaultLineMarkerSize',10);

panel = 'E'; % spectral scanning ICAs for A > B and B > A

ax = axgrid(6, 6, .5, .5, .5, .5, 2, 1, 'centimeters');

subplot(2,1,1);
hold on;
for i = 1:50
    scatter(sim2.GED.frex(i), sim2.ICA.nb_A(i).eucnorm(1), 'filled');
end
title('Narrowband ICA A > B');
xlim([1 21]);
xticks(0:2:20)
xlabel('Frequency (Hz)');
ylabel('Max. Euclidean norm');
set(gca,'box','off');

subplot(2,1,2);
hold on;
for i = 1:50
    scatter(sim2.GED.frex(i), sim2.ICA.nb_B(i).eucnorm(1), 'filled');
end
xticks(0:2:20)
title('Narrowband ICA B > A');
xlim([1 21]);
xlabel('Frequency (Hz)');
ylabel('Max. Euclidean norm');
set(gca,'box','off');

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.png']);

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%% Figure 2
fig = 2;
status = mkdir([dirs.sim_plots 'Figure ' num2str(fig)]);

setPlotDefaults();

panel = 'A1'; 
% Plot original topoplots
topo(:,1) = sim2.dipweight(:,sim2.dip.dip(1));
axis([axis * 1.2]); % prevent topoplot nose being cut off
topo(:,2) = sim2.dipweight(:,sim2.dip.dip(2));
axis([axis * 1.2]); % prevent topoplot nose being cut off
chanlocs = sim2.MEEG.chanlocs;

% Plot reconstructed topoplots
frexidx = dsearchn(sim2.GED.frex(:), sim2.dip.frex');

topo(:,3) = sim2.compAevec * squeeze(sim2.GED_A.covA(frexidx(1),:,:)) * - 1; % sign flip based on topography
topo(:,4) = sim2.compBevec * squeeze(sim2.GED_B.covB(frexidx(2),:,:)) * - 1; % sign flip based on topography

% Broadband ICAs
topo(:,5) = sim2.ICA.bb_A.icawinv(:,1) * - 1; % sign flip based on topography
topo(:,6) = sim2.ICA.bb_B.icawinv(:,1) * - 1; % sign flip based on topography

% % Narrowband ICAs - commented out after realizing ICA is not valid on
% narrowband-filtered data
% topo(:,7) = sim2.ICA.nb_A.icawinv(:,1);
% topo(:,8) = sim2.ICA.nb_B.icawinv(:,1) * - 1; % sign flip based on topography

% Compute spatial R^2 values to report in text
GED_A_R2 = corr(topo(:,1), topo(:,3))^2;
GED_B_R2 = corr(topo(:,2), topo(:,4))^2;
ICA_A_R2 = corr(topo(:,1), topo(:,5))^2;
ICA_B_R2 = corr(topo(:,1), topo(:,6))^2;

figure;
% Original
subplot(2,4,1);
topoplot(topo(:,1), chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
axis([axis * 1.2]); % prevent nose cutoff
fixaxis = axis; % used to set all topoplots to the same size, for easier overlaying
subplot(2,4,2);
topoplot(topo(:,2), chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
axis(fixaxis);
subplot(2,4,3); 
topoplot([], chanlocs(sim2.bestsensor(1)), 'style', 'blank', 'electrodes', 'on'); % plot best sensor separately, overlay in Illustrator
axis(fixaxis);
subplot(2,4,4); 
topoplot([], chanlocs(sim2.bestsensor(2)), 'style', 'blank', 'electrodes', 'on'); % plot best sensor separately, overlay in Illustrator
axis(fixaxis);

% GED reconstruction
subplot(2,4,5);
topoplot(topo(:,3), chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
axis(fixaxis);
subplot(2,4,6);
topoplot(topo(:,4), chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
axis(fixaxis);

% Broadband ICA
subplot(2,4,7); 
topoplot(topo(:,5), chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
axis(fixaxis);
subplot(2,4,8); 
topoplot(topo(:,6), chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
axis(fixaxis);

% % Narrowband ICA - commented out after realizing ICA is not valid on
% narrowband-filtered data
% subplot(2,4,7);
% topoplot(topo(:,7), chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
% axis(fixaxis);
% subplot(2,4,8);
% topoplot(topo(:,8), chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
% axis(fixaxis);

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
setPlotDefaults();

panel = 'B';
% Plot original time series (ERP); 6 and 10 Hz
erp(1,:) = mean(sim2.orig_signal(:,sim2.trialtype==1),2);
erp(2,:) = mean(sim2.orig_signal(:,sim2.trialtype==2),2);

% Plot reconstructed component time series (ERP) (note R2)
erp(3,:) = squeeze(mean(sim2.compAts(:,:,sim2.trialtype==1),3));
erp(4,:) = squeeze(mean(sim2.compBts(:,:,sim2.trialtype==2),3));

% Plot ERP at best sensors (note R2)
erp(5,:) = squeeze(mean(sim2.MEEG.data(sim2.bestsensor(1),:,sim2.trialtype==1),3));
erp(6,:) = squeeze(mean(sim2.MEEG.data(sim2.bestsensor(2),:,sim2.trialtype==2),3));

% Flip so correlation with ground truth is positive
gt1 = mean(sim2.orig_signal(:,sim2.trialtype==1),2);
gt2 = mean(sim2.orig_signal(:,sim2.trialtype==2),2);
for i = 1:3
    temp = corrcoef(erp(2*i-1,:), gt1);
    flip(2*i-1) = sign(temp(1,2));
    temp = corrcoef(erp(2*i,:), gt2);
    flip(2*i) = sign(temp(1,2));
end
erp = bsxfun(@times, erp, flip');

% Select window to depict
widx = dsearchn(sim2.MEEG.times(:), [0.8 2.2]'); % stim is 1 to 2 s

s = 0.5; % scale
d = 1; % distance

figure;
hold on;
for p = 1:6
   plot(sim2.MEEG.times(widx(1):widx(2)), erp(p,widx(1):widx(2)) * (s / max(abs(erp(p,:)))) - (p * d));
end
title(sprintf('SNR: %.3f, R2 comps: %.3f, %.3f; R2 sensors: %.3f, %.3f', sim2.SNR_value, sim2.R2_comps, sim2.R2_sensors))

set(gca, 'box', 'off');

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
setPlotDefaults();
set(0,'defaultAxesFontSize',10); 
close;

panel = 'C';
% Plot TF decompositions

% Initialize structure of data to simulate
MEG.trials = 200; % trials to simulate
MEG.srate = 1000; % sampling rate in Hz
MEG.times = 0:1/MEG.srate:3; % time vector
MEG.pnts = length(MEG.times); % time points per trial

% Set parameters for TF decomposition:
min_freq = 2;
max_freq = 20;
num_frex = 40;

spacing = 'lin';
method = 'means';

times2save = [0.5:0.01:2.5]; % set to right around stim
baseline = [0.2 0.8];
% DEBUG: also decompose baseline
% times2save = [0.0:0.01:2.5];

trialtype = sim2.trialtype;

orig_signal = reshape(sim2.orig_signal, 1, size(sim2.orig_signal,1), size(sim2.orig_signal,2)); % reshape input for tfdecomp

% TF decompose: ground truth, reconstruction, best sensor
f = dsearchn(sim2.GED.frex(:), dip.frex(:));

[gt, frex] = tfdecomp(orig_signal, MEG, min_freq, max_freq, num_frex, spacing, method, times2save, baseline, trialtype);
rec1 = tfdecomp(sim2.compAts, MEG, min_freq, max_freq, num_frex, spacing, method, times2save, baseline, trialtype);
rec2 = tfdecomp(sim2.compBts, MEG, min_freq, max_freq, num_frex, spacing, method, times2save, baseline, trialtype);
sens1 = tfdecomp(sim2.MEEG.data(sim2.bestsensor(1),:,:), MEG, min_freq, max_freq, num_frex, spacing, method, times2save, baseline, trialtype);
sens2 = tfdecomp(sim2.MEEG.data(sim2.bestsensor(2),:,:), MEG, min_freq, max_freq, num_frex, spacing, method, times2save, baseline, trialtype);
ICA_bb_1 = tfdecomp(sim2.ICA.bb_A.icaact, MEG, min_freq, max_freq, num_frex, spacing, method, times2save, baseline, trialtype);
ICA_bb_2 = tfdecomp(sim2.ICA.bb_B.icaact, MEG, min_freq, max_freq, num_frex, spacing, method, times2save, baseline, trialtype);
ICA_nb_1 = tfdecomp(sim2.ICA.nb_A(f(1)).icaact, MEG, min_freq, max_freq, num_frex, spacing, method, times2save, baseline, trialtype);
ICA_nb_2 = tfdecomp(sim2.ICA.nb_B(f(2)).icaact, MEG, min_freq, max_freq, num_frex, spacing, method, times2save, baseline, trialtype);

ax = axgrid(24, 16, .5, .5, .5, .5, 7, 2, 'centimeters');

% plot ground truth

axes(ax(1));
contourf(times2save, frex, squeeze(gt(1,:,:,2)), num_frex, 'linecolor', 'none'); % signal dipole 1 (condition-specific)
c(1) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(1));
ylabel('Frequency (Hz)');

axes(ax(2));
contourf(times2save, frex, squeeze(gt(1,:,:,3)), num_frex, 'linecolor', 'none'); % signal dipole 2 (condition-specific)
c(2) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(2));

axes(ax(3));
contourf(times2save, frex, squeeze(rec1(1,:,:,2)), num_frex, 'linecolor', 'none'); % reconstructed dipole 1 (condition-specific)
c(3) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(3));

axes(ax(4));
contourf(times2save, frex, squeeze(rec2(1,:,:,3)), num_frex, 'linecolor', 'none'); % reconstructed dipole 2 (condition-specific)
c(4) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(4));

axes(ax(5));
contourf(times2save, frex, squeeze(rec1(1,:,:,2) - rec1(1,:,:,3)), num_frex, 'linecolor', 'none'); % reconstructed dipole 1 (condition difference)
c(5) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(5));

axes(ax(6));
contourf(times2save, frex, squeeze(rec2(1,:,:,3) - rec2(1,:,:,2)), num_frex, 'linecolor', 'none'); % reconstructed dipole 2 (condition difference)
c(6) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(6));

axes(ax(7));
contourf(times2save, frex, squeeze(sens1(1,:,:,2)), num_frex, 'linecolor', 'none'); % best sensor for dipole 1 (condition-specific)
c(7) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-2 1] * c(7));

axes(ax(8));
contourf(times2save, frex, squeeze(sens2(1,:,:,3)), num_frex, 'linecolor', 'none'); % best sensor for dipole 2 (condition-specific)
c(8) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(8));

axes(ax(9));
contourf(times2save, frex, squeeze(sens1(1,:,:,2) - sens1(1,:,:,3)), num_frex, 'linecolor', 'none'); % best sensor for dipole 1 (condition difference)
c(9) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(9));

axes(ax(10));
contourf(times2save, frex, squeeze(sens2(1,:,:,3) - sens2(1,:,:,2)), num_frex, 'linecolor', 'none'); % best sensor for dipole 2 (condition difference)
c(10) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(10));

axes(ax(11));
contourf(times2save, frex, squeeze(ICA_bb_1(1,:,:,2)), num_frex, 'linecolor', 'none'); % broadband ICA dipole 1 (condition-specific)
c(11) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(11));

axes(ax(12));
contourf(times2save, frex, squeeze(ICA_bb_2(1,:,:,3)), num_frex, 'linecolor', 'none'); % broadband ICA dipole 2 (condition-specific)
c(12) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(12));

axes(ax(13));
contourf(times2save, frex, squeeze(ICA_bb_1(1,:,:,2) - ICA_bb_1(1,:,:,3)), num_frex, 'linecolor', 'none'); % broadband ICA dipole 1 condition difference
xlabel('Time (s)');
c(13) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(13));

axes(ax(14));
contourf(times2save, frex, squeeze(ICA_bb_2(1,:,:,3) - ICA_bb_2(1,:,:,2)), num_frex, 'linecolor', 'none'); % broadband ICA dipole 2 condition difference
xlabel('Time (s)');
c(14) = floor(max(abs(caxis)) * 2) / 2; % round clim down to nearest 0.5 dB
colormap jet; caxis([-1 1] * c(14));

suptitle(['dB: ' sprintf('%.1f ', c)]);

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%% Figure 3
setPlotDefaults();
fig = 3;
status = mkdir([dirs.sim_plots 'Figure ' num2str(fig)]);

panel = 'A';

% Make arrays of results, averaged over iterations
SNR_value = squeeze(mean(sim3(:,:,2),2)); % this doesn't actually vary over iterations
R2_sensors = squeeze(mean(sim3(:,:,3:4),2));
R2_comps = squeeze(mean(sim3(:,:,5:6),2));
R2_ICA_bb = squeeze(mean(sim3(:,:,7:8),2));
R2_ICA_nb = squeeze(mean(sim3(:,:,9:10),2));

% Plot scatterplot of SNR vs. reconstruction accuracy
% figure;
% scatter(log10(SNR_value), R2_comps(:,1), 50, 'filled', 'DisplayName', 'A > B component');
% xticklabels(10.^cellfun(@str2num, xticklabels)');
% hold on;
% scatter(log10(SNR_value), R2_comps(:,2), 50, 'filled', 'DisplayName', 'B > A component');
% scatter(log10(SNR_value), R2_sensors(:,1), 50, 'filled', 'DisplayName', 'A > B best sensor');
% scatter(log10(SNR_value), R2_sensors(:,2), 50, 'filled', 'DisplayName', 'B > A best sensor');
% scatter(log10(SNR_value), R2_ICA_bb(:,1), 50, 'filled', 'DisplayName', 'A > B broadband ICA');
% scatter(log10(SNR_value), R2_ICA_bb(:,2), 50, 'filled', 'DisplayName', 'B > A broadband ICA');
% scatter(log10(SNR_value), R2_ICA_nb(:,1), 50, 'filled', 'DisplayName', 'A > B narrowband ICA');
% scatter(log10(SNR_value), R2_ICA_nb(:,2), 50, 'filled', 'DisplayName', 'B > A narrowband ICA');

% Alternatively, plot averages over A > B and B > A
% figure;
% scatter(log10(SNR_value), mean(R2_comps(:,1:2),2), 50, 'filled', 'DisplayName', 'Component average');
% xticklabels(10.^cellfun(@str2num, xticklabels)');
% hold on;
% scatter(log10(SNR_value), mean(R2_sensors(:,1:2),2), 50, 'filled', 'DisplayName', 'Best sensor average');
% scatter(log10(SNR_value), mean(R2_ICA_bb(:,1:2),2), 50, 'filled', 'DisplayName', 'Broadband ICA average');
% scatter(log10(SNR_value), mean(R2_ICA_nb(:,1:2),2), 50, 'filled', 'DisplayName', 'Narrowband ICA average');

% Alternatively, plot line plots
% figure;
% plot(log10(SNR_value)', R2_comps(:,1)', '-o', 'DisplayName', 'A > B component');
% xticklabels(10.^cellfun(@str2num, xticklabels)');
% hold on;
% plot(log10(SNR_value)', R2_comps(:,2)', '-o', 'DisplayName', 'B > A component');
% plot(log10(SNR_value)', R2_sensors(:,1)', '-o', 'DisplayName', 'A > B best sensor');
% plot(log10(SNR_value)', R2_sensors(:,2)', '-o', 'DisplayName', 'B > A best sensor');
% plot(log10(SNR_value)', R2_ICA_bb(:,1)', '-o', 'DisplayName', 'A > B broadband ICA');
% plot(log10(SNR_value)', R2_ICA_bb(:,2)', '-o', 'DisplayName', 'B > A broadband ICA');
% plot(log10(SNR_value)', R2_ICA_nb(:,1)', '-o', 'DisplayName', 'A > B narrowband ICA');
% plot(log10(SNR_value)', R2_ICA_nb(:,2)', '-o', 'DisplayName', 'B > A narrowband ICA');

% Alternatively, plot line plots of averages
figure;
plot(log10(SNR_value)', mean(R2_comps(:,1:2),2)', '-o', 'DisplayName', 'Component average');
xticklabels(10.^cellfun(@str2num, xticklabels)');
hold on;
plot(log10(SNR_value)', mean(R2_sensors(:,1:2),2)', '-o', 'DisplayName', 'Best sensor average');
plot(log10(SNR_value)', mean(R2_ICA_bb(:,1:2),2)', '-o', 'DisplayName', 'Broadband ICA average');
plot(log10(SNR_value)', mean(R2_ICA_nb(:,1:2),2)', '-o', 'DisplayName', 'Narrowband ICA average');

ylim([-0.02, 1.02]);
ylabel('Reconstruction accuracy (R^{2})')
xlabel('Signal-to-noise ratio');
legend('Location', 'northwest');

saveas(gcf, [dirs.sim_plots 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
quit;
