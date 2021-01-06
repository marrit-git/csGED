%% Pool data across subjects and generate plots and figure panels
% Script produces plots (in dirs.emp_plots) and semi-polished panels for 
% manuscript figures (in dirs.figs).
%
% Author: Marrit Zuure

%% Preamble
close all; clear;

dirs = setdirs();

if(~exist(dirs.emp_plots, 'dir'))
    mkdir(dirs.emp_plots);
end

for subno = 1:30
    sublist{subno} = ['S' num2str(subno,'%02.f')];
end

% Reject subjects - here rejecting directly from sublist, to prevent gaps
% in multi-subject plots
sub2reject = [2 5 10 12 13]; % same as in Dijkstra et al., 2018
subidx = ones(size(sublist));
subidx(sub2reject) = 0;
sublist = sublist(logical(subidx));

%% Load data
disp('Loading data...');
data = load([dirs.data 'csGED_spectrum_peaks.mat']);
data = data.data;
disp('Done loading data.');

%% Loop around subjects
for subno = 1:length(sublist)
    disp(['Processing subject ' sublist{subno} '...']);
    
    % Plot individual subject data:
    %% Specify directory to save plots
    subplotdir = [dirs.emp_plots sublist{subno} '/'];
    if ~exist(subplotdir, 'dir')
        mkdir(subplotdir);
    end
    
    %% Set peakfinding parameters
    minprom = 0.01;
    mindist = 1.99; % in Hz
    
    %% Select subset of peaks (based on minprom and mindist parameters)
    
    % Use findpeaks on eigenvalue spectrum
    [~, Pfrex] = findpeaks(data(subno).GED_P.evals(:,1), data(subno).GED_P.frex, 'MinPeakProminence', minprom, 'MinPeakDistance', mindist);
    [~, IMfrex] = findpeaks(data(subno).GED_IM.evals(:,1), data(subno).GED_P.frex, 'MinPeakProminence', minprom, 'MinPeakDistance', mindist);
    
    % Compare found peak indices to set of original peak indices; construct
    % index
    Pmask = find(ismembertol(data(subno).peaks_P.frex, Pfrex));
    IMmask = find(ismembertol(data(subno).peaks_IM.frex, IMfrex));

    data(subno).Pmask = Pmask;
    data(subno).IMmask = IMmask;
    
    % Plot individual subject data:
    %% Plot spectral scanning spectra for P > IM and IM > P
    setPlotDefaults();
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(2,1,1);
    hold on;
    plot(data(subno).GED_P.frex, data(subno).GED_P.evals(:,1));
    plot(data(subno).peaks_P.frex(Pmask), data(subno).peaks_P.height(Pmask), 'rs');
    xlabel('Frequency (Hz)');
    ylabel('Power ratio (\lambda)');
    title('P > IM');
    xlim([min(data(subno).GED_P.frex) max(data(subno).GED_P.frex)]);
    
    subplot(2,1,2);
    hold on;
    plot(data(subno).GED_IM.frex, data(subno).GED_IM.evals(:,1));
    plot(data(subno).peaks_IM.frex(IMmask), data(subno).peaks_IM.height(IMmask), 'rs');
    xlabel('Frequency (Hz)');
    ylabel('Power ratio (\lambda)');
    title('IM > P');
    xlim([min(data(subno).GED_IM.frex) max(data(subno).GED_IM.frex)]);
    
    suptitle(['Eigenspectra for ' sublist{subno} ]);
    
    saveas(gca, [subplotdir '01_eigenspectra.png']);
    close;
    
    %% Plot topoplots at peaks
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    for i = 1:length(Pmask)
        subplot(2,max(length(Pmask), length(IMmask)), i);
        topoplot(data(subno).peaks_P.topo(Pmask(i),:), data(subno).chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
        title(['P, ' num2str(data(subno).peaks_P.frex(Pmask(i)),3) ' Hz, \lambda ' num2str(data(subno).GED_P.evals(data(subno).peaks_P.idx(Pmask(i)),1),3)]);
    end
    
    for i = 1:length(IMmask)
        subplot(2,max(length(Pmask), length(IMmask)), ...
            i+max(length(Pmask), length(IMmask)));
        topoplot(data(subno).peaks_IM.topo(IMmask(i),:), data(subno).chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
        title(['P, ' num2str(data(subno).peaks_IM.frex(IMmask(i)),3) ' Hz, \lambda ' num2str(data(subno).GED_IM.evals(data(subno).peaks_IM.idx(IMmask(i)),1),3)]);
    end
    
    suptitle(['Peak topoplots for ' sublist{subno} ]);
    
    saveas(gca, [subplotdir '02_peaktopoplots.png']);
    close;
    
    %% Plot component time series (ERPs) at peaks

    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    for i = 1:length(Pmask)
        subplot(length(Pmask), 1, i);
        plot(data(subno).peaks_P.tfd.metadata.times, mean(data(subno).peaks_P.compts(Pmask(i),:,:),3)); % mean over trials
        ylim([-3 3]*1e-4);
        hold on;
        % Add lines for FP window
        plot([1 1], ylim, 'k', 'LineWidth', 1.5);
        plot([1.8 1.8], ylim, 'r', 'LineWidth', 1.5);
        % Add lines for SP widnow
        plot([2.3 2.3], ylim, 'k', 'LineWidth', 1.5);
        plot([3.1 3.1], ylim, 'r', 'LineWidth', 1.5);
        % Add lines for IM window
        plot([4.6 4.6], ylim, 'k', 'LineWidth', 1.5);
        plot([8.1 8.1], ylim, 'r', 'LineWidth', 1.5);
        xlabel('Time (s)');
        title(['Peak comp at ' num2str(data(subno).peaks_P.frex(Pmask(i)), '%.2f') ' Hz (\lambda = ' num2str(data(subno).peaks_P.height(Pmask(i)), '%.2f') ')']);
    end
    suptitle(['P > IM peak component time series for ' sublist{subno}]);
    saveas(gca, [subplotdir '03_P_peak_ERP.png']);
    close;
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    for i = 1:length(IMmask)
        subplot(length(IMmask), 1, i);
        plot(data(subno).peaks_IM.tfd.metadata.times, mean(data(subno).peaks_IM.compts(IMmask(i),:,:),3)); % mean over trials
        ylim([-3 3]*1e-4);
        hold on;
        % Add lines for FP window
        plot([1 1], ylim, 'k', 'LineWidth', 1.5);
        plot([1.8 1.8], ylim, 'r', 'LineWidth', 1.5);
        % Add lines for SP widnow
        plot([2.3 2.3], ylim, 'k', 'LineWidth', 1.5);
        plot([3.1 3.1], ylim, 'r', 'LineWidth', 1.5);
        % Add lines for IM window
        plot([4.6 4.6], ylim, 'k', 'LineWidth', 1.5);
        plot([8.1 8.1], ylim, 'r', 'LineWidth', 1.5);
        xlabel('Time (s)');
        title(['Peak comp at ' num2str(data(subno).peaks_IM.frex(IMmask(i)), '%.2f') ' Hz (\lambda = ' num2str(data(subno).peaks_IM.height(IMmask(i)), '%.2f') ')']);
    end
    suptitle(['IM > P peak component time series for ' sublist{subno}]);
    saveas(gca, [subplotdir '04_IM_peak_ERP.png']);
    close;
    
    %% Plot TF decompositions at peaks
    
    clim = 6;
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    for i = 1:length(Pmask)
        subplot(max(3,length(Pmask)),1,i);
        contourf(data(subno).peaks_P.tfd.times2save, log10(data(subno).peaks_P.tfd.frex), squeeze(data(subno).peaks_P.tfd.tf(Pmask(i),:,:)), 40, 'linecolor', 'none');
        yticks(log10(logspace(log10(min(data(subno).peaks_P.tfd.frex)), log10(max(data(subno).peaks_P.tfd.frex)), 10)));
        yticklabels(compose('%.1f', 10.^yticks));
        hold on;
        ylabel('Frequency (Hz)');
        colorbar; colormap jet;
        caxis([-1 1]*clim);
        % Add lines for FP
        plot([1 1], ylim, 'k', 'LineWidth', 1.5);
        plot([1.8 1.8], ylim, 'r', 'LineWidth', 1.5);
        % Add lines for SP
        plot([2.3 2.3], ylim, 'k', 'LineWidth', 1.5);
        plot([3.1 3.1], ylim, 'r', 'LineWidth', 1.5);
        % Add lines for IM
        plot([4.6 4.6], ylim, 'k', 'LineWidth', 1.5);
        plot([8.1 8.1], ylim, 'r', 'LineWidth', 1.5);
        title(['Peak comp at ' num2str(data(subno).peaks_P.frex(Pmask(i)), '%.2f') ' Hz (\lambda = ' num2str(data(subno).peaks_P.height(Pmask(i)), '%.2f') ')']);
    end
    xlabel('Time (s)');
    suptitle(['P > IM TF decompositions for ' sublist{subno} '; FP, SP, IM onset black, offset red']);
    saveas(gca, [subplotdir '05_P_peak_TF.png']);
    close;
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    for i = 1:length(IMmask)
        subplot(max(3,length(IMmask)),1,i);
        contourf(data(subno).peaks_IM.tfd.times2save, log10(data(subno).peaks_IM.tfd.frex), squeeze(data(subno).peaks_IM.tfd.tf(IMmask(i),:,:)), 40, 'linecolor', 'none');
        yticks(log10(logspace(log10(min(data(subno).peaks_IM.tfd.frex)), log10(max(data(subno).peaks_IM.tfd.frex)), 10)));
        yticklabels(compose('%.1f', 10.^yticks));
        hold on;
        ylabel('Frequency (Hz)');
        colorbar; colormap jet;
        caxis([-1 1]*clim);
        % Add lines for FP
        plot([1 1], ylim, 'k', 'LineWidth', 1.5);
        plot([1.8 1.8], ylim, 'r', 'LineWidth', 1.5);
        % Add lines for SP
        plot([2.3 2.3], ylim, 'k', 'LineWidth', 1.5);
        plot([3.1 3.1], ylim, 'r', 'LineWidth', 1.5);
        % Add lines for IM
        plot([4.6 4.6], ylim, 'k', 'LineWidth', 1.5);
        plot([8.1 8.1], ylim, 'r', 'LineWidth', 1.5);
        title(['Peak comp at ' num2str(data(subno).peaks_IM.frex(IMmask(i)), '%.2f') ' Hz (\lambda = ' num2str(data(subno).peaks_IM.height(IMmask(i)), '%.2f') ')']);
    end
    xlabel('Time (s)');
    suptitle(['IM > P TF decompositions for ' sublist{subno} '; FP, SP, IM onset black, offset red']);
    saveas(gca, [subplotdir '06_IM_peak_TF.png']);
    close;

end

%% Compute KDF
% Compute kernel density plots based on individual peaks

allfrex = data(1).interp_P.frex;

% Restructure individual peak masks to apply to full set of peaks
Pmasks = [];
IMmasks = [];
for subno = 1:length(sublist)
    temp = zeros(length(data(subno).peaks_P.idx),1);
    temp(data(subno).Pmask) = 1;
    Pmasks = [Pmasks; temp];
    
    temp = zeros(length(data(subno).peaks_IM.idx),1);
    temp(data(subno).IMmask) = 1;
    IMmasks = [IMmasks; temp];
end
Pmasks = find(Pmasks);
IMmasks = find(IMmasks);

% Gather all peak data
Pfrex = cell2mat( arrayfun(@(c) c(:).peaks_P.frex, data, 'Uniform', 0));
IMfrex = cell2mat( arrayfun(@(c) c(:).peaks_IM.frex, data, 'Uniform', 0));

% Mask all peak data
Pfrex = Pfrex(Pmasks);
IMfrex = IMfrex(IMmasks);

% Apply kernel density smoothing:
GaussWidth = 2;
pdP = fitdist(Pfrex','Kernel','BandWidth',GaussWidth);
pdIM = fitdist(IMfrex','Kernel','BandWidth',GaussWidth);
kdfP = pdf(pdP,allfrex);
kdfIM = pdf(pdIM,allfrex);

% Apply peakfinding to kernel density:
peakthres = 0.01;
[kdfPpeakheight, kdfPpeakidx] = findpeaks(kdfP, 'MinPeakHeight', peakthres);
[kdfIMpeakheight, kdfIMpeakidx] = findpeaks(kdfIM, 'MinPeakHeight', peakthres);

%% Define peak widths

% Use Klimesch golden rule instead: ratio of 2.33 between center frequency
% and frequency band width (doi: 10.1111/ejn.14192)
ratio = 2.33;
dist_P = (allfrex(kdfPpeakidx) ./ ratio) ./ 2;
dist_IM = (allfrex(kdfIMpeakidx) ./ ratio) ./ 2;

lcrP = [];
for i = 1:length(kdfPpeakidx)
    left = allfrex(max(1, kdfPpeakidx(i))) - dist_P(i);
    center = allfrex(kdfPpeakidx(i));
    right = allfrex(kdfPpeakidx(i)) + dist_P(i);
    lcrP(i,:) = [left, center, right];
end

lcrIM = [];
for i = 1:length(kdfIMpeakidx)
    left = allfrex(max(1, kdfIMpeakidx(i))) - dist_IM(i);
    center = allfrex(kdfIMpeakidx(i));
    right = allfrex(kdfIMpeakidx(i)) + dist_IM(i);
    lcrIM(i,:) = [left, center, right];
end

%% Merge lcrIM for overlapping gamma peaks
gamma_merged = [min(lcrIM(5:6,1)), mean(lcrIM(5:6,2),1), max(lcrIM(5:6,3))];
lcrIM_orig = lcrIM;
lcrIM = [lcrIM(1:4,:); gamma_merged];

%% Compute average topoplots

disp('P > IM:');
topo2plot_P = {};
for i = 1:size(lcrP,1) % for each center frequency
    % Gather all components with peaks in range
    topo2plot = [];
    toposub = [];
    topofrex = [];
    for subno = 1:length(sublist)
        maskedfrex = data(subno).peaks_P.frex(data(subno).Pmask);
        maskedtopo = data(subno).peaks_P.topo(data(subno).Pmask,:);
        [~, closestpeakidx] = min(abs(maskedfrex - lcrP(i,2))); % grab peak that is closest to the center frequency
        closestpeak = maskedfrex(closestpeakidx); % convert index to frequency
        if closestpeak >= lcrP(i,1) && closestpeak <= lcrP(i,3) % check whether frequency is in range
            topo2plot = [topo2plot; maskedtopo(closestpeakidx,:)]; % if so, append to topography to topo2plot
            toposub = [toposub; subno];
            topofrex = [topofrex; closestpeak]; % append frequency for descriptive statistics
        end
    end
    % Compute grand average topoplot
    grandavg = mean(topo2plot,1);
    
    % Compute each topoplot's correlation with grand average, and flip
    for j = 1:size(topo2plot,1)
        topocorr = corrcoef(topo2plot(j,:), grandavg);
        flip = sign(topocorr(1,2));
        topo2plot(j,:) = topo2plot(j,:) * flip;
        % Store flip, for later flipping of component ERPs
        flipP(i,j) = 1;
    end
    
    % Save to a variable for later plotting next to TF
    topo2plot_P{i} = topo2plot;
    toposub_P{i} = toposub;
    
    % Debug: check how many comps are in range
    disp([ num2str(lcrP(i,2)) 'Hz, width ' num2str((lcrP(i,3) - lcrP(i,2)) * 2) 'Hz: ' num2str(size(topo2plot,1)) ' comps; freq mean = '  num2str(mean(topofrex)) ', SD = ' num2str(std(topofrex)) ]);
end

disp('IM > P');
topo2plot_IM = {};
for i = 1:size(lcrIM,1) % for each center frequency
    % Gather all components with peaks in range
    topo2plot = [];
    toposub = [];
    topofrex = [];
    for subno = 1:length(sublist)
        maskedfrex = data(subno).peaks_IM.frex(data(subno).IMmask);
        maskedtopo = data(subno).peaks_IM.topo(data(subno).IMmask,:);
        [~, closestpeakidx] = min(abs(maskedfrex - lcrIM(i,2))); % grab peak that is closest to the center frequency
        closestpeak = maskedfrex(closestpeakidx); % convert index to frequency
        if closestpeak >= lcrIM(i,1) && closestpeak <= lcrIM(i,3) % check whether frequency is in range
            topo2plot = [topo2plot; maskedtopo(closestpeakidx,:)]; % if so, append to topography to topo2plot
            toposub = [toposub; subno];
            topofrex = [topofrex; closestpeak]; % append frequency for descriptive statistics
        end        
    end
    % Compute grand average topoplot
    grandavg = mean(topo2plot,1);
    
    % Compute each topoplot's correlation with grand average, and flip
    for j = 1:size(topo2plot,1)
        topocorr = corrcoef(topo2plot(j,:), grandavg);
        flip = sign(topocorr(1,2));
        topo2plot(j,:) = topo2plot(j,:) * flip;
        % Store flip, for later flipping of component ERPs
        flipIM(i,j) = 1;
    end
    
    % Save to a variable for later plotting next to TF
    topo2plot_IM{i} = topo2plot;
    toposub_IM{i} = toposub;
    
    % check how many comps are in range
    disp([ num2str(lcrIM(i,2)) 'Hz, width ' num2str((lcrIM(i,3) - lcrIM(i,2)) * 2) 'Hz: ' num2str(size(topo2plot,1)) ' comps; freq mean = '  num2str(mean(topofrex)) ', SD = ' num2str(std(topofrex))  ]);
end

%% 
% Generating figure panels for the manuscript starts here.

%% Figure 4
% Spectral scanning averages and KDF plot
fig = 4;
status = mkdir([dirs.figs 'Figure ' num2str(fig)]);

setPlotDefaults();
set(0,'defaultLineMarkerSize',5);
set(0,'defaultAxesFontSize',10);

panel = 'A'; % Spectral scanning averages

ax = axgrid(10, 17.6, 1, .5, 1.5, .5, 2, 2, 'centimeters');

axes(ax(1));
hold on;
for subno = 1:length(sublist)
    plot(data(subno).GED_P.frex, data(subno).GED_P.evals(:,1), 'DisplayName', sublist{subno}, 'LineWidth', 1);
    plot(data(subno).peaks_P.frex(data(subno).Pmask), data(subno).peaks_P.height(data(subno).Pmask), '.k');
end
ylim([0 12]);
xticks(0:20:80);
xlim([2 80]);
xlabel('Frequency (Hz)');
ylabel('Highest \lambda');
legend show;
title('P > IM');
set(gca, 'box', 'off');

axes(ax(2));
hold on;
for subno = 1:length(sublist)
    plot(data(subno).GED_IM.frex, data(subno).GED_IM.evals(:,1), 'DisplayName', sublist{subno}, 'LineWidth', 1);
    plot(data(subno).peaks_IM.frex(data(subno).IMmask), data(subno).peaks_IM.height(data(subno).IMmask), '.k');
end
ylim([0 12]);
xticks(0:20:80);
xlim([2 80]);
ylabel('Highest \lambda');
set(gca, 'box', 'off');

axes(ax(3));
plot(allfrex, kdfP, 'color', 'k', 'LineWidth', 2); % full KDF
hold on;
plot(allfrex(kdfPpeakidx), kdfPpeakheight, 'sr'); % peaks on KDF
% Plot lines at x Hz to the left and right
for i = 1:length(kdfPpeakidx)
    ptA = [lcrP(i,1), 0];
    ptB = [lcrP(i,3), 0];
    ptC = [lcrP(i,1), kdfP(dsearchn(allfrex(:), lcrP(i,1)))];
    ptD = [lcrP(i,3), kdfP(dsearchn(allfrex(:), lcrP(i,3)))];
    
    % Plot lines
    plot([ptA(1) ptC(1)], [ptA(2) ptC(2)], '--k')
    plot([ptB(1) ptD(1)], [ptB(2) ptD(2)], '--k')

    % Plot patch instead of lines
%     X = [ptA(1) allfrex(dsearchn(allfrex(:), lcrP(i,1)):dsearchn(allfrex(:), lcrP(i,3))) ptB(1)];
%     Y = [ptA(2) kdfP(dsearchn(allfrex(:), lcrP(i,1)):dsearchn(allfrex(:), lcrP(i,3))) ptB(2)];
%     patch(X, Y, [0 0.4470 0.7410]);
    
    text(lcrP(i,2), kdfPpeakheight(i)+max(caxis)*0.003, 1, [num2str(lcrP(i,2), '%.1f') ' Hz'], 'FontSize', 10);
end
title('P > IM peaks per frequency');
xlabel('Frequency (Hz)');
ylabel('Kernel density estimate');
ylim([0 0.04]);
yticks([0:0.01:0.04]);
xticks(0:20:80);
xlim([2 80]);
set(gca, 'box', 'off');

% IM > P
axes(ax(4));
plot(allfrex, kdfIM, 'color', 'k', 'LineWidth', 2); % full KDF
hold on;
plot(allfrex(kdfIMpeakidx), kdfIMpeakheight, 'sr'); % peaks on KDF
% Plot lines at x Hz to the left and right
for i = 1:length(kdfIMpeakidx)
    ptA = [lcrIM_orig(i,1), 0];
    ptB = [lcrIM_orig(i,3), 0];
    ptC = [lcrIM_orig(i,1), kdfIM(dsearchn(allfrex(:), lcrIM_orig(i,1)))];
    ptD = [lcrIM_orig(i,3), kdfIM(dsearchn(allfrex(:), lcrIM_orig(i,3)))];
    
    plot([ptA(1) ptC(1)], [ptA(2) ptC(2)], '--k')
    plot([ptB(1) ptD(1)], [ptB(2) ptD(2)], '--k')
    
    % Plot patch instead of lines
%     X = [ptA(1) allfrex(dsearchn(allfrex(:), lcrIM(i,1)):dsearchn(allfrex(:), lcrIM(i,3))) ptB(1)];
%     Y = [ptA(2) kdfIM(dsearchn(allfrex(:), lcrIM(i,1)):dsearchn(allfrex(:), lcrIM(i,3))) ptB(2)];
%     patch(X, Y, [0 0.4470 0.7410]);
    
    text(lcrIM_orig(i,2), kdfIMpeakheight(i)+max(caxis)*0.003, 1, [num2str(lcrIM_orig(i,2), '%.1f') ' Hz'], 'FontSize', 10);
end
title('IM > P peaks per frequency');
xlabel('Frequency (Hz)');
ylabel('Kernel density estimate');
ylim([0 0.04]);
yticks([0:0.01:0.04]);
xticks(0:20:80);
xlim([2 80]);
set(gca, 'box', 'off');

saveas(gcf, [dirs.figs 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%% Figure 5 - topoplots and TF decompositions
fig = 5;
status = mkdir([dirs.figs 'Figure ' num2str(fig)]);

setPlotDefaults();
set(0,'defaultAxesFontSize',12);

%%
panel = 'A1'; % P > IM topoplots

figure;
for i = 1:length(lcrP)
    subplot(length(lcrP), 1, i);
    topoplot(mean(topo2plot_P{i},1), data(subno).chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
end

saveas(gcf, [dirs.figs 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
panel = 'A2'; % P > IM TF decomps

ax = axgrid(12, 8.5, .25, .25, .5, .5, length(lcrP), 1, 'centimeters');

dBs = {}; % collect clims for later outputting
for i = 1:length(lcrP)
    % Gather closest component TF decomp with peak in range
    comptf2plot = [];
    for subno = 1:length(sublist)
        maskedfrex = data(subno).peaks_P.frex(data(subno).Pmask);
        maskedtf = data(subno).peaks_P.tfd.tf(data(subno).Pmask,:,:);
        [~, closestpeakidx] = min(abs(maskedfrex - lcrP(i,2))); % grab peak that is closest to the center frequency
        closestpeak = maskedfrex(closestpeakidx); % convert index to frequency
        if closestpeak >= lcrP(i,1) && closestpeak <= lcrP(i,3) % check whether frequency is in range
            comptf2plot = [comptf2plot; maskedtf(closestpeakidx,:,:)]; % if so, append tf to tf2plot
        end      
    end
    
    % Average TF spectra
    grandavg = mean(comptf2plot,1);
    
    % Plot TF
    axes(ax(i));
    contourf(data(subno).peaks_P.tfd.times2save, log10(data(subno).peaks_P.tfd.frex), squeeze(grandavg), 100, 'linecolor', 'none');
    yticks(log10(logspace(log10(min(data(subno).peaks_P.tfd.frex)), log10(max(data(subno).peaks_P.tfd.frex)), 5)));
    yticklabels(compose('%.0f', 10.^yticks));
    hold on;
    ylabel('Frequency (Hz)');
    colormap jet;
    
    % Add horizontal line at peak frequency
    plot(xlim, [log10(lcrP(i,2)) log10(lcrP(i,2))], 'k--', 'LineWidth', 1);
    
    % Only show xticklabels for last plot
    if i ~= length(lcrP)
        xticklabels([]);
    end
    
    % Set color limits to match window (P, IM) with lowest max power; blow out
    % other window
    FPwidx = dsearchn(data(subno).peaks_P.tfd.times2save(:), [1 1.8]');
    SPwidx = dsearchn(data(subno).peaks_P.tfd.times2save(:), [2.3 3.1]');
    IMwidx = dsearchn(data(subno).peaks_IM.tfd.times2save(:), [4.8 7.9]'); % shortened a bit because of large deflections following imagery
    Pwin = [grandavg(:,:,FPwidx(1):FPwidx(2)), grandavg(:,:,SPwidx(1):SPwidx(2))];
    IMwin = grandavg(:,:,IMwidx(1):IMwidx(2));
    clim = min(max(abs(Pwin(:))), max(abs(IMwin(:)))) * 1.3; % extract maximum; multiply by scaling factor to look nicer
    clim = ceil(10*clim) / 10; % round up to nearest 10th dB
    caxis([-1 1]*clim);
    
    dBs{i} = num2str(clim);
    
    % Add lines for FP
    plot([1 1], ylim, 'k', 'LineWidth', 1.5);
    plot([1.8 1.8], ylim, 'k--', 'LineWidth', 1.5);
    % Add lines for SP
    plot([2.3 2.3], ylim, 'k', 'LineWidth', 1.5);
    plot([3.1 3.1], ylim, 'k--', 'LineWidth', 1.5);
    % Add lines for IM
    plot([4.6 4.6], ylim, 'k', 'LineWidth', 1.5);
    plot([8.1 8.1], ylim, 'k--', 'LineWidth', 1.5);
%     title([ num2str(lcrP(i,2)) ' Hz, n = ' num2str(size(comptf2plot,1)) ]);    
end

xlabel('Time (s)');
suptitle('P > IM average TF decompositions');
% suptitle(['dBs: ' dBs]);

saveas(gcf, [dirs.figs 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
panel = 'B1'; % IM > P topoplots

figure;
for i = 1:length(lcrIM)
    subplot(length(lcrIM), 1, i);
    topoplot(mean(topo2plot_IM{i},1), data(subno).chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
end

saveas(gcf, [dirs.figs 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;

%%
panel = 'B2'; % IM > P TFs

height = 3 * length(lcrIM); % changed to get TF subplots same vertical size as P > IM
ax = axgrid(height, 8.5, .25, .25, .5, .5, length(lcrIM), 1, 'centimeters');

dBs = {}; % collect clims for later outputting
for i = 1:length(lcrIM)
    % Gather closest component TF decomp with peak in range
    comptf2plot = [];
    for subno = 1:length(sublist)
        maskedfrex = data(subno).peaks_IM.frex(data(subno).IMmask);
        maskedtf = data(subno).peaks_IM.tfd.tf(data(subno).IMmask,:,:);
        [~, closestpeakidx] = min(abs(maskedfrex - lcrIM(i,2))); % grab peak that is closest to the center frequency
        closestpeak = maskedfrex(closestpeakidx); % convert index to frequency
        if closestpeak >= lcrIM(i,1) && closestpeak <= lcrIM(i,3) % check whether frequency is in range
            comptf2plot = [comptf2plot; maskedtf(closestpeakidx,:,:)]; % if so, append tf to tf2plot
        end
    end
    
    % Average TF spectra
    grandavg = mean(comptf2plot,1);
    
    % Plot TF
    axes(ax(i));
    contourf(data(subno).peaks_IM.tfd.times2save, log10(data(subno).peaks_IM.tfd.frex), squeeze(grandavg), 100, 'linecolor', 'none');
    yticks(log10(logspace(log10(min(data(subno).peaks_IM.tfd.frex)), log10(max(data(subno).peaks_IM.tfd.frex)), 5)));
    yticklabels(compose('%.0f', 10.^yticks));
    hold on;
    ylabel('Frequency (Hz)');
    colormap jet;
    
    % Add horizontal line at peak frequency
    plot(xlim, [log10(lcrIM(i,2)) log10(lcrIM(i,2))], 'k--', 'LineWidth', 1);
    
    % Only show xticklabels for last plot
    if i ~= length(lcrIM)
        xticklabels([]);
    end
    
    % Set color limits to match window (P, IM) with lowest max power; blow out
    % other window
    FPwidx = dsearchn(data(subno).peaks_P.tfd.times2save(:), [1 1.8]');
    SPwidx = dsearchn(data(subno).peaks_P.tfd.times2save(:), [2.3 3.1]');
    IMwidx = dsearchn(data(subno).peaks_IM.tfd.times2save(:), [4.8 7.9]'); % shortened a bit because of large deflections following perception
    Pwin = [grandavg(:,:,FPwidx(1):FPwidx(2)), grandavg(:,:,SPwidx(1):SPwidx(2))];
    IMwin = grandavg(:,:,IMwidx(1):IMwidx(2));
    clim = min(max(abs(Pwin(:))), max(abs(IMwin(:)))) * 1.3; % extract maximum; multiply by scaling factor to look nicer
    clim = ceil(10*clim) / 10; % round up to nearest 10th dB
    caxis([-1 1]*clim);
    
    dBs{i} = num2str(clim);
    
    % Add lines for FP
    plot([1 1], ylim, 'k', 'LineWidth', 1.5);
    plot([1.8 1.8], ylim, 'k--', 'LineWidth', 1.5);
    % Add lines for SP
    plot([2.3 2.3], ylim, 'k', 'LineWidth', 1.5);
    plot([3.1 3.1], ylim, 'k--', 'LineWidth', 1.5);
    % Add lines for IM
    plot([4.6 4.6], ylim, 'k', 'LineWidth', 1.5);
    plot([8.1 8.1], ylim, 'k--', 'LineWidth', 1.5);
%     title([ num2str(lcrIM(i,2)) ' Hz, n = ' num2str(size(comptf2plot,1)) ]);    
end

xlabel('Time (s)');
% suptitle(['dBs: ' dBs]);

saveas(gcf, [dirs.figs 'Figure ' num2str(fig) '/' num2str(fig) panel '.pdf']);
close;


%% Supplementary figures: plot all components going into alpha/gamma peak
setPlotDefaults();
set(0,'defaultAxesFontSize',12);

% Select alpha and gamma peak
Ppeaks2plot = [1 4];
IMpeaks2plot = [2 5];
figs = {'S2' 'S3'; 'S4' 'S5'}; % S2 = P > IM alpha, S3 = P > IM gamma, S4 = IM > P alpha, S5 = IM > P gamma

for i = 1:2 % alpha, gamma
    Ppeak = Ppeaks2plot(i);
    IMpeak = IMpeaks2plot(i);
    
    % Gather all component TF decomps with peaks in range
    Ptf2plot = [];
    Psub = [];
    PHz = [];
    for subno = 1:length(sublist)
        maskedfrex = data(subno).peaks_P.frex(data(subno).Pmask);
        maskedtf = data(subno).peaks_P.tfd.tf(data(subno).Pmask,:,:);
        [~, closestpeakidx] = min(abs(maskedfrex - lcrP(Ppeak,2))); % grab peak that is closest to the center frequency
        closestpeak = maskedfrex(closestpeakidx); % convert index to frequency
        if closestpeak >= lcrP(Ppeak,1) && closestpeak <= lcrP(Ppeak,3) % check whether frequency is in range
            Ptf2plot = [Ptf2plot; maskedtf(closestpeakidx,:,:)]; % if so, append tf to tf2plot
            Psub = [Psub subno];
            PHz = [PHz closestpeak];
        end
    end  
    
    IMtf2plot = [];
    IMsub = [];
    IMHz = [];
    for subno = 1:length(sublist)
        maskedfrex = data(subno).peaks_IM.frex(data(subno).IMmask);
        maskedtf = data(subno).peaks_IM.tfd.tf(data(subno).IMmask,:,:);
        [~, closestpeakidx] = min(abs(maskedfrex - lcrIM(IMpeak,2))); % grab peak that is closest to the center frequency
        closestpeak = maskedfrex(closestpeakidx); % convert index to frequency
        if closestpeak >= lcrIM(IMpeak,1) && closestpeak <= lcrIM(IMpeak,3) % check whether frequency is in range
            IMtf2plot = [IMtf2plot; maskedtf(closestpeakidx,:,:)]; % if so, append tf to tf2plot
            IMsub = [IMsub subno];
            IMHz = [IMHz closestpeak];
        end
    end 
    
    fig = figs{1,i};
    status = mkdir([dirs.figs 'Supplementary Figure ' num2str(fig)]);
    panel = 'A'; % P > IM topoplots
    
    figure;
    for p = 1:min(16, size(Ptf2plot,1)) % for each component to plot (max. 16, otherwise make new figure - preventing memory errors)
        subplot(4,4,p);
        topoplot(topo2plot_P{Ppeaks2plot(i)}(p,:), data(subno).chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
        title([sublist(Psub(p)) num2str(PHz(p), '%.1f') ' Hz']);
    end
    saveas(gcf, [dirs.figs 'Supplementary Figure ' num2str(fig) '/' num2str(fig) panel '_1.pdf']);
    close;
    if size(Ptf2plot,1) > 16
        figure;
        for p = 17:size(Ptf2plot,1)
            subplot(4,4,mod(p,16));
            topoplot(topo2plot_P{Ppeaks2plot(i)}(p,:), data(subno).chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
            title([sublist(Psub(p)) num2str(PHz(p), '%.1f') ' Hz']);
        end
        saveas(gcf, [dirs.figs 'Supplementary Figure ' num2str(fig) '/' num2str(fig) panel '_2.pdf']);
        close;
    end
    
    fig = figs{2,i};
    status = mkdir([dirs.figs 'Supplementary Figure ' num2str(fig)]);
    panel = 'A'; % IM > P topoplots
    
    figure;
    for p = 1:min(16, size(IMtf2plot,1)) % for each component to plot (max. 16, otherwise make new figure - preventing memory errors)
        subplot(4,4,p);
        topoplot(topo2plot_IM{IMpeaks2plot(i)}(p,:), data(subno).chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
        title([sublist(IMsub(p)) num2str(IMHz(p), '%.1f') ' Hz']);
    end
    saveas(gcf, [dirs.figs 'Supplementary Figure ' num2str(fig) '/' num2str(fig) panel '_1.pdf']);
    close;
    if size(IMtf2plot,1) > 16
        figure;
        for p = 17:size(IMtf2plot,1)
            subplot(4,4,mod(p,16));
            topoplot(topo2plot_IM{IMpeaks2plot(i)}(p,:), data(subno).chanlocs, 'numcontour',0,'electrodes','off','shading','interp');
            title([sublist(IMsub(p)) num2str(IMHz(p), '%.1f') ' Hz']);
        end
        saveas(gcf, [dirs.figs 'Supplementary Figure ' num2str(fig) '/' num2str(fig) panel '_2.pdf']);
        close;
    end
       
    fig = figs{1,i};
    panel = 'B'; % P > IM TF decomps
   
    for p = 1:size(Ptf2plot,1)
        if mod(p-1,5)+1 == 1
            figure;
            height = 3 * 5; % make all TF subplots same vertical size (2.5 cm)
            ax = axgrid(height, 8.5, .25, .25, .5, .5, 5, 1, 'centimeters');
            dBs = {}; % collect clims for later outputting
        end
        axes(ax(mod(p-1,5)+1));
        contourf(data(subno).peaks_P.tfd.times2save, log10(data(subno).peaks_P.tfd.frex), squeeze(Ptf2plot(p,:,:)), 100, 'linecolor', 'none');
        yticks(log10(logspace(log10(min(data(subno).peaks_P.tfd.frex)), log10(max(data(subno).peaks_P.tfd.frex)), 5)));
        yticklabels(compose('%.0f', 10.^yticks));
        hold on;
        ylabel('Frequency (Hz)');
        colormap jet;
        
        %     Add horizontal line at frequency
        plot(xlim, [log10(PHz(p)) log10(PHz(p))], 'k--', 'LineWidth', 1);
        
        % Only show xticklabels for last plot
        if p ~= size(Ptf2plot,1)
            xticklabels([]);
        end
        
        % Set color limits to match window (P, IM) with lowest max power; blow out
        % other window
        FPwidx = dsearchn(data(subno).peaks_P.tfd.times2save(:), [1 1.8]');
        SPwidx = dsearchn(data(subno).peaks_P.tfd.times2save(:), [2.3 3.1]');
        IMwidx = dsearchn(data(subno).peaks_IM.tfd.times2save(:), [4.8 7.9]'); % shortened a bit because of large deflections following perception
        Pwin = [Ptf2plot(p,:,FPwidx(1):FPwidx(2)), Ptf2plot(p,:,SPwidx(1):SPwidx(2))];
        IMwin = Ptf2plot(p,:,IMwidx(1):IMwidx(2));
        clim = min(max(abs(Pwin(:))), max(abs(IMwin(:)))) * 1.3; % extract maximum; multiply by scaling factor to look nicer
        clim = ceil(10*clim) / 10; % round up to nearest 10th dB
        caxis([-1 1]*clim);
        
        dBs{mod(p-1,5)+1} = num2str(clim);
       
         % Add lines for FP
        plot([1 1], ylim, 'k', 'LineWidth', 1.5);
        plot([1.8 1.8], ylim, 'k--', 'LineWidth', 1.5);
        % Add lines for SP
        plot([2.3 2.3], ylim, 'k', 'LineWidth', 1.5);
        plot([3.1 3.1], ylim, 'k--', 'LineWidth', 1.5);
        % Add lines for IM
        plot([4.6 4.6], ylim, 'k', 'LineWidth', 1.5);
        plot([8.1 8.1], ylim, 'k--', 'LineWidth', 1.5);
        %     title([ num2str(lcrIM(i,2)) ' Hz, n = ' num2str(size(comptf2plot,1)) ]);
        
        if mod(p-1,5)+1 == 5 || p == size(Ptf2plot,1)
            xlabel('Time (s)');
            % suptitle(['dBs: ' dBs]);
            
            saveas(gcf, [dirs.figs 'Supplementary Figure ' num2str(fig) '/' num2str(fig) panel '_' num2str(ceil(p/5)) '.pdf']);
            close;
        end
    end 
    
    fig = figs{2,i};
    panel = 'B'; % IM > P TF decomps
   
    for p = 1:size(IMtf2plot,1)
        if mod(p-1,5)+1 == 1
            figure;
            height = 3 * 5; % make all TF subplots same vertical size (2.5 cm)
            ax = axgrid(height, 8.5, .25, .25, .5, .5, 5, 1, 'centimeters');
            dBs = {}; % collect clims for later outputting
        end
        axes(ax(mod(p-1,5)+1));
        contourf(data(subno).peaks_IM.tfd.times2save, log10(data(subno).peaks_IM.tfd.frex), squeeze(IMtf2plot(p,:,:)), 100, 'linecolor', 'none');
        yticks(log10(logspace(log10(min(data(subno).peaks_IM.tfd.frex)), log10(max(data(subno).peaks_IM.tfd.frex)), 5)));
        yticklabels(compose('%.0f', 10.^yticks));
        hold on;
        ylabel('Frequency (Hz)');
        colormap jet;
        
        % Add horizontal line at frequency
        plot(xlim, [log10(IMHz(p)) log10(IMHz(p))], 'k--', 'LineWidth', 1);
        
        % Only show xticklabels for last plot
        if p ~= size(IMtf2plot,1)
            xticklabels([]);
        end
        
        % Set color limits to match window (P, IM) with lowest max power; blow out
        % other window
        FPwidx = dsearchn(data(subno).peaks_P.tfd.times2save(:), [1 1.8]');
        SPwidx = dsearchn(data(subno).peaks_P.tfd.times2save(:), [2.3 3.1]');
        IMwidx = dsearchn(data(subno).peaks_IM.tfd.times2save(:), [4.8 7.9]'); % shortened a bit because of large deflections following perception
        Pwin = [IMtf2plot(p,:,FPwidx(1):FPwidx(2)), IMtf2plot(p,:,SPwidx(1):SPwidx(2))];
        IMwin = IMtf2plot(p,:,IMwidx(1):IMwidx(2));
        clim = min(max(abs(Pwin(:))), max(abs(IMwin(:)))) * 1.3; % extract maximum; multiply by scaling factor to look nicer
        clim = ceil(10*clim) / 10; % round up to nearest 10th dB
        caxis([-1 1]*clim);
        
        dBs{mod(p-1,5)+1} = num2str(clim);
       
        % Add lines for FP
        plot([1 1], ylim, 'k', 'LineWidth', 1.5);
        plot([1.8 1.8], ylim, 'k--', 'LineWidth', 1.5);
        % Add lines for SP
        plot([2.3 2.3], ylim, 'k', 'LineWidth', 1.5);
        plot([3.1 3.1], ylim, 'k--', 'LineWidth', 1.5);
        % Add lines for IM
        plot([4.6 4.6], ylim, 'k', 'LineWidth', 1.5);
        plot([8.1 8.1], ylim, 'k--', 'LineWidth', 1.5);
        % title([ num2str(lcrIM(i,2)) ' Hz, n = ' num2str(size(comptf2plot,1)) ]);
        
        if mod(p-1,5)+1 == 5 || p == size(IMtf2plot,1)
            xlabel('Time (s)');
            % suptitle(['dBs: ' dBs]);
            
            saveas(gcf, [dirs.figs 'Supplementary Figure ' num2str(fig) '/' num2str(fig) panel '_' num2str(ceil(p/5)) '.pdf']);
            close;
        end
    end 
    
end % end loop over alpha, gamma

%% Debug: Plot individual component topoplots going into peak average
% for i = 1:length(lcrP) % for each peak
%     figure;
%     for j = 1:size(topo2plot_P{i},1) % for each topo
%         subplot(ceil(sqrt(size(topo2plot_P{i},1))), ceil(sqrt(size(topo2plot_P{i},1))), j);
%         topoplot(topo2plot_P{i}(j,:), data(subno).chanlocs,  'numcontour',0,'electrodes','off','shading','interp'); % plot that topo
%         title(sublist(toposub_P{i}(j))); % show subject name in title
%     end
%     suptitle(['Individual topographies in P > IM peak at ' num2str(lcrP(i,2)) ' Hz']);
%     saveas(gca, [dirs.emp_plots '06_allsubject_topoplots_P_' num2str(lcrP(i,2)) 'Hz.png']);
%     close;
% end
% 
% for i = 1:length(lcrIM) % for each peak
%     figure;
%     for j = 1:size(topo2plot_IM{i},1) % for each topo
%         subplot(ceil(sqrt(size(topo2plot_IM{i},1))), ceil(sqrt(size(topo2plot_IM{i},1))), j);
%         topoplot(topo2plot_IM{i}(j,:), data(subno).chanlocs,  'numcontour',0,'electrodes','off','shading','interp'); % plot that topo
%         title(sublist(toposub_IM{i}(j))); % show subject name in title
%     end
%     suptitle(['Individual topographies in IM > P peak at ' num2str(lcrIM(i,2)) ' Hz']);
%     saveas(gca, [dirs.emp_plots '07_allsubject_topoplots_IM_' num2str(lcrIM(i,2)) 'Hz.png']);
%     close;
% end
