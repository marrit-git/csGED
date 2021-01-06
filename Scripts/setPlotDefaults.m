function [] = setPlotDefaults()
%SETPLOTDEFAULTS Call this to set a bunch of MATLAB defaults for nicer
%plots.

    %% Set MATLAB plotting defaults
    
    lw = 1.5;   % line width
    msz = 15;    % marker size
    fsz = 16;   % font size
    alw = 0.75; % axes line width
    tlen = [0.01 0.01]; % tick length
    
    width = 8.5;  % width in cm
    height = 6; % height in cm
    % Note: J Neurosci requests 8.5 cm, 11.6 cm, or 17.6 cm (1, 1.5, 2 columns) as max width
    
    %% Propagate MATLAB plotting defaults
    set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
    set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
    set(0,'defaultAxesFontSize',fsz);
    set(0, 'defaultAxesLineWidth', alw);
    set(0, 'defaultAxesTickLength', tlen);
    
    % Set default renderer to vector format outputs. Note: doesn't support
    % alpha channels. Commented out now; generating for presentation, not pub-quality
    % plots.
    set(0, 'defaultFigureRenderer', 'painters')
    
    % Set the default Size for display
    defpos = get(0,'defaultFigurePosition');
    set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
    
    % Set the defaults for saving/printing to a file
    set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
    set(0,'defaultFigurePaperUnits','centimeters');
    defsize = get(gcf, 'PaperSize');
    left = (defsize(1)- width)/2;
    bottom = (defsize(2)- height)/2;
    defsize = [left, bottom, width, height];
    set(0, 'defaultFigurePaperPosition', defsize);
    
    close; % close ghost figure
end

