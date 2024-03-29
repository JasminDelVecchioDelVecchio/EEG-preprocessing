
function [EEGCleaned, EEGsteps] = MainCleaningPipeline(EEG, ID, pathRes)

%Input:
% - EEG: data to be cleaned in eeglab struct format
% - subj_name: string with subject's name 
% - pathRes: path to results' folder 

addCleaningPipelinePaths();

DefParam = DefaultParameters;

if nargin < 3
    pathRes = uigetdir(currentDir, 'Select folder to store results');
end

if nargin < 2
    ID = input('Please enter subject-ID: ', 's');
    pathRes = uigetdir(currentDir, 'Select folder to store results');
end

nfi = 1;

for iFile = 1:nfi

    %Exclude TENS artifact if present (used for synhcronization purposes)
    try
        d = 2*EEG.srate;
        ind_tens = find(ismember({EEG.event.type}, 'TENStrigger'));
        assert(length(ind_tens) == 2, 'need to have exactly two TENS triggers')
        EEG = pop_select(EEG, 'point', [EEG.event(ind_tens).latency] + [d -3*d]);
    catch
        warning('No TENS detected.')
    end


    EEG.badchan = false(1,EEG.nbchan);
    EEG.bad_segment = false(1,EEG.pnts*EEG.trials);
    EEG.bad_epoch = [];
    
    % Select only EEG channels
    iChannels =  find(strcmp({EEG.chanlocs.type},'EEG'));
    excludeChannels = [setdiff(1:EEG.nbchan,iChannels)];

    if iFile==1
        EEGCleaned = EEG;
        EEGmerged = pop_select(EEG, 'nochannel', excludeChannels);
        EEGmerged.badchan(excludeChannels) = [];
    end

    ecgchannels = find(cellfun(@(x)~isempty(x),strfind({EEG.chanlocs.type},'ECG')));
    ECG = pop_select(EEG,'channel', ecgchannels);
end

% Start cleaning pipeline
fprintf('Start cleaning pipeline...')

[EEGmerged, EEGsteps] = RunCleaningPipeline(EEGmerged, ECG);


if ~isempty(EEGmerged)
    % Update original EEG struct with cleaned EEG data
    EEGCleaned.data(iChannels,:) = EEGmerged.data(:,1:size(EEGCleaned.data,2));
    EEGCleaned.bad_segment = EEGmerged.bad_segment(1:size(EEGCleaned.data,2));
    EEGCleaned.badchan(iChannels) = EEGmerged.badchan;
    EEGCleaned.bad_epoch = [];
    EEGCleaned.icachansind = find(not(EEGmerged.badchan));
    EEGCleaned.icaweights = EEGmerged.icaweights;
    EEGCleaned.icasphere = EEGmerged.icasphere;
    EEGCleaned.icawinv    = EEGmerged.icawinv; % a priori same result as inv
    if ~isempty(EEGmerged.reject)
        EEGCleaned.reject.gcompreject = EEGmerged.reject.gcompreject;
    end

    delete(gcp('nocreate'))
end

if DefParam.save_fig == 1
    % Parameters
    load("cm17.mat") % custom colormaps
    frqs = 1:EEG.srate/2; % frequency range [1:F_Nyquist]
    bands = {[3 8], [8 12], [13 35],[35 80]}; % EEG bands of interest 
    bandnames = {'theta', 'alpha', 'beta','gamma'}; 

    % select only EEG channels for plotting - Raw data
    iChannels =  find(strcmp({EEG.chanlocs.type},'EEG')); 
    excludeChannels = [setdiff(1:EEG.nbchan,iChannels)];
    EEGplot = pop_select(EEG, 'nochannel', excludeChannels);
    EEGplot.badchan(excludeChannels) = [];
    
    % select only EEG channels for plotting - Cleaned data
    iChannels =  find(strcmp({EEGCleaned.chanlocs.type},'EEG')); 
    EEGplotcleaned = EEGCleaned;
    EEGplotcleaned.data = EEGCleaned.data(iChannels,:);
    EEGplotcleaned.chanlocs = EEGCleaned.chanlocs(1,iChannels);
    EEGplotcleaned.nbchan = length(iChannels);
    
    h = figure; pop_spectopo(EEGplot, 1, [0 EEGplot.xmax*1000], 'EEG' , 'percent', 100, 'freq', [10 20], ...
        'freqrange',[0 EEGplot.srate./2],'electrodes','on');
    title('Raw EEG','Position',[-42 48 1]);
    set(h,'windowstate','maximized');
    exportgraphics(h,[pathRes,ID,'_raw_spectra.png'],"Resolution",300)
    savefig([pathRes,ID,'raw_spectra']);
    close all;

    h = figure; pop_spectopo(EEGplotcleaned, 1, [0 EEGplotcleaned.xmax*1000], 'EEG' , 'percent', 100, 'freq', [10 20], ...
        'freqrange',[0 EEGplot.srate/2], 'electrodes','on');
    title('Cleaned EEG', 'Position', [-42 48 1]);
    exportgraphics(h, [pathRes, ID, '_cleaned_spectra.png'], "Resolution", 300)
    savefig([pathRes,ID,'cleaned_spectra']);
    close all;

    % Topoplot after each cleaning step
    prepro_steps = EEGsteps.Properties.VariableNames; % cell array containing preprocessing implemented steps
    
    sum_pow = []; psd = [];
    for istep = 1:length(prepro_steps)
        if iscell(EEGsteps.(istep))
            eeg = EEGsteps.(istep){1,1};
        else
            eeg = EEGsteps.(istep);
        end
        psd(:, :, istep) = pwelch(eeg.data', eeg.srate, eeg.srate/2, 1:eeg.srate/2, eeg.srate);
        for ib = 1:length(bands)
            frq_inds = frqs >= bands{ib}(1) & frqs <= bands{ib}(2);
            sum_pow(ib, :, istep) = squeeze(sum(psd(frq_inds, :, istep),1,'omitmissing'));
        end
    end
    
% Calculate the total power across all bands for each cleaning step
total_power = sum(sum_pow, 1);

% Calculate the relative power in each band for each cleaning step
relative_power = sum_pow ./ total_power;

idx = 1;
h = figure;
for istep = 1:size(relative_power, 3)
    for ib = 1:size(relative_power, 1)
        ax = subplot(7, 4, idx);
        power2plot = relative_power(ib, ~eeg.badchan, istep);
        topoplot(power2plot, eeg.chanlocs(~eeg.badchan));
        colormap(cm17), colorbar;
        %Using the 10th and 90th percentiles to set the color limits (clim) 
        % offers a good balance between including most of the data range while 
        % minimizing the influence of outliers.
        cmax = max(quantile(power2plot, 0.9));
        cmin = min(quantile(power2plot, 0.1));
        clim([cmin cmax]);
        % Add row title
        if ib == 1
            text(-2.8, 0.13, prepro_steps{istep}, 'FontSize', 18, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end

        % Add column title
        if istep == 1
            title(bandnames{ib}); % Use ib for column titles
        end
        set(ax, 'FontSize', 18);
        idx = idx + 1;
    end
end
set(h,'WindowState','Maximized');
exportgraphics(h,[pathRes,ID,'_topoplots.png'],"Resolution",300)
savefig([pathRes,ID,'topoplots']);
close all;
end

if DefParam.save_steps == 0
    EEGsteps = [];
end
% removeCleaningPipelinePaths();

end
