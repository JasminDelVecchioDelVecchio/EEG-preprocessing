
function [EEGCleaned, EEGsteps] = MainCleaningPipeline(EEG,subj_name)

%Input:
% - EEG: data to be cleaned in eeglab struct format

addCleaningPipelinePaths();

DefParam = DefaultParameters;

nfi = 1;

for iFile = 1:nfi

    %Exclude TENS
    try
        d = 2*EEG.srate;
        ind_tens = find(ismember({EEG.event.type}, 'TENStrigger'));
        assert(length(ind_tens) == 2, 'need to have exactly two TENS triggers')
        %         ind_tens = ind_tens + [1 -1]; %one event after and before TENS
        EEG = pop_select(EEG, 'point', [EEG.event(ind_tens).latency] + [d -3*d]);
    catch
        warning('No TENS detected.')
    end

    EEG.badchan = false(1,EEG.nbchan);
    EEG.bad_segment = false(1,EEG.pnts*EEG.trials);
    EEG.bad_epoch = [];

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
% [EEGmerged, EEGsteps] = RunCleaningPipeline(EEGmerged, ECG, 'EpochingParams',EpochingParams,'DetrendParams',DetrendParams,... %*
%     'FilterParams',FilterParams,'EOGParams',EOGParams,'ECGParams', ECGParams,'PrepParams',PrepParams,'CRDParams',CRDParams,...
%     'FastICAParams',FastICAParams,'prunedata',prunedata,'performICA',performICA,...
%     'subcompICA',subcompICA,'interpbadchannels',interpbadchannels,'save_fig',save_fig,'save_steps',save_steps);

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
    load("cm17.mat") % Low custom colormaps
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
    EEGplotcleaned.data = EEGCleaned.data(iChannels,:);
    EEGplotcleaned.chanlocs = EEGCleaned.chanlocs(1,iChannels);
    EEGplotcleaned.nbchan = length(iChannels);
    
    % PSDs (via pwelch)
    PSD_raw = pwelch(EEGplot.data', EEG.srate, EEG.srate/2, 1:EEG.srate/2, EEG.srate);
    PSD_cleaned = pwelch(EEGplotcleaned.data', EEG.srate, EEG.srate/2, 1:EEG.srate/2, EEG.srate);
    
    % Get psds averages over channels
    PSD_raw_mean = mean(PSD_raw,2,'omitnan');
    PSD_cleaned_mean = mean(PSD_cleaned,2,'omitnan');

    % get standard error of the mean
    PSD_raw_sem = std(PSD_raw, [], 2,'omitnan')./sqrt(size(PSD_raw,2));
    PSD_cleaned_sem = std(PSD_cleaned, [], 2,'omitnan')./sqrt(size(PSD_cleaned,2));

    % Calculate t-score for a 95% confidence interval
    t_score = tinv(.95, length(PSD_raw_mean)-1); % 6 degrees of freedom for 95% confidence

    % Plotting PSDs (mean +- sem)
    h = figure();
    % Plot raw PSD
    plot(frqs, PSD_raw_mean,'LineWidth',2.5,'Color', [1, .6, 0], 'DisplayName','Raw - MEAN PSD');
    hold on; grid on;
    upper_bound_raw = PSD_raw_mean(:) + t_score * PSD_raw_sem(:);
    lower_bound_raw = PSD_raw_mean(:) - t_score * PSD_raw_sem(:);
    patch([frqs(:); flipud(frqs(:))], [upper_bound_raw; flipud(lower_bound_raw)], [1, .6, 0], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName', 'MEAN ± SEM');

    % Plot cleaned PSDs
    plot(frqs, PSD_cleaned_mean,'LineWidth',2.5,'Color', 'b', 'DisplayName','Cleaned - MEAN PSD');
    upper_bound_cleaned = PSD_cleaned_mean(:) + t_score * PSD_cleaned_sem(:);
    lower_bound_cleaned = PSD_cleaned_mean(:) - t_score * PSD_cleaned_sem(:);
    patch([frqs(:); flipud(frqs(:))], [upper_bound_cleaned; flipud(lower_bound_cleaned)], 'b', ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName', 'MEAN ± SEM');
    xlim([1 125]),xlabel('Frequency [Hz]'); ylabel('PSD [au]');
    
    set(gca,'FontSize',22,'YScale','log');
    legend('show','interpreter','none');
    title('Averaged PSDs comparison before (raw) and after (cleaned) cleaning','FontSize',24);
    set(h,'WindowState','fullscreen');
    
    currentDir = pwd;
    matchingDir = dir(fullfile(currentDir, '**', ['*' 'CleaningPipeline\data_cleaned\plots\' '*']));
    plotsDir = matchingDir.folder; 

    exportgraphics(h,[plotsDir,subj_name,'_PSDs_rawcleaned.png'],"Resolution",300)
    savefig([plotsDir,subj_name,'_PSDs_rawcleaned'])

    % Topoplots for each cleaning step implemented
    % Cleaning steps are stored in EEGsteps table
    prepro_steps = EEGsteps.Properties.VariableNames; % cell array containing preprocessing implemented steps
    
    sum_pow = [];
    for istep = 1:length(prepro_steps)
        if iscell(EEGsteps.(istep))
            eeg = EEGsteps.(istep){1,1}; 
        else
            eeg = EEGsteps.(istep);
        end
        psd(:, :, istep) = pwelch(eeg.data',eeg.srate,eeg.srate/2,1:eeg.srate,eeg.srate);
        for ib = 1:length(bands)
            frq_inds = frqs >= bands{ib}(1) & frqs <= bands{ib}(2);
            sum_pow(ib, :, istep) = squeeze(sum(psd(frq_inds, :, :),1,'omitmissing'));
        end
    end


    % Power in frequency bands
    sum_pow_raw = []; sum_pow_cleaned = [];
    for ib = 1:length(bands)
        frq_inds = frqs >= bands{ib}(1) & frqs <= bands{ib}(2);
        sum_pow_raw(ib, :) = squeeze(sum(PSD_raw(frq_inds, :),1,'omitmissing'));
        sum_pow_cleaned(ib, :) = squeeze(sum(PSD_cleaned(frq_inds, :),1,'omitmissing'));
    end

    figure;
    fig_idx = 1;
    for ib = 1:length(bands)
        subplot(4,4,fig_idx);
        topoplot(zscore(sum_power_raw(ib,:)),EEG.chanlocs);
        colormap(cm17);
        colorbar
        fig_idx = fig_idx + 1;
    end
    for ib = 1:length(bands)
        subplot(4,4,fig_idx);
        topoplot(zscore(sum_power_cleaned(ib,:)),EEG.chanlocs);
        colormap(cm17);
        colorbar
        fig_idx = fig_idx + 1;
    end

end

if DefParam.save_steps == 0
    EEGsteps = [];
end
% removeCleaningPipelinePaths();

end
