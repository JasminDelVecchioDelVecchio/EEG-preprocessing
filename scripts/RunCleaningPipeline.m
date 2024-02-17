function [FullEEGCleaned, EEGsteps] = RunCleaningPipeline(EEG, ECG, varargin)
% preprocess the input EEG

%   [EEG, varargout] = preprocess(EEG, ECG)
% Inputs:
%   - EEG is an EEGLAB data structure with raw EEG data
%   - ECG is an optional parameter representing synchronized cardiac
%   activity. It is used only when the regression mode is specified in the
%   DefaultParameters.m file.
% Outputs:
%   - FullEEGCleaned is the output EEG after being preprocessed with additional fields.
%   - EEGsteps is a table containing EEG signals after each preprocessing step.
%
%   If ECG is ommited, ICA will be performed by default


% Contributions: Jasmin Del Vecchio Del Vecchio, Andrea Canessa, Chiara Palmisano, Franziska
% Pellegrini, Stefan Haufe


%% Get Default Parameters
Defaults = DefaultParameters;
p = inputParser;
addParameter(p,'FilterParams', Defaults.FilterParams, @isstruct);
addParameter(p,'PrepParams', Defaults.PrepParams, @isstruct);
addParameter(p,'CRDParams', Defaults.CRDParams, @isstruct);
addParameter(p,'FastICAParams', Defaults.FastICAParams, @isstruct);
addOptional(p,'EpochingParams', Defaults.EpochingParams, @isstruct);
addOptional(p,'EOGParams', Defaults.EOGParams, @isstruct);
addOptional(p,'ECGParams', Defaults.ECGParams, @isstruct);
addParameter(p,'prunedata', Defaults.prunedata, @islogical);
addParameter(p,'autoICA', Defaults.autoICA, @islogical);
addOptional(p,'perform_ICA', Defaults.perform_ICA, @islogical);
addOptional(p,'subcompICA', Defaults.subcompICA, @islogical);
addOptional(p,'interpbadchannels', Defaults.interpbadchannels, @islogical);
addOptional(p,'save_steps', Defaults.save_steps, @islogical);
addOptional(p,'save_fig', Defaults.interpbadchannels, @islogical);

parse(p, varargin{:});
params = p.Results;
FilterParams = p.Results.FilterParams;
CRDParams = p.Results.CRDParams;
PrepParams = p.Results.PrepParams;
FastICAParams = p.Results.FastICAParams;
EpochingParams = p.Results.EpochingParams;
EOGParams = p.Results.EOGParams;
ECGParams = p.Results.ECGParams;
prunedata = p.Results.prunedata;
autoICA = p.Results.autoICA;
perform_ICA = p.Results.perform_ICA;
subcompICA = p.Results.subcompICA;
interpbadchannels = p.Results.interpbadchannels;
save_steps = p.Results.save_steps;
save_fig = p.Results.save_fig;

clear p varargin;

if ~isfield(EEG,'badchan')
    EEG.badchan = false(1,EEG.nbchan);
end
if ~isfield(EEG,'bad_segment')
    EEG.bad_segment = false(1,EEG.pnts*EEG.trials);
    EEG.remove_range = [];
end
if ~isfield(EEG,'bad_epoch')
    if EEG.trials==1
        EEG.bad_epoch = [];
    else
        EEG.bad_epoch = sum(reshape(EEG.bad_segment,EEG.pnts,EEG.trials))>0;
    end
end
if ~isfield(EEG,'interpolatedchan')
    EEG.interpolatedchan = false(1,EEG.nbchan);
end

%% Remove biosemi before preprocessing steps to avoid potential conflicts
allPaths = path;
allPaths = strsplit(allPaths, pathsep);
idx = contains(allPaths, 'biosig');
allPaths(~idx) = [];
automagicPaths = strjoin(allPaths, pathsep);
rmpath(automagicPaths);

%% Start Preprocessing
EEG.remove_range = create_bad_range(EEG.bad_segment);

[~,filename] = fileparts(EEG.filename);

if ~isfield(EEG,'pipeline')
    EEG.pipeline.suffix = '_cleaned_';
    EEG.pipeline.detrending.performed = 'no';
    EEG.pipeline.filtering.performed = 'no';
    EEG.pipeline.filtering.highpass.performed = 'no';
    EEG.pipeline.filtering.lowpass.performed = 'no';
    EEG.pipeline.Zapline.performed = 'no';
    EEG.pipeline.EOGregression.performed = 'no';
    EEG.pipeline.RobustReference.performed = 'no';
    EEG.pipeline.CleanRawData.performed = 'no';
    EEG.pipeline.Epoching.performed = 'no';
    EEG.pipeline.FastICA.performed = 'no';
    EEG.pipeline.ICrej.performed = 'no';
    EEG.pipeline.interpolation.performed = 'no';
    EEG.pipeline.ECGremoval.performed = 'no';
end

if isfield(EEG.pipeline,'suffix')
    suffix = EEG.pipeline.suffix;
else
    suffix = '_cleaned_';
    EEG.pipeline.suffix = suffix;
end

EEG.pipeline.params = params;

if ~isempty(EEG.icaweights)
    recomputeICA = false;
else
    recomputeICA = true;
end


%% set the random seed to a fix value equal to 1
rng(1);

%% Speed up computation with Parallel Computing toolbox
if exist('gcp','file')
    running_pool = gcp('nocreate');
    if isempty(running_pool)
        pool = parpool('local');
    end
else
    pool = [];
end

%% EEGsteps initialization
% It will be used to store preprocessed EEG data after each cleaning step
EEGsteps = table();

% Convert the struct to a table
eegTable = table(EEG,VariableNames={'Raw'});
% Append the table as a new column to the main EEGsteps table
EEGsteps = [EEGsteps, eegTable];
clear eegTable;

%% Detect bad channels based on PSD energy in [1 45] Hz range
if ~any(EEG.badchan)
    [~,EEG2] = evalc('pop_select(EEG,''nopoint'',EEG.remove_range)');

    plt = PSDchannel_selector(EEG2);
    close all;
    clear EEG2;
    % if ~autoICA % user-based pipeline
        % waitfor(gcf)
    % end
    rmchan = [plt.rmchan];
    EEG.badchan = EEG.badchan | rmchan;

    EEG.badchan_fixed = EEG.badchan;

end

%% Detect flat portion and loss of connection and mark it as bad segment
if ~any(EEG.bad_segment)
    mask = abs(diff(EEG.data,1,2))<(10e-8);
    mask = [mask(:,1) mask];
    EEG.bad_segment_flat = sum(mask) > 0.9 * EEG.nbchan;
    EEG.bad_segment_flat = conv(double(EEG.bad_segment_flat),ones(1,EEG.srate),'same')>0;
    EEG.remove_range = create_bad_range(EEG.bad_segment_flat);

    EEG.bad_segment = EEG.bad_segment | EEG.bad_segment_flat;
    % Detect high amplitude segment and mark it as bad
    dataBp = wlb_bandstop_filter(EEG.data(~EEG.badchan, ~EEG.bad_segment)', EEG.srate, [60, EEG.srate/2 - 5], 4, 4);
    dataBp = abs(hilbert(wlb_bandpass_filter(dataBp, EEG.srate, [60, EEG.srate/2-5], 4, 1))');
    meandataBp = median(dataBp);
    meandataBp = wlb_detrend(meandataBp, round(5*EEG.srate), 0.5 * round(5 * EEG.srate), 3, [], 5);
    mad = (median(abs(meandataBp - median(meandataBp))));
    mask = meandataBp > (10 * 1.4826 * mad);
    mask = conv(double(mask), ones(1, round(0.5*EEG.srate)),'same')>0;
    EEG.bad_segment(~EEG.bad_segment) = mask;
    EEG.remove_range = create_bad_range(EEG.bad_segment);
end

%% Remove line-noise with Zapline algorithm (Cheveigné et al., 2020)
if strcmp(EEG.pipeline.filtering.performed,'no') || strcmp(EEG.pipeline.Zapline.performed,'no')
    EEG = performZapline(EEG, FilterParams.notch);
    suffix = [suffix 'z'];
end
EEG.filename = [filename suffix '.set'];
EEG.pipeline.suffix = suffix;
eegTable = table({EEG},'VariableNames',{'Zapline'});
% Append the table as a new column to the main EEGsteps table
EEGsteps = [EEGsteps, eegTable];
clear eegTable;

%% Robust Rereference
if ~isempty(PrepParams) && strcmp(EEG.pipeline.RobustReference.performed,'no')
    [EEGcleaned, referenceOut] = performReference(EEG, PrepParams);
    %%%remove reference signal also from original EEG variable for successive comparison
    EEG.data = bsxfun(@minus, EEG.data, referenceOut.referenceSignal);
    suffix = [suffix 'p'];
    recomputeICA = true;
    eegTable = table({EEGcleaned}, 'VariableNames', {'RobustReference'});
    % Append the table as a new column to the main EEGsteps table
    EEGsteps = [EEGsteps, eegTable];
    clear eegTable;
else
    EEGcleaned = EEG;
end

EEG.pipeline.suffix = suffix;
EEG.filename = [filename suffix '.set'];

%% EOG regression
if ~isempty(EOGParams) && strcmp(EEG.pipeline.EOGregression.performed,'no')
    % ---------------------------------------------------------------------- %
    % Two options will be available for EOG cleaning: 
    % 1) SVD + ICLabel classifier
    % 2) Linear regression
    if perform_ICA 
        % 1)
        EEGcleaned = performEOGregression(EEGcleaned, EOGParams);
        suffix = [suffix 'e'];
        recomputeICA = true;
    end
        % ---------------------------------------------------------------------- %
    % else
    %     % 2)
    %     fprintf('Regressing out EOG...\n');
    %     % Convert the comma-separated string to a cell array of channel labels
    %     eog_channels_cell = strsplit(EOGParams.channels, ',');
    % 
    %     % Remove leading and trailing whitespaces in each channel label
    %     eog_channels_cell = strtrim(eog_channels_cell);
    % 
    %     % Convert all channel labels to lowercase for case-insensitive comparison
    %     eog_channels_cell = lower(eog_channels_cell);
    % 
    %     % Find the intersection of lowercase channel labels
    %     [~, eogchannels] = intersect(lower({EEG.chanlocs.labels}'), eog_channels_cell);
    %     eogs = EEGcleaned.data(eogchannels, :)';
    %     data = EEGcleaned.data(:, :)';
    % 
    %     % regress out EOG activity
    %     kk_eog = eogs\data;
    %     eogs = EEGcleaned.data(eogchannels, :)';
    %     data = EEGcleaned.data - kk_eog'*eogs';
    %     EEGcleaned.data = data;
    % 
    %     % We take out the two empty channels (frontal ones) via amplitude
    %     % thresholding
    %     kIQDp =3;
    %     logpower = log(std(data,[], 2));
    %     Q=prctile(log(std(data,[],2)),[25 50 75]);
    %     EEGcleaned.badchan = EEG.badchan | (logpower'-Q(2) > kIQDp * (Q(3) - Q(1)));
    %     EEGcleaned.pipeline.EOGregression.performed = 'yes';
    %     EEGcleaned.pipeline.EOGregression.channels = EOGParams.channels;
    %     suffix = [suffix 'e'];
    %     recomputeICA = true;
    % end
    eegTable = table({EEGcleaned}, VariableNames={'EOGregression'});
    % Append the table as a new column to the main EEGsteps table
    EEGsteps = [EEGsteps, eegTable];
    clear eegTable;
end

%% High-Pass Filtering
if strcmp(EEG.pipeline.filtering.performed, 'no') || strcmp(EEG.pipeline.filtering.highpass.performed, 'no')
    EEGcleaned = performFilter(EEGcleaned, 'high', FilterParams.high);
    suffix = [suffix 'h'];
    EEGcleaned.filename = [filename suffix '.set'];
    EEGcleaned.pipeline.suffix = suffix;
end

%% Clean EEG using clean_rawdata()
if ~isempty(CRDParams) && strcmp(EEG.pipeline.CleanRawData.performed, 'no')

    [EEGcleaned] = performCleanrawdata(EEGcleaned, CRDParams);

    suffix = [suffix 'cr'];
    recomputeICA = true;

    % Detect high amplitude segment and mark it as bad
    dataBp = abs(hilbert(wlb_bandpass_filter(EEGcleaned.data(~EEGcleaned.badchan,:)', EEGcleaned.srate, [60, (EEGcleaned.srate/2 - 5)], 4, 1))'); %edited version compatible with percept sample frequency
    meandataBp = mean(dataBp);
    meandataBp = wlb_detrend(meandataBp, round(5 * EEGcleaned.srate), 0.5*round(5 * EEGcleaned.srate), 3, [], 5);
    mad =(median(abs(meandataBp - median(meandataBp))));
    mask = meandataBp > (10 * 1.4826 * mad);
    mask = conv(double(mask), ones(1,round(0.5*EEGcleaned.srate)), 'same') > 0;
    EEGcleaned.bad_segment = EEGcleaned.bad_segment | mask;
    EEGcleaned.remove_range = create_bad_range(EEGcleaned.bad_segment);


    eegTable = table({EEGcleaned},VariableNames={'HP+CleanRawData'});
    % Append the table as a new column to the main EEGsteps table
    EEGsteps = [EEGsteps, eegTable];
    clear eegTable;
end
EEGcleaned.pipeline.suffix = suffix;

FullEEGCleaned = EEGcleaned;
FullEEGCleaned.filename = [filename suffix '.set'];

%% Review Data Cleaning by eye
% Epoching dataset (required by fastplot2)
EpochDuration = EpochingParams.EpochDuration;
EventNames = EpochingParams.EventNames;
EvetnWdw = EpochingParams.EvetnWdw;
Baseline = EpochingParams.Baseline;


if ~isempty(EpochDuration)
    if EpochDuration > 0
        [EEGcleaned, original_timepoint] = mypop_epoch(EEGcleaned, 'EpochDuration', EpochDuration);
        EEGcleaned.pipeline.Epoching.performed = 'yes';
        EEGcleaned.pipeline.Epoching.params = EpochingParams;
    end
else
    [EEGcleaned,original_timepoint] = mypop_epoch(EEGcleaned, 'EventNames', EventNames, 'EvetnWdw', EvetnWdw);
    EEGcleaned.pipeline.Epoching.performed = 'yes';
    EEGcleaned.pipeline.Epoching.params = EpochingParams;

    if(not(isempty(Baseline)))
        EEGcleaned = pop_rmbase(EEGcleaned, Baseline);
        EEGcleaned.pipeline.Epoching.baseline.performed = 'yes';
        EEGcleaned.pipeline.Epoching.baseline.window = Baseline;
    end
end

    if ~autoICA
    plt = fastplot2(EEGcleaned);
        waitfor(gcf)
        ButtonName = questdlg('Are you satisfied with cleaning results', 'Verify channels removal', 'Yes', 'Cancel', 'Yes'); %ButtonName = 'Yes';%

        switch ButtonName
            case 'Yes'
                EEG = EEGcleaned;
                clear EEGcleaned
            case 'Cancel'
                return;
        end

        badChans = (plt.badchan{1});
        bad_segments = plt.bad_segment{1};
        bad_epoch = plt.bad_epoch{1};
        clear plt;

        if any(EEG.badchan~=badChans) || any(EEG.bad_segment~=bad_segments) || any(EEG.bad_epoch~=bad_epoch)
            recomputeICA = true;
        end

        % Update badchan and badsegment struct field in EEG struct
        EEG.badchan = badChans;
        EEG.bad_segment = bad_segments;
        EEG.bad_epoch = bad_epoch;

        % % Update badchan and badsegment struct field in pre-epoching FullEEGCleaned dataset
        FullEEGCleaned.badchan = EEG.badchan;
        FullEEGCleaned.bad_segment(original_timepoint) = EEG.bad_segment;
        FullEEGCleaned.bad_epoch = EEG.bad_epoch;

         FullEEGCleaned.remove_range = create_bad_range(FullEEGCleaned.bad_segment);
    else
        FullEEGCleaned = EEGcleaned;
        clear EEGCleaned
    end

%% Perform ICA
if perform_ICA
    if strcmp(EEG.pipeline.FastICA.performed,'no') || recomputeICA
        EEG.icaweights = [];
        EEG.icawinv = [];
        EEG.icasphere = [];
        EEG.icachansind=[];
        EEG.reject.gcompreject = [];
        EEG.icaact = [];
        epochparamstruct.EpochingParams = EpochingParams;
        EEG = performFastICA(EEG, FastICAParams);
        suffix = [suffix 'I'];
    end
    
    %Update ICA struct field in FullEEGCleaned struct
    FullEEGCleaned.icaweights = EEG.icaweights;
    FullEEGCleaned.icawinv = EEG.icawinv;
    FullEEGCleaned.icasphere = EEG.icasphere;
    FullEEGCleaned.icachansind=EEG.icachansind;
    FullEEGCleaned.pipeline.FastICA = EEG.pipeline.FastICA;
    FullEEGCleaned.filename = [filename suffix '.set'];
    FullEEGCleaned.pipeline.suffix = suffix;

    %ICLabel IC rejection
    if ~isempty(EEG.icaweights)
        % https://mne.tools/mne-icalabel/stable/generated/api/mne_icalabel.iclabel.iclabel_label_components.html
        epochparamstruct.EpochingParams = EpochingParams;
        [EEG] = performIClabel(FullEEGCleaned,epochparamstruct);

        %Visual IC rejection
        [EEG] = performICrejection(EEG);
    end

    % close all
    FullEEGCleaned.reject.gcompreject = EEG.reject.gcompreject;
    FullEEGCleaned.pipeline.ICrej.rejectComponents = EEG.pipeline.ICrej.rejectComponents;

%clear EEG

if prunedata

    if subcompICA

        artcomps = FullEEGCleaned.pipeline.ICrej.rejectComponents;

        % Reconstruct signal based on IC rejection
        if ~isempty(artcomps)
            FullEEGCleaned = pop_subcomp(FullEEGCleaned, artcomps);
        end
        suffix = [suffix 'j'];
        FullEEGCleaned.filename = [filename suffix '.set'];
        FullEEGCleaned.pipeline.suffix = suffix;
    end

    orig_icasphere=FullEEGCleaned.icasphere;
    orig_icachansind=FullEEGCleaned.icachansind;
    orig_icaweights=FullEEGCleaned.icaweights;
    orig_icawinv=FullEEGCleaned.icawinv;
end
    eegTable = table({FullEEGCleaned}, VariableNames={'ICA'});
    % Append the table as a new column to the main EEGsteps table
    EEGsteps = [EEGsteps, eegTable];
    clear eegTable;
end
% else
%     % ECG artifact regression
%     ecgs = ECG.data;
%     data = EEG.data(:,:)';
%     if size(ecgs,2) > size(ecgs,1)
%         ecgs = ecgs';
%     end
% 
%     kk_ecg = ecgs.\data;
%     data = EEG.data' - kk_ecg.*ecgs;
% 
%     % FullEEGCleaned.data = data';
%     suffix = [suffix 'r'];
%     FullEEGCleaned.filename = [filename suffix '.set'];
%     FullEEGCleaned.pipeline.suffix = suffix;
% 
%     orig_icasphere=[];
%     orig_icachansind=[];
%     orig_icaweights=[];
%     orig_icawinv=[];
% 
%     FullEEGCleaned.data = data';
% 
% 
%     eegTable = table({FullEEGCleaned}, VariableNames={'ECGregression'});
%     % Append the table as a new column to the main EEGsteps table
%     EEGsteps = [EEGsteps, eegTable];
%     clear eegTable;
% end


%% Lowpass filtering
% [a,b] = butter(5,params.FilterParams.low.freq(1)./FullEEGCleaned.srate);
% FullEEGCleaned.data = filtfilt(a,b,FullEEGCleaned.data);
% suffix = [suffix 'l'];
% FullEEGCleaned.filename = [filename suffix '.set'];
% FullEEGCleaned.pipeline.suffix = suffix;
FullEEGCleaned = performFilter(FullEEGCleaned, 'low', FilterParams.low);
suffix = [suffix 'l'];
FullEEGCleaned.filename = [filename suffix '.set'];
FullEEGCleaned.pipeline.suffix = suffix;
badchan = find(FullEEGCleaned.badchan);

%% Bad channel interpolation
    if interpbadchannels

        FullEEGCleaned = eeg_interp(FullEEGCleaned , sort(badchan) , 'spherical');
        FullEEGCleaned.pipeline.interpolation.performed = 'yes';
        suffix = [suffix 'i'];
        FullEEGCleaned.filename = [filename suffix '.set'];
        FullEEGCleaned.pipeline.suffix = suffix;

        % Put the original icadata back into the structure
        FullEEGCleaned.icasphere=orig_icasphere;
        FullEEGCleaned.icachansind=orig_icachansind;
        FullEEGCleaned.icaweights= orig_icaweights;
        FullEEGCleaned.icawinv=orig_icawinv;

        FullEEGCleaned.badchan = false(1, size(FullEEGCleaned.data,1));
        FullEEGCleaned.interpolatedchan(badchan) = true;
    end

    % Restore the seed generator to the previous state
    rng('default');
    eegTable = table({FullEEGCleaned}, VariableNames={'LP+Interp'});
    % Append the table as a new column to the main EEGsteps table
    EEGsteps = [EEGsteps, eegTable];
    clear eegTable;
end