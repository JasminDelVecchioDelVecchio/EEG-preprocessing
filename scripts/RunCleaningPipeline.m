function [FullEEGCleaned, EEGsteps] = RunCleaningPipeline(EEG, ECG, varargin)
% preprocess the input EEG

%   [EEG, varargout] = preprocess(data, params)
% Inputs:
%   - data is an EEGLAB data structure with raw EEG data
%   - params is an optional parameter specifying preprocessing parameters.
% Outputs:
%   - EEG is the output EEG after being preprocessed with additional fields.
%   - EEGsteps is a table (Ncol = N preprocessing step, each with EEGLAB data structs).
%
%   And example of the params where filtering parameters and clean_rawdata()
%   parameters are specified and other parameters are ommitted so that the
%   deafult value is used, is as shown below:
%   params = struct('FilterParams', struct('notch', struct('freq', 50), ...
%                                        'high',  struct('freq', 1,...
%                                                        'order', []),...
%                                        'low',   struct('freq', 80,...
%                                                         'order',
%                                                         []))),...
%                   'CRDParams',    struct('ChannelCriterion',   0.85,...
%                                          'LineNoiseCriterion', 4,...
%                                          'BurstCriterion',     5,...
%                                          'WindowCriterion',    0.25, ...
%                                          'Highpass',           [0.25 0.75])
%                    )
%
%   params must be a structure with optional fields
%   'FilterParams', 'NTParams','CRDParams',  'MARAParams', 'PrepParams',
%   'InterpolationParams',
%   to specify parameters for filtering, clean_rawdata(),
%   mara ica, prep robust average referencing, interpolation,
%   respectively.
%
%   If params is ommited, default values are used. If any of the fields
%   of params are ommited, corresponsing default values are used. If a
%   structure is given as 'struct([])' then the corresponding operation is
%   omitted and is not performed; for example, MARAParams = struct([])
%   skips the ICA and does not perform any ICA. Whereas if MARAParams =
%   struct() or if MARAParams is simply not given, then the default value
%   will be used.
%
%   'PrepParams' is an optional struture required by prep library. For more
%   information please see their documentation.
%
%   'InterpolationParams' is an optional structure with an optional field
%   'method' which can be on of the following chars: 'spherical',
%   'invdist' and 'spacetime'. To learn more about these
%   three methods please see eeg_interp.m of EEGLAB.
%
%   To learn more about 'FilterParams', 'CRDParams', 'MARAParams',
%   please see their corresponding functions performFilter.m,  performCleanrawdata.m,
%   performMaraICA.m.
%
%   A complete example of all possible parameters can be found in
%   RecommendedParameters.m (please do not change this file). All default
%   values are taken from DefaultParameters.m.
%
% Copyright (C) 2017  Amirreza Bahreini, methlabuzh@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Contributions: Andrea Canessa, Jasmin Del Vecchio Del Vecchio, Franziska
% Pellegrini, Stefan Haufe, Chiara Palmisano


    % Get Default Parameters
    Defaults = DefaultParameters;
    p = inputParser;
    % addParameter(p,'ResampleParams', Defaults.ResampleParams, @isstruct);
    addParameter(p,'FilterParams', Defaults.FilterParams, @isstruct);
    addParameter(p,'DetrendParams', Defaults.DetrendParams, @isstruct);
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
    DetrendParams = p.Results.DetrendParams;
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

% ---------------------------------------------------------------------- %
% Remove biosemi before preprocessing steps to avoid potential conflicts
allPaths = path;
allPaths = strsplit(allPaths, pathsep);
idx = contains(allPaths, 'biosig');
allPaths(~idx) = [];
automagicPaths = strjoin(allPaths, pathsep);
rmpath(automagicPaths);

% ---------------------------------------------------------------------- %
% Start Preprocessing
EEG.remove_range = create_bad_range(EEG.bad_segment);

[~,filename] = fileparts(EEG.filename);

if ~isfield(EEG,'pipeline')
    EEG.pipeline.suffix = '_ASR';
    EEG.pipeline.detrending.performed = 'no';
    EEG.pipeline.filtering.performed = 'no';
    EEG.pipeline.filtering.highpass.performed = 'no';
    EEG.pipeline.filtering.lowpass.performed = 'no';
    % EEG.pipeline.filtering.notch.performed = 'no';
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
    suffix = '_ASR';
    EEG.pipeline.suffix = suffix;
end

EEG.pipeline.params = params;

if ~isempty(EEG.icaweights)
    recomputeICA = false;
else
    recomputeICA = true;
end

% ---------------------------------------------------------------------- %

%set the random seed to a fix value equal to 1
% 
rng(1);

if exist('gcp','file')
    running_pool = gcp('nocreate');
    if isempty(running_pool)
        pool = parpool('local');
    end
else
    pool = [];
end

% EEGsteps will be used to store preprocessed EEG as EEGLAB structs after
% each cleaning step
EEGsteps = table();
% Convert the struct to a table
eegTable = table(EEG,VariableNames={'Raw'});
% Append the table as a new column to the main EEGsteps table
EEGsteps = [EEGsteps, eegTable];
clear eegTable;
% ---------------------------------------------------------------------- %
% Detect flat portion and loss of connection and mark it as bad segment
if ~any(EEG.bad_segment)
        mask = abs(diff(EEG.data,1,2))<(10e-8); % original value for thresholding: 10e-14
        mask = [mask(:,1) mask];
        EEG.bad_segment_flat = sum(mask)>0.9*EEG.nbchan;
        EEG.bad_segment_flat = conv(double(EEG.bad_segment_flat),ones(1,EEG.srate),'same')>0;
        EEG.remove_range = create_bad_range(EEG.bad_segment_flat);
end

% ---------------------------------------------------------------------- %
%Detrending
if ~isempty(DetrendParams) && strcmp(EEG.pipeline.detrending.performed,'no')
    EEG = performDetrend(EEG, DetrendParams);
    suffix = [suffix 'd'];
end

eegTable = table({EEG},VariableNames={'Detrend'});
% Append the table as a new column to the main EEGsteps table
EEGsteps = [EEGsteps, eegTable];
clear eegTable;
% ---------------------------------------------------------------------- %
% Detect bad channels based on PSD energy in [1 45] Hz range
if ~any(EEG.badchan)
    [~,EEG2] = evalc('pop_select(EEG,''nopoint'',EEG.remove_range)');

    plt = PSDchannel_selector(EEG2);
    clear EEG2;
    if ~autoICA % user-based pipeline
        waitfor(gcf)
    end
    rmchan = [plt.rmchan];
    EEG.badchan = EEG.badchan | rmchan;

    EEG.badchan_fixed = EEG.badchan;

end

% ---------------------------------------------------------------------- %
% Detect flat portion and loss of connection and mark it as bad segment
if ~any(EEG.bad_segment)    

    EEG.bad_segment = EEG.bad_segment | EEG.bad_segment_flat;
    % Detect high amplitude segment and mark it as bad
    dataBp = wlb_bandstop_filter(EEG.data(~EEG.badchan,~EEG.bad_segment)', EEG.srate,[100,120],4, 4); %edited version

    dataBp = abs(hilbert(wlb_bandpass_filter(dataBp, EEG.srate,[60, min(200,EEG.srate/2-10)], 4,1))');
    meandataBp = median(dataBp);
    meandataBp = wlb_detrend(meandataBp,round(5*EEG.srate),.5*round(5*EEG.srate),3,[],5);
    mad =(median(abs(meandataBp-median(meandataBp))));
    mask = meandataBp > (10*1.4826 * mad); % original value set to 10*...
    mask = conv(double(mask),ones(1,round(.5*EEG.srate)),'same')>0;
    EEG.bad_segment(~EEG.bad_segment) = mask;
    EEG.remove_range = create_bad_range(EEG.bad_segment);
end

eegTable = table({EEG},VariableNames={'ThrFreqAmp'});
% Append the table as a new column to the main EEGsteps table
EEGsteps = [EEGsteps, eegTable];
clear eegTable;

% ---------------------------------------------------------------------- %
% Remove line-noise with Zapline algorithm (Cheveigné et al., 2020)
if strcmp(EEG.pipeline.filtering.performed,'no') || strcmp(EEG.pipeline.Zapline.performed,'no')
    EEG = performZapline(EEG,FilterParams.notch);
    suffix = [suffix 'z'];
end
EEG.filename = [filename suffix '.set'];
EEG.pipeline.suffix = suffix;
eegTable = table({EEG},'VariableNames',{'Zapline'});
% Append the table as a new column to the main EEGsteps table
EEGsteps = [EEGsteps, eegTable];
clear eegTable;

% ---------------------------------------------------------------------- %
% Robust Rereference
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

if ~isempty(EOGParams) && strcmp(EEG.pipeline.EOGregression.performed,'no')
    % ---------------------------------------------------------------------- %
    % Two options are available for EOG cleaning: 
    % 1) ICA (ICLabel-based)
    % 2) Linear regression
    if perform_ICA 
        % 1)
        EEGcleaned = performICA(EEGcleaned, EOGParams, [], autoICA);
        suffix = [suffix 'e'];
        recomputeICA = true;
        
        % ---------------------------------------------------------------------- %
    else
        % 2)
        fprintf('Regressing out EOG...\n');
        % Convert the comma-separated string to a cell array of channel labels
        eog_channels_cell = strsplit(EOGParams.channels, ',');

        % Remove leading and trailing whitespaces in each channel label
        eog_channels_cell = strtrim(eog_channels_cell);

        % Convert all channel labels to lowercase for case-insensitive comparison
        eog_channels_cell = lower(eog_channels_cell);

        % Find the intersection of lowercase channel labels
        [~, eogchannels] = intersect(lower({EEG.chanlocs.labels}'), eog_channels_cell);
        eogs = EEGcleaned.data(eogchannels,:)';
        data = EEGcleaned.data(:,:)';

        % regress out EOG activity
        kk_eog = eogs\data;
        eogs = EEGcleaned.data(eogchannels,:)';
        data = EEGcleaned.data - kk_eog'*eogs';
        EEGcleaned.data = data;

        % We take out the two empty channels (frontal ones) via amplitude
        % thresholding
        kIQDp =3;
        logpower = log(std(data,[],2));
        Q=prctile(log(std(data,[],2)),[25 50 75]);
        EEGcleaned.badchan = EEG.badchan | (logpower'-Q(2)>kIQDp*(Q(3)-Q(1)));
        EEGcleaned.pipeline.EOGregression.performed = 'yes';
        EEGcleaned.pipeline.EOGregression.channels = EOGParams.channels;
        suffix = [suffix 'e'];
        recomputeICA = true;
    end
    eegTable = table({EEGcleaned},VariableNames={'EOGregression'});
    % Append the table as a new column to the main EEGsteps table
    EEGsteps = [EEGsteps, eegTable];
    clear eegTable;
end

% -------------------------------------------------------------------------------------------------------------%
% High-Pass Filtering
if strcmp(EEG.pipeline.filtering.performed,'no') || strcmp(EEG.pipeline.filtering.highpass.performed,'no')
    EEGcleaned = performFilter(EEGcleaned, 'high',FilterParams.high);
    suffix = [suffix 'h'];
    EEGcleaned.filename = [filename suffix '.set'];
    EEGcleaned.pipeline.suffix = suffix;
end

eegTable = table({EEGcleaned}, 'VariableNames',{'HP'});
% Append the table as a new column to the main EEGsteps table
EEGsteps = [EEGsteps, eegTable];
clear eegTable;

% ---------------------------------------------------------------------- %
% Clean EEG using clean_rawdata()
if ~isempty(CRDParams) && strcmp(EEG.pipeline.CleanRawData.performed,'no')

    [EEGcleaned] = performCleanrawdata(EEGcleaned, CRDParams);

    suffix = [suffix 'cr'];
    recomputeICA = true;

    % Detect high amplitude segment and mark it as bad
    dataBp = abs(hilbert(wlb_bandpass_filter(EEGcleaned.data(~EEGcleaned.badchan,:)', EEGcleaned.srate,[60, 120], 4,1))'); %edited version compatible with percept sample frequency
    meandataBp = mean(dataBp);
    meandataBp = wlb_detrend(meandataBp,round(5*EEGcleaned.srate),.5*round(5*EEGcleaned.srate),3,[],5);
    mad =(median(abs(meandataBp-median(meandataBp))));
    mask = meandataBp > (10*1.4826 * mad);
    mask = conv(double(mask),ones(1,round(.5*EEGcleaned.srate)),'same')>0;
    EEGcleaned.bad_segment = EEGcleaned.bad_segment | mask;
    EEGcleaned.remove_range = create_bad_range(EEGcleaned.bad_segment);


    eegTable = table({EEGcleaned},VariableNames={'CleanRawData'});
    % Append the table as a new column to the main EEGsteps table
    EEGsteps = [EEGsteps, eegTable];
    clear eegTable;
end
EEGcleaned.pipeline.suffix = suffix;

FullEEGCleaned = EEGcleaned;
FullEEGCleaned.filename = [filename suffix '.set'];

% ---------------------------------------------------------------------- %
% Epoching dataset
EpochDuration = EpochingParams.EpochDuration;
EventNames = EpochingParams.EventNames;
EvetnWdw = EpochingParams.EvetnWdw;
Baseline = EpochingParams.Baseline;


if ~isempty(EpochDuration)
    if EpochDuration > 0
        [EEGcleaned,original_timepoint] = mypop_epoch(EEGcleaned,'EpochDuration',EpochDuration);
        EEGcleaned.pipeline.Epoching.performed = 'yes';
        EEGcleaned.pipeline.Epoching.params = EpochingParams;
    end
else
    [EEGcleaned,original_timepoint] = mypop_epoch(EEGcleaned,'EventNames',EventNames,'EvetnWdw',EvetnWdw);
    EEGcleaned.pipeline.Epoching.performed = 'yes';
    EEGcleaned.pipeline.Epoching.params = EpochingParams;

    if(not(isempty(Baseline)))
        EEGcleaned = pop_rmbase( EEGcleaned, Baseline);
        EEGcleaned.pipeline.Epoching.baseline.performed = 'yes';
        EEGcleaned.pipeline.Epoching.baseline.window = Baseline;
    end
end
% -------------------------------------------------------------------------------------------------------------%
    % Review Data Cleaning by eye
    if ~autoICA
    plt = fastplot2(EEGcleaned);
        waitfor(gcf)
        ButtonName = questdlg('Are you satisfied with cleaning results', 'Verify channels removal', 'Yes', 'Cancel', 'Yes'); %ButtonName = 'Yes';%

        switch ButtonName
            case 'Yes'
                EEG = EEGcleaned;
                clear EEGcleaned
            case 'Cancel'
                EEG = [];
                return;
        end

        badChans = (plt.badchan{1});
        bad_segments = plt.bad_segment{1};
        bad_epoch = plt.bad_epoch{1};
        clear plt;

        if any(EEG.badchan~=badChans)||any(EEG.bad_segment~=bad_segments)||any(EEG.bad_epoch~=bad_epoch)
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

        % remove_range = [];
        % if any(FullEEGCleaned.bad_segment)
        %     firsts = find(diff(FullEEGCleaned.bad_segment) == 1) + 1;
        %     seconds = find(diff(FullEEGCleaned.bad_segment) == -1);
        %     if(firsts(1) > seconds(1))
        %         firsts = [1, firsts];
        %     end
        %     if(seconds(end) < firsts(end))
        %         seconds = [seconds, length(FullEEGCleaned.bad_segment)];
        %     end
        %     remove_range = [firsts;seconds]';
        % end
        FullEEGCleaned.remove_range = create_bad_range(FullEEGCleaned.bad_segment);
        % FullEEGCleaned.remove_range = remove_range;
        eegTable = table({FullEEGCleaned},VariableNames={'VisualInspection'});
        % Append the table as a new column to the main EEGsteps table
        EEGsteps = [EEGsteps, eegTable];
        clear eegTable;
    else
        FullEEGCleaned = EEGcleaned;
        clear EEGCleaned
    end
% -------------------------------------------------------------------------------------------------------------%
% Compute ICA
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
        % ICLabel is designed to classify ICs fitted with an extended infomax ICA decomposition algorithm on EEG 
        % datasets referenced to a common average and filtered between [1., 100.] Hz. It is possible to run ICLabel 
        % on datasets that do not meet those specification, but the classification performance might be negatively 
        % impacted.
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
    eegTable = table({FullEEGCleaned},VariableNames={'ICA'});
    % Append the table as a new column to the main EEGsteps table
    EEGsteps = [EEGsteps, eegTable];
    clear eegTable;
else
    % ECG artifact regression
    ecgs = ECG.data;
    data = EEG.data(:,:)';
    if size(ecgs,2)>size(ecgs,1)
        ecgs = ecgs';
    end

    kk_ecg = ecgs.\data;
    data = EEG.data' - kk_ecg.*ecgs;

    % FullEEGCleaned.data = data';
    suffix = [suffix 'r'];
    FullEEGCleaned.filename = [filename suffix '.set'];
    FullEEGCleaned.pipeline.suffix = suffix;

    orig_icasphere=[];
    orig_icachansind=[];
    orig_icaweights=[];
    orig_icawinv=[];

    FullEEGCleaned.data = data';


    eegTable = table({FullEEGCleaned},VariableNames={'ECGregression'});
    % Append the table as a new column to the main EEGsteps table
    EEGsteps = [EEGsteps, eegTable];
    clear eegTable;
end


% Lowpass filtering
% FullEEGCleaned = performFilter(FullEEGCleaned, 'low', FilterParams.low);
[a,b] = butter(5,params.FilterParams.low.freq(1)./FullEEGCleaned.srate);
FullEEGCleaned.data = filtfilt(a,b,FullEEGCleaned.data);
suffix = [suffix 'l'];
FullEEGCleaned.filename = [filename suffix '.set'];
FullEEGCleaned.pipeline.suffix = suffix;


eegTable = table({FullEEGCleaned},VariableNames={'LP'});
% Append the table as a new column to the main EEGsteps table
EEGsteps = [EEGsteps, eegTable];
clear eegTable;

badchan = find(FullEEGCleaned.badchan);

% Bad channel interpolation
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

        FullEEGCleaned.badchan = false(1,size(FullEEGCleaned.data,1));
        FullEEGCleaned.interpolatedchan(badchan) = true;
    end

    % Restore the seed generator to the previous state
    rng('default');

    % Plotting cleaning results 
    % 1) Topoplots for each cleaning step and for each frequency band of
    % interest

end