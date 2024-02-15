function [EEG,referenceOut] = performReference(EEG,params)
% Inputs:
%   - EEG is an EEGLAB data structure with EEG data

toRemove = find(EEG.badchan);
remove_range = [];%EEG.remove_range;
[~,EEGCleaned] = evalc('pop_select(EEG,''nochannel'',toRemove)');
oldremovedMask = EEG.badchan;
removedMask = EEG.badchan;
fprintf('Performing robust rereferencing...\n');
% [EEGCleaned,referenceOut] = Reference(EEG, params, remove_range);
[EEGCleaned,referenceOut] = Reference(EEGCleaned, params, remove_range);

if(isfield(referenceOut, 'noisyStatistics'))
    noisechan = referenceOut.noisyStatistics;
    if(isfield(noisechan, 'noisyChannels'))
        newmask = false(1,EEGCleaned.nbchan);
        newmask(noisechan.noisyChannels.all) = true;
        % removedMask(~removedMask) = newmask; 
        removedMask(~removedMask) = newmask(~removedMask); 
    end
end

EEG.badchan = removedMask;
EEG.data(~oldremovedMask,:) = EEGCleaned.data(~oldremovedMask,:);
% EEG.data(~oldremovedMask,:) = EEGCleaned.data;

clear EEGCleaned;

EEG.pipeline.RobustReference.performed = 'yes';
EEG.pipeline.RobustReference.params = params;
end

function [signal, referenceOut] =  Reference(signal, referenceIn, remove_range)
% Perform the specified reference
%
% Parameters:
%    signal       EEGLAB EEG or similar structure
%                 (must have .srate, .data, and .chanlocs fields)
%    referenceIn  (optional) Input reference structure --- or empty
%    signal       (output) Input signal after referencing and interpolation
%
% Default for the PREP pipeline is:
%    referenceType = 'robust'
%    interpolationOrder = 'post'
%    meanEstimateType = 'median'
%
% Other possible reference types: average, specific, and none
%
%% Check the input parameters
if nargin < 1
    error('performReference:NotEnoughArguments', 'requires at least 1 argument');
elseif ~isstruct(signal)
    error('performReference:SignalNotStructue', 'signal must be a structure');
elseif ~isfield(signal, 'data')
    error('performReference:NoDataField', 'signal structure requires a data field');
elseif ~isfield(signal, 'chanlocs') || ~isfield(signal.chanlocs, 'labels')
    error('performReference:NoChanLocs', 'signal structure requires a chanlocs field');
elseif size(signal.data, 3) ~= 1
    error('performReference:DataNotContinuous', 'signal.data must be a 2D array');
elseif size(signal.data, 2) < 2
    error('performReference:NoData', 'signal.data must have multiple points');
elseif ~exist('referenceIn', 'var') || isempty(referenceIn)
    referenceIn = struct();
end
if ~isstruct(referenceIn)
    error('performReference:NoData', 'second argument must be a structure')
end

%% Set the defaults and initialize as needed
referenceOut = getReferenceStructure();
defaults = getPrepDefaults(signal, 'reference');
[referenceOut, errors] = checkDefaults(referenceIn, referenceOut, defaults);
if ~isempty(errors)
    error('performReference:BadParameters', ['|' sprintf('%s|', errors{:})]);
end
defaults = getPrepDefaults(signal, 'detrend');
[referenceOut, errors] = checkDefaults(referenceIn, referenceOut, defaults);
if ~isempty(errors)
    error('performReference:BadParameters', ['|' sprintf('%s|', errors{:})]);
end
referenceOut.rereferencedChannels = sort(referenceOut.rereferencedChannels);
referenceOut.referenceChannels = sort(referenceOut.referenceChannels);
referenceOut.evaluationChannels = sort(referenceOut.evaluationChannels);

%% Calculate the reference for the original signal
if isempty(referenceOut.referenceChannels) || ...
        strcmpi(referenceOut.referenceType, 'none')
    referenceOut.referenceSignalOriginal = ...
        zeros(1, size(signal.data, 2));
else
    referenceOut.referenceSignalOriginal = mean(signal.data(referenceOut.referenceChannels, :), 1,'omitnan');
        %         nanmean(signal.data(referenceOut.referenceChannels, :), 1);
end
%% Make sure that reference channels have locations for interpolation
chanlocs = referenceOut.channelLocations(referenceOut.evaluationChannels);
if ~(length(cell2mat({chanlocs.X})) == length(chanlocs) && ...
        length(cell2mat({chanlocs.Y})) == length(chanlocs) && ...
        length(cell2mat({chanlocs.radius})) == length(chanlocs)) && ...
        ~(length(cell2mat({chanlocs.theta})) == length(chanlocs) && ...
        length(cell2mat({chanlocs.radius})) == length(chanlocs))
    error('performReference:NoChannelLocations', ...
        'evaluation channels must have locations');
end

%% Now perform the particular combinations
if  strcmpi(referenceOut.referenceType, 'robust') && ...
        strcmpi(referenceOut.interpolationOrder, 'post-reference')
    doRobustPost();
elseif  strcmpi(referenceOut.referenceType, 'robust') && ...
        strcmpi(referenceOut.interpolationOrder, 'pre-reference')
    doRobustPre();
elseif strcmpi(referenceOut.referenceType, 'average')
    doOther();
elseif strcmpi(referenceOut.referenceType, 'specific')
    doOther();
else
    doOther();
end


    function [] = doRobustPre()
        % Use the bad channels accumulated from reference search to robust
        sig2 = pop_select(signal,'nopoint',remove_range);
        
        referenceOut = robustReference(sig2, referenceOut);
        referenceOut.noisyStatisticsBeforeInterpolation = ...
            referenceOut.noisyStatistics;
        noisy = referenceOut.badChannels.all;
        if isempty(noisy)   %No noisy channels -- ordinary ref
            referenceOut.referenceSignal = ...
            mean(signal.data(referenceOut.referenceChannels, :), 1,'omitnan');
            %                 nanmean(signal.data(referenceOut.referenceChannels, :), 1);
        else
            bad = signal.data(noisy, :);
            sourceChannels = setdiff(referenceOut.evaluationChannels, noisy);
            signal = interpolateChannels(signal, noisy, sourceChannels);
            referenceOut.interpolatedChannelNumbers = referenceOut.interpolatedChannels.all;
            referenceOut.referenceSignal = ...
                mean(signal.data(referenceOut.referenceChannels, :), 1,'omitnan');
            %                 nanmean(signal.data(referenceOut.referenceChannels, :), 1);
            referenceOut.badSignalsUninterpolated = ...
                bad - repmat(referenceOut.referenceSignal, length(noisy), 1);
        end
        
        signal = removeReference(signal, referenceOut.referenceSignal, ...
            referenceOut.rereferencedChannels);
        
        sig2 = pop_select(signal,'nopoint',remove_range);
        
        referenceOut.noisyStatistics = ...
            findNoisyChannels(removeTrend(sig2, referenceOut), referenceOut);
    end

    function [] = doRobustPost()
        % Robust reference with interpolation afterwards
        [~,sig2] = evalc('pop_select(signal,''nopoint'',remove_range)');
        referenceOut = robustReference(sig2, referenceOut);
        noisy = referenceOut.badChannels.all;
        if isempty(noisy)   %No noisy channels -- ordinary ref
            referenceOut.referenceSignal = ...
                mean(signal.data(referenceOut.referenceChannels, :), 1,'omitnan');
        else
            sourceChannels = setdiff(referenceOut.evaluationChannels, noisy);
            [~,signalNew] = evalc('interpolateChannels(signal, noisy, sourceChannels)');
            referenceOut.referenceSignal = ...
                mean(signalNew.data(referenceOut.referenceChannels, :), 1,'omitnan');
            clear signalNew;
        end
        sig2 = removeReference(signal, referenceOut.referenceSignal, ...
            referenceOut.rereferencedChannels);
        
        [~,sig2] = evalc('pop_select(signal,''nopoint'',remove_range)');
        
        textbar = textprogressbar(sprintf('\tFinal iteration: '));
        
        referenceOut.noisyStatistics  = ...
            findNoisyChannels(removeTrend(sig2, referenceOut), referenceOut);
        
        referenceOut.noisyStatisticsBeforeInterpolation = ...
            referenceOut.noisyStatistics;
        %% Bring forward unusable channels from original data
        noisy = referenceOut.noisyStatisticsOriginal.noisyChannels;
        unusableChans = union(noisy.badChannelsFromNaNs, ...
            union(noisy.badChannelsFromNoData, ...
            noisy.badChannelsFromLowSNR));
        intChans = referenceOut.noisyStatistics.noisyChannels;
        chans = union(intChans.all, unusableChans);
        intChans.all = chans(:)';
        chans = union(intChans.badChannelsFromNaNs, noisy.badChannelsFromNaNs);
        intChans.badChannelsFromNaNs = chans(:)';
        chans = union(intChans.badChannelsFromNoData, noisy.badChannelsFromNoData);
        intChans.badChannelsFromNoData = chans(:)';
        chans = union(intChans.badChannelsFromLowSNR, noisy.badChannelsFromLowSNR);
        intChans.badChannelsFromLowSNR = chans(:)';
        referenceOut.badChannels = intChans;
        
        %% Now find the bad channels and interpolate
        %referenceOut.noisyStatisticsForReference = noisyStatistics;
        noisyChans = referenceOut.badChannels.all;
        if isempty(noisyChans)
            return;
        end
        bad = signal.data(noisyChans, :);
        sourceChannels = setdiff(referenceOut.evaluationChannels, noisyChans);
        [~,signalNew] = evalc('interpolateChannels(signal, noisyChans, sourceChannels)');
        referenceOut.interpolatedChannels = referenceOut.badChannels;
%         referenceOut.referenceSignal = nanmean(signalNew.data(referenceOut.referenceChannels, :), 1);
        referenceOut.referenceSignal = mean(signalNew.data(referenceOut.referenceChannels, :), 1,'omitnan');
        clear signalNew;
        signal = removeReference(signal, referenceOut.referenceSignal, ...
            referenceOut.rereferencedChannels);
        
        textbar.endprogressbar('done');
        %         sig2 = pop_select(signal,'nopoint',remove_range);
        %
        %         referenceOut.noisyStatistics  = ...
        %             findNoisyChannels(removeTrend(sig2, referenceOut), referenceOut);
    end


    function [] = doOther()
        % Do some type of non-robust referencing.
        if strcmpi(referenceOut.referenceType, 'average') && ...
                (~isempty( setdiff(referenceOut.evaluationChannels, ...
                referenceOut.referenceChannels)) ...
                || ~isempty( setdiff(referenceOut.referenceChannels, ...
                referenceOut.evaluationChannels)))
            warning('averageReference:EvaluationChannels', ...
                'Reference and evaluation channels should be same for average reference');
        elseif strcmpi(referenceOut.referenceType, 'specific') && ...
                (isempty( setdiff(referenceOut.evaluationChannels, ...
                referenceOut.referenceChannels))  && ...
                isempty( setdiff(referenceOut.referenceChannels, ...
                referenceOut.evaluationChannels)))
            warning('specificReference:EvaluationChannels', ...
                'Reference and evaluation channels should not be same for specific reference');
        elseif ~strcmpi(referenceOut.referenceType, 'average') && ...
                ~strcmpi(referenceOut.referenceType, 'specific') && ...
                ~strcmpi(referenceOut.referenceType, 'none')
            warning('specificReference:BadReferenceType', 'Unrecognized reference type');
        end;
        referenceOut.noisyStatisticsOriginal  = ...
            findNoisyChannels(removeTrend(signal, referenceOut), referenceOut);
        referenceOut.noisyStatisticsBeforeInterpolation = ...
            referenceOut.noisyStatisticsOriginal;
        referenceOut.referenceSignal = referenceOut.referenceSignalOriginal;
        signal = removeReference(signal, referenceOut.referenceSignal, ...
            referenceOut.rereferencedChannels);
        if  ~strcmpi(referenceOut.referenceType, 'none')
            referenceOut.noisyStatistics  = findNoisyChannels( ...
                removeTrend(signal, referenceOut), referenceOut);
        else
            referenceOut.noisyStatistics = referenceOut.noisyStatisticsOriginal;
        end
    end
end
