function [EEG] = performICA(EEG, EOGParams, isblink, autoICA)

fprintf('EOG artifacts removal through ICA...\n');

if nargin <3
    isblink = true;
    autoICA = false; % default: user-based IC selection
end

% EOG channels are specified in DefaultParameters.m as EOGParams struct
if isstruct(EOGParams)
    if isfield(EOGParams,'channels')
        EOGParams = EOGParams.channels;
    else
        fprintf('Param struct required the field ''channels''');
    end
end

if ischar(EOGParams)
    EOGParams = strsplit(EOGParams,',');
end


if isblink
    % Bad channels and bad segments are removed before performing ICA 
    [~,EEG2] = evalc('pop_select(EEG,''nochannel'',find(EEG.badchan))');
    EEG2.badchan(EEG.badchan) = [];
    [~,EEG2] = evalc('pop_select(EEG2,''nopoint'',EEG.remove_range)');
    EEG2.bad_segment(EEG.bad_segment) = [];

    [~,~,eogchannels] = intersect(lower(EOGParams),lower({EEG2.chanlocs.labels}), 'stable');

    if isempty(eogchannels)
        warning('No good channels EOG regression');
        return;
    end

    EOGParams = {EEG2.chanlocs(eogchannels).labels};

    if isempty(eogchannels)
        fprintf('No channels for EOG regression...no regression performed\n');
        return;
    end

    params.signalLabels = EOGParams;
    params.signalNumbers = eogchannels(:)';
    params.signalTypeIndicator = 'UseLabels';
    params.dumpBlinkerStructures = false;
    params.showMaxDistribution = false;
    params.dumpBlinkImages = false;
    params.dumpBlinkPositions = false;
    params.verbose = false;
    params.goodRatioThreshold = .7;
    params.minGoodBlinks = 10;
    params.blinkAmpRange = [3 50];

    ready = false;
    iteration = 0;
    while ~ready
        iteration = iteration + 1;
        fprintf('\tIteration %d\n',iteration);
        warning off
        [~, com, blinks, blinkFits, blinkProperties] = eval('pop_blinker(EEG2,params,true)');
        warning on
        if blinks.usedSignal>0
            ready = true;
        else
            if blinks.usedSignal == -1
                params.goodRatioThreshold = params.goodRatioThreshold - .1;
                if params.goodRatioThreshold < .4
                    fprintf('No blinks detected...no regression performed\n');
                    return;
                end
            elseif blinks.usedSignal == -2
                params.minGoodBlinks = params.minGoodBlinks - 1;
                if params.minGoodBlinks < 5
                    fprintf('No blinks detected...no regression performed\n');
                    return;
                end
            elseif blinks.usedSignal == -3
                params.blinkAmpRange(1) = params.blinkAmpRange(1) - .1;
                if params.blinkAmpRange(1) < 2
                    fprintf('No blinks detected...no regression performed\n');
                    return;
                end

            end
        end
    end
    [~,EEG2] = evalc('pop_eegfiltnew(EEG2, 1, 15)');

    leftbase = [blinkFits.leftBase];
    rightbase = [blinkFits.rightBase];
    points = [];
    for iblink = 1:length(leftbase)
        points = [points leftbase(iblink):rightbase(iblink)];
    end

    data = EEG2.data(:,points);

    fprintf("\tSVD decomposition...\n");
    [U,S,V] = svd(data,'econ');

    EEG2.icachansind = 1:EEG2.nbchan;
    EEG2.icaweights = pinv(U);
    EEG2.icawinv = (U);
    EEG2.icasphere = eye(size(U));
    EEG2.reject.gcompreject = false(1,size(EEG2.icawinv,2));

    fprintf("\tICLabel classification...\n");
    features = ICL_feature_extractor(EEG2, true);
    fprintf( sprintf('\t\tcalculating labels...\n'));
    iclabels = run_ICL('default', features{:});
    icclasses = {'Brain', 'Muscle', 'Eye', 'Heart', ...
        'Line Noise', 'Channel Noise', 'Other'};
    [~,idx] = max(iclabels,[],2);
    EEG2.etc.ic_classification.ICLabel.classifications = iclabels;
    EEG2.etc.ic_classification.ICLabel.classes = icclasses;

    maskU = false(1,size(U,2));

    [~,max_p_idx] = max(EEG2.etc.ic_classification.ICLabel.classifications,[],2);

    eye_classification = (max_p_idx == 3)';
    eye_probability = EEG2.etc.ic_classification.ICLabel.classifications(:,3)';

    if sum(eye_classification)>0
        fprintf("\tICLabel found %d EYE components with a probability > %d%%\n",sum(eye_classification),floor(100*min(eye_probability(eye_classification))));
    else
        fprintf("\tICLabel found no EYE components \n");
    end

    maskU = eye_classification;

    if sum(maskU)>0
        fprintf("\tOnly %d out of %d EYE components with a probability > 90%%\n\tPerform manual labelling\n",sum(maskU & eye_probability>.9),sum(maskU));
        tmp = U(:,maskU);
        U(:,maskU) = [];
        U = [tmp U];
        EEG2.icaweights = pinv(U);
        EEG2.icawinv = (U);
        EEG2.icasphere = eye(size(U));
        EEG2.reject.gcompreject(maskU) = true;
        tmp = EEG2.reject.gcompreject(maskU);
        EEG2.reject.gcompreject(maskU) = [];
        EEG2.reject.gcompreject = [tmp EEG2.reject.gcompreject];
        tmp = EEG2.etc.ic_classification.ICLabel.classifications(maskU,:);
        EEG2.etc.ic_classification.ICLabel.classifications(maskU,:) = [];
        EEG2.etc.ic_classification.ICLabel.classifications = [tmp; EEG2.etc.ic_classification.ICLabel.classifications];

        EEG2.data = wlb_bandpass_filter(EEG2.data',EEG2.srate, [0 10])';
        icplt = icaplotter(EEG2);
        % if ~autoICA 
        waitfor(gcf);
        % end
        maskU = icplt.selectionComp;
    else
        fprintf("\tNo EYE components detected \n\tPerform manual labelling\n");
        EEG2.data = wlb_bandpass_filter(EEG2.data',EEG2.srate, [0 10])';
        icplt = icaplotter(EEG2);
        % if ~autoICA 
        waitfor(gcf);
        % end
        maskU = icplt.selectionComp;
    end

    U = U(:,maskU);

    if ~isempty(U) && ~all(U(:) == 0)
        % Reorthogonalize the vectors
        [U,S,V] = svd(U,0);
        S = diag(S);
        % Enforce strict zero values (because Matlab SVD function can randomly return small values instead of zeros !!!)
        % If not, sometimes it adds small contributions of some channels that have nothing to do with the projectors (<1e-16).
        % If this added channel contains high values (eg. Stim channel in Volts) it can corrupt significantly the MEG values.
        iZero = find((abs(U) < 1e-15) & (U ~= 0));
        if ~isempty(iZero)
            U(iZero) = 0;
        end
        % Throw away the linearly dependent guys (threshold on singular values: 0.01 * the first one)
        iThresh = find(S < 0.01 * S(1),1);
        if ~isempty(iThresh)
            fprintf('SSP> %d linearly dependent vectors removed...', size(U,2)-iThresh+1);
            U = U(:, 1:iThresh-1);
        end
        % Compute projector in the form: I-UUt
        EEG.data(~EEG.badchan,:) = EEG.data(~EEG.badchan,:) - U*wlb_bandpass_hfilter(U'*EEG.data(~EEG.badchan,:),EEG.srate,0,10);
    end
end

EEG.pipeline.EOGregression.performed = 'yes';
EEG.pipeline.EOGregression.channels = EOGParams;
end