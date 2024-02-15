function EEG = performIClabel(EEG,varargin)
%
%% Parse and check parameters
p = inputParser;
addParameter(p,'EpochingParams', struct([]), @isstruct);
parse(p, varargin{:});
params = p.Results;
EpochingParams = params.EpochingParams;

if ~isempty(EpochingParams)
    EpochDuration = EpochingParams.EpochDuration;
    EventNames = EpochingParams.EventNames;
    EvetnWdw = EpochingParams.EvetnWdw;
    Baseline = EpochingParams.Baseline;
    if ~isempty(EpochDuration)
        if EpochDuration > 0
            [~,EEGICA] = evalc('mypop_epoch(EEG,''EpochDuration'',EpochDuration)');
        end
    else
        [~,EEGICA] = evalc('mypop_epoch(EEG,''EventNames'',EventNames,''EvetnWdw'',EvetnWdw)');        
        if(not(isempty(Baseline)))
            [~,EEGICA] = evalc('pop_rmbase( EEGICA, Baseline*1000)');
        end
    end
else
    EEGICA=EEG;
end

[~,EEGICA] = evalc('pop_select( EEGICA,''notrial'',find(EEGICA.bad_epoch))');
fprintf('Performing ICLabel detection...\n');

features = my_ICL_feature_extractor(EEGICA, true);
clear EEGICA;

fprintf( sprintf('\tcalculating labels...\n'));
iclabels = run_ICL('default', features{:});
icclasses = {'Brain', 'Muscle', 'Eye', 'Heart', ...
    'Line Noise', 'Channel Noise', 'Other'};
[~,idx] = max(iclabels,[],2);
EEG.reject.gcompreject = false(1,size(EEG.icawinv,2));
EEG.reject.gcompreject(idx>1 & idx <7) = true;
EEG.etc.ic_classification.ICLabel.classifications = iclabels;
EEG.etc.ic_classification.ICLabel.classes = icclasses;

EEG.pipeline.ICrej.performed = 'yes';
EEG.pipeline.ICrej.iclabels = iclabels;
EEG.pipeline.ICrej.icclasses = icclasses;
EEG.pipeline.ICrej.rejectComponents = find(EEG.reject.gcompreject);


end
