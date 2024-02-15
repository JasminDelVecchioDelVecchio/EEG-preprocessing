function EEG = performFastICA(EEG, varargin)
 
% Parse and check parameters
defaults = DefaultParameters.FastICAParams;
p = inputParser;
addParameter(p,'FilterParams', defaults.FilterParams, @isstruct);
parse(p, varargin{:});
params = p.Results;
FilterParams = params.FilterParams;

chanind = setdiff( 1:EEG.nbchan,find(EEG.badchan));

if( ~isempty(FilterParams) )
    EEGICA = performFilter(EEG, 'high',FilterParams.high,'low',FilterParams.low);
    EEGICA.pipeline.FastICA.filter.performed = 'yes';
    EEGICA.pipeline.FastICA.filter.filterparams = FilterParams;
else
    EEGICA = EEG;
    EEGICA.pipeline.FastICA.filter.performed = 'no';
end

if(EEGICA.trials>1)
    EEGICA = pop_select( EEGICA,'notrial',find(EEGICA.bad_epoch));
else
    EEGICA = pop_select( EEGICA,'nopoint',EEGICA.remove_range);
end

if isfield(EEGICA,'icaweights') && isempty(EEGICA.icaweights)
    EEGICA.data = EEGICA.data/mean(std(EEGICA.data(:,:),[],2));
    [EEGICA] = pop_runica(EEGICA, 'icatype','runica','chanind',chanind);
end

act = EEGICA.icaweights*EEGICA.icasphere*EEGICA.data(EEGICA.icachansind,:);
squaresig  = sum(sum(EEGICA.data(EEGICA.icachansind,:).^2));
varsPerc = zeros(1,size(EEGICA.icaweights,1));
icawinv    = pinv(EEGICA.icaweights*EEGICA.icasphere);

for i=1:size(EEGICA.icaweights,1)
    compproj = icawinv(:,i)*act(i,:);
    varsPerc(i) = 100*(sum(sum(compproj.^2))/squaresig);
end

[sortvar, windex] = sort(varsPerc,'descend');

EEGICA.icawinv = icawinv(:,windex); %alters inverse weights matrix based on component variance order
EEGICA.icaweights = EEGICA.icaweights(windex,:); %alters weights matrix based on component variance order

EEG.icaweights = EEGICA.icaweights;
EEG.icawinv = EEGICA.icawinv;
EEG.icasphere = EEGICA.icasphere;
EEG.icachansind=EEGICA.icachansind;

EEG.pipeline.FastICA.performed = 'yes';
EEG.pipeline.FastICA.params = params;


end
