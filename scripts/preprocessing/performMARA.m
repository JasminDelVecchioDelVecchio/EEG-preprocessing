  function EEG = performMARA(EEG, varargin)
% performMARA  perform Independent Component Analysis (ICA) on the high
%   passsed data and classifies bad components using MARA.
%   This function applies a high pass filter before the ICA. But the output
%   result is NOT high passed filtered, rather only cleaned with ICA. This
%   option allows to choose a separate high pass filter only for ICA from
%   the desired high pass filtered after the entire preprocessing. Please
%   note that at this stage of the preprocessing, another high pass filter
%   has been already applied on the data in the performFilter.m. Please use
%   accordingly.
%
%   data_out = performMARA(data, params) where data is the input EEGLAB data
%   structure and data_out is the output EEGLAB data structure after ICA.
%   params is an optional parameter which must be a structure with optional
%   fields 'chanlocMap' and 'high'. An example of params is given below:
%
%   params = = struct('chanlocMap', containers.Map, ...
%                     'largeMap',   0, ...
%                     'high',       struct('freq', 1.0, 'order', []))
%
%   params.chanlocMap must be a map (of type containers.Map) which maps all
%   "possible" current channel labels to the standard channel labels given
%   by FPz, F3, Fz, F4, Cz, Oz, ... as required by processMARA. Please note
%   that if the channel labels are already the same as the mentionned
%   standard, an empty map would be enough. However if the map is empty and
%   none of the labels has the same semantic as required, no ICA will be
%   applied. For more information please see processMARA. An example of
%   such a map is given in systemDependentParse.m where a containers.Map is
%   created for the MARAParams.chanlocMap in the case of EGI systems.
%   Sometimes this field may be ignored, but here then it get replaced with
%   a new empty mapping.
%
%   params.high is a structure indicating the high pass frequency
%   (params.high.freq) and order (params.high.order) of the high pass
%   filter applied on the data before ICA. For more information on this
%   parameter please see performFilter.m
%
%   If varargin is ommited, default values are used. If any fields of
%   varargin is ommited, corresponsing default value is used.
%
%   Default values are taken from DefaultParameters.m.
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

defaults = DefaultParameters.MARAParams;
% chanlocMap could be inexistant. This might be necessary so that in the
% SystemDependentParse.m the default mapping is not assigned. But here it
% does not matter if it is existent or not as nothing will happen depending
% on this. Existence and non existence of this field is only necessary for
% the above mentioned function.
if ~isfield(defaults, 'chanlocMap')
    defaults.chanlocMap = RecommendedParameters.MARAParams.chanlocMap;
end
%% Parse and check parameters
p = inputParser;
validate_param = @(x) isa(x, 'containers.Map');
addParameter(p,'chanlocMap', defaults.chanlocMap, validate_param);
addParameter(p,'largeMap', defaults.largeMap);
addParameter(p,'REQ_CHAN_LABELS', defaults.REQ_CHAN_LABELS);
parse(p, varargin{:});
chanlocMap = p.Results.chanlocMap;
REQ_CHAN_LABELS = p.Results.REQ_CHAN_LABELS;
% Change channel labels to their corresponding ones as required by
% processMARA. This is done only for those labels that are given in the map.
if( ~ isempty(chanlocMap))
    inverseChanlocMap = containers.Map(chanlocMap.values, ...
        chanlocMap.keys);
    idx = find(ismember({EEG.chanlocs.labels}, chanlocMap.keys));
    for i = idx
        EEG.chanlocs(1,i).labels = chanlocMap(EEG.chanlocs(1,i).labels);
    end
    
    % Temporarily change the name of all other labels to make sure they
    % don't create conflicts
    for i = 1:length(EEG.chanlocs)
        if(~ any(i == idx))
            EEG.chanlocs(1,i).labels = strcat(EEG.chanlocs(1,i).labels, ...
                '_automagiced');
        end
    end
end

% Check if the channel system is according to what Mara is expecting.
intersect_labels = intersect(cellstr(REQ_CHAN_LABELS), ...
    {EEG.chanlocs.labels});
if(length(intersect_labels) < 3)
    msg = ['The channel location system was very probably ',...
        'wrong and MARA ICA could not be used correctly.' '\n' 'MARA ICA for this ',...
        'file is skipped.'];
    ME = MException('Automagic:MARA:notEnoughChannels', msg);
    
    % Change back the labels to the original one
    if( ~ isempty(chanlocMap))
        for i = idx
            EEG.chanlocs(1,i).labels = inverseChanlocMap(...
                EEG.chanlocs(1,i).labels);
        end
        
        for i = 1:length(EEG.chanlocs)
            if(~ any(i == idx))
                EEG.chanlocs(1,i).labels = strtok(...
                    EEG.chanlocs(1,i).labels, '_automagiced');
            end
        end
    end
    EEG.pipeline.mara.performed = 'no';
    throw(ME)
end

%% Perform ICA

EEGMara = pop_select( EEG,'notrial',find(EEG.bad_epoch));

[artcomps,MARAinfo] = MARA(EEGMara);
EEG.reject.MARAinfo = MARAinfo;

EEG.reject.gcompreject = false(1,size(EEG.icawinv,2));
EEG.reject.gcompreject(artcomps) = true;


end
