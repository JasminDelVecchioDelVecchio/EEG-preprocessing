function EEG = performFilter(EEG, varargin)
% performFilter  perform a high pass, low pass and notch filter.
%
%   filtered = performFilter(EEG, params)
%   where EEG is the EEGLAB data structure. filtered is the resulting 
%   EEGLAB data structured after filtering. params is an optional
%   parameter which must be a structure with optional parameters
%   'notch', 'high' and 'low', each of which a struct. An example of this
%   parameter is given below:
%   params = struct('notch', struct('freq', 50),...
%                   'high',  struct('freq', 0.5, 'order', []),...
%                   'low',   struct('freq', 30,  'order', []))
%   
%   'notch.freq' is the frequency for the notch filter where from
%   (notch_freq - 3) to (notch_freq + 3) is attenued.
%
%   'high.freq' and 'high.order' are the frequency and filtering order for 
%   high pass filter respectively.
%
%   'low.freq' and 'low.order' are the frequency and filtering order for 
%   low pass filter respectively.
%
%   In the case of filtering ordering, if it is left to be high.order = []
%   (or low.order = []), then the default value of pop_eegfiltnew.m is
%   used.
%
%   If params is ommited default values are used. If any field of params
%   are ommited, corresponding default values are used. If 
%   'params.notch = struct([])', 'params.high = struct([])' or 
%   'params.low = struct([])' then notch filter, high pass filter or 
%   low pass filter are not perfomed respectively.
%
%   Default values are specified in DefaultParameters.m. If they are empty
%   then defaults of inexact_alm_rpca.m are used.
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

%% Parse parameters
p = inputParser;
addOptional(p,'notch', struct([]), @isstruct);
addOptional(p,'high', struct([]), @isstruct);
addOptional(p,'low', struct([]), @isstruct);
parse(p, varargin{:});
notch = p.Results.notch;
high = p.Results.high;
low = p.Results.low;

%% Perform filtering
if( ~isempty(high) || ~isempty(low) || ~isempty(notch))
    EEG.pipeline.filtering.performed = 'yes';
    if( ~isempty(high) )
%         [~, EEG, ~ , b] = evalc('pop_eegfiltnew(EEG, high.freq, 0, high.order)');
        [EEG, ~, b] = mypop_firws(EEG, high.freq, [], 0);
        EEG.pipeline.filtering.highpass.performed = 'yes';
        EEG.pipeline.filtering.highpass.freq = high.freq;
        EEG.pipeline.filtering.highpass.order = length(b)-1; 
        EEG.pipeline.filtering.highpass.transitionBandWidth = diff(high.freq);
    end

    if ( ~isempty(low) )
%         [~, EEG, ~ , b] = evalc('pop_eegfiltnew(EEG, 0, low.freq, low.order)');
        [EEG, ~, b] = mypop_firws(EEG, [], low.freq, 0);
        EEG.pipeline.filtering.lowpass.performed = 'yes';
        EEG.pipeline.filtering.lowpass.freq = low.freq;
        EEG.pipeline.filtering.lowpass.order = length(b)-1; 
        EEG.pipeline.filtering.lowpass.transitionBandWidth = diff(low.freq);
    end

    if( ~isempty(notch) )
%         notch.freq = str2double(strsplit(notch.freq,{' ',',',';'}));
        fline = notch.freq(1)/(EEG.srate);
        nHarmonics=floor((1/2)/fline);

        notch_freq = fline*EEG.srate*(1:nHarmonics);
        for ifreq = 1:length(notch_freq)
            if notch_freq(ifreq)+notch.freq(2)/2<EEG.srate/2;
                [EEG, ~, b] = mypop_firws(EEG, notch_freq(ifreq)-notch.freq(2)/2, notch_freq(ifreq)+notch.freq(2)/2, 0,1);
            else
                [EEG, ~, b] = mypop_firws(EEG, 0,[EEG.srate/2-5 EEG.srate/2], 0,0);
            end
        end
        EEG.pipeline.filtering.notch.performed = 'yes';
        EEG.pipeline.filtering.notch.freq = notch.freq;
        EEG.pipeline.filtering.notch.order = length(b)-1; 
        EEG.pipeline.filtering.notch.transitionBandWidth = 1;
    end
end

end