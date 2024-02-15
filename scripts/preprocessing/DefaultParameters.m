classdef DefaultParameters
    %DefaultParameters is a class containing default parameters for
    %   different preprocessing steps.
    %
    %   struct([]) desactivate the corresponding operation.
    %   struct() will be the default parameters used in the corresponding
    %   operation.
    %   If a field has [] as value, then the default value in the
    %   corresponding function is used.
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
    properties(Constant)

        ResampleParams = struct('freq',[],'factor',2);

        FilterParams = struct('notch',    struct('freq', 50),...
            'high',     struct('freq', [0.25 1]),...
            'low',      struct('freq', [80 85]))

        EOGParams = struct('channels','Fp1,Fpz,Fp2,AF7,F7,AF8,F8');

        ECGParams = struct('channels','ECG');

        PrepParams = struct('highFrequencyNoiseThreshold',   5,...
            'correlationThreshold',       0.4,...
            'detrendType' , 'none'); % Not needed because highpass is already implemented as first filtering step (see "FilterParams")

        CRDParams = struct('ChannelCriterion',     'off',...%already implemented in robustReference step
            'LineNoiseCriterion',   'off',...%already implemented in robustReference step
            'BurstCriterion',       25,...
            'WindowCriterion',      0.25, ...
            'Highpass',             'off'); %highpass already implemented as first filtering step (see "FilterParams")

        EpochingParams = struct('EpochDuration',2,'EventNames',[],'EvetnWdw',[],'Baseline',[]);

        FastICAParams = struct( 'FilterParams',...
            struct('notch',    struct([]),...
            'high',     struct('freq', [.25 1.0]),...
            'low',      struct('freq', [80 85])));

        prunedata = true; %if true the pipeline performs the final steps:
        %removes bad IC, reconstruct signal and interpolate bad channels

        autoICA = false; 

        perform_ICA = true; % if true performs ICA for artifacts removal, otherwise performs EOG/ECG regression

        subcompICA = true; %if is true removes bad IC, reconstruct signal (valid only if "prunedata" is true too)

        interpbadchannels = true; %if is true interpolate bad channels (valid only if "prunedata" is true too)

        save_steps = true;

        save_fig = true;

    end
end