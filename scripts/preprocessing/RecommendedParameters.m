classdef RecommendedParameters
    %RecommendedParameters is a class containing recommended parameters for
    %   different preprocessing steps. This is different from the
    %   DefaultParameters.m in the sense that there is a recommended value
    %   even for preprocessing steps that are not used by default: These
    %   structures must have all the possible fields for each parameter
    %   (with the exception of PrepParams), so that they are used later in
    %   the GUI or other places where no default parameter is given but
    %   another recommendation is needed. 
    %
    %   Please do not change anything in this file unless you are sure what
    %   you are doing.
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
        
        ResampleParams = struct('freq',250);

        FilterParams = struct('notch',    struct('freq', []),...
                             'high',     struct('freq', [0.25 1],...
                                                 'order', []),... % Default
                             'low',      struct('freq', [80 85],...
                                                 'order', []))    % Default
       
         EOGParams = struct('channels','Fp1,Fpz,Fp2');

         PrepParams = struct('highFrequencyNoiseThreshold',   5,...
             'correlationThreshold',       0.4,...
             'detrendType' , 'none');

         CRDParams = struct('ChannelCriterion',     0.85,...
                               'LineNoiseCriterion',   4,...
                               'BurstCriterion',       25,...
                               'WindowCriterion',      0.25, ...
                               'Highpass',             []);
                
        EpochingParams = struct('EpochDuration',2,'EventNames',[],'EvetnWdw',[],'Baseline',[]);
        
        FastICAParams = struct( 'FilterParams',...
                                        struct('notch',    struct([]),...
                                               'high',     struct('freq', [.25 1.0],'order', []),...
                                               'low',      struct('freq', [80 85],'order', [])));
        
        prunedata = true;
        
        subcompICA = true;
        
        interpbadchannels = true
    end
end