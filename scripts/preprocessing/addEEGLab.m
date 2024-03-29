function parts = addEEGLab(libraryPath)
% addEEGLab  Add the eeglab package to the path
%   eeglab package is assumed to be in folder matlab_scripts located in the
%   parent directory. The need for this function is to remove the path to
%   few folders of eeglab which make conflicts with other MATLAB functions.
%   This functions is strongly recommended over the use of
%   `addpath(genpath('path/to/eeglab'))`.
%
%   parts = addEEGLab()
%
%   parts: all the paths in EEGLAB folder. This is returned for the cases
%   where it is needed.
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


folderName = 'eeglab14_1_2b';

folderName = fullfile(libraryPath,folderName);

addpath(genpath(folderName));

    
end