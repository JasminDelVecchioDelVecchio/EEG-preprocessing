function [v] = markers_velocity(tracks, h,  frequency)
% Computes markers' velocity using the second order forward differentiation
% (approx with finite differences)

% tracks - vector having #samples x #coordinates
% h - incremental step
% frequency - sample frequency
    
    flag = 0;
    
    % Most likely the number of coordinates is lower than the samples
    if size(tracks,1) < size(tracks,2)
        tracks = tracks';
        flag = 1;
    end
    
    tracks = [tracks; zeros(h,size(tracks,2))+tracks(end,:)];
    
    nFRM    = size(tracks,1);
    nTracks = size(tracks,2);
    
    v = zeros(nFRM-h, nTracks);
    
    for ii = 2:nFRM-h
        for jj = 1:nTracks
%             v(ii,jj) = (tracks(ii+h,jj) - tracks(ii,jj))*frequency; % 1st order forward differentiation (error ~ 1/fs)
            v(ii,jj) = (tracks(ii+h,jj) - tracks(ii-h,jj))*0.5*frequency; % 2nd order forward differentiation (error ~ (1/fs).^2)
        end
    end 
    
    if flag 
        v = v';
    end
end

