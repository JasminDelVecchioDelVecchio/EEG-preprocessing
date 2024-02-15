function [x,FiltSpec] = wlb_bandstop_filter(x, sfreq, FreqList,FreqWidth,order,isMirror)
if (nargin < 6) || isempty(isMirror)
    isMirror = 0;
end


% Nyqist frequency
Fnyq = sfreq/2;

if(not(isstruct(FreqList)))
    FiltSpec = [];
    % Check list of freq to remove
    if isempty(FreqList) || isequal(FreqList, 0)
        return;
    end
    
    % Remove all the frequencies sequencially
    for ifreq = 1:length(FreqList)
        % Frequency band to remove
        FreqBand = [FreqList(ifreq) - FreqWidth/2, FreqList(ifreq) + FreqWidth/2];
        % Filter order
        % Butterworth filter
        [FiltSpec(ifreq).B,FiltSpec(ifreq).A] = butter(order, FreqBand ./ Fnyq, 'stop');
        FiltSpec(ifreq).mirror = isMirror;
        FiltSpec(ifreq).order = order;
        
    end
else
    FiltSpec = FreqList;
end

if ~isempty(x)
    filter_data = true;
    [ nTime ,nChan] = size(x);
    xmean = mean(x);
    x = bsxfun(@minus, x, xmean);
    for ifreq = 1:length(FreqList)
        M = length(FiltSpec(ifreq).B);
        if isMirror
            x = [flipud(x(1:M,:)); x; flipud(x(end-M+1:end,:))];
            % Zero-padding
        else
            x = [zeros(M,nChan); x; zeros(M,nChan)] ;
        end
        % Filter signal
        x = filtfilt(FiltSpec(ifreq).B, FiltSpec(ifreq).A, x);
        x = x(M+1:end-M,:);

    end
    x = bsxfun(@plus, x, xmean);
    
end
% Restore the mean of the signal

