function [b] = myfirws(Fs, HighPass, LowPass, isRelax, isstop)
if nargin<5 || isempty(isstop)
    isstop=false;
end

Nyquist = Fs/2;
% High-pass filter
if ~isempty(HighPass) && all(HighPass ~= 0)
    if numel(HighPass)==2
        HighPass = sort(HighPass);
        highTranBand = diff(HighPass);
        HighPass = HighPass(2);
    else
        highTranBand = [];
    end
    f_highpass = HighPass / Nyquist;    % Change frequency from Hz to normalized scale (0-1)
    if isempty(highTranBand) || highTranBand==0
        if (HighPass <= 5)
            highTranBand = .5 ; %Hz
        else
            highTranBand = 1 ; %Hz
        end
        f_highstop = f_highpass - highTranBand/Nyquist;
    else
        f_highstop = max(0, HighPass - highTranBand) / Nyquist;
        % f_highstop = max(0.2, HighPass - TranBand) / Nyquist;
        highTranBand   = (f_highpass - f_highstop)*Nyquist ;  % Adjusted Transition band
    end
else
    f_highpass = 0;
    f_highstop = 0;
    highTranBand = 1 ;
end
% Low-pass filter
if ~isempty(LowPass) && all(LowPass ~= 0)
    if numel(LowPass)==2
        LowPass = sort(LowPass);
        lowTranBand = diff(LowPass);
        LowPass = LowPass(1);
    else
        lowTranBand = [];
    end
    f_lowpass = LowPass / Nyquist;
    if isempty(lowTranBand) || lowTranBand==0
        lowTranBand = 1;
        f_lowstop = f_lowpass + lowTranBand/Nyquist;
    else
        f_lowstop = f_lowpass + lowTranBand/Nyquist;
    end
else
    f_lowpass  = 0;
    f_lowstop  = 0;
end
% If both high-pass and low-pass are zero
if (f_highpass == 0) && (f_lowpass == 0)
    return;
    % Input frequencies are too high
elseif (f_highpass >= 1) || (f_lowpass >= 1)
    return;
end
% Transition parameters
if isRelax
    Ripple = 10^(-2);
    Atten  = 10^(-2);   % Equals 40db
else
    Ripple = 10^(-3);   % pass band ripple
    Atten  = 10^(-3);   % Equals 60db
end

% ===== DESIGN FILTER =====
% Build the general case first
fcuts = [f_highstop, f_highpass, f_lowpass, f_lowstop];
if isstop
    mags = [1 0 1];
    devs = [Ripple 10^-4 Ripple];
else
    mags  = [0 1 0];               % filter magnitudes
    devs  = [Atten Ripple Atten];  % deviations
end
% Now adjust for desired properties
fcuts = max(0,fcuts);      % Can't go below zero
fcuts = min(1-eps, fcuts); % Can't go above or equal to 1

% We have implicitly created a bandpass, but now adjust for desired filter
if (f_lowpass == 0)  % User didn't want a lowpass
    fcuts(3:4) = [];
    mags(3) = [];
    devs(3) = [];
end
if (f_highpass == 0)  % User didn't want a highpass
    fcuts(1:2) = [];
    mags(1) = [];
    devs(1) = [];
end

% Generate FIR filter
% Using Matlab's Signal Processing toolbox
[n,Wn,beta,ftype] = kaiserord(fcuts, mags, devs, 2);
n = n + rem(n,2);  % ensure even order
b = fir1(n, Wn, ftype, kaiser(n+1,beta), 'noscale');
end