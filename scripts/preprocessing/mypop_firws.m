% pop_firws() - Filter data using windowed sinc FIR filter
function [EEG, com, b] = mypop_firws(EEG, HighPass, LowPass, isRelax, isstop,minphase,plotresp)

if isstruct(EEG)
    if isempty(EEG.data)
        error('Cannot process empty dataset');
    else
        Fs = EEG.srate;
    end
else
    Fs = EEG;
    EEG = [];
    if isempty(Fs)
        error('Cannot define filter');
    end
end
if nargin<4 || isempty(isRelax)
    isRelax = false;
end
if nargin<5 || isempty(isstop)
    isstop = false;
end
if nargin<6 || isempty(minphase)
    minphase = false;
end
if nargin<7 || isempty(plotresp)
    plotresp = false;
end

b = myfirws(Fs, HighPass, LowPass, isRelax, isstop);

% Filter
if isstruct(EEG)
    m = mean(EEG.data,2);
    EEG.data = bsxfun(@minus,EEG.data,m);

    if minphase
        b = minphaserceps(b);
        EEG = firfiltsplit(EEG, b, 1);
    else
        EEG = firfilt(EEG, b);
    end
    if isempty(HighPass) || any(HighPass==0) || isstop
        EEG.data = bsxfun(@plus,EEG.data,m);
    end
    % History string
    com = sprintf('%s = pop_firws(%s', inputname(1), inputname(1));
%     for c = fieldnames(args)'
%         if ischar(args.(c{:}))
%             com = [com sprintf(', ''%s'', ''%s''', c{:}, args.(c{:}))];
%         else
%             com = [com sprintf(', ''%s'', %s', c{:}, mat2str(args.(c{:})))];
%         end
%     end
    com = [com ');'];
end
if plotresp
    plotfresp(b, 1, [], EEG.srate, minphase);
end