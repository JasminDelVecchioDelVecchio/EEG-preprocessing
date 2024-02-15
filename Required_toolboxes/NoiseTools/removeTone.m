function [r] = removeTone(x,Fs,Fl,bandwidth)


N=length(x);            %no of data

xMean = mean(x);%calculate dc
Fl(Fl>(Fs/2)-2)=[];
v = 0;
nChan = size(x,2);
sigma = 1/bandwidth*Fs;

for i=1:length(Fl)
    xMean = mean(x);
    x = bsxfun(@minus,x,xMean);
    if Fl(i)+bandwidth<Fs/2
        b = myfirws(Fs, Fl(i)-bandwidth, Fl(i)+bandwidth, [], 1);
    else
        b = myfirws(Fs, 0,[Fl(i)-2*bandwidth, Fl(i)], [], 0);
    end
    M = (length(b)-1)/2;
    x1 = [zeros(M,nChan); x; zeros(M,nChan)] ;
    x1=filter(b,1,x1);
    x = x1(2*M+1:end,:);

%     [B,A] = butter(3,[Fl(i)-bandwidth, min(Fl(i)+bandwidth,0.5*Fs)]/(0.5*Fs),'stop');
%     x = filtfilt(B,A,x);

    x = bsxfun(@plus,x,xMean);

end

r=x-v;    %I wasn't sure what was meant by removing, so i just subtracted

function [f,findx]=getfgrid(Fs,nfft,fpass)
% Helper function that gets the frequency grid associated with a given fft based computation
% Called by spectral estimation routines to generate the frequency axes 
% Usage: [f,findx]=getfgrid(Fs,nfft,fpass)
% Inputs:
% Fs        (sampling frequency associated with the data)-required
% nfft      (number of points in fft)-required
% fpass     (band of frequencies at which the fft is being calculated [fmin fmax] in Hz)-required
% Outputs:
% f         (frequencies)
% findx     (index of the frequencies in the full frequency grid). e.g.: If
% Fs=1000, and nfft=1048, an fft calculation generates 512 frequencies
% between 0 and 500 (i.e. Fs/2) Hz. Now if fpass=[0 100], findx will
% contain the indices in the frequency grid corresponding to frequencies <
% 100 Hz. In the case fpass=[0 500], findx=[1 512].
if nargin < 3; error('Need all arguments'); end;
df=Fs/nfft;
f=0:df:Fs; % all possible frequencies
f=f(1:nfft);
if length(fpass)~=1;
   findx=find(f>=fpass(1) & f<=fpass(end));
else
   [fmin,findx]=min(abs(f-fpass));
   clear fmin
end;
f=f(findx);
