function [varargout]=nt_zapline (x,srate,fline,nfft,nkeep,dsskeep,plotflag,textbar,mask)
%    x: data matrix Chans x Time
%    srate: sample rate
%    fline:line frequency principal from which the harmonics will be computed
%    nfft: size of FFT [default:1024]
%    nkeep: number of components to keep in DSS [default: all]
%    dsskeep: number of components to keep to describe the line noise in dss if empty try an automatic selection (default: [])


%testcode
if nargin ==0
    sr=400;
    nsamples=100000; nchans=100;
    signal=randn(nchans,nsamples);
    artifact=sin((1:nsamples)/sr*2*pi*50);
    artifact=max(artifact,0).^3; % introduce harmonics
    artifact=3*nt_demean(randn(nchans,1)*artifact);
    disp(nt_wpwr(artifact)/nt_wpwr(signal+artifact));
    [varargout] = nt_zapline(signal+artifact,sr,50,[],[],3);
    return;
end

assert(nargin>=3, '!');

if nargin<4||isempty(nfft); nfft=1024; end
if nargin<5||isempty(nkeep); nkeep=[]; end
if nargin<6||isempty(dsskeep); dsskeep=[]; end

if nargin<7||isempty(plotflag)
    plotflag=0;
end
if nargin<8||isempty(textbar)
    textbar=[];
end
if nargin<9
    mask=[];
end
if isempty(x); error('!'); end
if nkeep>=size(x,1); error('!'); end
fline = fline.freq;
if fline>srate/2; error('fline should be less than Nyquist'); end
if size(x,2)<nfft
    warning(['reducing nfft to ',num2str(size(x,2))]); 
    nfft=size(x,2);
end

fline = fline/(srate);

nHarmonics=floor((1/2)/fline);

if(~isempty(textbar))
    textbar.updatebar(0);
end

oldx = x;
emdexpansion = false;

if(size(x,1)<32)
    emdexpansion = true;
end

if emdexpansion
    yy=[];
    for ich = 1:size(x,1)
        IMF = emd(x(ich,:));
        numIMF(ich) = size(IMF,2);
        yy = [yy IMF];
    end
    x = yy';
end
% xx=removeTone(x',srate,fline*srate*(1:nHarmonics),.5);
xx=nt_smooth(x',1/fline,1,1); % cancels line_frequency and harmonics, light lowpass
% xx=nt_smooth_filt(x',1/fline,fline,1); % cancels line_frequency and harmonics, light lowpass

if(~isempty(textbar))
    textbar.updatebar(16.66);
end

if isempty(nkeep); nkeep=size(x,1); end

kk = x'-xx;
clear x; 

xxxx=nt_pca(kk,[],nkeep,[],mask); % reduce dimensionality to avoid overfitting

if(~isempty(textbar))
    textbar.updatebar(33.33);
end

% xxxx = x-xx;
% DSS to isolate line components from residual:
[c0,c1]=nt_bias_fft(xxxx,fline*(1:nHarmonics), nfft);

if(~isempty(textbar))
    textbar.updatebar(50);
end

[todss,pwr0,pwr1]=nt_dss0(c0,c1);

if(~isempty(textbar))
    textbar.updatebar(66.66);
end

assert(size(todss,2)>0, '!'); 
 
score = (pwr1)./(pwr0);
score = (score-min(score))/(max(score)-min(score));

if isempty(dsskeep)
%     fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
%     p = sigm_fit(1:length(score),score,[],[],0);
%     score = 1-fsigm(p,1:length(score))/p(2);
    idx = find(score<.01,1,'first');
    if (isempty(idx))
        idx = length(score);
    end
else
    idx = dsskeep;
end
if isempty(idx)
    idx = 1;
end

xxxxx=nt_mmat(xxxx,todss(:,1:idx)); % line-dominated components

if(~isempty(textbar))
    textbar.updatebar(83.33);
end

[xxx]=nt_tsr(kk,xxxxx,[],mask,mask); % project them out

if(~isempty(textbar))
    textbar.updatebar(100);
end

clear xxxx

% reconstruct clean signal
y=xx+xxx; clear xx xxx

if emdexpansion
    
    yy=zeros(size(oldx'));
    
    for ich = 1:length(numIMF)
        yy(:,ich) = sum(y(:,1:numIMF(ich)),2);
        y(:,1:numIMF(ich)) = [];
    end
    y = yy;
end
yy=oldx'-y;

% yy = yy - mean(yy);
% 
% xxx=nt_smooth(yy,1/fline,10,1); % cancels line_frequency and harmonics, light lowpass
% xxx=nt_smooth_filt(yy,1/fline,10,1); % cancels line_frequency and harmonics, light lowpass
% % 
% yy = yy-xxx;
% 
% y=oldx'-yy;

varargout{1} = y';
varargout{2} = yy';
varargout{3} = score;
varargout{4} = idx;
varargout{5} = todss;

if ~nargout || plotflag
    % print result and display spectra    
%     nfft = size(x,2);
    figure;clf;subplot(223);
    plot((score), '.-'); xlabel('component'); ylabel('score'); title('DSS to enhance line frequencies');
    hold on
    plot(idx,score(idx),'og','MarkerFaceColor','g')
    h = subplot(221);
    [pxx,f]=nt_spect_plot(x',20*srate,10*srate,[],1/fline);
    divisor=mean(pxx(:));
    semilogy(f,bsxfun(@rdivide,abs(pxx),divisor));
    title('original'); %legend boxoff
    set(gca,'ygrid','on','xgrid','on');
    xlabel('frequency (relative to line)');
    ylabel('relative power');
    yl1=get(gca,'ylim');
    hh=get(gca,'children');
    h(2)=subplot(222);
    [pyy,f]=nt_spect_plot(y,20*srate,10*srate,[],1/fline);
    semilogy(f,bsxfun(@rdivide,abs(abs(pyy)),divisor));
    title('clean')
    set(gca,'ygrid','on','xgrid','on');
    set(gca,'yticklabel',[]); ylabel([]);
    xlabel('frequency (relative to line)');
    yl2=get(gca,'ylim');
    
    h(3)=subplot(224);
    [pxx,f]=nt_spect_plot(yy,20*srate,10*srate,[],1/fline);
    semilogy(f,bsxfun(@rdivide,abs(pxx),divisor));
    title('removed'); %legend boxoff
    set(gca,'ygrid','on','xgrid','on');
    set(gca,'yticklabel',[]); ylabel([]);
    xlabel('frequency (relative to line)');
    yl3=get(gca,'ylim');
    hh=get(gca,'children');
    
   	yl(1)=min([yl1(1),yl2(1),yl3(1)]); yl(2)=max([yl1(2),yl2(2),yl3(2)]);
    subplot 221; ylim(yl); subplot 222; ylim(yl); subplot 224; ylim(yl);
    linkaxes(h)
    drawnow;

end