function [EEG,yy]=pop_nt_zapline(EEG,p,plotflag)
%[y,yy]=nt_zapline(x,fline,nremove,p,plotflag) - remove power line artifact
%
%  y: denoised data
%  yy: artifact
%
%  x: data
% %  p: additional parameters:
%    p.lineFrequency: power line frequency 
%    p.nfft: size of FFT [default:1024]
%    p.nkeep: number of components to keep in DSS [default: all]
%    p.dsskeep: number of component to keep to describe the line noise in dss
%    p.fig1: figure to use for DSS score [default: 100]
%    p.fig2: figure to use for results [default: 101]
%  plotflag: plot
%
%Examples:
%  nt_zapline(x,60/1000) 
%    apply to x, assuming line frequency=60Hz and sampling rate=1000Hz, plot results
%  nt_zapline(x,60/1000,4)
%    same, removing 4 line-dominated components 
%  p=[];p.nkeep=30; nt_zapline(x,60/1000,4,p);
%    same, truncating PCs beyond the 30th to avoid overfitting
%  [y,yy]=nt_zapline(x,60/1000)
%    return cleaned data in y, noise in yy, don't plot
%

% NoiseTools
try, nt_greetings; catch, disp('You must download NoiseNools from http://audition.ens.fr/adc/NoiseTools/'); return; end

assert(nargin>=1, '!'); 
if nargin<2||isempty(p); p=struct(); end
if ~isfield(p,'fline'); p.lineFrequency=50; end
if ~isfield(p,'nfft'); p.nfft=1024; end
if ~isfield(p,'nkeep'); p.nkeep=[]; end
if ~isfield(p,'dsskeep'); p.dsskeep=[]; end
if nargin<3||isempty(plotflag); plotflag=0; end

if isempty(EEG.data); error('!'); end
if p.nkeep>=size(EEG.data,1); error('!'); end
if p.lineFrequency>EEG.srate/2; error('lineFrequency should be less than Nyquist'); end
if size(EEG.data,2)<p.nfft
    warning(['reducing nfft to ',num2str(size(EEG.data,2))]); 
    p.nfft=size(EEG.data,2);
end

toRemove = find(EEG.badchan);

EEGCleaned = pop_select(EEG,'nochannel',toRemove);
EEGCleaned = pop_select(EEGCleaned,'nopoint',EEG.remove_range);

fprintf('Removing Line noise...\n');

p.nfft = EEGCleaned.pnts*EEGCleaned.trials;
 
if rem(p.nfft,2) == 1
    p.nfft = p.nfft-1;
end

p.nfft = min(p.nfft ,2^nextpow2(round(20*EEGCleaned.srate)));
p.dsskeep = 20;

[~,EEG.data(~EEG.badchan,~EEG.bad_segment)] = evalc('nt_zapline (EEGCleaned.data,EEGCleaned.srate,p.lineFrequency,p.nfft,p.nkeep,p.dsskeep,1)');
