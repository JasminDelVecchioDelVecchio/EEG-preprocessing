function [y,w]=wlb_detrend(x,wsize,overlap,order,thresh,niter,w)
%[y,w,r]=nt_detrend(x,order,w,basis,thresh,niter,wsize) - robustly remove trend
% 
%  y: detrended data
%  w: updated weights
%  r: basis matrix used
%
%  x: raw data
%  order: order of polynomial or number of sin/cosine pairs
%  w: weights
%  basis: 'polynomials' 
%  thresh: threshold for outliers [default: 3 sd]
%  niter: number of iterations [default: 3]
%  wsize: window size for local detrending [default: all]
%
% This NEW (circa Oct 2019) version of detrend allows detrending to be
% applied to smaller overlapping windows, which are then stitched together
% using overlap-add.


%% arguments
if nargin<2; wsize = []; end
if nargin<3; overlap=[]; end
basis='polynomials';
if nargin<4||isempty(order); order=3; end
if nargin<5||isempty(thresh); thresh=3; end
if nargin<6||isempty(niter); niter=3; end
if nargin<7; w=[]; end

dims=size(x); 
flipsignal = false;

if dims(1)<dims(2)
    x = x';
    dims=size(x);
    flipsignal = true;
end

nchans=dims(2);
x=x(:,:); % concatenates dims >= 2
w=w(:,:);

x = demean(x);

if isempty(wsize) || ~wsize
    % standard detrending (trend fit to entire data)
    [y,w,r]=single_detrend(x,order,w,basis,thresh,niter);
elseif wsize==dims(1)
    [y,w,r]=single_detrend(x,order,w,basis,thresh,niter);
else 
    wsize=2*floor(wsize/2);
    if isempty(overlap)
        overlap = round(0.5*wsize);
    end
    
    if overlap>=wsize
        error('overlap should be less than window length')
    end

    winstep = wsize - overlap;
    
    % do some sanity checks because many parameters
    if numel(order)>1; error('!'); end
    if ~isempty(w) && ~(size(w,1)==size(x,1)) ; disp(size(w)); error('!'); end
    if thresh==0; error('thresh=0 is not what you want...'); end % common mistake
    if numel(thresh)>1; error('!'); end
    if numel(niter)>1; error('!'); end

    if isempty(w); w=ones(size(x)); end
    if size(w,2)==1; w=repmat(w,1,size(x,2)); end

    % (1) divide into windows, (2) detrend each, (3) stitch together, (4)
    % estimate w

    for iIter=1:niter

        y=zeros(size(x));
        trend=zeros(size(x));
        a=zeros(size(x,1),1);

    %     figure(1); clf
%         overlap = wsize/2;
        offset=0;
        while true
            start=offset+1;
            stop=min(size(x,1),offset+wsize);

            % if not enough valid samples grow window:
            counter=5;
            while any (sum(min(w(start:stop),2)) <wsize)
                if counter <= 0 ; break; end 
                start=max(1,start-wsize/2);
                stop=min(size(x,1),stop+wsize/2);
                counter=counter-1;
            end
            if rem(stop-start+1,2)==1; stop=stop-1; end
            wsize2=stop-start+1;

            % detrend this window
            NITER=1;

            yy=single_detrend(x(start:stop,:),order,w(start:stop,:),basis,thresh,NITER);

            % triangular weighting
            if start==1
                wt=mytukeywin(wsize2,2*overlap/wsize2);
                b=[ones(floor(wsize2/2),1); wt((floor(wsize2/2)+1):end)];
            elseif stop==size(x,1)
                wt=mytukeywin(wsize2,2*overlap/wsize2);
                b=[wt(1:(floor(wsize2/2))); ones(wsize2/2,1)];
            else
                b=mytukeywin(wsize2,2*overlap/wsize2);
            end

            % overlap-add to output
            y(start:stop,:)=y(start:stop,:)+bsxfun(@times,yy,b);
            trend(start:stop,:)=trend(start:stop,:)+bsxfun(@times,x(start:stop,:)-yy,b);

            a(start:stop,1)=a(start:stop)+b;

            offset=offset+winstep;
            if offset>size(x,1)-winstep; break; end
        end
        
        if stop<size(x,1); y(end,:)=y(end-1,:); a(end,:)=a(end-1,:); end; % last sample can be missed
        
        y=bsxfun(@times,y,1./a); y((isnan(y)))=0;
        trend=bsxfun(@times,trend,1./a);  trend((isnan(trend)))=0;

        % find outliers
        d=x-trend; 


        if ~isempty(w); d=bsxfun(@times,d,w); end
        ww=ones(size(x));
        ww(find(abs(d)>thresh*repmat(std(d),size(x,1),1))) = 0;
        clear d

        % update weights
        if isempty(w); 
            w=ww;
        else
            w=min(w,ww);
        end
        clear ww;

    end % for iIter...
end % if isempty(wsize)...

if flipsignal
    y = y';
end

if ~nargout
    % don't return, just plot
    subplot 411; plot(x); title('raw'); xlim([1,size(x,1)])
    subplot 412; plot(y); title('detrended'); xlim([1,size(x,1)])
    subplot 413; plot(x-y'); title('trend'); xlim([1,size(x,1)])
    subplot 414; nt_imagescc(w'); title('weight');
    clear y w r
end
end

%% Original version of detrend (no windows) is called by new version (windows)
function [y,w,r]=single_detrend(x,order,w,basis,thresh,niter)
%[y,w,r]=nt_detrend(x,order,w,basis,thresh,niter) - robustly remove trend
% 
%  y: detrended data
%  w: updated weights
%  r: basis matrix used
%
%  x: raw data
%  order: order of polynomial or number of sin/cosine pairs
%  w: weights
%  basis: 'polynomials' [default] or 'sinusoids', or user-provided matrix
%  thresh: threshold for outliers [default: 3 sd]
%  niter: number of iterations [default: 3]
%
% Noise tools
% See nt_regw().
%
% The data are fit to the basis using weighted least squares. The weight is
% updated by setting samples for which the residual is greater than 'thresh' 
% times its std to zero, and the fit is repeated at most 'niter'-1 times.
%
% The choice of order (and basis) determines what complexity of the trend
% that can be removed.  It may be useful to first detrend with a low order
% to avoid fitting outliers, and then increase the order.
%
% Examples:
% Fit linear trend, ignoring samples > 3*sd from it, and remove:
%   y=nt_detrend(x,1); 
% Fit/remove polynomial order=5 with initial weighting w, threshold = 4*sd:
%   y=nt_detrend(x,5,w,[],4);
% Fit/remove linear then 3rd order polynomial:
%   [y,w]=nt_detrend(x,1);
%   [yy,ww]=nt_detrend(y,3);
%

% arguments
if nargin<2; error('!'); end
if nargin<3; w=[]; end
if nargin<4||isempty(basis); basis='polynomials'; end
if nargin<5||isempty(thresh); thresh=3; end
if nargin<6||isempty(niter); niter=3; end

if thresh==0; error('thresh=0 is not what you want...'); end % common mistake

dims=size(x);
x=x(:,:); % concatenates dims >= 2
w=w(:,:);

if size(w,2)==1; w=repmat(w,1,size(x,2)); end

%% regressors
r=zeros(size(x,1),numel(order));
lin=linspace(-1,1,size(x,1));
for k=1:order
    r(:,k)=lin.^k;
end

% remove trends
% The tricky bit is to ensure that weighted means are removed before
% calculating the regression (see nt_regw).

for iIter=1:niter
    
    % weighted regression on basis
    [~,y]=weigthed_regression(x,r,w);
    
    % find outliers
    d=x-y; 
    if ~isempty(w); d=bsxfun(@times,d,w); end
    ww=ones(size(x));
    ww(find(abs(d)>thresh*repmat(std(d),size(x,1),1))) = 0;
     
    % update weights
    if isempty(w); 
        w=ww;
    else
        w=min(w,ww);
    end
    clear ww;    
end
y=x-y;
y=reshape(y,dims);
w=reshape(w,dims);

end

function [b,z]=weigthed_regression(y,x,w)
%[b,z]=nt_regw(y,x,w) - weighted regression
%
%  b: regression matrix (apply to x to approximate y)
%  z: regression (x*b)
%
%  y: data
%  x: regressor
%  w: weight to apply to y
%
%  w is either a matrix of same size as y, or a column vector to be applied
%  to each column of y
%
% NoiseTools

PCA_THRESH=0.0000001; % discard dimensions of x with eigenvalue lower than this

if nargin<3; w=[]; end
if nargin<2; error('!'); end

%% check/fix sizes
mn=y-demean(y,w);
%%
if isempty(w) 
    %% simple regression
    xx=demean(x);
    yy=demean(y);
    [V,D]=eig(xx'*xx); V=real(V); D=real(D);
    topcs=V(:,find(D/max(D) > PCA_THRESH)); % discard weak dims
    xxx=xx*topcs;
    b=(yy'*xxx) / (xxx'*xxx); b=b';
    if nargout>1; z=demean(x,w)*topcs*b; z=z+mn; end
else
    %% weighted regression
    if size(w,1)~=size(x,1); error('!'); end
    if size(w,2)==1; 
        %% same weight for all channels
        if sum(w(:))==0; 
            warning('weights all zero');
            b=zeros(size(x,2),1);
            topcs = zeros(size(x,2));
        else
            yy=demean(y,w).*repmat(w,1,size(y,2)); 
            xx=demean(x,w).*repmat(w,1,size(x,2));  
            [V,D]=eig(xx'*xx); V=real(V); D=real(D); D=diag(D);
            topcs=V(:,find(D/max(D) > PCA_THRESH)); % discard weak dims
            xxx=xx*topcs;
            b=(yy'*xxx) / (xxx'*xxx); b=b';
        end
        if nargout>1; z=demean(x,w)*topcs*b; z=z+mn; end
    else
        %% each channel has own weight
        if size(w,2) ~= size(y,2); error('!'); end 
        if nargout; z=zeros(size(y)); end
        for iChan=1:size(y,2)
            if sum(w(:,iChan))==0; %disp(iChan); 
                %warning('weights all zero'); 
                b=zeros(size(y,2));
                topcs = 0;
            else
                yy=demean(y(:,iChan),w(:,iChan)) .* w(:,iChan); 
                x=demean(x,w(:,iChan)); % remove channel-specific-weighted mean from regressor
                xx=x.*repmat(w(:,iChan),1,size(x,2)); 
                [V,D]=eig(xx'*xx); V=real(V); D=real(diag(D));
                topcs=V(:,find(D/max(D) > PCA_THRESH)); % discard weak dims
                xxx=xx*topcs;
                b(iChan,1:size(topcs,2))=(yy'*xxx) / (xxx'*xxx); 
            end
            if nargout>1; z(:,iChan)=x*(topcs*b(iChan,1:size(topcs,2))') + mn(:,iChan); end
        end
    end             
end
end

function [x,mn]=demean(x,w)
%[y,mn]=nt_demean(x,w) - remove weighted mean over cols
% 
%  w is optional
%
%  if w is a vector with fewer samples than size(x,1), it is interpreted as
%  a vector of indices to be set to 1, the others being set to 0.
%
% NoiseTools

if nargin<2; w=[]; end

[m,n,o]=size(x);

if isempty(w)
    
    mn=mean(double(x),1);
    x=bsxfun(@minus,x,mn);
    
else
        
    if size(w,2)==1
        mn=sum(bsxfun(@times,double(x),w),1) ./ (sum(w,1)+eps);
    elseif size(w,2)==n
        mn=sum(x.*w) ./ (sum(w,1)+eps);
    else
        error('W should have same ncols as X, or else 1');
    end

    x=bsxfun(@minus,x,mn);
    
end
end
