classdef ZaplinePlotter < handle
    properties
        fig
        h
        score_ax
        data
        data_cleaned
        divisor
        srate
        fline
        nfft
        dsskeep
        nkeep
        todss
        textbar
        mask
        numiter
    end
    methods
        
        function obj = ZaplinePlotter(data,srate,fline,nfft,dsskeep,nkeep,txtbar,mask,numiter)
            
            if nargin==0
                srate=400;
                nsamples=100000; 
                nchans=10;
                signal=randn(nchans,nsamples);
                artifact=sin((1:nsamples)/srate*2*pi*50);
                artifact=max(artifact,0).^3; % introduce harmonics
                artifact=3*nt_demean(randn(nchans,1)*artifact);
                disp(nt_wpwr(artifact)/nt_wpwr(signal+artifact));
                nfft = 2^nextpow2(1*srate);
                obj = ZaplinePlotter(signal+artifact,srate,50,nfft);
                return
            end
            
            if nargin<4||isempty(nfft); nfft=1024; end
            if nargin<5||isempty(dsskeep); dsskeep=[]; end
            if nargin<6||isempty(nkeep); nkeep=[]; end
            if nargin<7||isempty(txtbar)
                txtbar=[];
            end
            if nargin<8 || isempty(mask)
                mask=true(size(data,2),1);
            end
            if nargin<9 || isempty(numiter)
                numiter=1;
            end
            % Plot
            obj.fig = figure;
            obj.h = subplot(221); % cleaned signal
            obj.data = data;
            obj.srate = srate;
            obj.nfft = nfft;
            obj.dsskeep = dsskeep;
            obj.nkeep = nkeep;
            obj.todss = [];
            obj.fline = fline;
            obj.mask = mask;
            obj.textbar = txtbar;
            obj.numiter = numiter;
            
            if ~isempty(obj.textbar)
                obj.textbar.updatebar(0);
            end
            
            [pxx,f]=nt_spect_plot(obj.data(:,obj.mask)',2*obj.srate,obj.srate,[],obj.srate);
            obj.divisor=mean(pxx(:));
            semilogy(f,bsxfun(@rdivide,abs(pxx),obj.divisor));
            title('original'); %legend boxoff
            set(obj.h(1),'ygrid','on','xgrid','on');
            xlabel('frequency (Hz)');
            ylabel('relative power');
            
            obj.h(2)=subplot(222);cla(obj.h(2))
            obj.h(3)=subplot(224);cla(obj.h(3))
            obj.score_ax = subplot(223);cla(obj.score_ax)

            if ~isempty(obj.textbar)
                obj.textbar.updatebar(5);
            end
            
            obj.compute
            
        end
        
        function compute(obj)
            x = obj.data;
            
            if ~isempty(obj.textbar)
                obj.textbar.min = 5;
                obj.textbar.max = 90;
            end
            
            for iter = 1:obj.numiter
                [obj.data_cleaned,artefact,score,obj.dsskeep,obj.todss] = nt_zapline (x,obj.srate,obj.fline,obj.nfft,obj.nkeep,obj.dsskeep,0,obj.textbar,obj.mask);
                x = obj.data_cleaned;
            end
            cla(obj.score_ax)
            plot(score, '.-','HitTest', 'on', ...
                'Parent', obj.score_ax, 'LineWidth', 1,'ButtonDownFcn',@(o,e) obj.onClickedScore(o,e));
            xlabel('component'); ylabel('score'); title(obj.score_ax,'DSS to enhance line frequencies');
            hold(obj.score_ax,'on')
            plot(obj.dsskeep,score(obj.dsskeep),'og','MarkerFaceColor','g','Parent', obj.score_ax)
            
            [pyy,f]=nt_spect_plot(obj.data_cleaned(:,obj.mask)',2*obj.srate,obj.srate,[],obj.srate);
            semilogy(f,bsxfun(@rdivide,abs(abs(pyy)),obj.divisor),'Parent', obj.h(2));
            title(obj.h(2),'clean')
            set(obj.h(2),'ygrid','on','xgrid','on');
            set(obj.h(2),'yticklabel',[]); ylabel([]);
            xlabel('frequency (Hz)');
            
            if ~isempty(obj.textbar)
                obj.textbar.min = 0;
                obj.textbar.max = 100;
                
                obj.textbar.updatebar(95);
            end
            
            [pxx,f]=nt_spect_plot(artefact(:,obj.mask)',2*obj.srate,obj.srate,[],obj.srate);
            semilogy(f,bsxfun(@rdivide,abs(pxx),obj.divisor),'Parent', obj.h(3));
            title(obj.h(3),'removed'); %legend boxoff
            set(obj.h(3),'ygrid','on','xgrid','on');
            xlabel('frequency (Hz)');

            if ~isempty(obj.textbar)
                obj.textbar.min = 0;
                obj.textbar.max = 100;
                
                obj.textbar.updatebar(100);
            end
            yl1=get(obj.h(1),'ylim');
            yl2=get(obj.h(2),'ylim');
            yl3=get(obj.h(3),'ylim');
            
            yl(1)=min([yl1(1),yl2(1),yl3(1)]); yl(2)=max([yl1(2),yl2(2),yl3(2)]);
            
            set(obj.h(1),'ylim',yl);
            set(obj.h(2),'ylim',yl);
            set(obj.h(3),'ylim',yl);
            
            linkaxes(obj.h)

        end
        
        function onClickedScore(obj,o,e)
                        
            obj.dsskeep = round(e.IntersectionPoint(1));
                        
            obj.compute
            
        end
        
        
    end
    
    
end