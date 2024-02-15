classdef CleaningPipelineParameterDef < handle
    
    properties
        handles
        results
    end
    
    methods
       
        function obj = CleaningPipelineParameterDef(varargin)
            
            Defaults = DefaultParameters;
            p = inputParser;
            addParameter(p,'FilterParams', Defaults.FilterParams, @isstruct);
            addParameter(p,'NTParams', Defaults.NTParams, @isstruct);
            addParameter(p,'CRDParams', Defaults.CRDParams, @isstruct);
            addParameter(p,'PrepParams', Defaults.PrepParams, @isstruct);
            parse(p, varargin{:});
            FilterParams = p.Results.FilterParams;
            CRDParams = p.Results.CRDParams;
            NTParams = p.Results.NTParams;
            PrepParams = p.Results.PrepParams;

            obj.results.NTParams = NTParams;
            obj.results.CRDParams = CRDParams;
            obj.results.PrepParams = PrepParams;
            obj.results.FilterParams = FilterParams;
            
            obj.handles.f = figure;

            obj.handles.NoiseToolFrequencies_label = uicontrol('Style','text','String','Line noise [Hz]:',...
                'units','normalize','FontUnits','normalize','position',[.1 .9 .25 .05],'HorizontalAlignment', 'left');
            obj.handles.NoiseToolFrequencies_text = uicontrol('Style','edit','Min',1,'Max',1,'String',num2str(NTParams.lineFrequency),...
                'units','normalize','FontUnits','normalize','position',[.1 .85 .25 .05]);
            obj.handles.NoiseToolnfft_label = uicontrol('Style','text','String','nfft:',...
                'units','normalize','FontUnits','normalize','position',[.375 .9 .2 .05],'HorizontalAlignment', 'left');
            obj.handles.NoiseToolnfft_text = uicontrol('Style','edit','Min',1,'Max',1,'String',num2str(NTParams.nfft),...
                'units','normalize','FontUnits','normalize','position',[.375 .85 .2 .05]);
            obj.handles.NoiseTooldsskeep_label = uicontrol('Style','text','String','N° component [empty = auto]:',...
                'units','normalize','FontUnits','normalize','position',[.605 .9 .35 .05],'HorizontalAlignment', 'left');
            obj.handles.NoiseTooldsskeep_text = uicontrol('Style','edit','Min',1,'Max',1,'String',num2str(NTParams.dsskeep),...
                'units','normalize','FontUnits','normalize','position',[.605 .85 .3 .05]);
            
            obj.handles.correlationThreshold_label = uicontrol('Style','text','String','Minimum channel correlation:',...
                'units','normalize','FontUnits','normalize','position',[.1 .8-.075 .9 .05],'HorizontalAlignment', 'left');
            obj.handles.correlationThreshold_text = uicontrol('Style','edit','Min',1,'Max',1,'String',num2str(PrepParams.correlationThreshold),...
                'units','normalize','FontUnits','normalize','position',[.1 .75-0.075 .2 .05]);
            obj.handles.highFrequencyNoiseThreshold_label = uicontrol('Style','text','String','HF Noise Criterion standard deviations:',...
                'units','normalize','FontUnits','normalize','position',[.1 .7-0.075 .9 .05],'HorizontalAlignment', 'left');
            obj.handles.highFrequencyNoiseThreshold_text = uicontrol('Style','edit','Min',1,'Max',1,'String',num2str(PrepParams.highFrequencyNoiseThreshold),...
                'units','normalize','FontUnits','normalize','position',[.1 .65-0.075 .2 .05]);
            
            obj.handles.BurstCriterion_label1= uicontrol('Style','text','String','Burst Criterion Standard deviation cutoff: Default: 5',...
                'units','normalize','FontUnits','normalize','position',[.1 .6-0.075 .9 .05],'HorizontalAlignment', 'left');
            obj.handles.BurstCriterion_label2 = uicontrol('Style','text','String','Reasonable range: 3 (mimum,very aggressive) to 40 (very lax)',...
                'units','normalize','FontUnits','normalize','position',[.1 .55-0.075 .9 .05],'HorizontalAlignment', 'left');
            obj.handles.BurstCriterion_text = uicontrol('Style','edit','Min',1,'Max',1,'String',num2str(CRDParams.BurstCriterion),...
                'units','normalize','FontUnits','normalize','position',[.1 .5-0.075 .2 .05]);
            obj.handles.WindowCriterion_label1 = uicontrol('Style','text','String','Window Criterion: Default: 0.25.',...
                'units','normalize','FontUnits','normalize','position',[.1 .45-0.075 .9 .05],'HorizontalAlignment', 'left');
            obj.handles.WindowCriterion_label2 = uicontrol('Style','text','String','Reasonable range: 0.05 (very aggressive) to 0.3 (very lax)',...
                'units','normalize','FontUnits','normalize','position',[.1 .4-0.075 .9 .05],'HorizontalAlignment', 'left');
            obj.handles.WindowCriterion_text = uicontrol('Style','edit','Min',1,'Max',1,'String',num2str(CRDParams.WindowCriterion),...
                'units','normalize','FontUnits','normalize','position',[.1 .35-0.075 .2 .05]);

            obj.handles.HighPass_label = uicontrol('Style','text','String','High Pass [Hz]:',...
                'units','normalize','FontUnits','normalize','position',[.1 .3-0.12 .3 .05],'HorizontalAlignment', 'left');
            obj.handles.HighPass_text = uicontrol('Style','edit','Min',1,'Max',1,'String',num2str(FilterParams.high.freq),...
                'units','normalize','FontUnits','normalize','position',[.1 .25-0.12 .2 .05]);
            obj.handles.LowPass_label = uicontrol('Style','text','String','Low Pass [Hz]:',...
                'units','normalize','FontUnits','normalize','position',[.1 .2-0.12 .3 .05],'HorizontalAlignment', 'left');
            obj.handles.LowPass_text = uicontrol('Style','edit','Min',1,'Max',1,'String',num2str(FilterParams.low.freq),...
                'units','normalize','FontUnits','normalize','position',[.1 .15-0.12 .2 .05]);
            
            obj.handles.ok_button = uicontrol('Style','pushbutton','String','Ok','units','normalize',...
                'FontUnits','normalize','position',[.65 .15 .1 .05],'Callback',@(h,e)obj.ok_button_Callback(h,e));
            obj.handles.cancel_button = uicontrol('Style','pushbutton','String','Cancel','units','normalize',...
                'FontUnits','normalize','position',[.5 .15 .1 .05],'Callback',@(h,e)obj.cancel_button_Callback(h,e));

        end
        
        function cancel_button_Callback(obj,h,e)
            obj.results = [];
            close(obj.handles.f);
        end
        
        function ok_button_Callback(obj,h,e)

            obj.results.NTParams.lineFrequency = str2double((get(obj.handles.NoiseToolFrequencies_text,'String')));
            obj.results.NTParams.nfft = str2double((get(obj.handles.NoiseToolnfft_text,'String')));
            obj.results.NTParams.dsskeep = str2double((get(obj.handles.NoiseTooldsskeep_text,'String')));
            if(isnan(obj.results.NTParams.dsskeep))
                obj.results.NTParams.dsskeep = [];
            end
            obj.results.PrepParams.correlationThreshold = str2double((get(obj.handles.correlationThreshold_text,'String')));
            obj.results.PrepParams.highFrequencyNoiseThreshold = str2double((get(obj.handles.highFrequencyNoiseThreshold_text,'String')));
            
            obj.results.CRDParams.BurstCriterion = str2double((get(obj.handles.BurstCriterion_text,'String')));
            obj.results.CRDParams.WindowCriterion = str2double((get(obj.handles.WindowCriterion_text,'String')));
            obj.results.CRDParams.Highpass = [];
            
            obj.results.FilterParams.high.freq = str2double((get(obj.handles.HighPass_text,'String')));
            obj.results.FilterParams.low.freq = str2double((get(obj.handles.LowPass_text,'String')));
            
            close(obj.handles.f)
        end
    end
    
    
end