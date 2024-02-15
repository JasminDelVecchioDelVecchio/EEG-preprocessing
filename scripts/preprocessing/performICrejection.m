function [EEG] = performICrejection(EEG)

options.figSize = 'small';
if EEG.trials > 1
    options.plotTimeX = [EEG.xmin EEG.xmax];
    
else
    options.plotTimeX = [0 5];
end

options.plotFreqX = [1 100];
plt = icaplotter(EEG,'figSize',options.figSize,'plotTimeX',options.plotTimeX,'plotFreqX',options.plotFreqX);

waitfor(gcf);

EEG.reject.gcompreject = plt.selectionComp; 

EEG.pipeline.ICrej.performed = 'yes';
EEG.pipeline.ICrej.rejectComponents = find(plt.selectionComp);
