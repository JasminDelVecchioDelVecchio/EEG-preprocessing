function EEG = performZapline(EEG,fline,dsskeep)

if nargin<3
    dsskeep = [];
end
fprintf('Removing power line with Zapline...\n');

toRemove = find(EEG.badchan);
bad_segments = EEG.bad_segment;

[~,EEGCleaned] = evalc('pop_select(EEG,''nochannel'',toRemove)');
oldremovedMask = EEG.badchan;

nfft = 2^nextpow2(5*EEG.srate);
nkeep = [];

if ~isempty(dsskeep)
    textbar = textprogressbar(sprintf('\tzapline : '));
    [data_cleaned,artefact,score,dsskeep,todss] = eval(['nt_zapline(EEGCleaned.data(:,:), EEGCleaned.srate, ' ...
        'double(fline.freq), nfft, nkeep, dsskeep, 0, textbar, ~bad_segments'')']);
    textbar.endprogressbar('done');
else
    textbar = textprogressbar(sprintf('\tzapline : '));
    [plt]  = eval(['ZaplinePlotter(EEGCleaned.data(:,:), EEGCleaned.srate, fline, nfft, dsskeep,' ...
        'nkeep,textbar, ~bad_segments'')']);
    % waitfor(gcf)
    data_cleaned = plt.data_cleaned;

    clear plt
    textbar.endprogressbar('done');

end


EEG.data(~oldremovedMask,:) = data_cleaned;
EEG.pipeline.Zapline.performed = 'yes';
EEG.pipeline.Zapline.fline = fline;
end


