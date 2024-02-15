function [EEG,events] = loadDataFromBrainstorm(subj_name)
% Load EEG as EEGLAB data struct from brainstorm

[sub,sub_idx] = bst_get('Subject', subj_name);

[studies,istudies] = bst_get('StudyWithSubject', sub.FileName);

[study_idx,~] = listdlg('ListString',{studies.Name}.');

resting_istudy = istudies(study_idx);
resting_data   = studies(study_idx).Data;

resting_ChannelMat = in_bst_channel(studies(study_idx).Channel.FileName);

LoadOptions.IgnoreBad      = 0;  % From raw files: ignore the bad segments
LoadOptions.ProcessName    = 'process_unknown';
LoadOptions.RemoveBaseline = 'no';
LoadOptions.UseSsp         = 1;


resting_DataMat    = in_bst_data(resting_data.FileName);

[sMat, nSignals, iRows] = bst_process('LoadInputFile', resting_data.FileName, [], [], LoadOptions);
[DataMat, matName] = in_bst(resting_data.FileName, [], 0);
events=DataMat.F.events;

EEG = bst2eeglab(sMat.Data,sMat.Time,events,1:size(sMat.Data,1),resting_ChannelMat,DataMat.ChannelFlag);

end