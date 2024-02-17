% Clear workspace
clc; clear; close all;

% addCleaningPipelinePaths()

% Check if brainstorm is installed
if exist('bst_get') == 2
    disp('Brainstorm is installed and available.');
else
    disp('Brainstorm is not installed or not available');
end

% Specify required paths
pathData = 'C:\Users\jasmi\Lab Isaias Dropbox\Team\WUMIGE_shared\WUMIGE_shared_projects\Dystonia_jdvdv\EEG_analysis\CleaningPipeline\EEG-preprocessing\data\';
pathRes = 'C:\Users\jasmi\Lab Isaias Dropbox\Team\WUMIGE_shared\WUMIGE_shared_projects\Dystonia_jdvdv\EEG_analysis\CleaningPipeline\EEG-preprocessing\data\plots\';

% Load data
subject_ID = 'CR78';
fprintf("Loading data...\n")
fprintf(subject_ID)
[EEG,events] = loadDataFromBrainstorm(subject_ID);

% Preprocessing EEG data
[EEGCleaned, EEG_steps] = MainCleaningPipeline(EEG,subject_ID,pathRes);
save([pathRes,subject_ID,'_cleaned_pipeline.mat'],'EEGCleaned','-v7.3');

Fs = EEGCleaned.srate;
iChannels = find(strcmp({EEGCleaned.chanlocs.type},'EEG')); % EEG channels
iChannels_LFP = find(strcmp({EEGCleaned.chanlocs.type}, 'STN')); % LFPs channels
iChannels_KIN = find(strcmp({EEG.chanlocs.type}, 'MISC')); % kinematic data channels
iChannel_ECG = find(strcmp({EEGCleaned.chanlocs.type}, 'ECG')); % ECG channel

% Remove LFP-ECG artifact
ECG = EEGCleaned.data(iChannel_ECG,:);
LFP = EEGCleaned.data(iChannels_LFP,:);
[eventQRS] = detect_qrs(ECG, 1, length(ECG), Fs); % identify qrs peaks
time = (1:length(LFP))./Fs;
[LFP_cleaned] = my_ecg_svd(LFP, Fs, 1, eventQRS, time);

% Remove tremor artifact
if ~isempty(iChannels_KIN)
    % Marker tracks identification
    loi={'HEAD_MX_x_kin';'HEAD_RX_x_kin';'HEAD_LX_x_kin'};
    kin_labels = {EEGCleaned.chanlocs(iChannels_KIN).labels}';
    markers_tracks = EEGCleaned.data(iChannels_KIN,:)';

    [index_loi]=extract_markers(kin_labels',loi)';
    head=markers_tracks(:,index_loi(1):index_loi(1)+2);
    head_rx=markers_tracks(:,index_loi(2):index_loi(2)+2);
    head_lx=markers_tracks(:,index_loi(3):index_loi(3)+2);

    head = FilterData(head, 8, Fs, 30);
    head_rx = FilterData(head_rx, 8, Fs, 30);
    head_lx = FilterData(head_lx, 8, Fs, 30);

    % Ref system of the head
    med=(head_rx+head_lx)/2;
    AP=(med-head)'; % antero-posterios axis
    temp=(head_lx-head)';
    VERT=(cross(AP,temp)); % vertical axis
    ML=(cross(AP,VERT)); % medio-lateral axis: points towards the right

    %head angular position
    tilt=[(head(:,2)-med(:,2)),(med(:,1)-head(:,1))*-1];
    obl=[(head_rx(:,2)-head_lx(:,2)),(head_rx(:,3)-head_lx(:,3))*-1];
    rot=(head_rx(:,[1,3])-head_lx(:,[1,3]));
    ang_tilt=(atand(tilt(:,1)./tilt(:,2)));
    ang_rot=(atand(rot(:,1)./rot(:,2)));
    ang_obl=(atand(obl(:,1)./obl(:,2)));
    head_ang=[ang_tilt ang_rot ang_obl];

    % Computes markers' velocity using the first order forward differentiation
    h = 1; % incremental step
    vel = MarkersVelocity(head_ang, h, Fs);

    % High-pass filtering at 1Hz (Butterworth, 2nd order)
    ord = 2; fc = 1; fn = Fs/2;
    [b,a] = butter(ord,fc/fn,'high');
    vel_filt = filtfilt(b,a,vel);

    lags = linspace(-Fs,Fs,51); %2sec window

    % Head-tremor artifact removal
    sigsToReg = vel_filt';
    data = EEGCleaned.data([iChannels_LFP,iChannels],:);
    vel_embed = embedIMU(sigsToReg', lags)';
    data_embed = embedLFP(data', lags)';
    if size(vel_embed,2)>size(data_embed,2)
        vel_embed(:,end)=[];
    elseif size(data_embed,2)>size(vel_embed,2)
        data_embed(:,end)=[];
    end
    % Regress out of LFP+EEG signals multiple time-lagged copies of head
    % angular velocity
    data_cleaned = padIMU((data_embed - (data_embed/vel_embed)*vel_embed)', lags, 0)';
else
    data_cleaned = EEGCleaned.data;
end

save([pathRes,subject_ID,'_cleaned_tremor.mat'],'EEGCleaned','-v7.3');
