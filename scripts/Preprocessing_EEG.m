% Clear workspace
clc; clear; close all;

% addCleaningPipelinePaths()

% Check if brainstorm is installed
if exist('bst_get') == 2
    disp('Brainstorm is installed and available.');
else
    disp('Brainstorm is not installed or not available');
end

% Set data directory
datadir = 'C:\Users\jasmi\Lab Isaias Dropbox\Team\WUMIGE_shared\WUMIGE_shared_projects\Dystonia_jdvdv\EEG_analysis\CleaningPipeline\EEG-preprocessing\data\';

% Load data
subj_name = 'SR72';
fprintf("Loading data...\n")
fprintf(subj_name)
[EEG,events] = loadDataFromBrainstorm(subj_name);

% Preprocessing EEG data
% EEGCleaned is an EEGLAB struct containing cleaned EEG + remaining
% recorded signals not preprocessed (e.g. LFPs, EMGs, etc.)
[EEGCleaned, EEG_steps] = MainCleaningPipeline(EEG,subj_name);


iChannels = find(strcmp({EEGCleaned.chanlocs.type},'EEG'));
iChannels_LFP = find(strcmp({EEGCleaned.chanlocs.type}, 'STN'));
iChannels_KIN = find(strcmp({EEG.chanlocs.type}, 'MISC'));
iChannel_ECG = find(strcmp({EEGCleaned.chanlocs.type}, 'ECG'));

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

% Compute PSD for cleaned LFP and EEG
% if strcmp(subj_name, 'ER48')
%     data_cleaned(3:end,:) = [];
% end

EEGCleaned.data = data_cleaned;
save([subj_name,'_data_rest_off.mat'],'EEGCleaned','EEG','EEG_steps');


subj_name = 'CR78';
datadir = 'C:\Users\jasmi\Lab Isaias Dropbox\Team\WUMIGE_shared\WUMIGE_shared_projects\Dystonia_jdvdv\EEG_analysis\CleaningPipeline\data_cleaned\fig';
power_analysis(EEG, EEGCleaned, subj_name, datadir)


% 
% epochLength = 2; % length of epochs in sec for a frequency resolution of .5 Hz
% [data_epoched] = epoch_data(data,Fs,epochLength); % epoching data
% 
% morder = 20; % model order for the spectral/connectivity analysis
% nbootstrap = 100;
% tic
% conn_uni = data2spwctrgc(data_epoched, Fs, morder, 0, nbootstrap, [], {'COH','TRGC'});
% toc %data2spwctrgc(data, fres, nlags, cond, nboot, maxfreq, output, verbose, varargin)
% iCOH_off = abs(imag(conn_uni.COH));
% TRGC_net_off = conn_uni.TRGC - permute(conn_uni.TRGC, [1 3 2 4]); % NET CONNECTIVITY i->j and j->i
% toc



psdFig = figure();
semilogy(frqs,PSDs,'LineWidth',1.5);
grid on, xlabel('Frequency [Hz]'), ylabel('10*log10(PSD) [au]');
xlim([1 45]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute PSD for raw EEG
PSD_eeg = pwelch(EEG.data(iChannels, :)', Fs, Fs/2, frqs, Fs);

% Compute PSD for cleaned EEG using the first cleaning pipeline
PSD_eeg_cleaned1 = pwelch(EEGCleaned.data(iChannels, :)', Fs, Fs/2, frqs, Fs);

% Compute PSD for cleaned EEG using the second cleaning pipeline
PSD_eeg_cleaned2 = pwelch(EEGCleaned_2.data(iChannels, :)', Fs, Fs/2, frqs, Fs);

% Normalize PSD in the range 4-80 Hz (Canessa 2020)
PSD_eeg = PSD_eeg./sum(PSD_eeg(4:80,:));
PSD_eeg_cleaned1 = PSD_eeg_cleaned1./sum(PSD_eeg_cleaned1(4:80,:));
PSD_eeg_cleaned2 = PSD_eeg_cleaned2./sum(PSD_eeg_cleaned2(4:80,:));

% Compute percentage changes in PSD
percent_change_cleaned1 = ((PSD_eeg_cleaned1 - PSD_eeg) ./ PSD_eeg) * 100;
percent_change_cleaned2 = ((PSD_eeg_cleaned2 - PSD_eeg) ./ PSD_eeg) * 100;

for ib = 1:length(bands)
    frq_inds{ib} = frqs >= bands{ib}(1) & frqs <= bands{ib}(2);
end


for iba = 1:length(bands)
    % subplot(3, 2, iba);
    h = figure;
    % Extract data for the current frequency band
    frq_ind = frq_inds{iba};
    raw_data = PSD_eeg(frq_ind, :);
    cleaned1_data = PSD_eeg_cleaned1(frq_ind, :);
    cleaned2_data = PSD_eeg_cleaned2(frq_ind, :);

    % Combine data across channels for each frequency band
    raw_data_combined = raw_data(:);
    cleaned1_data_combined = cleaned1_data(:);
    cleaned2_data_combined = cleaned2_data(:);

    % Combine data for boxplot
    combined_data = [raw_data_combined, cleaned1_data_combined, cleaned2_data_combined];

    % Generate labels for boxplot
    labels = repelem({'Raw', 'Cleaned1', 'Cleaned2'}, numel(raw_data_combined));

    % Create the boxplot
    boxplot(combined_data, 'Labels', labels);

    title(['Power in ', bandnames{iba}, ' band']);
    ylabel('Power');
    set(gca,'FontSize',22); 
    set(gcf,'windowstate','maximized');
    saveas(h, ['boxplots_comparison_',bandnames{iba},'.jpg'])
end

PSD_zapline = pwelch(EEG_steps.zapline.data', Fs, Fs/2, frqs, Fs);
PSD_rr = pwelch(EEG_steps.robust_reference.data', Fs, Fs/2, frqs, Fs);
PSD_ica1 = pwelch(EEG_steps.ica_eog_ecg.data', Fs, Fs/2, frqs, Fs);
PSD_cleanrawdata = pwelch(EEG_steps.clean_rawdata.data', Fs, Fs/2, frqs, Fs);
PSD_ica2 = pwelch(EEG_steps.ica.data', Fs, Fs/2, frqs, Fs);

PSD_eeg = PSD_eeg./sum(PSD_eeg(4:80,:));
PSD_eeg_cleaned = PSD_eeg_cleaned./sum(PSD_eeg_cleaned(4:80,:));
PSD_zapline = PSD_zapline./sum(PSD_zapline(4:80,:));
PSD_rr = PSD_rr./sum(PSD_rr(4:80,:));
PSD_ica1 = PSD_ica1./sum(PSD_ica1(4:80,:));
PSD_cleanrawdata = PSD_cleanrawdata./sum(PSD_cleanrawdata(4:80,:));
PSD_ica2 = PSD_ica2./sum(PSD_ica2(4:80,:));

% PSDs and topomaps comparison
figure; pop_spectopo(EEG, 1, [0 EEG.xmax*1000], 'EEG' , 'percent', 100, 'freq', [10 20 ], 'freqrange',[0 50],'electrodes','on');
% figure; pop_spectopo(EEG_steps.zapline, 1, [0 EEG_steps.zapline.xmax*1000], 'EEG' , 'percent', 100, 'freq', [10 20], 'freqrange',[0 50],'electrodes','on');
% figure; pop_spectopo(EEG_steps.zapline, 1, [0 EEG_steps.zapline.xmax*1000], 'EEG' , 'percent', 100, 'freq', [10 20], 'freqrange',[0 50],'electrodes','on');
% figure; pop_spectopo(FullEEGCleaned, 1, [0 FullEEGCleaned.xmax*1000], 'EEG' , 'percent', 100, 'freq', [10 20 ], 'freqrange',[0 50],'electrodes','on');



% Power per frequency band
EEG.chanlocs = EEG
Powband = [];
Powband_cleaned = [];
fig = figure();
for iba = 1:length(bands)
    Powband(iba, :) = sum(PSD_eeg(frq_inds{iba}, :));
    ax(iba) = subplot(1,4,iba);
    if iba == 1
        topoplot(Powband(iba,:),EEG.chanlocs(1:end-2),'electrodes','on','chaninfo',EEG.chaninfo,'maplimits',[-3.82,3.82]);
    else
        topoplot(Powband(iba,:),EEG.chanlocs,'electrodes','on','chaninfo',EEG.chaninfo,'maplimits',[-0.7,0.7]);
    end
    title(bandnames{iba});
    set(gca,'FontSize',18);
end

sgtitle('Raw EEG','FontSize',20)
set(fig,'WindowState', 'maximized')
saveas(gca,'rawEEG_topo.jpg')
close all;

Powband_zapline = [];
fig = figure();
for iba = 1:length(bands)
    Powband_zapline(iba, :) = sum(PSD_zapline(frq_inds{iba}, :));
    ax(iba) = subplot(1,4,iba);
    if iba == 1
        topoplot(Powband_zapline(iba,:),EEG_steps.zapline.chanlocs,'electrodes','on','chaninfo',EEG_steps.zapline.chaninfo,'maplimits',[-3.82,3.82]);
    else
        topoplot(Powband_zapline(iba,:),EEG_steps.zapline.chanlocs,'electrodes','on','chaninfo',EEG_steps.zapline.chaninfo,'maplimits',[-.7,.7]);
    end
    title(bandnames{iba});
    set(gca,'FontSize',18);
end
sgtitle('After High-pass filter and Zapline','FontSize',20)
set(fig,'WindowState', 'maximized')
saveas(gca,'zaplineEEG_topo.jpg')
close all;

Powband_rr = [];
fig = figure();
for iba = 1:length(bands)
    Powband_rr(iba, :) = sum(PSD_rr(frq_inds{iba}, :));
    ax(iba) = subplot(1,4,iba);
    if iba == 1
        topoplot(Powband_rr(iba,:),EEG_steps.robust_reference.chanlocs,'electrodes','on','chaninfo',EEG_steps.robust_reference.chaninfo,'maplimits',[-3.82,3.82]);
    else
        topoplot(Powband_rr(iba,:),EEG_steps.robust_reference.chanlocs,'electrodes','on','chaninfo',EEG_steps.robust_reference.chaninfo,'maplimits',[-.7,.7]);
    end
    title(bandnames{iba});
    set(gca,'FontSize',18);
end
sgtitle('After robust reference','FontSize',20)
set(fig,'WindowState', 'maximized')
saveas(gca,'rrEEG_topo.jpg')
close all;

Powband_ica1 = [];
fig = figure();
for iba = 1:length(bands)
    Powband_ica1(iba, :) = sum(PSD_ica1(frq_inds{iba}, :));
    ax(iba) = subplot(1,4,iba);
    if iba == 1
        topoplot(Powband_ica1(iba,:),EEG_steps.ica_eog_ecg.chanlocs,'electrodes','on','chaninfo',EEG_steps.ica_eog_ecg.chaninfo,'maplimits',[-3.82,3.82]);
    else
        topoplot(Powband_ica1(iba,:),EEG_steps.ica_eog_ecg.chanlocs,'electrodes','on','chaninfo',EEG_steps.ica_eog_ecg.chaninfo,'maplimits',[-.7,.7]);
    end
    title(bandnames{iba});
    set(gca,'FontSize',18);
end
sgtitle('After EOG removal','FontSize',20)
set(fig,'WindowState', 'maximized')
saveas(gca,'ica1EEG_topo.jpg')
close all;

Powband_crd = [];
fig = figure();
for iba = 1:length(bands)
    Powband_crd(iba, :) = sum(PSD_cleanrawdata(frq_inds{iba}, :));
    ax(iba) = subplot(1,4,iba);
    if iba == 1
        topoplot(Powband_crd(iba,:),EEG_steps.clean_rawdata.chanlocs,'electrodes','on','chaninfo',EEG_steps.clean_rawdata.chaninfo,'maplimits',[-3.82, 3.82]);
    else
        topoplot(Powband_crd(iba,:),EEG_steps.clean_rawdata.chanlocs,'electrodes','on','chaninfo',EEG_steps.clean_rawdata.chaninfo,'maplimits',[-.7,.7]);
    end
    title(bandnames{iba});
    set(gca,'FontSize',18);
end
sgtitle('After clean_rawdata','FontSize',20)
set(fig,'WindowState', 'maximized')
saveas(gca,'clean_rawdataEEG_topo.jpg')
close all;

Powband_ica2 = [];
fig = figure();
for iba = 1:length(bands)
    Powband_ica2(iba, :) = sum(PSD_ica2(frq_inds{iba}, :));
    ax(iba) = subplot(1,4,iba);
    if iba == 1
        topoplot(Powband_ica2(iba,:),EEG_steps.ica.chanlocs,'electrodes','on','chaninfo',EEG_steps.ica.chaninfo,'maplimits',[-3.82,3.82]);
    else
        topoplot(Powband_ica2(iba,:),EEG_steps.ica.chanlocs,'electrodes','on','chaninfo',EEG_steps.ica.chaninfo,'maplimits',[-.7,.7]);
    end
    title(bandnames{iba});
    set(gca,'FontSize',18);
end
sgtitle('After ICA','FontSize',20)
set(fig,'WindowState', 'maximized')
saveas(gca,'icaEEG_topo.jpg')
close all;


fig = figure();
iplot = 1;
for iba = 1:length(bands)
    Powband_cleaned(iba, :) = sum(PSD_eeg_cleaned1(frq_inds{iba}, :));
    ax(iplot) = subplot(1,5,iplot);
    if iba == 1
        topoplot(Powband_cleaned(iba,:),EEG.chanlocs,'electrodes','on','chaninfo',EEG.chaninfo,'maplimits',[-3.82, 3.82]);
    else
        topoplot(Powband_cleaned(iba,:),EEG.chanlocs,'electrodes','on','chaninfo',EEG.chaninfo,'maplimits',[-.7, .7]);
    end
    title(bandnames{iba});
    set(gca,'FontSize',18);
    iplot = iplot + 1;
end
sgtitle('Cleaned EEG','FontSize',20)
set(fig,'WindowState', 'maximized')
saveas(gca,'cleanedEEG_topo.jpg')
close all;

% Plot the frequency analysis
% Frequency Vs Power for Channel
fig = figure;
% ax(1) = subplot(121);
semilogy(frqs, PSD_eeg);
grid on;
xlim([0 80]);
xlabel('Frequency (Hz)');
ylabel('Power [log10]');
% legend({EEG.chanlocs.labels}');
title('PSDs raw');
xlim([1 80]), ylim([10^-4 10]);
set(fig,'WindowState', 'maximized')
saveas(gca,'rawEEG_psd.jpg')
close all;

fig = figure;
% ax(1) = subplot(121);
semilogy(frqs, PSD_zapline);
grid on;
xlim([0 80]);
xlabel('Frequency (Hz)');
ylabel('Power [log10]');
% legend({EEG.chanlocs.labels}');
title('PSDs after zapline');
xlim([1 80]), ylim([10^-4 10]);
set(fig,'WindowState', 'maximized')
saveas(gca,'zaplineEEG_psd.jpg')
close all;

fig = figure;
% ax(1) = subplot(121);
semilogy(frqs, PSD_rr);
grid on;
xlim([0 80]);
xlabel('Frequency (Hz)');
ylabel('Power [log10]');
% legend({EEG.chanlocs.labels}');
title('PSDs after robust reference');
xlim([1 80]), ylim([10^-4 10]);
set(fig,'WindowState', 'maximized')
saveas(gca,'rrEEG_psd.jpg')
close all;

fig = figure;
% ax(1) = subplot(121);
semilogy(frqs, PSD_ica1);
grid on;
xlim([0 80]);
xlabel('Frequency (Hz)');
ylabel('Power [log10]');
% legend({EEG.chanlocs.labels}');
title('PSDs after EOG/ECG removal');
xlim([1 80]), ylim([10^-4 10]);
set(fig,'WindowState', 'maximized')
saveas(gca,'ica1EEG_psd.jpg')
close all;

fig = figure;
% ax(1) = subplot(121);
semilogy(frqs, PSD_cleanrawdata);
grid on;
xlim([0 80]);
xlabel('Frequency [Hz]');
ylabel('PSD uV^2/Hz [dB]');
% legend({EEG.chanlocs.labels}');
title('PSDs after clean_rawdata()');
xlim([1 80]), ylim([10^-4 10]);
set(fig,'WindowState', 'maximized')
saveas(gca,'clean_rawdataEEG_psd.jpg')
close all;

fig = figure;
% ax(1) = subplot(121);
semilogy(frqs, PSD_ica2);
grid on;
xlim([0 80]);
xlabel('Frequency (Hz)');
ylabel('Power [log10]');
% legend({EEG.chanlocs.labels}');
title('PSDs after final ICA');
xlim([1 80]), ylim([10^-4 10]);
set(fig,'WindowState', 'maximized')
saveas(gca,'ica2EEG_psd.jpg')
close all;

% ax(2) = subplot(122);
fig = figure();
semilogy(frqs, PSD_eeg_cleaned);
grid on;
xlim([0 80]);
xlabel('Frequency (Hz)');
ylabel('Power [log10]');
% legend({EEG.chanlocs.labels}');
xlim([1 80]), ylim([10^-4 10]);
title('PSDs cleaned EEG');
set(gca,'FontSize',22);
set(fig,'WindowState', 'maximized')
saveas(gca,'cleanedEEG_psd.jpg')
close all;

%save plot
saveas(psdfig, ['C:\Users\jasmi\OneDrive\Desktop\EEG_analysis\Preprocessing\psd_', subj_name, '.png'])
saveas(powfig, ['C:\Users\jasmi\OneDrive\Desktop\EEG_analysis\Preprocessing\pow_topo_', subj_name, '.png'])

