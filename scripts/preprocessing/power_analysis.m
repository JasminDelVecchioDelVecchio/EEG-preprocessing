function power_analysis(EEG, EEGCleaned, subj_name, datadir)
% Input EEG: EEGlab data struct
% frequency bin indices for some EEG rhythms

% Selecting only EEG and LFP channels
iChs = find((strcmp({EEG.chanlocs.type}, 'STN') + strcmp({EEG.chanlocs.type},'EEG')) == 1);
EEG.chanlocs = EEG.chanlocs(:,iChs);
EEGCleaned.chanlocs = EEGCleaned.chanlocs(:,iChs);

% Sample rate
Fs = EEG.srate;
frqs = 1:Fs/2;
% Power spectral density via pwelch (1sec window, 50% overlap)
PSD = pwelch(EEG.data', Fs, Fs/2, frqs, Fs);
PSDCleaned = pwelch(EEGCleaned.data', Fs, Fs/2, frqs, Fs);

% Power normalization (Canessa 2020)
PSD = PSD./sum(PSD(4:80,:));
PSDCleaned = PSDCleaned./sum(PSDCleaned(4:80,:));

% Visualize EEG (FD) before and after cleaning pipeline
% Split EEG and LFP PSDs
lfp_psd = PSD(:,strcmp({EEG.chanlocs.type}, 'STN'));
eeg_psd = PSD(:,strcmp({EEG.chanlocs.type},'EEG'));

lfp_psdcleaned = PSDCleaned(:,strcmp({EEG.chanlocs.type}, 'STN'));
eeg_psdcleaned = PSDCleaned(:,strcmp({EEG.chanlocs.type},'EEG'));

% Remove 1/f component (fooof algorithm)
figure();
hraw = plot(frqs, 10*log10(mean(eeg_psd,2)),'Color',[.5 .5 .5],'LineWidth',2,'DisplayName','meanEEGraw');
plot(frqs, 10*log10(eeg_psd),'Color',[.5 .5 .5],'LineWidth',2,'LineStyle',':');
hold on; grid on; 
hcleaned = plot(frqs, 10*log10(mean(eeg_psdcleaned,2)),'Color','blue','LineWidth',2,'DisplayName','meanEEGcleaned');
plot(frqs, 10*log10(eeg_psdcleaned),'Color','blue','LineWidth',2,'LineStyle',':');
xlabel('Frequency (Hz)'); xlim([1 45]);
ylabel('Power/Frequency (dB/Hz)');
title(['PSDs EEG Comparison pre/post pre-processing for ',subj_name]);
grid on;
hold off;
set(gca,'FontSize',24);
set(gcf,'windowstate','maximized')

% Plot and compare PSDs
plot_psd_comparison(eeg_psd, lfp_psd, frqs, subj_name);
plot_power_map(eeg_psd, lfp_psd, EEG, subj_name)
plot_topomap(EEG, frqs, subj_name)

end

function plot_psd_comparison(eeg_psd, lfp_psd, frqs,subj_name)
% Plot and compare PSDs
figure;
mean_eeg_line = plot(frqs, 10*log10(mean(eeg_psd,2,'omitnan')), 'b', 'LineWidth', 2,'DisplayName','meanEEG');
hold on;
plot(frqs, 10*log10(eeg_psd), 'b', 'LineWidth', 2,'LineStyle',':');
mean_lfp_line = plot(frqs, 10*log10(mean(lfp_psd,2,'omitnan')), 'm', 'LineWidth', 2,'DisplayName','meanLFP');
plot(frqs, 10*log10(lfp_psd), 'm', 'LineWidth', 2,'LineStyle',':');
xlabel('Frequency (Hz)'); xlim([1 45]);
ylabel('Power/Frequency (dB/Hz)');
legend([mean_eeg_line, mean_lfp_line]);
title(['PSDs EEG-LFP Comparison for ',subj_name]);
grid on;
hold off;
set(gca,'FontSize',24);
set(gcf,'windowstate','maximized')
end

function plot_power_map(eeg_psd,lfp_psd,EEG,subj_name)
% Input: eeg_data - EEG data struct
bands = {[3 8], [8 12], [13 35],[13 20],[20 35]};
bandnames = {'theta', 'alpha', 'beta','lowbeta','highbeta'};
iChannels = strcmp({EEG.chanlocs.type},'EEG');
iChannels_LFP = strcmp({EEG.chanlocs.type}, 'STN');
frqs = 1:125;

% Extract EEG data for different regions
SMA_L_data = mean(eeg_psd(:,ismember({EEG.chanlocs.labels},{'F1','F3','FC1','FC3'})),2,'omitnan');
SMA_R_data = mean(eeg_psd(:,ismember({EEG.chanlocs.labels},{'F2','F4','FC2','FC4'})),2,'omitnan');
M1_L_data = mean(eeg_psd(:,ismember({EEG.chanlocs.labels},{'C1','C3','CP1','CP3'})),2,'omitnan');
M1_R_data = mean(eeg_psd(:,ismember({EEG.chanlocs.labels},{'C2','C4','CP2','CP4'})),2,'omitnan');
PC_L_data = mean(eeg_psd(:,ismember({EEG.chanlocs.labels},{'P1','P3'})),2,'omitnan');
PC_R_data = mean(eeg_psd(:,ismember({EEG.chanlocs.labels},{'P2','P4'})),2,'omitnan');
TC_L_data = mean(eeg_psd(:,ismember({EEG.chanlocs.labels},{'FT9','T7','TP9'})),2,'omitnan');
TC_R_data = mean(eeg_psd(:,ismember({EEG.chanlocs.labels},{'FT10','T8','TP10'})),2,'omitnan');

GPi_L_data = mean(lfp_psd(:,contains({EEG.chanlocs(iChannels_LFP).labels}, '_L_')),2,'omitnan');
GPi_R_data = mean(lfp_psd(:,contains({EEG.chanlocs(iChannels_LFP).labels}, '_R_')),2,'omitnan');

% Create a power heatmap
power_matrix = [SMA_L_data, SMA_R_data, M1_L_data, M1_R_data, PC_L_data, PC_R_data, TC_L_data, TC_R_data,...
    GPi_L_data, GPi_R_data];

% Sum power over each frequency band for each region
sum_power_matrix = zeros(length(bands), size(power_matrix, 2));
for ib = 1:length(bands)
    frq_inds = frqs >= bands{ib}(1) & frqs <= bands{ib}(2);
    sum_power_matrix(ib, :) = sum(power_matrix(frq_inds, :));
end

% Normalize the power matrix over frequency bands
normalized_power_matrix = sum_power_matrix ./ sum(sum_power_matrix, 2);
numColors = 256;
% blue_white_red_map = zeros(numColors, 3);
%
% % Blue to white (first half)
% blue_white_red_map(1:numColors/2, 1) = linspace(0, 1, numColors/2); % R
% blue_white_red_map(1:numColors/2, 2) = linspace(0, 1, numColors/2); % G
% blue_white_red_map(1:numColors/2, 3) = 1; % B
%
% % White (middle)
% blue_white_red_map(numColors/2+1, :) = [1 1 1];
%
% % White to red (second half)
% blue_white_red_map(numColors/2+2:end, 1) = 1; % R
% blue_white_red_map(numColors/2+2:end, 2) = linspace(1, 0, numColors/2-1); % G
% blue_white_red_map(numColors/2+2:end, 3) = linspace(1, 0, numColors/2-1); % B
% Set the maximum value color
min_value_color = [1, 1, 1]; % White color
max_value_color = [0 0 1];

% Create a custom colormap transitioning from the specified color to blue or red
custom_colormap = [linspace(min_value_color(1), max_value_color(1), numColors)', ...
    linspace(min_value_color(2), max_value_color(2), numColors)', ...
    linspace(min_value_color(3), max_value_color(3), numColors)'];

% Plot the normalized matrix using imagesc
figure;
imagesc(normalized_power_matrix);
colormap(custom_colormap)
% colormap(blue_white_red_map);
title(['Normalized Power Matrix for ',subj_name]);
xlabel('Cortical and Subcortical Regions');
ylabel('Frequency Bands');

% Add color bar
c = colorbar;
c.Label.String = 'Relative Power';

% Adjust axis ticks and labels
xticks(1:size(normalized_power_matrix, 2));
xticklabels({'SMA_L','SMA_R','M1_L','M1_R','PC_L','PC_R','TC_L','TC_R','GPi_L','GPi_R'});
yticks(1:length(bands));
yticklabels(bandnames);

% Set axis labels font size
set(gca, 'FontSize', 24);
set(gcf,'windowstate','maximized')
end

function plot_topomap(eeg_psd,frqs,EEG)
    bands = {[3 8], [8 12], [13 35],[13 20],[20 35],[35 80]};
    bandnames = {'theta', 'alpha', 'beta','lowbeta','highbeta','gamma'};
    fig = figure();
    for ib = 1:length(bands)
        frq_inds = frqs >= bands{ib}(1) & frqs <= bands{ib}(2);
        eeg_pow(iba, :) = sum(eeg_psd(frq_inds, :));
        ax(iba) = subplot(1,4,iba);
        if iba == 1
            topoplot(eeg_pow(iba,:),EEG.chanlocs,'electrodes','on','chaninfo',EEG.chaninfo,'maplimits',[-3.82,3.82]);
        else
            topoplot(eeg_pow(iba,:),EEG.chanlocs,'electrodes','on','chaninfo',EEG.chaninfo,'maplimits',[-0.7,0.7]);
        end
        title(bandnames{iba});
        set(gca,'FontSize',18);
    end
    linkaxes(ax,'x')
    sgtitle('Raw EEG','FontSize',20)
    set(fig,'WindowState', 'maximized')
    saveas(gca,'rawEEG_topo.jpg')
    close all;
end
