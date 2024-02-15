function [eventSamples]=detect_qrs(EMG_ekg,on,off,EMGfs)

EMG_ekg = EMG_ekg(on:off);

flip=0;
if abs(min(EMG_ekg))>abs(max(EMG_ekg))
    EMG_ekg=EMG_ekg*(-1);
    flip=1;
end

[pks,locs] = findpeaks(zscore(EMG_ekg),'MinPeakHeight',3,'MinPeakDistance',20);
if flip==0 && length(pks)<((off-on)/EMGfs)
    EMG_ekg=EMG_ekg*(-1);
    [pks,locs] = findpeaks(zscore(EMG_ekg),'MinPeakHeight',3,'MinPeakDistance',20);
end

% figure;
% plot(zscore(EMG_ekg))
% hold on
% plot(locs,pks,'*r')
% xlabel('samples')
% ylabel('Voltage [mV]')
% title('QRS EKG EMG')
eventSamples = locs;
