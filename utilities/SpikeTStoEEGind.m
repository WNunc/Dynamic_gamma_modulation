%% 找到离spike发生最近的eeg时间点
function [spikesEeg] = SpikeTStoEEGind(spikeTS, eegTS)
spikesEeg = [];
%eliminate spikes that occur earlier than the first time stamp of the EEG
for k = 1:length(spikeTS)
    % Find closest eeg timestamp to the current spike timestamp
    tDiff = (eegTS-spikeTS(k)).^2;
    [~,eegTS_ind] = min(tDiff);
    spikesEeg = [spikesEeg,eegTS_ind];
end
spikesEeg(spikesEeg<1) = 1;
end