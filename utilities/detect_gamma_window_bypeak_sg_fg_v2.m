% function [peak_ind_unique,non_overlap_peak_ind_unique,start2,stop2,non_overlap_gamma_windows_CA1,ts_gamma_window] = detect_gamma_window(ts,EEG, f1, f2, Fs);
function [non_overlap_slow_gamma_windows_CA1,non_overlap_slow_gamma_windows_ts,non_overlap_slow_gamma_windows_TFR_z,...
    non_overlap_slow_gamma_windows_TFR_z0,non_overlap_slow_gamma_windows_bp,start2_slow,stop2_slow,...
    non_overlap_fast_gamma_windows_CA1,non_overlap_fast_gamma_windows_ts,non_overlap_fast_gamma_windows_TFR_z,...
    non_overlap_fast_gamma_windows_TFR_z0,non_overlap_fast_gamma_windows_bp,start2_fast,stop2_fast]...
    = detect_gamma_window_bypeak_sg_fg_v2(ts,EEG, sgsig, fgsig, s1,s2,f1, f2, Fs,delta_f,win_len)
%EEG = ch3_rev_begin1;
% detect both slow and fast gamma

%f1 = 100;    %  frequency band of interest
%f2 = 140;

if nargin<9
    win_len=0.2;  % default window length is 0.2s
end

[TFR_s,temp_s] = TFR_frequency_band_cz(EEG,Fs,5, s1, s2,delta_f);
[TFR_f,temp_f] = TFR_frequency_band_cz(EEG,Fs,5, f1, f2,delta_f);
% TFR should be time-length vector, each value should be mean power on each time point

TFR_f_z = zscore(TFR_f);
TFR_s_z = zscore(TFR_s);

TFR_threshold=3;  % can be change to 3

%% detect fast gamma windows
frequency='fast';
TFR_z=TFR_f_z;
TFR_z0=TFR_s_z;
temp=temp_f;
temp0=temp_s;
bp = fgsig;
% bp = fftbandpass(EEG,Fs,f1-2,f1,f2,f2+2);
high = find(TFR_z > TFR_threshold & TFR_z0 < TFR_threshold);
% high = find(TFR_z > TFR_threshold);

%bp = bandpass(EEG,f1,f2,Fs);
detect_length = round(0.16 * Fs); %   length of the window of interest (time * sampling frequency)
detect_start = round(0.08 * Fs);  %   start of the window of interest
kc = 0;
peak_ind=[];
for k = 1:size(high,2)              %   for all samples with a high power
    start = high(k) - detect_start;    %   determine start of window of interest
    stop = start + detect_length;               %   determine stop of window of interest
    if 0 < start && length(bp) >= stop        %   window has to be inside the whole EEG
        [max_val max_ind]= max(bp(start:stop));  %   determine the index of the maximal value in the bandpassed EEG
        ind0=start + max_ind - 1;
        if TFR_z(ind0)>TFR_threshold
            kc = kc + 1;                            %   counter for windows of interest inside the whole EEG
            peak_ind(kc) = start + max_ind - 1;      %   add the offset, in peak_ind are now all the sample indices with high power
        end
    end
end
if(isempty(peak_ind))
    start2_fast=[];
    stop2_fast=[];
    non_overlap_fast_gamma_windows_CA1=[];
    non_overlap_fast_gamma_windows_ts=[];
    non_overlap_fast_gamma_windows_TFR_z=[];
    non_overlap_fast_gamma_windows_TFR_z0=[];
    non_overlap_fast_gamma_windows_bp=[];
    start2_slow=[];
    stop2_slow=[];
    non_overlap_slow_gamma_windows_CA1=[];
    non_overlap_slow_gamma_windows_ts=[];
    non_overlap_slow_gamma_windows_TFR_z=[];
    non_overlap_slow_gamma_windows_TFR_z0=[];
    non_overlap_slow_gamma_windows_bp=[];
    return;
end
p = 0;
counter = 0;
peak_ind_unique=[];
for k = 1:size(peak_ind,2)
    if p < peak_ind(k)                      %   we don't want duplicates
        p = peak_ind(k);
        start = p - round((win_len/2) * Fs);      %   start of the windows that become averaged
        stop = start + round(win_len * Fs);   %   stop of the windows that become averaged
        %start = p - (round(Fs/10));      %   start of the windows that become averaged
        %stop = p + (round(Fs/10));   %   stop of the windows that become averaged
        if 0 < start & length(EEG) >= stop %   windows have to be inside the whole EEG
            counter = counter + 1;
            peak_ind_unique(counter) = p;   %   here are the indices with eliminiated duplicates
            %gamma_windows_25_140_z2(counter,:) = EEG(start:stop);
            %fast_gamma_windows_ca1(counter,:) = egf_27110711(start:stop);
            %fast_gamma_windows_ec(counter,:) = t_ec.EEG_rev_uvolts(start:stop);
            %fast_gamma_windows_ca3(counter,:) = t_ca3.EEG_rev_uvolts(start:stop);
        end
    end
end

non_overlap_peak_ind_unique = [];
for n=2:length(peak_ind_unique)
    if peak_ind_unique(n-1) + (round(Fs/10)) < peak_ind_unique(n);
        non_overlap_peak_ind_unique = [non_overlap_peak_ind_unique,peak_ind_unique(n-1)];
    else
        non_overlap_peak_ind_unique = non_overlap_peak_ind_unique;
    end
end

start2_fast=[];
stop2_fast=[];
non_overlap_fast_gamma_windows_CA1=[];
non_overlap_fast_gamma_windows_ts=[];
non_overlap_fast_gamma_windows_TFR_z=[];
non_overlap_fast_gamma_windows_TFR_z0=[];
non_overlap_fast_gamma_windows_bp=[];
ts_gamma_window=[];
p = 0;
counter = 0;
for k = 1:size(non_overlap_peak_ind_unique,2)
    if p < non_overlap_peak_ind_unique(k)                      %   we don't want duplicates
        p = non_overlap_peak_ind_unique(k);
        start = p - round((win_len/2) * Fs);
        stop = start + round(win_len * Fs);
        if 0 < start && length(EEG) >= stop
            counter = counter + 1;
            % ts_gamma_window(counter)=ts(p);
            non_overlap_fast_gamma_windows_CA1(counter,:) = EEG(start:stop);
            non_overlap_fast_gamma_windows_TFR_z(counter,:) = TFR_z(start:stop);
            non_overlap_fast_gamma_windows_TFR_z0(counter,:) = TFR_z0(start:stop);
            non_overlap_fast_gamma_windows_ts(counter,:) = ts(start:stop);
            non_overlap_fast_gamma_windows_bp(counter,:) = bp(start:stop);
            %non_overlap_gamma_windows_EC(counter,:) = EEG_EC(start:stop);
            %non_overlap_gamma_windows_CA3(counter,:) = EEG_CA3(start:stop);
            
            start2_fast(counter) = start;
            stop2_fast(counter) = stop;
            
        end
    end
end



%% detect slow gamma windows
frequency='slow';
TFR_z=TFR_s_z;
TFR_z0=TFR_f_z;
temp=temp_s;
temp0=temp_f;
bp = sgsig;
% bp = fftbandpass(EEG,Fs,s1-2,s1,s2,s2+2);
high = find(TFR_z > TFR_threshold & TFR_z0 < TFR_threshold);
% high = find(TFR_z > TFR_threshold);
%bp = bandpass(EEG,f1,f2,Fs);
detect_length = round(0.16 * Fs); %   length of the window of interest (time * sampling frequency)
detect_start = round(0.08 * Fs);  %   start of the window of interest
kc = 0;
peak_ind=[];
for k = 1:size(high,2)              %   for all samples with a high power
    start = high(k) - detect_start;    %   determine start of window of interest
    stop = start + detect_length;               %   determine stop of window of interest
    if 0 < start && length(bp) >= stop        %   window has to be inside the whole EEG
        [max_val max_ind]= max(bp(start:stop));  %   determine the index of the maximal value in the bandpassed EEG
        ind0=start + max_ind - 1;
        if TFR_z(ind0)>TFR_threshold
            kc = kc + 1;                            %   counter for windows of interest inside the whole EEG
            peak_ind(kc) = start + max_ind - 1;      %   add the offset, in peak_ind are now all the sample indices with high power
        end
    end
end
if(isempty(peak_ind))
    start2_fast=[];
    stop2_fast=[];
    non_overlap_fast_gamma_windows_CA1=[];
    non_overlap_fast_gamma_windows_ts=[];
    non_overlap_fast_gamma_windows_TFR_z=[];
    non_overlap_fast_gamma_windows_TFR_z0=[];
    non_overlap_fast_gamma_windows_bp=[];
    start2_slow=[];
    stop2_slow=[];
    non_overlap_slow_gamma_windows_CA1=[];
    non_overlap_slow_gamma_windows_ts=[];
    non_overlap_slow_gamma_windows_TFR_z=[];
    non_overlap_slow_gamma_windows_TFR_z0=[];
    non_overlap_slow_gamma_windows_bp=[];
    return;
end
p = 0;
counter = 0;
peak_ind_unique=[];
for k = 1:size(peak_ind,2)
    if p < peak_ind(k)                      %   we don't want duplicates
        p = peak_ind(k);
        start = p - round((win_len/2) * Fs);      %   start of the windows that become averaged
        stop = start + round(win_len * Fs);   %   stop of the windows that become averaged
        %start = p - (round(Fs/10));      %   start of the windows that become averaged
        %stop = p + (round(Fs/10));   %   stop of the windows that become averaged
        if 0 < start & length(EEG) >= stop %   windows have to be inside the whole EEG
            counter = counter + 1;
            peak_ind_unique(counter) = p;   %   here are the indices with eliminiated duplicates
            %gamma_windows_25_140_z2(counter,:) = EEG(start:stop);
            %slow_gamma_windows_ca1(counter,:) = egf_27110711(start:stop);
            %fast_gamma_windows_ec(counter,:) = t_ec.EEG_rev_uvolts(start:stop);
            %fast_gamma_windows_ca3(counter,:) = t_ca3.EEG_rev_uvolts(start:stop);
        end
    end
end

non_overlap_peak_ind_unique = [];
for n=2:length(peak_ind_unique)
    if peak_ind_unique(n-1) + (round(Fs/10)) < peak_ind_unique(n);
        non_overlap_peak_ind_unique = [non_overlap_peak_ind_unique,peak_ind_unique(n-1)];
    else
        non_overlap_peak_ind_unique = non_overlap_peak_ind_unique;
    end
end

start2_slow=[];
stop2_slow=[];
non_overlap_slow_gamma_windows_CA1=[];
non_overlap_slow_gamma_windows_ts=[];
non_overlap_slow_gamma_windows_TFR_z=[];
non_overlap_slow_gamma_windows_TFR_z0=[];
non_overlap_slow_gamma_windows_bp=[];
ts_gamma_window=[];
p = 0;
counter = 0;
for k = 1:size(non_overlap_peak_ind_unique,2)
    if p < non_overlap_peak_ind_unique(k)                      %   we don't want duplicates
        p = non_overlap_peak_ind_unique(k);
        start = p - round((win_len/2) * Fs);
        stop = start + round(win_len * Fs);
        if 0 < start && length(EEG) >= stop
            counter = counter + 1;
            % ts_gamma_window(counter)=ts(p);
            non_overlap_slow_gamma_windows_CA1(counter,:) = EEG(start:stop);
            non_overlap_slow_gamma_windows_TFR_z(counter,:) = TFR_z(start:stop);
            non_overlap_slow_gamma_windows_TFR_z0(counter,:) = TFR_z0(start:stop);
            non_overlap_slow_gamma_windows_ts(counter,:) = ts(start:stop);
            non_overlap_slow_gamma_windows_bp(counter,:) = bp(start:stop);
            %non_overlap_gamma_windows_EC(counter,:) = EEG_EC(start:stop);
            %non_overlap_gamma_windows_CA3(counter,:) = EEG_CA3(start:stop);
            
            start2_slow(counter) = start;
            stop2_slow(counter) = stop;
            
        end
    end
end




%{
figure;
plot(mean(non_overlap_gamma_windows_CA1));
%hold on;plot(mean(non_overlap_gamma_windows_CA3),'g');
title('num2str(f1)');
%}
%figure
%plot((1:size(non_overlap_gamma_windows_CA1,2))/Fs*1000,mean(non_overlap_gamma_windows_CA1));

%hold on;plot(mean(non_overlap_gamma_windows_CA3),'g');