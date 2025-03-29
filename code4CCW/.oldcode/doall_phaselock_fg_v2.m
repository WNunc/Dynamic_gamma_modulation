% To get phase-lock value of fast gamma
clear
clc
directories_allData_v1

%%
for ns = 1:isession
    path_ns = path{ns};
    date_p = date{ns};
    disp(date_p);
    dir_C=strcat(path_ns,'Cells_singleLap_v2_vel_5.mat');
    TTlist=strcat(path_ns,'TTList_dCA1_pyr.txt');
    F = ReadFileList(TTlist);
    load(dir_C,'spikes','fieldProp_singlelap','mapAxis');
    B=[];
    for nt = 1:3
        NC=length(spikes);
        switch nt
            case 1
                laps=6;
            case 2
                laps=8;
            case 3
                laps=8;
            case 4
                laps=6;
        end
        for nl = 1:laps
            ts_sp=[];
            old_name=0;
            for nc = 1:NC
                ts = spikes{nc,1};
                file_input = strcat(date_p,'scores',num2str(nt),'-',num2str(nl),'.mat');
                load(file_input,'Ts_start_stop')
                ind_ts = find(ts >= Ts_start_stop(nl,1) & ts <= Ts_start_stop(nl,2));
                ts = ts(ind_ts);
                ts_sp{nc,1} = ts;
                % find eeg file
                if F{nc,1}(4) ~= '_';
                    eeg_name = str2double(F{nc,1}(3:4));
                elseif F{nc,1}(4) == '_';
                    eeg_name = str2double(F{nc,1}(3));
                else
                    disp('problem with something');
                end
                if eeg_name~=old_name;
                    filnam = strcat(path_ns,'\CSC',num2str(eeg_name),'.ncs');
                    % read EEG
                    [sample,tt,tt_raw] = loadCSC_new_cz(filnam); % sample unit:uv, tt and tt_raw unit: s
                    f1 = 65;
                    f2 = 100;
                    bpry = fftbandpass(sample,2000, f1-2,f1,f2,f2+2); % fast gamma
                    phsry = angle(hilbert(-bpry(1,:)))*180/pi+180;
                    old_name = eeg_name;
                end
                % find the spike phase
                spikecells{nc,1} = phsry(SpikeTStoEegIndex(ts_sp{nc,1},tt_raw,2000));
                if ~isnan(spikecells{nc,1})
                    [meanAngle,angDev,vectorLength,rayleighP,rayleighZ,cirDev] = circularStat(spikecells{nc,1});
                    result(nc,:) = [Ind_Rat(ns),nt,nl,nc,size(spikecells{nc,1},2),rayleighP,vectorLength,meanAngle];
                else
                    result(nc,:) = [Ind_Rat(ns),nt,nl,nc,0,0,0,0];
                end
%                 % find the spike phase in case you need downsample
%                 spikecells{nc,file_id} = phsry(SpikeTStoEegIndex(spikes{nc,file_id},tt_raw,2000));
%                 spikecells_ds{nc,file_id} = spikecells{nc,file_id}((randperm(numel(spikecells{nc,file_id}),spkmin(j,2))));
%                 if ~isnan(spikecells_ds{nc,file_id})
%                     [meanAngle,angDev,vectorLength,rayleighP,rayleighZ,cirDev] = circularStat(spikecells_ds{nc,file_id});
%                     result(j,:) = [Ind_Rat(ns),nt,nl,nc,size(spikecells_ds{nc,file_id},2),rayleighP,vectorLength,meanAngle];
%                 else
%                     result(j,:) = [Ind_Rat(ns),nt,nl,nc,0,0,0,0];
%                 end
            end
            B{nt,nl} = result;
        end
    end
    file_output=strcat(date_p,'spk_phs_fg.mat');
    save(file_output,'B')
    clc
end