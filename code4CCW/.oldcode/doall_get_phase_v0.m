% Get phase of error peaks PHASE
clc
clear
directories_allData_v1
%%
for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
    disp(pwd)
    date_p=date{ns};
    for nt = 1:3
        switch nt
            case 1
                laps = 6;
            case 2
                laps = 8;
            case 3
                laps = 8;
        end
        PHASE=[];
        for nl = 1:laps
            dir1=strcat(path_ns,'scores',num2str(nt),'-',num2str(nl),'.mat');
            dir_T=trackdata{ns,1};
            load(dir1);
            load(dir_T,'Ang_RewardLoc_ontrack');
            Pxn=scores{nl,3};
            xbin=scores{nl,4}; % actual bin
            mapAxis=scores{nl,5};
            tbin = scores{nl,6};
            xpos=scores{nl,9}; % actual position
            [max_a,index]=max(Pxn);
            ind_pos1 = find(xpos <= Ang_RewardLoc_ontrack(1,1),1,'last');
            a=ind_pos1;
            b=length(xbin)-1;
            %find error peak
            xdc=mapAxis(index);
            error=xpos-xdc';
            [yp,idxpE]=findpeaks(abs(error(a:b)),'MinPeakHeight',0.1);% different from chronux
            L=length(idxpE);
            idxpE=idxpE+a;
            % get theta EEG in sample lap
            t_lap_sample = scores{nl,6};
            ind_t_lap_sample = find(t_lap_sample >= Ts_start_stop(nl,1)...
                & t_lap_sample <= Ts_start_stop(nl,2));
            t_lap_sample = t_lap_sample(ind_t_lap_sample);
            %
            t_EEG_sample = scores{nl,10};
            ind_t_EEG_sample = find(t_EEG_sample >= Ts_start_stop(nl,1)...
                & t_EEG_sample <= Ts_start_stop(nl,2));
            t_EEG_sample = t_EEG_sample(ind_t_EEG_sample);
            t_EEG_proto = t_EEG_sample;
            t_EEG_sample = t_EEG_sample - t_lap_sample(1);
            EEG_sample = scores{nl,11}(:,ind_t_EEG_sample);
            smo = 1000; %in index points
            thetadelta = nan(size(EEG_sample,1),1);
            for ntt = 1:size(EEG_sample,1)
                TFRt = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,6,10); %theta
                TFRd = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,2,4);%delta
                thetadelta(ntt,1) = mean(smooth(TFRt./TFRd,smo));
            end
            [~,ind_EEG] = max(thetadelta);
            phase=scores{nl,21}(ind_EEG,:);
            Difa=abs(t_EEG_proto-tbin(a));
            [~,Ia]=min(Difa);
            Difb=abs(t_EEG_proto-tbin(b));
            [~,Ib]=min(Difb);
            t_phase=t_EEG_proto(:,Ia:Ib);
            phase=phase(:,Ia:Ib);
            %
            ind_peakE=zeros(1,L);
            for s=1:L
                dis_t=abs(t_phase-tbin(idxpE(1,s)));
                [~,ind_peakE(1,s)]=min(dis_t);
            end
            phase_peaks=phase(ind_peakE);
            L=length(phase_peaks);
            for i=1:L
                PHASE(nl,i)=phase_peaks(1,i);
            end
        end
        dir2=strcat(path_ns,'phase',num2str(nt),'.mat');
        save(dir2,'PHASE');
    end
end

% cd 'F:/MATLAB/MATLAB Production Server/R2015a/bin'