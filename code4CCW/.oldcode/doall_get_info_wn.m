% Get info.mat
% updated on 2021/7/9
% clc
% clear
% directories_allData_v1

% 临时组会用

%%
for ns = 1:isession
    path_ns = path{ns};
    %     date_s = date{ns};
    phase_peak = max_phase;
    %     cd(path_ns);
    disp(path_ns)
    for nseg = 1%:3
        switch nseg
            case 1
                laps = 6;
            case 2
                laps = 8;
            case 3
                laps = 8;
        end
        for nl = 1:laps
            file_input1=strcat(path_ns,'scores',num2str(nseg),'-',num2str(nl),'-wn.mat');
            load(file_input1);
            scores_sl = scores;
            load([path_ns 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v2.mat']);
            %             file_input2=strcat(path_ns,'Cells_singleLap_v2_vel_5.mat');
            %             load(file_input2,'spikes');
            Pxn=[];
            xbin=[];
            tbin=[];
            speed=[];
            Pxn=scores_sl{nl,3};
            xbin=scores_sl{nl,4}; % actual bin
            tbin = scores_sl{nl,6};
            speed=scores_sl{nl,8}(1,:);
            a=find(speed>20,1);
            if isempty(a)
                a=1;
            end
            b=length(xbin)-1;
            % get theta EEG
            t_lap_sample = [];
            t_lap_sample = scores_sl{nl,6};
            ind_t_lap_sample = find(t_lap_sample >= Ts_start_stop(nl,1)...
                & t_lap_sample <= Ts_start_stop(nl,2));
            t_lap_sample = t_lap_sample(ind_t_lap_sample);
            t_EEG_sample = [];
            t_EEG_sample = scores{nseg}{nl,10};
            ind_t_EEG_sample = find(t_EEG_sample >= Ts_start_stop(nl,1)...
                & t_EEG_sample <= Ts_start_stop(nl,2));
            t_EEG_sample = t_EEG_sample(ind_t_EEG_sample);
            t_EEG_proto = t_EEG_sample;
            t_EEG_sample = t_EEG_sample - t_lap_sample(1);
            EEG_sample = [];
            EEG_sample = scores{nseg}{nl,11}(:,ind_t_EEG_sample);
            smo = 1000; % in index points
            thetadelta = nan(size(EEG_sample,1),1);
            for ntt = 1:size(EEG_sample,1)
                TFRt = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,6,10); % theta
                TFRd = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,2,4);% delta
                thetadelta(ntt,1) = mean(smooth(TFRt./TFRd,smo));
            end
            [~,ind_EEG] = max(thetadelta);
            phase=[];
            phase=scores{nseg}{nl,21}(ind_EEG,:);
            Difa=abs(t_EEG_proto-tbin(a));
            [~,Ia]=min(Difa);
            Difb=abs(t_EEG_proto-tbin(b));
            [~,Ib]=min(Difb);
            t_phase=t_EEG_proto(:,Ia:Ib);
            phase=phase(:,Ia:Ib);
            % find index of phase_peak-idxp
            ind_peak=[];
            idxp=[];
            dis_p=phase-phase_peak;
            c=1;
            for j=2:length(dis_p)-1
                if dis_p(j-1)<0 & dis_p(j+1)>0
                    ind_peak(1,c)=j;
                    c=c+1;
                end
            end
            tt=t_phase(ind_peak);
            for i=1:length(ind_peak)
                Dif=abs(tbin-tt(1,i));
                [~,idxp(1,i)]=min(Dif);
            end
            idxp=unique(idxp);
            DifI=diff(idxp);
            while find(DifI<=7)
                for j=1:length(idxp)-1
                    if DifI(j)<=7
                        idxp(j)=idxp(j+1);
                    end
                end
                idxp=unique(idxp);
                DifI=diff(idxp);
            end
            %
            t_thc=tbin(idxp);
            mid=[];
            for m = 1:length(t_thc)-1
                % spikes and active cells in single lap
                for n=1:length(spikes)
                    ind_sp=find(spikes{n,1}>=t_thc(1,m) & spikes{n,1}<=t_thc(1,m+1));
                    sp(n,m)=length(ind_sp);
                end
                spk(1,m)=sum(sp(:,m));
                ind_ac=find(sp(:,m)~=0);
                ac(1,m)=length(ind_ac);
                mid(1,m)=(idxp(m)+idxp(m+1))/2;
         centemp   end
            % get part of the Pxn
            mid=round(mid);
            for m = 1:length(t_thc)-1
                if mid(1,m)+20>=b
                    ind_th(m,1)=mid(1,m)-20;
                    ind_th(m,2)=b;
                elseif mid(1,m)-20<=1
                    ind_th(m,1)=1;
                    ind_th(m,2)=mid(1,m)+20;
                else
                    ind_th(m,1)=mid(1,m)-20;
                    ind_th(m,2)=mid(1,m)+20;
                end
            end
            pos_mid=xbin(mid);
            for m = 1:length(t_thc)-1
                if pos_mid(1,m)<=6
                    pos(m,1)=1;
                else
                    pos(m,1)=pos_mid(1,m)-6;
                end
                if pos_mid(1,m)+6>=90
                    pos(m,2)=90;
                else
                    pos(m,2)=pos_mid(1,m)+6;
                end
            end
            spd=[];
            for m = 1:length(t_thc)-1
                spd(m,1)=mean(speed(idxp(m):idxp(m+1)));
            end
            
            
            % power and its zscore
            for m = 1:length(t_thc)-1
                ind_s=[idxp(m),idxp(m+1)];
                for n = 1:2
                    dis_t=abs(t_EEG_proto-tbin(ind_s(1,n)));
                    [~,ind_tc(m,n)]=min(dis_t);
                end
            end
            TFRsg=[];
            TFRfg=[];
            for ntt = 1:size(EEG_sample,1)
                TFRsg(ntt,:) = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,25,55); % slow gamma
                TFRfg(ntt,:) = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,65,100); % fast gamma
            end
            Zsg=zscore(TFRsg,0,2); % zscore between tetrodes
            Zfg=zscore(TFRfg,0,2);
            for m = 1:length(t_thc)-1
                for ntt = 1:size(EEG_sample,1)
                    Zsg_s(ntt,:)=Zsg(ntt,ind_tc(m,1):ind_tc(m,2));
                    Zfg_s(ntt,:)=Zfg(ntt,ind_tc(m,1):ind_tc(m,2));
                end
                Psg=TFRsg(:,ind_tc(m,1):ind_tc(m,2)
                Zsg_m=mean(Zsg_s(:)););
                Pfg=TFRfg(:,ind_tc(m,1):ind_tc(m,2));
                PowSG(m,1)=mean(Psg(:));
                PowFG(m,1)=mean(Pfg(:));
                Zfg_m=mean(Zfg_s(:));
                Mzsg(m,1)=Zsg_m;
                Mzfg(m,1)=Zfg_m;
                Zsg_s=[];
                Zsg_m=[];
                Zfg_s=[];
                Zfg_m=[];
            end
            % output the information
            info=[];
            for m = 1:length(t_thc)-1
                info{m,1}=[t_thc(1,m),t_thc(1,m+1)]; % timestamp of thetacycle
                info{m,2}=[idxp(m),idxp(m+1)]; % index of thetacycle
                info{m,3}=Pxn(pos(m,1):pos(m,2),ind_th(m,1):ind_th(m,2)); % Pxn of thetacycle
                info{m,4}=ac(1,m); % active cells
                info{m,5}=spk(1,m); % spikes
                info{m,6}=NXP(:,ind_th(m,1):ind_th(m,2));
                info{m,7}=PowFG(m,1); % power of fast gamma
                info{m,8}=Mzfg(m,1); % normalized fast gamma
                info{m,9}=PowSG(m,1); % power of slow gamma
                info{m,10}=Mzsg(m,1); % normalized of slow gamma
                info{m,11}=spd(m,1); % mean speed of each cycle
            end
            file_output=strcat(path_ns,'info',num2str(nseg),'-',num2str(nl),'.mat');
            save(file_output,'info')
        end
    end
end