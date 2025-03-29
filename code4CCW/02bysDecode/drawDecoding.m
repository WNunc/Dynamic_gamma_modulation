clear
clc
directories_allData_v1

for ns = 1
    path_ns = path{ns};
%     date_p=date{ns};
%     d=date_p(1:length(date_p)-1);
    for nt = 1
        Laps=6;
        phase_peak = peak(ns);
        d=0;
        for nl = 1
            dir1 = strcat(path_ns,'scores',num2str(nt),'-',num2str(nl),'.mat');
            dir_T=trackdata{ns,1};
            load(dir1);
            load(dir_T,'Ang_RewardLoc_ontrack');
            Pxn=scores{nl,3};
            xbin=scores{nl,4}; % actual bin
            mapAxis=scores{nl,5};
            tbin = scores{nl,6};
            xpos=scores{nl,9}; % actual position
            ind_pos1 = find(xpos <= Ang_RewardLoc_ontrack(1,1),1,'last');
            speed=scores{nl,8}(1,:);
            a=ind_pos1;
            b=length(xbin);
            % get theta EEG
            t_EEG_sample = [];
            t_EEG_sample = scores{nl,10};
            ind_t_EEG_sample = find(t_EEG_sample >= Ts_start_stop(nl,1)...
                & t_EEG_sample <= Ts_start_stop(nl,2));
            t_EEG_proto = t_EEG_sample;
            EEG_sample = [];
            EEG_sample = scores{nl,11}(:,ind_t_EEG_sample);
            smo = 1000; % in index points
            thetadelta = nan(size(EEG_sample,1),1);
            for ntt = 1:size(EEG_sample,1)
                TFRt = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,6,10); % theta
                TFRd = TFR_frequency_band(EEG_sample(ntt,:)',2000,5,2,4);% delta
                thetadelta(ntt,1) = mean(smooth(TFRt./TFRd,smo));
            end
            [~,ind_EEG] = max(thetadelta);
            phase=[];
            phase=scores{nl,21}(ind_EEG,:);
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
            idxp=idxp-ind_pos1;
            idxp(1,1)=idxp(1,1)+1;
            
            ind_tbin = ind_pos1:1:length(tbin);
            tbin = tbin(ind_tbin);
            t_lap0_test = tbin - tbin(1);
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
            %%
            clf
%             uimagesc(t_lap0_test,scores{nl,5},Pxn(:,ind_tbin))
            imagesc(t_lap0_test,scores{nl,5},Pxn(:,ind_tbin))
            L=length(idxp);
            y=zeros(1,L);
            y(1,:)=2*pi;
            PosX=scores{nl,9}(ind_tbin);
            fig=figure(1);
            hold on
            plot(t_lap0_test,PosX,'w','LineWidth',2) % »­ÔË¶¯¹ì¼£
            stem(t_lap0_test(idxp),y,'Marker','none','Color','m','LineWidth',1.5) % »­thetaÐòÁÐÇÐ¸îÏàÎ»
            hold off
            axis xy
            set(gca,'YTick',0:pi/4:2*pi);
            set(gca,'YTickLabel',{'0','1/4 pi','1/2 pi','3/4 pi','pi','5/4 pi','3/2 pi','7/4 pi','2pi'});
            set(gca,'fontsize',14);
            xlabel('Time / s','fontsize',14)
            ylabel('Location on track / rad','fontsize',14)
            h=colorbar;
            h.Label.String = 'Probability';
            % xlim([10,30])
            % ylim([pi,3/2*pi])
            % caxis([0 0.18])
            % ±£´æÍ¼Æ¬
            frame=getframe(fig);
            img=frame2im(frame);
            % output=strcat(num2str(d),'-',num2str(nt),'-',num2str(nl),'.png');
            d=d+1;
            output=strcat(num2str(d),'-',num2str(nt),'-',num2str(nl),'.png');
            imwrite(img,output);
        end
    end
end