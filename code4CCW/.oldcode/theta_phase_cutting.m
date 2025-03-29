%% cut theta phase on maxspkdensity

nseg = 1;%哪个部分？ 1=pre-running, 2=sample, 3=test, 4=post
% nl = 2;%哪一圈？
for nl = 1:5
    lap_ts = scores{nseg}{nl,6};%行为学时间
    lap_ang = scores{nseg}{nl,5};%位置（角度）
    lap_pxn = scores{nseg}{nl,3};%解码位置（Pxn）
    lap_pos = scores{nseg}{nl,5}(scores{nseg}{nl,4});%实际位置
    lap_EEG = scores{nseg}{nl,11};%EEG数据
    lap_EEG_ts = scores{nseg}{nl,10};%EEG时间
    lap_theta = scores{nseg}{nl,18};%scores中的theta
    lap_thetaphs = scores{nseg}{nl,21};%scores中的theta相位
    lap_timepoint_track = size(lap_pxn,2);%路径总时间点
    lap_timepoint_EEG = size(lap_EEG,2);%EEG总时间点
    lap_speed = scores{nseg}{nl,8}(1,:);
    %% 运动时的时间
    
    limit = TS_running(nl,:);
    
    %% 用最大spike密度处的相位切割theta
    % 对应的时间点
    
    % [~,I] = min(abs(lap_thetaphs(4,:)-max_phase),[],1);
    [~,I]=findpeaks(-abs(lap_thetaphs(mth,:)-max_phase),'MinPeakHeight',-5);
    
    % 对应的位置
    IND = [];
    for i= 1:length(I)
        [~,IND(i)] = min(abs(lap_ts - lap_EEG_ts(I(i))));
    end
    
    thetacycle.ind{nl} = IND;% 全部的index
    thetacycle.ts{nl} = lap_ts(IND); % 全部的时间点
    thetacycle.EEG_ts{nl} = lap_EEG_ts(I); % 全部的EEG时间点
    for i = 1:length(I)-1
        thetacycle.cycle_ind{nl}(i,:) = [IND(i),IND(i+1)]; %每个theta周期开始和结束的index
        thetacycle.cycle_ts{nl}(i,:) = [lap_ts(IND(i)),lap_ts(IND(i+1))];%每个theta周期开始和结束的时间点
    end
    
    figure
    set(gcf,'Position',[2 42 1438 952]);
    subplot(2,1,1)
    imagesc(lap_ts,lap_ang,lap_pxn)
    hold on
    plot(lap_ts,lap_pos,'--w','LineWidth',1.5) % 画运动轨迹
    stem(lap_ts(IND),repmat( [lap_ang(90)],1,length(IND)),'w','LineWidth',1,'Marker','none')
    hold off
    colormap jet
    colorbar('Position',[0.91 0.7093,0.015,0.21])
    caxis([0,0.25])
    axis xy
    xlim(limit)
    
    title('decoding')
    subplot(2,1,2)
    plot(lap_EEG_ts,lap_theta(mth,:)','k')
    hold on
    scatter(lap_EEG_ts(I),lap_theta(mth,I)')
    hold off
    xlim(limit)
    saveas(gcf,[subfolder 'decoding_maxspk_phase_cutting_lap' num2str(nl) '_' num2str(nbars),'bars-min.png'])
    saveas(gcf,[subfolder 'decoding_maxspk_phase_cutting_lap' num2str(nl) '_' num2str(nbars),'bars-min.fig'])
end
save([subfolder, 'data_AllLap_thetaphase_',num2str(nbars),'bars-min.mat'],'thetacycle','-append');