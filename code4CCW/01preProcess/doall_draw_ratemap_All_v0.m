% edited on 2021/7/9
% only for CW experiments
clc
clear

directories_allData_v1
figure('Units','normalized','Position',[0 0.2 1 0.5]);
%%
for ns = 1:isession
    path_ns = path{ns};
    date_p=date{ns};
    d=date_p(1:length(date_p)-1);
    trackdata_n = trackdata{ns};
    file_input1 = strcat(path_ns,'Cells_ALLLaps_v2_vel_5.mat');
      
    % Load ratemaps of "ALLLaps" from the input1 file
    % Sort the cells by using COM from ratemaps of "ALLLaps"
    load(file_input1)
    Ncell = size(Ratemap_AllLaps,2);
    COM_allcell = nan(Ncell,1);
    peakFR_allcell = nan(Ncell,1);
    for nc = 1:Ncell
        if ~isempty(fieldProp_AllLaps{nc})
            COM_allcell(nc,1) = fieldProp_AllLaps{nc}(1,1).x_COM;  % x_COM of the 1st place field
            peakFR_allcell(nc,1) = fieldProp_AllLaps{nc}(1,1).peakRate;  % x_COM of the 1st place field
        end
    end
    [COM_cell_sort,ind_cell_sort] = sort(COM_allcell);
    ind = find(~isnan(COM_cell_sort));
    ind_cell_sort = ind_cell_sort(ind);
    Ncell_sort = length(ind);
    
    %% plot the ratemaps for CW laps
    ratemap0 = Ratemap_AllLaps;
    peakrate = repmat(peakFR_allcell',size(ratemap0,1),1);
    ratemap0 = ratemap0./peakrate;
    Ratemap_norm = ratemap0;
    clf
    fig = figure(1);
    hold on
    ratemap0_sort = ratemap0(:,ind_cell_sort);
    imagesc(mapAxis,1:Ncell_sort,ratemap0_sort');
    %     ang_reward = Ang_RewardLoc_ontrack(Ind_rewardloc);
    %     line([ang_reward,ang_reward],[1,Ncell_sort],'color','w')
    hold off
    xlim([0,2*pi]);
    set(gca,'XTick',0:pi:2*pi);
    set(gca,'XTickLabel',{'0','pi','2pi'});
    ylim([1,Ncell_sort]);
    set(gca,'YTick',[1,Ncell_sort]);
    set(gca,'fontsize',20);
    xlabel('Angle on the track (rad)')
    ylabel('Cell ID')
    title('All Laps')
    
    h=colorbar;
    h.Label.String = 'Normalized Firing Rate';
    
    frame=getframe(fig);
    img=frame2im(frame);
    output=strcat(num2str(d),'-RatemapAll.png');
    imwrite(img,output);
    
end