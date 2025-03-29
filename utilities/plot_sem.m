function plot_sem(data,xLabel,yLabel,color)
%plot_sem 画数据的折线图，error是sem
%   lineplot with error bar
if nargin < 4 && nargin~= 2
    color = 'flat';
end
if size(data,2)>1
    warning('sem计算可能有误');
end
data_plot = mean(data,'omitnan');
data_sd = std(data,'omitnan');
data_num = size(data(~isnan(data)),1);
data_sem = data_sd/sqrt(data_num);
errorbar(data_plot,data_sem,'LineWidth',2,'Color',color)
% xlim([0.5,5.5])
% ylim([min(data_plot)*0.8,max(data_plot)*1.2])
if nargin >= 4
    xlabel(xLabel)
    ylabel(yLabel)
end

set(gca,'FontSize',15,'Color','none','LineWidth',1)
box off

