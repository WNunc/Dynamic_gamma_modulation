function [r2,calphase,xaxis,p,slope,xspan,tspan] = imageseq_slope(tt,AA,loc)
numbins = 90; % Number of bins
bin_ang = 2*pi/numbins;
mapAxis = bin_ang/2:bin_ang:(2*pi-bin_ang/2);
[~, n] = size(AA);
tt0 = 0:(tt(2)-tt(1))/(n-1):tt(2)-tt(1);
bins2use = find(~isnan(sum(AA)));
if length(bins2use)/size(AA,2)<0.1
    r2 = [];
    calphase = [];
    xaxis = [];
    p = [];
    slope = [];
    xspan = [];
    tspan = [];
    return
end

[r2,calphase,xaxis,p,slope,xspan,tspan] = Cir_reg(AA,mapAxis,tt0,bins2use);

% sequence fit
xFit = tt0;
yFit = calphase(1)+slope*(xFit-tt0(1));

imagesc(tt0,mapAxis,AA)
hold on
plot(tt0,mapAxis(loc),'w--','LineWidth',1.5)
plot(xFit,yFit,'m','LineWidth',2.5)
hold off
% 获取图像坐标范围
slope = round(slope,3);
xrange = xlim;
yrange = ylim;
% 计算文本框的位置
xpos = xrange(2) - 0.1 * diff(xrange);
ypos = yrange(1) + 0.05 * diff(yrange);
text(xpos, ypos, [ 'Slope = ' num2str(slope)] , ...
    'Color', 'w',...
    'HorizontalAlignment', 'right',...
    'VerticalAlignment', 'bottom', 'FontSize', 15);

xticks([tt0(1),tt0(n)])
yticks([mapAxis(1),mapAxis(numbins)])
yticklabels({'0','2π'})
axis xy
axis square
axis tight
caxis([0,0.15])
set(gca,'FontSize',12)
colormap(parula(128));
end