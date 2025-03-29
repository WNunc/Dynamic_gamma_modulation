function [p1,P1s,b,bint,r,rint,stats] = scatimg_XY(x,y,xbin,ybin,xl,yl,flag)
% 画图散点图和拟合直线以及概率热图
% xbin = -pi:pi/20:pi;
% ybin = -1:1/20:1;

X = [ones(length(x),1),x];

% linear regress
[b,bint,r,rint,stats] = regress(y,X);
% 平滑后的概率
p1=hist2d(x,y,xbin,ybin);
P1=p1/sum(sum(p1));
P1s = smooth2a(normz(P1,0),3,3);

if flag==1
    %% 画散点图和拟合直线(可以注释掉)
    subplot(1,4,1)
    scatter(x,y)
    xFit = -pi:0.1*pi:pi;
    yFit = b(1) + b(2)*xFit;
    hold on
    plot(xFit,yFit,'k--','LineWidth',1.5)
    hold off
    title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',11)
    axis xy; axis square;
    xlabel(xl,'Fontsize',11)
    set(gca,'XTick',-pi:pi:pi);
    set(gca,'XTickLabel',{'-π','0','π'},'Fontsize',11)
    ylabel(yl,'Fontsize',11)
    ylim([-1,1])
    %% 画平滑后的概率热图(可以注释掉)
    subplot(1,4,2)
    imagesc(xbin,ybin,P1s)
    axis xy; axis square;
    xlabel(xl,'Fontsize',11)
    set(gca,'XTick',-pi:pi:pi);
    set(gca,'XTickLabel',{'-π','0','π'},'Fontsize',11)
    ylabel(yl,'Fontsize',11)
    colormap(jet(128));
end
end