function plot_NPX(NPX,mid_tbins,cycind,abc)
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明

imagesc(NPX);
hold on
stem(cycind,abc,'r','Marker','none','LineWidth',1)
stem(mid_tbins,[abc 36],'w--','Marker','none','LineWidth',1)
hold off
caxis([0.005,0.12])
colorbar
colormap jet
axis xy
yticks([1  18.5  36])
yticklabels({-63,0,63})
ylabel('distance(cm)')
end

