function  plot_TSS(TSS)
% theta_seqence_averaging 专用画图函数
%   此处显示详细说明
imagesc(TSS);

x=[27 42 42 27];%四个顶点
y=[10 10 27 27];%四个顶点
X=[x x(1)];%首尾相接
Y=[y y(1)];
xl = zeros(1,68)+18.5;
yl = zeros(1,36)+34.5;
hold on
plot(X,Y,'w','LineWidth',1.5);
plot(1:68,xl,'w','LineWidth',1.5)
plot(yl,1:36,'w','LineWidth',1.5)
hold off

caxis([0.005,0.12])
colormap jet
axis xy
yticks([1  18.5  36])
yticklabels({-63,0,63})
ylabel('distance(cm)')
xticks([1 34.5 68])
xticklabels({-170, 0, 170})
xlabel('time(ms)')
set(gca,'FontSize',12)
end

