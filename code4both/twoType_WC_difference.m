% scatter plot of fgcell ratio with WC difference between the two type
% decoded sequence
% after running the doall_fgcell_percentage.m

cd('H:\neuralynx\gamma in sequence result\theta sequence indicator\0829')
load('WC-f2lapup.mat')
% 所有session所有sequence的WC的平均
X1 = GPdataB1_nfgcell;
X2 = GPdataB1_fgcell;

% cd('H:\neuralynx\gamma in sequence result\theta sequence indicator\230214')
% load('WC-f2lapup.mat')
% % 所有session的平均sequence的WC
% X1 = GPdataB1_nfgcell0;
% X2 = GPdataB1_fgcell0;

close all
X = [ones(length(ratio),1) ratio'*100];
diff_fgvsnfg = mean(X1(1:5,:)-X2(1:5,:));
y = diff_fgvsnfg';
% linear regress
[b,bint,r,rint,stats] = regress(y,X);
ffa = figure('unit','centimeters','position',[20 10 12 10]);
scatter(ratio*100,diff_fgvsnfg,80,'LineWidth',2)
hold on
xFit = min(ratio'*100)-3:1:max(ratio'*100)+5;
yFit = b(1) + b(2)*xFit;
plot(xFit,yFit,'k--','LineWidth',1.5)
hold off
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Fraction of excluded cell /%','FontSize',18)
ylabel({'Weight correlation difference ','(nfgcell - fgcell)'},'FontSize',18)
ylim([-0.05 0.15])
ylim([-0.2 0.2])
title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',18)
saveas(ffa,'all_lap_diff.png')
saveas(ffa,'all_lap_diff','epsc')
saveas(ffa,'all_lap_diff.fig')

diff_fgvsnfg = mean(X1(1:2,:)-X2(1:2,:),1);
y = diff_fgvsnfg()';
% linear regress
[b,bint,r,rint,stats] = regress(y,X);
ffa = figure('unit','centimeters','position',[20 10 12 10]);
scatter(ratio*100,diff_fgvsnfg,80,'LineWidth',2)
hold on
xFit = min(ratio'*100)-3:1:max(ratio'*100)+5;
yFit = b(1) + b(2)*xFit;
plot(xFit,yFit,'k--','LineWidth',1.5)
hold off
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Fraction of excluded cell /%','FontSize',18)
ylabel({'Weight correlation difference ','(nfgcell - fgcell)'},'FontSize',18)
ylim([-0.12 0.12])
ylim([-0.2 0.2])
title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',18)
saveas(ffa,'f2_lap_diff.png')
saveas(ffa,'f2_lap_diff.eps','psc2')
saveas(ffa,'f2_lap_diff.fig')

diff_fgvsnfg = mean(X1(4:5,:)-X2(4:5,:));
y = diff_fgvsnfg()';
% linear regress
[b,bint,r,rint,stats] = regress(y,X);
ffa = figure('unit','centimeters','position',[20 10 12 10]);
scatter(ratio*100,diff_fgvsnfg,80,'LineWidth',2)
hold on
xFit = min(ratio'*100)-3:1:max(ratio'*100)+5;
yFit = b(1) + b(2)*xFit;
plot(xFit,yFit,'k--','LineWidth',1.5)
hold off
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Fraction of excluded cell /%','FontSize',18)
ylabel({'Weight correlation difference ','(nfgcell - fgcell)'},'FontSize',18)
ylim([-0.05 0.15])
ylim([-0.2 0.2])
title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',18)
saveas(ffa,'l2_lap_diff.png')
saveas(ffa,'l2_lap_diff.eps','psc2')
saveas(ffa,'l2_lap_diff.fig')