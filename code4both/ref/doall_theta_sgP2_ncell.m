%Part2
% 按cell
clear
close all
clc
directories_allData_v0_allgood
outputFolder = ['H:\neuralynx\phase precession rseult\'];
%%
CNEG = []; CZERO = []; CPOS = [];
nt = 1;
for ns = [1,2,4:6,8,10:13,15,16]%1:isession
    path_ns = path{ns};
    outdir = fullfile(outputFolder,path_ns(13:end));
    cd(outdir)

    load('ThetaSlowR_5_cell.mat');
    CNEG = [CNEG;cNeg{nt}];
    CZERO = [CZERO;cZero{nt}];
    CPOS = [CPOS;cPos{nt}];
end
ffa = figure('unit','centimeters','position',[20 10 12 10]);
for nl = 1:5%5圈
    
[N,Z,P,na,zb,pc] = scatterphs(CNEG(:,nl),CZERO(:,nl),CPOS(:,nl));
    
X = [ones(length([na,zb,pc]),1) [na,zb,pc]'];
cellphase = [N;Z;P];
y = cellphase;
% linear regress
[b,bint,r,rint,stats] = regress(y,X);

subplot(5,1,nl)
scatter([na,zb,pc],cellphase,50,'LineWidth',1.5,'Marker','o',...
'MarkerEdgeAlpha',0.8)%'MarkerEdgeColor',[0.8,0.8,0.8],
hold on

xFit = -1.5:0.1:1.5;
yFit = b(1) + b(2)*xFit;
plot(xFit,yFit,'k--','LineWidth',1.5)
hold off
xticks([-1,0,1])
yticks([0,360])
ylim([0,360])
title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',10)
ylabel('sg phase')
end
xlabel('sg cycle')

ffb = figure('unit','centimeters','position',[20 10 12 10]);
CNEGx(:,1) = [CNEG(:,1);CNEG(:,2)];
CNEGx(:,2) = [CNEG(:,4);CNEG(:,5)];
CZEROx(:,1) = [CZERO(:,1);CZERO(:,2)];
CZEROx(:,2) = [CZERO(:,4);CZERO(:,5)];
CPOSx(:,1) = [CPOS(:,1);CPOS(:,2)];
CPOSx(:,2) = [CPOS(:,4);CPOS(:,5)];
AAA = [CNEGx(:,1),ones(308,1)*-1;
 CZEROx(:,1),ones(308,1)*0;
 CPOSx(:,1),ones(308,1)*1];

BBB = [CNEGx(:,2),ones(308,1)*-1;
 CZEROx(:,2),ones(308,1)*0;
 CPOSx(:,2),ones(308,1)*1];
for nl = 1:2%前期后期
    
[N,Z,P,na,zb,pc] = scatterphs(CNEGx(:,nl),CZEROx(:,nl),CPOSx(:,nl));
    
X = [ones(length([na,zb,pc]),1) [na,zb,pc]'];
cellphase = [N;Z;P];
y = cellphase;
% linear regress
[b,bint,r,rint,stats] = regress(y,X);

subplot(2,1,nl)
scatter([na,zb,pc],cellphase,50,'LineWidth',1.5,'Marker','o',...
'MarkerEdgeAlpha',0.8)%'MarkerEdgeColor',[0.8,0.8,0.8],
hold on

xFit = -1.5:0.1:1.5;
yFit = b(1) + b(2)*xFit;
plot(xFit,yFit,'k--','LineWidth',1.5)
hold off
xticks([-1,0,1])
yticks([0,360])
ylim([0,360])
title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',10)
ylabel('sg phase')
end
xlabel('sg cycle')
%%

function [neg,zero,pos,Na,Nb,Nc] = scatterphs(cNeg,cZero,cPos)
neg = cNeg(~isnan(cNeg));
zero= cZero(~isnan(cZero));
pos = cPos(~isnan(cPos));

Na = ones(1,length(neg))*-1;
Nb = ones(1,length(zero))*0;
Nc = ones(1,length(pos))*1;

% figure
% scatter(Na,neg)
% hold on
% scatter(Nb,zero)
% scatter(Nc,pos)
% hold off
% xlabel('sg cycle')
% ylabel('sg phase')
% xlim([-1.5,1.5])
end

