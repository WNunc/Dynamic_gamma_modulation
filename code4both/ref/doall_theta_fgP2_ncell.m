%Part2
% 按cell
clear
% close all
clc
directories_allData_v0_allgood
outputFolder = ['H:\neuralynx\phase precession rseult\'];
%%
CNEG2 = [];CNEG1 = []; CZERO = []; CPOS1 = [];CPOS2 = [];
nt = 1;
for ns = [1,2,4:6,8,10:13,15,16]%1:isession
    path_ns = path{ns};
    outdir = fullfile(outputFolder,path_ns(13:end));
    cd(outdir)

    load('ThetaFastR_5_cell.mat');
    CNEG2 = [CNEG2;cNeg2{nt}];
    CNEG1 = [CNEG1;cNeg1{nt}];
    CZERO = [CZERO;cZero{nt}];
    CPOS1 = [CPOS1;cPos1{nt}];
    CPOS2 = [CPOS2;cPos2{nt}];
end
ffa = figure('unit','centimeters','position',[20 10 12 10]);
for nl = 1:5%5圈
    
[N2,N1,Z,P1,P2,n2a,n1b,zc,p1d,p2e] = scatterphs(CNEG2(:,nl),CNEG1(:,nl),CZERO(:,nl),CPOS1(:,nl),CPOS2(:,nl));
    
X = [ones(length([n2a,n1b,zc,p1d,p2e]),1) [n2a,n1b,zc,p1d,p2e]'];
cellphase = [N2;N1;Z;P1;P2];
y = cellphase;
% linear regress
[b,bint,r,rint,stats] = regress(y,X);

subplot(5,1,nl)
scatter([n2a,n1b,zc,p1d,p2e],cellphase,50,'LineWidth',1.5,'Marker','o',...
'MarkerEdgeAlpha',0.8)%'MarkerEdgeColor',[0.8,0.8,0.8],
hold on

xFit = -2.5:0.1:2.5;
yFit = b(1) + b(2)*xFit;
plot(xFit,yFit,'k--','LineWidth',1.5)
hold off
xticks([-1,0,1])
yticks([0,360])
ylim([0,360])
title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',10)
ylabel('fg phase')
end
xlabel('fg cycle')

ffb = figure('unit','centimeters','position',[20 10 12 10]);
CNEG2x(:,1) = [CNEG2(:,1);CNEG2(:,2)];
CNEG2x(:,2) = [CNEG2(:,4);CNEG2(:,5)];
CNEG1x(:,1) = [CNEG1(:,1);CNEG1(:,2)];
CNEG1x(:,2) = [CNEG1(:,4);CNEG1(:,5)];
CZEROx(:,1) = [CZERO(:,1);CZERO(:,2)];
CZEROx(:,2) = [CZERO(:,4);CZERO(:,5)];
CPOS1x(:,1) = [CPOS1(:,1);CPOS1(:,2)];
CPOS1x(:,2) = [CPOS1(:,4);CPOS1(:,5)];
CPOS2x(:,1) = [CPOS2(:,1);CPOS2(:,2)];
CPOS2x(:,2) = [CPOS2(:,4);CPOS2(:,5)];
for nl = 1:2%前期后期
    
[N2,N1,Z,P1,P2,n2a,n1b,zc,p1d,p2e] = scatterphs(CNEG2x(:,nl),CNEG1x(:,nl),CZEROx(:,nl),CPOS1x(:,nl),CPOS2x(:,nl));
    
X = [ones(length([n2a,n1b,zc,p1d,p2e]),1) [n2a,n1b,zc,p1d,p2e]'];
cellphase = [N2;N1;Z;P1;P2];
y = cellphase;
% linear regress
[b,bint,r,rint,stats] = regress(y,X);

subplot(2,1,nl)
scatter([n2a,n1b,zc,p1d,p2e],cellphase,50,'LineWidth',1.5,'Marker','o',...
'MarkerEdgeAlpha',0.8)%'MarkerEdgeColor',[0.8,0.8,0.8],
hold on

xFit = -2.5:0.1:2.5;
yFit = b(1) + b(2)*xFit;
plot(xFit,yFit,'k--','LineWidth',1.5)
hold off
xticks([-2,-1,0,1,2])
yticks([0,360])
ylim([0,360])
title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',10)
ylabel('fg phase')
end
xlabel('fg cycle')
%%

function [neg2,neg1,zero,pos1,pos2,Na,Nb,Nc,Nd,Ne] = scatterphs(cNeg2,cNeg1,cZero,cPos1,cPos2)
neg2 = cNeg2(~isnan(cNeg2));
neg1 = cNeg1(~isnan(cNeg1));
zero= cZero(~isnan(cZero));
pos1 = cPos1(~isnan(cPos1));
pos2 = cPos2(~isnan(cPos2));

Na = ones(1,length(neg2))*-2;
Nb = ones(1,length(neg1))*-1;
Nc = ones(1,length(zero))*0;
Nd = ones(1,length(pos1))*1;
Ne = ones(1,length(pos2))*2;

% figure
% scatter(Na,neg2)
% hold on
% scatter(Nb,neg1)
% scatter(Nc,zero)
% scatter(Nd,pos1)
% scatter(Ne,pos2)
% hold off
% xlabel('sg cycle')
% ylabel('sg phase')
% xlim([-2.5,2.5])
end

