 %Part2
clear
close all
clc
directories_allData_v0_allgood
outputFolder = ['H:\neuralynx\phase precession rseult\'];
lapN = [];%['early'];%
for ns = [1,2,4:6,8,10:13,15,16]%1:isession
    path_ns = path{ns};
    outdir = fullfile(outputFolder,path_ns(13:end));
    disp(outdir)
    cd(outdir)
    file_input=strcat(outdir,'ThetaFastR_5',lapN,'3.mat'); % choose different files
    load(file_input);
    for nt = 1 % 1pre 2sample 
        if ~isempty(Neg2)%~isempty(Neg2{nt,1}) & ~isempty(Neg2)
            Neg2X{1,ns}=Neg2{nt,1}(1,:);
            Neg2Y{1,ns}=Neg2{nt,1}(2,:);
        end
        if ~isempty(Neg1)
            Neg1X{1,ns}=Neg1{nt,1}(1,:);
            Neg1Y{1,ns}=Neg1{nt,1}(2,:);
        end
        if ~isempty(Zero)
            ZeroX{1,ns}=Zero{nt,1}(1,:);
            ZeroY{1,ns}=Zero{nt,1}(2,:);
        end
        if ~isempty(Pos1)
            Pos1X{1,ns}=Pos1{nt,1}(1,:);
            Pos1Y{1,ns}=Pos1{nt,1}(2,:);
        end
        if ~isempty(Pos2)
            Pos2X{1,ns}=Pos2{nt,1}(1,:);
            Pos2Y{1,ns}=Pos2{nt,1}(2,:);
        end
    end
end
cd(outputFolder)
Neg2X=cell2mat(Neg2X);
Neg2Y=cell2mat(Neg2Y);
Neg1X=cell2mat(Neg1X);
Neg1Y=cell2mat(Neg1Y);
ZeroX=cell2mat(ZeroX);
ZeroY=cell2mat(ZeroY);
Pos1X=cell2mat(Pos1X);
Pos1Y=cell2mat(Pos1Y);
Pos2X=cell2mat(Pos2X);
Pos2Y=cell2mat(Pos2Y);
clearvars -except Neg1X Neg1Y Neg2X Neg2Y ZeroX ZeroY Pos1X Pos1Y Pos2X Pos2Y lapN

%%
x=0:9:360;
y=0:9:360;
X=0:9:720;
sbin = 10;
p1=hist2d(Neg2X,Neg2Y,x,y);
P1=p1/sum(sum(p1));
P1P=[P1 P1];
cyc_1 = smooth2a(normz(P1P,0),sbin,sbin);
[r,c] = find(cyc_1 == max(max(cyc_1(:,20:60))));
r = r(c>20&c<60);c = c(c>20&c<60);
r = mean(r*9);c = mean(c*9);
phase_max(1) = r;
figure('Position', [2 454 1916 450])
subplot(1,5,1)
uimagesc(X,y,cyc_1);
hold on
plot(c,r,'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
hold off
axis xy;axis square;
xlabel('Theta phase (degrees)','Fontsize',15)
set(gca,'XTick',0:90:720);
set(gca,'XTickLabel',{'0','90','180','270','360','450','540','630','720'},'Fontsize',13);
ylabel('Fast gamma phase (degrees)','Fontsize',15)
set(gca,'YTick',0:90:360);
set(gca,'YTickLabel',{'0','90','180','270','360'},'Fontsize',13);
colormap(jet(128));
colorbar
caxis([0.0005 0.0009])
xlim([180 540])

%%
p2=hist2d(Neg1X,Neg1Y,x,y);
P2=p2/sum(sum(p2));
P2P=[P2 P2];
cyc_2 = smooth2a(normz(P2P,0),sbin,sbin);
[r,c] = find(cyc_2 == max(max(cyc_2(:,20:60))));
r = r(c>20&c<60);c = c(c>20&c<60);
r = mean(r*9);c = mean(c*9);
phase_max(2) = r;
subplot(1,5,2)
uimagesc(X,y,smooth2a(normz(P2P,0),sbin,sbin));
hold on
plot(c,r,'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
hold off
axis xy;axis square;
xlabel('Theta phase (degrees)','Fontsize',15)
set(gca,'XTick',0:90:720);
set(gca,'XTickLabel',{'0','90','180','270','360','450','540','630','720'},'Fontsize',13);
ylabel('Fast gamma phase (degrees)','Fontsize',15)
set(gca,'YTick',0:90:360);
set(gca,'YTickLabel',{'0','90','180','270','360'},'Fontsize',13);
colormap(jet(128));
colorbar
caxis([0.0005 0.0009])
xlim([180 540])

%%
p3=hist2d(ZeroX,ZeroY,x,y);
P3=p3/sum(sum(p3));
P3P=[P3 P3];
cyc_3 = smooth2a(normz(P3P,0),sbin,sbin);
[r,c] = find(cyc_3 == max(max(cyc_3(:,20:60))));
r = r(c>20&c<60);c = c(c>20&c<60);
r = mean(r*9);c = mean(c*9);
phase_max(3) = r;
subplot(1,5,3)
uimagesc(X,y,cyc_3);
hold on
plot(c,r,'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
hold off
axis xy;axis square;
xlabel('Theta phase (degrees)','Fontsize',15)
set(gca,'XTick',0:90:720);
set(gca,'XTickLabel',{'0','90','180','270','360','450','540','630','720'},'Fontsize',13);
ylabel('Fast gamma phase (degrees)','Fontsize',15)
set(gca,'YTick',0:90:360);
set(gca,'YTickLabel',{'0','90','180','270','360'},'Fontsize',13);
colormap(jet(128));
colorbar
caxis([0.0005 0.0009])
xlim([180 540])

%%
p4=hist2d(Pos1X,Pos1Y,x,y);
P4=p4/sum(sum(p4));
P4P=[P4 P4];
cyc_4 = smooth2a(normz(P4P,0),sbin,sbin);
[r,c] = find(cyc_4 == max(max(cyc_4(:,20:60))));
r = r(c>20&c<60);c = c(c>20&c<60);
r = mean(r*9);c = mean(c*9);
phase_max(4) = r;
subplot(1,5,4)
uimagesc(X,y,cyc_4);
hold on
plot(c,r,'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
hold off
axis xy;axis square;
xlabel('Theta phase (degrees)','Fontsize',15)
set(gca,'XTick',0:90:720);
set(gca,'XTickLabel',{'0','90','180','270','360','450','540','630','720'},'Fontsize',13);
ylabel('Fast gamma phase (degrees)','Fontsize',15)
set(gca,'YTick',0:90:360);
set(gca,'YTickLabel',{'0','90','180','270','360'},'Fontsize',13);
colormap(jet(128));
colorbar
caxis([0.0005 0.0009])
xlim([180 540])

%%
p5=hist2d(Pos2X,Pos2Y,x,y);
P5=p5/sum(sum(p5));
P5P=[P5 P5];
cyc_5 = smooth2a(normz(P5P,0),sbin,sbin);
[r,c] = find(cyc_5 == max(max(cyc_5(:,20:60))));
r = r(c>20&c<60);c = c(c>20&c<60);
r = mean(r*9);c = mean(c*9);
phase_max(5) = r;
subplot(1,5,5)
uimagesc(X,y,cyc_5);
hold on
plot(c,r,'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
hold off
axis xy;axis square;
xlabel('Theta phase (degrees)','Fontsize',15)
set(gca,'XTick',0:90:720);
set(gca,'XTickLabel',{'0','90','180','270','360','450','540','630','720'},'Fontsize',13);
ylabel('Fast gamma phase (degrees)','Fontsize',15)
set(gca,'YTick',0:90:360);
set(gca,'YTickLabel',{'0','90','180','270','360'},'Fontsize',13);
colormap(jet(128));
colorbar
caxis([0.0005 0.0009])
xlim([180 540])
% saveas(gcf,['fgcyc ' lapN '2.png'])
% saveas(gcf,['fgcyc ' lapN '2'],'epsc')