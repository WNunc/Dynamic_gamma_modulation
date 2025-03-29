%Part2
clear
% close all
clc
directories_allData_v0_allgood
outputFolder = ['H:\neuralynx\phase precession rseult\'];
lapN = [];%['-late'];%['-early'];%['-first']%
for ns = [1,2,4:6,8,10:13,15,16]%16%[3,6,8,9,11,14,16:20]%isession%[1,2,4:6,8,10:13,15,16]%[4,6,8,11,15]%1:isession%[1,5,6,8,11,16]%
    path_ns = path{ns};
    outdir = fullfile(outputFolder,path_ns(13:end));
    disp(outdir)
    cd(outdir)
    file_input=strcat(outdir,'ThetaSlowR_5',lapN,'2.mat'); % choose different files
    load(file_input);
    for nt = 1 % 1pre 2sample 3test
        if ~isempty(Neg)
            NegX{1,ns}=Neg{nt,1}(1,:);
            NegY{1,ns}=Neg{nt,1}(2,:);
        end
        if ~isempty(Zero)
            ZeroX{1,ns}=Zero{nt,1}(1,:);
            ZeroY{1,ns}=Zero{nt,1}(2,:);
        end
        if ~isempty(Pos)
            PosX{1,ns}=Pos{nt,1}(1,:);
            PosY{1,ns}=Pos{nt,1}(2,:);
        end
    end
end
cd(outputFolder)
NegX=cell2mat(NegX);
NegY=cell2mat(NegY);
ZeroX=cell2mat(ZeroX);
ZeroY=cell2mat(ZeroY);
PosX=cell2mat(PosX);
PosY=cell2mat(PosY);
clearvars -except NegX NegY ZeroX ZeroY PosX PosY lapN
% %% do downsample if you need
% ind1=randperm(length(NegX),5442);
% NegX=NegX(ind1);
% NegY=NegY(ind1);
% ind2=randperm(length(ZeroX),5442);
% ZeroX=ZeroX(ind2);
% ZeroY=ZeroY(ind2);
% ind3=randperm(length(PosX),5442);
% PosX=PosX(ind3);
% PosY=PosY(ind3);
% P1=circ_wwtest(NegY,PosY);

%%
x=0:9:360;
y=0:9:360;
X=0:9:720;
Y=0:9:720;
sbin = 9;


p1=hist2d(NegX,NegY,x,y);
P1=p1/sum(sum(p1));
P1P=[P1 P1;P1 P1];
cyc_1 = smooth2a(normz(P1P,0),sbin,sbin);
[r,c] = find(cyc_1 == max(max(cyc_1(:,20:60))));
r = r(c>20&c<60);c = c(c>20&c<60);
% r = mean(r*9);c = mean(c*9);
phase_max(:,1) = r;
figure('Position',[131 537 1637 390])
subplot(1,5,2)
uimagesc(X,Y,[cyc_1]);
hold on
% scatter(X(c),Y(r),'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
scatter(X([34,34]),Y([35,75]),'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
hold off
axis xy;axis equal;
xlabel('Theta phase (degrees)','Fontsize',15)
set(gca,'XTick',0:90:720);
% set(gca,'XTickLabel',{'0','90','180','270','360','450','540','630','720'},'Fontsize',13);
set(gca,'XTickLabel',{'-360','-270','-180','-90','0','90','180','270','360'},'Fontsize',13);
ylabel('Slow gamma phase (degrees)','Fontsize',15)
set(gca,'YTick',0:180:720);
set(gca,'YTickLabel',{'0','180','360','540','720'},'Fontsize',13);
colormap(jet(128));
colorbar
caxis([0.0005 0.0009])
xlim([180 540])
ylim([0,720])
%%
p2=hist2d(ZeroX,ZeroY,x,y);
P2=p2/sum(sum(p2));
P2P=[P2 P2;P2 P2];
cyc_2 = smooth2a(normz(P2P,0),sbin,sbin);
[r,c] = find(cyc_2 == max(max(cyc_2(:,20:60))));
r = r(c>20&c<60);c = c(c>20&c<60);
% r = mean(r*9);c = mean(c*9);
% phase_max(:,2) = r;
% figure
subplot(1,5,3)
uimagesc(X,Y,cyc_2);
hold on
% scatter(X(c),Y(r),'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
scatter(X([40 40]),Y([19,59]),'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])

hold off
axis xy;axis equal;
xlabel('Theta phase (degrees)','Fontsize',15)
set(gca,'XTick',0:90:720);
% set(gca,'XTickLabel',{'0','90','180','270','360','450','540','630','720'},'Fontsize',13);
set(gca,'XTickLabel',{'-360','-270','-180','-90','0','90','180','270','360'},'Fontsize',13);
ylabel('Slow gamma phase (degrees)','Fontsize',15)
set(gca,'YTick',0:180:720);
set(gca,'YTickLabel',{'0','180','360','540','720'},'Fontsize',13);
colormap(jet(128));
colorbar
caxis([0.0005 0.0009])
xlim([180 540])
ylim([0,720])
%%
p3=hist2d(PosX,PosY,x,y);
P3=p3/sum(sum(p3));
P3P=[P3 P3;P3 P3];
% P3P=[P3 P3];
cyc_3 = smooth2a(normz(P3P,0),sbin,sbin);
cyc_31 = smooth2a(normz([P3 P3],0),sbin,sbin);
[r,c] = find(cyc_31 == max(max(cyc_31(:,20:60))));
r = r(c>20&c<60);c = c(c>20&c<60);
% r = mean(r*9);c = mean(c*9);
phase_max(:,3) = r;
% figure
subplot(1,5,4)
uimagesc(X,Y,cyc_3);
hold on
% scatter([351,351],[27,387],'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
% scatter(X(c),Y(r),'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])
scatter(X([48,48]),Y([17,57]),'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1])

hold off
axis xy;axis equal;
xlabel('Theta phase (degrees)','Fontsize',15)
set(gca,'XTick',0:90:720);
% set(gca,'XTickLabel',{'0','90','180','270','360','450','540','630','720'},'Fontsize',13);
set(gca,'XTickLabel',{'-360','-270','-180','-90','0','90','180','270','360'},'Fontsize',13)

ylabel('Slow gamma phase (degrees)','Fontsize',15)
set(gca,'YTick',0:180:720);
set(gca,'YTickLabel',{'0','180','360','540','720'},'Fontsize',13);
colormap(jet(128));
colorbar
caxis([0.0005 0.0009])
xlim([180,540])
ylim([0,720])
% saveas(gcf,['sgcyc ' lapN '2.png'])
% saveas(gcf,['sgcyc ' lapN '2'],'epsc')
