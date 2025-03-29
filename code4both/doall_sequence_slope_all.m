% investigating the relationship between the sequence slope and the phase shift
% step 1 calculate all sequence slope
% load all sequence
% load all phase diff 
% calculate the slope and correlation between slope and shift

% clear
% close all
n_std = 1.5;
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
Dx = {'-cw','-ccw'};
case1 = {'-ontrack','-ontrack_exfg','-ontrack_dsfg'};
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';

fileinput2 = [];% good theta cycle in each lap
fileinput3 = [];% all theta cycle in each lap

lockat = {'firstlap','alllap','f2lap'};
L = 2;
PHSDvsSLOPE = [];
for ns = 16%[1,2,4:6,8,10:13,15,16]%1:isession%
    path_ns = path{ns};
    disp(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    cd(outFolder);
    mkdir('eachSequence')
    pause(0.1)
    goodphase_ns = Seq_cutPhase{ns,1};%***CW和CCW分割相位一致，用1和2都一样
    case3 = num2str(goodphase_ns);
    % all theta cycle in each lap both Direction
    TGa = load([outFolder,'data_theta_seq_info','_AllLap',case1{1},case2,case3,midmod,'_v5-both.mat'],...
        'theta_INFO');
    TGa = TGa.theta_INFO;
    % phase difference
    pd_Folder = ['H:\neuralynx\sgamma result v3 std-' num2str(n_std) '\' path_ns(13:end) '\'];
    pd_file = 'sgphase_difference_v3_3spk.mat';
    load([pd_Folder pd_file]);
    
    slope_good = cell(5,1);
    slope_all = cell(5,1);
    phsvsslope = [];
    sn = 0;
    for nl = 5
        disp(['======== plot lap ' num2str(nl) ' ========'])
%         % load good theta cycle in each lap both Direction
%         TG = load([outFolder,'data_theta_seq_info',case1{1},'_lap',num2str(nl),case2,case3,midmod,'_v5-both.mat'],'ThetaGood');
%         TG = TG.ThetaGood;
%         % good sequence slope
%         Num_cyc = size(TG,1);
%         for t = 1:Num_cyc
%             ffname = ['good lap ' num2str(nl) ' cyc-' num2str(TG{t,1}) '-' num2str(t)];
%             % sequence & wc
%             all = TG{t,8};
%             % real pos, time
%             loc = TG{t,7};
%             tt = TG{t,3};
%             ffa=figure('Position', [2 42 478 952],'Visible','off');
%             [r2,calphase,xaxis,p,slope,xspan,tspan] = imageseq_slope(tt,all,loc);
%             xlabel('time(s)')
%             ylabel('location(rad)')
%             saveas(ffa,['eachSequence\slope-',ffname,'.png']);
%             saveas(ffa,['eachSequence\slope-',ffname],'epsc');
%             disp([ffname,'...figure saved'])
%             clear ffa
%             close all
%             slope_good{nl}{t,1} = r2;
%             slope_good{nl}{t,2} = calphase;
%             slope_good{nl}{t,3} = xaxis;
%             slope_good{nl}{t,4} = p;
%             slope_good{nl}{t,5} = slope;%斜率
%             slope_good{nl}{t,6} = xspan;
%             slope_good{nl}{t,7} = tspan;
%             slope_good{nl}{t,8} = TG{t,1};%theta ID
%         end
        
        % all sequence slope
        Num_cyc = size(TGa{nl},1);
        for t = 1:Num_cyc
            
            ffname = ['all lap ' num2str(nl) ' cyc-' num2str(TGa{nl}{t,1}) '-' num2str(t)];
            % sequence & wc
            all = TGa{nl}{t,8};
            % real pos, time
            loc = TGa{nl}{t,7};
            tt = TGa{nl}{t,3};
            ffa=figure('Position', [2 42 478 952],'Visible','off');
            [r2,calphase,xaxis,p,slope,xspan,tspan] = imageseq_slope(tt,all,loc);
            xlabel('time(s)')
            ylabel('location(rad)')
            saveas(ffa,['eachSequence\',ffname,'.png']);
            saveas(ffa,['eachSequence\',ffname],'epsc');
            disp([ffname,'...figure saved'])
            clear ffa
            close all
            slope_all{nl}{t,1} = r2;
            slope_all{nl}{t,2} = calphase;
            slope_all{nl}{t,3} = xaxis;
            slope_all{nl}{t,4} = p;
            slope_all{nl}{t,5} = slope;
            slope_all{nl}{t,6} = xspan;
            slope_all{nl}{t,7} = tspan;
            slope_all{nl}{t,8} = TGa{nl}{t,1};%theta ID
            inds = find(phasediff1(:,2)==nl & phasediff1(:,3)==t);
            if ~isempty(inds)
                sn = sn+1;
                phsvsslope(sn,1) = phasediff1(inds,1);
                phsvsslope(sn,2) = slope;
                phsvsslope(sn,3) = nl;
                phsvsslope(sn,4) = t;
                phsvsslope(sn,5) = phasediff1(inds,4);%good seq or not
            end
        end
    end
    save([outFolder 'Slope_eachgoodseq.mat'],'slope_good','slope_all','phsvsslope')
    PHSDvsSLOPE = [PHSDvsSLOPE;phsvsslope];
end


X = [ones(length(PHSDvsSLOPE(:,1)),1) PHSDvsSLOPE(:,1)];
slopeall = PHSDvsSLOPE(:,2);
y = slopeall;
% linear regress
[b,bint,r,rint,stats] = regress(y,X);
figure
scatter(PHSDvsSLOPE(:,1),PHSDvsSLOPE(:,2))

xFit = -pi:0.1*pi:pi;
yFit = b(1) + b(2)*xFit;
hold on
plot(xFit,yFit,'k--','LineWidth',1.5)
hold off
% xticks([-2,-1,0,1,2])
% yticks([0,360])
% ylim([0,360])
title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',10)
ylabel('slope')
xlabel('phasediff')

ind_good = PHSDvsSLOPE(:,end)==1;
X = [ones(length(PHSDvsSLOPE(ind_good,1)),1) PHSDvsSLOPE(ind_good,1)];
slopeall = PHSDvsSLOPE(ind_good,2);
y = slopeall;
% linear regress
[b,bint,r,rint,stats] = regress(y,X);
figure
scatter(PHSDvsSLOPE(ind_good,1),PHSDvsSLOPE(ind_good,2))

xFit = -pi:0.1*pi:pi;
yFit = b(1) + b(2)*xFit;
hold on
plot(xFit,yFit,'k--','LineWidth',1.5)
hold off
% xticks([-2,-1,0,1,2])
% yticks([0,360])
% ylim([0,360])
title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',10)
ylabel('slope')
xlabel('phasediff')



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