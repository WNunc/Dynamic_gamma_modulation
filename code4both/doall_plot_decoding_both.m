% 执行画图的操作画出所有解码出来的sequence
% 并显示出weight correlation
clear
close all

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

for ns = isession
    path_ns = path{ns};
    disp(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    cd(outFolder);
    mkdir('eachSequence')
    pause(0.1)
    goodphase_ns = Seq_cutPhase{ns,1};%***CW和CCW分割相位一致，用1和2都一样
    case3 = num2str(goodphase_ns);
    
    % all theta cycle in each lap both Direction
    
    % all theta cycle in each lap exfg both Direction
    
    % all theta cycle in each lap dsfg both Direction
    
    % weight correlation
    WC = load([outFolder 'weight_corr_eachseq.mat']);
    
    for nl = 1:5
        disp(['======== plot lap ' num2str(nl) ' ========'])
        % load good theta cycle in each lap both Direction
        TG = load([outFolder,'data_theta_seq_info',case1{1},'_lap',num2str(nl),case2,case3,midmod,'_v5-both.mat'],'ThetaGood');
        TG = TG.ThetaGood;
        % load good theta cycle in each lap exfg both Direction
        TG_exfg = load([outFolder,'data_theta_seq_info',case1{2},'_lap',num2str(nl),case2,case3,midmod,lockat{2},'_v5-bothv2.mat'],'ThetaGood');
        TG_exfg = TG_exfg.ThetaGood;
        % load good theta cycle in each lap dsfg both Direction
        TG_dsfg = load([outFolder,'data_theta_seq_info',case1{3},'_lap',num2str(nl),case2,case3,midmod,lockat{2},'_v5-bothv3.mat'],'ThetaGood');
        TG_dsfg = TG_dsfg.ThetaGood;
        
        Num_cyc = size(TG,1);
        for t = 1:Num_cyc
            ffname = ['lap ' num2str(nl) ' cyc-' num2str(TG{t,1}) 'slope'];
            % sequence & wc
            all = TG{t,8};wc1 = WC.wc1{nl}(t,1);
            exfg = TG_exfg{t,8};wc2 = WC.wc1{nl}(t,2);
            dsfg = TG_dsfg{t,8};wc3 = WC.wc1{nl}(t,3);
            % real pos, time
            loc = TG{t,7};
            tt = TG{t,3};
            ffa=figure('Position', [2 42 478 952],'Visible','off');
            imageseq(tt,all,loc,wc1)
            title(ffname)
            xlabel('time(s)')
            saveas(ffa,['eachSequence\',ffname,'.png']);
            saveas(ffa,['eachSequence\',ffname],'epsc');
            disp([ffname,'...figure saved'])
            clear ffa
            close all
        end
    end
end






function imageseq(tt,AA,loc,wc)
numbins = 90; % Number of bins
bin_ang = 2*pi/numbins;
mapAxis = bin_ang/2:bin_ang:(2*pi-bin_ang/2);
[~, n] = size(AA);
tt0 = 0:(tt(2)-tt(1))/(n-1):tt(2)-tt(1);
imagesc(tt0,mapAxis,AA)
hold on
plot(tt0,mapAxis(loc),'w--','LineWidth',2)
hold off
[wc(2,1),~] = nanweightcorr(AA,fix((n-1)/2),fix((numbins-1)/2));

[wc(2,1),~] = Cir_reg(AA,fix((n-1)/2),fix((numbins-1)/2));

% 获取图像坐标范围
wc = round(wc,3);
xrange = xlim;
yrange = ylim;
% 计算文本框的位置
xpos = xrange(2) - 0.1 * diff(xrange);
ypos = yrange(1) + 0.05 * diff(yrange);
text(xpos, ypos, num2str(wc), ...
    'Color', 'w',...
    'HorizontalAlignment', 'right',...
    'VerticalAlignment', 'bottom', 'FontSize', 10);

xticks([tt0(1),tt0(n)])
yticks([mapAxis(1),mapAxis(numbins)])
yticklabels({'0','2π'})
axis xy
axis square
axis tight
caxis([0,0.15])
set(gca,'FontSize',12)
colormap(jet(128));
end
