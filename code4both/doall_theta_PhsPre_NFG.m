% theta 相位进动
% FG-CELL




clear
close all
% directories_allData_v2
directories_allData_v0_allgood
peak0 = 1; % firing rate threshold to remove place cells
nsn = 0;
TTList0 = 'TTList_dCA1_pyr.txt';
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
plot_sign = 0;
start_phase=[];PP_feature=[];

for ns = [1,2,4:6,8,10:13,15,16]%[1:4,6,7]%[1;2;4;5;6;8;10;11;16;20;22;25]'%[1:11,16,20:22,26,27] %1:isession%[1;2;5;6;8;10;11;13;20;22;23;24;27;28]'%
    nsn = nsn+1;
    path_ns = path{ns};
    outFolder = [resuletFolder path_ns(13:end)];
    cd(path_ns);
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};load(trackdata_ns,'Ang_RewardLoc_ontrack')
    fgcell = 0;fgcellinfo = {};
    for D = 1:2
        
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        
        file_input1 = [subfolder1 'Cells_allsegment_v1_vel_0.mat'];  % use to get all spikes
        file_input2 = [subfolder1 'Cells_allsegment_v1_vel_5.mat'];  % use to get ratemap
        file_input3 = [subfolder1 'Cells_singleLap_v2_vel_5.mat'];  % use itself as the decoder
        file_input4 = [subfolder1 'data_phaselocking_TSlap_vel0_new.mat'];  % ** phase locking feature trial TS 用vel0 spike 的结果,不带new也可以
        file_input5 = [subfolder1 'data_phaselocking_TSlap_vel0_f2lap.mat'];  % ** phase locking feature trial TS 用vel0 spike 的前两圈相锁的结果
        % file_input5 = [subfolder1 'data_phaselocking_spkmin.mat'];  % phase locking feature running TS
        fprintf(1,'all spike:\t%s\n%s\ndecoder:\t%s\nphase lock:\t%s\n%s\n',file_input1,file_input2,file_input3,file_input4,file_input5)
        load(file_input4)
        load(file_input5)
        
        if Ncell < 1
            continue
        end
        
        disp(pwd)
        for nt = 1% this file for pre_running only
            switch nt
                case 1
                    laps=5;
                case 2
                    laps=8;
                case 3
                    laps=8;
                case 4
                    laps=6;
            end
            
            % Load Spikes data
            S1 = load(file_input1);  % used to get all spikes
            spikes = S1.spikes;
            S2 = load(file_input2);  % used to get ratemap
            Ratemap_seg = S2.Ratemap_seg{nt};
            mapAxis = S2.mapAxis;
            % Remove the place cells whose peak firing rate<1Hz in all laps
            peak_all = max(Ratemap_seg); % use S2 for detecting
            ind = find(peak_all >= peak0);
            spikes = spikes(ind,:);
            % find cell on track as detecting
            X = cell2struct(S2.fieldProp_seg{nt}(ind),'placefield',1);
            TT0t64 = ReadFileList(S2.TTList0);
            TT0t64 = TT0t64(ind);
            COM = [];
            
            for ncell  = 1:length(ind)
                COM(ncell) = X(ncell).placefield(1).x_COM;
                BIN_st(ncell) = X(ncell).placefield(1).startBin;
                BIN_ed(ncell) = X(ncell).placefield(1).stopBin;
            end
            
            cind_ot = find(COM>Ang_RewardLoc_ontrack(1) & COM<Ang_RewardLoc_ontrack(18));
            cellnumontrack = length(cind_ot);
            BIN_st = BIN_st(cind_ot);
            BIN_ed = BIN_ed(cind_ot);
            fgamma_rayleighP = [plFeature_alllap{1,3}.rayleighP];% all lap 算出来的相锁
            %             fgamma_rayleighP = [plFeature_f2lap{1,3}.rayleighP];% 前两圈 lap 算出来的相锁
            fgamma_rayleighP_ontrack = fgamma_rayleighP(ind);
            fgamma_rayleighP_ontrack = fgamma_rayleighP_ontrack(cind_ot);
            
            fgPhase_singlap_ot = spikePhase_singlap(ind,3,:);
            fgPhase_singlap_ot = fgPhase_singlap_ot(cind_ot,1,:);
            fgPhase_singlap_ot = squeeze(fgPhase_singlap_ot);
            [spkNum_singlap,spkNum_ind] = spknumfinder(fgPhase_singlap_ot,10);
            spkNum_alllap = sum(spkNum_singlap,2);
            cind_ot_fg = find(fgamma_rayleighP_ontrack<0.05 & spkNum_ind'==1)%find(fgamma_rayleighP_ontrack<0.05 & ) % ontrack的ind
            cellnum_fg(nsn,D) = length(cind_ot_fg);
            cind_ot_nfg = setdiff(1:cellnumontrack,cind_ot_fg)
            cellnum_nfg(nsn,D) = length(cind_ot_nfg);
            
            if cellnum_fg(nsn,D)==0
                disp([path_ns '没有fastgamma相锁的神经元'])
                continue
            end
            %% 组合 生成计算theta 相位进动的数据
            icell = cellnum_nfg(nsn,D);
            thPhs_alllap = spikePhase_alllap(ind, 1);
            thPhs_alllap = thPhs_alllap(cind_ot,:);
            Pos_alllap =  spikePos_alllap(ind, 1);
            Pos_alllap = Pos_alllap(cind_ot,:);
            
            for ifgcell = 1:icell
                cellid = cind_ot_nfg(ifgcell);
                x = Pos_alllap{cellid,nt}; %第nc个cell，在第nseg任务中spike的position
                y = thPhs_alllap{cellid,nt}'; %第nc个cell，在第nseg任务中spike的phase
                Xstart = mapAxis(BIN_st(cellid)); % 位置域的边界
                Xend = mapAxis(BIN_ed(cellid));
                
                if Xstart<Xend
                    ind_field = find(x>=Xstart & x<=Xend); %位置域内的点
                    x=x(ind_field);
                    y=y(ind_field);
                    sign_turn = '';
                else % 这是位置域越过2pi的情况
                    cirshift = 2*pi-Xstart;
                    x = mod(x+cirshift,2*pi);
                    Xend = mod(Xend+cirshift,2*pi);
                    ind_field = find(x>=0 & x<=Xend); %位置域内的点
                    x=x(ind_field);
                    y=y(ind_field);
                    sign_turn = '   从2pi位置转到从0开始！！';
                end
                [~,ind]=sort(x);
                temp=y(ind);
                if ~isnan(temp(1))
                    start_phase0 = temp(1) ; % 非nan的起始相位
                else
                    start_phase0 = temp(2) ;% 非nan的起始相位
                end
                
                if length(x)<10
                    continue
                end
                [correlation,pval,slope, phi, RR] = cl_corr(x, y, -2.5,2.5); % slope 单位是deg/m，这里的范围[-2.5,2.5 ]是slope得出的值中间变量的区间，实际范围为[-800 800] for deg/rad, ATTENTION!! 不要随便改！！！
                SLA=[slope,phi,pval];
                phia = phi-slope*median(x); % 进行拟合直线画图的时候，截距是多少
                SP2=[SLA(1),phia];
                
                % ==========================数据存储==============================%
                start_phase = [start_phase;start_phase0];
                PP_feature0 = [ns,... %1-第几天
                    cellid,...            %2-第几个cell
                    nt,...          %3-这一天内的第几个session
                    length(x),...     %4-多少个spike
                    SLA(1),...        %5-拟合的斜率
                    SLA(2),...        %6-截距
                    SLA(3),...        %7-显著性
                    correlation,...   %8-corr-coeff
                    RR];              %9-权重
                PP_feature = [PP_feature;PP_feature0];
                % 画图
                % ==========================画图==================================%
                if plot_sign==1
                    fig = figure;
                    set(fig,'unit','centimeters','position',[25 10 20 12])
                    hold on
                    %plot([x x],[y/180*pi y/180*pi+2*pi],'r.','MarkerSize',9); % 实际位置为横坐标，纵坐标两个周期叠加
                    plot([x x],[y y+360],'r.','MarkerSize',9); % 实际位置为横坐标，纵坐标两个周期叠加
                    plot([0,2*pi],[start_phase0,start_phase0],'r');% start phase
                    axis([0 2*pi,0 4*pi]);
                    tickx=0:pi/2:2*pi;
                    ticky=0:pi:4*pi;
                    axis([0 2*pi,0 720]);
                    tickx=0:pi/2:2*pi;
                    ticky=0:180:720;
                    set(gca, 'XTick',tickx);
                    set(gca, 'XTickLabel',{'0','1/2pi','pi','3/2pi','2pi'});
                    xlabel('Position (rad)')
                    set(gca, 'YTick',ticky);
                    set(gca, 'YTickLabel',{'0','180','360','540','720'});
                    ylabel('Phase (degree)')
                    set(gca,'FontSize',20);
                    
                    h2line = refline(SP2);
                    h2line.Color = 'b';
                    hold off
                    TL1=TT0t64{cellid}; lt=length(TL1); TL1=TL1(1:lt-4);
                    cell_name=strrep(TL1,'_','-');
                    title0=strcat(cell_name,'  spike =',num2str(length(x)),sign_turn);%strcat(cell_name,' lap =',num2str(nl),'  spike =',num2str(length(x)),sign_turn);
                    title1=strcat('slope=',num2str(slope),'   p=',num2str(pval),'   r=',num2str(correlation),'   R=',num2str(RR)); %r是coefficient correlation值，p是显著性，phi是初始相位
                    if SLA(3)<0.05
                        title({title0;title1},'FontSize',15,'color','r');
                    else
                        title({title0;title1},'FontSize',15);
                    end
                    ylim([0 720]);
                    drawnow
                    %     output=strcat(num2str(Groupname),'_Rat',num2str(Rat_id),'_day',num2str(ns),'_seg',num2str(nseg),'_nc',num2str(nc));
                    %     set(gcf,'PaperPositionMode','auto') %把figure按照显示的大小保存
                    %     saveas(gcf,[fig_out_dir,output],'png');
                    %     close all
                end
            end
        end
    end
end



function [N,L] = spknumfinder(spike, threshold)
% 找出来每一圈数量不超过threshold的spike
% 输出每一圈spike数量和每一圈满足要求的cell的逻辑值
[i,j] = size(spike);
for I = 1:i
    for J = 1:j
        N(I,J) = length(spike{I,J});
        
        if N(I,J)>=threshold
            L(I,J) = 1;
        else
            L(I,J) = 0;
        end
    end
end
L = sum(L,2);
L = L == j;
end



function plotPP

end