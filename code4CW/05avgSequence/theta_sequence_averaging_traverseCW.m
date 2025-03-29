% crate by WN on 2022/02/18
% find good theta cycle included theta seqence and average them
% 注意解码的参数,此处使用dt40ms_5ms
% use decodingv5 data


clear
close all
% directories_allData_v2
directories_allData_v0
subfolder = 'Tseq\';
case1 = '-ontrack'; %
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
bars = 12;
nbars = 3*bars;
% case3 = sprintf('_%ubars',nbars);
midmod = '_cycmid';
TPsweep = 7;
SPsweep = 10;

peak0 = 1;

for ns = 1:isession
    path_ns = path{ns};
    disp(path_ns)
    cd(path_ns)
    trackdata_ns = trackdata{ns};
    traverseIND = 0;
    phasestep = 10;
    value_stat = [];
    while traverseIND < 36
        cut_phase = traverseIND * phasestep;
        traverseIND = traverseIND+1;
        case3 = num2str(cut_phase);
        mkdir(['CuttingPhaseCW\phase_' case3])
        close all
        theta_info = {};
        nseg = 1;
        m_vel = zeros(1,5);
        sem_vel = zeros(1,5);
        
        ProbDiff = nan(1,5);WeigCorr = nan(1,5);FlipCorr = nan(1,5);
        for nl = 1:5
            file_input1 = strcat(path_ns,subfolder,'scores',num2str(nseg),'-',num2str(nl),case1,'_v5.mat');
            load(file_input1)
            % file_input2 = strcat(path_ns,subfolder,'data_AllLap_thetaphase_',num2str(nbars),'bars-min.mat');
            file_input2 = strcat(path_ns,'CuttingPhaseCW\data_AllLap_thetaphasecut_',case3,'.mat');
            load(file_input2)
            file_input3 = strcat(path_ns,subfolder,'Cells_allsegment_v1_vel_5.mat');
            S = load(file_input3);
            file_input4 = strcat(path_ns,subfolder,'Data_angle_ontrack.mat');
            load(file_input4,'data_angle');
            fprintf(1,'input file:\n%s\n%s\n%s\n%s\n',file_input1,file_input2,file_input3,file_input4)
            lap_ts = scores{nl,6};
            lap_pxn = scores{nl,3};
            lap_bin = scores{nl,4};
            lap_linvel = data_angle{nseg}{nl,1}(:,3);
            lap_angvel = data_angle{nseg}{nl,1}(:,4);
            lap_ts_vel = data_angle{nseg}{nl,1}(:,1);
            lg = movdirection(lap_bin);
            %%
            % find theta cycle on running
            ind_ok = find(thetacycle.cycle_ts{nl}(:,1)>TS_running(nl,1) &...
                thetacycle.cycle_ts{nl}(:,2)<TS_running(nl,2));
            % find speed on running
            ind_tsv = find(lap_ts_vel>=TS_running(nl,1) & ...
                lap_ts_vel<=TS_running(nl,2));
            m_runspeed(nl) = mean(lap_linvel(ind_tsv));
            sem_runspeed(nl) = std(lap_linvel(ind_tsv))/sqrt(length(ind_tsv));
            % 实际用到的运动过程中theta
            theta_info{nl}(:,1) = num2cell(ind_ok);
            % theta cycle持续时间
            theta_info{nl}(:,2) = num2cell(thetacycle.cycle_ts{nl}(ind_ok,2)...
                -thetacycle.cycle_ts{nl}(ind_ok,1));
            % theta cycle 开始和结束对应的时间
            theta_info{nl}(:,3) = num2cell(thetacycle.cycle_ts{nl}(ind_ok,:),2);
            % theta cycle 开始和结束对应的时间的bin
            theta_info{nl}(:,4) = num2cell(thetacycle.cycle_ind{nl}(ind_ok,:),2); % 对不同解码参数使用 ....*10-1,2)
            % 计算每个theta周期的spike个数和cell个数
            spikes = S.spikes;
            spk = [];sp = [];
            ac = [];
            peak_all = max(S.Ratemap_seg{nseg}); % use S2 for detecting
            ind = find(peak_all >= peak0);
            spikes = spikes(ind,:);
            % find cell on track
            spikes = spikes(cind_ot,:);
            cell_ontrack = length(cind_ot)
            for cyclenum = 1:length(ind_ok)
                t_thc = theta_info{nl}{cyclenum,3};
                for n=1:size(spikes,1)
                    ind_sp=find(spikes{n,1}>=t_thc(1) & spikes{n,1}<=t_thc(2));
                    sp(n,cyclenum)=length(ind_sp);
                end
                spk(1,cyclenum)=sum(sp(:,cyclenum));
                ind_ac=find(sp(:,cyclenum)~=0);
                ac(1,cyclenum)=length(ind_ac);
            end
            theta_info{nl}(:,5) = num2cell(ac);% theta cycle的cell数
            theta_info{nl}(:,6) = num2cell(spk);% theta cycle的spike数
            
            % 找到每个theta cycle动物的实际位置bin和每个cycle的pxn，以及sequence重心
            xbins = scores{nl,5};
            diff_comreal = {};
            nonancol = [];
            for ncycle = 1:length(ind_ok)
                ind_thc = theta_info{nl}{ncycle,4};% theta周期的索引
                ts_thc = theta_info{nl}{ncycle,3};% theta周期的开始和结束的时间戳
                [ts_v,cycle_lvel,cycle_avel] = get_thetacycle_vel(lap_ts_vel,lap_linvel,lap_angvel,ts_thc);
                theta_info{nl}{ncycle,11} = cycle_lvel;% 线速度
                theta_info{nl}{ncycle,12} = cycle_avel;% 角速度
                % 每个theta cycle动物的实际位置bin
                theta_info{nl}{ncycle,7} = lap_bin(ind_thc(1):ind_thc(2));
                % 每个cycle的运动方向
                if isempty(find(lg(ind_thc(1):ind_thc(2))==-1))
                    theta_info{nl}{ncycle,10} = 1;
                else
                    theta_info{nl}{ncycle,10} = -1;
                end
                % 每个cycle的pxn
                Pxn_ncycle = lap_pxn(:,ind_thc(1):ind_thc(2));
                theta_info{nl}{ncycle,8} = Pxn_ncycle;
                % 每个cycle的pxn的COM
                bin2use = find(~isnan(sum(Pxn_ncycle)));
                COM_ncycle = nan(1,size(Pxn_ncycle,2));
                for ibins = bin2use
                    COM_ncycle(ibins) = sum(Pxn_ncycle(:,ibins).* xbins,1) ./ sum(Pxn_ncycle(:,ibins),1);
                end
                theta_info{nl}{ncycle,9} = COM_ncycle;
                % 计算每个cycle长度和nan列数的差
                nonancol(ncycle) = length(COM_ncycle)-length(find(isnan(COM_ncycle)));
                % 算com和实际位置的差
                diff_comreal{ncycle,1} = theta_info{nl}{ncycle,9} - xbins(theta_info{nl}{ncycle,7})';
                diff_comreal{ncycle,2} = round(diff_comreal{ncycle,1}./0.0698);% 这个差是几个角度bin
                % 差最小所在的时间bin
                [diff_comreal{ncycle,3}, diff_comreal{ncycle,4}]= min(abs(diff_comreal{ncycle,1}),[],'omitnan');
            end
            
            ac = [theta_info{nl}{:,5}];
            spk =  [theta_info{nl}{:,6}];
            % 找到满足要求的theta cycle
            thetaind = find([theta_info{nl}{:,2}]>0.1&[theta_info{nl}{:,2}]<0.2...时间长度在100ms到200ms之间
                & ac>=ncells&spk>=nspkikes ...3个cells 5个spikes
                & nonancol>=12 ...非nan的值大于等于12个
                & [theta_info{nl}{:,10}]>0 ...正方向的seq
                & ...
                & [diff_comreal{:,3}]<pi);
            thetagood = theta_info{nl}(thetaind,:);
            diffgood = diff_comreal(thetaind,:);
            if isempty(thetaind) || length(thetaind)<3
                disp(['no enough sequence in lap-' num2str(nl)])
                continue
            end
            %%
            pos = [];
            centre_pos = [];
            NPX1 = [];%nan(90,68);% 解码位置-轨迹中点
            NPX2 = [];%nan(90,68);% 解码位置-运动轨迹
            
            mid_cycle = [];mid_cyc = [];% theta cycle 中点
            mid_sequence = [];mid_seq = [];% theta sequence 中点
            ind_pxn = [];
            midsq_pxn = [];%在pxn中的theta cycle 中点
            midcy_pxn = [];%在pxn中的theta sequence 中点
            alignpxn = {};
            for m = 1:length(thetaind)
                pxn = thetagood{m,8};
                realpos = thetagood{m,7};
                realCOM = thetagood{m,9};
                [realPEAK,Pi] = max(pxn);
                T = size(pxn,2);% theta length
                mid_cyc = round(0.5*(T));
                mid_cycle(m) = mid_cyc;
                
                pxn1 = circshift(pxn,-realpos(mid_cyc));
                pxn1 = circshift(pxn1,45);
                pxn2 = [];
                for tn = 1:T
                    pxn2(:,tn) = circshift(pxn(:,tn),-realpos(tn));
                    pxn2(:,tn) = circshift(pxn2(:,tn),45);
                end
                ipart = find(~isnan(realPEAK));
                seq_part = find_continue_part(ipart);
                seq_part = seq_part{1};
                mid_seq = median (seq_part);
                mid_sequence(m) = mid_seq;
                
                %                 bin2use = find(~isnan(sum(pxn2)));
                %                 COM_pxn2 = nan(1,size(pxn2,2));
                %                 xbins_pxn2 = 1:90;
                %                 for ibins = bin2use
                %                     COM_pxn2(ibins) = sum(pxn2(:,ibins).* xbins_pxn2',1) ./ sum(pxn2(:,ibins),1);
                %                 end
                %                 [fitresult, ~] = createseqFit(seq_part, COM_pxn2(seq_part));
                %                 difftemp = round(45-fitresult(mid_seq));
                %                 pxn3 = circshift(pxn2,difftemp);
                
                if m == 1
                    ind_pxn(m,1:2) = [1   size(pxn,2)];
                    midsq_pxn(m) = mid_seq;
                    midcy_pxn(m) = mid_cyc;
                else
                    ind_pxn(m,1:2) = [1   size(pxn,2)]+ ind_pxn(m-1,2);
                    midsq_pxn(m) = mid_seq + ind_pxn(m-1,2);
                    midcy_pxn(m) = mid_cyc + ind_pxn(m-1,2);
                end
                NPX1 = [NPX1 pxn1];
                NPX2 = [NPX2 pxn2];
            end
            
            alignpxn{nl,1} = NPX1;alignpxn{nl,2} = NPX2;
            if size(NPX1,2)<68
                bu = 68 - size(NPX1,2);
                NPX1 = [NPX1 nan(90,bu)];
                NPX2 = [NPX2 nan(90,bu)];
            end
            % 截取theta序列
            TSS = {};%theta sequence structure
            TSS_temp1 = [];
            TSS_temp2 = [];
            iii = 0;
            %         mid_tbins = round(midsq_pxn);
            mid_tbins = round(midcy_pxn);
            for imid = mid_tbins
                iii = iii+1;
                TSS{iii,1}= nan(36,68);
                TSS{iii,2}= nan(36,68);
                
                if imid == mid_tbins(1) || imid-34<1
                    th_a = 1;
                    th_b = imid+33;%34个时间bin，5ms/bin
                    TSS{iii,1}(:,end-th_b+1:end) = NPX1(28:63,th_a:th_b);
                    TSS{iii,2}(:,end-th_b+1:end) = NPX2(28:63,th_a:th_b);
                    
                elseif imid == mid_tbins(end) || imid+34>size(NPX1,2)
                    th_a = imid-34;
                    th_b = size(NPX1,2);
                    TSS{iii,1}(:,1:th_b-th_a+1) = NPX1(28:63,th_a:th_b);
                    TSS{iii,2}(:,1:th_b-th_a+1) = NPX2(28:63,th_a:th_b);
                else
                    th_a = imid-34;
                    th_b = imid+33;
                    TSS{iii,1} = NPX1(28:63,th_a:th_b);
                    TSS{iii,2} = NPX2(28:63,th_a:th_b);
                end
                TSS_temp1(:,:,iii) =  TSS{iii,1};
                TSS_temp2(:,:,iii) =  TSS{iii,2};
            end
            
            TSS_mean1 = mean(TSS_temp1,3,'omitnan');
            TSS_mean2 = mean(TSS_temp2,3,'omitnan');
            %%
            figure(1)
            plot_decodeingandtheta_CW(scores,thetacycle,nl,TS_running,mth)
            saveas(figure(1),[path_ns, 'CuttingPhaseCW\phase_', case3, '\decoding_single_lap',num2str(nl),midmod,'.png'])
            close(figure(1))
            figure(2); %对齐位置
            %         set(gcf, 'position', get(0,'ScreenSize'));
            set(gcf, 'position', [0 666 1920 320]);
            subplot(1,5,nl)
            plot_TSS(TSS_mean1)
            title(['lap-' num2str(nl) ', ThSeq_N_u_m = ' num2str(iii)],'FontSize',14)
            if nl==5
                colorbar('Position',[0.92,0.25,0.01,0.6])
            end
            figure(3);
            set(gcf, 'position', [0 166 1920 320]);
            subplot(1,5,nl)
            plot_TSS(TSS_mean2)
            title(['lap-' num2str(nl) ', ThSeq_N_u_m = ' num2str(iii)],'FontSize',14)
            if nl==5
                colorbar('Position',[0.92,0.25,0.01,0.6])
            end
            
            cycind = [];
            for ii = 1:length(ind_pxn)-1
                cycind(ii) = 0.5*(ind_pxn(ii,2)+ind_pxn(ii+1,1));
            end
            abc = zeros(1,ii)+36;
            %% 速度统计
            m_vel(nl) = mean([thetagood{:,11}]);
            sem_vel(nl) = std([thetagood{:,11}])./sqrt(length(thetaind));
            %%
            save([path_ns,'CuttingPhaseCW\phase_', case3, '\data_theta_seq_info',case1,'_lap',num2str(nl),case2,case3,midmod,'_v5.mat'],...
                'thetaind','ind_ok',...
                'mid_cycle','mid_sequence','midcy_pxn','midsq_pxn',...
                'diff_comreal',...
                'TSS','TSS_mean1','TSS_mean2','iii')
            [ProbDiff(1,nl) ,~] = probdiff(TSS_mean1,TPsweep,SPsweep);
            [WeigCorr(1,nl),~] = weightcorr(TSS_mean1,TPsweep,SPsweep);
            [FlipCorr(1,nl),~] = flipcorr(TSS_mean1,TPsweep,SPsweep);
        end
        if isempty(sem_vel)&&isempty(m_vel)
            disp([path_ns '  ignored'])
            continue
        end
        saveas(figure(2),[path_ns,'CuttingPhaseCW\phase_', case3,'\thetaseq_strcture',case1,'1',case2,case3,midmod,'.png'])
        saveas(figure(3),[path_ns,'CuttingPhaseCW\phase_', case3,'\thetaseq_strcture',case1,'2',case2,case3,midmod,'.png'])
        %     saveas(figure(3),[path_ns,subfolder,'thseq-lap-',subfix,'3-x.png'])
        %     saveas(figure(6),[path_ns,subfolder,'thseq-lap-pxn',subfix,'3-x.png'])
        save([path_ns,'CuttingPhaseCW\phase_', case3,'\data_theta_seq_info','_AllLap',case1,case2,case3,midmod,'_v5.mat'],...
            'theta_info','alignpxn','m_vel','sem_vel')
        m_vel = [];% 平均速度
        sem_vel = [];% 速度标准误
        
        ProbDiff_stat(traverseIND,:) = ProbDiff;
        WeightCorr_stat(traverseIND,:) = WeigCorr;
        FlipCorr_stat(traverseIND,:) = FlipCorr;
    end
    [~,goodphase1]= max(sum(ProbDiff_stat,2,'omitnan'),[],1);
    goodphase1 = (goodphase1-1)*phasestep;
    [~,goodphase2]= max(sum(WeightCorr_stat,2,'omitnan'),[],1);
    goodphase2 = (goodphase2-1)*phasestep;
    [~,goodphase3]= max(sum(FlipCorr_stat,2,'omitnan'),[],1);
    goodphase3 = (goodphase3-1)*phasestep;
    save('CuttingPhaseCW\value_stat.mat','ProbDiff_stat','WeightCorr_stat','FlipCorr_stat',...
        'goodphase1','goodphase2','goodphase3')
end