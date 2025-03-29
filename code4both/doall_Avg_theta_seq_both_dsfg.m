% crate by WN on 2022/06/28
% 找到好的相位后使用
% 用原始条件的窗口截取排除相同数量的fastgamma相锁神经元解码出来的TTS
% 注意解码的参数,此处使用dt40ms_5ms
% use decodingv5 data
%%
clear
close all
% directories_allData_v0
directories_allData_v0_allgood

resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
case1 = '-ontrack_dsfg'; %exclude non-fgamma
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';

peak0 = 1;
TSSall = cell(5,1);
lockat = {'firstlap','alllap','f2lap'};

for L = 2%1:3
    for ns = [1,2,4:6,8,10:13,15,16]%1:isession
        path_ns = path{ns};
        disp(path_ns);cd(path_ns);
        trackdata_ns = trackdata{ns};
        nseg = 1;
        outFolder = [resuletFolder path_ns(13:end)];
        alignpxn = {}; % 对齐后的pxn
        theta_INFO  =cell(1,5);% CW 和 CCW 总的theta infomation
        
        
        for nl = 1:5
            ThetaGood = {};
            for D = 1:2
                theta_info = {};
                
                goodphase_ns = Seq_cutPhase{ns,D};%***
                subfolder1 = Directionfolder{D};
                subfolder2 = Phasecutfolder{D};
                case3 = num2str(goodphase_ns);
                
                file_input1 = strcat(path_ns,subfolder1,...
                    'scores',num2str(nseg),'-',num2str(nl),case1,'_TSlap_vel0',lockat{L},'v3.mat');
                if ~exist(file_input1,'file')
                    disp([ lockat{L},': lap' num2str(nl) 'have no phase lock cell'])
                    continue
                end
                load(file_input1)
                
                file_input2 = strcat(path_ns,subfolder2,'data_AllLap_thetaphasecut_',case3,'.mat');
                load(file_input2)
                
                file_input3 = strcat(path_ns,subfolder1,'Cells_allsegment_v1_vel_0.mat');% load spike
                S1 = load(file_input3);
                
                file_input4 = strcat(path_ns,subfolder1,'Cells_allsegment_v1_vel_5.mat');% load rate map
                S2 = load(file_input4);
                
                file_input5 = strcat(path_ns,subfolder1,'Data_angle_ontrack.mat');
                load(file_input5,'data_angle');
                
                file_input6 = strcat(path_ns,subfolder2,...
                    'phase_',case3,'\data_theta_seq_info-ontrack_lap',num2str(nl),case2,case3,midmod,'_v5.mat');
                if ~exist(file_input6,'file')
                    disp([ Directionfolder{D},lockat{L},': lap' num2str(nl) 'have no enough sequence'])
                    continue
                end
                load(file_input6)
                fprintf(1,'input file:\n%s\n%s\n%s\n%s\n',file_input1,file_input2,file_input3,file_input5)
                fprintf(1,'theta windows file:\n%s\n',file_input6)
                
                lap_ts = scores{nl,6};
                if D == 1
                    lap_pxn = scores{nl,3};
                    lap_bin = scores{nl,4};
                else
                    lap_pxn = flipud(scores{nl,3});
                    lap_bin = 90-scores{nl,4};
                end
                lap_linvel = data_angle{nseg}{nl,1}(:,3);
                lap_angvel = data_angle{nseg}{nl,1}(:,4);
                lap_ts_vel = data_angle{nseg}{nl,1}(:,1);
                lg = movdirection(lap_bin);
                %% 画原始解码图
                figure(1)
                if D == 1
%                     plot_decodeingandtheta_CW(scores,thetacycle,nl,TS_running,mth)
%                     saveas(figure(1),[outFolder, '1-CWphase_', case3, '_decoding_single_lap',num2str(nl),case1,lockat{L},'v3.png'])
                else
%                     plot_decodeingandtheta_CCW(scores,thetacycle,nl,TS_running,mth)
%                     saveas(figure(1),[outFolder, '2-CCWphase_', case3, '_decoding_single_lap',num2str(nl),case1,lockat{L},'v3.png'])
                end
                clf
                %%
                % find theta cycle on running
                %         ind_ok = find(thetacycle.cycle_ts{nl}(:,1)>TS_running(nl,1) &...
                %             thetacycle.cycle_ts{nl}(:,2)<TS_running(nl,2));
                
                % find speed on running
                ind_tsv = find(lap_ts_vel>=TS_running(nl,1) & ...
                    lap_ts_vel<=TS_running(nl,2));
                
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
                spikes = S1.spikes; % use S1 spike
                spk = [];sp = [];
                ac = [];
                peak_all = max(S2.Ratemap_seg{nseg}); % use S2 ratemap
                ind = find(peak_all >= peak0);
                spikes = spikes(ind,:);
                % find cell on track and exclude fgamma phase locking cell
                spikes = spikes(cind_ot_exnonfg ,:);
                cell_ontrack = length(cind_ot_exnonfg)
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
                    ind_thc = theta_info{nl}{ncycle,4};
                    ts_thc = theta_info{nl}{ncycle,3};
                    [ts_v,cycle_lvel,cycle_avel] = get_thetacycle_vel(lap_ts_vel,lap_linvel,lap_angvel,ts_thc);
                    theta_info{nl}{ncycle,11} = cycle_lvel;
                    theta_info{nl}{ncycle,12} = cycle_avel;
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
                    [diff_comreal{ncycle,3}, diff_comreal{ncycle,4}]= min(abs(diff_comreal{ncycle,1}./0.0698),[],'omitnan');
                end
                
                ac = [theta_info{nl}{:,5}];
                spk =  [theta_info{nl}{:,6}];
                % 用与全部神经元结果相同的theta cycle
                
                %thetaind = find([theta_info{nl}{:,2}]>0.1&[theta_info{nl}{:,2}]<0.2...时间长度在100ms到200ms之间
                %    & ac>=ncells&spk>=nspkikes ...3个cells 6个spikes
                %    & nonancol>=16 ...非nan的值大于等于16个
                %    & [theta_info{nl}{:,10}]>0 ...正方向的seq
                %    & [diff_comreal{:,3}]<pi);
                
                thetagood = theta_info{nl}(thetaind,:);
                diffgood = diff_comreal(thetaind,:);
                theta_INFO{nl} = [theta_INFO{nl};theta_info{nl}];
                ThetaGood = [ThetaGood; thetagood];
            end
            %%
            NPX1 = [];% 解码位置-轨迹中点
            NPX2 = [];% 解码位置-运动轨迹
            mid_cycle = [];mid_cyc = [];% theta cycle 中点
            %         mid_sequence = [];mid_seq = [];% theta sequence 中点
            ind_pxn = [];
            midsq_pxn = [];%在pxn中的theta cycle 中点
            midcy_pxn = [];%在pxn中的theta sequence 中点
            for m = 1:size(ThetaGood,1)
                pxn = ThetaGood{m,8};
                realpos = ThetaGood{m,7};
                realCOM = ThetaGood{m,9};
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
                %             ipart = find(~isnan(realPEAK));
                %             seq_part = find_continue_part(ipart);
                %             seq_part = seq_part{1};
                %             mid_seq = median (seq_part);
                %             mid_sequence(m) = mid_seq;
                if m == 1
                    ind_pxn(m,1:2) = [1   size(pxn,2)];
                    %                 midsq_pxn(m) = mid_seq;
                    midcy_pxn(m) = mid_cyc;
                else
                    ind_pxn(m,1:2) = [1   size(pxn,2)]+ ind_pxn(m-1,2);
                    %                 midsq_pxn(m) = mid_seq + ind_pxn(m-1,2);
                    midcy_pxn(m) = mid_cyc + ind_pxn(m-1,2);
                end
                NPX1 = [NPX1 pxn1];
                NPX2 = [NPX2 pxn2];
            end
            alignpxn{nl,1} = NPX1;alignpxn{nl,2} = NPX2;
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
            end
            if isempty(TSS)
                continue
            end
            
            A = size(TSS{1,1});
            B = size(TSS,1);
            TSS_temp1 = reshape(cell2mat(TSS(:,1)),[A(1) B A(2)]);
            TSS_temp2 = reshape(cell2mat(TSS(:,2)),[A(1) B A(2)]);
            TSS_mean1 = squeeze( mean(TSS_temp1,2,'omitnan'));
            TSS_mean2 = squeeze( mean(TSS_temp2,2,'omitnan'));
            %% 画出 平均后的前扫结构
%             figure(2); %对齐位置
%             set(gcf, 'position', [0 666 1920 320]);
%             subplot(1,5,nl)
%             plot_TSS(TSS_mean1)
%             title(['lap-' num2str(nl) ', ThSeq_N_u_m = ' num2str(iii)],'FontSize',14)
%             colorbar('Position',[0.92,0.586,0.01,0.33])
%             figure(3);
%             set(gcf, 'position', [0 166 1920 320]);
%             subplot(1,5,nl)
%             plot_TSS(TSS_mean2)
%             title(['lap-' num2str(nl) ', ThSeq_N_u_m = ' num2str(iii)],'FontSize',14)
%             colorbar('Position',[0.92,0.586,0.01,0.33])
            
            cycind = [];
            for ii = 1:length(ind_pxn)-1
                cycind(ii) = 0.5*(ind_pxn(ii,2)+ind_pxn(ii+1,1));
            end
            abc = zeros(1,ii)+36;
            
%             figure(4)
%             set(gcf, 'position', get(0,'ScreenSize'));
%             subplot(5,1,nl)
%             plot_NPX(NPX1(28:63,:),mid_tbins,cycind,abc)
%             figure(5)
%             set(gcf, 'position', get(0,'ScreenSize'));
%             subplot(5,1,nl)
%             plot_NPX(NPX2(28:63,:),mid_tbins,cycind,abc)
            
            %%
%             save([outFolder,'data_theta_seq_info',case1,'_lap',num2str(nl),case2,case3,midmod,lockat{L},'_v5-bothv3.mat'],...
%                 'ThetaGood','mid_cycle','midcy_pxn',...
%                 'mid_tbins','cycind','abc',...
%                 'TSS','TSS_mean1','TSS_mean2','iii')
            
            TSSall{nl} = [TSSall{nl};TSS]; 
        end
        
%         saveas(figure(2),[outFolder,'thetaseq_strcture',case1,'1',case2,case3,midmod,lockat{L},'-bothv3.png']);clf
%         saveas(figure(3),[outFolder,'thetaseq_strcture',case1,'2',case2,case3,midmod,lockat{L},'-bothv3.png']);clf
%         
%         saveas(figure(4),[outFolder,'thetaseq_pxn',case1,'1',case2,case3,midmod,lockat{L},'-bothv3.png']);clf
%         saveas(figure(5),[outFolder,'thetaseq_pxn',case1,'2',case2,case3,midmod,lockat{L},'-bothv3.png']);clf
%         
%         save([outFolder,'data_theta_seq_info','_AllLap',case1,case2,case3,midmod,lockat{L},'_v5-bothv3.mat'],...
%             'theta_INFO','alignpxn')
    end
end
