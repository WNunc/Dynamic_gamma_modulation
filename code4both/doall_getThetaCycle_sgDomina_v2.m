% 找到slow gamma能量高于n_std倍标准差的theta cycle
% 用全部的thetacycle----------v2

clear
close all

n_std = 1.5;%0;%1;%2;%

directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
Dx = {'-cw','-ccw'};
case1 = '-ontrack';
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
nsn = 0;
fileoutput = ['sgamma_dominant_thetacyc_IND_std' num2str(n_std) '_v2.mat'];

for ns = 1:isession
    path_ns = path{ns};cd(path_ns)
    disp(path_ns)
    trackdata_ns = trackdata{ns};
    power_fg = cell(1,5);power_sg = cell(1,5);
    power_fgZ = cell(1,5);power_sgZ = cell(1,5);
    ind_th_sl = cell(1,5);ind_th_al = cell(1,5);
    for D = 1:2 % 1 = CW; 2 = CCW;
        goodphase_ns = Seq_cutPhase{ns,D};% CW是第1列
        case3 = num2str(goodphase_ns);
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        x = Dx{D};
        outFolder = [resuletFolder path_ns(13:end)];
        nseg = 1;
        % 输入1：gamma power
        file_input1 = [outFolder,'data_gamma_pow_info_AllLap',case1,case2,case3,midmod,'_v5' ,x,'.mat'];
        load(file_input1,'gpower_seq')
        % % 输入2：theta序列的信息
        % file_input2 = [path_ns,subfolder2,'phase_', case3,'\data_theta_seq_info_AllLap',case1,case2,case3,midmod,'_v5.mat'];
        % load(file_input2)
        % power and its zscore
        for nl = 1:5
            %                 ptemp2_fg = gpower_seq{nl}(thetaind,2);
            %                 ptemp4_fg = gpower_seq{nl}(thetaind,4);
            power_sg{nl} = [power_sg{nl},gpower_seq{nl}(:,1)'];
            power_sgZ{nl} = [power_sgZ{nl},gpower_seq{nl}(:,3)'];
        end
    end
    
    threshold_sl = [];% sgamma-detected theshold single lap
    threshold_al = [];% sgamma-detected theshold all lap
    
    threshold_al = mean([power_sg{:}])+n_std*std([power_sg{:}]);% 计算所有圈阈值
    thr_altemp = threshold_al;
    
    for nl = 1:5
        threshold_sl(nl) = mean([power_sg{nl}])+n_std*std([power_sg{nl}]);% 计算单圈阈值
        thr_sltemp = threshold_sl(nl);
        % 用于Theta_Info变量的ind
        if n_std>0
            ind_th_sl{nl} = find([power_sg{nl}]>thr_sltemp);% sg_domina thetacyc ind single lap
            ind_th_al{nl} = find([power_sg{nl}]>thr_altemp);% sg_domina thetacyc ind all lap
        else
            ind_th_sl{nl} = 1:size(power_sg{nl},2);% sg_domina thetacyc ind single lap
            ind_th_al{nl} = 1:size(power_sg{nl},2);% sg_domina thetacyc ind all lap
        end
    end
    length([ind_th_sl{:}])
    length([ind_th_al{:}])
    nsn = nsn + 1;
    cd(outFolder)
    save([outFolder fileoutput],'ind_th_sl','ind_th_al','threshold_sl','threshold_al',...
        'power_sg','power_sgZ');
end


cd('E:\code\theta_precession_gamma\code4both')