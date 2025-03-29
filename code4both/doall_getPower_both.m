clear
close all
% directories_allData_v2
directories_allData_v0
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
Dx = {'-cw','-ccw'};
case1 = '-ontrack'; %
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';

for ns = 28:30 %1:isession
    for D = 1:2
        path_ns = path{ns};
        cd(path_ns)
        disp(path_ns)
        trackdata_ns = trackdata{ns};
        close all
        outFolder = [resuletFolder path_ns(13:end)];
        inFolder = outFolder;
        nseg = 1;
        goodphase_ns = Seq_cutPhase{ns,2};%CW是第2列
        case3 = num2str(goodphase_ns);
        file_input1 = [path_ns,Phasecutfolder{D},'phase_', case3,'\data_theta_seq_info_AllLap',case1,case2,case3,midmod,'_v5.mat'];
        load(file_input1)
        file_output = [outFolder,'data_gamma_pow_info_AllLap',case1,case2,case3,midmod,'_v5' ,Dx{D},'.mat'];
        fprintf(1,'output file:\n%s\n',file_output)
        gpower_seq = cell(1,5);
        % power and its zscore
        for nl = 1:5
            file_input2 = strcat(path_ns,Directionfolder{D},'scores',num2str(nseg),'-',num2str(nl),case1,'_v5.mat');
            load(file_input2)
            fprintf(1,'input file:\n%s\n%s\n',file_input1,file_input2)
            
            lap_EEG_ts = scores{nl,10};
            lap_EEG_sample = scores{nl,11};
            % find theta cycle ts on running
            t_thc = cell2mat(theta_info{nl}(:,3));
            ind_tc = [];
            for m = 1:size(t_thc,1)
                for n = 1:2
                    dis_t=abs(lap_EEG_ts-t_thc(m,n));
                    [~,ind_tc(m,n)]=min(dis_t);
                end
            end
            
            TFRsg=[];
            TFRfg=[];
            for ntt = 1:size(lap_EEG_sample,1)
                TFRsg(ntt,:) = TFR_frequency_band(lap_EEG_sample(ntt,:)',2000,5,25,45); % slow gamma
                TFRfg(ntt,:) = TFR_frequency_band(lap_EEG_sample(ntt,:)',2000,5,65,100); % fast gamma
            end
            Zsg=zscore(TFRsg,0,2); % zscore between tetrodes
            Zfg=zscore(TFRfg,0,2);
            PowSG = []; PowFG = []; Mzsg = []; Mzfg = [];
            Zsg_s=[]; Zsg_m=[]; Zfg_s=[]; Zfg_m=[];
            for m = 1:size(t_thc,1)
                Zsg_s=Zsg(:,ind_tc(m,1):ind_tc(m,2));
                Zfg_s=Zfg(:,ind_tc(m,1):ind_tc(m,2));
                Psg=TFRsg(:,ind_tc(m,1):ind_tc(m,2));
                Pfg=TFRfg(:,ind_tc(m,1):ind_tc(m,2));
                PowSG(m,1)=mean(Psg(:));%**
                PowFG(m,1)=mean(Pfg(:));%**
                Zsg_m=mean(Zsg_s(:));
                Zfg_m=mean(Zfg_s(:));
                Mzsg(m,1)=Zsg_m;%**
                Mzfg(m,1)=Zfg_m;%**
            end
            
            gpower_seq{nl} = [PowSG,PowFG,Mzsg,Mzfg];
        end
        save(file_output,'gpower_seq','TFRsg','TFRfg');
    end
end
cd E:\code\theta_precession_gamma\