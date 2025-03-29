%%
clear
close all
% directories_allData_v0
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
Dx = {'-cw','-ccw'};
case1 = '-ontrack'; %
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
nsn = 0;
% 
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%[2;4;5;6;8;10;11;16;20;22;25]'%[1,2,4,10,21,24]%[1:11,16,20:22,25,26,28,29,30]%1:isession
    path_ns = path{ns};
    disp(path_ns)
    trackdata_ns = trackdata{ns};
    power_fg = cell(1,5);power_sg = cell(1,5);
    power_fgZ = cell(1,5);power_sgZ = cell(1,5);
    for D = 1:2 % 1 = CW; 2 = CCW;
        goodphase_ns = Seq_cutPhase{ns,D};%CW是第1列
        case3 = num2str(goodphase_ns);
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        x = Dx{D};
        outFolder = [resuletFolder path_ns(13:end)];
        nseg = 1;
        file_input1 = [outFolder,'data_gamma_pow_info_AllLap',case1,case2,case3,midmod,'_v5' ,x,'.mat'];
        load(file_input1,'gpower_seq')
        
        % power and its zscore
        for nl = 1:5
            file_input2 = [path_ns,subfolder2,'phase_', case3,'\data_theta_seq_info',case1,'_lap',num2str(nl),case2,case3,midmod,'_v5.mat'];

            if exist(file_input2,'file')
            load(file_input2,'thetaind')
            fprintf(1,'input file:\n%s\n',file_input2)
            %         sgpower_lap(ns,nl) = mean(gpower_seq(thetaind,1));
            %         fgpower_lap(ns,nl) = mean(gpower_seq{nl}(thetaind,2));
            %         sgpowerZ_lap(ns,nl) = mean(gpower_seq(thetaind,3));
            %         fgpowerZ_lap(ns,nl) = mean(gpower_seq{nl}(thetaind,4));
            ptemp2_fg = gpower_seq{nl}(thetaind,2);
            ptemp4_fg = gpower_seq{nl}(thetaind,4);
            ptemp1_sg = gpower_seq{nl}(thetaind,1);
            ptemp3_sg = gpower_seq{nl}(thetaind,3);
            else
                ptemp2_fg = [];
                ptemp4_fg = [];
            end
            power_fg{nl} = [power_fg{nl};ptemp2_fg];
            power_fgZ{nl} = [power_fgZ{nl};ptemp4_fg];
            
            power_sg{nl} = [power_sg{nl};ptemp1_sg];
            power_sgZ{nl} = [power_sgZ{nl};ptemp3_sg];
        end
    end
    
    nsn = nsn + 1;
    
    for nl = 1:5
            fgpower_lap(nsn,nl) = mean(power_fg{nl});
            fgpowerZ_lap(nsn,nl) = mean(power_fgZ{nl});
            sgpower_lap(nsn,nl) = mean(power_sg{nl});
            sgpowerZ_lap(nsn,nl) = mean(power_sgZ{nl});
    end
end


plotg(fgpower_lap)
ylabel('Fast Gamma Power(μV^2)')
axis square
plotg(sgpower_lap)
ylabel('Slow Gamma Power(μV^2)')
axis square
function plotg(power)
n = size(power,1);
mean_power = mean(power,1);
sem_power = std(power,[],1)/sqrt(n);
lap = size(power,2);
figure('Color',[1 1 1])
errorbar(1:lap,mean_power,sem_power,'LineWidth',1.5)
xlim([0.5,5.5])
ylim([min(mean_power)*0.8,max(mean_power)*1.2])
xlabel('Lap Num')
set(gca,'FontSize',15,'Color','none','LineWidth',1)
box off

end



