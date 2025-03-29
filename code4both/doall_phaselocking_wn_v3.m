
% origin file is doall_phaselocking_wn_v2.m
%
% this file can read those v2 result *.mat files
% concatenate first 2 lap spike arrays and
% calculate PLV on each cell
%
% input CW & CCW
%
% created by WN on 20220725
%
% ==========use timestamp at trial start and end==========
%%
clear
close all
% directories_allData_v2
directories_allData_v0
Directionfolder = {'Tseq\','Tseq-CCW\'};
for D = 1:2
    subfolder = Directionfolder{D};
    for ns = 1:isession
        
        path_ns = path{ns};
        cd(path_ns);
        disp(path_ns)
        
        % step 1 Load spike firing data
        file_input1 = [subfolder 'Cells_singleLap_v2_vel_0.mat'];  % use to get all spikes
        load(file_input1)
        % step 2 load trackdata
        file_input2 = trackdata{ns};
        load(file_input2,'Ts_prerunning','Ts_sample','Ts_test','Ts_posttest')
        % read timestamps
        Ts_start_stop{1,1} = Ts_prerunning{1,1}./ 1000000;
        % Ts_start_stop{1,2} = Ts_sample{1,1}(:,[1,3])./ 1000000;
        % Ts_start_stop{1,3} = Ts_test{1,1}(:,[1,3])./ 1000000;
        % Ts_start_stop{1,4} = Ts_posttest{1,1}(:,[1,3])./ 1000000;
        
        % step 3 load plv result form doall_phaselocking_wn_v2.m
        file_input3 = [subfolder 'data_phaselocking_TSlap_vel0.mat'];
        load(file_input3)
        spikePhase_f2lap = cell(1,3,1);plFeature_f2lap = cell(1,3);
        fprintf(['=====session %.3g start=====' '\n'], ns)
        
        for nc = 1:Ncell
            for nseg = 1%%%only pre-running
                spikePhase_f2lap{nc,1,nseg} = [spikePhase_singlap{nc,1,1},spikePhase_singlap{nc,1,2}];
                spikePhase_f2lap{nc,2,nseg} = [spikePhase_singlap{nc,2,1},spikePhase_singlap{nc,2,2}];
                spikePhase_f2lap{nc,3,nseg} = [spikePhase_singlap{nc,3,1},spikePhase_singlap{nc,3,2}];
                plFeature_f2lap{nseg,1}(nc) = circularStat_WN(spikePhase_f2lap{nc,1,nseg});
                plFeature_f2lap{nseg,2}(nc) = circularStat_WN(spikePhase_f2lap{nc,2,nseg});
                plFeature_f2lap{nseg,3}(nc) = circularStat_WN(spikePhase_f2lap{nc,3,nseg});
            end
        end
        fprintf('\ndone!!\n');
        
        file_output = [subfolder 'data_phaselocking_TSlap_vel0_f2lap.mat'];
        save(file_output,'spikePhase_f2lap','plFeature_f2lap');
    end
end