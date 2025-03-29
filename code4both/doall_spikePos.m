% created by WN on 20230117
% original code: doall_phaselocking_wn_v2.m
%%
clear
close all
% directories_allData_v2
directories_allData_v0
subfolder = {'Tseq\','Tseq-CCW\'};

for ns = 1:isession
    
    path_ns = path{ns};
    cd(path_ns);
    disp(path_ns)
    for D = 1:2
        % step 1 Load spike firing data
        file_input1 = [subfolder{D} 'Cells_singleLap_v2_vel_0.mat'];  % use to get all spikes
        load(file_input1)
        
        file_input2 = trackdata{ns};
        load(file_input2,'Ts_prerunning','Ts_sample','Ts_test','Ts_posttest')
        % read timestamps
        Ts_start_stop{1,1} = Ts_prerunning{1,D}./ 1000000;%1=CW；2=CCW
        % Ts_start_stop{1,1} = Ts_prerunning{1,1}./ 1000000;%CW
        % Ts_start_stop{1,2} = Ts_sample{1,1}(:,[1,3])./ 1000000;
        % Ts_start_stop{1,3} = Ts_test{1,1}(:,[1,3])./ 1000000;
        % Ts_start_stop{1,4} = Ts_posttest{1,1}(:,[1,3])./ 1000000;
        
        % step 2 load trackdata
        % load(trackdata{ns})
        
        
        % open TTlist
        celllist = TTList0;
        fid=fopen(celllist);
        if (fid == -1)
            warning([ 'Could not open tfile ' celllist]);
        else
            % read the file names from the t-file list
            TT0 = ReadFileList(celllist);
            numCells0 = length(TT0);
            if numCells0==1 && max(TT0{1}==-1)
                % no cells in ttlist
                numcells=0;
            else
                numcells=numCells0;
            end
        end
        
        % create the CSC file list
        CSC0 = cell(Ncell,1);
        for ii = 1:Ncell
            TTind = strfind(  TT0{ii}, '_' );
            CSC0{ii,1} = ['CSC' TT0{ii}(3:TTind-1) '.ncs'] ;
        end
        nseg = 1;
        
        spikePos = cell(Ncell,1);
        spikePos_alllap = cell(Ncell,1,1);
        spikePos_singlap = cell(Ncell,1,1);
        
        fprintf(['=====session %.3g D=%.3g start=====' '\n'], ns,D)
        for nc = 1:Ncell
            
            spikePos{nc,1} = spikes{nc,2};% spike 所在的posbin
            
            %% find the feature
            slap = 1;
            for nseg=1%%%only pre-running
                Ts_start_stop0 = Ts_start_stop{nseg};
                nlap = size(Ts_start_stop0,1);
                IND = [];
                for nl = 1:nlap
                    ind = find(spikes{nc,1}>=Ts_start_stop0(nl,1) & spikes{nc,1}<=Ts_start_stop0(nl,2));
                    spikePos_singlap{nc,1,slap} = spikePos{nc,1}(ind);
                    IND = [IND; ind];
                    slap = slap + 1;
                end
                spikePos_alllap{nc,1,nseg} = spikePos{nc,1}(IND);
            end
            
            %         nb = length(sprintf([repmat('>' , 1, round((nc-1)/Ncell*10)) '%.2f'],(nc-1)/Ncell*100));
            %         per = sprintf([repmat('>' , 1, round(nc/Ncell*10)) '%.2f'],nc/Ncell*100);
            %         fprintf(1,[repmat('\b',1,nb+1) '%1$s%2$s'],    per, '%');
        end
        fprintf('\ndone!!\n');
        % save the waveform feature
        %     [currPath, name, ext] = fileparts(file_input2);
        
        file_output = [subfolder{D} 'data_phaselocking_TSlap_vel0_new.mat'];
        save(file_output,'spikePos_alllap','spikePos_singlap','spikePos','-append');
    end
end