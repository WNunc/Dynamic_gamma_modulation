
% origin file is phaselocking_feature_catch.m
%
% this file can read those result *.mat files
% extract and filt the LFP signal from *.ncs
% calculate the phase of LFP the spike occured
% and PLV
%
% output each folder PLV result and LFP phase when spikes occured
% spikePhase  cell numbers × frequency bands
%                       every cell array is 1 × the number of spikes
% spikePhase_alllap  cell numbers × frequency bands × segment numbers
%                       every cell array is 1 × the number of spikes
% spikePhase_singlap  cell numbers × frequency bands × lap numbers
%                       every cell array is 1 × the number of spikes
% plFeature  1 × frequency bands
% plFeature_alllap  segment numbers × frequency bands
% plFeature_singlap  lap numbers × frequency bands
%
% created by WN on 20220412
%
% ==========use timestamp at trial start and end==========

%%
clear
close all
% directories_allData_v2
directories_allData_v0
subfolder = 'Tseq\';

f1 = 65;
f2 = 100;
s1 = 25;
s2 = 45;
th1 = 4;
th2 = 12;

for ns = 1:isession
    
    path_ns = path{ns};
    cd(path_ns);
    disp(path_ns)
    
    % step 1 Load spike firing data
    file_input1 = [subfolder 'Cells_singleLap_v2_vel_0.mat'];  % use to get all spikes
    load(file_input1)
    
    file_input2 = trackdata{ns};
    load(file_input2,'Ts_prerunning','Ts_sample','Ts_test','Ts_posttest')
    % read timestamps
    Ts_start_stop{1,1} = Ts_prerunning{1,1}./ 1000000;
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
    
    spikePhase = cell(Ncell,1);
    spikePhase_alllap = cell(Ncell,1,1);
    spikePhase_singlap = cell(Ncell,1,1);
    plFeature = cell(1,3);
    plFeature_alllap = cell(1,3);
    plFeature_singlap = cell(1,3);
    
    Current = [];
    
    fprintf(['=====session %.3g start=====' '\n'], ns)
    for nc = 1:Ncell
        
        if nc == 1 || (nc>=2 && (~strcmp(Current, CSC0{nc,1})))
            Current = CSC0{nc,1};
            disp(Current)
            [sample,tt,tt_raw] = loadCSC_new_WN(Current); % sample unit:uv, tt and tt_raw unit: s
            
            [sigTheta,sigTheta_filtwts]= eegfilt(sample', 2000, th1, th2,length(sample'),[],0,'fir1',0); % pop-up window mode
            [sigGammaS,sigGammaS_filtwts]= eegfilt(sample', 2000, s1, s2,length(sample'),[],0,'fir1',0);
            [sigGammaF,sigGammaF_filtwts]= eegfilt(sample', 2000, f1, f2,length(sample'),[],0,'fir1',0);
            fprintf('\n')
            % 相位与 bayesGamma_CT_v6.m 保持一致
            phaseTheta = angle(hilbert(-sigTheta(1,:)))*180/pi+180;
            phaseGammaS = angle(hilbert(sigGammaS(1,:)))*180/pi+180;
            phaseGammaF = angle(hilbert(sigGammaF(1,:)))*180/pi+180;
            
            % other slow way to filt eeg
            % xf = fftbandpass(samples1',2000,63,65,100,102);
            
            spikePhase{nc,1} = phaseTheta(SpikeTStoEegIndex(spikes{nc,1},tt_raw,2000));% theta
            spikePhase{nc,2} = phaseGammaS(SpikeTStoEegIndex(spikes{nc,1},tt_raw,2000));% slow gamma
            spikePhase{nc,3} = phaseGammaF(SpikeTStoEegIndex(spikes{nc,1},tt_raw,2000));% fast gamma
            
        else if nc>=2 && strcmp(Current, CSC0{nc,1})
                disp(Current)
                spikePhase{nc,1} = phaseTheta(SpikeTStoEegIndex(spikes{nc,1},tt_raw,2000));% theta
                spikePhase{nc,2} = phaseGammaS(SpikeTStoEegIndex(spikes{nc,1},tt_raw,2000));% slow gamma
                spikePhase{nc,3} = phaseGammaF(SpikeTStoEegIndex(spikes{nc,1},tt_raw,2000));% fast gamma
            end
        end
        plFeature{1}(nc) = circularStat_WN(spikePhase{nc,1});
        plFeature{2}(nc) = circularStat_WN(spikePhase{nc,2});
        plFeature{3}(nc) = circularStat_WN(spikePhase{nc,3});
        
        
        %% find the feature
        slap = 1;
        for nseg=1%%%only pre-running
            Ts_start_stop0 = Ts_start_stop{nseg};
            nlap = size(Ts_start_stop0,1);
            IND = [];
            for nl = 1:nlap
                ind = find(spikes{nc,1}>=Ts_start_stop0(nl,1) & spikes{nc,1}<=Ts_start_stop0(nl,2));
                spikePhase_singlap{nc,1,slap} = spikePhase{nc,1}(ind);
                spikePhase_singlap{nc,2,slap} = spikePhase{nc,2}(ind);
                spikePhase_singlap{nc,3,slap} = spikePhase{nc,3}(ind);
                plFeature_singlap{slap,1}(nc) = circularStat_WN(spikePhase_singlap{nc,1,slap});
                plFeature_singlap{slap,2}(nc) = circularStat_WN(spikePhase_singlap{nc,2,slap});
                plFeature_singlap{slap,3}(nc) = circularStat_WN(spikePhase_singlap{nc,3,slap});
                IND = [IND; ind];
                slap = slap + 1;
            end
            
            spikePhase_alllap{nc,1,nseg} = spikePhase{nc,1}(IND);
            spikePhase_alllap{nc,2,nseg} = spikePhase{nc,2}(IND);
            spikePhase_alllap{nc,3,nseg} = spikePhase{nc,3}(IND);
            plFeature_alllap{nseg,1}(nc) = circularStat_WN(spikePhase_alllap{nc,1,nseg} );
            plFeature_alllap{nseg,2}(nc) = circularStat_WN(spikePhase_alllap{nc,2,nseg} );
            plFeature_alllap{nseg,3}(nc) = circularStat_WN(spikePhase_alllap{nc,3,nseg});
            
        end
        
        %         nb = length(sprintf([repmat('>' , 1, round((nc-1)/Ncell*10)) '%.2f'],(nc-1)/Ncell*100));
        %         per = sprintf([repmat('>' , 1, round(nc/Ncell*10)) '%.2f'],nc/Ncell*100);
        %         fprintf(1,[repmat('\b',1,nb+1) '%1$s%2$s'],    per, '%');
    end
    fprintf('\ndone!!\n');
    % save the waveform feature
    %     [currPath, name, ext] = fileparts(file_input2);
    
    file_output = [subfolder 'data_phaselocking_TSlap_vel0_new.mat'];
    save(file_output,'Ts_start_stop','spikePhase','spikePhase_alllap','spikePhase_singlap', 'plFeature', 'plFeature_alllap', 'plFeature_singlap');
end