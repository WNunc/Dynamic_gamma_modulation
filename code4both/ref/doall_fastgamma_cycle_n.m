% get the timepoint of each slow gamma cycle
clear
close all
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
tcolor = {'black','red'};
peak0 = 1; % firing rate threshold to remove place cells

% fileinput1 = ['sgamma_dominant_thetacyc_IND_std' num2str(n_std) '_v2.mat'];% sg dominant theta cycle ind
% fileinput2 = [];% good theta cycle in each lap
% fileinput3 = 'data_phaselocking_TSlap_vel0.mat';% theta cycle
% fileinput4 = 'Cells_allsegment_v1_vel_0.mat';% spike info
% fileinput5 = 'Cells_allsegment_v1_vel_5.mat';% ratemap info
fileinput6 = [];% decoded info
fileinput7 = [];% all theta cycle in each lap
fileinput8 = 'data_info_allphase_both.mat';
fileoutput1 = 'data_CycleList_fg.mat';% sg cycle timestamp
outputFolder = 'H:\neuralynx\phase precession rseult\';
%%
for ns = 1:isession
    path_ns = path{ns};
    disp(path_ns);
    cd(path_ns)
    dir_TT=strcat(path_ns,'TTList_dCA1_pyr.txt');
    F = ReadFileList(dir_TT);
    ff8dir = fullfile(outputFolder,path_ns(13:end));
    load([ff8dir,fileinput8])
    
    
    for nt = 1%1:3
        switch nt
            case 1
                laps=5;
            case 2
                laps=8;
            case 3
                laps=8;
        end
        
        ListC = cell(nt,1);
        
        for D = 1:2 % CW��CCW��������
            goodphase_ns = Seq_cutPhase{ns,D};%***
            case3 = num2str(goodphase_ns);% directories_allData_v0_allgood ��������λһ��
            % input folder
            ff6dir = [path_ns Directionfolder{D}];
            ff7dir = [Phasecutfolder{D} 'phase_' num2str(case3) '\'];
            fileinput7 = ['data_theta_seq_info_AllLap-ontrack_3cell_5spk' num2str(case3) '_cycmid_v5.mat'];
            load([ff7dir,fileinput7],'theta_info')% load the theta sequence Infomation
            ListC0=[];
            for nl = 1:laps
                fileinput6 = strcat('scores',num2str(nt)','-',num2str(nl),'-ontrack_v5.mat');
                load([ff6dir,fileinput6])
                
                numberCell = length(cind_ot);
                switch D
                    case 1
                        CSCd = CSCnum(1:numberCell);
                    case 2
                        CSCd = CSCnum(end-numberCell+1:end);
                end
                %%
                tbin = scores{nl,6};
                t_EEG_sample = scores{nl,10};
                info = theta_info{nl};
                L=size(info,1);
                Map=zeros(NC,L);
                ABC=cell(NC,L);
                for nc = 1:numberCell
                    EEG_sample_fg=[];
                    phase_fg=[];
                    yp_fg=[];
                    ind_fg=[];
                    tt_fg=[];
                    List=[];
                    EEG_sample_fg = scores{nl,16}(CSCd(nc),:);
                    phase_fg=scores{nl,19}(CSCd(nc),:);
                    ind_fg =findpeaks(EEG_sample_fg);
                    ind_fg = ind_fg.loc;
                    ind_fg = ind_fg';
                    yp_fg = EEG_sample_fg(ind_fg);
                    tt_fg=t_EEG_sample(ind_fg);
                    for i=1:L % theta cycle
                        List{i,1}=find(tt_fg > info{i,3}(1) & tt_fg < info{i,3}(2));
                        Peak{i,1}=yp_fg(List{i,1});
                        List{i,2}=length(Peak{i,1});
                        [~,List{i,3}]=max(Peak{i,1});
                        if (List{i,2}-List{i,3})>=2 & List{i,3}>2 & List{i,1}(1,1)~=1
                            ListC0{nt,1}{nc,nl}{1,i}=1;
                            List{i,5}=[List{i,1}(List{i,3}-2),List{i,1}(List{i,3}-1),List{i,1}(List{i,3}),List{i,1}(List{i,3}+1),List{i,1}(List{i,3}+2)];
                            ind_gc=ind_fg(1,List{i,5});
                            gc1=ind_gc(1);
                            while EEG_sample_fg(gc1-1)<EEG_sample_fg(gc1)
                                gc1=gc1-1;
                            end
                            gc2=ind_gc(2);
                            while EEG_sample_fg(gc2-1)<EEG_sample_fg(gc2)
                                gc2=gc2-1;
                            end
                            gc3=ind_gc(3);
                            while EEG_sample_fg(gc3-1)<EEG_sample_fg(gc3)
                                gc3=gc3-1;
                            end
                            gc4=ind_gc(4);
                            while EEG_sample_fg(gc4-1)<EEG_sample_fg(gc4)
                                gc4=gc4-1;
                            end
                            gc5=ind_gc(5);
                            while EEG_sample_fg(gc5-1)<EEG_sample_fg(gc5)
                                gc5=gc5-1;
                            end
                            gc6=ind_gc(5);
                            while EEG_sample_fg(gc6+1)<EEG_sample_fg(gc6)
                                gc6=gc6+1;
                            end
                            List{i,6}=[gc1,gc2,gc3,gc4,gc5,gc6];
                            for z=1:6
                                Df=[];
                                Df=abs(tbin-t_EEG_sample(List{i,6}(z)));
                                [~,List{i,8}(1,z)]=min(Df);
                            end
                            ListC0{nt,1}{nc,nl}{2,i}=tbin(List{i,8});
                        elseif List{i,3}==1 & length(List{i,1})>=3 & List{i,1}(1,1)~=1
                            ListC0{nt,1}{nc,nl}{1,i}=2;
                            List{i,5}=[List{i,1}(List{i,3}),List{i,1}(List{i,3}+1),List{i,1}(List{i,3}+2)];
                            ind_gc=ind_fg(1,List{i,5});
                            gc1=ind_gc(1);
                            gc2=ind_gc(2);
                            while EEG_sample_fg(gc2-1)<EEG_sample_fg(gc2)
                                gc2=gc2-1;
                            end
                            gc3=ind_gc(3);
                            while EEG_sample_fg(gc3-1)<EEG_sample_fg(gc3)
                                gc3=gc3-1;
                            end
                            gc4=ind_gc(3);
                            while EEG_sample_fg(gc4+1)<EEG_sample_fg(gc4)
                                gc4=gc4+1;
                            end
                            List{i,6}=[gc1,gc2,gc3,gc4];
                            for z=1:4
                                Df=[];
                                Df=abs(tbin-t_EEG_sample(List{i,6}(z)));
                                [~,List{i,8}(1,z)]=min(Df);
                            end
                            ListC0{nt,1}{nc,nl}{2,i}=tbin(List{i,8});
                        elseif (List{i,2}-List{i,3})==0 & length(List{i,1})>=3 & i~=L
                            ListC0{nt,1}{nc,nl}{1,i}=3;
                            List{i,5}=[List{i,1}(List{i,3}-2),List{i,1}(List{i,3}-1),List{i,1}(List{i,3})];
                            ind_gc=ind_fg(1,List{i,5});
                            gc1=ind_gc(1);
                            while EEG_sample_fg(gc1-1)<EEG_sample_fg(gc1)
                                gc1=gc1-1;
                            end
                            gc2=ind_gc(2);
                            while EEG_sample_fg(gc2-1)<EEG_sample_fg(gc2)
                                gc2=gc2-1;
                            end
                            gc3=ind_gc(3);
                            while EEG_sample_fg(gc3-1)<EEG_sample_fg(gc3)
                                gc3=gc3-1;
                            end
                            gc4=ind_gc(3);
                            List{i,6}=[gc1,gc2,gc3,gc4];
                            for z=1:4
                                Df=[];
                                Df=abs(tbin-t_EEG_sample(List{i,6}(z)));
                                [~,List{i,8}(1,z)]=min(Df);
                            end
                            ListC0{nt,1}{nc,nl}{2,i}=tbin(List{i,8});
                        else
                            ListC0{nt,1}{nc,nl}{1,i}=0;
                        end
                    end
                end
            end
            ListC{nt} = [ListC{nt};ListC0{nt}];
        end
    end
    file_output=fullfile(ff8dir,fileoutput1);
    save(file_output,'ListC');
end
