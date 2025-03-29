% edited on 07/17/2017
% modified on 12/13/2017
% get all the theta windows
% modified from doall_add_thetawin_cz_v1
% use EEG from CSClist, but not TTList
% 可在调用程序中插入ind_tet= 任意值 替代本程序的输出
clear

directories_allData_v2

TTList0 = 'TTList_dCA1_pyr.txt';

for ns = 1
    path_ns = path{ns};
    cd(path_ns);
    
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};
    csclist_ns = CSClist_CA1{ns};
    
    tets = gettetnumbers(TTList0);
    if(min(ismember(csclist_ns,tets)) == 0)
        error('A tetrode in CSCList does not pick up cells!')
        break
    end
    
    if ~isempty(csclist_ns)
        %% cut theta cycles on theta peaks
        disp('Getting theta windows');
        
        % use EEG from csclist_ns
        [windows,ind_tet,name_tet] = thetawindows_all_CT_v2_EEG(csclist_ns);
        Code = 'thetawindows_all_CT_v2_EEG.m';
        
        file_output = 'thetawindows_thetapeak.mat';
        cd(path_ns);
        save(file_output,'windows','ind_tet','name_tet','Code');
        
        %% cut theta cycles on minimum spike number
        disp('Getting theta windows');
        
        % use EEG from csclist_ns
        [windows,ind_tet,name_tet]= thetawindows_all_CT_v3_EEG(csclist_ns,TTList0,path_ns);
        Code = 'thetawindows_all_CT_v3_EEG';
        
        file_output = 'thetawindows_minspk.mat';
        cd(path_ns);
        save(file_output,'windows','ind_tet','name_tet','Code');
        
        %% cut theta cycles on minimum spike number,
        % by using the tetrode with most number of cells
        disp('Getting theta windows');
        
        % use EEG from csclist_ns
        [windows,ind_tet,name_tet]= thetawindows_all_CT_v4_EEG(csclist_ns,TTList0,path_ns);
        Code = 'thetawindows_all_CT_v4_EEG';
        
        file_output = 'thetawindows_minspk_v2.mat';
        cd(path_ns);
        save(file_output,'windows','ind_tet','name_tet','Code');
        
        %% cut theta cycles on global maximum spike number,
        % by using the tetrode with most number of cells
        disp('Getting theta windows');
        
        % use EEG from csclist_ns
        load('G:\Colgin DATA\Circular Track\RawThetaPhs_AllRunSpks.mat','cutph_highFR','binsize')
        [windows,ind_tet,name_tet]= thetawindows_all_CT_v5_EEG(csclist_ns,TTList0,path_ns,cutph_highFR,binsize);
        Code = 'thetawindows_all_CT_v5_EEG';
        
        file_output = 'thetawindows_global_maxspk_v2.mat';
        cd(path_ns);
        save(file_output,'windows','ind_tet','name_tet','Code');
        
        %% cut theta cycles on maximum spike number,
        % by using the tetrode with most number of cells
        disp('Getting theta windows');
        
        % use EEG from csclist_ns
        [windows,ind_tet,name_tet]= thetawindows_all_CT_v6_EEG(csclist_ns,TTList0,path_ns);
        Code = 'thetawindows_all_CT_v6_EEG';
        
        file_output = 'thetawindows_maxspk_v2.mat';
        cd(path_ns);
        save(file_output,'windows','ind_tet','name_tet','Code');
    end
    %%
    cd ../
    clearvars windows ind_tet name_tet Code
    close all
end