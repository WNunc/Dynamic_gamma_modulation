%%
clear
close all
directories_allData_v2
subfolder = 'Tseq-CCW\';
case1 = '-ontrack'; %
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
for ns = 1:isession
    path_ns = path{ns};
    disp(path_ns)
    trackdata_ns = trackdata{ns};
    goodphase_ns = Seq_cutPhase{ns,2};%CCW是第2列
    case3 = num2str(goodphase_ns);
    close all
    nseg = 1;
    file_input1 = [path_ns,subfolder,'data_gamma_pow_info_AllLap',case1, case2, case3, midmod,'_v5.mat'];
    load(file_input1,'gpower_seq')
    % power and its zscore
    for nl = 1:5
        file_input2 = [path_ns,subfolder,'data_theta_seq_info',case1,'_lap',num2str(nl),case2,case3,midmod,'_v5.mat'];
        load(file_input2,'thetaind')
        fprintf(1,'input file:\n%s\n',file_input2)
%         sgpower_lap(ns,nl) = mean(gpower_seq(thetaind,1));
        fgpower_lap(ns,nl) = mean(gpower_seq{nl}(thetaind,2));
%         sgpowerZ_lap(ns,nl) = mean(gpower_seq(thetaind,3));
        fgpowerZ_lap(ns,nl) = mean(gpower_seq{nl}(thetaind,4));
    end
end