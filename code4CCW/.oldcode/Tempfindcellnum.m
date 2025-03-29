clear

close all
% directories_allData_v2
directories_allData_control
for ns = 1:isession
    path_ns = path{ns};
    
    TTlist0 = [path_ns 'TTList_dCA1_pyr.txt'];
    cellresult{ns,1} = path_ns;
    ncell = getnumberofcells_cz_v1(TTlist0);
    cellresult{ns,2} = ncell;
end