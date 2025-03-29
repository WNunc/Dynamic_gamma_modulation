% 使用统一的预处理程序获取数据
clear

directories_allData_v0
for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
    trackdata_ns = trackdata{ns};
    pathT = trackdata_ns;
    Data_getvideo2(pathT)
    Data_getspikes
end