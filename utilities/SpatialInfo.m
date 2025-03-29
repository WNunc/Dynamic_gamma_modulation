function [SI1,SI2] = SpatialInfo(Ratemap,posang_ontrack,mapAxis)
%SpatialInfo 计算空间信息率——圆轨数据
%   input
%   ratemap 放电率变化
%   posang_ontrack 在轨道上上的位置
%   mapAxis 空间bin的edge
%   output
%   SI1 = bit/spike
%   SI2 = bit/second
Ratemap = Ratemap';
bin_num = length(Ratemap);
%% probablity on each position bin
count_pos(1) = length(find(posang_ontrack<=mapAxis(1) | posang_ontrack>=mapAxis(bin_num)));
for i = 2:bin_num
    count_pos(i) = length(find(posang_ontrack>=mapAxis(i-1) & posang_ontrack<=mapAxis(i)));
end
prob_pos= count_pos/length(posang_ontrack);
%% mean firerate
lambda = mean(Ratemap);

%% infomation per spike
info_per_spk = 0;
for i = 1:bin_num
    info = prob_pos(i)*(Ratemap(i)/lambda)*log2((Ratemap(i)/lambda));
    info_per_spk = info_per_spk+info;
end
SI1 = info_per_spk;
%% information per second
info_per_sec = 0;
for i = 1:bin_num
    info = prob_pos(i)*Ratemap(i)*log2((Ratemap(i)/lambda));
    info_per_sec = info_per_sec+info;
end
SI2 = info_per_sec;
end

