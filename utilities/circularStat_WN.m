function plv_struct = circularStat_WN(spikePhase)
% created by 王宁
% input: spikePhase (每个spike对应的相位)
% output: plv_struct (相锁值结构体数组)
%               'meanAngle' 平均角度
%               'angDev' 角度方差
%               'vectorLength' 相锁值
%               'rayleighP' p值
%               'rayleighZ'
%               'cirDev1'
%
if ~isempty(spikePhase)
    
    [meanAngle,angDev,vectorLength,rayleighP,rayleighZ,cirDev] = circularStat(spikePhase);
    plv_struct = struct('meanAngle', meanAngle, 'angDev', angDev, 'vectorLength',vectorLength,'rayleighP',rayleighP,'rayleighZ',rayleighZ,'cirDev',cirDev);
else
%     plv_struct = struct('meanAngle', 0, 'angDev', 0, 'vectorLength',0,'rayleighP',0,'rayleighZ',0,'cirDev',0);% 可能改成nan比较好
    plv_struct = struct('meanAngle', nan, 'angDev', nan, 'vectorLength',nan,'rayleighP',nan,'rayleighZ',nan,'cirDev',nan);
end


end