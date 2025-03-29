function plv_struct = circularStat_WN(spikePhase)
% created by ����
% input: spikePhase (ÿ��spike��Ӧ����λ)
% output: plv_struct (����ֵ�ṹ������)
%               'meanAngle' ƽ���Ƕ�
%               'angDev' �Ƕȷ���
%               'vectorLength' ����ֵ
%               'rayleighP' pֵ
%               'rayleighZ'
%               'cirDev1'
%
if ~isempty(spikePhase)
    
    [meanAngle,angDev,vectorLength,rayleighP,rayleighZ,cirDev] = circularStat(spikePhase);
    plv_struct = struct('meanAngle', meanAngle, 'angDev', angDev, 'vectorLength',vectorLength,'rayleighP',rayleighP,'rayleighZ',rayleighZ,'cirDev',cirDev);
else
%     plv_struct = struct('meanAngle', 0, 'angDev', 0, 'vectorLength',0,'rayleighP',0,'rayleighZ',0,'cirDev',0);% ���ܸĳ�nan�ȽϺ�
    plv_struct = struct('meanAngle', nan, 'angDev', nan, 'vectorLength',nan,'rayleighP',nan,'rayleighZ',nan,'cirDev',nan);
end


end