function [sgcycind] = get_sgphase(sgts, sgphase,thetacyc,spk)
%%
% input
% sgts = slow gamma 时间点 一维数组
% sgphase = slow gamma 相位 一维数组
% thetacyc = theta 周期开始和结束的时间
% spk = spike
% output
% sgcycind = 每个spike所在的sg周期

t_sgcyc = findpeaks(sgphase);
t_sgcyc = t_sgcyc.loc;
t_sgcyc = [thetacyc(1) sgts(t_sgcyc(1:end-1)) thetacyc(2)];
sgcycind = [];
for  cy = 1:length(t_sgcyc)-1
    ii = find(spk>t_sgcyc(cy) & spk<=t_sgcyc(cy+1));
    if ~isempty(ii)
        cind = ones(length(ii),1)*cy;
        sgcycind = [sgcycind;cind];
    end
end
end

