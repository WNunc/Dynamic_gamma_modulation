function [FlipCorr,Qsize ]= flipcorr(ThetaSquenceStructure,TPsweep,SPsweep)
% 反转后做相关系数 - 不反转做相关系数
% 此处显示详细说明
%%
TP = size(ThetaSquenceStructure,2);TPmid = 0.5*TP;
SP = size(ThetaSquenceStructure,1);SPmid = 0.5*SP;
Qsize = sprintf('_%utbins_%uxbins',TPsweep,SPsweep);
tp1 = [TPmid-TPsweep:TPmid];
tp2 = [TPmid+1:TPmid+1+TPsweep];
sp1 = [SPmid-SPsweep:SPmid];
sp2 = [SPmid+1:SPmid+1+SPsweep];
ThetaSquenceStructure = ThetaSquenceStructure([sp1,sp2],[tp1,tp2]);

xbin = size(ThetaSquenceStructure,2);
nFC = nan(1,xbin/2);FC = nan(1,xbin/2);
for x = 1:xbin/2
    a = ThetaSquenceStructure(:,x);
    b = flip(ThetaSquenceStructure(:,x),1);
    c = ThetaSquenceStructure(:,end-x+1);
nFC(x) = corr(a,c);
FC(x) = corr(b,c);
end
mean_nFC = mean(nFC,'omitnan');mean_FC = mean(FC,'omitnan');
FlipCorr = mean_FC - mean_nFC;
end