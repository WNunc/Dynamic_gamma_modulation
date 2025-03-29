function [WeightCorr, Qsize]= weightcorr(ThetaSquenceStructure,TPsweep,SPsweep)
% Feng T, et al.2015 weighted_correlation
% 参考了Guo

%   此处显示详细说明
%%
Qsize = sprintf('_%utbins_%uxbins',TPsweep,SPsweep);
TP = size(ThetaSquenceStructure,2);TPmid = 0.5*TP;
SP = size(ThetaSquenceStructure,1);SPmid = 0.5*SP;
tp1 = [TPmid-TPsweep:TPmid];
tp2 = [TPmid+1:TPmid+1+TPsweep];
sp1 = [SPmid-SPsweep:SPmid];
sp2 = [SPmid+1:SPmid+1+SPsweep];
ThetaSquenceStructure = ThetaSquenceStructure([sp1,sp2],[tp1,tp2]);

t = -TPsweep-0.5:1:TPsweep+0.5;t = t';
p = -SPsweep-0.5:1:SPsweep+0.5;

a = sum(ThetaSquenceStructure(:),'omitnan');
b = ThetaSquenceStructure * t ;
c = p * ThetaSquenceStructure;
mT = sum(b)/a;
mP = sum(c)/a;
T = (t-mT)';
P = (p-mP)';
for i = 1:(2*SPsweep+2)
    for j = 1:(2*TPsweep+2)
        x2(i,j) = ThetaSquenceStructure(i,j) * P(i,1) * T(1,j);
        y2(i,j) = ThetaSquenceStructure(i,j) * P(i,1)^2;
        z2(i,j) = ThetaSquenceStructure(i,j) * T(1,j)^2;
    end
end
X2 = sum(x2(:),'omitnan')/a;
Y2 = sum(y2(:),'omitnan')/a;
Z2 = sum(z2(:),'omitnan')/a;

WeightCorr = X2/((Y2*Z2)^(1/2));
end

