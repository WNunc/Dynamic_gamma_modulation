function [WeightCorr, Qsize]= nanweightcorr(ThetaSquenceStructure,TPsweep,SPsweep)
% Feng T, et al.2015 weighted_correlation
% 参考了Guo

%   此处显示详细说明
%%
Qsize = sprintf('_%utbins_%uxbins',TPsweep,SPsweep);
TP = size(ThetaSquenceStructure,2);TPmid = fix(0.5*TP);
SP = size(ThetaSquenceStructure,1);SPmid = fix(0.5*SP);
tp1 = [TPmid-TPsweep:TPmid];
tp2 = [TPmid+1:TPmid+1+TPsweep];
if tp1(1)==0
    tp1 = tp1(2:end);
end
sp1 = [SPmid-SPsweep:SPmid];
sp2 = [SPmid+1:SPmid+1+SPsweep];
ThetaSquenceStructure = ThetaSquenceStructure([sp1,sp2],[tp1,tp2]);
TSS = sum(ThetaSquenceStructure);
t_ind = ~isnan(TSS);
if isempty(find(t_ind,1)) || length(t_ind(t_ind==1))==1
    WeightCorr = 0;
    return
end
ThetaSquenceStructure = ThetaSquenceStructure(:,t_ind);
t = -TPsweep-0.5:1:TPsweep+0.5;t = t';t = t(t_ind);
p = -SPsweep-0.5:1:SPsweep+0.5;

a = sum(ThetaSquenceStructure(:));
b = ThetaSquenceStructure * t ;
c = p * ThetaSquenceStructure;
mT = sum(b)/a;
mP = sum(c)/a;
T = (t-mT)';
P = (p-mP)';
for i = 1:length(P)
    for j = 1:length(T)
        x2(i,j) = ThetaSquenceStructure(i,j) * P(i,1) * T(1,j);
        y2(i,j) = ThetaSquenceStructure(i,j) * P(i,1)^2;
        z2(i,j) = ThetaSquenceStructure(i,j) * T(1,j)^2;
    end
end
X2 = sum(x2(:))/a;
Y2 = sum(y2(:))/a;
Z2 = sum(z2(:))/a;

WeightCorr = X2/((Y2*Z2)^(1/2));
end

