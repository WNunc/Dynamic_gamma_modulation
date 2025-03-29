function [ProbDiff ,Qsize] = probdiff(ThetaSquenceStructure,TPsweep,SPsweep)
% Probability Differences
%   此处显示详细说明
%%
    TP = size(ThetaSquenceStructure,2);TPmid = 0.5*TP;
    SP = size(ThetaSquenceStructure,1);SPmid = 0.5*SP;
    Qsize = sprintf('_%utbins_%uxbins',TPsweep,SPsweep);
    tp1 = [TPmid-TPsweep:TPmid];
    tp2 = [TPmid+1:TPmid+1+TPsweep];
    sp1 = [SPmid-SPsweep:SPmid];
    sp2 = [SPmid+1:SPmid+1+SPsweep];
    
    q24 = sum(sum(ThetaSquenceStructure(sp1,tp1),'omitnan')...
        +sum(ThetaSquenceStructure(sp2,tp2),'omitnan'),'omitnan');
    q13 = sum(sum(ThetaSquenceStructure(sp1,tp2),'omitnan')...
        +sum(ThetaSquenceStructure(sp2,tp1),'omitnan'),'omitnan');
    ProbDiff  = (q24 - q13)/(q24 + q13);
end

