%Part1
% 按cell
clear
clc
directories_allData_v0_allgood
outputFolder = ['H:\neuralynx\phase precession rseult\'];
%%
CNEG = []; CZERO = []; CPOS = [];
for ns = [1,2,4:6,8,10:13,15,16]%1:isession
    path_ns = path{ns};
    outdir = fullfile(outputFolder,path_ns(13:end));
    cd(outdir)
    file_input1='data_info_allphase_both.mat';% phase data
    file_input2='data_CycleList_sg.mat';% slow gamma cycle
    file_input3='data_PeakFRate_singlelap.mat';% pk fire rate
    load(file_input1);
    load(file_input2);
    load(file_input3);
    M0{1} = M;
    Neg=[];
    Zero=[];
    Pos=[];
    cNeg = [];cZero = [];cPos = [];
    K=[];
    for nt = 1%1:3
        switch nt
            case 1
                laps=5;
            case 2
                laps=8;
            case 3
                laps=8;
        end
        Sg=[];
        a=1;
        b=1;
        c=1;
        
        Ind_c=find(MNP>=5);
        Lc=length(Ind_c);
        for nl = 1:laps
            CC = 1;%%cell计数
            OP{1} = SPK_in{nl}(:,[1,4,5]);
            OPG{1} = SPK_in{nl}(:,[1,3,5]);
            for nlc = 1:Lc
                nc = Ind_c(nlc);
                Tspk=[];
                if M0{nt,1}(nc,nl)~=0
                    for m = 1:M0{nt,1}(nc,nl)
                        Tspk=OPG{nt}{nc,3}(m);
                        L=size(ListC{nt,1}{nc,nl},2);
                        Sg{nc,nl}(nl,m)=0;
                        for n = 1:L
                            Tc=ListC{nt,1}{nc,nl}{2,n};
                            switch ListC{nt,1}{nc,nl}{1,n}
                                case 1
                                    if Tc(1)<Tspk & Tc(2)>Tspk
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m)+1;
                                        Neg{nt,1}(1,a) = OP{nt}{nc,2}(m);
                                        Neg{nt,1}(2,a) = OPG{nt}{nc,2}(m);
                                        a=a+1;
                                    elseif Tc(2)<Tspk & Tc(3)>Tspk
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m)+2;
                                        Zero{nt,1}(1,b) = OP{nt}{nc,2}(m);
                                        Zero{nt,1}(2,b) = OPG{nt}{nc,2}(m);
                                        b=b+1;
                                    elseif Tc(3)<Tspk & Tc(4)>Tspk
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m)+3;
                                        Pos{nt,1}(1,c) = OP{nt}{nc,2}(m);
                                        Pos{nt,1}(2,c) = OPG{nt}{nc,2}(m);
                                        c=c+1;
                                    else
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m);
                                    end
                                case 2
                                    if Tc(1)<Tspk & Tc(2)>Tspk
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m)+2;
                                        Zero{nt,1}(1,b) = OP{nt}{nc,2}(m);
                                        Zero{nt,1}(2,b) = OPG{nt}{nc,2}(m);
                                        b=b+1;
                                    elseif Tc(2)<Tspk & Tc(3)>Tspk
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m)+3;
                                        Pos{nt,1}(1,c) = OP{nt}{nc,2}(m);
                                        Pos{nt,1}(2,c) = OPG{nt}{nc,2}(m);
                                        c=c+1;
                                    else
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m);
                                    end
                                case 3
                                    if Tc(1)<Tspk & Tc(2)>Tspk
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m)+1;
                                        Neg{nt,1}(1,a) = OP{nt}{nc,2}(m);
                                        Neg{nt,1}(2,a) = OPG{nt}{nc,2}(m);
                                        a=a+1;
                                    elseif Tc(2)<Tspk & Tc(3)>Tspk
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m)+2;
                                        Zero{nt,1}(1,b) = OP{nt}{nc,2}(m);
                                        Zero{nt,1}(2,b) = OPG{nt}{nc,2}(m);
                                        b=b+1;
                                    else
                                        Sg{nc,nl}(nl,m)=Sg{nc,nl}(nl,m);
                                    end
                            end
                        end
                    end
                end
                if ~isempty(Neg)
                    tempN = Neg{nt,1}(2,:);
                else
                    tempN = Neg;
                end
                if ~isempty(Zero)
                    tempZ = Zero{nt,1}(2,:);
                else
                    tempZ = Zero;
                end
                if ~isempty(Pos)
                    tempP = Pos{nt,1}(2,:);
                else
                    tempP = Pos;
                end
                
                if CC == 1
                    cNeg{nt}(CC,nl) = ifmean(tempN);
                    cZero{nt}(CC,nl) = ifmean(tempZ);
                    cPos{nt}(CC,nl) = ifmean(tempP);
                else
                    cNeg{nt}(CC,nl) = ifmean(tempN,aa,a);
                    cZero{nt}(CC,nl) = ifmean(tempZ,bb,b);
                    cPos{nt}(CC,nl) = ifmean(tempP,cc,c);
                end
                aa = a;% neg ind
                bb = b;% zero ind
                cc = c;% pos ind
                CC = CC+1;
            end
            %             K{nt,1}(nl,1)=size(Neg{nt,1},2);
            %             K{nt,2}(nl,1)=size(Zero{nt,1},2);
            %             K{nt,3}(nl,1)=size(Pos{nt,1},2);
        end
    end
%     CNEG = [CNEG;cNeg{nt}];
%     CZERO = [CZERO;cZero{nt}];
%     CPOS = [CPOS;cPos{nt}];
%     scatterphs(cNeg{nt},cZero{nt},cPos{nt});
    save([outdir,'ThetaSlowR_5_cell.mat'],'cNeg','cZero','cPos');

end

%%
% [N,Z,P,na,zb,pc] = scatterphs(CNEG,CZERO,CPOS);



function m = ifmean(NZP,aabbcc,abc)
% NZP means Neg Zero Pos
% abc means a b c
% aabbcc means aa bb cc
if nargin == 1
    if isempty(NZP)
        m = nan;
    else
        m = mean(NZP);
        %m = circ_mean(NZP);
    end
else
    if aabbcc == abc
        m = nan;
    else
        m = mean(NZP(aabbcc:abc-1));
        %m = circ_mean(NZP(bb+1:b));
    end
end
end


function [neg,zero,pos,Na,Nb,Nc] = scatterphs(cNeg,cZero,cPos)
neg = cNeg(~isnan(cNeg));
zero= cZero(~isnan(cZero));
pos = cPos(~isnan(cPos));

Na = ones(1,length(neg))*-1;
Nb = ones(1,length(zero))*0;
Nc = ones(1,length(pos))*1;

figure
scatter(Na,neg)
hold on
scatter(Nb,zero)
scatter(Nc,pos)
hold off
xlabel('sg cycle')
ylabel('sg phase')
xlim([-1.5,1.5])
end

