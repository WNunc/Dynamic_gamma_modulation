%Part1
% 按cell
clear
clc
directories_allData_v0_allgood
outputFolder = ['H:\neuralynx\phase precession rseult\'];
%%
CNEG2 = [];CNEG1 = [];CZERO = [];CPOS1 = [];CPOS2 = [];
for ns = [1,2,4:6,8,10:13,15,16]%1:isession
    cNeg2 = [];cNeg1 = [];cZero = [];cPos1 = [];cPos2 = [];
    path_ns = path{ns};
    outdir = fullfile(outputFolder,path_ns(13:end));
    cd(outdir)
    file_input1='data_info_allphase_both.mat';% phase data
    file_input2='data_CycleList_fg.mat';% slow gamma cycle
    file_input3='data_PeakFRate_singlelap.mat';% pk fire rate
    load(file_input1);
    load(file_input2);
    load(file_input3);
    M0{1} = M;
    Neg2=[];
    Neg1=[];
    Zero=[];
    Pos1=[];
    Pos2=[];
    cNeg2 = [];cNeg1 = [];cZero = [];cPos1 = [];cPos2 = [];
    for nt = 1%1:3
        switch nt
            case 1
                laps=5;
            case 2
                laps=8;
            case 3
                laps=8;
        end
        % part1
        Fg=[];
        a=1;
        b=1;
        c=1;
        d=1;
        e=1;
        Ind_c=find(MNP>=5);
        Lc=length(Ind_c);
        for nl = 1:laps
            CC = 1;%%cell计数
            OP{1} = SPK_in{nl}(:,[1,4,5]);
            OPG{1} = SPK_in{nl}(:,[1,2,5]);
            for nlc = 1:Lc
                nc = Ind_c(nlc);
                Tspk=[];
                if M0{nt,1}(nc,nl)~=0
                    for m = 1:M0{nt,1}(nc,nl)
                        Tspk=OPG{nt}{nc,3}(m);
                        L=size(ListC{nt,1}{nc,nl},2);
                        Fg{nc,nl}(nl,m)=0;
                        for n = 1:L
                            Tc=ListC{nt,1}{nc,nl}{2,n};
                            switch ListC{nt,1}{nc,nl}{1,n}
                                case 1
                                    if Tc(1)<Tspk & Tc(2)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+1;
                                        Neg2{nt,1}(1,a) = OP{nt}{nc,2}(m);
                                        Neg2{nt,1}(2,a) = OPG{nt}{nc,2}(m);
                                        a=a+1;
                                    elseif Tc(2)<Tspk & Tc(3)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+1.5;
                                        Neg1{nt,1}(1,b) = OP{nt}{nc,2}(m);
                                        Neg1{nt,1}(2,b) = OPG{nt}{nc,2}(m);
                                        b=b+1;
                                    elseif Tc(3)<Tspk & Tc(4)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+2;
                                        Zero{nt,1}(1,c) = OP{nt}{nc,2}(m);
                                        Zero{nt,1}(2,c) = OPG{nt}{nc,2}(m);
                                        c=c+1;
                                    elseif Tc(4)<Tspk & Tc(5)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+2.5;
                                        Pos1{nt,1}(1,d) = OP{nt}{nc,2}(m);
                                        Pos1{nt,1}(2,d) = OPG{nt}{nc,2}(m);
                                        d=d+1;
                                    elseif Tc(5)<Tspk & Tc(6)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+3;
                                        Pos2{nt,1}(1,e) = OP{nt}{nc,2}(m);
                                        Pos2{nt,1}(2,e) = OPG{nt}{nc,2}(m);
                                        e=e+1;
                                    else
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m);
                                    end
                                case 2
                                    if Tc(1)<Tspk & Tc(2)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+2;
                                        Zero{nt,1}(1,c) = OP{nt}{nc,2}(m);
                                        Zero{nt,1}(2,c) = OPG{nt}{nc,2}(m);
                                        c=c+1;
                                    elseif Tc(2)<Tspk & Tc(3)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+2.5;
                                        Pos1{nt,1}(1,d) = OP{nt}{nc,2}(m);
                                        Pos1{nt,1}(2,d) = OPG{nt}{nc,2}(m);
                                        d=d+1;
                                    elseif Tc(3)<Tspk & Tc(4)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+3;
                                        Pos2{nt,1}(1,e) = OP{nt}{nc,2}(m);
                                        Pos2{nt,1}(2,e) = OPG{nt}{nc,2}(m);
                                        e=e+1;
                                    else
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m);
                                    end
                                case 3
                                    if Tc(1)<Tspk & Tc(2)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+1;
                                        Neg2{nt,1}(1,a) = OP{nt}{nc,2}(m);
                                        Neg2{nt,1}(2,a) = OPG{nt}{nc,2}(m);
                                        a=a+1;
                                    elseif Tc(2)<Tspk & Tc(3)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+1.5;
                                        Neg1{nt,1}(1,b) = OP{nt}{nc,2}(m);
                                        Neg1{nt,1}(2,b) = OPG{nt}{nc,2}(m);
                                        b=b+1;
                                    elseif Tc(3)<Tspk & Tc(4)>Tspk
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m)+2;
                                        Zero{nt,1}(1,c) = OP{nt}{nc,2}(m);
                                        Zero{nt,1}(2,c) = OPG{nt}{nc,2}(m);
                                        c=c+1;
                                    else
                                        Fg{nc,nl}(nl,m)=Fg{nc,nl}(nl,m);
                                    end
                            end
                        end
                    end
                end
                
                if ~isempty(Neg2)
                    tempN2 = Neg2{nt,1}(2,:);
                else
                    tempN2 = Neg2;
                end
                
                if ~isempty(Neg1)
                    tempN1 = Neg1{nt,1}(2,:);
                else
                    tempN1 = Neg1;
                end
                
                if ~isempty(Zero)
                    tempZ = Zero{nt,1}(2,:);
                else
                    tempZ = Zero;
                end
                
                if ~isempty(Pos1)
                    tempP1 = Pos1{nt,1}(2,:);
                else
                    tempP1 = Pos1;
                end
                
                if ~isempty(Pos2)
                    tempP2 = Pos2{nt,1}(2,:);
                else
                    tempP2 = Pos2;
                end
                
                if CC == 1
                    cNeg2{nt}(CC,nl) = ifmean(tempN2);
                    cNeg1{nt}(CC,nl) = ifmean(tempN1);
                    cZero{nt}(CC,nl) = ifmean(tempZ);
                    cPos1{nt}(CC,nl) = ifmean(tempP1);
                    cPos2{nt}(CC,nl) = ifmean(tempP2);
                else
                    cNeg2{nt}(CC,nl) = ifmean(tempN2,aa,a);
                    cNeg1{nt}(CC,nl) = ifmean(tempN1,bb,b);
                    cZero{nt}(CC,nl) = ifmean(tempZ,cc,c);
                    cPos1{nt}(CC,nl) = ifmean(tempP1,dd,d);
                    cPos2{nt}(CC,nl) = ifmean(tempP2,ee,e);
                end
                aa = a;% neg2 ind
                bb = b;% neg1 ind
                cc = c;% zero ind
                dd = d;% pos1 ind
                ee = e;% pos2 ind
                CC = CC+1;
            end
%             scatterphs(cNeg2{nt}(:,nl),cNeg1{nt}(:,nl),cZero{nt}(:,nl),cPos1{nt}(:,nl),cPos2{nt}(:,nl))
        end
    end

    save([outdir,'ThetaFastR_5_cell.mat'],'cNeg2','cNeg1','cZero','cPos1','cPos2');
end

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


function [neg2,neg1,zero,pos1,pos2,Na,Nb,Nc,Nd,Ne] = scatterphs(cNeg2,cNeg1,cZero,cPos1,cPos2)
neg2 = cNeg2(~isnan(cNeg2));
neg1 = cNeg1(~isnan(cNeg1));
zero= cZero(~isnan(cZero));
pos1 = cPos1(~isnan(cPos1));
pos2 = cPos2(~isnan(cPos2));

Na = ones(1,length(neg2))*-2;
Nb = ones(1,length(neg1))*-1;
Nc = ones(1,length(zero))*0;
Nd = ones(1,length(pos1))*1;
Ne = ones(1,length(pos2))*2;

figure
scatter(Na,neg2)
hold on
scatter(Nb,neg1)
scatter(Nc,zero)
scatter(Nd,pos1)
scatter(Ne,pos2)
hold off
xlabel('sg cycle')
ylabel('sg phase')
xlim([-2.5,2.5])
end
