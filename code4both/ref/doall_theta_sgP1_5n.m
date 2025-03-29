%Part1
clear
clc
directories_allData_v0_allgood
outputFolder = ['H:\neuralynx\phase precession rseult\'];
%%
for ns = 1:isession
    path_ns = path{ns};
    outdir = fullfile(outputFolder,path_ns(13:end));
    cd(outdir)
    file_input1='data_info_allphase_both.mat';% phase data
%     file_input2='data_CycleList_sg.mat';% slow gamma cycle
%     file_input2='data_CycleList_sg_v2.mat';% slow gamma cycle
    file_input2='data_CycleList_sg_v3.mat';% slow gamma cycle
    file_input3='data_PeakFRate_singlelap.mat';% pk fire rate
    load(file_input1);
    load(file_input2);
    load(file_input3);
    M0{1} = M;
    Neg=[];
    Zero=[];
    Pos=[];
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
        Ind_c=find(MNP>=5);% 最低放电率的阈值
        Lc=length(Ind_c);
        for nl = 1%4:5%1:2%1:laps%
            OP{1} = SPK_in{nl}(:,[1,4,5]);% theta 相位
            OPG{1} = SPK_in{nl}(:,[1,3,5]);% slow gamma 相位
            for nlc = 1:Lc
                nc = Ind_c(nlc);
                Tspk=[];
                if M0{nt,1}(nc,nl)~=0
                    for m = 1:M0{nt,1}(nc,nl)
                        Tspk=OPG{nt}{nc,3}(m);
                        L=size(ListC{nt,1}{nc,nl},2);
                        Sg{nc,nl}(nl,m)=0;
                        for n = 1:L
                            % Tc=ListC{nt,1}{nc,nl}{2,n};
                            switch ListC{nt,1}{nc,nl}{1,n}
                                case 1
                                    Tc=ListC{nt,1}{nc,nl}{2,n};
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
                                    Tc=ListC{nt,1}{nc,nl}{2,n};
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
                                    Tc=ListC{nt,1}{nc,nl}{2,n};
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
            end
            %             K{nt,1}(nl,1)=size(Neg{nt,1},2);
            %             K{nt,2}(nl,1)=size(Zero{nt,1},2);
            %             K{nt,3}(nl,1)=size(Pos{nt,1},2);
        end
    end
    lapN = ['-first'];%['-late'];%['-early'];%[];%
    save([outdir,'ThetaSlowR_5' lapN '3_v3.mat'],'Neg','Zero','Pos','K');
end