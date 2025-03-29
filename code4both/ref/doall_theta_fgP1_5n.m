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
    file_input2='data_CycleList_fg.mat';% fast gamma cycle
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
        Ind_c=find(MNP>=1);
        Lc=length(Ind_c);
        for nl = 1:laps%4:5%
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
            end
        end
    end
    lapN = [];%['late'];%
    save([outdir,'ThetaFastR_5' lapN '3.mat'],'Neg2','Neg1','Zero','Pos1','Pos2');
end