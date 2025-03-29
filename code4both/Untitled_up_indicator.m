% 确定weight correlation


clear
close all
directories_allData_v0
resuletFolder = 'H:\neuralynx\gamma in sequence result';

ncells = 3;nspkikes = 5;
ncases = sprintf('_%ucell_%uspk', ncells, nspkikes);
midmod = '_cycmid';
lockat = {'firstlap','alllap','f2lap'};L = 2;
D = 1;
TPsweep = 7;SPsweep = 10;
Qsize = sprintf('_%utbins_%uxbins',TPsweep,SPsweep);
peak0 = 1;

nlaps = 5;
nsn = 0;
for ns = [1:22,24:29,31]%1:isession
    nsn = nsn + 1;
    path_ns = path{ns};
    disp(path_ns);cd(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    inFolder = outFolder;
    goodphase_ns = Seq_cutPhase{ns,D};%***
    case3 = num2str(goodphase_ns);
    
    file_input1 = [inFolder 'weight_corr_eachseq.mat'];
    load(file_input1)
    file_input2 =[inFolder 'Weight_Correlation' ncases Qsize case3 midmod lockat{L} '.mat'];
    load(file_input2)
    
    a = [];b = [];wc1_m = [];wc2_m = [];
    for nl = 1:nlaps
        S1 = size(wc1{nl});
        a = [a;[wc1{nl}(:,1),ones(S1(1),1)*nl]];
        b = [b;[wc2{nl}(:,1),ones(S1(1),1)*nl]];
        wc1_m(nl,:) = mean(wc1{nl});%WC1_sem(nl,:) = std(wc1{nl})./sqrt(S1(1));
        wc2_m(nl,:) = mean(wc2{nl});%WC2_sem(nl,:) = std(wc2{nl})./sqrt(S1(1));
    end
    %     %%%↓不shuffle直接算↓%%%
    %     A(:,1) = wc1_m(:,1);
    %     A(:,2) = [1:nlap]';
    %     [R_A,P_A] = corrcoef(A)
    %     B(:,1) = wc2_m(:,1);
    %     B(:,2) = [1:nlap]';
    %     [R_B,P_B] = corrcoef(B)
    %
    %     [R_a,P_a] = corrcoef(a)
    %     [R_b,P_b] = corrcoef(b)
    %
    %     R(ns,1) = R_A(1,2);R(ns,2) = R_B(1,2);
    %     r(ns,1) = R_a(1,2);r(ns,2) = R_b(1,2);
    %     P(ns,1) = P_A(1,2);P(ns,2) = P_B(1,2);
    %     p(ns,1) = P_a(1,2);p(ns,2) = P_b(1,2);
    %     %%%↑不shuffle直接算↑%%%
    
    %%%↓shuffle↓%%%
    nlap = 2;%nlaps
    x = 1:nlap;
    X = [ ones(nlap,1), x' ];
    %     [r1_real,~] = corr(wc1_m(:,1),x','type','Pearson');
    %     [r2_real,~] = corr(wc2_m(:,1),x','type','Pearson');
    %     [R1_real,~] = corr(WC1(:,1),x','type','Pearson');
    %     [R2_real,~] = corr(WC2(:,1),x','type','Pearson');
    
    [b1_real(nsn,:),~,~,~,stat1_real(nsn,:)]= regress(wc1_m(x,1),X);
    [b2_real(nsn,:),~,~,~,stat2_real(nsn,:)]= regress(wc2_m(x,1),X);
    [B1_real(nsn,:),~,~,~,Stat1_real(nsn,:)]= regress(WC1(x,1),X);
    [B2_real(nsn,:),~,~,~,Stat2_real(nsn,:)]= regress(WC2(x,1),X);
    
    if nlap == nlaps
    shuffletimes=100;
    for ish = 1:shuffletimes % 混洗100次
        sh = randperm(nlap);%生成打乱的随机序列
        % 单个sequence的wc平均
        wc1_mshuff = wc1_m(sh,1);
        wc2_mshuff = wc2_m(sh,1);
        % 平均sequence的wc
        WC1_shuff = WC1(sh,1);
        WC2_shuff = WC2(sh,1);
        
        %         [r1(ish),p1(ish)] = corr(wc1_mshuff,x','type','Pearson');
        %         [r2(ish),p1(ish)] = corr(wc2_mshuff,x','type','Pearson');
        %         [R1(ish),P1(ish)] = corr(WC1_shuff,x','type','Pearson');
        %         [R2(ish),P1(ish)] = corr(WC2_shuff,x','type','Pearson');
        
        [b1(ish,:),~,~,~,stat1(ish,:)]= regress(wc1_mshuff,X);
        [b2(ish,:),~,~,~,stat2(ish,:)] = regress(wc2_mshuff,X);
        [B1(ish,:),~,~,~,Stat1(ish,:)] = regress(WC1_shuff,X);
        [B2(ish,:),~,~,~,Stat2(ish,:)] = regress(WC2_shuff,X);
    end
    rr(:,1) = sort(stat1(:,1));rr(:,2) = sort(stat2(:,1));
    RR(:,1) = sort(Stat1(:,1));RR(:,2) = sort(Stat2(:,1));
    bb(:,1) = sort(b1(:,2));bb(:,2) = sort(b2(:,2));
    BB(:,1) = sort(B1(:,2));BB(:,2) = sort(B2(:,2));
    
    
    if stat1_real(1)>rr(shuffletimes*0.95,1)
        sig_r1(ns) = 1
    else
        sig_r1(ns) = 0
    end
    
    if stat2_real(1)>rr(shuffletimes*0.95,2)
        sig_r2(ns) = 1
    else
        sig_r2(ns) = 0
    end
    
    if Stat1_real(1)>RR(shuffletimes*0.95,1)
        sig_R1(ns) = 1
    else
        sig_R1(ns) = 0
    end
    
    if Stat2_real(1)>RR(shuffletimes*0.95,2)
        sig_R2(ns) = 1
    else
        sig_R2(ns) = 0
    end
    
    if b1_real(2)>bb(shuffletimes*0.95,1)
        sig_b1(ns) = 1
    else
        sig_b1(ns) = 0
    end
    
    if b2_real(2)>bb(shuffletimes*0.95,2)
        sig_b2(ns) = 1
    else
        sig_b2(ns) = 0
    end
    
    if B1_real(2)>BB(shuffletimes*0.95,1)
        sig_B1(ns) = 1
    else
        sig_B1(ns) = 0
    end
    
    if B2_real(2)>BB(shuffletimes*0.95,2)
        sig_B2(ns) = 1
    else
        sig_B2(ns) = 0
    end
    %%%↑shuffle↑%%%
    elseif nlap == 2
        aa = b1_real(:,2)>0
        bb = b2_real(:,2)>0
        cc = B1_real(:,2)>0
        dd = B2_real(:,2)>0
    end
end
% SIG = sig_r1 + sig_r2  + sig_b1 + sig_b2 %+ sig_R1 + sig_R2+ sig_B1 + sig_B2

SIG = aa + bb  + cc + dd %+ sig_R1 + sig_R2+ sig_B1 + sig_B2