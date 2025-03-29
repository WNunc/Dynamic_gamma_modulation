% 对应一下sg找出来的cell和fg相锁的cell
% 把doall_sgPhaseDiff_distribution得到的结果中的cellID与解码时找到的fgcellID
% 对应一下

clear
close all


n_std = 1.5;
n_spk = 3;
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
Dx = {'-cw','-ccw'};
case1 = '-ontrack';
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
nsn = 0;
tcolor = {'black','red'};
peak0 = 1; % firing rate threshold to remove place cells

fileinput1 = ['sgphase_difference_v3_' num2str(n_spk) 'spk.mat'];%
fileinput2 = [];% 排除掉的fgcell
fileinput3 = [];% 排除掉的nonfgcell
inputFolder1 = ['H:\neuralynx\sgamma result v3 std-' num2str(n_std) '\'];
inputFolder2 = [];% 根目录
fileoutput = fileinput1;
nlap = 5;
lockat = {'firstlap','alllap','f2lap'};L = 2;
PhaseDiff2 = [];LABEL2 = [];
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    path_ns = path{ns};
    disp(path_ns)
    outputFolder = [inputFolder1 path_ns(13:end)];
    cd(outputFolder)
    load([outputFolder fileinput1])
    % load fg & nonfg cellID form decoding
    for D = 1:2 % CW和CCW两个方向
        goodphase_ns = Seq_cutPhase{ns,D};%***
        case3 = num2str(goodphase_ns);% directories_allData_v0_allgood 俩方向相位一样
        % 排除神经元解码的 input folder
        inputFolder = [path_ns Directionfolder{D}];
        
        file_exfg = ['scores1-1-ontrack_exfg_TSlap_vel0',lockat{L},'v2.mat'];
        load([inputFolder file_exfg],'cind_ot','cind_ot_fg')
        if D == 1
            cind_ot1 = cind_ot;
            L_D1 = length(cind_ot);
        else
            cind_ot2 = cind_ot;
            L_D2 = length(cind_ot);
        end
        file_exnfg = ['scores1-1-ontrack_dsfg_TSlap_vel0',lockat{L},'v3.mat'];
        load([inputFolder file_exnfg],'cind_ot_nonfg')
        Pool_cellID{1,D} = cind_ot_fg;
        Pool_cellID{2,D} = cind_ot_nonfg;
    end
    
    
    % 遍历每个theta cycle
    for nth = 1:size(phasediff2,1)% phasediff2是按cell统计的结果
        cellID = phasediff2(nth,2);
        if cellID <= L_D1
            dd = 1;% CW
            cellotID = cind_ot1(cellID);
        else
            cellID = cellID - L_D1;
            dd = 2;% CCW
            cellotID = cind_ot2(cellID);
        end
        phasediff2(nth,7) = ismember(cellotID,Pool_cellID{1,dd});
        phasediff2(nth,8) = ismember(cellotID,Pool_cellID{2,dd});
    end
    
    PhaseDiff2 = [PhaseDiff2;phasediff2];
    if isempty(phasediff2)
        continue
    end
    [C,ia,~] = unique(phasediff2(:,2));
    label_fg = phasediff2(ia,:);
    label_fg = label_fg(label_fg(:,1)<0,7);% phasediff为负数的cell的标签，1=fgcell 0=nfgcell
    LABEL2 = [LABEL2;label_fg];
    fgcell_L = phasediff2(:,7);% 能算出来phasediff的cell的标签，1=fgcell 0=nfgcell
    save([outputFolder fileinput1],'fgcell_L','Pool_cellID',...
        'L_D1','cind_ot1','cind_ot2','phasediff2','-append')
end


ff1 = figure('Position', [487 253 890 420]);
subplot(1,2,1)
total_c1 = length(PhaseDiff2(:,end));
fg_per1 = length(find(PhaseDiff2(:,end)==1))/total_c1;
p = pie([fg_per1,1-fg_per1]);
p(1).FaceColor = '#F24444';
p(3).FaceColor = '#F2CA50';
% histogram(PhaseDiff2(:,end));
title(['包含重复的神经元 total cell = ' num2str(total_c1)],'FontSize',15);
% ylabel('Cell Num')
subplot(1,2,2)
total_c2 = length(LABEL2);
fg_per2 = length(find(LABEL2==1))/total_c2;
p = pie([fg_per2,1-fg_per2]);
p(1).FaceColor = '#F24444';
p(3).FaceColor = '#F2CA50';
% histogram(LABEL);
title(['不包含重复的神经元 total cell = ' num2str(total_c2)],'FontSize',15)
% ylabel('Cell Num')
legend({'fg cell','nfg cell'},'Location','northwest','FontSize',13)