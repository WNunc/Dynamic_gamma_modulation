% 检测theta sequence的显著性
shuffletimes=1500;
signifi=[]; dis=[];
mkdir(['H:\Lap-' num2str(nl)]);
tic;
for Nseq = 1:length(ind_ok)
    ind_thc = theta_info{nl}{Nseq,4};
    realpos_nseq = xbins(theta_info{nl}{Nseq,7});
    Pxn_nseq = lap_pxn(:,ind_thc(1):ind_thc(2));
    bin2use = find(~isnan(sum(Pxn_nseq)));
    tbins = lap_ts(ind_thc(1):ind_thc(2));
    tbins = tbins - tbins(1);
    if isempty(bin2use)
        r2 = 0;
        slope = 0;
        dis_nseq = NaN;
        signifi_nseq = 0;
    else
        tbins = lap_ts(ind_thc(1):ind_thc(2));
        tbins = tbins - tbins(1);
        [r2,~,~,~,slope] = Cir_reg(Pxn_nseq,xbins,tbins,bin2use);
        %shuffle
        pxn = Pxn_nseq;
        parfor ii = 1:shuffletimes %混洗1500次
            shuff = randperm(size(pxn,2));%生成打乱的随机序列
            pxn_shuff = pxn(:,shuff);%生成时间尺度打乱的pxn矩阵
            bin2use_shuff = find(~isnan(sum(pxn_shuff,1)));
            [r_shuffle(ii),~,~,~,~] = Cir_reg(pxn_shuff,xbins,tbins,bin2use_shuff);
        end
        r_shuffle(isnan(r_shuffle))=0;%将NaN置0
        r_distribution = sort(r_shuffle);%对shuffle后的值进行排序以便后面取分位数
        [~,dis_nseq]=min(abs(r_distribution-r2));
        if r2>r_distribution(shuffletimes*0.95)
            signifi_nseq=1;
        else
            signifi_nseq = 0;
        end
    end
    signifi(Nseq,1)= signifi_nseq;
    dis(Nseq,1)=dis_nseq;
    
    figure(1)
    imagesc(tbins,xbins,Pxn_nseq);
    hold on
    plot(tbins,realpos_nseq,'w')
    hold off
    caxis([0.005,0.12])
    colorbar
    colormap jet
    axis xy
    set(gca,'FontSize',18)
    title(['ThSeq_N_u_m = ' num2str(Nseq)])
    saveas(gcf,['H:\Lap-' num2str(nl) '\' num2str(Nseq) '.png'])
    
end
toc;