% doall_getInfo_sgamma_v3_2 画图子程序
% 单个trial中的单个cell在所有sg高的theta画图
FN = 'sgphase_difference_v3_3spk.mat';
mkdir([outputFolder path_ns(13:end)])
if exist([outputFolder path_ns(13:end) FN],'file')
    load([outputFolder path_ns(13:end) FN])
else
    phasediff2 = [0 0 0 0 0 0];
    fgcell_L = zeros(100,1)+2;
end

tcolor = {'black','red'};
theN = [];%每个cell在每圈中有几个thetacyc
theN_g = [];%每个cell在每圈中有几个goodthetacyc
[NCell,NLap] = size(Cell_TheSG);
for nl = 1:NLap
    for nc = 1:NCell
        if ~isempty(Cell_TheSG{nc,nl})
            Cell_TheSG{nc,nl}(cellfun(@isempty,Cell_TheSG{nc,nl}(:,1)),:)=[];
        end
        theN(nc,nl) = size(Cell_TheSG{nc,nl},1);
        
    end
end
mtheN = max(theN,[],2);
for nc = 1:NCell
    close all
    mt = mtheN(nc);
    pw = 0.98/mt;%plotwidth 
    ph = 0.96/(2*NLap);%plothight
    if mt>0
        cellname = [' cell-' num2str(nc)];
        FFA = figure('Name',cellname,'Color','w');%,'Visible','off');
        set(FFA,'Position',[0 50 3000 900])
    else
        continue
    end
    if exist('Pool_cellID','var')
        
        cellID = nc;
        if cellID <= L_D1
            dd = 1;% CW
            cellotID = cind_ot1(cellID);
        else
            cellID = cellID - L_D1;
            dd = 2;% CCW
            cellotID = cind_ot2(cellID);
        end
        if ismember(cellotID,Pool_cellID{1,dd})
            cellname = [cellname '-fg'];
        end
    end
    for nl = 1:NLap
        if theN(nc,nl)==0
            continue
        end
        
        for nth = 1:theN(nc,nl)
            % single cell SGphase in each lap
            indth = Cell_TheSG{nc,nl}{nth,1};
            thseq = Cell_TheSG{nc,nl}{nth,2};
            thseq_s = Cell_TheSG{nc,nl}{nth,3};
            thseq_t = Cell_TheSG{nc,nl}{nth,4};
            thcyc_onset = Cell_TheSG{nc,nl}{nth,5}(1);
            thcyc_offset = Cell_TheSG{nc,nl}{nth,5}(end);
            % spike time/phase scatter size
            t = Cell_TheSG{nc,nl}{nth,6}(1,:);
            p = Cell_TheSG{nc,nl}{nth,6}(2,:);
            sz = Cell_TheSG{nc,nl}{nth,6}(3,:);
            % slow gamma time; sg+theta; sgphase;
            SGTT = Cell_TheSG{nc,nl}{nth,7}(1,:);
            SGwave = Cell_TheSG{nc,nl}{nth,7}(2,:);
            SGphase = Cell_TheSG{nc,nl}{nth,7}(3,:);
            % sg; theta
            LFP_sgth = Cell_TheSG{nc,nl}{nth,7}(4:5,:);
            % sg cycle
            t_sgcyc = Cell_TheSG{nc,nl}{nth,8};
            sgcycind = Cell_TheSG{nc,nl}{nth,9};
            c = Cell_TheSG{nc,nl}{nth,10};
            figname = Cell_TheSG{nc,nl}{nth,11};
            
            
            %subplot(15,mt,[nth+3*(nl-1)*mt,nth+mt+3*(nl-1)*mt])
            %subplot(10,mt,nth+2*(nl-1)*mt,'Position',[0.02+(nth-1)*0.04,0.98-(2*nl-1)*0.04,0.03,0.04])
            axes('Position',[0.02+(nth-1)*pw,0.99-(2*nl-1)*ph,pw*0.75,ph*0.9-0.03])
            imagesc(thseq_t,thseq_s,thseq);axis xy
            hold on
            cyc_x = Cell_TheSG{nc,nl}{nth,5}(2:end-1);
            cyc_y = repmat(thseq_s(end),1,length(cyc_x));
            stem(cyc_x,cyc_y,'--w','Linewidth',2,'Marker','none')
            hold off
            
            % ylabel('position(rad)')
            colormap(jet)
            xticklabels([])
            xlim([thcyc_onset,thcyc_offset])
            set(gca,'FontSize',8);
            yticks([0,6.28])
            yticklabels({0,'2π'})
            ff = find(phasediff2(:,2)==nc & phasediff2(:,3)==nl & phasediff2(:,5)==indth);
            if isempty(ff)
                TT = sprintf('L-%u tCyc-%u',nl, indth);
            else
                TT = sprintf('L-%u tCyc-%u\nPhsDiff=%.3f',nl, indth,phasediff2(ff,1));
            end
            title(TT,'color',tcolor{c})
            %subplot(15,mt,nth+2*mt+3*(nl-1)*mt)
            %subplot(10,mt,nth+1*mt+2*(nl-1)*mt,'Position',[0.02+(nth-1)*0.04,0.98-(2*nl)*0.04,0.03,0.033])
            axes('Position',[0.02+(nth-1)*pw,0.99-(2*nl)*ph,pw*0.75,ph*0.9])
            yyaxis left
            scatter(t-thcyc_onset,p,sz, 'r|','LineWidth',2);%SPIKE
            %             hold on
            %             stem(t_sgcyc,720*ones(1,length(t_sgcyc)),...
            %                 'r','LineWidth',1,'Marker','none','ShowBaseLine','off')
            %             hold off
            xlim([0,thcyc_offset-thcyc_onset])
            ylim([0 720])
            yticks([0,360,720])
            yticklabels({0,'2π','4π'})
            %             ylim([0 360])
            %             yticks([0,360])
            %             yticklabels([0,360])
            %             xticklabels([])
            % ylabel('gamma_s phase(deg)')
            set(gca,'FontSize',8);
            set(gca,'Ycolor','r')
            % title(num2str(sgcycind'))
            yyaxis right
            plot(SGTT-thcyc_onset,LFP_sgth(1,:),'LineWidth',1.2,'Color',[0.3,0.3,0.3])% 只画sg
            % plot(SGTT,LFP_sgth(2,:))% 只画th
            % plot(SGTT,SGwave)% 画sg+th
            maxy = get(gca,'YLim');
            yticklabels([])
            hold on
            stem(t_sgcyc-thcyc_onset,maxy(2)*ones(1,length(t_sgcyc)),...
                '-k','LineWidth',1,'Marker','none','ShowBaseLine','off')
            stem(t_sgcyc-thcyc_onset,maxy(1)*ones(1,length(t_sgcyc)),...
                '-k','LineWidth',1,'Marker','none','ShowBaseLine','off')
            hold off
            %xlim([thcyc_onset,thcyc_offset])
            xlim([0,thcyc_offset-thcyc_onset])
            % xlabel('time(s)')
            axis tight
            box off
            % ylabel('amp(μV)')
            set(gca,'FontSize',8,'LineWidth',1,'Ycolor','k');
            
            % pause
        end
    end
%     saveas(FFA,[outputFolder path_ns(13:end-1) cellname '.png'])
%     saveas(FFA,[outputFolder path_ns(13:end-1) cellname],'epsc')

    saveas(FFA,[outputFolder path_ns(13:end-1) cellname '-3cyc.png'])
    saveas(FFA,[outputFolder path_ns(13:end-1) cellname '-3cyc'],'epsc')
    clear FFA
end
clear Pool_cellID
