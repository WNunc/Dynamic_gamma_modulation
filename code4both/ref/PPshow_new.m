clear
clc
directories_Data_Ratall % ����directories�ļ�

% ����ʵ������ѭ��
for ns = 1:isession
    dir_T=trackdata{ns,1};
    date_p=date{ns};
    path_p=path{ns};
    d=date_p(1:length(date_p)-1);
    file_input=strcat(path_p,'PP2017_v5.mat'); % ����theta��λ�������
    load(dir_T,'ang_test_reward_ontrack','ang_sample_reward_ontrack','sign_correct_test') % ����sample��test�Ľ���λ��
    load(file_input)
    % ����cellѭ����NCΪϸ����
    for n = 1:NC
        m=1;
        x=[];
        y=[];
        % Part1 ������ͼ
        for nt = 1
            % ����pre/sample/test�޸�Ȧ��
            switch nt
                case 1
                    laps=10;
                case 2
                    laps=8;
                case 3
                    laps=8;
            end
            % ѭ��Ŀ�ģ���������Ȧ�ŵ絽x��y
            for nl=1:laps
                if ~isnan(OP{nt,1}{n,nl})
                    L=size(OP{nt,1}{n,nl},2);
                    for l=1:L
                        x(1,m)=OP{nt,1}{n,nl}(nl,l);
                        y(1,m)=OP{nt,2}{n,nl}(nl,l);
                        m=m+1;
                    end
                else
                    continue
                end
            end
        end
        % ��ͼ����
        if ~isnan(x)
            clf % ���
            fig=figure(1);
            [a,~,p] = circ_lin_regress(x, y, 2); % �������ֱ��
            if a(1,1)>0 % �ж�ʹ�����ֱ���ں��ʵ�λ��
                k=0;
            else
                k=2;
            end
            SLA=[a,p];
            switch nt % ����Ȧ�����Ų���
                case 1
                    subplot(3,5,3);
                case 2
                    subplot(3,4,2);
                case 3
                    subplot(3,4,2);
            end   
            plot([x x],[y/180*pi y/180*pi+2*pi],'r.','MarkerSize',9); % ��ͼ
            axis([0 2*pi,0 4*pi]);
            tickx=0:pi/2:2*pi;
            ticky=0:pi:4*pi;
            set(gca, 'XTick',tickx);
            set(gca, 'XTickLabel',{'0','1/2pi','pi','3/2pi','2pi'});
            set(gca, 'YTick',ticky);
            set(gca, 'YTickLabel',{'0','pi','2pi','3pi','4pi'});
            xlabel('Position / rad')
            ylabel('Phase / ��')
            hold on
            SP=[SLA(1)*2*pi,SLA(2)+k*2*pi]; % б�ʡ���ʼ��λ
            refline(SP) % �����ֱ��
            hold off
            % get R
            sum_c=0;
            sum_s=0;
            for i=1:length(y)
                y1=SP(1)*x(i);
                sum_c=sum_c+cos(y(i)/180*pi-y1);
                sum_s=sum_s+sin(y(i)/180*pi-y1);
            end
            R=sqrt((sum_c/length(y))^2+(sum_s/length(y))^2); % �������ϵ��R
            % titles
            TL1=F{n,1};
            lt=length(TL1);
            TL1=TL1(1:lt-2);
            TL1=strrep(TL1,'_','-');
            TL2=strcat('a0=',num2str(SLA(1)));
            TL=strcat(TL1,'��',TL2);
            TL3=strcat('p=',num2str(SLA(3)),'��R=',num2str(R));
            title({TL;TL3},'FontSize',10);
            ylim([0 4*pi]);
        else
            continue
        end
        % ntӦ������һ�� ��ÿһȦ����λ����
        for nt = 1
            switch nt
                case 1
                    laps=10;
                case 2
                    laps=8;
                case 3
                    laps=8;
            end
            for nl = 1:laps
                switch nt
                    case 1
                        subplot(3,5,5+nl)
                    case 2
                        subplot(3,4,4+nl)
                    case 3
                        subplot(3,4,4+nl)
                end
                X=[];
                Y=[];
                if ~isnan(OP{nt,1}{n,nl})
                    X=OP{nt,1}{n,nl}(nl,:);
                    Y=OP{nt,2}{n,nl}(nl,:);
                    plot([X X],[Y/180*pi (Y+360)/180*pi],'r.','MarkerSize',9);
                    axis([0 2*pi,0 4*pi]);
                    tickx=0:pi/2:2*pi;
                    ticky=0:pi:4*pi;
                    set(gca, 'XTick',tickx);
                    set(gca, 'XTickLabel',{'0','1/2pi','pi','3/2pi','2pi'});
                    set(gca, 'YTick',ticky);
                    set(gca, 'YTickLabel',{'0','pi','2pi','3pi','4pi'});
                    xlabel('Position / rad')
                    ylabel('Phase / ��')
                    hold on
                    refline_dash(SP)
                    switch nt
                        case 2
                            stem(ang_sample_reward_ontrack(nl,1),4*pi,'Marker','none','Color','m')
                        case 3
                            stem(ang_test_reward_ontrack(nl,1),4*pi,'Marker','none','Color','m')
                    end
                    hold off
                    % get R
                    sum_c=0;
                    sum_s=0;
                    for i=1:size(OP{nt,2}{n,nl},2)
                        y1=SP(1)*OP{nt,1}{n,nl}(nl,i);
                        sum_c=sum_c+cos(OP{nt,2}{n,nl}(nl,i)/180*pi-y1);
                        sum_s=sum_s+sin(OP{nt,2}{n,nl}(nl,i)/180*pi-y1);
                    end
                    R=sqrt((sum_c/size(OP{nt,2}{n,nl},2))^2+(sum_s/size(OP{nt,2}{n,nl},2))^2);
                else
                    continue
                end
                TL1=strcat('lap',num2str(nl));
                TL3=strcat('R=',num2str(R));
                % �ж�testʱ������
                switch nt
                    case 3
                        if sign_correct_test(nl,1)==1
                            TL3=strcat('R=',num2str(R),',correct');
                        else
                            TL3=strcat('R=',num2str(R),',false');
                        end
                end
                title({TL1;TL3},'FontSize',10);
                ylim([0 4*pi]);
            end
        end
%         % �������ͼƬ
%         frame=getframe(fig);
%         img=frame2im(frame);
%         file_output=strcat(num2str(d),'-',num2str(n),'.png');
%         imwrite(img,file_output);
    end
end