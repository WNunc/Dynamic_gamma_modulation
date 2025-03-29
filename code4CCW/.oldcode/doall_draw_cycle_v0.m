clear
clc
directories_allData_v1
%%
nt = 1;
for nl = 3
    for ns = 1:isession
        dir0 = path{ns};
        dird = date{ns};
        dir=strcat(dir0,'info',num2str(nt),'-',num2str(nl),'.mat');
        load(dir);
        % select and find the Pxn we need
        L=length(info);
        for m=1:L
            ac(m,1)=info{m,4};
            spk(m,1)=info{m,5};
            x(m,1)=size(info{m,6},1);
            y(m,1)=size(info{m,6},2);
        end
        ind=find(ac>=3 & spk>=8 & x==13 & y==41);
        L=length(ind);
        if L==0
            continue
        else
            Pxn=zeros(13,41);
            %
            for x = 1:13
                for y = 1:41
                    for n = 1:L
                        k(1,n)=info{ind(n,1),6}(x,y);
                    end
                    p=mean(k,2,'omitnan');
                    Pxn(x,y)=p;
                    y=y+1;
                end
                x=x+1;
            end
            % Weighted Correlation
            l1=21;
            l2=7;
            r1=5;
            r2=3;
            pxn = Pxn(l2-r2:l2+r2,l1-r1:l1+r1);
            t=-r1:1:r1;
            p=-r2:1:r2;
            t = t';
            a = sum(pxn(:),'omitnan');
            b = pxn * t;
            c = p * pxn;
            mT = sum(b)/a;
            mP = sum(c)/a;
            T = (t-mT)';
            P = (p-mP)';
            for i = 1:(2*r2+1)
                for j = 1:(2*r1+1)
                    x1(i,j) = pxn(i,j) * P(i,1) * T(1,j);
                    y1(i,j) = pxn(i,j) * P(i,1)^2;
                    z1(i,j) = pxn(i,j) * T(1,j)^2;
                end
            end
            X = sum(x1(:),'omitnan')/a;
            Y = sum(y1(:),'omitnan')/a;
            Z = sum(z1(:),'omitnan')/a;
            Corr = X/((Y*Z)^(1/2));
            clearvars i j
            % Probability Differences
            Q3=sum(pxn(1:r2,1:r1));
            Q1=sum(pxn(r2+2:2*r2+1,r1+2:2*r1+1));
            Q2=sum(pxn(r2+2:2*r2+1,1:r1));
            Q4=sum(pxn(1:r2,r1+2:2*r1+1));
            Q1=sum(Q1);
            Q2=sum(Q2);
            Q3=sum(Q3);
            Q4=sum(Q4);
            PD=(Q1+Q3-Q2-Q4)/(Q1+Q2+Q3+Q4);
            % Theta Sequence Slope
            [~,indP]=max(pxn,[],1);
            for i = 1:length(indP)
                x(1,i)=t(i,1)*0.01;
                y(1,i)=p(indP(i))*pi/45;
            end
            [Slope,~,~] = circ_lin_regress(x, y, 2);
            G{ns,1} = Pxn;
            G{ns,2} = Corr;
            G{ns,3} = PD;
            G{ns,4} = Slope(1,1);
        end
        clearvars ac spk x y k spd x1 y1 z1
    end
    %
    for x = 1:13
        for y = 1:41
            for n = 1:size(G,1) % change
                if ~isempty(G{n,1})
                    k(1,n)=G{n,1}(x,y);
                else
                    continue
                end
            end
            p=mean(k,2,'omitnan');
            Pxn(x,y)=p;
            y=y+1;
        end
        x=x+1;
    end
end
%%
p=-6:1:6;
t=-20:1:20;
fig=figure(1);
uimagesc(t,p,smooth2a(normz(Pxn,0),1,1))
axis xy
xlabel('Time / ms','Fontsize',15)
set(gca,'XTick',-20:5:20);
set(gca,'XTickLabel',{'-200','-150','-100','-50','0','50','100','150','200'},'Fontsize',13);
ylabel('Location on track / rad','Fontsize',15)
set(gca,'YTick',-6:2:6);
set(gca,'YTickLabel',{'-6/45 pi','-4/45 pi','-2/45 pi','0','2/45 pi','4/45 pi','6/45 pi'},'Fontsize',13);
colormap(jet(128));
h=colorbar;
h.Label.String = 'Nomarlized Probability';
% caxis([0.01 0.13])
% %±£´æÍ¼Æ¬
% frame=getframe(fig);
% img=frame2im(frame);
% output=strcat('ThetaCycle-',num2str(nt),'-',num2str(nl),'.png');
% imwrite(img,output);
%% Weighted Correlation
l1=21;
l2=7;
r1=5;
r2=3;
pxn = Pxn(l2-r2:l2+r2,l1-r1:l1+r1);
t=-r1:1:r1;
p=-r2:1:r2;
t = t';
a = sum(pxn(:),'omitnan');
b = pxn * t;
c = p * pxn;
mT = sum(b)/a;
mP = sum(c)/a;
T = (t-mT)';
P = (p-mP)';
for i = 1:(2*r2+1)
    for j = 1:(2*r1+1)
        x2(i,j) = pxn(i,j) * P(i,1) * T(1,j);
        y2(i,j) = pxn(i,j) * P(i,1)^2;
        z2(i,j) = pxn(i,j) * T(1,j)^2;
    end
end
X2 = sum(x2(:),'omitnan')/a;
Y2 = sum(y2(:),'omitnan')/a;
Z2 = sum(z2(:),'omitnan')/a;
Corr = X2/((Y2*Z2)^(1/2));

%% Probability Differences
Q3=sum(pxn(1:r2,1:r1));
Q1=sum(pxn(r2+2:2*r2+1,r1+2:2*r1+1));
Q2=sum(pxn(r2+2:2*r2+1,1:r1));
Q4=sum(pxn(1:r2,r1+2:2*r1+1));
Q1=sum(Q1);
Q2=sum(Q2);
Q3=sum(Q3);
Q4=sum(Q4);
PD=(Q1+Q3-Q2-Q4)/(Q1+Q2+Q3+Q4);