clear
clc
directories_Data_Ratall % change according to experiment
%%
for ns = 1
    path_ns=path{ns};
    nt = 1;
    nl = 1;
    file_input=strcat(path_ns,'info',num2str(nt),'-',num2str(nl),'.mat');
    load(file_input);
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
    Pxn=zeros(13,41);
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
    
    SG=0;
    FG=0;
    for n = 1:L
        k=ind(n,1);
        SG=SG+info{k,9};
        FG=FG+info{k,7};
    end
    SG=SG/L;
    FG=FG/L;
    %
    SGno=0;
    FGno=0;
    L=length(info);
    for m=1:L
        fgno(m,1)=info{m,8};
        sgno(m,1)=info{m,10};
    end
    L=length(ind);
    for n = 1:L
        k=ind(n,1);
        SGno=SGno+sgno(k,1);
        FGno=FGno+fgno(k,1);
    end
    SN=SGno/L;
    FN=FGno/L;
    dir1=strcat('Pow',num2str(nt),'-',num2str(nl),'.mat');
    save(dir1,'FG','FN','SG','SN')
end