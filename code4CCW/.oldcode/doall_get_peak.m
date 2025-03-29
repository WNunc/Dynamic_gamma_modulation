clear
clc
directories_allData_v1

%%
for ns = 1:isession
    path_ns=path{ns};
    PA=[];
    phase_peak=[];
    for nt = 1:3
        clf
        file_input=strcat(path_ns,'phase',num2str(nt),'.mat');
        load(file_input)
        PHASE=PHASE(:);
        a=1;
        P=[];
        for i=1:length(PHASE)
            if PHASE(i,1)~=0
                P(a,1)=PHASE(i,1);
                a=a+1;
            end
        end
        pd = fitdist(P,'Kernel','Kernel','epanechnikov');
        xgrid = linspace(0,360,72)';
        pdfEst = pdf(pd,xgrid);
        L=line(xgrid,pdfEst);
        [~,z]=max(L.YData);
        phase_peak(nt,1)=round(L.XData(z));
        PA=[PA;P];
    end
    pda = fitdist(P,'Kernel','Kernel','epanechnikov');
    pdfEst = pdf(pda,xgrid);
    figure(1)
    La=line(xgrid,pdfEst);
    xlabel('Angle / бу','fontsize',12)
    ylabel('Probability','fontsize',12)
    [~,z]=max(La.YData);
    phase_peak(4,1)=round(La.XData(z));
    file_output=strcat(path_ns,'phasepeak.mat');
    save(file_output,'phase_peak');
end

%%
figure(2)
h=histogram(P,[1:12:372]);
xlim([0 400]);
ylim([0 110]);
xlabel('Angle / бу','fontsize',12)
ylabel('Counts','fontsize',12)