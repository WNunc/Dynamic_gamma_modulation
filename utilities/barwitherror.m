function [ba,er] = barwitherror(x,data)

datanum = size(data(data~=nan),1);

datamean = mean(data,'omitnan');
sem = std(data)/sqrt(datanum);

ba = bar(x,datamean);
hold on
er = errorbar(x,datamean,sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
box off
end

