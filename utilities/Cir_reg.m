function [r2,calphase,xaxis,p,slope,xspan,tspan] = Cir_reg(Pxn,xbins,tbins,bins2use)
% CZ: shuffle the position bins instead of time bins, which was uesd in
% Davison et al Neuron 2009 paper

%changed to use actual time and space dimensions rather than
%just bin numbers
%also added in com measeure and its difference from the predicted line as
%an additional measure of accuracy of prediction


%how many surrogate r2 values to get:
drawnum = 1000;

all = [];
s = RandStream('mt19937ar','Seed',1); 
for i = 1:size(bins2use,2) %passing through time bins
   bin = bins2use(i);
   pdtemp = Pxn(:,bin)';
   if isempty(find(isnan(pdtemp), 1))
       randlocs = xbins(  randsample(s,length(pdtemp),drawnum,true,pdtemp)  ); 
       randlocs(:,2) = tbins(bin);
       all = [all;randlocs];
   end
end

if ~isnan(all) & ~isempty(all)
    %% calculate the original sequence
    x = all(:,2);%ʱ��
    y = all(:,1);%λ�ã�0��2*pi��Բ��λ�ã�
    xaxis = sort(unique(x));
    [para,~,p] = circ_lin_regress_v2(x, y, 8);
    calphase = 2*pi*para(1,1)*xaxis+para(1,2);
    slope = para(1,1)*2*pi;
    
%     [phi0,para,R] =cl_regression(x,y,0,15) % wenbo
%     tang��ƪ��������ĺ����������������õı�����һ����
%      slope = para*2*pi;
%     calphase = slope*xaxis+phi0;
    
    
    xspan = slope * (x(end)-x(1));
    tspan = x(end)-x(1);
    
    over = calphase>=2*pi;
    under = calphase<0;
    calphase(over) = calphase(over)-2*pi;
    calphase(under) = calphase(under)+2*pi;
    
    % compute the residual(�вʵ��ֵ�����ֵ֮��Ĳ�)
    yresidual = zeros(drawnum,size(xaxis,1));
    ytotal = zeros(drawnum,size(xaxis,1));
    for ii = 1:size(xaxis,1)
        yfit = repmat(calphase(ii),drawnum,1);
         yresidual(:,ii) = circdistance(y(x==xaxis(ii)), yfit, 1);% residual(�вʵ��ֵ�����ֵ֮��Ĳ�)
         ytotal(:,ii) = circdistance(y(x==xaxis(ii)), circ_mean(y), 1); % ʵ��ֵ�� ʵ��ֵ��ƽ��ֵ֮��Ĳ�
    end
    SSresid = sum(sum(yresidual));
    SStotal = sum(sum(ytotal));
    r2 = 1 - SSresid/SStotal;
    if r2<0
        r2 = NaN;
    end
        
end    
end
