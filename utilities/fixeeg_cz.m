function [removinds_all,windows] = fixeeg_cz(e,gap,Fs)

if nargin < 2
    gap =1; % default window width is 1 second
end

% Fs=2000;

e =e';%×ªÖÃ
inds = find(e > 1.99e+03);
df = diff(inds);
inddf = [find(df>gap*Fs),length(df)+1]; % connect the windows if two windows are close (<1s in between)
if ~isempty(df)
    windows1 = [[inds(1),inds(inddf(1:end-1) +1)]' , inds(inddf)'];
else 
    windows1 = [];
end
if ~isempty(windows1)
    windows1(:,1) = windows1(:,1)-round(gap*Fs/2);
    windows1(:,2) = windows1(:,2)+round(gap*Fs/2);
    windows1(windows1(:,1)<1,1) = 1;
    windows1(windows1(:,2)>length(e),2) = length(e);
end

removinds1 = [];
for i = 1:size(windows1,1)
      removinds1 = [removinds1, windows1(i,1):windows1(i,2)];
end
removinds1 = sort(unique(removinds1));



% 
inds = find(e < -1.99e+03);
df = diff(inds);
inddf = [find(df>gap*Fs),length(df)+1]; % connect the windows if two windows are close (<1s in between)
% 
if ~isempty(df)
    windows2 = [[inds(1),inds(inddf(1:end-1) +1)]' , inds(inddf)'];
else 
    windows2 = [];
end
% 
if ~isempty(windows2) 
    windows2(:,1) = windows2(:,1)-round(gap*Fs/2);
    windows2(:,2) = windows2(:,2)+round(gap*Fs/2);
    windows2(windows2(:,1)<1,1) = 1;
    windows2(windows2(:,2)>length(e),2) = length(e);
end
% 
removinds2 = [];
for i = 1:size(windows2,1)
      removinds2 = [removinds2, windows2(i,1):windows2(i,2)];
end
removinds2 = sort(unique([removinds2,removinds1]));

% 

% combine he high cut-off and low cut-off windows together
removinds_all=sort(unique([removinds1,removinds2]));
df = diff(removinds_all);
inddf = [find(df>gap*Fs),length(df)+1]; % connect the windows if two windows are close (<1s in between)
if ~isempty(df)
    windows = [[removinds_all(1),removinds_all(inddf(1:end-1) +1)]' , removinds_all(inddf)'];
else 
    windows = [];
end

