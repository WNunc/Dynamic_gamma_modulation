function Xz = nanzscore_cz(X,flag,dim)

if dim == 2
    X = X';
end

Xz = nan(size(X));
for col = 1:size(X,2)
    x = X(:,col);
    x_mean = nanmean(x);
    if flag == 0
        x_std = nanstd(x);
    elseif flag == 1
        x_std = nanstd(x,1);
    end
    xz = (x-x_mean)/x_std;
    Xz(:,col) = xz;
end

if dim == 2
    Xz = Xz';
end