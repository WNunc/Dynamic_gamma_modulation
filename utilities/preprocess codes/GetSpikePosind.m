function Ind = GetSpikePosind(ts,post)

N = length(ts);
Ind = zeros(N,1);
for ii = 1:N
    [~,ind] = min(abs(post-ts(ii)));
    Ind(ii) = ind(1);
end