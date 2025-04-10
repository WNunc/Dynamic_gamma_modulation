% Calculates the rate map.
function map = ratemap_ang_cz(spk_ang,posang,post,h,mapAxis,vfs)
invh = 1/h; % h > 0 is a smoothing parameter called the bandwidth
map = zeros(length(mapAxis),1);
idx = 0;
for ang0 = mapAxis   
    idx = idx + 1;
    map(idx) = rate_estimator_ang_cz(spk_ang,ang0,invh,posang,post,vfs);    
end

% Calculate the rate for one position value
function r = rate_estimator_ang_cz(spk_ang,ang0,invh,posang,post,vfs)
% edge-corrected kernel density estimator
n_post=length(post);
post0=1/vfs:1/vfs:n_post/vfs;

delta = min([abs(spk_ang - ang0),2*pi-abs(spk_ang - ang0)],[],2);
conv_sum = sum(gaussian_kernel(delta*invh));
delta = min([abs(posang - ang0);2*pi-abs(posang - ang0)],[],1);
edge_corrector =  trapz(post0,gaussian_kernel((delta*invh)));
r = (conv_sum / (edge_corrector + 0.01)) + 0.01; % regularised firing rate for "wellbehavedness"
                                                       % i.e. no division by zero or log of zero
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x)
% k(u) = ((2*pi)^(-length(u)/2)) * exp(u'*u)
r = 0.15915494309190 * exp(-0.5*(x.*x));