% ####################### circular firing rate map #########################
% Calculates the rate map.
function map = ratemap_ang_circ_cz(spk_ang,posang,post,h,mapAxis,vfs)
invh = 1/h; % h > 0 is a smoothing parameter called the bandwidth
map = zeros(length(mapAxis),1);
idx = 0;
for ang0 = mapAxis   
    idx = idx + 1;
    map(idx) = rate_estimator_ang_circ_cz(spk_ang,ang0,invh,posang,post,vfs);    
end

% Calculate the rate for one position value
function r = rate_estimator_ang_circ_cz(spk_ang,ang0,invh,posang,post,vfs)
% edge-corrected kernel density estimator
n_post=length(post);
post0=1/vfs:1/vfs:n_post/vfs;
conv_sum = sum(gaussian_kernel_circ(spk_ang-ang0,invh));
edge_corrector =  trapz(post0,gaussian_kernel_circ(posang-ang0,invh));
r = (conv_sum / (edge_corrector + 0.01)) + 0.01; % regularised firing rate for "wellbehavedness"
                                                       % i.e. no division by zero or log of zero
% Gaussian kernel for the rate calculation
function r = gaussian_kernel_circ(x,invh)
% k(u) = ((2*pi)^(-length(u)/2)) * exp(u'*u)
r = exp(10*invh*cos(x));
% Large values of invh lead to highly
% variable estimations whereas small values provide
% oversmoothed circular densities (Oliveira, et al. 2012)
% Oliveira, et al. 2012: A plug-in rule for bandwidth selection in circular density estimation