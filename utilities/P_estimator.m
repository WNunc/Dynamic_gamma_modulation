% Estimate the mean P for one position value
function p = P_estimator(P,vel,vbin,invh)
% edge-corrected kernel density estimator
num = P*gaussian_kernel((vel-vbin),invh);
den = sum(gaussian_kernel((vel-vbin),invh));
p = num/den;
                                                    
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,invh)
r =  invh/sqrt(2*pi)*exp(-0.5*x.*x*invh^2);