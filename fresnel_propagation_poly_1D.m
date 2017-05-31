%% fresnel_propagation_poly_1D.m
%
% DESCRIPTION: polychromatic fresnel propagation for 1D waves
% 
% CALL: Do = fresnel_propagation_poly_1D(Di,L,lambda,z)
%   - Di: intput wave, DIM 2 spatial dimension, DIM 1 energy dimension  
%   - L : field of view
%   - lambda: vector containing all the wavelegths for which the
%   propagation is performed
%
%   - z : propagation distance
%   - Do: output wave
%
%
% UPDATES:
%   30.09.2013 (Matias) : fixed dimesion problem
%   30.09.2013 (Matias) : first version
%
%%
function Do = fresnel_propagation_poly_1D(Di,L,lambda,z) 
     
     M = size(Di,2); %get input field array size
     dx = L/M; %sample interval
     fx = -1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords
     fx = repmat(fx,length(lambda),1);
     lambda = repmat(lambda',1,length(Di));
     H = exp(-1i*pi*lambda*z.*fx.^2); %trans func
     H = fftshift(H,2); %shift trans func
     di = fft(fftshift(Di,2),[],2); %shift, fft src fields
     do = H.*di; %multiply
     Do = ifftshift(ifft(do,[],2),2); %inv fft, center obs field
     
 end